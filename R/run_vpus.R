# A build that finishes is not a build that passed.
#
# run_vpus only ever recorded a VPU as failed when run_pipeline threw. The
# invariant suite runs inside the build and writes its verdict to
# invariants.json, but nothing read it back, so a VPU whose gates went red was
# reported "ok" and fed straight into the merge, which is how an Alaska build
# carrying 2429 non-integer poi_id values came to be counted among "3 ok".
.rv_invariant_verdict <- function(ngen_gpkg) {
  if (is.null(ngen_gpkg) || is.na(ngen_gpkg)) return(NULL)
  f <- file.path(dirname(ngen_gpkg), "work", "invariants.json")
  if (!file.exists(f)) return(NULL)
  j <- tryCatch(jsonlite::fromJSON(f, simplifyVector = FALSE),
    error = function(e) NULL)
  if (is.null(j) || is.null(j$checks)) return(NULL)
  # An entry with no `kind` is a gate, not something to skip. The merge suite
  # writes its structural checks as bare list(ok = ...), so filtering on
  # kind == "check" silently ignored 12 of 21 entries: is_dag and the
  # duplicate-id checks among them. Absence of a level is not permission to
  # disregard a result, for the same reason absence of an input is not a pass.
  failed <- names(Filter(function(e)
    !identical(e$kind, "info") && identical(e$ok, FALSE), j$checks))
  list(ok = isTRUE(j$ok) && !length(failed), failed = failed)
}

#' Build a set of VPUs, then merge them
#'
#' @description
#' A loop over [run_pipeline()] with per-VPU inputs, so a national build is a
#' list of paths rather than a script. One VPU failing does not stop the rest by
#' default; the return value says which succeeded, and [merge_vpus()] can be run
#' from it once they all have.
#'
#' `inputs` is a named list, one entry per `vpuid`, each holding the arguments
#' that differ for that VPU. Anything passed in `...` applies to every VPU, so
#' shared paths are given once.
#'
#' @param inputs Named list of per-VPU argument lists, keyed by `vpuid`.
#' @param out_dir Parent output directory; each VPU writes `vpu_{vpuid}/`.
#' @param ... Arguments common to every VPU, forwarded to [run_pipeline()].
#' @param continue_on_error If `TRUE` (default), record a failure and carry on.
#'   `FALSE` stops at the first one.
#' @param merge If `TRUE`, merge the VPUs that succeeded once the loop finishes.
#'   Default `FALSE`.
#' @param merge_dir Directory for the merged domain products when `merge = TRUE`.
#'   Defaults to `file.path(out_dir, "global")`.
#'
#' @returns A data frame with one row per VPU: `vpuid`, `status`
#'   (`"ok"`/`"error"`/`"invariant"`/`"skipped"`), `ngen_gpkg`, `minutes`, and
#'   `message`. `"invariant"` means the build finished but its gates went red;
#'   it counts as a failure and blocks the merge, because a product that does
#'   not pass its own suite is not one to publish.
#'   Returned invisibly when every VPU succeeded.
#'
#' @examples
#' \dontrun{
#' conus <- c("01","02","03N","03S","03W","04","05","06","07","08","09",
#'            "10L","10U","11","12","13","14","15","16","17","18")
#'
#' # Shared inputs once, per-VPU overrides where they differ.
#' inputs <- stats::setNames(vector("list", length(conus)), conus)
#' inputs[["20"]] <- list(reference_gpkg = "hi_reference.gpkg",
#'                        fdr_path = "hi/hi_fdr.vrt", fac_path = "hi/hi_fac.vrt")
#'
#' res <- run_vpus(
#'   inputs,
#'   out_dir        = "out",
#'   reference_gpkg = "conus_reference.gpkg",
#'   fdr_path       = "conus/fdr.vrt",
#'   fac_path       = "conus/fac.vrt",
#'   skip_if_exists = TRUE,
#'   merge          = TRUE)
#'
#' subset(res, status != "ok")
#' }
#' @export
run_vpus <- function(inputs, out_dir, ..., continue_on_error = TRUE,
                     merge = FALSE, merge_dir = NULL) {
  if (is.character(inputs))
    inputs <- stats::setNames(vector("list", length(inputs)), inputs)
  stopifnot(is.list(inputs), length(inputs) > 0L, !is.null(names(inputs)),
    all(nzchar(names(inputs))))
  common <- list(...)

  rows <- vector("list", length(inputs))
  for (i in seq_along(inputs)) {
    v  <- names(inputs)[i]
    t0 <- proc.time()[["elapsed"]]
    message(sprintf("[run_vpus] %s (%d of %d)", v, i, length(inputs)))
    # Per-VPU values win over the shared ones, so a domain overrides only what
    # differs for it.
    args <- utils::modifyList(common, inputs[[i]] %||% list())
    args$vpuid   <- v
    args$out_dir <- out_dir
    res <- tryCatch(do.call(run_pipeline, args),
      error = function(e) structure(list(message = conditionMessage(e)),
        class = "hf_vpu_error"))
    mins <- round((proc.time()[["elapsed"]] - t0) / 60, 1)
    if (inherits(res, "hf_vpu_error")) {
      .hf_log_warn(sprintf("[run_vpus] %s failed: %s", v, res$message))
      if (!continue_on_error) stop(sprintf("run_vpus: %s failed: %s", v, res$message),
        call. = FALSE)
      rows[[i]] <- data.frame(vpuid = v, status = "error", ngen_gpkg = NA_character_,
        minutes = mins, message = res$message, stringsAsFactors = FALSE)
    } else {
      g   <- res$ngen_gpkg %||% NA_character_
      inv <- .rv_invariant_verdict(g)
      st  <- if (isTRUE(res$skipped)) "skipped" else "ok"
      msg <- NA_character_
      if (!is.null(inv) && !inv$ok) {
        st  <- "invariant"
        msg <- sprintf("%d failed gate(s): %s", length(inv$failed),
          paste(inv$failed, collapse = ", "))
        .hf_log_warn(sprintf("[run_vpus] %s built but did not pass: %s", v, msg))
      } else if (is.null(inv) && st == "ok") {
        .hf_log_warn(sprintf(
          "[run_vpus] %s wrote no invariants.json; its gates are unverified", v))
      }
      rows[[i]] <- data.frame(vpuid = v, status = st, ngen_gpkg = g,
        minutes = mins, message = msg, stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)

  n_bad <- sum(out$status %in% c("error", "invariant"))
  message(sprintf(
    "[run_vpus] %d ok, %d skipped, %d failed, %d built but did not pass (%.1f min total)",
    sum(out$status == "ok"), sum(out$status == "skipped"),
    sum(out$status == "error"), sum(out$status == "invariant"), sum(out$minutes)))

  if (isTRUE(merge)) {
    if (n_bad > 0L)
      stop(sprintf(paste("run_vpus: refusing to merge with %d failed VPU(s) (%s).",
        "The global id pool is assigned across the whole set, so merging a",
        "partial set produces ids that a later full run will not reproduce."),
        n_bad, paste(out$vpuid[out$status %in% c("error", "invariant")],
          collapse = ", ")),
        call. = FALSE)
    built <- out[!is.na(out$ngen_gpkg) & out$status != "invariant", ]
    merge_vpus(gpkg_map = stats::setNames(as.list(built$ngen_gpkg), built$vpuid),
      out_dir = merge_dir %||% file.path(out_dir, "global"))
  }
  if (n_bad == 0L) invisible(out) else out
}
