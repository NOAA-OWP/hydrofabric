# Recording what a build produced, next to what it consumed.
#
# hf_input_manifest() pins the inputs; nothing pinned the outputs. A recipient
# could not verify a product arrived intact, and after the fact there was no
# record tying a delivered GeoPackage to the commit and the input set that made
# it. The provenance stamped inside each GeoPackage names the commit but not the
# inputs, and says nothing a consumer can check a download against.

.om_layer_counts <- function(gpkg) {
  con <- DBI::dbConnect(RSQLite::SQLite(), gpkg)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  lyr <- tryCatch(DBI::dbGetQuery(con,
    "SELECT table_name FROM gpkg_contents ORDER BY table_name")$table_name,
    error = function(e) character(0))
  stats::setNames(lapply(lyr, function(l) tryCatch(
    DBI::dbGetQuery(con, sprintf('SELECT COUNT(*) n FROM "%s"', l))$n,
    error = function(e) NA_integer_)), lyr)
}

.om_invariants <- function(gpkg) {
  d <- dirname(gpkg)
  dom <- sub("_ngen\\.gpkg$", "", basename(gpkg))
  cand <- c(file.path(d, sprintf("invariants_%s.json", dom)),
            file.path(d, "work", "invariants.json"))
  f <- cand[file.exists(cand)]
  if (!length(f)) return(list(record = NA_character_, ok = NA, failing = character(0)))
  j <- tryCatch(jsonlite::fromJSON(f[1], simplifyVector = FALSE),
    error = function(e) NULL)
  if (is.null(j)) return(list(record = basename(f[1]), ok = NA, failing = character(0)))
  # A missing `kind` is a gate; only an explicit "info" is exempt.
  bad <- names(Filter(function(e)
    !identical(e$kind, "info") && identical(e$ok, FALSE), j$checks %||% list()))
  list(record = basename(f[1]), ok = isTRUE(j$ok) && !length(bad), failing = bad)
}

#' Record what a build produced
#'
#' @description
#' Writes a manifest of finished products: the md5, size, layer counts and
#' invariant verdict of each, together with the commit that built them and the
#' digest of the input manifest they were built from. Verifiable with
#' [hf_verify_manifest()], the same way inputs are.
#'
#' @details
#' This is the other half of [hf_input_manifest()]. Without it a product cannot
#' be checked against what was shipped, and the tie between a delivered file, the
#' code that made it and the inputs it consumed lives only in a build log.
#'
#' The input manifest is recorded by digest rather than copied, so one line
#' establishes which input set was used without duplicating it.
#'
#' @param gpkgs Paths to the finished products.
#' @param out Path to write the manifest JSON to.
#' @param input_manifest Optional path to the input manifest used for the build;
#'   its md5 is recorded so the two can be tied together.
#' @param git_sha Commit that produced the products. Defaults to `HF_GIT_SHA`.
#' @returns Invisibly, the manifest as a list. Written to `out` as JSON.
#' @examples
#' \dontrun{
#' hf_output_manifest(
#'   Sys.glob("global/*_ngen.gpkg"),
#'   out = "global/output_manifest.json",
#'   input_manifest = "owp-inputs/national_inputs_manifest.json")
#' }
#' @export
hf_output_manifest <- function(gpkgs, out, input_manifest = NULL,
                               git_sha = Sys.getenv("HF_GIT_SHA", "")) {
  gpkgs <- gpkgs[file.exists(gpkgs)]
  if (!length(gpkgs)) stop("no existing product paths given", call. = FALSE)

  products <- lapply(gpkgs, function(g) {
    inv <- .om_invariants(g)
    message(sprintf("[manifest:out] %-24s %s  invariants=%s", basename(g),
      format(structure(file.size(g), class = "object_size"), units = "auto"),
      if (is.na(inv$ok)) "no record" else if (inv$ok) "ok" else
        paste0("FAIL(", paste(inv$failing, collapse = ","), ")")))
    list(name = basename(g), path = normalizePath(g),
         bytes = file.size(g), md5 = unname(tools::md5sum(g)),
         layers = .om_layer_counts(g),
         invariants_record = inv$record, invariants_ok = inv$ok,
         invariants_failing = inv$failing)
  })

  man <- list(
    created    = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    git_sha    = if (nzchar(git_sha)) git_sha else NULL,
    software   = paste("hydrofabric", as.character(utils::packageVersion("hydrofabric"))),
    input_manifest = if (!is.null(input_manifest) && file.exists(input_manifest))
      list(path = normalizePath(input_manifest),
           md5  = unname(tools::md5sum(input_manifest))) else NULL,
    products   = stats::setNames(products, vapply(products, `[[`, "", "name")))

  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  writeLines(jsonlite::toJSON(man, auto_unbox = TRUE, pretty = TRUE), out)
  message(sprintf("[manifest:out] %d product(s) -> %s", length(products), out))
  invisible(man)
}

#' Verify products against an output manifest
#'
#' @param manifest Path to a manifest written by [hf_output_manifest()].
#' @param strict If `TRUE`, raise an error when a product is missing or altered.
#' @returns Invisibly, a data frame with `name`, `exists`, `md5_ok`.
#' @examples
#' \dontrun{
#' hf_verify_output_manifest("global/output_manifest.json", strict = TRUE)
#' }
#' @export
hf_verify_output_manifest <- function(manifest, strict = FALSE) {
  m <- jsonlite::fromJSON(manifest, simplifyVector = FALSE)
  rows <- lapply(m$products, function(p) {
    ex <- file.exists(p$path)
    ok <- ex && identical(unname(tools::md5sum(p$path)), p$md5)
    message(sprintf("[verify:out] %-8s %-24s %s", if (ok) "ok" else "MISMATCH",
      p$name, if (!ex) "missing" else if (ok) "" else "content changed"))
    data.frame(name = p$name, exists = ex, md5_ok = ok, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  bad <- sum(!out$md5_ok)
  msg <- sprintf("hf_verify_output_manifest: %d of %d product(s) do not match",
    bad, nrow(out))
  if (bad && strict) stop(msg, call. = FALSE)
  if (bad) .hf_log_warn(msg) else
    message(sprintf("[verify:out] all %d product(s) match", nrow(out)))
  invisible(out)
}
