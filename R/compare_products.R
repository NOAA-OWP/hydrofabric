# Comparing a rebuild against the product it replaces.
#
# Every invariant is intra-build: it asks whether one fabric is self-consistent,
# not whether it still says what the last one said. That leaves a class of
# regression invisible. A flowpath that loses its divide and is relabelled a
# connector satisfies `flowpath_divide_assigned`, because connectors are
# legitimately divide-less and the count of "explained" gaps simply rises. The
# orphaned area is absorbed by neighbours, so conservation holds to four decimal
# places. Both times this happened during the OWP delivery the suite passed and
# only a direct diff of the two GeoPackages showed it.
#
# `connector_claim_valid` does detect it, 36 of 109 became 44 of 117, but it
# is informational, and it cannot be gated at zero while a clean build already
# reports 36. Comparing against the previous product is what turns that signal
# into a decision.

.cp_counts <- function(con, layers) {
  vapply(layers, function(l) tryCatch(
    DBI::dbGetQuery(con, sprintf('SELECT COUNT(*) n FROM "%s"', l))$n,
    error = function(e) NA_integer_), numeric(1L))
}

.cp_divide_owners <- function(f) {
  con <- DBI::dbConnect(RSQLite::SQLite(), f)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  tryCatch(DBI::dbGetQuery(con,
    "SELECT flowpath_id, SUM(areasqkm) AS areasqkm FROM divides
      WHERE flowpath_id IS NOT NULL GROUP BY flowpath_id"),
    error = function(e) NULL)
}

.cp_flowpath_ids <- function(f) {
  con <- DBI::dbConnect(RSQLite::SQLite(), f)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  tryCatch(as.character(DBI::dbGetQuery(con,
    "SELECT flowpath_id FROM flowpaths")$flowpath_id), error = function(e) NULL)
}

#' Compare a rebuilt product against the one it replaces
#'
#' @description
#' Diffs two NextGen GeoPackages on the things an intra-build invariant cannot
#' see: whether features were lost, and whether a flowpath that owned a divide
#' still owns one. A flowpath silently reclassified as a connector passes every
#' check in the suite, so a rebuild that drops divides can otherwise look clean.
#'
#' @section What this can and cannot decide:
#' A difference is only evidence of a regression when the inputs were the same.
#' Hold them constant and any loss is attributable to the code: that is the case
#' this was written for, and the one where `strict = TRUE` belongs.
#'
#' Across an input change it reports rather than decides. Hydrolocations force
#' reach splits, so a different hydrolocation set legitimately adds, removes and
#' resizes divides; a divide whose POI is no longer in the set is *supposed* to
#' disappear. Used as a gate there it fails a correct build, so read the output
#' and account for the differences instead.
#'
#' With both changed at once nothing is attributable. When a new input set and new
#' code both need to land, build one unit on the new inputs with the current code
#' first and compare that, one VPU is a couple of hours against a national run's
#' tens of hours, and it buys an interpretable answer.
#'
#' @param old_gpkg,new_gpkg Paths to the previous and rebuilt products.
#' @param area_tol Relative change in a flowpath's divide area to report.
#'   Default `0.01` (1%).
#' @param strict If `TRUE`, raise an error when a flowpath loses its divide or a
#'   layer loses features. Default `FALSE` reports and returns.
#' @returns Invisibly, a list with `ok`, `counts`, `lost_divide`,
#'   `gained_divide` and `area_changed`.
#' @examples
#' \dontrun{
#' # Same inputs, changed code: a loss is a regression, so gate on it.
#' hf_compare_products("before/ngen_vpu04.gpkg", "after/ngen_vpu04.gpkg",
#'                     strict = TRUE)
#'
#' # New hydrolocation set: report and account for the differences.
#' hf_compare_products("baseline/ngen_vpu04.gpkg", "rebuild/ngen_vpu04.gpkg")
#' }
#' @export
hf_compare_products <- function(old_gpkg, new_gpkg, area_tol = 0.01,
                                strict = FALSE) {
  stopifnot(file.exists(old_gpkg), file.exists(new_gpkg))
  layers <- c("flowpaths", "divides", "nexus", "flowlines", "pois",
              "hydrolocations", "lakes")

  co <- DBI::dbConnect(RSQLite::SQLite(), old_gpkg)
  cn <- DBI::dbConnect(RSQLite::SQLite(), new_gpkg)
  counts <- data.frame(layer = layers,
    old = .cp_counts(co, layers), new = .cp_counts(cn, layers),
    stringsAsFactors = FALSE)
  DBI::dbDisconnect(co); DBI::dbDisconnect(cn)
  counts$delta <- counts$new - counts$old
  counts <- counts[!is.na(counts$old) | !is.na(counts$new), ]
  rownames(counts) <- NULL

  o <- .cp_divide_owners(old_gpkg)
  n <- .cp_divide_owners(new_gpkg)
  # Only flowpaths present in both builds can be compared. Identifiers are
  # minted per build from the aggregation, so a changed POI set re-mints most of
  # them: on VPU 01 only 7,010 of 15,345 ids survived a hydrolocation change.
  # Joining on the raw identifier counted every re-minted reach as a divide
  # lost, reporting 8,351 where 46 flowpaths had actually lost one.
  common <- intersect(.cp_flowpath_ids(old_gpkg) %||% character(0),
                      .cp_flowpath_ids(new_gpkg) %||% character(0))
  churn  <- length(unique(c(.cp_flowpath_ids(old_gpkg),
                            .cp_flowpath_ids(new_gpkg)))) - length(common)
  lost <- gained <- changed <- data.frame()
  if (!is.null(o) && !is.null(n)) {
    o <- o[as.character(o$flowpath_id) %in% common, , drop = FALSE]
    n <- n[as.character(n$flowpath_id) %in% common, , drop = FALSE]
    m <- merge(o, n, by = "flowpath_id", all = TRUE,
               suffixes = c("_old", "_new"))
    lost   <- m[!is.na(m$areasqkm_old) & is.na(m$areasqkm_new), ]
    gained <- m[is.na(m$areasqkm_old) & !is.na(m$areasqkm_new), ]
    both   <- m[!is.na(m$areasqkm_old) & !is.na(m$areasqkm_new), ]
    if (nrow(both)) {
      rel <- abs(both$areasqkm_new - both$areasqkm_old) /
        pmax(both$areasqkm_old, .Machine$double.eps)
      changed <- both[rel > area_tol, ]
      changed <- changed[order(-abs(changed$areasqkm_new - changed$areasqkm_old)), ]
    }
  }

  for (i in seq_len(nrow(counts)))
    message(sprintf("[compare] %-15s %8s -> %-8s %+d",
      counts$layer[i], counts$old[i], counts$new[i], counts$delta[i]))

  message(sprintf("[compare] %d flowpath id(s) common to both; %d re-minted or %s",
    length(common), churn, "added"))
  if (nrow(lost))
    message(sprintf("[compare] %d common flowpath(s) LOST their divide: %s", nrow(lost),
      paste(utils::head(lost$flowpath_id, 8), collapse = ", ")))
  else message("[compare] no common flowpath lost its divide")
  if (nrow(gained))
    message(sprintf("[compare] %d flowpath(s) gained a divide", nrow(gained)))
  if (nrow(changed))
    message(sprintf("[compare] %d divide area(s) changed by more than %.1f%% (largest: %s)",
      nrow(changed), 100 * area_tol, changed$flowpath_id[1L]))

  shrunk <- counts$layer[!is.na(counts$delta) & counts$delta < 0]
  ok <- nrow(lost) == 0L && length(shrunk) == 0L
  if (!ok) {
    msg <- sprintf(paste("hf_compare_products: %d flowpath(s) lost a divide;",
      "layer(s) with fewer features: %s"), nrow(lost),
      if (length(shrunk)) paste(shrunk, collapse = ", ") else "none")
    if (strict) stop(msg, call. = FALSE) else .hf_log_warn(msg)
  } else {
    message("[compare] no regression: nothing lost")
  }
  invisible(list(ok = ok, counts = counts, lost_divide = lost,
                 gained_divide = gained, area_changed = changed))
}
