#' Add a nested-set upstream index to a hydrofabric GeoPackage
#'
#' Stamps `upstream_id` (depth-first pre-order ordinal from each outlet) and
#' `num_upstreams` (strict upstream count) on the flowpaths, flowlines, divides,
#' and network layers, so the flowpaths upstream of X are those with
#' `upstream_id` in `(upstream_id[X], upstream_id[X] + num_upstreams[X]]`, an
#' O(1) range filter instead of a traversal.
#'
#' The values are specific to this file's network (a per-VPU, merged, or subset
#' gpkg each gets a correct-scope index) and are not a persistent key. The
#' compute (fp -> nexus -> fp resolution + nested set) is
#' [hfutils::hf_upstream_index()]; reads and writes are attribute-only via DBI,
#' so it is safe at merged/national scale.
#'
#' @param ngen_gpkg Path to the NextGen GeoPackage.
#' @param flowpaths_layer Flowpath layer name. Default `"flowpaths"`.
#' @param nexus_layer Nexus layer name. Default `"nexus"`.
#' @param propagate_to Layers to stamp by `flowpath_id`. Default
#'   `c("flowlines", "divides", "network")`.
#'
#' @return Invisibly, a data.frame of `flowpath_id`, `upstream_id`, and
#'   `num_upstreams`. Adds those integer columns to the layers in place.
#' @examples
#' \dontrun{
#' add_upstream_index("ngen_vpu01.gpkg")
#' # flowpaths strictly upstream of fp-8:
#' library(DBI)
#' con <- dbConnect(RSQLite::SQLite(), "ngen_vpu01.gpkg")
#' anchor <- dbGetQuery(con, "SELECT upstream_id u, num_upstreams n
#'                            FROM flowpaths WHERE flowpath_id = '8'")
#' dbGetQuery(con, sprintf(
#'   "SELECT flowpath_id FROM flowpaths WHERE upstream_id > %d AND upstream_id <= %d",
#'   anchor$u, anchor$u + anchor$n))
#' }
#' @export
add_upstream_index <- function(ngen_gpkg,
                               flowpaths_layer = "flowpaths",
                               nexus_layer     = "nexus",
                               propagate_to    = c("flowlines", "divides", "network")) {
  if (!hfutils::layer_exists(ngen_gpkg, flowpaths_layer))
    stop("add_upstream_index: no '", flowpaths_layer, "' layer in ", ngen_gpkg,
      call. = FALSE)

  con <- DBI::dbConnect(RSQLite::SQLite(), ngen_gpkg)
  on.exit(if (DBI::dbIsValid(con)) DBI::dbDisconnect(con), add = TRUE)

  # Attribute-only reads (no geometry): safe at national/merged scale.
  fp  <- DBI::dbGetQuery(con, sprintf(
    'SELECT flowpath_id, flowpath_toid FROM "%s"', flowpaths_layer))
  nex <- if (hfutils::layer_exists(ngen_gpkg, nexus_layer))
    DBI::dbGetQuery(con, sprintf(
      'SELECT nexus_id, nexus_toid FROM "%s"', nexus_layer)) else NULL

  # The index is only defined on an acyclic network; degrade to a skip (rather
  # than aborting finalize) if a cycle is ever present.
  idx <- tryCatch(hfutils::hf_upstream_index(fp, nex), error = function(e) {
    .hf_log_warn("add_upstream_index: skipped -- ", conditionMessage(e))
    NULL
  })
  if (is.null(idx) || nrow(idx) == 0L) return(invisible(idx))

  # Write back via DBI UPDATE to every layer keyed by flowpath_id (geometry
  # untouched). Triggers are dropped/restored around each update, matching the
  # drainage-area finalize, so gpkg rtree triggers never fire on the write.
  bf <- data.frame(flowpath_id = as.character(idx$flowpath_id),
    upstream_id = as.integer(idx$upstream_id),
    num_upstreams = as.integer(idx$num_upstreams), stringsAsFactors = FALSE)
  DBI::dbWriteTable(con, "ui_tmp", bf, temporary = TRUE, overwrite = TRUE)
  DBI::dbExecute(con, "CREATE INDEX ui_idx ON ui_tmp(flowpath_id)")

  for (lyr in unique(c(flowpaths_layer, propagate_to))) {
    if (!DBI::dbExistsTable(con, lyr)) next
    cols <- DBI::dbListFields(con, lyr)
    if (!("flowpath_id" %in% cols)) next
    for (c_ in c("upstream_id", "num_upstreams"))
      if (!(c_ %in% cols))
        DBI::dbExecute(con, sprintf('ALTER TABLE "%s" ADD COLUMN "%s" INTEGER', lyr, c_))
    trg <- DBI::dbGetQuery(con, sprintf(
      "SELECT name, sql FROM sqlite_master WHERE type='trigger' AND tbl_name='%s'", lyr))
    for (nm in trg$name) DBI::dbExecute(con, sprintf('DROP TRIGGER IF EXISTS "%s"', nm))
    DBI::dbExecute(con, sprintf(
      'UPDATE "%s" SET upstream_id=(SELECT u.upstream_id FROM ui_tmp u WHERE u.flowpath_id="%s".flowpath_id),
                       num_upstreams=(SELECT u.num_upstreams FROM ui_tmp u WHERE u.flowpath_id="%s".flowpath_id)
       WHERE flowpath_id IN (SELECT flowpath_id FROM ui_tmp)', lyr, lyr, lyr))
    for (s in trg$sql) if (!is.na(s) && nzchar(s))
      tryCatch(DBI::dbExecute(con, s), error = function(e) NULL)
  }
  DBI::dbExecute(con, "DROP TABLE ui_tmp")
  try(DBI::dbExecute(con, sprintf(
    'CREATE INDEX IF NOT EXISTS idx_%s_upstream_id ON "%s"(upstream_id)',
    flowpaths_layer, flowpaths_layer)), silent = TRUE)

  message(sprintf(
    "add_upstream_index: indexed %d flowpath(s) across %d outlet(s); max num_upstreams = %d",
    nrow(idx), attr(idx, "n_outlets"), max(idx$num_upstreams)))
  invisible(idx)
}
