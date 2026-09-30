#' Pipeline ID lineage table
#'
#' @description
#' A normalized lineage table: one row per leaf NHD flowline (including
#' post-refactor split children), with columns tracking the ID progression
#' through each pipeline stage. Replaces comma-separated `member_comid`
#' string parsing with data.table joins.
#'
#' Schema:
#' - `nhd_comid`      (character): original NHD COMID, no suffix
#' - `split_part`     (integer)  : 0 if unsplit; 1, 2, ... for post-split children
#' - `refactored_id`  (character): e.g. `"85000100004325.1"` or `"85000100004325"`
#' - `reconciled_id`  (integer)  : reconciled group ID after refactor_flowlines
#' - `aggregated_id`  (character): set by `lineage_attach_aggregated()`
#' - `final_divide_id`(character): set by `lineage_attach_ngen()`
#' - `final_flowpath_id`(character): set by `lineage_attach_ngen()`
#'
#' @name lineage
#' @keywords internal
NULL

#' Build the initial lineage table from refactor_flowlines output
#'
#' @param refactored An sf (or data frame) of refactored flowlines with a
#'   `flowpath_id` column containing base COMIDs and `.N` split suffixes.
#' @param reconciled An sf with columns `ID` (integer) and
#'   `member_flowpath_id` (comma-separated refactored IDs).
#'
#' @return A data.table with one row per refactored flowline, sorted by
#'   `reconciled_id` then `refactored_id`.
#'
#' @examples
#' \dontrun{
#' refactored <- sf::st_read("refactored.gpkg", "refactored")
#' reconciled <- sf::st_read("reconciled.gpkg", "reconciled")
#' lin <- lineage_build(refactored, reconciled)
#' }
#'
#' @export
lineage_build <- function(refactored, reconciled) {
  rf <- data.table::as.data.table(sf::st_drop_geometry(refactored))
  rf <- rf[, .(refactored_id = as.character(flowpath_id))]
  rf[, nhd_comid  := sub("\\.[0-9]+$", "", refactored_id)]
  has_part <- grepl("\\.[0-9]+$", rf$refactored_id)
  rf[, split_part := 0L]
  if (any(has_part))
    rf[has_part, split_part := suppressWarnings(as.integer(
      sub(".*\\.", "", refactored_id)))]

  # Expand reconciled members: one row per (reconciled_id, refactored member)
  rec_dt <- data.table::as.data.table(sf::st_drop_geometry(reconciled))
  .id_col <- if ("rec_id" %in% names(rec_dt)) "rec_id" else "ID"
  rec_dt <- rec_dt[, .(reconciled_id = as.integer(get(.id_col)),
    mem = as.character(member_flowpath_id))]
  splits <- strsplit(rec_dt$mem, ",", fixed = TRUE)
  rec_long <- data.table::data.table(
    reconciled_id = rep(rec_dt$reconciled_id, lengths(splits)),
    refactored_id = trimws(unlist(splits))
  )
  rec_long <- rec_long[nzchar(refactored_id)]

  # Join: every refactored row gets its reconciled_id; detect any orphan
  out <- merge(rf, rec_long, by = "refactored_id", all.x = TRUE, sort = FALSE)
  n_orphan <- sum(is.na(out$reconciled_id))
  if (n_orphan > 0L) .hf_log_warn(sprintf(
    "lineage_build: %d refactored flowline(s) not in any reconciled group",
    n_orphan))

  out[, aggregated_id     := NA_character_]
  out[, final_divide_id   := NA_character_]
  out[, final_flowpath_id := NA_character_]

  data.table::setkeyv(out, c("reconciled_id", "refactored_id"))
  data.table::setcolorder(out, c("nhd_comid", "split_part", "refactored_id",
    "reconciled_id", "aggregated_id",
    "final_divide_id", "final_flowpath_id"))
  out
}

#' Attach aggregated flowpath IDs to a lineage table
#'
#' @param lineage A lineage `data.table` from `lineage_build()`.
#' @param network A data.frame or sf object with columns `flowpath_id`
#'   (aggregated, integer or character) and `rec_id_list` (comma-separated
#'   reconciled IDs): typically `sf::st_read(out_gpkg, "network")`.
#'
#' @return The lineage table with `aggregated_id` populated (in place).
#'
#' @examples
#' \dontrun{
#' lin     <- lineage_build(refactored, reconciled)
#' network <- sf::st_read("out_vpu01.gpkg", "network")
#' lin <- lineage_attach_aggregated(lin, network)
#' }
#'
#' @export
lineage_attach_aggregated <- function(lineage, network) {
  stopifnot(data.table::is.data.table(lineage))
  nw <- data.table::as.data.table(
    as.data.frame(network)[, c("flowpath_id", "rec_id_list")])
  nw[, flowpath_id := as.character(flowpath_id)]
  nw[, rec_id_list := as.character(rec_id_list)]
  nw <- unique(nw, by = c("flowpath_id", "rec_id_list"))
  splits <- strsplit(nw$rec_id_list, ",", fixed = TRUE)
  long <- data.table::data.table(
    aggregated_id = rep(nw$flowpath_id, lengths(splits)),
    reconciled_id = as.integer(trimws(unlist(splits)))
  )
  long <- unique(long[!is.na(reconciled_id)], by = "reconciled_id")
  lineage[, aggregated_id := NULL]
  out <- long[lineage, on = "reconciled_id"]
  data.table::setnames(out, "aggregated_id", "aggregated_id")
  # Reorder columns back
  data.table::setcolorder(out, c("nhd_comid", "split_part", "refactored_id",
    "reconciled_id", "aggregated_id",
    "final_divide_id", "final_flowpath_id"))
  out
}

#' Attach NextGen IDs (`fp-N`, `cat-N`) to a lineage table
#'
#' `nextgen_from_gpkg` prefixes aggregated IDs with `"fp-"` and `"cat-"` deterministically,
#' so this attach is a trivial string paste, no join to the network table
#' required.
#'
#' @param lineage A lineage `data.table` that already has `aggregated_id`.
#' @param flowpath_prefix Prefix for final flowpath IDs. Default `"fp-"`.
#' @param divide_prefix   Prefix for final divide IDs. Default `"cat-"`.
#'
#' @return Lineage with `final_flowpath_id` and `final_divide_id` populated.
#'
#' @examples
#' \dontrun{
#' lin <- lineage_attach_aggregated(lin, network)
#' lin <- lineage_attach_ngen(lin)
#' }
#'
#' @export
lineage_attach_ngen <- function(lineage,
                                flowpath_prefix = "fp-",
                                divide_prefix   = "cat-") {
  stopifnot(data.table::is.data.table(lineage))
  has_agg <- !is.na(lineage$aggregated_id)
  lineage[, final_flowpath_id := NA_character_]
  lineage[, final_divide_id   := NA_character_]
  lineage[has_agg, final_flowpath_id := paste0(flowpath_prefix, aggregated_id)]
  lineage[has_agg, final_divide_id   := paste0(divide_prefix,   aggregated_id)]
  lineage
}

#' Recompute `member_comid` CSV strings from a lineage table
#'
#' Utility to validate that the lineage table carries the same information
#' as the existing `member_comid` strings.
#'
#' @param lineage A lineage data.table.
#' @param group_by Column to group by. Default `"aggregated_id"` (NextGen /
#'   aggregated-flowpath member lists).
#'
#' @return data.table with columns `group_by` and `member_comid` (CSV).
#'
#' @examples
#' lin <- data.table::data.table(
#'   refactored_id = c("101", "102", "103"),
#'   aggregated_id = c("fp-1", "fp-1", "fp-2"))
#' lineage_member_csv(lin, group_by = "aggregated_id")
#'
#' @export
lineage_member_csv <- function(lineage, group_by = "aggregated_id") {
  stopifnot(data.table::is.data.table(lineage))
  if (!(group_by %in% names(lineage)))
    stop(sprintf("column '%s' not in lineage", group_by))
  lineage[!is.na(get(group_by)),
    .(member_comid = paste(refactored_id, collapse = ",")),
    by = c(group_by)]
}
