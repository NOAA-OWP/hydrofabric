# Internal helpers shared across partition_reconcile_divides and
# reconcile_divides. These are package-private (dot-prefixed); not exported.

# Coerce named columns to character without touching others.
.hf_as_character_cols <- function(df, cols) {
  for (nm in cols) {
    if (nm %in% names(df))
      df[[nm]] <- as.character(df[[nm]])
  }
  df
}

# Drop list-columns from an sf object before writing. GeoPackage (and most OGR
# drivers) cannot serialize list-columns, so they must be removed: geometry
# preserved: prior to sf::st_write. Mirrors the inline idiom it replaces
# (`names(lc)[lc]`), so the result is byte-identical.
.drop_list_cols <- function(x) {
  lc <- vapply(sf::st_drop_geometry(x), is.list, logical(1L))
  for (nm in names(lc)[lc]) x[[nm]] <- NULL
  x
}

# Split a delimited character column into one row per token, recycling the
# other columns. Base-R replacement for tidyr::separate_longer_delim() /
# separate_rows(); `delim` is a fixed string. NA / empty cells yield a single
# row carrying that value (callers trim/filter as needed).
.separate_longer <- function(df, col, delim = ",") {
  parts <- strsplit(as.character(df[[col]]), delim, fixed = TRUE)
  parts[lengths(parts) == 0L] <- list(NA_character_)
  out <- df[rep.int(seq_len(nrow(df)), lengths(parts)), , drop = FALSE]
  out[[col]] <- unlist(parts, use.names = FALSE)
  rownames(out) <- NULL
  out
}

# Read a lake-membership mask: a parquet mapping a reference flowline id
# (`comid`) to the waterbody/lake id it belongs to (`wbareacomi`). This is the
# materialized, versioned replacement for a live nhdplusTools::get_vaa() call.
# The CONUS mask is built from NHDPlusV2 wbareacomi;
# oCONUS / superconus domains supply an equivalent parquet with the same two
# columns. Errors clearly when the path is required but absent.
.read_lake_mask <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("lake_mask_path is required for lake processing (Step 6h) but was ",
      "not supplied or does not exist: ",
      if (is.null(path)) "NULL" else path, call. = FALSE)
  }
  if (!requireNamespace("arrow", quietly = TRUE))
    stop("reading the lake mask requires the 'arrow' package: ",
      "install.packages('arrow')", call. = FALSE)
  m <- as.data.frame(arrow::read_parquet(path))
  if (!all(c("comid", "wbareacomi") %in% names(m))) {
    stop("lake mask parquet must have columns `comid` and `wbareacomi`; got: ",
      paste(names(m), collapse = ", "), call. = FALSE)
  }
  m
}

# .hf_fmt_id (scientific-notation-safe id stringifier) moved to the shared base:
# hfutils::hf_fmt_id (hfutils >= 0.3.5). Call sites use hfutils::hf_fmt_id().

# Standardise flowline column names by scanning candidate lists.
#
# Picks the first existing column from each candidate vector and renames it
# to the canonical name (flowpath_id / flowpath_toid / levelpath). Terminal
# toid values that are NA or empty-string are replaced with default_toid.
#
# Returns a list:
#   $flowline : normalised sf/data.frame
#   $id_col   : canonical id column name after renaming (always "flowpath_id")
#   $toid_col : canonical toid column name ("flowpath_toid") or NA
#   $level_col: canonical level column name ("levelpath") or NA
.hf_normalize_flowline <- function(flowline,
                                   id_candidates    = c("flowpath_id", "COMID", "id"),
                                   toid_candidates  = c("flowpath_toid", "toCOMID", "toid"),
                                   level_candidates = c("levelpath", "levelpathid", "mainstem_id"),
                                   default_toid     = "0",
                                   ensure_alias     = FALSE) {
  nms <- names(flowline)

  # -- id column ------------------------------------------------------
  id_found <- id_candidates[id_candidates %in% nms][1]
  if (is.na(id_found)) stop("No id column found; tried: ", paste(id_candidates, collapse = ", "))
  if (id_found != "flowpath_id") names(flowline)[names(flowline) == id_found] <- "flowpath_id"
  flowline$flowpath_id <- as.character(flowline$flowpath_id)

  # -- toid column -----------------------------------------------------
  toid_found <- toid_candidates[toid_candidates %in% names(flowline)][1]
  if (!is.na(toid_found)) {
    if (toid_found != "flowpath_toid") names(flowline)[names(flowline) == toid_found] <- "flowpath_toid"
    flowline$flowpath_toid <- as.character(flowline$flowpath_toid)
    flowline$flowpath_toid[is.na(flowline$flowpath_toid) | !nzchar(flowline$flowpath_toid)] <- default_toid
    toid_col <- "flowpath_toid"
  } else {
    flowline$flowpath_toid <- default_toid
    toid_col <- NA_character_
  }

  # -- level column ----------------------------------------------------
  level_found <- level_candidates[level_candidates %in% names(flowline)][1]
  if (!is.na(level_found)) {
    if (level_found != "levelpath") names(flowline)[names(flowline) == level_found] <- "levelpath"
    level_col <- "levelpath"
  } else {
    level_col <- NA_character_
  }

  list(flowline  = flowline,
    id_col    = "flowpath_id",
    toid_col  = toid_col,
    level_col = level_col)
}

# Standardise a lookup table data.frame.
#
# Returns NULL if lookup_table is NULL / empty or lacks required_id_col.
# Coerces required_id_col to character, deduplicates, and optionally
# fills original_flowpath_id from the id column by stripping the ".suffix"
# appended to split flowpath ids (e.g. "123.1" -> "123").
.hf_normalize_lookup_table <- function(lookup_table,
                                       required_id_col       = "refactored_flowpath_id",
                                       ensure_original_from_id = FALSE) {
  if (is.null(lookup_table) || nrow(lookup_table) == 0) return(NULL)
  if (!required_id_col %in% names(lookup_table)) return(NULL)

  lookup_table[[required_id_col]] <- as.character(lookup_table[[required_id_col]])
  lookup_table <- lookup_table[!is.na(lookup_table[[required_id_col]]) &
    nzchar(lookup_table[[required_id_col]]), , drop = FALSE]
  lookup_table <- lookup_table[!duplicated(lookup_table[[required_id_col]]), , drop = FALSE]

  if (isTRUE(ensure_original_from_id) && !"original_flowpath_id" %in% names(lookup_table)) {
    lookup_table$original_flowpath_id <-
      sub("\\..*$", "", as.character(lookup_table[[required_id_col]]))
  }

  if (nrow(lookup_table) == 0) return(NULL)
  lookup_table
}

# Read a non-spatial table from a GeoPackage (SQLite).
#
# Returns a data.frame (possibly with original_flowpath_id filled from
# required_id_col when fill_original_from_id = TRUE), or NULL when the
# table is absent or the file does not exist.
.hf_read_lookup_table <- function(gpkg,
                                  table_name          = "lookup_table",
                                  required_id_col     = "refactored_flowpath_id",
                                  fill_original_from_id = FALSE,
                                  warn_on_missing_db  = FALSE,
                                  warn_prefix         = "") {
  if (!file.exists(gpkg)) {
    if (warn_on_missing_db)
      .hf_log_warn(warn_prefix, "gpkg not found: ", gpkg)
    return(NULL)
  }

  con <- tryCatch(DBI::dbConnect(RSQLite::SQLite(), gpkg, flags = RSQLite::SQLITE_RO),
    error = function(e) NULL)
  if (is.null(con)) return(NULL)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  tables <- DBI::dbListTables(con)
  if (!table_name %in% tables) {
    if (warn_on_missing_db)
      .hf_log_warn(warn_prefix, "table '", table_name, "' not found in ", gpkg)
    return(NULL)
  }

  tbl <- tryCatch(DBI::dbReadTable(con, table_name), error = function(e) NULL)
  if (is.null(tbl) || nrow(tbl) == 0) return(NULL)

  if (!required_id_col %in% names(tbl)) return(NULL)
  tbl[[required_id_col]] <- as.character(tbl[[required_id_col]])

  if (isTRUE(fill_original_from_id) && !"original_flowpath_id" %in% names(tbl)) {
    tbl$original_flowpath_id <-
      sub("\\..*$", "", as.character(tbl[[required_id_col]]))
  }

  tbl
}

# hfutils::hf_assert_network_dag and hfutils::hf_recompute_hydroseq moved to the shared base:
# hfutils::hf_assert_network_dag / hf_recompute_hydroseq (hfutils >= 0.3.6).

# Apply DAG enforcement and hydroseq recompute to a flowpaths layer in a
# GeoPackage, writing the corrected layer back in place.  Call this after any
# operation that may have modified flowpath_id / flowpath_toid topology.
# Returns invisibly TRUE if cycles were found (and broken), FALSE otherwise.
.hf_enforce_dag_gpkg <- function(gpkg, layer = "flowpaths",
                                 id_col   = "flowpath_id",
                                 toid_col = "flowpath_toid",
                                 vpuid    = NULL) {
  if (!file.exists(gpkg)) return(invisible(FALSE))
  if (!hfutils::layer_exists(gpkg, layer)) return(invisible(FALSE))

  fp    <- sf::read_sf(gpkg, layer, quiet = TRUE)
  dag_v <- hfutils::hf_assert_network_dag(fp, id_col = id_col, toid_col = toid_col)
  tag   <- if (!is.null(vpuid)) sprintf("[%s][%s]", vpuid, layer) else sprintf("[%s]", layer)

  had_cycles <- FALSE
  if (!dag_v$is_dag) {
    .hf_log_warn(sprintf("%s BUG: cycles detected (%s) -- topology propagation is incorrect; breaking as fallback",
      tag, dag_v$message))
    fp         <- hfutils::hf_break_cycles(fp, id_col = id_col, toid_col = toid_col)
    had_cycles <- TRUE
  }

  fp <- hfutils::hf_recompute_hydroseq(fp, id_col = id_col, toid_col = toid_col)

  fp <- .drop_list_cols(fp)
  sf::st_write(fp, gpkg, layer, delete_layer = TRUE, quiet = TRUE)

  if (!had_cycles)
    message(sprintf("%s DAG OK; hydroseq recomputed (%d fps)", tag, nrow(fp)))
  invisible(had_cycles)
}

# `.hf_break_cycles` now lives in hfutils as `hf_break_cycles` (hfutils >= 0.3.6).

# GPKG safe-write helpers moved to the shared base (hfutils >= 0.3.5):
#    hfutils::gpkg_update_col / gpkg_update_geom / gpkg_exec
#    (with the internal trigger save/restore helpers). Call sites use hfutils::.

# Write the partition/reconcile output layers to a GeoPackage.
# Used by partition_reconcile_divides(_from_objects) when write_outputs = TRUE.
# Each non-empty spatial layer is written under the same name the gpkg-based
# entry point reads (flowline / divides / events / split_events / pois); the
# non-spatial lookup_table is written as a plain table. require_lookup_db stops
# if a lookup_table was expected but is empty.
.hf_write_outputs_bundle <- function(output_gpkg, flowline, events = NULL,
                                     split_events = NULL, divides = NULL,
                                     lookup_table = NULL, pois = NULL,
                                     require_lookup_db = FALSE) {
  .write_layer <- function(obj, layer) {
    if (is.null(obj) || !inherits(obj, "sf") || nrow(obj) == 0L)
      return(invisible(FALSE))
    sf::st_write(obj, output_gpkg, layer, delete_layer = TRUE, quiet = TRUE)
    invisible(TRUE)
  }
  .write_layer(flowline,     "flowline")
  .write_layer(divides,      "divides")
  .write_layer(events,       "events")
  .write_layer(split_events, "split_events")
  .write_layer(pois,         "pois")

  has_lookup <- !is.null(lookup_table) && nrow(lookup_table) > 0L
  if (has_lookup) {
    con <- DBI::dbConnect(RSQLite::SQLite(), output_gpkg)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    DBI::dbWriteTable(con, "lookup_table",
      as.data.frame(lookup_table), overwrite = TRUE)
  } else if (isTRUE(require_lookup_db)) {
    stop(".hf_write_outputs_bundle: lookup_table required but empty/NULL")
  }
  invisible(output_gpkg)
}

# Pipeline diagnostics go through logger, not R's condition system.
#
# R buffers a warning raised inside a pipeline stage and the driver prints only
# "There were 50 or more warnings" at the end, with no content, a finding that
# never reaches the log cannot be acted on.
#
# Takes `...` and concatenates, so a call passing several fragments keeps all of
# them. `skip_formatter` stops a literal brace being read as glue interpolation.
#
# Exported functions keep the condition system, where an interactive caller
# actually sees it.
#' @keywords internal
#' @noRd
.hf_log_warn <- function(...) {
  msg <- paste0(...)
  logger::log_warn(logger::skip_formatter(paste(msg, collapse = "")))
  invisible(NULL)
}
