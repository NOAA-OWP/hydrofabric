#' Finalize NextGen output to the Lynker hydrofabric data model spec
#'
#' @description
#' Applies the final column renames, adds VPU/streamorder/slope fields, rolls
#' up hydrolocation attributes on POIs, and promotes a dedicated
#' `hydrolocations` layer. Mutates the provided ngen GeoPackage in place.
#'
#' See `docs/data_model.md` for the full mapping against the spec.
#'
#' @param ngen_gpkg Path to the NextGen GeoPackage produced by `nextgen_from_gpkg` +
#'   cleanup steps.
#' @param reference_gpkg Path to the reference fabric (used to source `slope`
#'   per member reference COMID). If `NULL`, slope stays `NA`.
#' @param vpuid Single VPU identifier to stamp on every layer (character).
#' @param events_layer Source layer for hydrolocation attributes. Default
#'   `"pois"`, the pre-ngen pois layer carries the full hl_* set (hl_class,
#'   hl_type, hl_source, hl_uri, etc). The `events` layer in the intermediate
#'   gpkg only has a subset, so prefer `pois`.
#' @param events_gpkg Optional: if the hydrolocation source layer isn't in
#'   `ngen_gpkg`, read it from this gpkg (typically the pre-ngen intermediate).
#' @param slope_parquet Optional path to a parquet of per-reference-id slope
#'   used to populate `slope`: either a slope table (`id`, `slope_m_per_m`) or
#'   a channel table (`flowpath_id`, `slope_m_m`). `NULL` leaves `slope` `NA`
#'   (default).
#'
#' @return Invisibly, the input `ngen_gpkg` path.
#' @examples
#' \dontrun{
#' finalize_data_model(ngen_gpkg      = "ngen_vpu01.gpkg",
#'                     reference_gpkg = "reference_01.gpkg",
#'                     vpuid          = "01",
#'                     events_gpkg    = "out_vpu01.gpkg")
#' }
#' @export
finalize_data_model <- function(ngen_gpkg,
                                reference_gpkg = NULL,
                                vpuid,
                                events_layer   = "pois",
                                events_gpkg    = NULL,
                                slope_parquet  = NULL) {
  stopifnot(is.character(vpuid), length(vpuid) == 1L, nzchar(vpuid))

  # Length-weighted slope is computed by a single vectorised data.table
  # aggregation (fast), so it is ON by default, a populated, positive `slope`
  # is part of the NextGen data model and must appear consistently across
  # flowpaths, flowlines, and network.
  # Disable with HF_NGEN_SLOPE=false or options(hydrofabric.ngen_slope = FALSE).
  compute_slope <- isTRUE(as.logical(getOption("hydrofabric.ngen_slope",
    Sys.getenv("HF_NGEN_SLOPE", "true"))))
  if (is.na(compute_slope)) compute_slope <- FALSE

  # sub-step timer (diagnostic: per-step breakdown of where finalize spends time)
  .ft <- function(label, expr) {
    t0 <- proc.time()[["elapsed"]]
    force(expr)
    message(sprintf("  finalize[%-22s] %6.1fs", label, proc.time()[["elapsed"]] - t0))
  }

  # -- 1. Rename columns to match the spec ----------------------------------
  .ft("rename_cols", {
    .rename_layer_cols(ngen_gpkg, "flowpaths", list(
      mainstem              = "mainstem_id",
      tot_drainage_areasqkm = "total_dasqkm"))
    .rename_layer_cols(ngen_gpkg, "flowlines", list(
      levelpathid           = "mainstem_id",
      flowline_incr         = "incremental_areasqkm",
      tot_drainage_areasqkm = "total_dasqkm"))
    # network is aspatial; read/rename/write via DBI to avoid sf round-trips
    .rename_aspatial_cols(ngen_gpkg, "network", list(
      flowline_area = "flowline_total_dasqkm"))
    # The flowlines rename is a no-op when total_dasqkm already exists (a clean,
    # divide-based value is set upstream), leaving the fac-derived
    # tot_drainage_areasqkm orphaned alongside it. That column is unreliable at
    # flowline grain: FAC is sampled at fragment endpoints, which land
    # off-channel for ~16% of fragments (near-zero drainage), and is fully
    # superseded by the authoritative divide-based total_dasqkm. Drop it so only
    # the correct drainage area ships (guarded to never remove it unless
    # total_dasqkm is present).
    .drop_layer_col(ngen_gpkg, "flowlines", "tot_drainage_areasqkm",
      require = "total_dasqkm")
  })

  # -- 2. Stamp vpuid on every layer ----------------------------------------
  .ft("stamp_vpuid", {
    for (lyr in c("flowpaths", "divides", "flowlines", "pois", "nexus")) {
      if (hfutils::layer_exists(ngen_gpkg, lyr))
        .add_static_col(ngen_gpkg, lyr, "vpuid", vpuid)
    }
    .add_aspatial_col(ngen_gpkg, "network", "vpuid", vpuid)
  })

  # -- 3. Add Strahler streamorder (native, fp -> fp flattened via nexus) ----
  .ft("streamorder", .add_streamorder(ngen_gpkg))

  # -- 4. Length-weighted slope from reference per member reference COMID ---------
  .ft("slope", {
    if (compute_slope && !is.null(reference_gpkg) && file.exists(reference_gpkg)) {
      .add_slope(ngen_gpkg, reference_gpkg, vpuid, slope_parquet = slope_parquet)
    } else {
      if (compute_slope)
        message("finalize_data_model: slope requested but reference_gpkg missing -- writing NA")
      else
        message("finalize_data_model: slope skipped (set HF_NGEN_SLOPE=true to enable) -- writing NA")
      # Ensure column exists as NA so downstream consumers see the field
      .add_static_col(ngen_gpkg, "flowpaths", "slope", NA_real_)
    }
  })

  # -- 5. Propagate flowpath attrs (streamorder, slope, flowpath_hydroseq)
  .ft("propagate_fp_to_fl", .propagate_fp_to_fl(ngen_gpkg))

  # -- 6. Promote hydrolocations into a dedicated layer ---------------------
  .ft("promote_hydrolocations", {
    ev_source <- if (!is.null(events_gpkg)) events_gpkg else ngen_gpkg
    if (hfutils::layer_exists(ev_source, events_layer)) {
      .promote_hydrolocations(ngen_gpkg, ev_source, events_layer, vpuid)
    }
  })

  # -- 6b. Withdraw POIs the network cannot place at an outlet --------------
  # Runs after the hydrolocations layer exists, so the hydrolocation and its
  # reason survive even though the POI does not.
  .ft("withdraw_unplaceable_pois", .withdraw_unplaceable_pois(ngen_gpkg))

  # -- 7. POI rollup fields (hl_classes / hl_types / hl_count) --------------
  .ft("poi_rollups", {
    if (hfutils::layer_exists(ngen_gpkg, "hydrolocations")) {
      .add_poi_rollups(ngen_gpkg)
    }
  })

  # -- 8. Denormalize network with flowpath-scale + hf_* + hl_reference -----
  .ft("finalize_network", .finalize_network(ngen_gpkg, reference_gpkg, vpuid))

  # -- 9. Populate incremental area + total drainage area ------------------
  .ft("drainage_area", .finalize_drainage_area(ngen_gpkg))

  # -- 10. Flag sub-threshold divides with the reason they are small --------
  .ft("flag_small_divides", .flag_small_divides(ngen_gpkg, min_area_sqkm = 3))

  # -- 11. Gage-POI drainage-area consistency gate --------------------------
  # Publish, per gage POI, the drainage area it is delineated to alongside the
  # reference drainage area of the reaches it is built from. Reported, not
  # judged: the two numbers and their ratio are written to the POI layer so a
  # consumer can decide what divergence matters for their use.
  .ft("report_poi_drainage_area", .report_poi_drainage_area(ngen_gpkg))

  # -- 11b. Nested-set upstream index (O(1) upstream range queries) ---------
  # Placed after all layer-rewriting steps so the index is not clobbered, and
  # before enforce_integer_ids so upstream_id / num_upstreams are typed INTEGER.
  .ft("upstream_index", add_upstream_index(ngen_gpkg))

  # -- 12. Type integer-valued id columns as INTEGER + index join keys ------
  # Value-guarded: prefixed string ids (fp-/cat-/nex-) stay TEXT; only genuine
  # integer columns (hydroseq, mainstem_id, streamorder, ...) are retyped. Keeps
  # ids off floating point and indexes the toid joins consumers rely on.
  .ft("enforce_integer_ids", enforce_integer_ids(ngen_gpkg))

  # -- 13. National mainstem ids --------------------------------------------
  # Aggregation keys level paths by a number local to the VPU; the reference's
  # own level-path id travels alongside as orig_levelpathID and becomes
  # mainstem_id here, so mainstem_id means the same river in every VPU. It runs
  # after the columns are typed 64-bit INTEGER and changes them in SQL only:
  # reference ids exceed 32 bits, and a read through GDAL of a 32-bit column
  # would wrap them.
  .ft("national_mainstem", .national_mainstem_ids(ngen_gpkg))

  invisible(ngen_gpkg)
}

# Annotate every sub-`min_area_sqkm` divide with WHY it is small, in a
# `small_reason` column (NA for divides at/above the threshold). After the
# ideal-size cross-levelpath collapse + mainstem upstream-absorb, the only
# legitimately-small divides are gage POIs (`poi_outlet`), domain terminals
# (`terminal_outlet`), and real junctions (`confluence`). The `headwater` and
# `mid_reach` buckets SHOULD be empty, any non-zero count is emitted as a
# warning, i.e. a built-in tripwire that the aggregation left a reducible unit.
# Purely additive: rewrites the divides layer with one extra column, no geometry
# or topology change. Classification reuses the final, reconciled data model
# (poi_id, terminal toid, and inflow count over the fp->nexus->fp graph).
.flag_small_divides <- function(gpkg, min_area_sqkm = 3) {
  if (!hfutils::layer_exists(gpkg, "divides") ||
    !hfutils::layer_exists(gpkg, "flowpaths")) return(invisible(NULL))
  dv <- sf::st_read(gpkg, "divides", quiet = TRUE)
  if (!all(c("divide_id", "areasqkm") %in% names(dv))) return(invisible(NULL))
  fp  <- sf::st_drop_geometry(sf::st_read(gpkg, "flowpaths", quiet = TRUE))
  nex <- if (hfutils::layer_exists(gpkg, "nexus"))
    sf::st_drop_geometry(sf::st_read(gpkg, "nexus", quiet = TRUE)) else NULL
  fp$flowpath_id <- as.character(fp$flowpath_id)

  # inflow count per flowpath over the resolved fp -> nexus -> fp topology
  ds <- if (!is.null(nex)) {
    nm <- stats::setNames(as.character(nex$nexus_toid), as.character(nex$nexus_id))
    unname(nm[as.character(fp$flowpath_toid)])
  } else rep(NA_character_, nrow(fp))
  inflow <- table(ds[!is.na(ds)])
  nin <- as.integer(inflow[fp$flowpath_id])
  nin[is.na(nin)] <- 0L
  is_term <- is.na(fp$flowpath_toid) | as.character(fp$flowpath_toid) == "0" |
    grepl("^tnx-", as.character(fp$flowpath_toid))
  haspoi  <- "poi_id" %in% names(fp) & !is.na(fp$poi_id)
  # does this flowpath drain INTO a connector (divide-less waterbody flowpath)?
  # such a small divide is a lake-edge tributary: genuinely small and unable to
  # merge across the waterbody, NOT a reducible aggregation miss.
  hd_map  <- if ("has_divide" %in% names(fp))
    stats::setNames(as.logical(fp$has_divide), fp$flowpath_id) else NULL
  ds_is_connector <- if (!is.null(hd_map)) {
    v <- hd_map[ds]
    !is.na(v) & !v
  } else rep(FALSE, nrow(fp))

  # precedence: a fixed outlet (POI/terminal) is the reason even at a junction;
  # waterbody-isolation outranks the headwater/mid_reach tripwire categories
  reason <- rep(NA_character_, nrow(fp))
  reason[haspoi]                    <- "poi_outlet"
  reason[is.na(reason) & is_term]   <- "terminal_outlet"
  reason[is.na(reason) & nin >= 2L] <- "confluence"
  reason[is.na(reason) & ds_is_connector] <- "waterbody"  # drains to a connector
  reason[is.na(reason) & nin == 1L] <- "mid_reach"   # tripwire: expected empty
  reason[is.na(reason) & nin == 0L] <- "headwater"   # tripwire: expected empty
  rmap <- stats::setNames(reason, fp$flowpath_id)

  dv$flowpath_id <- as.character(dv$flowpath_id)
  small <- !is.na(dv$areasqkm) & dv$areasqkm < min_area_sqkm
  dv$small_reason <- NA_character_
  dv$small_reason[small] <- unname(rmap[dv$flowpath_id[small]])

  bad <- sum(dv$small_reason %in% c("headwater", "mid_reach"), na.rm = TRUE)
  if (bad > 0L)
    .hf_log_warn(sprintf(paste0(
      ".flag_small_divides: %d small divide(s) are headwater/mid_reach -- ",
      "aggregation left a reducible unit"), bad))
  tb <- table(dv$small_reason)
  message(sprintf("  small_reason: %d/%d divides < %g km2 (%s)",
    sum(small), nrow(dv), min_area_sqkm,
    if (length(tb)) paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = ", ") else "none"))

  dv <- .drop_list_cols(dv)
  sf::st_write(dv, gpkg, "divides", delete_layer = TRUE, quiet = TRUE)
  invisible(NULL)
}

# Populate incremental catchment area (`areasqkm`) and total upstream drainage
# area (`total_dasqkm`) on flowpaths + network from the final divides, keyed by
# the divide->flowpath link (`flowpath_id`). This is the authoritative
# computation: the aggregate-stage `add_measures` can silently zero these when
# the divides at that point lack the `flowpath_id` link (divide_id != flowpath_id
# in the cat-*/fp-* schema), so we (re)compute on the final, reconciled data
# model where the link is guaranteed present (step 9). Total DA accumulates over
# the flowpath -> nexus -> flowpath resolved topology (ngen flowpath_toid points
# at a nexus). Idempotent; trigger-safe; rolls back cleanly on flaky storage.
.finalize_drainage_area <- function(gpkg) {
  con <- DBI::dbConnect(RSQLite::SQLite(), gpkg)
  on.exit(if (DBI::dbIsValid(con)) try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)
  DBI::dbExecute(con, "PRAGMA synchronous = NORMAL")
  DBI::dbExecute(con, "PRAGMA journal_mode = DELETE")
  tabs <- DBI::dbListTables(con)
  if (!all(c("divides", "flowpaths", "nexus") %in% tabs)) return(invisible(NULL))

  div <- DBI::dbGetQuery(con, "SELECT flowpath_id, areasqkm FROM divides WHERE flowpath_id IS NOT NULL")
  fp  <- DBI::dbGetQuery(con, "SELECT flowpath_id, flowpath_toid FROM flowpaths")
  nex <- DBI::dbGetQuery(con, "SELECT nexus_id, nexus_toid FROM nexus")
  fp$flowpath_id <- as.character(fp$flowpath_id)
  nexmap <- stats::setNames(as.character(nex$nexus_toid), as.character(nex$nexus_id))
  fp$ds_fp <- unname(nexmap[as.character(fp$flowpath_toid)])           # resolve fp -> nex -> fp
  area <- stats::setNames(as.numeric(div$areasqkm), as.character(div$flowpath_id))
  fp$areasqkm <- unname(ifelse(is.na(area[fp$flowpath_id]), 0, area[fp$flowpath_id]))
  fp$total_dasqkm <- hfutils::accumulate_downstream(fp, id = "flowpath_id", toid = "ds_fp", attr = "areasqkm")

  DBI::dbWriteTable(con, "bf_tmp", fp[, c("flowpath_id", "areasqkm", "total_dasqkm")],
    temporary = TRUE, overwrite = TRUE)
  DBI::dbExecute(con, "CREATE INDEX bf_idx ON bf_tmp(flowpath_id)")
  DBI::dbExecute(con, "BEGIN")
  tryCatch(
    {
      for (tab in intersect(c("flowpaths", "network"), tabs)) {
        cols <- DBI::dbListFields(con, tab)
        if (!all(c("flowpath_id", "areasqkm", "total_dasqkm") %in% cols)) next  # defensive
        trg <- DBI::dbGetQuery(con, sprintf(
          "SELECT name, sql FROM sqlite_master WHERE type='trigger' AND tbl_name='%s'", tab))
        for (nm in trg$name) DBI::dbExecute(con, sprintf('DROP TRIGGER IF EXISTS "%s"', nm))
        DBI::dbExecute(con, sprintf(
          "UPDATE %s SET areasqkm=(SELECT b.areasqkm FROM bf_tmp b WHERE b.flowpath_id=%s.flowpath_id),
                      total_dasqkm=(SELECT b.total_dasqkm FROM bf_tmp b WHERE b.flowpath_id=%s.flowpath_id)
         WHERE flowpath_id IN (SELECT flowpath_id FROM bf_tmp)", tab, tab, tab))
        for (s in trg$sql[!is.na(trg$sql)]) DBI::dbExecute(con, s)         # recreate triggers
      }
    },
    error = function(e) {
      try(DBI::dbExecute(con, "ROLLBACK"), silent = TRUE)
      stop(e)
    })
  DBI::dbExecute(con, "COMMIT")
  DBI::dbExecute(con, "DROP TABLE bf_tmp")
  invisible(NULL)
}

# --------------------------- internals -----------------------------------

# Join streamorder, slope, and hydroseq (as flowpath_hydroseq) from flowpaths
# onto flowlines. All three are flowpath-scale attributes; sub-flowline segments
# inherit from their parent flowpath via the flowpath_id foreign key.
.propagate_fp_to_fl <- function(gpkg) {
  if (!hfutils::layer_exists(gpkg, "flowlines") ||
    !hfutils::layer_exists(gpkg, "flowpaths")) return(invisible(NULL))

  fl <- sf::st_read(gpkg, "flowlines", quiet = TRUE)
  already <- all(c("streamorder", "slope", "flowpath_hydroseq") %in% names(fl))
  if (already) return(invisible(NULL))

  fp <- sf::st_drop_geometry(sf::st_read(gpkg, "flowpaths", quiet = TRUE))
  lu_cols <- intersect(c("flowpath_id", "streamorder", "slope", "hydroseq"), names(fp))
  fp_lu <- fp[, lu_cols, drop = FALSE]
  fp_lu$flowpath_id <- as.character(fp_lu$flowpath_id)
  fl$flowpath_id    <- as.character(fl$flowpath_id)

  ix <- match(fl$flowpath_id, fp_lu$flowpath_id)
  if (!("streamorder"       %in% names(fl)) && "streamorder" %in% names(fp_lu))
    fl$streamorder       <- as.integer(fp_lu$streamorder[ix])
  if (!("slope"             %in% names(fl)) && "slope" %in% names(fp_lu)) {
    fl$slope             <- as.numeric(fp_lu$slope[ix])
    # Inherit the flowpath floor; guard the unmatched (NA) case so no flowline
    # carries a non-positive/null slope into routing. Skip if slope was never
    # computed (flowpath slopes all NA -> leave NA, "slope not computed").
    if (any(!is.na(fl$slope))) {
      bad <- is.na(fl$slope) | fl$slope < 1e-5
      if (any(bad)) fl$slope[bad] <- 1e-5
    }
  }
  if (!("flowpath_hydroseq" %in% names(fl)) && "hydroseq" %in% names(fp_lu))
    fl$flowpath_hydroseq <- as.integer(fp_lu$hydroseq[ix])

  fl <- .drop_list_cols(fl)
  sf::st_write(fl, gpkg, "flowlines", delete_layer = TRUE, quiet = TRUE)
}

.rename_layer_cols <- function(gpkg, layer, rename_map) {
  if (!hfutils::layer_exists(gpkg, layer)) return(invisible(NULL))
  x <- sf::st_read(gpkg, layer, quiet = TRUE)
  renamed_any <- FALSE
  for (old in names(rename_map)) {
    nw <- rename_map[[old]]
    if (old %in% names(x) && !(nw %in% names(x))) {
      names(x)[names(x) == old] <- nw
      renamed_any <- TRUE
    }
  }
  if (renamed_any) {
    x <- .drop_list_cols(x)
    sf::st_write(x, gpkg, layer, delete_layer = TRUE, quiet = TRUE)
  }
}

.drop_layer_col <- function(gpkg, layer, col, require = NULL) {
  if (!hfutils::layer_exists(gpkg, layer)) return(invisible(NULL))
  x <- sf::st_read(gpkg, layer, quiet = TRUE)
  if (!(col %in% names(x))) return(invisible(NULL))
  # Never strip the column unless the authoritative replacement is present, so
  # the layer is never left without a drainage-area attribute.
  if (!is.null(require) && !(require %in% names(x))) return(invisible(NULL))
  x[[col]] <- NULL
  x <- .drop_list_cols(x)
  sf::st_write(x, gpkg, layer, delete_layer = TRUE, quiet = TRUE)
}

.rename_aspatial_cols <- function(gpkg, table, rename_map) {
  con <- DBI::dbConnect(RSQLite::SQLite(), gpkg)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  if (!DBI::dbExistsTable(con, table)) return(invisible(NULL))
  cols <- DBI::dbListFields(con, table)
  need <- intersect(names(rename_map), cols)
  if (length(need) == 0L) return(invisible(NULL))
  df <- DBI::dbReadTable(con, table)
  for (old in need) {
    nw <- rename_map[[old]]
    names(df)[names(df) == old] <- nw
  }
  DBI::dbWriteTable(con, table, df, overwrite = TRUE)
}

.add_static_col <- function(gpkg, layer, col, value) {
  if (!hfutils::layer_exists(gpkg, layer)) return(invisible(NULL))
  x <- sf::st_read(gpkg, layer, quiet = TRUE)
  if (col %in% names(x)) return(invisible(NULL))
  x[[col]] <- value
  x <- .drop_list_cols(x)
  sf::st_write(x, gpkg, layer, delete_layer = TRUE, quiet = TRUE)
}

.add_aspatial_col <- function(gpkg, table, col, value) {
  con <- DBI::dbConnect(RSQLite::SQLite(), gpkg)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  if (!DBI::dbExistsTable(con, table)) return(invisible(NULL))
  cols <- DBI::dbListFields(con, table)
  if (col %in% cols) return(invisible(NULL))
  DBI::dbExecute(con, sprintf('ALTER TABLE "%s" ADD COLUMN "%s" TEXT', table, col))
  DBI::dbExecute(con,
    sprintf('UPDATE "%s" SET "%s" = ?', table, col),
    params = list(value))
}

# Flatten fp -> nexus -> fp into direct fp -> fp edges, then native Strahler
# stream order via hfutils::get_streamorder (no hydroloom dependency).
.add_streamorder <- function(gpkg) {
  fp  <- sf::st_read(gpkg, "flowpaths", quiet = TRUE)
  nex <- if (hfutils::layer_exists(gpkg, "nexus"))
    sf::st_drop_geometry(sf::st_read(gpkg, "nexus", quiet = TRUE)) else NULL
  if ("streamorder" %in% names(fp)) return(invisible(NULL))

  nex_lu <- if (!is.null(nex))
    stats::setNames(as.character(nex$nexus_toid), as.character(nex$nexus_id)) else
    character(0)
  next_fp <- nex_lu[as.character(fp$flowpath_toid)]
  real_ids <- as.character(fp$flowpath_id)
  # Terminal nexuses in our schema point to synthetic "fp-1000..." ids that
  # aren't real flowpaths. Map those to the hydroloom outlet sentinel "".
  next_fp <- ifelse(is.na(next_fp) | !(next_fp %in% real_ids), "",
    as.character(next_fp))

  # Native Strahler order via hfutils::get_streamorder (igraph topo-sort, same
  # approach as get_hydroseq). Replaces hydroloom::add_streamorder, which
  # transiently flagged this fp->nexus->fp edge set "non-dendritic" and silently
  # dropped streamorder; the native path has no such failure mode and matches
  # hydroloom output exactly where hydroloom succeeds.
  next_fp[next_fp == ""] <- "0"
  edges <- data.frame(flowpath_id = real_ids, flowpath_toid = next_fp,
    stringsAsFactors = FALSE)
  so <- tryCatch(
    hfutils::get_streamorder(edges),
    error = function(e) {
      message("finalize_data_model: get_streamorder failed: ", conditionMessage(e))
      NULL
    })
  if (is.null(so)) return(invisible(NULL))
  fp$streamorder <- as.integer(so)   # aligned to fp row order

  fp <- .drop_list_cols(fp)
  sf::st_write(fp, gpkg, "flowpaths", delete_layer = TRUE, quiet = TRUE)
}

# Length-weighted mean slope per flowpath, pulled from reference by reference COMID.
# Sentinel 1e-5 (NHD "unknown" slope) is excluded IF the flowpath has at least
# one non-sentinel member: otherwise the sentinel is kept so the field is
# never silently dropped.
# slope_parquet: optional path to a parquet keyed by reference id, in either of
# two shapes: a slope table (`id`, `slope_m_per_m`) or a channel table
# (`flowpath_id` or `feature_id`, `slope_m_m`). When provided and the reference
# lacks slope, the parquet is used as the slope source; lengthkm still comes
# from the reference.
.slope_table_cols <- function(nms) {
  if (all(c("id", "slope_m_per_m") %in% nms))
    return(c(id = "id", slope = "slope_m_per_m"))
  key <- intersect(c("feature_id", "flowpath_id"), nms)
  sl  <- intersect(.HF_DINGMAN_COLS$slope, nms)
  if (length(key) && length(sl)) c(id = key[[1]], slope = sl[[1]]) else NULL
}

.add_slope <- function(gpkg, reference_gpkg, vpuid, slope_parquet = NULL) {
  fp <- sf::st_read(gpkg, "flowpaths", quiet = TRUE)
  if ("slope" %in% names(fp) && any(!is.na(fp$slope))) return(invisible(NULL))

  con <- DBI::dbConnect(RSQLite::SQLite(), reference_gpkg)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  ref_cols <- DBI::dbListFields(con, "flowpaths")

  # Pull lengthkm (always from reference) + slope (from parquet or reference)
  ref_lkm <- DBI::dbGetQuery(con, sprintf(
    "SELECT flowpath_id, lengthkm FROM flowpaths WHERE vpuid = '%s'", vpuid))
  ref_lkm$flowpath_id <- as.character(ref_lkm$flowpath_id)
  length_lu <- stats::setNames(as.numeric(ref_lkm$lengthkm), ref_lkm$flowpath_id)

  use_parquet <- !is.null(slope_parquet) && file.exists(slope_parquet)
  if (use_parquet) {
    if (!requireNamespace("arrow", quietly = TRUE))
      stop("reading `slope_parquet` requires the 'arrow' package: ",
        "install.packages('arrow')", call. = FALSE)
    slope_cols <- .slope_table_cols(
      names(arrow::read_parquet(slope_parquet, as_data_frame = FALSE)$schema))
    if (is.null(slope_cols)) {
      message("  .add_slope: ", basename(slope_parquet),
        " carries no slope column; not used as the slope source")
      use_parquet <- FALSE
    }
  }
  ref_has_slope <- "slope" %in% ref_cols

  if (!use_parquet && !ref_has_slope) {
    fp$slope <- NA_real_
    fp <- .drop_list_cols(fp)
    sf::st_write(fp, gpkg, "flowpaths", delete_layer = TRUE, quiet = TRUE)
    return(invisible(NULL))
  }

  if (use_parquet) {
    sp <- arrow::read_parquet(slope_parquet, col_select = unname(slope_cols))
    ids <- sp[[slope_cols[["id"]]]]
    ids <- if (is.character(ids)) ids else hfutils::hf_fmt_id(ids)
    slope_lu <- stats::setNames(as.numeric(sp[[slope_cols[["slope"]]]]), ids)
  } else {
    # Pull slope + lengthkm for all ref flowpaths in this VPU
    ref <- DBI::dbGetQuery(con, sprintf(
      "SELECT flowpath_id, slope FROM flowpaths WHERE vpuid = '%s'", vpuid))
    ref$flowpath_id <- as.character(ref$flowpath_id)
    slope_lu <- stats::setNames(as.numeric(ref$slope), ref$flowpath_id)
  }

  # NHD slope sentinels: 1e-5 ("unknown flat") and negative values in the
  # -9999 family ("unknown slope"). Always drop invalid negatives (no such
  # thing as a negative reach slope in the hydrofabric). Drop the 1e-5
  # sentinel only when at least one non-sentinel positive survives (per
  # user spec: "excluding 1e-5 if more than one value exists"). If nothing
  # usable survives, return NA rather than a bogus sentinel-averaged value.
  is_flat_sentinel <- function(s) abs(s - 1e-5) < 1e-9

  # Vectorised length-weighted slope. Flatten all members once, strip the NHD
  # split suffix in a single regex pass, look up slope/length once, then a
  # grouped weighted mean in one O(total_members) data.table aggregation.
  # Sentinel rule preserved: drop the 1e-5 "flat" sentinel within a group only
  # when a non-sentinel positive survives.
  members <- strsplit(as.character(fp$member_comid), ",", fixed = TRUE)
  grp     <- rep.int(seq_along(members), lengths(members))
  m_flat  <- sub("[.][0-9]+$", "", trimws(unlist(members, use.names = FALSE)))
  sl      <- unname(slope_lu[m_flat])
  ln      <- unname(length_lu[m_flat])
  ok      <- is.finite(sl) & is.finite(ln) & ln > 0 & sl >= 0

  fp$slope <- NA_real_
  if (any(ok)) {
    dt <- data.table::data.table(grp = grp[ok], sl = sl[ok], ln = ln[ok])
    dt[, sent := is_flat_sentinel(sl)]
    dt[, keep := if (.N > 1L && any(!sent)) !sent else rep(TRUE, .N), by = grp]
    res <- dt[keep == TRUE, .(slope = sum(sl * ln) / sum(ln)), by = grp]
    fp$slope[res$grp] <- res$slope
  }

  # Guarantee a routable slope: floor at min_slope and fill any residual NA
  # (flowpaths whose members carried no usable DEM slope). Channel routing
  # (Manning / t-route) requires a strictly positive slope, no negative, zero,
  # or null may reach the model. The source parquet already floors at 1e-5
  # ("flat/unknown"); we apply the same floor here so aggregation, an alternate
  # slope source, or missing members can never emit a non-positive/NA slope.
  min_slope <- 1e-5
  bad <- is.na(fp$slope) | fp$slope < min_slope
  if (any(bad)) {
    message(sprintf("  .add_slope: floored %d non-positive/NA slope(s) to %g",
      sum(bad), min_slope))
    fp$slope[bad] <- min_slope
  }

  fp <- .drop_list_cols(fp)
  sf::st_write(fp, gpkg, "flowpaths", delete_layer = TRUE, quiet = TRUE)
}

# Read events, pull hl_* fields, write a dedicated `hydrolocations` layer.
.promote_hydrolocations <- function(ngen_gpkg, ev_gpkg, ev_layer, vpuid) {
  ev <- sf::st_read(ev_gpkg, ev_layer, quiet = TRUE)
  # Some pipelines keep event_kind; keep only POI-anchored rows.
  if ("event_kind" %in% names(ev)) {
    ev <- ev[is.na(ev$event_kind) | ev$event_kind == "poi", , drop = FALSE]
  }
  # hl_flowpath_id_source / hl_reach_change_reason carry the reach-assignment
  # lineage: the id the source prescribed, and the rule that re-pointed it. Without
  # them a consumer joining on the supplied flowpath_id sees a silent mismatch.
  wanted <- c("poi_id", "flowpath_id", "hl_type", "hl_link", "hl_class",
    "hl_source", "hl_reference", "hl_uri", "hl_id",
    "hl_flowpath_id_source", "hl_reach_change_reason")
  keep   <- intersect(wanted, names(ev))
  if (length(keep) == 0L) return(invisible(NULL))
  # Zero-POI fabric: no poi events / no ngen pois layer -> nothing to promote.
  if (nrow(ev) == 0L || !hfutils::layer_exists(ngen_gpkg, "pois")) return(invisible(NULL))

  # Re-resolve flowpath_id from the ngen flowpath layer if events carry the
  # pre-ngen integer id, the pois ngen layer already has the fp-prefixed id.
  ngen_pois <- sf::st_read(ngen_gpkg, "pois", quiet = TRUE)
  ngen_poi_to_fp <- stats::setNames(as.character(ngen_pois$flowpath_id),
    as.character(ngen_pois$poi_id))
  ev$flowpath_id <- ngen_poi_to_fp[as.character(ev$poi_id)]

  # Bring geometry from pre-ngen events (they have snap geom)
  geom_col <- attr(ev, "sf_column")
  if (sf::st_crs(ev) != sf::st_crs(ngen_pois))
    ev <- sf::st_transform(ev, sf::st_crs(ngen_pois))

  out_cols <- unique(c("poi_id", "flowpath_id", keep))
  out <- sf::st_sf(ev[, out_cols, drop = FALSE])
  out$vpuid <- vpuid
  # Drop rows with no hl_reference (not a real hydrolocation)
  if ("hl_reference" %in% names(out)) {
    out <- out[!is.na(out$hl_reference) & nzchar(out$hl_reference), , drop = FALSE]
  }
  if (nrow(out) == 0L) return(invisible(NULL))

  out <- .drop_list_cols(out)
  sf::st_write(out, ngen_gpkg, "hydrolocations",
    delete_layer = TRUE, quiet = TRUE)
}

# Withdraw POIs that the reconciled network cannot place at a reach outlet.
#
# A hydrolocation at the head of a headwater reach has no upstream reach to be the
# outlet of, so no placement satisfies the outlet convention: see "Locations that
# cannot be placed" in vignette("poi-contract").
#
# Nothing is lost: the hydrolocation keeps its identity, position and reach in the
# `hydrolocations` layer and carries `poi_drop_reason`, and the poi_id is cleared
# from `flowpaths` and `network` so no layer references a POI that is not there.
.withdraw_unplaceable_pois <- function(gpkg, tol_m = 30, min_split_m = 250) {
  if (!hfutils::layer_exists(gpkg, "pois") ||
      !hfutils::layer_exists(gpkg, "flowlines")) return(invisible(NULL))
  pois <- sf::st_read(gpkg, "pois", quiet = TRUE)
  fl   <- sf::st_read(gpkg, "flowlines", quiet = TRUE)
  if (!nrow(pois) || !nrow(fl) || !"poi_id" %in% names(pois)) return(invisible(NULL))
  if (!all(c("flowline_id", "flowline_toid") %in% names(fl))) return(invisible(NULL))

  if (sf::st_crs(pois) != sf::st_crs(fl)) pois <- sf::st_transform(pois, sf::st_crs(fl))
  ends <- function(first) sf::st_sfc(lapply(sf::st_geometry(fl), function(g) {
    m <- sf::st_coordinates(g)
    sf::st_point(unname(m[if (first) 1L else nrow(m), 1:2]))
  }), crs = sf::st_crs(fl))
  outs <- ends(FALSE); ins <- ends(TRUE)

  pg    <- sf::st_geometry(pois)
  d_out <- as.numeric(sf::st_distance(pg, outs[sf::st_nearest_feature(pg, outs)],
    by_element = TRUE))
  ni    <- sf::st_nearest_feature(pg, ins)
  d_in  <- as.numeric(sf::st_distance(pg, ins[ni], by_element = TRUE))

  # Nothing drains into the host reach -> no upstream outlet this POI could occupy.
  n_up  <- table(as.character(fl$flowline_toid))
  host  <- as.character(fl$flowline_id)[ni]
  ups   <- unname(n_up[host]); ups[is.na(ups)] <- 0L

  # Unplaceable means BOTH routes to an outlet are closed.
  #
  # The placement rule reassigns a POI within snap_dist_m of a reach's inlet to the
  # outlet of an incoming reach, and that tolerance is the same 250 m as the split
  # floor. On a headwater there is no incoming reach, so that route is closed; and a
  # POI closer to the inlet than min_split_m cannot be split into an outlet either,
  # because the upstream piece would fall under the floor. A POI further than
  # min_split_m from the inlet of a headwater IS splittable, so it is deliberately
  # not withdrawn, that would be a split failure to investigate, not a limitation
  # of the network.
  headwater <- ups == 0L
  at_inlet  <- d_out > tol_m & d_in <= 1 & headwater
  below_fl  <- d_out > tol_m & d_in > 1 & d_in <= min_split_m & headwater
  at_inlet[is.na(at_inlet)] <- FALSE
  below_fl[is.na(below_fl)] <- FALSE
  drop <- at_inlet | below_fl
  if (!any(drop)) {
    message("  finalize: no unplaceable POIs to withdraw")
    return(invisible(NULL))
  }
  dropped_ids <- unique(as.character(pois$poi_id)[drop])
  reason_of <- stats::setNames(
    ifelse(at_inlet[drop], "headwater_inlet_no_upstream_reach",
           "headwater_inlet_below_split_floor"),
    as.character(pois$poi_id)[drop])
  reason <- paste(sort(unique(unname(reason_of))), collapse = " / ")

  types <- if ("hl_types" %in% names(pois)) {
    tt <- sort(table(as.character(pois$hl_types)[drop]), decreasing = TRUE)
    sprintf(" [%s]", paste(sprintf("%s=%d", names(tt), as.integer(tt)), collapse = ", "))
  } else ""
  message(sprintf(paste0(
    "  finalize: withdrew %d POI(s) the network cannot place at an outlet%s; ",
    "hydrolocation(s) kept with poi_drop_reason = '%s'"),
    length(dropped_ids), types, reason))

  # Record the reason on the hydrolocation lineage before the POI goes.
  if (hfutils::layer_exists(gpkg, "hydrolocations")) {
    hl <- sf::st_read(gpkg, "hydrolocations", quiet = TRUE)
    if ("poi_id" %in% names(hl)) {
      hl$poi_drop_reason <- unname(reason_of[as.character(hl$poi_id)])
      sf::st_write(.drop_list_cols(hl), gpkg, "hydrolocations",
        delete_layer = TRUE, quiet = TRUE)
    }
  }

  sf::st_write(.drop_list_cols(pois[!drop, , drop = FALSE]), gpkg, "pois",
    delete_layer = TRUE, quiet = TRUE)

  # Clear the now-dangling references so no layer names a POI that is not shipped.
  #
  # The rtree triggers on a spatial table call Spatialite functions (ST_IsEmpty, ...)
  # that plain RSQLite does not have, and SQLite compiles a table's triggers when
  # preparing any UPDATE on it, so they are dropped for the statement and restored
  # verbatim from sqlite_master afterwards, as .retype_cols_integer does.
  con <- DBI::dbConnect(RSQLite::SQLite(), gpkg)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  qids <- paste(sprintf("'%s'", dropped_ids), collapse = ",")
  for (tb in c("flowpaths", "network")) {
    if (!DBI::dbExistsTable(con, tb)) next
    if (!"poi_id" %in% DBI::dbGetQuery(con,
      sprintf('PRAGMA table_info("%s")', tb))$name) next
    trg <- DBI::dbGetQuery(con, sprintf(paste0(
      "SELECT name, sql FROM sqlite_master WHERE type='trigger' AND tbl_name='%s' ",
      "AND sql IS NOT NULL"), tb))
    DBI::dbBegin(con)
    n <- tryCatch({
      for (nm in trg$name) DBI::dbExecute(con, sprintf('DROP TRIGGER "%s"', nm))
      k <- DBI::dbExecute(con, sprintf(
        'UPDATE "%s" SET poi_id = NULL WHERE CAST(poi_id AS TEXT) IN (%s)', tb, qids))
      for (s in trg$sql) DBI::dbExecute(con, s)
      DBI::dbCommit(con)
      k
    }, error = function(e) {
      DBI::dbRollback(con)
      .hf_log_warn(sprintf("withdraw_unplaceable_pois: could not clear poi_id on %s: %s",
        tb, conditionMessage(e)))
      0L
    })
    if (n > 0L) message(sprintf("  finalize: cleared poi_id on %d %s row(s)", n, tb))
  }
  invisible(dropped_ids)
}

# Denormalize the network table per the spec: add flowpath-scale attributes
# (mainstem_id, hydroseq, areasqkm, lengthkm, total_dasqkm), NHD-source
# crosswalks (hf_source, hf_toid, hf_hydroseq), flowline-scale attrs
# (flowline_hydroseq, flowline_lengthkm, flowline_areasqkm), and
# hl_reference from the hydrolocations layer.
.finalize_network <- function(ngen_gpkg, reference_gpkg, vpuid) {
  con <- DBI::dbConnect(RSQLite::SQLite(), ngen_gpkg)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  if (!DBI::dbExistsTable(con, "network")) return(invisible(NULL))

  net <- DBI::dbReadTable(con, "network")

  # Flowpath-scale mirror
  fp <- sf::st_drop_geometry(sf::st_read(ngen_gpkg, "flowpaths", quiet = TRUE))
  fp_cols <- c("flowpath_id", "mainstem_id", "hydroseq", "areasqkm",
    "lengthkm", "total_dasqkm", "slope")
  fp_cols <- intersect(fp_cols, names(fp))
  fp_mirror <- fp[, fp_cols, drop = FALSE]
  fp_mirror$flowpath_id <- as.character(fp_mirror$flowpath_id)

  net$flowpath_id <- as.character(net$flowpath_id)
  ix <- match(net$flowpath_id, fp_mirror$flowpath_id)
  for (col in setdiff(fp_cols, "flowpath_id")) {
    if (!(col %in% names(net))) net[[col]] <- fp_mirror[[col]][ix]
  }

  # Flowline-scale mirror
  fl <- sf::st_drop_geometry(sf::st_read(ngen_gpkg, "flowlines", quiet = TRUE))
  fl$flowline_id <- as.character(fl$flowline_id)
  net$flowline_id <- as.character(net$flowline_id)
  fl_ix <- match(net$flowline_id, fl$flowline_id)
  if ("flowline_hydroseq" %in% names(fl) && !("flowline_hydroseq" %in% names(net)))
    net$flowline_hydroseq <- fl$flowline_hydroseq[fl_ix]
  if ("lengthkm" %in% names(fl) && !("flowline_lengthkm" %in% names(net)))
    net$flowline_lengthkm <- fl$lengthkm[fl_ix]
  # Spec's `flowline_areasqkm` maps to our incremental catchment on flowlines.
  if ("incremental_areasqkm" %in% names(fl) &&
    !("flowline_areasqkm" %in% names(net)))
    net$flowline_areasqkm <- fl$incremental_areasqkm[fl_ix]
  # Flowline-scale slope (mirrors flowline_lengthkm), so slope is present and
  # consistent across flowpaths (`slope`), flowlines (`slope`), and network
  # (flowpath-scale `slope` + flowline-scale `flowline_slope`).
  if ("slope" %in% names(fl) && !("flowline_slope" %in% names(net)))
    net$flowline_slope <- fl$slope[fl_ix]

  # hl_reference from hydrolocations: expand pipe-delimited values so each
  # distinct hl_reference gets its own row. A poi_id with N pipe-delimited
  # references produces N copies of every network row carrying that poi_id.
  if (hfutils::layer_exists(ngen_gpkg, "hydrolocations")) {
    hl <- sf::st_drop_geometry(
      sf::st_read(ngen_gpkg, "hydrolocations", quiet = TRUE))
    if (all(c("poi_id", "hl_reference") %in% names(hl))) {
      # Reconcile net$poi_id against the pois layer (authoritative source).
      # build_aspatial_network derives poi_id from fp$poi_id, which may differ
      # from the finalized pois layer. Use pois' flowpath_id -> poi_id mapping
      # to backfill NA poi_ids in the network. The pois layer must have exactly
      # one poi_id per flowpath_id (enforced upstream); duplicates are dropped
      # here with last-wins so the code stays safe even if that invariant fails.
      if (hfutils::layer_exists(ngen_gpkg, "pois")) {
        pois_tbl <- sf::st_drop_geometry(
          sf::st_read(ngen_gpkg, "pois", quiet = TRUE))
        if (all(c("poi_id", "flowpath_id") %in% names(pois_tbl))) {
          pois_tbl$poi_id      <- as.character(pois_tbl$poi_id)
          pois_tbl$flowpath_id <- as.character(pois_tbl$flowpath_id)
          # Warn loudly if invariant is violated
          dup_fps <- pois_tbl$flowpath_id[duplicated(pois_tbl$flowpath_id)]
          if (length(dup_fps) > 0L) {
            message(sprintf(
              ".finalize_network: WARNING -- pois layer has %d flowpath_id(s) with >1 poi_id: %s",
              length(dup_fps),
              paste(utils::head(unique(dup_fps), 5), collapse = ", ")))
          }
          # One poi_id per flowpath (last-wins on duplicates)
          pois_lu <- stats::setNames(pois_tbl$poi_id, pois_tbl$flowpath_id)
          net$flowpath_id <- as.character(net$flowpath_id)
          net$poi_id      <- as.character(net$poi_id)

          # Backfill only NA poi_ids from pois
          na_mask <- is.na(net$poi_id)
          if (any(na_mask)) {
            filled <- unname(pois_lu[net$flowpath_id[na_mask]])
            net$poi_id[na_mask] <- filled
            n_filled <- sum(!is.na(filled))
            if (n_filled > 0L)
              message(sprintf(
                ".finalize_network: backfilled poi_id for %d network rows from pois layer",
                n_filled))
          }
        }
      }

      # Build poi_id -> character vector of individual hl_reference tokens
      hl_clean <- hl[!is.na(hl$hl_reference) & nzchar(hl$hl_reference), ,
        drop = FALSE]
      if (nrow(hl_clean) > 0L) {
        # Split pipe-delimited hl_reference values per hydrolocation row, then
        # aggregate all distinct tokens per poi_id.
        tokens_by_poi <- tapply(
          as.character(hl_clean$hl_reference),
          as.character(hl_clean$poi_id),
          function(v) sort(unique(trimws(unlist(strsplit(v, "|", fixed = TRUE))))))
        # tokens_by_poi: named list; names = poi_id, values = character vector

        if (!("hl_reference" %in% names(net))) net$hl_reference <- NA_character_

        # For each network row, look up its poi_id's token list and expand.
        # Rows with no poi_id or no tokens stay as-is (hl_reference = NA).
        poi_vec <- as.character(net$poi_id)

        n_tokens <- vapply(poi_vec, function(p) {
          if (is.na(p) || !nzchar(p)) return(1L)
          toks <- tokens_by_poi[[p]]
          if (is.null(toks)) 1L else length(toks)
        }, integer(1L))

        if (any(n_tokens > 1L)) {
          # Expand: replicate rows, then fill hl_reference per replicated row
          rep_idx <- rep(seq_len(nrow(net)), n_tokens)
          net_exp <- net[rep_idx, , drop = FALSE]
          rownames(net_exp) <- NULL

          flat_refs <- unlist(lapply(poi_vec, function(p) {
            if (is.na(p) || !nzchar(p)) return(NA_character_)
            toks <- tokens_by_poi[[p]]
            if (is.null(toks)) NA_character_ else toks
          }), use.names = FALSE)
          net_exp$hl_reference <- flat_refs
          net <- net_exp
        } else {
          # No expansion needed: just fill hl_reference in place
          net$hl_reference <- unlist(lapply(poi_vec, function(p) {
            if (is.na(p) || !nzchar(p)) return(NA_character_)
            toks <- tokens_by_poi[[p]]
            if (is.null(toks)) NA_character_ else toks[1L]
          }), use.names = FALSE)
        }
      }
    }
  }

  # hf_* crosswalks from the reference fabric (join on reference_id)
  if (!is.null(reference_gpkg) && file.exists(reference_gpkg) &&
    "reference_id" %in% names(net)) {
    rcon <- DBI::dbConnect(RSQLite::SQLite(), reference_gpkg)
    on.exit(DBI::dbDisconnect(rcon), add = TRUE)
    ref_cols <- DBI::dbListFields(rcon, "flowpaths")

    # Pull only what we need
    wanted <- c("flowpath_id", "flowpath_toid", "hydroseq")
    if ("source" %in% ref_cols) wanted <- c(wanted, "source")
    wanted <- intersect(wanted, ref_cols)
    q <- sprintf("SELECT %s FROM flowpaths WHERE vpuid = '%s'",
      paste(sprintf('"%s"', wanted), collapse = ","), vpuid)
    ref <- DBI::dbGetQuery(rcon, q)
    ref$flowpath_id <- as.character(ref$flowpath_id)

    r_ix <- match(as.character(net$reference_id), ref$flowpath_id)
    if (!("hf_toid" %in% names(net)) && "flowpath_toid" %in% ref_cols)
      net$hf_toid <- as.character(ref$flowpath_toid)[r_ix]
    if (!("hf_hydroseq" %in% names(net)) && "hydroseq" %in% ref_cols)
      net$hf_hydroseq <- ref$hydroseq[r_ix]
    if (!("hf_source" %in% names(net))) {
      if ("source" %in% names(ref)) {
        # per-row source if available (superconus has this)
        net$hf_source <- as.character(ref$source)[r_ix]
      } else {
        # PRVI reference lacks `source` -> fall back to the NWMv4 default
        net$hf_source <- "nhdplusv2"
      }
    }
  } else if (!("hf_source" %in% names(net))) {
    net$hf_source <- "nhdplusv2"
  }

  # Reorder columns into the spec's canonical order where possible
  spec_order <- c("flowpath_id", "flowpath_toid", "mainstem_id", "hydroseq",
    "areasqkm", "lengthkm", "total_dasqkm", "vpuid",
    "divide_id", "poi_id", "hl_reference",
    "reference_id", "hf_toid", "hf_source", "hf_hydroseq",
    "flowline_id", "flowline_toid", "is_mainstem",
    "flowline_areasqkm", "flowline_hydroseq",
    "flowline_lengthkm", "flowline_total_dasqkm")
  have <- intersect(spec_order, names(net))
  others <- setdiff(names(net), spec_order)
  net <- net[, c(have, others), drop = FALSE]

  DBI::dbWriteTable(con, "network", net, overwrite = TRUE)
}

# Aggregate hydrolocations per poi_id into summary rollup fields on pois.
.add_poi_rollups <- function(gpkg) {
  if (!hfutils::layer_exists(gpkg, "pois") ||
      !hfutils::layer_exists(gpkg, "hydrolocations")) return(invisible(NULL))
  hl   <- sf::st_drop_geometry(sf::st_read(gpkg, "hydrolocations", quiet = TRUE))
  pois <- sf::st_read(gpkg, "pois", quiet = TRUE)
  # Either layer can exist and still be empty, which the rollup join cannot use.
  if (nrow(hl) == 0L || nrow(pois) == 0L) return(invisible(NULL))

  uniq_join <- function(v) paste(sort(unique(v[!is.na(v) & nzchar(v)])),
    collapse = ",")

  # hl_class / hl_type may be absent from the hydrolocations layer. Fall back
  # to NA vectors of the right length so tapply doesn't complain about
  # length mismatch.
  hl_cls <- if ("hl_class" %in% names(hl)) as.character(hl$hl_class)
  else rep(NA_character_, nrow(hl))
  hl_typ <- if ("hl_type"  %in% names(hl)) as.character(hl$hl_type)
  else rep(NA_character_, nrow(hl))

  # tapply returns a 1D array; as.character strips the class so the column
  # writes as plain character (GeoPackage can't serialize arrays). Compute each
  # grouped reduction once and read both its values and (poi_id) names off it.
  r_cls <- tapply(hl_cls, hl$poi_id, uniq_join)
  r_typ <- tapply(hl_typ, hl$poi_id, uniq_join)
  classes <- as.character(r_cls)
  types   <- as.character(r_typ)
  counts  <- table(hl$poi_id)
  names(classes) <- names(r_cls)
  names(types)   <- names(r_typ)

  m <- match(as.character(pois$poi_id), names(counts))
  pois$hl_count   <- as.integer(counts[m])
  pois$hl_classes <- unname(classes[as.character(pois$poi_id)])
  pois$hl_types   <- unname(types[as.character(pois$poi_id)])
  pois$hl_classes[is.na(pois$hl_classes) | !nzchar(pois$hl_classes)] <- NA_character_
  pois$hl_types[is.na(pois$hl_types)     | !nzchar(pois$hl_types)]   <- NA_character_

  pois <- .drop_list_cols(pois)
  sf::st_write(pois, gpkg, "pois", delete_layer = TRUE, quiet = TRUE)
}

# Report the drainage area a gage POI is delineated to, against the reference
# drainage area of the reaches it is built from.
#
# Measurement only, with no tolerance and no flag. `ref` is the largest
# reference DA among the flowpath's member reaches, so it depends on how many
# reaches the aggregation bundled rather than on where the gage sits: the same
# gage passes or fails a fixed threshold according to unit size, with its
# delineated area unchanged. A threshold that moves with unit size does not
# measure placement, so the areas are published and the judgement is left to
# the reader.
#
# `poi_da_ref_sqkm`  largest reference DA among the POI's member reaches
# `poi_da_ngen_sqkm` the delineated DA of the flowpath the POI is an outlet of
# `poi_da_ratio`     |ref - ngen| / max(ref, 1), NA when ref is below one DEM
#                    cell (9e-4 km2), where a reference DA is absent rather than
#                    small
.report_poi_drainage_area <- function(gpkg) {
  if (!hfutils::layer_exists(gpkg, "network")) return(invisible(NULL))
  net <- tryCatch(sf::st_drop_geometry(sf::st_read(gpkg, "network", quiet = TRUE)),
    error = function(e) NULL)
  if (is.null(net)) return(invisible(NULL))
  need <- c("hl_reference", "flowpath_id", "total_dasqkm", "flowline_total_dasqkm")
  if (!all(need %in% names(net))) return(invisible(NULL))

  # CONUS references a streamgage as `nwis-01027200`; the oCONUS domains use
  # `gages-15239800`. A row can carry several references, so match the prefix at
  # a boundary rather than anchoring the whole string.
  is_gage <- !is.na(net$hl_reference) &
    grepl("(^|\\|)(nwis|gages)-", net$hl_reference)
  if (!any(is_gage)) return(invisible(NULL))
  g <- net[is_gage, , drop = FALSE]

  fmax <- function(v, by) tapply(suppressWarnings(as.numeric(v)), by,
    function(x) suppressWarnings(max(x, na.rm = TRUE)))
  ref_by  <- fmax(g$flowline_total_dasqkm, g$hl_reference)
  ngen_by <- fmax(g$total_dasqkm,          g$hl_reference)
  hlr  <- names(ref_by)
  ref  <- as.numeric(ref_by)
  ngen <- as.numeric(ngen_by[hlr])

  .DA_FLOOR_SQKM <- 1e-3
  usable <- is.finite(ref) & ref >= .DA_FLOOR_SQKM
  ratio <- abs(ref - ngen) / pmax(ref, 1)
  ratio[!is.finite(ratio) | !usable] <- NA_real_

  qc <- data.frame(hl_reference = hlr,
    poi_da_ref_sqkm = round(ref, 2), poi_da_ngen_sqkm = round(ngen, 2),
    poi_da_ratio = round(ratio, 3), stringsAsFactors = FALSE)
  attached <- .attach_poi_da_qc(gpkg, qc)
  message(sprintf(paste0("finalize_data_model: drainage-area comparison written ",
    "for %d gage POI(s)%s; %d had no usable reference area (< %.0e km2)"),
    sum(usable), if (length(attached)) sprintf(" to %s",
      paste(attached, collapse = ", ")) else "", sum(!usable), .DA_FLOOR_SQKM))
  invisible(qc)
}

# Attach the drainage-area QC to whichever POI layer can actually be joined.
#
# The join key is `hl_reference`, which `hydrolocations` carries and `pois` does
# not, so assuming a layer silently writes nothing. Pick the layers that hold the
# key, and say so when none can be joined, a QC computed and dropped is worse
# than one never run, because the caller reports it as done.
.attach_poi_da_qc <- function(gpkg, qc) {
  cand <- c("hydrolocations", "pois")
  cand <- cand[vapply(cand, function(l) hfutils::layer_exists(gpkg, l), logical(1))]
  attached <- character(0)
  for (poi_layer in cand) {
    pl <- sf::st_read(gpkg, poi_layer, quiet = TRUE)
    if (!"hl_reference" %in% names(pl)) next
    m <- match(as.character(pl$hl_reference), qc$hl_reference)
    if (!any(!is.na(m))) next
    for (cc in setdiff(names(qc), "hl_reference")) pl[[cc]] <- qc[[cc]][m]
    sf::st_write(.drop_list_cols(pl), gpkg, poi_layer, delete_layer = TRUE, quiet = TRUE)
    attached <- c(attached, poi_layer)
  }
  if (!length(attached))
    message(sprintf(paste0(
      "finalize_data_model: drainage-area QC computed for %d POI(s) but no layer ",
      "carries a matching `hl_reference` to join it to (%s) -- flags NOT written"),
      nrow(qc), if (length(cand)) paste(cand, collapse = ", ") else "no POI layer"))
  attached
}

# Replace the VPU-local mainstem_id with the reference level-path id recorded
# beside it in orig_levelpathID, then drop that column. The local ids stay if the
# record is absent or does not map one-to-one, so no river is merged or split.
.national_mainstem_ids <- function(gpkg) {
  if (!hfutils::layer_exists(gpkg, "flowlines")) return(invisible(NULL))
  con <- DBI::dbConnect(RSQLite::SQLite(), gpkg)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  if (!all(c("mainstem_id", "orig_levelpathID") %in% DBI::dbListFields(con, "flowlines"))) {
    .hf_log_warn("national mainstem ids: flowlines carry no orig_levelpathID; keeping local ids")
    return(invisible(NULL))
  }
  nat <- DBI::dbGetQuery(con, 'SELECT DISTINCT mainstem_id AS local,
    CAST("orig_levelpathID" AS INTEGER) AS lp FROM flowlines
    WHERE mainstem_id IS NOT NULL AND "orig_levelpathID" IS NOT NULL')
  n_local <- DBI::dbGetQuery(con,
    "SELECT COUNT(DISTINCT mainstem_id) n FROM flowlines WHERE mainstem_id IS NOT NULL")$n
  lyrs <- Filter(function(l) DBI::dbExistsTable(con, l) &&
    "mainstem_id" %in% DBI::dbListFields(con, l), c("flowpaths", "flowlines", "network"))
  DBI::dbDisconnect(con); on.exit()

  if (nrow(nat) != n_local || anyDuplicated(nat$local) || anyDuplicated(nat$lp)) {
    .hf_log_warn(sprintf(paste("national mainstem ids: %d local ids do not map",
      "one-to-one onto orig_levelpathID; keeping local ids"), n_local))
  } else {
    for (lyr in lyrs)
      hfutils::gpkg_update_col(gpkg, lyr, "mainstem_id", nat$local, "mainstem_id",
                               format(nat$lp, scientific = FALSE, trim = TRUE))
  }
  con <- DBI::dbConnect(RSQLite::SQLite(), gpkg)
  DBI::dbExecute(con, 'ALTER TABLE flowlines DROP COLUMN "orig_levelpathID"')
  DBI::dbDisconnect(con)
  invisible(nat)
}
