#' Build an aspatial long-form network crosswalk
#'
#' @description
#' Replaces the spatial `network` layer in a hydrofabric GeoPackage with an
#' aspatial attribute table keyed on `(flowline_id, reference_id)`, one row
#' per sub-flowline x reference-COMID member. The table is rebuilt from the
#' flowlines layer and the reconciled-group memberships so it always matches
#' the trimmed topology produced by `extend_flowpaths_to_mainstem`:
#' flowlines on a flowpath's mainstem levelpathid keep the mainstem's
#' `flowpath_id`; flowlines on absorbed tributaries (different levelpathid)
#' report a synthetic `fp-t<levelpathid>` identity and inherit the mainstem's
#' downstream nexus as their `flowpath_toid`.
#'
#' Output columns:
#' - `flowpath_id`, `flowpath_toid`        : aggregated flowpath topology
#' - `mainstem_id`                         : levelpathid of the flowpath's
#'                                            mainstem; group all flowlines
#'                                            on the same mainstem channel
#' - `is_mainstem`                         : TRUE when this flowline is on
#'                                            the flowpath's mainstem
#'                                            (levelpathid == mainstem_id);
#'                                            FALSE for absorbed tributaries
#' - `flowline_id`,  `flowline_toid`       : split-segment topology
#' - `reference_id`                        : originating reference COMID
#' - `reference_percent`                   : share of reference NHD length
#'                                            represented by this flowline's
#'                                            reconciled group
#' - `poi_id`                              : aggregation POI (nullable)
#' - `hydrolocation_id`                    : `hl_id` from hydrolocations
#' - `divide_id`, `divide_areasqkm`        : divide this flowpath drains
#'                                            (all flowlines in a flowpath,
#'                                            mainstem and tributary: carry
#'                                            the same `divide_id`)
#' - `flowline_area`                       : drainage area at flowline
#'                                            outlet (from FAC extraction)
#'
#' @param gpkg Path to the ngen GeoPackage produced by `nextgen_from_gpkg`.
#' @param refactored `sf` of refactored flowlines from `refactor_flowlines`
#'   (`out$refactored`). Used to compute `reference_percent` per reference COMID.
#' @param reconciled `sf` of reconciled flowlines (`out$reconciled`): must
#'   contain `ID` and `member_flowpath_id`. Used to map each flowline's
#'   levelpathid back to the reference COMIDs it represents.
#' @param hydrolocations `sf` of hydrolocations with `poi_id` and `hl_id`,
#'   or path to a gpkg with a `hydrolocations` layer. If `NULL`,
#'   `hydrolocation_id` is sourced from the `events` layer if present.
#' @param verbose If TRUE, emit a rich breakdown: row class counts,
#'   tributary-flowline share, reference_percent distribution, hydrolocation
#'   match rate, divide-area / flowline-area summaries. Default FALSE.
#'
#' @return Invisibly, the path to `gpkg`.
#'
#' @examples
#' \dontrun{
#' build_aspatial_network(
#'   gpkg           = "ngen_vpu01.gpkg",
#'   refactored     = sf::st_read("refactored.gpkg", "refactored"),
#'   reconciled     = sf::st_read("reconciled.gpkg", "reconciled"),
#'   hydrolocations = "hydrolocations.gpkg")
#' }
#'
#' @export
build_aspatial_network <- function(gpkg,
                                   refactored     = NULL,
                                   reconciled     = NULL,
                                   hydrolocations = NULL,
                                   verbose        = FALSE) {
  fl <- sf::st_drop_geometry(sf::st_read(gpkg, "flowlines",  quiet = TRUE))
  fp <- sf::st_drop_geometry(sf::st_read(gpkg, "flowpaths", quiet = TRUE))
  div <- sf::st_drop_geometry(sf::st_read(gpkg, "divides",  quiet = TRUE))

  # Per-flowpath mainstem levelpathid (used for mainstem_id + is_mainstem flag)
  mainstem_col   <- if ("mainstem" %in% names(fp)) "mainstem" else                    if ("mainstem_id" %in% names(fp)) "mainstem_id" else NULL
  mainstem_by_fp <- if (!is.null(mainstem_col))
    stats::setNames(as.character(fp[[mainstem_col]]), as.character(fp$flowpath_id))
  else
    stats::setNames(rep(NA_character_, nrow(fp)), as.character(fp$flowpath_id))

  # Reconciled-line membership: map each reconciled ID to ITS OWN reference
  # COMIDs (member_flowpath_id). Key by the reconciled ID, NOT levelpathid.
  # A length-split flowline fragment must inherit only the COMIDs of the
  # reconciled line it was split from. Keying by levelpathid unions every
  # reconciled row on the levelpath (i.e. the whole mainstem channel), so a
  # sub-km fragment ends up attributed the entire mainstem's COMID set
  # (e.g. a 0.49 km fragment on a 155 km levelpath inheriting all 149 COMIDs).
  # Each fragment carries its parent reconciled ID via `reconciled_id` (set by
  # nextgen_from_gpkg); it equals reconciled$ID even for lines absorbed during
  # Stage-2 (add_flowlines preserves orig_fl_id through the toid-walk).
  rec_to_members <- list()
  if (!is.null(reconciled) && nrow(reconciled) > 0L) {
    rec <- sf::st_drop_geometry(reconciled)
    mem_col <- if ("member_flowpath_id" %in% names(rec))
      "member_flowpath_id" else "member_comid"
    rec_ids <- as.character(rec[["ID"]])
    mems    <- strsplit(as.character(rec[[mem_col]]), ",", fixed = TRUE)
    rec_to_members <- lapply(
      split(seq_along(rec_ids), rec_ids),
      function(i) unique(trimws(unlist(mems[i])))
    )
  }

  # Helper: normalize a refactored_id ("NHD" or "NHD.N") -> base reference COMID
  base_of <- function(x) format(suppressWarnings(as.numeric(
    sub("\\..*$", "", as.character(x)))), scientific = FALSE, trim = TRUE)

  # Refactored lengths for reference_percent
  pct_lookup <- NULL
  if (!is.null(refactored) && nrow(refactored) > 0L) {
    rfc <- sf::st_drop_geometry(refactored)
    rfc_base <- base_of(rfc$flowpath_id)
    nhd_tot  <- tapply(as.numeric(rfc$lengthkm), rfc_base, sum)
    # Per refactored part length (for the part->reference percent)
    # Used when we can tell the specific part; otherwise we sum parts later.
    rfc_key <- paste0(rfc_base,
      ifelse(grepl("\\.", as.character(rfc$flowpath_id)),
        paste0(".", sub(".*\\.", "", as.character(rfc$flowpath_id))),
        ""))
    rfc_len <- stats::setNames(as.numeric(rfc$lengthkm), rfc_key)
    pct_lookup <- list(nhd_tot = nhd_tot, rfc_len = rfc_len)
  }

  # Normalize flowlines cols we need
  fl$flowline_id   <- as.integer(fl$flowline_id)
  fl$flowline_toid <- suppressWarnings(as.integer(fl$flowline_toid))
  fl$levelpathid   <- as.character(fl$levelpathid)
  fl$flowpath_id   <- as.character(fl$flowpath_id)
  fl$flowpath_toid <- as.character(fl$flowpath_toid)
  fl$reconciled_id <- if ("reconciled_id" %in% names(fl)) {
    as.character(fl$reconciled_id)
  } else if ("rec_id" %in% names(fl)) {
    as.character(fl$rec_id)
  } else {
    NA_character_
  }

  # For each flowline, build its row(s), one per reference COMID of the
  # reconciled line the fragment was split from (keyed by reconciled_id, so a
  # fragment gets only its own line's members, not the whole levelpath's).
  members_per_row <- rec_to_members[fl$reconciled_id]
  # Flowlines whose levelpathid isn't in the lookup get a single NA-member
  # row so they're not silently dropped from the lookup table.
  n_members <- vapply(members_per_row, length, integer(1L))
  n_members[n_members == 0L] <- 1L

  # divide_id and mainstem_id are flowpath-level; derive per flowline via
  # the flowpath mapping. Both are set before any is_mainstem logic so that
  # tributary flowlines carry the same divide_id as their parent flowpath.
  div_by_fp  <- stats::setNames(as.character(fp$divide_id),
    as.character(fp$flowpath_id))
  fl$divide_id   <- unname(div_by_fp[fl$flowpath_id])
  fl$mainstem_id <- unname(mainstem_by_fp[fl$flowpath_id])

  # is_mainstem: a flowline is on the mainstem of its flowpath when its
  # flowline_id appears in fp$rec_id_list (the CSV list of mainstem rec_ids
  # written by Stage-2 aggregation). This is the only reliable signal,
  # levelpathid on flowlines uses sequential IDs from add_flowlines while
  # mainstem_id on flowpaths carries NHDPlus levelpathids; they cannot be
  # compared directly.
  fp_fl_col <- if ("rec_id_list" %in% names(fp)) "rec_id_list"
  else if ("flowline_id" %in% names(fp)) "flowline_id"
  else NULL
  if (!is.null(fp_fl_col)) {
    ms_ids_per_fp <- strsplit(as.character(fp[[fp_fl_col]]), ",", fixed = TRUE)
    ms_set <- as.integer(unlist(ms_ids_per_fp))
    fl$is_mainstem <- as.integer(fl$flowline_id) %in% ms_set
  } else {
    fl$is_mainstem <- NA
  }

  long <- data.frame(
    flowpath_id   = rep(fl$flowpath_id,       n_members),
    flowpath_toid = rep(fl$flowpath_toid,     n_members),
    mainstem_id   = rep(fl$mainstem_id,       n_members),
    is_mainstem   = rep(fl$is_mainstem,       n_members),
    flowline_id   = rep(fl$flowline_id,       n_members),
    flowline_toid = rep(fl$flowline_toid,     n_members),
    flowline_area = rep(as.numeric(fl$tot_drainage_areasqkm), n_members),
    divide_id     = rep(fl$divide_id,         n_members),
    levelpathid   = rep(fl$levelpathid,       n_members),
    stringsAsFactors = FALSE)
  long$ref_raw <- unlist(lapply(seq_along(members_per_row), function(i) {
    m <- members_per_row[[i]]
    if (length(m) == 0L) NA_character_ else trimws(m)
  }))

  long$reference_id <- base_of(long$ref_raw)
  long$reference_id[is.na(long$ref_raw)] <- NA_character_

  # reference_percent: per (reference_id): use refactored_id length if the
  # member lists split parts (e.g. "NHD.1"); otherwise fall back to 100%.
  if (!is.null(pct_lookup)) {
    # Build key per row: base + optional ".N"
    row_part <- ifelse(grepl("\\.", long$ref_raw),
      paste0(".", sub(".*\\.", "", long$ref_raw)), "")
    row_key  <- ifelse(is.na(long$reference_id), NA_character_,
      paste0(long$reference_id, row_part))
    this_len <- unname(pct_lookup$rfc_len[row_key])
    nhd_tot  <- unname(pct_lookup$nhd_tot[long$reference_id])
    long$reference_percent <- ifelse(
      !is.na(this_len) & !is.na(nhd_tot) & nhd_tot > 0,
      round(100 * this_len / nhd_tot, 4),
      ifelse(!is.na(long$reference_id) & !is.na(nhd_tot), 100, NA_real_))
  } else {
    long$reference_percent <- NA_real_
  }

  # Divide area (per divide_id)
  div_area <- stats::setNames(as.numeric(div$areasqkm), as.character(div$divide_id))
  long$divide_areasqkm <- unname(div_area[long$divide_id])

  # POI per flowpath (from the flowpaths layer)
  poi_lookup <- stats::setNames(as.character(fp$poi_id),
    as.character(fp$flowpath_id))
  long$poi_id <- unname(poi_lookup[long$flowpath_id])

  # hydrolocation_id from hydrolocations (or events fallback)
  hl_lookup <- NULL
  if (!is.null(hydrolocations)) {
    if (is.character(hydrolocations) && file.exists(hydrolocations)) {
      hydrolocations <- sf::st_read(hydrolocations, "hydrolocations",
        quiet = TRUE)
    }
    if (inherits(hydrolocations, "sf") &&
      all(c("poi_id", "hl_id") %in% names(hydrolocations))) {
      hl_lookup <- stats::setNames(
        as.character(hydrolocations$hl_id),
        as.character(hydrolocations$poi_id))
      hl_lookup <- hl_lookup[!duplicated(names(hl_lookup))]
    }
  }
  if (is.null(hl_lookup) && hfutils::layer_exists(gpkg, "events")) {
    ev <- sf::st_drop_geometry(sf::st_read(gpkg, "events", quiet = TRUE))
    if (all(c("poi_id", "hl_id") %in% names(ev))) {
      ev <- ev[!is.na(ev$poi_id) & !is.na(ev$hl_id), c("poi_id", "hl_id"),
        drop = FALSE]
      hl_lookup <- stats::setNames(as.character(ev$hl_id),
        as.character(ev$poi_id))
      hl_lookup <- hl_lookup[!duplicated(names(hl_lookup))]
    }
  }
  long$hydrolocation_id <- if (!is.null(hl_lookup))
    unname(hl_lookup[long$poi_id]) else NA_character_

  n_trib <- sum(!long$is_mainstem & !is.na(long$flowpath_id), na.rm = TRUE)

  # Before dropping NA-flowpath rows, remap any valid-fp flowlines whose
  # flowline_toid pointed to a NA-fp row by following the actual NHD
  # flowline_toid chain until reaching a valid-fp flowline.
  valid_fp     <- !is.na(long$flowpath_id)
  valid_fl_ids <- as.character(long$flowline_id[valid_fp])
  na_rows_lkp  <- long[!valid_fp & !is.na(long$flowline_toid),
    c("flowline_id", "flowline_toid"), drop = FALSE]
  if (nrow(na_rows_lkp) > 0L) {
    na_set <- as.character(na_rows_lkp$flowline_id)

    na_lkp <- stats::setNames(as.character(na_rows_lkp$flowline_toid), na_set)
    for (iter in seq_len(50L)) {
      still_na <- na_lkp %in% na_set
      if (!any(still_na, na.rm = TRUE)) break
      na_lkp[still_na] <- na_lkp[na_lkp[still_na]]
    }
    na_lkp <- na_lkp[na_lkp %in% valid_fl_ids]

    needs_remap <- which(valid_fp & (as.character(long$flowline_toid) %in% na_set))
    if (length(needs_remap) > 0L) {
      cur_toids_r <- as.character(long$flowline_toid[needs_remap])
      fb_vals     <- na_lkp[cur_toids_r]
      long$flowline_toid[needs_remap] <- ifelse(!is.na(fb_vals), as.integer(fb_vals), 0L)
    }

    # Attribute NA-fp flowlines' COMIDs to their nearest downstream valid-fp.
    # This recovers absorbed tributary COMIDs that lost their flowpath_id during
    # Stage-2 collapse but can be traced downstream to a surviving fp.
    if (length(na_lkp) > 0L) {
      na_long_idx <- which(!valid_fp &
        as.character(long$flowline_id) %in% names(na_lkp) &
        !is.na(long$ref_raw))
      if (length(na_long_idx) > 0L) {
        ds_fl <- na_lkp[as.character(long$flowline_id[na_long_idx])]
        ok    <- !is.na(ds_fl)
        if (any(ok)) {
          vfp  <- long[valid_fp, , drop = FALSE]
          ufp  <- vfp[!duplicated(as.character(vfp$flowline_id)), ]
          fp_lkp   <- stats::setNames(as.character(ufp$flowpath_id),   as.character(ufp$flowline_id))
          toid_lkp <- stats::setNames(as.character(ufp$flowpath_toid), as.character(ufp$flowline_id))
          ms_lkp   <- stats::setNames(as.character(ufp$mainstem_id),   as.character(ufp$flowline_id))
          div_lkp  <- stats::setNames(as.character(ufp$divide_id),     as.character(ufp$flowline_id))
          idx <- na_long_idx[ok]
          ds  <- ds_fl[ok]
          long$flowpath_id[idx]   <- fp_lkp[ds]
          long$flowpath_toid[idx] <- toid_lkp[ds]
          long$mainstem_id[idx]   <- ms_lkp[ds]
          long$divide_id[idx]     <- div_lkp[ds]
          long$is_mainstem[idx]   <- FALSE
        }
      }
    }
  }
  long <- long[!is.na(long$flowpath_id), , drop = FALSE]

  out <- unique(data.frame(
    flowpath_id       = long$flowpath_id,
    flowpath_toid     = long$flowpath_toid,
    mainstem_id       = long$mainstem_id,
    flowline_id       = long$flowline_id,
    is_mainstem       = long$is_mainstem,
    flowline_toid     = long$flowline_toid,
    reference_id      = long$reference_id,
    reference_percent = long$reference_percent,
    poi_id            = long$poi_id,
    hydrolocation_id  = long$hydrolocation_id,
    divide_id         = long$divide_id,
    divide_areasqkm   = long$divide_areasqkm,
    flowline_area     = long$flowline_area,
    stringsAsFactors  = FALSE))

  # Delete the spatial network layer, write aspatial via DBI
  con <- DBI::dbConnect(RSQLite::SQLite(), gpkg)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  tryCatch(
    {
      DBI::dbExecute(con,
        "DELETE FROM gpkg_geometry_columns WHERE table_name = 'network'")
      DBI::dbExecute(con,
        "DELETE FROM gpkg_contents WHERE table_name = 'network'")
      DBI::dbExecute(con,
        "DELETE FROM gpkg_ogr_contents WHERE table_name = 'network'")
    },
    error = function(e) invisible(NULL))
  for (nm in c("rtree_network_geom",
    "rtree_network_geom_node",
    "rtree_network_geom_parent",
    "rtree_network_geom_rowid")) {
    if (DBI::dbExistsTable(con, nm)) DBI::dbRemoveTable(con, nm)
  }
  if (DBI::dbExistsTable(con, "network")) {
    DBI::dbRemoveTable(con, "network")
  }
  DBI::dbWriteTable(con, "network", out, overwrite = TRUE)
  DBI::dbExecute(con, paste0(
    "INSERT OR REPLACE INTO gpkg_contents ",
    "(table_name, data_type, identifier) ",
    "VALUES ('network', 'attributes', 'network')"))

  n_ms  <- if ("is_mainstem" %in% names(out)) sum(out$is_mainstem, na.rm = TRUE) else nrow(out) - n_trib
  n_trib_out <- nrow(out) - n_ms
  message(sprintf(
    "build_aspatial_network: wrote %d rows (%d mainstem, %d tributary)",
    nrow(out), n_ms, n_trib_out))

  if (isTRUE(verbose)) {
    message("-- build_aspatial_network detail --------------------------")
    message(sprintf("  rows: %d  (mainstem: %d, tributary: %d)",
      nrow(out), n_ms, n_trib_out))
    message(sprintf("  distinct flowpath_ids: %d   distinct flowline_ids: %d   distinct reference NHDs: %d",
      length(unique(stats::na.omit(out$flowpath_id))),
      length(unique(out$flowline_id)),
      length(unique(stats::na.omit(out$reference_id)))))

    rp <- out$reference_percent
    rp <- rp[!is.na(rp)]
    if (length(rp) > 0L) {
      buckets <- cut(rp, breaks = c(-0.01, 25, 50, 75, 99.99, 100.01),
        labels = c("0-25%", "25-50%", "50-75%", "75-100% (split)", "100% (unsplit)"))
      tab <- table(buckets)
      message("  reference_percent distribution:")
      for (nm in names(tab)) message(sprintf("    %-18s %6d", nm, tab[[nm]]))
    }

    if ("hydrolocation_id" %in% names(out)) {
      hl_rows <- sum(!is.na(out$hydrolocation_id))
      poi_rows <- sum(!is.na(out$poi_id))
      message(sprintf("  hydrolocation_id populated on %d / %d POI-bearing rows (%.1f%%)",
        hl_rows, poi_rows,
        if (poi_rows > 0) 100 * hl_rows / poi_rows else 0))
    }

    if ("divide_areasqkm" %in% names(out)) {
      da <- stats::na.omit(unique(data.frame(
        divide_id = out$divide_id, a = out$divide_areasqkm,
        stringsAsFactors = FALSE)))$a
      fa <- stats::na.omit(unique(data.frame(
        flowline_id = out$flowline_id, a = out$flowline_area,
        stringsAsFactors = FALSE)))$a
      if (length(da) > 0L)
        message(sprintf("  divide_areasqkm -- min: %.3f  median: %.3f  max: %.3f  (n=%d)",
          min(da), stats::median(da), max(da), length(da)))
      if (length(fa) > 0L)
        message(sprintf("  flowline_area   -- min: %.3f  median: %.3f  max: %.3f  (n=%d)",
          min(fa), stats::median(fa), max(fa), length(fa)))
    }
    message("------------------------------------------------------------")
  }
  invisible(gpkg)
}
