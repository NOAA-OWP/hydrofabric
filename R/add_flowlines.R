# Which grid a fallback FAC value came from, for the `fac_source` column.
#
# Matched on the distributed file and directory names, the same way the grid
# credits are, so a build records the product it read rather than a fixed label.
# An unrecognised path is reported as "fallback": accurate, if less specific,
# and it keeps the column honest about the value not being primary.
.fac_source_label <- function(path) {
  if (!is.character(path) || !length(path) || !nzchar(path[1L])) return("fallback")
  p <- tolower(path[1L])
  if (grepl("gedtm", p, fixed = TRUE)) return("gedtm30")
  if (grepl("hydrosheds", p, fixed = TRUE) || grepl("_acc_1s_|_dir_1s_", p))
    return("hydrosheds")
  "fallback"
}

# Flow accumulation is a non-negative cell count, so a negative reading is a void
# the source encoded with a sentinel instead of its declared NoData. Treat it as
# missing so it is filled or left NA, never carried through as data.
.fac_drop_sentinels <- function(x) {
  bad <- !is.na(x) & x < 0
  n <- sum(bad)
  if (n > 0L) {
    x[bad] <- NA_real_
    .hf_log_warn(sprintf(
      "add_flowlines: %d negative FAC value(s) treated as NoData", n))
  }
  x
}

#' Add Sub-Flowline Segments to Aggregated Hydrofabric
#'
#' @description
#' Takes the full-resolution reconciled flowlines (as passed to
#' `aggregate_to_distribution`), splits any that exceed `split_flines_meters`,
#' assigns each segment its parent aggregated `flowpath_id` and `divide_id`
#' via the `network` lookup table in `input_gpkg`, remaps topology so the
#' resulting segments form a valid directed graph, and attributes each segment
#' with total and incremental contributing area derived from a
#' flow-accumulation raster (FAC).
#'
#' If `events` is supplied (the full event set from
#' `rf_generate_events_from_hydrolocations_fast`), each event is assigned the
#' aggregated `flowpath_id` and `divide_id` of the reconciled flowline it
#' belongs to, and written as an `events` layer to `input_gpkg`.
#'
#' @details
#' Workflow:
#' 1. Read the `network` table from `input_gpkg` and expand its
#'    comma-separated `flowline_id` column to build a per-ID lookup:
#'    reconciled flowline ID -> aggregated `flowpath_id`, `flowpath_toid`,
#'    `divide_id`.
#' 2. Split `flowlines` segments that exceed `split_flines_meters` using
#'    `split_flowlines`.
#' 3. Assign `flowpath_id` / `divide_id` to each segment via the lookup.
#' 4. Remap topology: segments within the same parent flowpath chain
#'    upstream-to-downstream; the terminal segment of each chain points to
#'    the first segment of the downstream chain (by `flowpath_toid`).
#' 5. Extract the FAC value at the downstream endpoint of each segment to
#'    derive `tot_drainage_areasqkm` (FAC * cell area in km2).
#' 6. Compute `areasqkm` (incremental) as `tot_drainage_areasqkm` of this
#'    segment minus the maximum `tot_drainage_areasqkm` of immediate upstream
#'    segments.  Head segments with no upstream neighbour use their own total.
#' 7. If `poi_events` and/or `junction_events` are provided, each is
#'    attributed with the aggregated `flowpath_id` and `divide_id`, then
#'    the combined set is written as an `events` layer in `input_gpkg`.
#'    POI events are mapped via their `reconciled_id` (integer reconciled
#'    flowline ID from `refactor_flowlines`) directly through the network lookup.
#'    Junction events are mapped via their NHD `flowpath_id` column through
#'    `flowlines$member_flowpath_id` to the reconciled integer ID.
#' 8. Write flowline segments as a `flowlines` layer to `input_gpkg` and
#'    return the `sf` object invisibly.
#'
#' @param flowlines  An `sf` object of full-resolution reconciled flowlines
#'   (the data passed to `aggregate_to_distribution`).  Must have integer
#'   columns `ID` and `toID`, and a character column `member_flowpath_id`
#'   (comma-separated NHD IDs per reconciled flowline).
#' @param input_gpkg  Path to a GeoPackage produced by
#'   `aggregate_to_distribution` (must contain a `network` table).  The
#'   `flowlines` (and optionally `events`) layer will be written here.
#' @param split_flines_meters Numeric.  Maximum individual segment length in
#'   metres.  Default `2000`.
#' @param fac Character path to a flow-accumulation GeoTIFF, or a
#'   [terra::SpatRaster] already in memory.
#' @param fallback_fac Optional character path (or `NULL`) to a fallback FAC
#'   raster used when `fac` returns NA after the standard
#'   2-cell buffer fill.  Values filled from the fallback are tagged
#'   `fac_source` naming that grid in the output.  When `fallback_fac` resolves to
#'   the same file as `fac` (path equality, i.e. a domain configured with the
#'   fallback grid as its primary), the fallback lookup is skipped
#'   automatically.  Default `NULL`.
#' @param poi_events Optional `sf` object of snapped POI events returned by
#'   `refactor_flowlines` (i.e. `out$events`).  Must contain a `reconciled_id`
#'   column (integer reconciled flowline ID) and `poi_id`.  These events
#'   carry geometries repositioned onto the reconciled network.
#' @param junction_events Optional `sf` object of junction events from
#'   `rf_generate_events_from_hydrolocations_fast` (i.e.
#'   `filter(events, event_kind == "junction")`).  Must contain a
#'   `flowpath_id` column with original NHD flowline IDs.
#'
#' @return Invisibly, the `sf` data frame of flowline segments.  The layer
#'   `flowlines` is written/overwritten in `input_gpkg`.  If `events` is
#'   provided, an `events` layer is also written/overwritten.
#'
#' @examples
#' \dontrun{
#' reconciled <- sf::st_read("reconciled.gpkg", "reconciled")
#' add_flowlines(flowlines = reconciled,
#'               input_gpkg = "out_vpu01.gpkg",
#'               split_flines_meters = 2000,
#'               fac = "fac.tif")
#' }
#'
#' @export
#'
#' @importFrom sf st_read st_drop_geometry st_transform st_crs st_length st_write
#' @importFrom dplyr select any_of
#' @importFrom terra rast res vect extract
#' @importFrom lwgeom st_endpoint


add_flowlines <- function(flowlines,
                          input_gpkg,
                          split_flines_meters = 2000,
                          fac,
                          fallback_fac    = NULL,
                          poi_events      = NULL,
                          junction_events = NULL) {
  # ---- 1. Build lookup from network table ------------------------------------
  # network$rec_id_list is comma-separated reconciled IDs, e.g. "1,2" or "9,10"
  net <- sf::st_read(input_gpkg, layer = "network", quiet = TRUE)
  net$flowpath_id   <- as.character(net$flowpath_id)
  net$flowpath_toid <- as.character(net$flowpath_toid)
  net$divide_id     <- as.character(net$divide_id)

  # Expand: one row per individual flowline ID. Single-allocation build via
  # vectorized strsplit + rep: avoids an O(n^2) iterative rbind.
  splits_fl <- strsplit(as.character(net$rec_id_list), ",", fixed = TRUE)
  lens_fl   <- lengths(splits_fl)
  net_expanded <- data.frame(
    flowline_id   = as.integer(trimws(unlist(splits_fl, use.names = FALSE))),
    flowpath_id   = rep(net$flowpath_id,   lens_fl),
    flowpath_toid = rep(net$flowpath_toid, lens_fl),
    divide_id     = rep(net$divide_id,     lens_fl),
    stringsAsFactors = FALSE)
  net_expanded <- net_expanded[!duplicated(net_expanded$flowline_id), ]

  fl_to_fp     <- stats::setNames(net_expanded$flowpath_id,   net_expanded$flowline_id)
  fl_to_fptoid <- stats::setNames(net_expanded$flowpath_toid, net_expanded$flowline_id)
  fl_to_div    <- stats::setNames(net_expanded$divide_id,     net_expanded$flowline_id)

  # ---- 2. Load FAC -----------------------------------------------------------
  fac_r         <- if (is.character(fac)) terra::rast(fac) else fac
  fac_res_m     <- terra::res(fac_r)[1]
  fac_cell_sqkm <- (fac_res_m^2) / 1e6

  # ---- 3. Prepare flowlines for splitting ------------------------------------
  # split_flowlines needs: ID (id), toid/flowpath_toid/toCOMID (downstream),
  # lengthkm.  Reconciled flowlines use toID: rename to toid.
  flines        <- flowlines
  # Normalize: rec_id/rec_toid (on-disk names) -> ID/toID (internal convention).
  if ("rec_id"   %in% names(flines)) names(flines)[names(flines) == "rec_id"]   <- "ID"
  if ("rec_toid" %in% names(flines)) names(flines)[names(flines) == "rec_toid"] <- "toID"
  flines$ID     <- as.integer(flines$ID)
  names(flines)[names(flines) == "toID"] <- "toid"
  flines$toid   <- as.integer(flines$toid)

  # lwgeom::st_linesubstring (inside split_flowlines) requires LINESTRING, so
  # coerce MULTILINESTRING reconciled reaches. Two steps, both needed to keep
  # downstream connectivity:
  #   1. line_merge FIRST: many reaches are stored as multiple raw parts that
  #      merge into a single LINESTRING; merging keeps the whole reach (a raw
  #      keep-longest would drop a connecting part and open a flowline_toid gap).
  #   2. For reaches still multi-part after merge (genuinely disconnected),
  #      keep the component that TOUCHES the downstream reach (the connector) so
  #      topology is preserved; fall back to the longest when none connects.
  is_multi <- sf::st_geometry_type(flines) == "MULTILINESTRING"
  if (any(is_multi)) {
    message(sprintf("  Coercing %d MULTILINESTRING reconciled flowlines to LINESTRING (merge; connector-preserving)",
      sum(is_multi)))
    target_crs <- sf::st_crs(flines)
    geoms      <- sf::st_geometry(flines)
    mi_idx     <- which(is_multi)
    geom_list  <- as.list(geoms)
    fl_ix      <- stats::setNames(seq_len(nrow(flines)), as.character(flines$ID))
    di_all     <- unname(fl_ix[as.character(flines$toid)])   # downstream reach row
    chunk_sz <- 25000L
    for (start in seq(1L, length(mi_idx), by = chunk_sz)) {
      idx    <- mi_idx[start:min(start + chunk_sz - 1L, length(mi_idx))]
      merged <- sf::st_line_merge(sf::st_sfc(geoms[idx], crs = target_crs))
      for (k in seq_along(idx)) {
        g <- merged[[k]]
        if (inherits(g, "LINESTRING")) { geom_list[[idx[k]]] <- g; next }
        parts <- lapply(g, sf::st_linestring)
        psfc  <- sf::st_sfc(parts, crs = target_crs)
        lens  <- suppressWarnings(as.numeric(geos::geos_length(geos::as_geos_geometry(psfc))))
        pick  <- which.max(lens)
        dj    <- di_all[idx[k]]
        if (!is.na(dj)) {
          dd   <- as.numeric(sf::st_distance(psfc, geoms[dj]))
          near <- which(dd <= 30)
          if (length(near)) pick <- near[which.max(lens[near])]  # longest connecting component
        }
        geom_list[[idx[k]]] <- parts[[pick]]
      }
      rm(merged)
    }
    sf::st_geometry(flines) <- sf::st_sfc(geom_list, crs = target_crs)
  }

  # Drop zero-length flowlines (degenerate source geometries, e.g. AK TDX).
  # lwgeom::st_endpoint segfaults on start==end linestrings.
  zero_len <- as.numeric(sf::st_length(flines)) == 0
  if (any(zero_len, na.rm = TRUE)) {
    message(sprintf("add_flowlines: dropped %d zero-length flowline(s) (degenerate geometry)",
      sum(zero_len, na.rm = TRUE)))
    flines <- flines[!zero_len, ]
  }

  .t_split <- Sys.time()
  segs <- split_flowlines(flines, max_length = split_flines_meters)
  message(sprintf(
    "  [add_flowlines timing] split_flowlines: %.1fs (%d flowline(s) -> %d sub-segment(s))",
    as.numeric(difftime(Sys.time(), .t_split, units = "secs")),
    nrow(flines), nrow(segs)))

  # ---- 4. Assign parent flowline ID -----------------------------------------
  # split_flowlines sets the `flowpath_id` column to the ID value for unsplit
  # rows (e.g. "5") and to "5.1", "5.2" for split children.
  # Strip any decimal suffix to recover the original integer reconciled ID.
  segs$orig_fl_id <- as.integer(sub("\\.[0-9]+$", "", as.character(segs$flowpath_id)))

  segs$flowpath_id   <- fl_to_fp[as.character(segs$orig_fl_id)]
  segs$flowpath_toid <- fl_to_fptoid[as.character(segs$orig_fl_id)]
  segs$divide_id     <- fl_to_div[as.character(segs$orig_fl_id)]

  # ---- 4b. Fallback: absorbed reconciled IDs via toID chain walk ---------------
  # reconciled.gpkg may contain features that were merged (absorbed) into a
  # larger representative unit in rec_fp_agg during Stage-1 post-processing.
  # Those absorbed IDs are absent from network$flowline_id -> fl_to_fp lookup,
  # leaving their split sub-segments with flowpath_id = NA.
  # Fix: for each unmatched reconciled ID, walk its toID chain until reaching
  # an ID that IS in the lookup, that ID is the rec_fp_agg representative
  # that absorbed it, so its assigned fp is the correct one.
  na_fp_mask <- is.na(segs$flowpath_id)
  if (any(na_fp_mask)) {
    unmatched_rec_ids <- unique(as.character(segs$orig_fl_id[na_fp_mask]))
    rec_toid_lu       <- stats::setNames(as.character(flines$toid), as.character(flines$ID))

    walk_res <- lapply(unmatched_rec_ids, function(rid) {
      cur  <- rid
      seen <- character(0L)
      for (hop in seq_len(500L)) {
        if (is.na(cur) || cur == "0" || cur %in% seen)
          return(list(fp = NA_character_, div = NA_character_, fptoid = NA_character_))
        seen <- c(seen, cur)
        fp <- fl_to_fp[cur]
        if (!is.na(fp))
          return(list(fp = fp, div = fl_to_div[cur], fptoid = fl_to_fptoid[cur]))
        nxt <- rec_toid_lu[cur]
        if (is.null(nxt) || length(nxt) == 0L || is.na(nxt)) break
        cur <- nxt
      }
      list(fp = NA_character_, div = NA_character_, fptoid = NA_character_)
    })
    names(walk_res) <- unmatched_rec_ids

    fix_orig <- as.character(segs$orig_fl_id[na_fp_mask])
    segs$flowpath_id[na_fp_mask]   <- vapply(fix_orig, function(r) walk_res[[r]]$fp,
      character(1L), USE.NAMES = FALSE)
    segs$flowpath_toid[na_fp_mask] <- vapply(fix_orig, function(r) walk_res[[r]]$fptoid,
      character(1L), USE.NAMES = FALSE)
    segs$divide_id[na_fp_mask]     <- vapply(fix_orig, function(r) walk_res[[r]]$div,
      character(1L), USE.NAMES = FALSE)

    n_filled   <- sum(!is.na(segs$flowpath_id[na_fp_mask]))
    n_still_na <- sum(is.na(segs$flowpath_id))
    message(sprintf(
      "add_flowlines: toid-walk filled %d sub-segment(s) for %d absorbed reconciled ID(s); %d remain unmatched",
      n_filled, length(unmatched_rec_ids), n_still_na))
  }

  # Propagate levelpathid / hydroseq to split children. split_flowlines may
  # leave these NA on the .N suffix pieces; look up from the reconciled
  # parent flowline by orig_fl_id. Downstream navigation (main-stem picking
  # in extend_flowpaths_to_mainstem) relies on levelpathid: if split pieces
  # are lp=NA, every confluence falls back to drainage-area tie-breaking and
  # can follow a tributary instead of the NHD-canonical mainstem.
  flines_df <- sf::st_drop_geometry(flines)
  for (col in c("levelpathid", "orig_levelpathID", "hydroseq")) {
    if (col %in% names(flines_df)) {
      lp_lu <- stats::setNames(flines_df[[col]], as.character(flines_df$ID))
      if (!col %in% names(segs)) segs[[col]] <- NA
      na_mask <- is.na(segs[[col]])
      if (any(na_mask))
        segs[[col]][na_mask] <- lp_lu[as.character(segs$orig_fl_id[na_mask])]
    }
  }

  # ---- 5. Unique sub-segment IDs and topology ---------------------------------
  # id   = seg_id, a sequential integer unique per sub-segment row.  This is
  #        used by extend_flowpaths_to_mainstem and other steps that walk the
  #        sub-flowline chain and require a unique ID per row.
  # toid = next seg_id in the chain (seg_id-space), derived from the reconciled
  #        topology.  For the last sub of reconciled R, toid = first sub of
  #        R.toid; for internal subs, toid = next sub of same R.
  #
  # rec_id   = orig_fl_id, the parent reconciled integer ID (may repeat for
  #            split sub-segments).  This is the externally-visible flowline_id
  #            in the ngen network table, keeping the lineage chain
  #            S0-COMID -> rec_id -> flowpath_id traceable at every stage.
  # rec_toid = next reconciled ID downstream (rec_id-space), for topology
  #            expressed at the reconciled-flowline granularity in the network.
  segs$seg_id <- seq_len(nrow(segs))

  rec_toid_map <- stats::setNames(as.integer(flines$toid), as.character(flines$ID))

  # rec_id / rec_toid columns (rec_id space, for network table)
  segs$rec_id   <- segs$orig_fl_id
  segs$rec_toid <- {
    raw <- rec_toid_map[as.character(segs$orig_fl_id)]
    as.integer(ifelse(is.na(raw), 0L, raw))
  }

  # seg_id topology (for sub-flowline walk algorithms)
  segs$toid <- NA_integer_

  .t_topo <- Sys.time()
  rec_to_subs   <- split(segs$seg_id, segs$orig_fl_id)
  seg_geom      <- sf::st_geometry(segs)

  # Order each reconciled reach's sub-segments in true flow order (inlet ->
  # outlet). split() returns them in split (line-fraction .1..N) order, which
  # equals flow order ONLY when the reconciled line runs upstream->downstream;
  # some lines are stored reversed. Anchor the order on the reach's own
  # downstream (rec_toid) connection: the chain-end sub touching the downstream
  # reach is the OUTLET (last). If instead the FIRST split sub touches the
  # downstream reach, the line is reversed: flip the order. Reaches that are
  # terminal, or whose downstream reach is not adjacent (a genuinely wrong
  # toid), keep split order. This fixes both the internal sub chain and the
  # boundary wiring below at once, so `rec_first_sub[downstream]` is the true
  # inlet sub the upstream outlet should connect to.
  multi <- which(lengths(rec_to_subs) > 1L)
  if (length(multi)) {
    dn_of    <- rec_toid_map[names(rec_to_subs)]
    for (i in multi) {
      dn <- dn_of[i]
      if (is.na(dn) || dn == 0L) next
      dn_subs <- rec_to_subs[[as.character(dn)]]
      if (is.null(dn_subs)) next
      x  <- rec_to_subs[[i]]
      df <- min(as.numeric(sf::st_distance(seg_geom[x[1L]],        seg_geom[dn_subs])))
      dl <- min(as.numeric(sf::st_distance(seg_geom[x[length(x)]], seg_geom[dn_subs])))
      # Reverse ONLY on confident evidence of reversal: the FIRST split sub
      # touches the downstream reach (it is really the outlet) AND the LAST does
      # not. When the downstream reach is not adjacent to either end (a wrong
      # toid), both are far and we must NOT guess: leaving split order avoids
      # spuriously reversing correctly-ordered reaches.
      if (is.finite(df) && is.finite(dl) && df <= 30 && dl > 30) rec_to_subs[[i]] <- rev(x)
    }
  }
  rec_first_sub <- vapply(rec_to_subs, `[`, integer(1L), 1L)
  rec_last_sub  <- vapply(rec_to_subs, function(x) x[length(x)], integer(1L))

  # Fully vectorized sub-segment topology. `segs$seg_id <- seq_len(nrow(segs))`
  # above, so a seg_id value IS its own row index: assign into segs$toid by
  # seg_id directly. Two edge classes:
  #   * internal: within a reconciled group, each sub -> the next sub.
  #   * boundary: each group's LAST sub -> the FIRST sub of its downstream
  #     reconciled group (rec_toid_map), else 0 (terminal).
  # Every lookup is a single vectorized match(): O(total sub-segments), versus
  # a per-rid loop doing an O(N) named-vector lookup per iteration (-> O(N^2)),
  # which dominates the add_flowlines tail on large VPUs.
  segs$toid <- 0L
  has_multi <- lengths(rec_to_subs) > 1L
  if (any(has_multi)) {
    froms <- unlist(lapply(rec_to_subs[has_multi], function(x) x[-length(x)]), use.names = FALSE)
    tos   <- unlist(lapply(rec_to_subs[has_multi], function(x) x[-1L]),        use.names = FALSE)
    segs$toid[froms] <- tos
  }
  # Boundary: each reach's OUTLET sub (rec_last_sub, now flow-ordered) -> the
  # downstream reach's INLET sub. Start from the downstream reach's flow-ordered
  # first sub, then, when it has >1 sub, refine to the sub geometrically nearest
  # this reach's outlet, the true confluence sub, robust to any residual
  # ordering imperfection in the downstream reach.
  next_rid  <- rec_toid_map[names(rec_to_subs)]        # group -> downstream rec id
  nxt_first <- rec_first_sub[as.character(next_rid)]   # downstream reach's inlet sub
  dn_nsub   <- lengths(rec_to_subs)[as.character(next_rid)]
  refine    <- which(!is.na(next_rid) & next_rid != 0L & !is.na(dn_nsub) & dn_nsub > 1L)
  if (length(refine)) {
    for (gi in refine) {
      from_sub <- rec_last_sub[gi]
      dn_subs  <- rec_to_subs[[as.character(next_rid[gi])]]
      d        <- as.numeric(sf::st_distance(seg_geom[from_sub], seg_geom[dn_subs]))
      nxt_first[gi] <- dn_subs[which.min(d)]
    }
  }
  nxt_first[is.na(next_rid) | next_rid == 0L | is.na(nxt_first)] <- 0L
  segs$toid[rec_last_sub] <- as.integer(nxt_first)

  # -- Residual boundary repair -----------------------------------------------
  # A handful of reaches carry a reconciled toid that points to a NON-adjacent
  # reach (a wrong pointer, not an ordering issue the flow-order pass can fix),
  # so the wired boundary sub is still far. Last-resort repair: re-point the
  # reach's outlet sub to the nearest sub it actually TOUCHES, rejecting any
  # choice that would create a cycle (which would reverse flow: e.g. an
  # upstream inflow is the only neighbour). Bounded to the few far edges; the
  # DAG-clean 99.98% are never examined. Reaches with no cycle-safe touch (a
  # genuinely missing/dropped downstream, or an isolated dangle) are left as-is.
  bnd_from <- rec_last_sub
  bnd_to   <- segs$toid[bnd_from]
  vb       <- !is.na(bnd_to) & bnd_to != 0L
  if (any(vb)) {
    fdist <- rep(0, length(bnd_from))
    fdist[vb] <- as.numeric(sf::st_distance(
      seg_geom[bnd_from[vb]], seg_geom[bnd_to[vb]], by_element = TRUE))
    far <- which(vb & fdist > 30)
    if (length(far)) {
      dn_now <- segs$toid                       # seg_id-indexed downstream (row == seg_id)
      reaches_seg <- function(start, target) {
        cur <- start
        for (h in seq_len(5000L)) {
          if (is.na(cur) || cur == 0L) return(FALSE)
          if (cur == target) return(TRUE)
          cur <- dn_now[cur]
        }
        FALSE
      }
      nbf  <- sf::st_intersects(sf::st_buffer(seg_geom[bnd_from[far]], 30), seg_geom)
      nfix <- 0L
      for (j in seq_along(far)) {
        s    <- bnd_from[far[j]]
        cand <- setdiff(nbf[[j]], s)
        if (!length(cand)) next
        cand <- cand[order(as.numeric(sf::st_distance(seg_geom[s], seg_geom[cand])))]
        # A CONFLUENCE IS NOT A CHAIN. Two reaches that share a reconciled
        # downstream are siblings joining at the same node, so the nearest touch
        # to an outlet is very often the sibling arriving beside it rather than
        # the continuation below it. Accepting one turns the confluence into a
        # chain: the sibling's flowpath is spliced between two members of this
        # reach's flowpath, giving that flowpath a second outlet and a
        # flowline-grain path that leaves it and comes back. Cycle-checking
        # cannot see this, the result is still acyclic.
        #
        # `rec_toid` is the reference-derived reconciled topology and is
        # authoritative on the distinction, so use it rather than a tighter
        # tolerance: reject a candidate that drains where I drain. A candidate
        # lying IN my own downstream reach is my successor, not a sibling, and
        # stays eligible. Rare (4 in 284k cross-reach edges over VPUs 18/01/09/13,
        # 2 of them crossing a flowpath boundary), and where it fires the reach
        # keeps its far pointer, a geometric gap already reported as INFO by
        # flowline_toid_geometry_consistent, which is the better of the two.
        s_rec_toid <- segs$rec_toid[s]
        for (pk in cand) {
          is_sibling <- !is.na(s_rec_toid) && s_rec_toid != 0L &&
            !is.na(segs$rec_toid[pk]) && segs$rec_toid[pk] == s_rec_toid &&
            !is.na(segs$rec_id[pk]) && segs$rec_id[pk] != s_rec_toid
          if (is_sibling) next
          if (!reaches_seg(pk, s)) {
            segs$toid[s] <- pk; dn_now[s] <- pk; nfix <- nfix + 1L; break
          }
        }
      }
      if (nfix) message(sprintf(
        "  [add_flowlines] residual boundary repair: re-pointed %d far edge(s) to nearest cycle-safe touch",
        nfix))
    }
  }
  message(sprintf(
    "  [add_flowlines timing] sub-segment topology: %.1fs (%d reconciled id(s), %d sub-segment(s))",
    as.numeric(difftime(Sys.time(), .t_topo, units = "secs")),
    length(rec_to_subs), nrow(segs)))

  # -- 5b. Verify sub-flowline topology is a DAG (seg_id level) ---------------
  {
    ids      <- segs$seg_id
    toids    <- segs$toid
    active   <- !is.na(toids) & toids != 0L
    in_deg   <- tabulate(match(toids[active], ids), nbins = length(ids))
    queue    <- which(in_deg == 0L)
    visited  <- 0L
    while (length(queue) > 0L) {
      n       <- queue[1L]
      queue <- queue[-1L]
      visited <- visited + 1L
      dn      <- toids[n]
      if (!is.na(dn) && dn != 0L) {
        dn_idx <- match(dn, ids)
        if (!is.na(dn_idx)) {
          in_deg[dn_idx] <- in_deg[dn_idx] - 1L
          if (in_deg[dn_idx] == 0L) queue <- c(queue, dn_idx)
        }
      }
    }
    if (visited < nrow(segs)) {
      cycle_ids <- ids[in_deg > 0L]
      stop("add_flowlines: sub-flowline topology contains ", length(cycle_ids),
        " cycle node(s) -- fix upstream in refactor/reconcile, not here.",
        " Example seg_id(s): ", paste(utils::head(cycle_ids, 5), collapse = ", "))
    }
  }

  # ---- 7. FAC at downstream endpoint -----------------------------------------
  # Use terra::crs() (WKT string) rather than sf::st_crs(fac_r): unnamed
  # rasters return NA from sf::st_crs(), leaving transformed points CRS-less
  # and causing terra::extract to fail on the CRS mismatch check.
  fac_wkt   <- terra::crs(fac_r)
  segs_fac  <- sf::st_transform(segs, fac_wkt)
  endpoints <- lwgeom::st_endpoint(sf::st_geometry(segs_fac))
  # Extract FAC by UNIQUE raster cell, in cell (row-major) order, instead of one
  # extract over all endpoints in network order. The endpoints are spatially
  # scattered, so a naive extract thrashes the VRT block cache (re-reading tiles
  # per point), which is I/O-bound on a continental VRT. cellFromXY is pure
  # arithmetic (no I/O): dedup to unique cells (adjacent sub-segments/confluences
  # share cells) and sort so consecutive reads hit the same/adjacent tiles.
  # Result is byte-identical (terra::extract at a point reads the value of the
  # cell containing it = the cell center we read).
  ep_xy     <- sf::st_coordinates(endpoints)
  cells     <- terra::cellFromXY(fac_r, ep_xy)
  fac_vals  <- rep(NA_real_, length(cells))
  uc        <- sort(unique(cells[!is.na(cells)]))
  if (length(uc) > 0L) {
    uc_df   <- terra::extract(fac_r, terra::xyFromCell(fac_r, uc))
    uc_vals <- uc_df[[ncol(uc_df)]]
    fac_vals <- uc_vals[match(cells, uc)]
  }

  fac_vals <- .fac_drop_sentinels(fac_vals)

  # ---- 7b. Fallback FAC setup -----------------------------------------------
  fac_source   <- rep("primary", length(fac_vals))
  fb_fac_path  <- if (is.character(fallback_fac)) fallback_fac else NULL
  primary_path <- if (is.character(fac)) fac else NULL
  use_fallback <- !is.null(fb_fac_path) && nzchar(fb_fac_path) &&
    !identical(normalizePath(if (is.null(primary_path)) "" else primary_path, mustWork = FALSE),
      normalizePath(fb_fac_path, mustWork = FALSE))

  # Buffer fills against VRT-backed rasters (fun=max over polygon buffers) are
  # catastrophically slow. With only a handful of primary NAs per VPU, skip the
  # buffer fill, a direct point extract suffices.
  still_na <- which(is.na(fac_vals))
  use_fallback <- use_fallback && length(still_na) > 0L

  if (use_fallback) {
    fb_r_full <- tryCatch(terra::rast(fb_fac_path), error = function(e) NULL)
    if (is.null(fb_r_full))
      # A configured-but-unreadable fallback is a broken configuration, not an
      # absent option. Degrading quietly here produced a build that skipped the
      # fallback at every call site, passed every invariant, and silently
      # differed from the product it was meant to reproduce.
      stop("add_flowlines: fallback FAC is configured but cannot be read: ",
           fb_fac_path,
           "\n  pass fallback_fac_path = NULL to build without a fallback.",
           call. = FALSE)
  }

  if (use_fallback) {
    fb_wkt    <- terra::crs(fb_r_full)
    ep_fb     <- sf::st_transform(
      sf::st_sfc(endpoints[still_na], crs = fac_wkt),
      fb_wkt)

    # Extract directly from the VRT. terra::extract on point geometries reads
    # only the raster blocks under the points (windowed I/O), so there is no
    # full-VRT tile-scan to avoid. Pre-cropping to the NA-endpoint bbox is
    # catastrophic when many endpoints take the fallback (the bbox can span
    # billions of cells), whereas the windowed extract returns identical values
    # quickly and is just as fast for the typical handful of NA points per VPU.
    fb_vals  <- tryCatch(
      terra::extract(fb_r_full, terra::vect(sf::st_as_sf(ep_fb)))[[2]],
      error = function(e) {
        .hf_log_warn("add_flowlines: fallback FAC extract failed (missing tiles?): ",
          conditionMessage(e))
        rep(NA_real_, length(still_na))
      })

    # Name the grid the value actually came from. This was hard-coded to
    # "gedtm30" for every fallback, so a CONUS build filling from HydroSHEDS
    # shipped rows crediting a product it never read.
    fb_label <- .fac_source_label(fb_fac_path)

    filled_fb <- !is.na(fb_vals)
    fac_vals[still_na[filled_fb]]    <- fb_vals[filled_fb]
    fac_source[still_na[filled_fb]]  <- fb_label

    n_fb <- sum(filled_fb)
    if (n_fb > 0L)
      message(sprintf("add_flowlines: filled %d NA FAC value(s) from the %s fallback",
        n_fb, fb_label))
  }

  segs$fac_value             <- fac_vals
  segs$fac_source            <- fac_source
  segs$tot_drainage_areasqkm <- fac_vals * fac_cell_sqkm

  # ---- 9. Segment length from geometry ---------------------------------------
  # Skip the round-trip transform when segs is already in a projected
  # (metres) CRS, the reconciled flowlines handed in by reconcile_divides are
  # always in a projected CRS. Transform only for longlat input, where
  # st_length would otherwise return spherical metres.
  if (sf::st_is_longlat(segs)) {
    segs$lengthkm <- as.numeric(
      sf::st_length(sf::st_transform(segs, 5072))) / 1000
  } else {
    segs$lengthkm <- as.numeric(sf::st_length(segs)) / 1000
  }

  # ---- 10. Select final columns and write flowlines --------------------------
  # id       = seg_id  (unique per row; used by topology walk algorithms)
  # toid     = next seg_id in chain (seg_id-space)
  # rec_id   = orig_fl_id (parent reconciled ID; becomes flowline_id in ngen)
  # rec_toid = next reconciled ID downstream (rec_id-space; becomes flowline_toid)
  out <- dplyr::select(segs,
    id    = seg_id,   # unique sub-segment ID (for walk algorithms)
    toid,
    rec_id,
    rec_toid,
    flowpath_id,
    flowpath_toid,
    divide_id,
    lengthkm,
    fac_value,
    fac_source,
    tot_drainage_areasqkm,
    dplyr::any_of(c("levelpathid", "hydroseq", "orig_levelpathID"))
  )

  sf::st_write(out, input_gpkg, layer = "flowlines",
    delete_layer = TRUE, quiet = TRUE)

  # ---- 11. Carry events through ----------------------------------------------
  has_poi  <- !is.null(poi_events)
  has_junc <- !is.null(junction_events)

  if (has_poi || has_junc) {
    ev_parts <- list()

    # -- POI events: snapped output of refactor_flowlines, carry reconciled_id -----
    if (has_poi) {
      ev_poi <- poi_events
      # reconciled_id IS the integer reconciled flowline ID -> look up directly
      ev_poi$flowpath_id <- fl_to_fp[as.character(ev_poi$reconciled_id)]
      ev_poi$divide_id   <- fl_to_div[as.character(ev_poi$reconciled_id)]
      ev_poi$event_kind  <- "poi"
      ev_parts[["poi"]]  <- ev_poi
    }

    # -- Junction events: raw NHD flowpath_id, map via member_flowpath_id -----
    if (has_junc) {
      # Build NHD flowpath_id -> reconciled integer ID lookup. Single-alloc
      # construction via vectorized strsplit + rep (was O(n2) rbind).
      flines_df <- sf::st_drop_geometry(flowlines)
      splits_nhd <- strsplit(as.character(flines_df$member_flowpath_id),
        ",", fixed = TRUE)
      lens_nhd   <- lengths(splits_nhd)
      nhd_to_rec <- data.frame(
        nhd_id = sub("\\.[0-9]+$", "",
          trimws(unlist(splits_nhd, use.names = FALSE))),
        rec_id = rep(as.integer(flines_df$ID), lens_nhd),
        stringsAsFactors = FALSE)
      nhd_to_rec    <- nhd_to_rec[!duplicated(nhd_to_rec$nhd_id), ]
      nhd_to_rec_id <- stats::setNames(nhd_to_rec$rec_id, nhd_to_rec$nhd_id)

      ev_junc <- junction_events
      ev_junc$flowpath_id <- as.character(ev_junc$flowpath_id)
      rec_ids             <- nhd_to_rec_id[ev_junc$flowpath_id]
      ev_junc$flowpath_id <- fl_to_fp[as.character(rec_ids)]
      ev_junc$divide_id   <- fl_to_div[as.character(rec_ids)]
      ev_junc$event_kind  <- "junction"
      ev_parts[["junction"]] <- ev_junc
    }

    # Combine and write: use dplyr::bind_rows which fills missing cols with NA
    ev_out <- dplyr::bind_rows(ev_parts)

    sf::st_write(sf::st_as_sf(ev_out), input_gpkg, layer = "events",
      delete_layer = TRUE, quiet = TRUE)
  }

  invisible(out)
}
