#' Generate events from hydrolocations (data.table + sf nearest-feature, GEOS)
#'
#' Reads only the necessary layers and uses vectorized nearest-feature queries
#' to snap hydrolocations onto reference flowpaths and emit split events.
#' This is the sole event-generation entry point used by the refactor phase.
#'
#' @param hydrolocations_gpkg Path to the GeoPackage of hydrolocations (POIs) to
#'   snap onto flowpaths.
#' @param flowpaths_gpkg Optional path to the reference flowpaths GeoPackage;
#'   defaults to `hydrolocations_gpkg` when `NULL`.
#' @param hydroloc_layer Hydrolocations layer name; auto-detected when `NULL`.
#' @param flowpath_layer Reference flowpaths layer name. Default
#'   `"reference_flowpaths"`.
#' @param out_gpkg Optional output path; when `NULL` the events `sf` is returned
#'   instead of written.
#' @param force_lr Logical; request left/right bank attribution where available.
#'   Default `TRUE`.
#' @param tolerance_m numeric (meters). Maximum distance from a flowpath for a
#'   hydrolocation to match; farther candidates are dropped. Default 50.
#' @param collapse_threshold_m numeric (meters). An event whose split would leave
#'   a piece shorter than this collapses onto the nearest flowpath endpoint.
#'   Default 250.
#' @param min_area_m numeric. Minimum drainage-area filter for events; `0`
#'   disables it (default).
#' @param snap_threshold_m numeric (meters). Distance within which an event snaps
#'   to a flowpath outlet/inlet node. Default 250.
#' @return An `sf` of events; when `out_gpkg` is supplied it is written there and
#'   returned invisibly.
#' @keywords internal
rf_generate_events_from_hydrolocations_fast <- function(
    hydrolocations_gpkg,
    flowpaths_gpkg = NULL,
    hydroloc_layer = NULL,
    flowpath_layer = "reference_flowpaths",
    out_gpkg = NULL,
    force_lr = TRUE,
    tolerance_m = 50,
    collapse_threshold_m = 250,
    min_area_m = 0,
    snap_threshold_m = 250
    ) {
  if (!requireNamespace("sf", quietly = TRUE)) stop("sf required")
  if (!suppressWarnings(requireNamespace("data.table", quietly = TRUE))) stop("data.table required for fast generation")

  if (is.null(flowpaths_gpkg)) flowpaths_gpkg <- hydrolocations_gpkg

  # Read hydrolocations (filter by layer if provided)
  if (is.null(hydroloc_layer)) {
    events <- tryCatch(sf::st_read(hydrolocations_gpkg, quiet = TRUE), error = function(e) NULL)
  } else {
    events <- tryCatch(sf::st_read(hydrolocations_gpkg, hydroloc_layer, quiet = TRUE), error = function(e) NULL)
  }
  if (is.null(events) || nrow(events) == 0) {
    if (!is.null(out_gpkg) && !is.null(events)) sf::st_write(events, out_gpkg, layer = "events", delete_layer = TRUE, quiet = TRUE)
    return(events)
  }

  # Optionally filter by vpuid if events have that column (caller may pass vpuid filter separately)
  # Keep geometry and minimal attributes
  events <- sf::st_as_sf(events)

  # Read flowpaths but keep only id/toid/length + geometry
  available_layers <- sf::st_layers(flowpaths_gpkg)$name
  if (!flowpath_layer %in% available_layers) {
    fp_candidates <- c("flowpaths", "reference_flowpaths", "reference_flowlines", "flowlines", "reference_flowpath")
    match_layer <- intersect(tolower(available_layers), fp_candidates)
    if (length(match_layer) == 0) stop("Could not find flowpath layer in ", flowpaths_gpkg)
    flowpath_layer <- available_layers[tolower(available_layers) == match_layer[1]][1]
  }
  fp <- sf::st_read(flowpaths_gpkg, layer = flowpath_layer, quiet = TRUE)
  # ensure crs matches events for spatial ops
  fp <- sf::st_transform(fp, sf::st_crs(events))

  # Normalize columns
  if (!"flowpath_id" %in% names(fp)) stop("flowpath_id required in flowpath layer")
  if (!"flowpath_toid" %in% names(fp)) stop("flowpath_toid required in flowpath layer")
  fp$flowpath_id <- as.character(fp$flowpath_id)
  fp$flowpath_toid <- as.character(fp$flowpath_toid)
  if (!"lengthkm" %in% tolower(names(fp))) {
    fp$length_m <- as.numeric(sf::st_length(sf::st_geometry(fp)))
  } else {
    len_col <- names(fp)[tolower(names(fp)) == "lengthkm"][1]
    fp$length_m <- as.numeric(fp[[len_col]]) * 1000
  }

  # Build endpoints: inlet (first coord) and outlet (last coord)
  coords_list <- sf::st_coordinates(sf::st_geometry(fp))
  # Use sf functions to extract endpoints robustly
  get_endpoints <- function(g) {
    cr <- sf::st_coordinates(g)
    if (nrow(cr) < 1) return(list(inlet = NULL, outlet = NULL))
    list(inlet = sf::st_point(cr[1, c("X", "Y")]), outlet = sf::st_point(cr[nrow(cr), c("X", "Y")]))
  }
  inlet_pts <- vector("list", nrow(fp))
  outlet_pts <- vector("list", nrow(fp))
  for (i in seq_len(nrow(fp))) {
    eps <- get_endpoints(sf::st_geometry(fp)[[i]])
    inlet_pts[[i]] <- eps$inlet
    outlet_pts[[i]] <- eps$outlet
  }

  # Build endpoint sfc combining type and index
  inlet_sfc <- sf::st_sfc(inlet_pts, crs = sf::st_crs(fp))
  outlet_sfc <- sf::st_sfc(outlet_pts, crs = sf::st_crs(fp))

  # Build data.tables for fast attribute ops
  dt_fp <- data.table::as.data.table(sf::st_drop_geometry(fp))
  dt_fp[, fp_idx := .I]
  dt_fp[, flowpath_id := as.character(flowpath_id)]
  dt_fp[, flowpath_toid := as.character(flowpath_toid)]
  dt_fp[, length_m := as.numeric(length_m)]

  # Build a fast lookup: flowpath_id -> row index in fp
  fp_row_by_id <- stats::setNames(seq_len(nrow(fp)), as.character(dt_fp$flowpath_id))

  # Precompute upstream list mapping: flowpath_id -> fp row indices of upstreams
  upstream_map <- split(dt_fp$fp_idx, dt_fp$flowpath_toid)

  # Initialize fields
  events$FromMeas_m <- ifelse("FromMeas" %in% names(events), as.numeric(events$FromMeas), NA_real_)
  events$was_snapped <- FALSE

  # Preserve the hydrolocation's original flowpath_id (reference comid) before any
  # spatial snap overwrites it. .snap_events_to_reconciled() uses this for the
  # member-list lookup that correctly maps reference comids to reconciled fp groups.
  if ("flowpath_id" %in% names(events))
    events$orig_flowpath_id <- as.character(events$flowpath_id)

  # Lineage: why a hydrolocation's reach assignment changed.
  #
  # Two blocks below re-point a hydrolocation from the reach it is the INLET of to
  # the upstream reach it is the OUTLET of: the node-type endpoint snap, and the
  # near-endpoint cleanup that keeps a split from leaving a sub-threshold sliver.
  # Both are correct and both change flowpath_id, so recording the reason keeps the
  # delivered layer able to explain a mismatch against the supplied id.
  events$hl_reach_change_reason <- NA_character_

  # Initial flowpath assignment:
  # The reference fabric already carries flowpath_id on each hydrolocation (1:1
  # from linear referencing during fabric build). Use that directly. Only fall
  # back to st_nearest_feature for hydrolocations that lack a pre-assigned id.
  {
    has_id <- "flowpath_id" %in% names(events) &&
      !all(is.na(events$flowpath_id) | events$flowpath_id == "")
    if (has_id) {
      events$flowpath_id <- as.character(events$flowpath_id)
      # Validate: any pre-assigned id not in the flowpath layer?
      bad <- !is.na(events$flowpath_id) &
        !(events$flowpath_id %in% dt_fp$flowpath_id)
      if (any(bad)) {
        message(sprintf(
          "rf_generate_events: %d hydrolocation(s) have flowpath_id absent from network -- snapping spatially",
          sum(bad)))
        fallback_idx <- sf::st_nearest_feature(
          sf::st_transform(sf::st_geometry(events)[bad], sf::st_crs(fp)),
          sf::st_geometry(fp))
        events$flowpath_id[bad] <- dt_fp$flowpath_id[fallback_idx]
      }
      # NA/empty ids also get spatial fallback
      missing <- is.na(events$flowpath_id) | !nzchar(events$flowpath_id)
      if (any(missing)) {
        fallback_idx <- sf::st_nearest_feature(
          sf::st_transform(sf::st_geometry(events)[missing], sf::st_crs(fp)),
          sf::st_geometry(fp))
        events$flowpath_id[missing] <- dt_fp$flowpath_id[fallback_idx]
      }
    } else {
      nearest_idx <- sf::st_nearest_feature(
        sf::st_transform(events, sf::st_crs(fp)),
        sf::st_geometry(fp))
      events$flowpath_id <- dt_fp$flowpath_id[nearest_idx]
    }
    events$flowpath_id <- as.character(events$flowpath_id)
  }

  # Assignment-aware endpoint snapping:
  # For each event, check the inlet and outlet of its **assigned** flowpath only.
  # This prevents teleportation to a different flowpath whose endpoint happens to
  # be the globally nearest point.
  if ("flowpath_id" %in% names(events)) {
    events$flowpath_id <- as.character(events$flowpath_id)

    # Precompute coordinate matrices once: per-event XY, per-fp inlet XY,
    # per-fp outlet XY. Replaces 2xN per-row sf::st_distance calls with a
    # single vectorized rowSums.
    ev_xy <- sf::st_coordinates(events)[, c("X", "Y"), drop = FALSE]
    in_xy <- do.call(rbind,
      lapply(inlet_pts, function(p) {
        if (is.null(p) || length(p) == 0L) c(NA_real_, NA_real_) else
          as.numeric(p)[1:2]
      }))
    out_xy <- do.call(rbind,
      lapply(outlet_pts, function(p) {
        if (is.null(p) || length(p) == 0L) c(NA_real_, NA_real_) else
          as.numeric(p)[1:2]
      }))
    fp_is_all <- fp_row_by_id[events$flowpath_id]
    d_in_all  <- rep(NA_real_, nrow(events))
    d_out_all <- rep(NA_real_, nrow(events))
    v <- !is.na(fp_is_all)
    if (any(v)) {
      d_in_all[v]  <- sqrt(rowSums((ev_xy[v, , drop = FALSE] -
        in_xy[fp_is_all[v], , drop = FALSE])^2))
      d_out_all[v] <- sqrt(rowSums((ev_xy[v, , drop = FALSE] -
        out_xy[fp_is_all[v], , drop = FALSE])^2))
    }

    # hl_types that must live at network nodes (inlet/outlet): eligible for
    # endpoint-snap. Gages and similar point observations snap to their nearest
    # position on the line instead (handled in the tolerance_m block below).
    node_types <- c("huc12", "wbin", "wbout", "nwm_waterbody", "coastal")
    ev_hl_type <- if ("hl_type" %in% names(events)) as.character(events$hl_type)
    else rep(NA_character_, nrow(events))

    for (ri in seq_len(nrow(events))) {
      fid <- events$flowpath_id[ri]
      if (is.na(fid) || fid == "") next
      fp_i <- fp_row_by_id[fid]
      if (is.na(fp_i)) next

      # Skip endpoint-snap for non-node types (e.g. gages); let them fall
      # through to the tolerance_m nearest-point snap below.
      hl_t <- ev_hl_type[ri]
      if (!is.na(hl_t) && !hl_t %in% node_types) next

      d_inlet  <- d_in_all[ri]
      d_outlet <- d_out_all[ri]

      near_outlet <- !is.na(d_outlet) && d_outlet <= snap_threshold_m
      near_inlet  <- !is.na(d_inlet)  && d_inlet  <= snap_threshold_m

      if (near_outlet && (!near_inlet || d_outlet <= d_inlet)) {
        events$FromMeas_m[ri] <- dt_fp$length_m[fp_i]
        sf::st_geometry(events)[[ri]] <- outlet_pts[[fp_i]]
        events$was_snapped[ri] <- TRUE
        events$hl_reach_change_reason[ri] <- "node_at_outlet_moved_to_own_outlet"
      } else if (near_inlet) {
        # Redirect to the outlet of the upstream flowpath continuing the same
        # mainstem. Where none does, take the largest-drainage upstream rather
        # than whichever happens to be first: at a confluence, list order can put
        # a small tributary in place of the channel the POI belongs to.
        upstream_idxs <- upstream_map[[fid]]
        upstream_idxs <- upstream_idxs[!is.na(upstream_idxs) & dt_fp$flowpath_id[upstream_idxs] != fid]
        chosen <- .pick_upstream(upstream_idxs, fp_i, dt_fp)
        if (!is.na(chosen)) {
          events$flowpath_id[ri] <- dt_fp$flowpath_id[chosen]
          events$FromMeas_m[ri]  <- dt_fp$length_m[chosen]
          sf::st_geometry(events)[[ri]] <- outlet_pts[[chosen]]
          events$was_snapped[ri] <- TRUE
          events$hl_reach_change_reason[ri] <- "node_at_inlet_moved_to_upstream_outlet"
        } else {
          # Headwater: no upstream reach exists, so the POI keeps its own.
          sf::st_geometry(events)[[ri]] <- inlet_pts[[fp_i]]
          events$FromMeas_m[ri] <- 0
          events$was_snapped[ri] <- TRUE
          events$hl_reach_change_reason[ri] <- "at_inlet_of_headwater_kept_own_reach"
        }
      }
    }
  }

  # Snap remaining events (not already endpoint-snapped) to nearest point on their
  # assigned flowline, if they are within tolerance_m. This is the core geometric
  # move that places events onto the line rather than leaving them hanging off it.
  if (is.numeric(tolerance_m) && tolerance_m > 0) {
    unsnapped <- which(!events$was_snapped)
    if (length(unsnapped) > 0 && "flowpath_id" %in% names(events)) {
      fp_geoms <- sf::st_geometry(fp)
      # Batch the (event -> its assigned fp-line) distance and nearest-point
      # calls across all unsnapped events with single pairwise calls plus a
      # lookup, keeping the per-event measure assignment inline.
      fp_is_uns <- fp_row_by_id[as.character(events$flowpath_id[unsnapped])]
      valid_mask <- !is.na(fp_is_uns)
      if (any(valid_mask)) {
        u         <- unsnapped[valid_mask]
        fp_v      <- fp_is_uns[valid_mask]
        ev_pts    <- sf::st_geometry(events)[u]
        line_geom <- fp_geoms[fp_v]
        d_all     <- suppressWarnings(as.numeric(
          sf::st_distance(ev_pts, line_geom, by_element = TRUE)))
        within    <- is.finite(d_all) & d_all <= tolerance_m
        n_far <- sum(!within, na.rm = TRUE)
        if (n_far > 0L)
          message(sprintf(paste0("rf_generate_events: %d hydrolocation(s) > %g m ",
            "from their flowpath (kept on nearest, not precisely snapped) -- ",
            "raise tolerance_m if these should snap"), n_far, tolerance_m))

        snap_coords <- matrix(NA_real_, nrow = length(u), ncol = 2L,
          dimnames = list(NULL, c("X", "Y")))
        if (any(within)) {
          # A multipart flowpath can make GEOS return MULTIPOINT for one
          # pair and LINESTRING for another.  Combining those results creates
          # sfc_GEOMETRY, for which st_coordinates() is not implemented.
          # The event count is small, so extract the second endpoint per pair
          # while retaining the vectorized distance calculation above.
          within_rows <- which(within)
          snap_pts <- do.call(rbind, lapply(within_rows, function(k) {
            np <- sf::st_nearest_points(
              ev_pts[k], line_geom[k], pairwise = TRUE)
            npc <- sf::st_coordinates(np)
            npc[nrow(npc), c("X", "Y"), drop = FALSE]
          }))
          snap_coords[within, ] <- snap_pts
        }

        # Map unsnapped->(within, snap_coords) for inline per-event use below
        d_cache     <- stats::setNames(d_all,           as.character(u))
        within_ix   <- stats::setNames(within,          as.character(u))
        snap_x_all  <- stats::setNames(snap_coords[, "X"], as.character(u))
        snap_y_all  <- stats::setNames(snap_coords[, "Y"], as.character(u))
      } else {
        d_cache <- within_ix <- snap_x_all <- snap_y_all <- stats::setNames(
          numeric(0), character(0))
      }

      for (ri in unsnapped) {
        fid <- as.character(events$flowpath_id[ri])
        if (is.na(fid) || fid == "") next
        fp_i <- fp_row_by_id[fid]
        if (is.na(fp_i)) next
        key <- as.character(ri)
        if (is.na(d_cache[key]) || !within_ix[key]) next
        line <- fp_geoms[[fp_i]]
        snapped_coord <- c(X = unname(snap_x_all[key]),
          Y = unname(snap_y_all[key]))
        snapped_pt <- sf::st_point(snapped_coord)
        sf::st_geometry(events)[[ri]] <- snapped_pt
        # Recompute FromMeas_m as cumulative distance from inlet to snapped point
        coords <- sf::st_coordinates(line)
        seg_lens <- sqrt(diff(coords[, "X"])^2 + diff(coords[, "Y"])^2)
        dists_to_vertices <- sqrt((coords[, "X"] - snapped_coord["X"])^2 +
          (coords[, "Y"] - snapped_coord["Y"])^2)
        nearest_v <- which.min(dists_to_vertices)
        if (nearest_v == 1L) {
          events$FromMeas_m[ri] <- 0
        } else {
          events$FromMeas_m[ri] <- sum(seg_lens[seq_len(nearest_v - 1L)])
        }
        events$was_snapped[ri] <- TRUE
      }
    }
  }

  # -- Near-endpoint POI cleanup ---------------------------------------------
  # A POI at position P on a flowpath of length L forces split pieces [0,P]
  # and [P,L]. If either piece < collapse_threshold_m it becomes an
  # uncollapseable sliver. Fix universally (all poi types) after tolerance snap:
  #   near outlet (L-P < threshold) -> snap to this flowpath's outlet
  #   near inlet  (P   < threshold) -> reassign to upstream outlet, same levelpath
  # This mirrors the existing node_type endpoint-snap logic and respects the
  # invariant that a POI at an inlet IS the outlet of its upstream segment.
  if (is.numeric(collapse_threshold_m) && collapse_threshold_m > 0 &&
    "FromMeas_m" %in% names(events) && "flowpath_id" %in% names(events)) {
    ev_fp_ne <- as.character(events$flowpath_id)
    ev_mea_ne <- as.numeric(events$FromMeas_m)
    n_near_out <- 0L
    n_near_in  <- 0L
    for (ri in seq_len(nrow(events))) {
      fid  <- ev_fp_ne[ri]
      if (is.na(fid) || !nzchar(fid)) next
      fp_i <- fp_row_by_id[fid]
      if (is.na(fp_i)) next
      fl <- dt_fp$length_m[fp_i]
      m  <- ev_mea_ne[ri]
      if (!is.finite(m) || !is.finite(fl) || fl <= 0) next
      down_piece <- fl - m
      up_piece   <- m

      # A POI already AT the outlet needs no action, but it must also not fall
      # through to the inlet branch below. On a reach shorter than the threshold
      # its up_piece is the whole reach length, which is under threshold too, so
      # the inlet branch fired and walked a correctly-placed POI up to its
      # upstream neighbour: 71 on VPU 01, 61 of them wbout, up to 239 m the wrong
      # way, and recorded as an inlet move.
      if (down_piece == 0) next

      # `down_piece <= up_piece`: the outlet branch is evaluated first, so without
      # it a POI nearer the INLET still snapped to the outlet whenever its
      # downstream piece happened to be under threshold. Gage 01022840 sat 157 m
      # from its inlet and 222 m from its outlet on a 394 m reach and was sent to
      # the outlet, the farther end, and the wrong drainage area.
      if (m > fl || (down_piece > 0 && down_piece < collapse_threshold_m &&
                     down_piece <= up_piece)) {
        # Past outlet (float precision) or near outlet: snap to outlet
        ev_mea_ne[ri]                 <- fl
        events$FromMeas_m[ri]         <- fl
        sf::st_geometry(events)[[ri]] <- outlet_pts[[fp_i]]
        events$was_snapped[ri]        <- TRUE
        events$hl_reach_change_reason[ri] <- "within_split_floor_of_outlet_moved_to_own_outlet"
        n_near_out <- n_near_out + 1L
        next
      }

      # `>= 0`, not `> 0`: a POI exactly at the inlet has up_piece == 0, and that
      # is the case this branch exists for. (The outlet branch above is correctly
      # `> 0`, a POI already at the outlet needs no action.)
      #
      # Both ends can be under threshold on a short reach, so the inlet branch
      # wins only when the inlet is actually nearer. The rule is "move to the
      # nearer end", not "prefer upstream".
      if (up_piece >= 0 && up_piece < collapse_threshold_m && up_piece <= down_piece) {
        # Near inlet: reassign to the outlet of the upstream segment continuing
        # the same mainstem.
        ups <- upstream_map[[fid]]
        ups <- ups[!is.na(ups) & dt_fp$flowpath_id[ups] != fid]
        chosen <- .pick_upstream(ups, fp_i, dt_fp)
        if (!is.na(chosen)) {
          ev_fp_ne[ri]                  <- dt_fp$flowpath_id[chosen]
          ev_mea_ne[ri]                 <- dt_fp$length_m[chosen]
          events$flowpath_id[ri]        <- dt_fp$flowpath_id[chosen]
          events$FromMeas_m[ri]         <- dt_fp$length_m[chosen]
          sf::st_geometry(events)[[ri]] <- outlet_pts[[chosen]]
          events$was_snapped[ri]        <- TRUE
          events$hl_reach_change_reason[ri] <- "within_split_floor_of_inlet_moved_to_upstream_outlet"
          n_near_in <- n_near_in + 1L
        }
        # No upstream (headwater): the POI stays on its own reach, because there
        # is no upstream outlet for it to belong to. Recorded rather than left
        # silent: 33 of VPU 01's POIs take this path and every other departure
        # from the source assignment carries a reason.
        if (is.na(chosen))
          events$hl_reach_change_reason[ri] <- "at_inlet_of_headwater_kept_own_reach"
      }
    }
    if (n_near_out + n_near_in > 0L)
      message(sprintf(
        "rf_generate_events: near-endpoint cleanup: %d snapped to outlet, %d reassigned to upstream outlet",
        n_near_out, n_near_in))
  }

  # -- Collapse co-located POIs within collapse_threshold_m ------------------
  # Multiple POIs snapped near each other (e.g. two USGS gauges ~16 m apart on
  # the same reach) would each force a refactor split, creating pathological
  # sub-cell divide pieces that PRD cannot partition. Collapse such clusters
  # down to ONE representative and aggregate the dropped rows' hl_link /
  # hl_reference values as pipe-delimited alias lists on the keeper so the
  # gauges remain discoverable via the hydrolocation lookup.
  #
  # Keeper priority:
  #   1. An event whose snapped coord is at a JUNCTION (a flowpath inlet
  #      that is the target of >=2 flowpath_toid pointers: NHD confluence).
  #   2. Otherwise, the event on the most DOWNSTREAM flowpath (min
  #      hydroseq): NHD convention: lower hydroseq = further downstream.
  if (is.numeric(collapse_threshold_m) && collapse_threshold_m > 0 &&
    nrow(events) > 1L) {
    ev_crs <- sf::st_crs(events)
    fp_t   <- if (sf::st_crs(fp) != ev_crs) sf::st_transform(fp, ev_crs) else fp

    # Junction points: inlet of any flowpath that is itself the target of
    # >=2 flowpath_toid pointers (2+ upstream flowpaths converge there).
    toid_tbl  <- table(as.character(fp_t$flowpath_toid)[
      !is.na(fp_t$flowpath_toid) & nzchar(fp_t$flowpath_toid) &
        as.character(fp_t$flowpath_toid) != "0"])
    junc_toids <- names(toid_tbl)[toid_tbl > 1L]
    junction_pts <- if (length(junc_toids)) {
      junc_fp <- fp_t[as.character(fp_t$flowpath_id) %in% junc_toids, ]
      if (nrow(junc_fp) > 0L)
        lwgeom::st_startpoint(sf::st_geometry(junc_fp))
      else sf::st_sfc(crs = ev_crs)
    } else sf::st_sfc(crs = ev_crs)

    # Per-event: at_junction flag (<= 5 m from any junction point)
    at_junction <- rep(FALSE, nrow(events))
    if (length(junction_pts) > 0L) {
      dmin_j <- suppressWarnings(as.numeric(sf::st_distance(
        sf::st_geometry(events), junction_pts)))
      dmin_j <- matrix(dmin_j, nrow = nrow(events),
        ncol = length(junction_pts))
      at_junction <- apply(dmin_j, 1L, function(r) any(r < 5, na.rm = TRUE))
    }

    # Hydroseq lookup (for keeper selection when no junction member)
    hyd_by_fp <- if ("hydroseq" %in% names(fp_t))
      stats::setNames(as.numeric(fp_t$hydroseq),
        as.character(fp_t$flowpath_id))
    else stats::setNames(rep(NA_real_, nrow(fp_t)),
      as.character(fp_t$flowpath_id))
    ev_hyd <- unname(hyd_by_fp[as.character(events$flowpath_id)])

    # Cluster events that are co-located ON THE SAME FLOWPATH within
    # collapse_threshold_m along-line distance. Using Euclidean/spatial
    # distance via st_is_within_distance causes chain-collapse across the
    # network: A within 250m of B, B within 250m of C -> all merged even
    # when A and C are hundreds of km apart. Instead, cluster only within
    # the same flowpath using |FromMeas_m_a - FromMeas_m_b| <= threshold.
    ev_fp  <- as.character(events$flowpath_id)
    ev_mea <- as.numeric(events$FromMeas_m)
    cluster_id <- seq_len(nrow(events))
    fp_groups <- split(seq_len(nrow(events)), ev_fp)
    for (idxs in fp_groups) {
      if (length(idxs) < 2L) next
      mea <- ev_mea[idxs]
      # Pairwise collapse within this flowpath (small group, O(n2) is fine)
      parent_local <- seq_along(idxs)
      find_local <- function(x) {
        while (parent_local[x] != x) x <- parent_local[x]
        x
      }
      for (a in seq_along(idxs)) {
        for (b in seq_along(idxs)) {
          if (b <= a) next
          if (!is.na(mea[a]) && !is.na(mea[b]) &&
            abs(mea[a] - mea[b]) <= collapse_threshold_m) {
            ra <- find_local(a)
            rb <- find_local(b)
            if (ra != rb) parent_local[rb] <- ra
          }
        }
      }
      local_cluster <- vapply(seq_along(idxs), find_local, integer(1L))
      cluster_id[idxs] <- idxs[local_cluster]
    }

    # Aggregate clusters: pick keeper + accumulate alias lists
    keep_mask <- rep(TRUE, nrow(events))
    alias_cols <- intersect(
      c("hl_link", "hl_reference", "hl_id", "hl_type", "hl_source", "hl_uri"),
      names(events))
    orig_alias <- lapply(alias_cols, function(cc) as.character(events[[cc]]))
    names(orig_alias) <- alias_cols

    cl_groups <- split(seq_len(nrow(events)), cluster_id)
    n_collapsed <- 0L
    for (g in cl_groups) {
      if (length(g) < 2L) next
      # Prefer junction member(s); break ties by most-downstream hydroseq
      cand <- if (any(at_junction[g])) g[at_junction[g]] else g
      cand_hyd <- ev_hyd[cand]
      keeper <- if (all(is.na(cand_hyd))) cand[1L] else
        cand[which.min(replace(cand_hyd, is.na(cand_hyd), Inf))]
      # Concatenate dropped members' alias values into keeper (pipe-separated,
      # unique, dropping NAs; keep keeper's own value first).
      for (cc in alias_cols) {
        vals <- unique(stats::na.omit(orig_alias[[cc]][g]))
        vals <- vals[nzchar(vals)]
        if (length(vals) > 1L) {
          events[[cc]][keeper] <- paste(vals, collapse = "|")
        }
      }
      keep_mask[setdiff(g, keeper)] <- FALSE
      n_collapsed <- n_collapsed + (length(g) - 1L)
    }
    if (n_collapsed > 0L) {
      events <- events[keep_mask, , drop = FALSE]
      message(sprintf(
        "rf_generate_events: collapsed %d co-located POI(s) (<%g m); kept %d",
        n_collapsed, collapse_threshold_m, nrow(events)))
    }
  }

  # Label original hydrolocation rows as "poi" before junction rows are added
  if (!"event_kind" %in% names(events)) {
    events$event_kind <- "poi"
  } else {
    events$event_kind <- as.character(events$event_kind)
    events$event_kind[is.na(events$event_kind)] <- "poi"
  }
  if (!"is_junction" %in% names(events)) events$is_junction <- FALSE

  # local helpers for endpoints used in junction construction
  get_outlet_point_local <- function(geom) {
    coords <- sf::st_coordinates(geom)
    if (nrow(coords) < 1) return(NULL)
    sf::st_point(coords[nrow(coords), c("X", "Y")])
  }
  # Add junction-generated events similar to original if needed
  junction_counts <- dt_fp[
    !is.na(flowpath_toid) & nzchar(flowpath_toid) & flowpath_toid != "0",
    .(upstream_count = .N),
    by = flowpath_toid
  ][upstream_count > 1L]

  if (nrow(junction_counts) > 0) {
    inflow_idx <- dt_fp[flowpath_toid %in% junction_counts$flowpath_toid, fp_idx]
    if (length(inflow_idx) > 0) {
      inflow_rows <- fp[inflow_idx, , drop = FALSE]
      inflow_rows$upstream_count <- junction_counts$upstream_count[match(
        as.character(inflow_rows$flowpath_toid),
        as.character(junction_counts$flowpath_toid)
      )]
      junction_rows <- events[rep(1, nrow(inflow_rows)), , drop = FALSE]
      junction_rows$flowpath_id <- as.character(inflow_rows$flowpath_id)
      # Synthetic junction rows are topology-derived; clear inherited HL metadata.
      for (nm in intersect(c(
        "poi_id", "hl_id", "hl_link", "hl_type", "hl_reference", "hl_class",
        "hl_source", "hl_uri", "identifier", "FromMeas_raw", "REACH_meas_raw",
        "reachmeas"
      ), names(junction_rows))) {
        junction_rows[[nm]] <- NA
      }
      if ("FromMeas_m" %in% names(junction_rows)) junction_rows$FromMeas_m <- as.numeric(inflow_rows$length_m)
      inflow_geoms <- sf::st_geometry(inflow_rows)
      junction_geoms <- sf::st_geometry(junction_rows)
      for (i in seq_len(nrow(junction_rows))) {
        out_pt <- get_outlet_point_local(inflow_geoms[[i]])
        if (!is.null(out_pt)) junction_geoms[[i]] <- out_pt
      }
      sf::st_geometry(junction_rows) <- sf::st_sfc(junction_geoms, crs = sf::st_crs(events))
      junction_rows$is_junction <- rep(TRUE, nrow(junction_rows))
      junction_rows$event_kind <- rep("junction", nrow(junction_rows))
      if ("was_snapped" %in% names(junction_rows)) junction_rows$was_snapped <- FALSE
      events <- dplyr::bind_rows(events, junction_rows)
    }
  }

  # -- Output invariants -----------------------------------------------------
  # Warn (not stop) so a bad VPU doesn't abort a batch run, but every
  # violation is surfaced before the events reach refactor_flowlines.
  poi_ev  <- events[events$event_kind == "poi", , drop = FALSE]
  junc_ev <- events[events$event_kind == "junction", , drop = FALSE]

  # 1. Every poi has a flowpath_id that exists in the reference network
  poi_fp <- as.character(poi_ev$flowpath_id)
  na_fp  <- is.na(poi_fp) | !nzchar(poi_fp)
  if (any(na_fp))
    .hf_log_warn(sprintf(
      "[events] %d poi event(s) have NA/empty flowpath_id after all snapping",
      sum(na_fp)))
  bad_fp <- !na_fp & !(poi_fp %in% dt_fp$flowpath_id)
  if (any(bad_fp))
    .hf_log_warn(sprintf(
      "[events] %d poi event(s) have flowpath_id absent from VPU network",
      sum(bad_fp)))
  if (sum(na_fp) + sum(bad_fp) == 0L)
    message(sprintf(
      "[events] invariant OK: all %d poi flowpath_id values in network",
      nrow(poi_ev)))

  # 2. FromMeas_m within [0, flowpath_length] for all poi events
  if ("FromMeas_m" %in% names(poi_ev)) {
    mea  <- as.numeric(poi_ev$FromMeas_m)
    fp_m <- dt_fp$length_m[match(poi_fp, dt_fp$flowpath_id)]
    oor  <- !is.na(mea) & !is.na(fp_m) & (mea < -1 | mea > fp_m + 1)
    if (any(oor))
      .hf_log_warn(sprintf(
        "[events] %d poi event(s) have FromMeas_m outside [0, flowpath_length]",
        sum(oor)))
    else
      message("[events] invariant OK: all poi FromMeas_m values in range")
  }

  # 3. Duplicate poi events on the same flowpath (informational, each forces
  #    a split; pairs far apart are valid, near-coincident pairs are not)
  dup_fp <- unique(poi_fp[!na_fp & duplicated(poi_fp[!na_fp])])
  if (length(dup_fp) > 0L) {
    message(sprintf(
      "[events] %d flowpath(s) carry >1 poi event (each forces a split): %s%s",
      length(dup_fp),
      paste(utils::head(dup_fp, 5), collapse = ", "),
      if (length(dup_fp) > 5L) " ..." else ""))
  } else {
    message("[events] invariant OK: no flowpath has >1 poi event")
  }

  # 4. Junction count vs. network confluences
  #    Expected = one junction row per inflow segment at each confluence.
  tc_tbl <- dt_fp[
    !is.na(flowpath_toid) & nzchar(flowpath_toid) & flowpath_toid != "0",
    .N, by = flowpath_toid]
  expected_junc <- sum(tc_tbl$N[tc_tbl$N > 1L])
  n_junc <- nrow(junc_ev)
  if (n_junc == 0L && expected_junc > 0L) {
    .hf_log_warn(sprintf(
      "[events] 0 junction events generated; network has %d expected inflow segs",
      expected_junc))
  } else {
    message(sprintf(
      "[events] invariant OK: %d junction event(s) (~%d expected inflow segs)",
      n_junc, expected_junc))
  }

  # 5. Summary line always emitted
  message(sprintf(
    "[events] output: %d poi, %d junction (total %d events)",
    nrow(poi_ev), n_junc, nrow(events)))

  if (!is.null(out_gpkg)) sf::st_write(events, out_gpkg, layer = "events", delete_layer = TRUE, quiet = TRUE)
  events
}

# Choose the upstream flowpath a POI at an inlet belongs to.
#
# The mainstem continuation is the right answer: the POI marks a point on a
# channel, and the reach carrying that channel upstream is the one whose outlet it
# is. Matching on levelpathid finds it, with two things the earlier form got
# wrong.
#
# A candidate with no levelpathid must be excluded before comparing, or `== lp`
# yields NA and NA-subscripting injects empty rows into the match set.
#
# Where no upstream shares the mainstem the choice fell to list order, which at a
# confluence can hand a POI a small tributary in place of the channel it belongs
# to: measured on VPU 01, two POIs landed on a 4 km2 and a 9 km2 reach whose
# mainstem carried 903 km2 and 44 km2. Largest drainage area is the mainstem by
# definition, so that is the fallback.
.pick_upstream <- function(idxs, fp_i, dt_fp) {
  if (!length(idxs)) return(NA_integer_)
  if ("levelpathid" %in% names(dt_fp)) {
    lp <- dt_fp$levelpathid[fp_i]
    if (!is.na(lp)) {
      same <- idxs[!is.na(dt_fp$levelpathid[idxs]) & dt_fp$levelpathid[idxs] == lp]
      if (length(same)) idxs <- same
    }
  }
  if (length(idxs) == 1L) return(as.integer(idxs))
  da_col <- intersect(c("tot_drainage_area", "total_dasqkm", "totdasqkm"), names(dt_fp))
  if (length(da_col)) {
    da <- suppressWarnings(as.numeric(dt_fp[[da_col[1]]][idxs]))
    if (any(is.finite(da))) return(as.integer(idxs[which.max(replace(da, !is.finite(da), -Inf))]))
  }
  as.integer(idxs[1L])
}
