# fast_reactor.R
# Refactor engine: reference accuracy guards + O(n) data.table speed.
#
# Two targeted corrections over fast_refactor_nhdplus / fast_collapse_flowlines:
#
#  1. Problem-headwater guard (Phase 2):
#     A headwater whose downstream is a confluence (n_up > 1) OR already has a
#     joined_from marker is flagged joined_to = "-9999" instead of being forcibly
#     donated downstream.  Mirrors the reference
#     reconcile_downstream(remove_problem_headwaters = TRUE).
#
#  2. Stage 4 / Stage 5 threshold split (Phase 3b):
#     Single-upstream short segments (mainstem runs, n_up_Y == 1) are collapsed
#     with mainstem_thresh; multi-upstream (confluence) short segments use thresh.
#     Mirrors the reference reroute_mainstem_set (Stage 4) /
#     reroute_confluence_set (Stage 5).
#
# Everything else: stale-ref resolution, integer-index graph, incremental n_up,
# event protection, parallel split, three-pass, reconcile: is unchanged from
# fast_refactor_nhdplus.

# -- refactor_flowlines -------------------------------------------------------------

#' @title Refactor reference flowlines into reconciled reaches
#' @description Refactors an NHDPlus flowline network by splitting long segments,
#'   collapsing short inter-confluence segments, and reconciling collapsed groups
#'   into a clean topology. Combines the reference accuracy guards with data.table
#'   O(n) performance. Two corrections over fast_collapse_flowlines: (1) a
#'   problem-headwater guard that marks confluence-adjacent headwaters as
#'   uncollapsible ("-9999") rather than forcing them downstream; and (2) a
#'   Stage 4/5 threshold split that applies mainstem_thresh to single-upstream
#'   (mainstem-run) segments and thresh to multi-upstream (confluence) segments.
#' @param flines sf data.frame with columns flowpath_id, flowpath_toid, lengthkm,
#'   hydroseq, levelpathid. totdasqkm optional.
#' @param split_flines_meters numeric max segment length (meters).
#' @param collapse_flines_meters numeric minimum inter-confluence segment length (meters).
#' @param collapse_flines_main_meters numeric minimum mainstem segment length (meters).
#'   Defaults to collapse_flines_meters.
#' @param out_refactored character path for split+collapsed gpkg output.
#' @param out_reconciled character path for reconciled gpkg output.
#' @param three_pass logical run collapse at 10%, 50%, 100% of threshold.
#' @param purge_non_dendritic logical (unused, kept for API compatibility).
#' @param exclude_cats integer/character IDs to protect from collapse.
#' @param events sf point event locations (passed to split_flowlines). Their
#'   geometry also determines which reach outlets a collapse may not merge away.
#' @param cores integer cores for parallel split (uses future/furrr when >1).
#' @param warn logical print progress messages.
#' @return invisibly returns list(refactored, reconciled) as sf data.frames.
#' @examples
#' \dontrun{
#' flines <- sf::st_read("reference_01.gpkg", "flowpaths")
#' out <- refactor_flowlines(flines,
#'                           split_flines_meters    = 2000,
#'                           collapse_flines_meters = 1000,
#'                           out_refactored = "refactored.gpkg",
#'                           out_reconciled = "reconciled.gpkg")
#' }
#' @export
refactor_flowlines <- function(flines,
                               split_flines_meters,
                               collapse_flines_meters,
                               collapse_flines_main_meters = NULL,
                               out_refactored,
                               out_reconciled,
                               three_pass  = TRUE,
                               purge_non_dendritic = TRUE,
                               exclude_cats = NULL,
                               events       = NULL,
                               cores        = 1L,
                               warn         = FALSE) {

  if (is.null(collapse_flines_main_meters))
    collapse_flines_main_meters <- collapse_flines_meters

  # Normalize level-path column. Reference fabrics vary: PRVI uses
  # `levelpathid`, superconus uses `mainstemlp`, some use `levelpath` or
  # `mainstem_id`. Downstream code expects `levelpathid`. Fail loudly if none
  # is present: silently defaulting to NA/1 collapses ALL flowlines into a
  # single level-path group, catastrophically mis-aggregating unrelated reaches.
  lp_src <- intersect(c("levelpathid", "mainstemlp", "levelpath", "mainstem_id"),
    names(flines))
  if (length(lp_src) == 0L) {
    stop("refactor_flowlines:flines must include one of: levelpathid, mainstemlp, ",
      "levelpath, mainstem_id. No level-path column found in: ",
      paste(names(flines), collapse = ", "))
  }
  if (lp_src[1] != "levelpathid") flines$levelpathid <- flines[[lp_src[1]]]
  if (any(is.na(flines$levelpathid))) {
    stop("refactor_flowlines:levelpathid has ",
      sum(is.na(flines$levelpathid)), " NA value(s). NHD input must have a ",
      "valid level path for every flowline.")
  }

  in_crs <- sf::st_crs(flines)
  if (sf::st_is_longlat(flines)) flines <- sf::st_transform(flines, 5070)

  # Reference NHD topology (COMID -> toCOMID, hydroseq) captured before split /
  # collapse, for the reference-faithful reconciled-toid derivation after
  # reconcile. Scientific-notation-safe ids (as.character(1e5) is "1e+05").
  .idc       <- function(x) format(as.numeric(x), scientific = FALSE, trim = TRUE)
  ref_toCOMID <- stats::setNames(.idc(flines$flowpath_toid), .idc(flines$flowpath_id))
  ref_hydseq  <- stats::setNames(as.numeric(flines$hydroseq), .idc(flines$flowpath_id))

  # Drop zero-length flowlines (degenerate source geometries, e.g. AK TDX).
  # These survive st_line_merge and cause downstream segfaults / GDAL write
  # failures. They carry no hydrologic information so removal is safe.
  zero_len_in <- as.numeric(sf::st_length(flines)) == 0
  if (any(zero_len_in, na.rm = TRUE)) {
    message(sprintf("refactor_flowlines: dropped %d zero-length input flowline(s)",
      sum(zero_len_in, na.rm = TRUE)))
    flines <- flines[!zero_len_in, ]
  }

  # -- Step 1: split ----------------------------------------------------------
  if (warn) message("Splitting flowlines to max ", split_flines_meters, " m...")

  if (!is.null(events) && "flowpath_id" %in% names(events)) {
    ev_protect   <- unique(c(events$flowpath_id, events$flowpath_toid))
    exclude_cats <- unique(c(exclude_cats, ev_protect[!is.na(ev_protect)]))
  }

  if (cores > 1L && requireNamespace("furrr", quietly = TRUE)) {
    future::plan(future::multisession, workers = cores)
    on.exit(future::plan(future::sequential), add = TRUE)
  }

  # Force events that would create a sub-collapse-threshold sliver (< collapse_flines_meters)
  # through the "don't split; donate downstream" path (handled just below).
  # Without this guard a parent can split into one long sibling plus a tiny
  # sliver; both land in exclude_cats and neither can collapse, leaving a
  # reconciled-flowline-length floor.
  flines_split <- .fast_split(flines,
    max_length  = split_flines_meters,
    events      = events,
    avoid       = exclude_cats,
    cores       = cores,
    min_split_m = collapse_flines_meters)

  dbg_split_rds <- file.path(dirname(out_refactored), "debug_split_output.rds")
  if (!file.exists(dbg_split_rds))
    saveRDS(flines_split, dbg_split_rds)

  # Prune exclude_cats: if an event flowpath was NOT split (event near an
  # endpoint, suppressed by the 10-m tolerance), it has no interior event and
  # should not be protected from collapse.  Instead, it should donate
  # DOWNSTREAM so the event sits at the divide boundary.
  donate_downstream_cats <- NULL
  if (!is.null(events) && "flowpath_id" %in% names(events) &&
    length(exclude_cats) > 0L) {
    split_ids_post <- as.character(
      sf::st_drop_geometry(flines_split)[["flowpath_id"]])
    dotted_all <- split_ids_post[grepl("\\.", split_ids_post)]
    bases_actually_split <- unique(sub("\\.[0-9]+$", "", dotted_all))
    ev_all_ids <- unique(as.character(
      c(events$flowpath_id, events$flowpath_toid)))
    ev_all_ids <- ev_all_ids[!is.na(ev_all_ids) & nzchar(ev_all_ids)]
    ev_unsplit <- ev_all_ids[!ev_all_ids %in% bases_actually_split]
    if (length(ev_unsplit) > 0L) {
      exclude_cats <- exclude_cats[!exclude_cats %in% ev_unsplit]
      ev_direct <- unique(as.character(events$flowpath_id[!is.na(events$flowpath_id)]))
      donate_downstream_cats <- ev_direct[!ev_direct %in% bases_actually_split]
    }
  }

  # Reaches whose OUTLET carries an event, which is to say the set of nodes the
  # events are indexed to. The collapse cuts its graph at these, so an event node
  # survives every merge by construction rather than by each phase testing for it.
  #
  # Measured by geometry rather than by parsing split suffixes or comparing
  # measures, so one test covers every way an event reaches an endpoint: the
  # near-endpoint cleanup in prepare_hydrolocations snaps an event to its reach's
  # outlet (or reassigns an inlet event to the upstream reach whose outlet that
  # inlet is), and an interior event becomes the outlet of the split child above
  # the split point.
  #
  # An event still sitting at an inlet has no reach whose outlet it is, so it is
  # absent here: that is the headwater case with no upstream to reassign to, and
  # donate_downstream_cats already keeps it on the merged reach's inlet boundary.
  event_outlet_ids <- character(0)
  if (!is.null(events) && nrow(events) > 0L &&
    inherits(events, "sf") && nrow(flines_split) > 0L) {
    ev_geom <- sf::st_geometry(events)
    ev_geom <- ev_geom[!sf::st_is_empty(ev_geom)]
    if (length(ev_geom) > 0L) {
      if (sf::st_crs(ev_geom) != sf::st_crs(flines_split))
        ev_geom <- sf::st_transform(ev_geom, sf::st_crs(flines_split))
      outs <- sf::st_sfc(lapply(sf::st_geometry(flines_split), function(g) {
        m <- sf::st_coordinates(g)
        sf::st_point(unname(m[nrow(m), 1:2]))
      }), crs = sf::st_crs(flines_split))
      near <- sf::st_nearest_feature(outs, ev_geom)
      d    <- as.numeric(sf::st_distance(outs, ev_geom[near], by_element = TRUE))
      # 1 m: the snap placed these exactly on an endpoint, so this is a
      # coordinate-precision tolerance, not a search radius.
      event_outlet_ids <- unique(as.character(
        sf::st_drop_geometry(flines_split)[["flowpath_id"]][d <= 1]))
    }
  }
  if (length(event_outlet_ids) > 0L)
    message(sprintf(paste0("refactor_flowlines: collapse graph cut at %d event ",
      "node(s)"), length(event_outlet_ids)))

  # Expand exclude_cats to include split children of event segments
  if (!is.null(exclude_cats) && length(exclude_cats) > 0L) {
    split_ids_all <- as.character(sf::st_drop_geometry(flines_split)[["flowpath_id"]])
    dotted_set    <- split_ids_all[grepl(".", split_ids_all, fixed = TRUE)]
    if (length(dotted_set) > 0L) {
      dotted_bases <- sub("\\.[0-9]+$", "", dotted_set)
      new_exc      <- dotted_set[dotted_bases %in% exclude_cats]
      if (length(new_exc) > 0L) exclude_cats <- unique(c(exclude_cats, new_exc))
    }
  }

  # -- Step 2: collapse (three-pass) -----------------------------------------
  if (warn) message("Collapsing short segments...")
  thresh_km <- collapse_flines_meters      / 1000
  mstem_km  <- collapse_flines_main_meters / 1000
  flines_attr <- sf::st_drop_geometry(flines_split)

  dbg_collapse_rds <- file.path(dirname(out_refactored), "debug_collapse_input.rds")
  if (!file.exists(dbg_collapse_rds))
    saveRDS(flines_attr, dbg_collapse_rds)

  flines_attr[["flowpath_id"]]   <- as.character(flines_attr[["flowpath_id"]])
  flines_attr[["flowpath_toid"]] <- as.character(flines_attr[["flowpath_toid"]])
  if (!"lengthkm" %in% names(flines_attr)) stop("No lengthkm column found after split.")

  # Propagate parent attributes (levelpathid, hydroseq, totdasqkm) to split children.
  # CRITICAL: look up from the ORIGINAL input (pre-split), not post-split data.
  # split_flowlines replaces parent rows with children.
  if ("flowpath_id" %in% names(flines_attr)) {
    parent_id      <- sub("[.][0-9]+$", "", flines_attr[["flowpath_id"]])
    is_split_child <- parent_id != flines_attr[["flowpath_id"]]
    if (any(is_split_child)) {
      orig_attr <- sf::st_drop_geometry(flines)
      orig_attr[["flowpath_id"]] <- as.character(orig_attr[["flowpath_id"]])
      parent_lu <- stats::setNames(seq_len(nrow(orig_attr)), orig_attr[["flowpath_id"]])
      child_idx <- parent_lu[parent_id[is_split_child]]
      ok        <- !is.na(child_idx)
      for (col in c("levelpathid", "hydroseq", "totdasqkm")) {
        if (any(ok) && col %in% names(flines_attr) && col %in% names(orig_attr)) {
          flines_attr[[col]][is_split_child][ok] <-
            dplyr::coalesce(flines_attr[[col]][is_split_child][ok],
              orig_attr[[col]][child_idx[ok]])
        }
      }
    }
  }

  if (three_pass) {
    ct <- .reactor_collapse(flines_attr, 0.1 * thresh_km, add_category = FALSE,
      mainstem_thresh = 0.1 * mstem_km,
      exclude_cats = exclude_cats,
      donate_downstream_cats = donate_downstream_cats,
      event_outlet_cats = event_outlet_ids,
      warn = FALSE)
    ct <- .reactor_collapse(ct,          0.5 * thresh_km, add_category = FALSE,
      mainstem_thresh = 0.5 * mstem_km,
      exclude_cats = exclude_cats,
      donate_downstream_cats = donate_downstream_cats,
      event_outlet_cats = event_outlet_ids,
      warn = FALSE)
    ct <- .reactor_collapse(ct,                thresh_km, add_category = TRUE,
      mainstem_thresh = mstem_km,
      exclude_cats = exclude_cats,
      donate_downstream_cats = donate_downstream_cats,
      event_outlet_cats = event_outlet_ids,
      warn = warn)
  } else {
    ct <- .reactor_collapse(flines_attr, thresh_km, add_category = TRUE,
      mainstem_thresh = mstem_km,
      exclude_cats = exclude_cats,
      donate_downstream_cats = donate_downstream_cats,
      event_outlet_cats = event_outlet_ids,
      warn = warn)
  }

  flines_split_chr <- flines_split |>
    dplyr::mutate(flowpath_id = hfutils::hf_fmt_id(flowpath_id)) |>
    dplyr::select(flowpath_id)

  refactored_sf <- flines_split_chr |>
    dplyr::inner_join(ct, by = "flowpath_id") |>
    sf::st_as_sf() |>
    sf::st_transform(in_crs)
  refactored_sf$lengthkm <- as.numeric(sf::st_length(refactored_sf)) / 1000

  sf::st_write(refactored_sf, out_refactored,
    layer_options = "OVERWRITE=YES", quiet = !warn, append = FALSE)

  if (warn) message("Collapse complete. Reconciling...")

  # -- Step 3: reconcile -----------------------------------------------------
  dbg_rds <- file.path(dirname(out_refactored), "debug_reconcile_input.rds")
  if (!file.exists(dbg_rds))
    saveRDS(list(ct = ct, geom = flines_split_chr), dbg_rds)
  reconciled <- reconcile_collapsed_flowlines(ct,
    geom = flines_split_chr,
    id   = "flowpath_id")
  if (!identical(sf::st_crs(reconciled), in_crs))
    reconciled <- sf::st_transform(reconciled, in_crs)

  # Flow-direction fix: hfutils::union_linestrings / st_line_merge can emit
  # reversed geometry when merging multiple collapsed pieces whose original
  # orientations differ. Detect and repair at source so downstream code can
  # rely on the upstream->downstream convention (add_flowlines chain building,
  # extend_flowpaths_to_mainstem, etc.).
  # Reference-faithful toid first (authoritative: follows reference COMID->toCOMID at
  # confluences, where geometry alone can't tell mainstem from tributary), then
  # orient geometry to match.
  reconciled <- .fix_reconciled_toid_from_reference(reconciled, ref_toCOMID, ref_hydseq)
  reconciled <- .fix_reconciled_flowdir(reconciled, warn = warn)

  # Stringify member_flowpath_id (list -> CSV) without injecting the unsplit
  # parent alongside its dotted children. Adding the bare NHD ID next to its
  # "N.1/N.2/..." siblings double-counts the parent during the aggregation
  # member_comid concat: a parent landing in two reconciled groups appears in
  # both aggregated flowpaths' member lists, inflating drainage lineage and
  # producing fp<->div join artifacts.
  #
  # Downstream consumers that genuinely need base-form lookup (event snap,
  # add_flowlines junction mapping, build_aspatial_network reference_id
  # derivation) already strip ".N" suffixes themselves.
  reconciled[["member_flowpath_id"]] <- vapply(
    reconciled$member_flowpath_id,
    function(x) paste(x, collapse = ","),
    character(1L))

  # Homogenise geometry type before writing: mixed LINESTRING/MULTILINESTRING
  # in a single sfc causes GDAL/GEOS to segfault on Apple Silicon (GEOS 3.14).
  # Cast all to MULTILINESTRING so the layer has a single geometry type.
  # Downstream readers that need LINESTRING cast back with st_cast("LINESTRING").
  rec_geom <- sf::st_geometry(reconciled)
  if (length(unique(as.character(sf::st_geometry_type(rec_geom)))) > 1L) {
    sf::st_geometry(reconciled) <- sf::st_cast(rec_geom, "MULTILINESTRING",
      warn = FALSE)
  }

  sf::st_write(reconciled, out_reconciled,
    layer_options = "OVERWRITE=YES", quiet = !warn, append = FALSE)

  # -- Step 4: snap events to reconciled outlets ----------------------------
  reconciled_events <- NULL
  if (!is.null(events) && nrow(events) > 0L) {
    reconciled_events <- .snap_events_to_reconciled(
      events      = events,
      reconciled  = reconciled,
      snap_dist_m = collapse_flines_meters,
      in_crs      = in_crs
    )
    if (!is.null(reconciled_events) && nrow(reconciled_events) > 0L) {
      tryCatch(
        sf::st_write(reconciled_events, out_reconciled, layer = "events",
          append = FALSE),
        error = function(e) {
          sf::st_write(reconciled_events, out_reconciled, layer = "events",
            delete_layer = TRUE)
        }
      )
    }
  }

  invisible(list(refactored = refactored_sf, reconciled = reconciled,
    events = reconciled_events))
}

# -- .reactor_collapse ---------------------------------------------------------
# Internal collapse engine. Identical to fast_collapse_flowlines except:
#
#  Change 1: Phase 2 problem-headwater guard:
#    Before donating headwater X to its downstream Y, check:
#      (a) n_up[Y] > 1  (Y is a confluence), or
#      (b) Y already has a joined_from marker (prior or current pass).
#    If either holds, set absorber_to[X] = -9999L (sentinel) and skip,
#    do NOT donate or remove X.  is.na(absorber_to) then excludes X from
#    any_hw() on subsequent iterations, exactly as the reference -9999 guard.
#
#  Change 2: Phase 3b Stage 4/5 threshold split:
#    any_conf() now emits two candidate sets:
#      Stage 4 (mainstem run): n_up[Y] == 1, len[Y] < mainstem_thresh
#      Stage 5 (confluence):   n_up[Y] >  1, len[Y] < thresh
#    The loop body is unchanged; category tracking separates ms_run_tracker
#    from conf_tracker and both roll into the "mainstem"/"confluence" fcase.
#
#  Write-back adds -9999 handling for absorber_to (problem headwater sentinels).

.reactor_collapse <- function(flines, thresh,
                              add_category    = FALSE,
                              mainstem_thresh = NULL,
                              exclude_cats    = NULL,
                              donate_downstream_cats = NULL,
                              event_outlet_cats      = NULL,
                              warn            = TRUE) {
  # -- Setup ------------------------------------------------------------------
  dt <- data.table::as.data.table(flines)

  keep_cols  <- c("flowpath_id", "flowpath_toid", "levelpathid", "hydroseq",
    "lengthkm", "totdasqkm",
    "joined_from_flowpath_id", "joined_to_flowpath_id")
  extra_cols <- setdiff(names(dt), keep_cols)
  if (length(extra_cols) > 0L) dt[, (extra_cols) := NULL]

  if (!"joined_from_flowpath_id" %in% names(dt)) dt[, joined_from_flowpath_id := NA_character_]
  if (!"joined_to_flowpath_id"   %in% names(dt)) dt[, joined_to_flowpath_id   := NA_character_]
  dt[, flowpath_id   := hfutils::hf_fmt_id(flowpath_id)]
  dt[, flowpath_toid := hfutils::hf_fmt_id(flowpath_toid)]
  dt[flowpath_toid == "0" | flowpath_toid == "NA", flowpath_toid := NA_character_]

  # TRUE outlets: original toCOMID was NA/0. Stale refs (point to non-existent IDs)
  # are NOT true outlets: resolve them first.
  true_outlet_orig <- is.na(dt$flowpath_toid)

  comid_set  <- dt$flowpath_id
  stale_mask <- !is.na(dt$flowpath_toid) & !dt$flowpath_toid %in% comid_set
  if (any(stale_mask)) {
    has_dot <- grepl(".", comid_set, fixed = TRUE)
    if (any(has_dot)) {
      split_dt <- data.table::data.table(
        flowpath_id = comid_set[has_dot],
        base        = sub("\\..*$", "", comid_set[has_dot]),
        part        = as.integer(sub("^[^.]+\\.", "", comid_set[has_dot]))
      )
      split_ends <- split_dt[, .(downstream_COMID = flowpath_id[which.min(part)]), by = base]
      stale_map  <- stats::setNames(split_ends$downstream_COMID, split_ends$base)
      stale_tc   <- dt$flowpath_toid[stale_mask]
      resolved   <- stale_map[sub("\\..*$", "", stale_tc)]
      dt$flowpath_toid[stale_mask] <- ifelse(!is.na(resolved), resolved, stale_tc)
    }
    still_stale <- !is.na(dt$flowpath_toid) & !dt$flowpath_toid %in% comid_set
    dt[still_stale, flowpath_toid := NA_character_]
  }

  if (!"levelpathid" %in% names(dt)) dt[, levelpathid := NA_real_]
  if (!"hydroseq"   %in% names(dt)) dt[, hydroseq    := seq_len(.N)]

  if (is.null(mainstem_thresh)) mainstem_thresh <- max(dt$lengthkm, na.rm = TRUE)
  if (!is.null(exclude_cats)) exclude_cats    <- as.character(exclude_cats)

  # -- Integer-index graph ----------------------------------------------------
  n       <- nrow(dt)
  pos_map <- stats::setNames(seq_len(n), dt$flowpath_id)
  to_pos  <- .lookup_pos(dt$flowpath_toid, pos_map)   # integer, 0 = outlet/NA
  len     <- dt$lengthkm
  lp      <- as.numeric(dt$levelpathid)
  lp[is.na(lp)] <- Inf
  excl    <- rep(FALSE, n)
  if (!is.null(exclude_cats)) excl <- dt$flowpath_id %in% exclude_cats
  donate_ds <- rep(FALSE, n)
  if (!is.null(donate_downstream_cats))
    donate_ds <- dt$flowpath_id %in% as.character(donate_downstream_cats)

  # Cut the graph at every node an event is indexed to.
  #
  # A collapse merges two reaches, which destroys the node between them. A POI on
  # that node ends up inside the merged reach, so it is no longer an outlet and
  # gets no divide. What has to survive is therefore a NODE, while membership of
  # a set like exclude_cats is a property of a REACH, which is why per-site
  # reach tests could not express it consistently and disagreed between phases.
  #
  # Rather than have each phase ask permission before merging, the edges that
  # could violate the property are removed from the graph the phases run on: a
  # reach whose outlet carries an event has no downstream, so the reach below the
  # node has one fewer upstream. No phase can merge across an edge that is not
  # there, and no phase needs to know that events exist. The collapse then does
  # its ordinary job independently within each chain the cuts leave behind.
  #
  # A node is cut once for every reach that ends on it, so a confluence carrying
  # an event is fully severed: not just the inflow the event was indexed to, but
  # every sibling inflow, which would otherwise be free to absorb the reach below
  # the node and destroy it from the side.
  #
  # `dn_cut` keeps the real downstream so the topology can be restored on the way
  # out. Cutting after `true_outlet_orig` is read above matters: these reaches are
  # terminal only for the duration of the collapse, not genuinely.
  dn_cut <- rep(0L, n)
  if (!is.null(event_outlet_cats)) {
    cut_i <- which(dt$flowpath_id %in% as.character(event_outlet_cats) &
      to_pos > 0L)
    if (length(cut_i) > 0L) {
      dn_cut[cut_i] <- to_pos[cut_i]
      to_pos[cut_i] <- 0L
    }
  }

  removed <- rep(FALSE, n)

  already_absorbed <- dt$lengthkm == 0 &
    (!is.na(dt$joined_from_flowpath_id) | !is.na(dt$joined_to_flowpath_id))
  removed[already_absorbed] <- TRUE

  absorber_from <- rep(NA_integer_, n)
  absorber_to   <- rep(NA_integer_, n)
  n_up          <- tabulate(to_pos[to_pos > 0L], nbins = n)

  # -- Helpers ----------------------------------------------------------------
  best_lp_idx <- function(idxs) idxs[which.min(lp[idxs])]

  update_n_up <- function(i, old_t, new_t) {
    if (old_t > 0L) n_up[old_t] <<- n_up[old_t] - 1L
    if (new_t > 0L) n_up[new_t] <<- n_up[new_t] + 1L
  }

  absorb_upstream <- function(Y_i, X_i) {
    len[X_i]     <<- len[X_i] + len[Y_i]
    len[Y_i]     <<- 0
    removed[Y_i] <<- TRUE
    absorber_from[Y_i] <<- X_i
    old_t <- to_pos[X_i]
    if (old_t == Y_i) {
      update_n_up(X_i, old_t, 0L)
      to_pos[X_i] <<- 0L
    }
    if (to_pos[Y_i] > 0L) n_up[to_pos[Y_i]] <<- n_up[to_pos[Y_i]] - 1L
    n_up[Y_i] <<- 0L
  }

  absorb_downstream <- function(X_i, Y_i) {
    y_to <- to_pos[Y_i]
    len[X_i]     <<- len[X_i] + len[Y_i]
    len[Y_i]     <<- 0
    removed[Y_i] <<- TRUE
    absorber_to[Y_i] <<- X_i
    old_x_to <- to_pos[X_i]
    update_n_up(X_i, old_x_to, y_to)
    to_pos[X_i] <<- y_to
    n_up[Y_i]   <<- 0L
    if (y_to > 0L) n_up[y_to] <<- n_up[y_to] - 1L
  }

  skip_past <- function(X_i, Y_i) {
    y_to     <- to_pos[Y_i]
    old_x_to <- to_pos[X_i]
    update_n_up(X_i, old_x_to, y_to)
    to_pos[X_i] <<- y_to
  }

  donate_to_next <- function(X_i, Y_i) {
    len[Y_i]     <<- len[Y_i] + len[X_i]
    len[X_i]     <<- 0
    removed[X_i] <<- TRUE
    absorber_from[X_i] <<- Y_i
    n_up[Y_i]    <<- n_up[Y_i] - 1L
  }

  # -- Phase 1: Short Outlets -------------------------------------------------
  # Handles TRUE original outlets AND stranded "pseudo-outlets" whose
  # downstream was removed in a prior pass. Both are safe to absorb into
  # upstream mainstem: the absorbing X becomes the new terminal, which matches
  # the actual post-collapse network state. Without pseudo-outlets, short
  # terminals whose downstream was already removed survive with no downstream
  # and no collapse path.
  {
    any_short_outlets <- function()
      which(!removed & len < thresh & to_pos == 0L & !excl)

    so_tracker <- character()
    prev_so    <- -1L
    so         <- any_short_outlets()
    iter       <- 0L

    while (length(so) > 0L) {
      if (length(so) == prev_so) break
      prev_so <- length(so)
      for (Y_i in so) {
        if (removed[Y_i]) next
        ups <- which(!removed & to_pos == Y_i)
        if (length(ups) == 0L) {
          absorber_from[Y_i] <- -9999L
          removed[Y_i]       <- TRUE
          n_up[Y_i]          <- 0L
          so_tracker <- c(so_tracker, dt$flowpath_id[Y_i])
          next
        }
        X_i <- best_lp_idx(ups)
        absorb_upstream(Y_i, X_i)
        for (o in setdiff(ups, X_i)) {
          to_pos[o]  <- X_i
          n_up[X_i]  <- n_up[X_i] + 1L
        }
        so_tracker <- c(so_tracker, dt$flowpath_id[Y_i])
      }
      iter <- iter + 1L
      if (iter > 200L) {
        .hf_log_warn("Phase 1 stuck after 200 iterations")
        break
      }
      so <- any_short_outlets()
    }
  }

  # -- Phase 2: Short Headwaters WITH problem-headwater guard -----------------
  # Change 1: before donating X -> Y, check whether Y is a confluence (n_up > 1)
  # or already carries a joined_from marker from a prior or current pass.
  # If so, mark absorber_to[X] = -9999L (sentinel) and skip: do NOT donate/remove.
  # is.na(absorber_to) in any_hw() then prevents X from being re-selected.
  # Matches the reference reconcile_downstream(remove_problem_headwaters = TRUE).
  {
    hw_tracker <- character()

    any_hw <- function()
      which(!removed & n_up == 0L & len < thresh & to_pos > 0L & !excl &
        is.na(absorber_from) & is.na(absorber_to))

    prev_hw <- -1L
    hws     <- any_hw()
    iter    <- 0L

    while (length(hws) > 0L) {
      if (length(hws) == prev_hw) break
      prev_hw <- length(hws)

      for (X_i in hws) {
        if (removed[X_i]) next
        Y_i <- to_pos[X_i]
        if (Y_i == 0L || removed[Y_i]) next

        # Problem-headwater guard (reference accuracy guard):
        # Y is "problematic" if it is a confluence OR already has a
        # joined_from marker.  But the concern is burying an EVENT boundary;
        # for plain NHD confluences with no event association (Y not in
        # exclude_cats), donating a short raw-NHD headwater only grows Y's
        # length: Y's outlet and other upstreams are untouched, and no event
        # point exists at the old inlet.  Without this relaxation many unsplit
        # short NHD headwaters never collapse.
        already_jf <- (!is.na(dt$joined_from_flowpath_id[Y_i]) &&
          dt$joined_from_flowpath_id[Y_i] != "-9999") ||
          !is.na(absorber_from[Y_i])

        if ((n_up[Y_i] > 1L || already_jf) && excl[Y_i]) {
          # Event-adjacent downstream: preserve confluence as a node.
          absorber_to[X_i] <- -9999L
          next
        }

        donate_to_next(X_i, Y_i)
        hw_tracker <- c(hw_tracker, dt$flowpath_id[X_i])
      }

      iter <- iter + 1L
      if (iter > 200L) {
        .hf_log_warn("Phase 2 stuck after 200 iterations")
        break
      }
      hws <- any_hw()
    }
  }

  # -- Phase 3a: Short Mainstem Tops ------------------------------------------
  # Identical to fast_collapse_flowlines.
  mst_tracker <- character()
  if (!is.null(mainstem_thresh) && mainstem_thresh < max(dt$lengthkm, na.rm = TRUE)) {
    any_mst <- function() {
      ai <- which(!removed & !excl & len < mainstem_thresh & n_up > 1L & to_pos > 0L)
      ai[!removed[to_pos[ai]] & n_up[to_pos[ai]] == 1L]
    }
    prev_mst <- -1L
    msts     <- any_mst()
    iter     <- 0L

    while (length(msts) > 0L) {
      if (length(msts) == prev_mst) break
      prev_mst <- length(msts)
      for (X_i in msts) {
        if (removed[X_i]) next
        Y_i <- to_pos[X_i]
        if (Y_i == 0L || removed[Y_i] || n_up[Y_i] != 1L) next
        if (excl[Y_i]) next
        absorb_downstream(X_i, Y_i)
        mst_tracker <- c(mst_tracker, dt$flowpath_id[Y_i])
      }
      iter <- iter + 1L
      if (iter > 200L) {
        .hf_log_warn("Phase 3a stuck")
        break
      }
      msts <- any_mst()
    }
  }

  # -- Phase 3b: Stage 4 (mainstem run) + Stage 5 (confluence) ---------------
  # Change 2: any_conf() emits two candidate sets distinguished by n_up[Y]:
  #   Stage 4: n_up[Y] == 1 (straight mainstem run): use mainstem_thresh
  #   Stage 5: n_up[Y] >  1 (confluence):            use thresh
  # Matches the reference reroute_mainstem_set / reroute_confluence_set.
  # Loop body and event-exclusion special case are unchanged.
  {
    ms_run_tracker <- character()
    conf_tracker   <- character()

    any_conf <- function() {
      # Stage 4: single-upstream Y short by mainstem_thresh (mainstem run)
      # !excl: event-protected segments must never be collapse candidates, else
      # a short event-bearing segment gets absorbed into an adjacent segment,
      # merging its geometry and divide into the wrong reconciled group.
      ai_ms <- which(!removed & !excl & n_up == 1L & len < mainstem_thresh & to_pos > 0L)
      ds_ms <- to_pos[ai_ms]
      ai_ms <- ai_ms[ds_ms > 0L & !removed[ds_ms]]

      # Stage 5: multi-upstream Y short by thresh (confluence)
      ai_c  <- which(!removed & !excl & n_up > 1L & len < thresh & to_pos > 0L)
      ds_c  <- to_pos[ai_c]
      ai_c  <- ai_c[ds_c > 0L & !removed[ds_c]]

      unique(c(ai_ms, ai_c))
    }

    prev_conf <- -1L
    cands     <- any_conf()
    iter      <- 0L

    while (length(cands) > 0L) {
      if (length(cands) == prev_conf) break
      prev_conf <- length(cands)

      for (Y_i in cands) {
        if (removed[Y_i]) next
        ds_Y <- to_pos[Y_i]
        if (ds_Y == 0L || removed[ds_Y]) next

        # Event-protection: short segment between a junction and an excluded
        # (event-bearing) downstream -> donate Y into that excluded segment.
        if (excl[ds_Y]) {
          ups_Y_all <- which(!removed & to_pos == Y_i)
          len[ds_Y]        <- len[ds_Y] + len[Y_i]
          len[Y_i]         <- 0
          removed[Y_i]     <- TRUE
          absorber_from[Y_i] <- ds_Y
          n_up[ds_Y]       <- n_up[ds_Y] - 1L
          for (o in ups_Y_all) {
            if (!removed[o]) skip_past(o, Y_i)
          }
          n_up[Y_i]        <- 0L
          conf_tracker     <- c(conf_tracker, dt$flowpath_id[Y_i])
          next
        }

        # Endpoint-event segments: the event is at the inlet (upstream end)
        # so donating Y downstream keeps the event at a divide boundary.
        if (donate_ds[Y_i]) {
          ups_Y_all <- which(!removed & to_pos == Y_i)
          len[ds_Y]        <- len[ds_Y] + len[Y_i]
          len[Y_i]         <- 0
          removed[Y_i]     <- TRUE
          absorber_from[Y_i] <- ds_Y
          n_up[ds_Y]       <- n_up[ds_Y] - 1L
          for (o in ups_Y_all) {
            if (!removed[o]) skip_past(o, Y_i)
          }
          n_up[Y_i]        <- 0L
          conf_tracker     <- c(conf_tracker, dt$flowpath_id[Y_i])
          next
        }

        ups_Y <- which(!removed & to_pos == Y_i)
        if (length(ups_Y) == 0L) next

        if (length(ups_Y) == 1L) {
          X_i <- ups_Y[1L]
          absorb_downstream(X_i, Y_i)
          ms_run_tracker <- c(ms_run_tracker, dt$flowpath_id[Y_i])
        } else {
          X_i <- best_lp_idx(ups_Y)
          absorb_downstream(X_i, Y_i)
          # Y's geometry and length have just moved INTO X_i, so X_i's feature
          # now spans [X_i .. Y]. A sibling upstream reach joined the network at
          # Y's inlet, which is a point in the middle of that merged feature,
          # so it must enter X_i. Routing it past to Y's old downstream (the
          # skip_past used when Y is donated DOWNSTREAM, where Y's geometry does
          # end up) would make its water bypass a real channel segment and move
          # its confluence downstream. This mirrors Phase 1's absorb_upstream
          # handling of the same situation; n_up[Y_i] is already 0 there, so
          # only the new target is incremented.
          for (o in setdiff(ups_Y, X_i)) {
            if (!removed[o]) {
              to_pos[o] <- X_i
              n_up[X_i] <- n_up[X_i] + 1L
            }
          }
          conf_tracker <- c(conf_tracker, dt$flowpath_id[Y_i])
        }
      }

      iter <- iter + 1L
      if (iter > 500L) {
        .hf_log_warn("Phase 3b stuck after 500 iterations")
        break
      }
      cands <- any_conf()
    }
  }

  # -- Write back -------------------------------------------------------------
  dt[, lengthkm := len]

  # Restore the edges cut at event nodes. The reach below a cut may itself have
  # been absorbed while the cut was in force, so follow the removal markers to
  # whichever reach its geometry ended up in; absorber_to and absorber_from both
  # name the survivor, and -9999 marks a reach deleted outright, which leaves a
  # genuine terminal behind. The de-alias pass below cannot do this job: it only
  # follows joined_to, and a severed upstream is never rerouted by skip_past
  # because the phases could not see it.
  if (any(dn_cut > 0L)) {
    survivor_of <- function(i) {
      for (k in seq_len(200L)) {
        if (is.na(i) || i <= 0L || !removed[i]) return(if (is.na(i)) 0L else i)
        nxt <- if (!is.na(absorber_to[i])) absorber_to[i] else absorber_from[i]
        if (is.na(nxt) || nxt <= 0L) return(0L)
        i <- nxt
      }
      0L
    }
    for (i in which(dn_cut > 0L)) {
      s <- survivor_of(dn_cut[i])
      to_pos[i] <- if (s == i) 0L else s
    }
  }

  safe_pos    <- pmax(to_pos, 1L)
  new_toCOMID <- dt$flowpath_id[safe_pos]
  new_toCOMID[to_pos == 0L &  true_outlet_orig] <- "0"           # genuine terminal
  new_toCOMID[to_pos == 0L & !true_outlet_orig] <- NA_character_  # stale/pseudo-outlet
  dt[, flowpath_toid := new_toCOMID]

  # joined_from_flowpath_id: preserve prior-pass markers, add new this pass.
  new_jfc          <- dt$joined_from_flowpath_id
  is_minus9_from   <- !is.na(absorber_from) & absorber_from == -9999L
  is_valid_from    <- !is.na(absorber_from) & absorber_from >  0L
  new_jfc[is_minus9_from] <- "-9999"
  new_jfc[is_valid_from]  <- dt$flowpath_id[absorber_from[is_valid_from]]
  dt[, joined_from_flowpath_id := new_jfc]

  # joined_to_flowpath_id: preserve prior-pass markers; add new + "-9999" sentinels.
  new_jtc          <- dt$joined_to_flowpath_id
  is_minus9_to     <- !is.na(absorber_to) & absorber_to == -9999L  # problem headwaters
  is_valid_to      <- !is.na(absorber_to) & absorber_to >  0L
  new_jtc[is_minus9_to] <- "-9999"
  new_jtc[is_valid_to]  <- dt$flowpath_id[absorber_to[is_valid_to]]
  dt[, joined_to_flowpath_id := new_jtc]

  # Chain-chase stale joined_to pointers (A->B->C becomes A->C).
  for (i in seq_len(30L)) {
    j_to  <- dt$joined_to_flowpath_id
    src   <- which(!is.na(j_to) & j_to %in% dt$flowpath_id)
    if (length(src) == 0L) break
    new_to  <- dt$joined_to_flowpath_id[match(j_to[src], dt$flowpath_id)]
    updated <- !is.na(new_to) & new_to != "-9999"
    if (!any(updated, na.rm = TRUE)) break
    dt$joined_to_flowpath_id[src[updated]] <- new_to[updated]
  }

  # Fix reciprocal bad pairs (identical to fast_collapse_flowlines cleanup).
  potential_bad_rows <- dt$flowpath_id[
    !is.na(dt$joined_to_flowpath_id) & dt$joined_to_flowpath_id != "-9999"]
  bad_jfc <- potential_bad_rows[potential_bad_rows %in% dt$joined_from_flowpath_id]
  if (length(bad_jfc) > 0L) {
    bfc_idx <- match(bad_jfc, dt$joined_from_flowpath_id)
    new_jtc <- dt$joined_to_flowpath_id[match(bad_jfc, dt$flowpath_id)]
    dt$joined_to_flowpath_id[bfc_idx]  <- new_jtc
    dt$joined_from_flowpath_id[bfc_idx] <- NA_character_
  }

  # De-alias flowpath_toid through joined_to pointers.
  alias_map <- stats::setNames(dt$joined_to_flowpath_id, dt$flowpath_id)
  new_tc    <- alias_map[dt$flowpath_toid]
  valid     <- !is.na(new_tc) & new_tc != "-9999"
  dt[valid, flowpath_toid := new_tc[valid]]

  dt <- unique(dt, by = "flowpath_id")

  if (add_category) {
    ms_set_all <- unique(c(mst_tracker, ms_run_tracker))
    dt[, join_category := data.table::fcase(
      flowpath_id %in% so_tracker,  "outlet",
      flowpath_id %in% ms_set_all,  "mainstem",
      flowpath_id %in% conf_tracker, "confluence",
      flowpath_id %in% hw_tracker,  "headwater",
      default = NA_character_
    )]
  }

  as.data.frame(dt)
}

# Flow-direction repair for reconciled flowlines. For each line with toID in
# the set, reverse its geometry if its start-point is closer to the toID
# line's start-point than its end-point is (i.e. line is drawn
# downstream-first). Runs a bounded number of passes so that reversing a line
# can still cascade to fix lines that depended on its orientation.
.fix_reconciled_flowdir <- function(rec_sf, tol_m = 1, max_passes = 5L,
                                    warn = FALSE) {
  if (!inherits(rec_sf, "sf") || nrow(rec_sf) == 0L) return(rec_sf)
  if (!all(c("ID", "toID") %in% names(rec_sf))) return(rec_sf)

  # Reconcile stores each reach as a single-path MULTILINESTRING and downstream
  # steps (generate_events) require that type, so keep `geom` as the raw
  # MULTILINESTRING and return it. lwgeom's endpoint/startpoint segfault on
  # MULTILINESTRING, so detection and the orientation/re-derivation passes run on
  # a merged LINESTRING copy `mgeom`; every flip is applied to BOTH so they stay
  # in lockstep. Branching multi-part geometry does not merge to a single
  # LINESTRING and stays excluded (no single orientation is meaningful for a
  # branch); zero-length lines are excluded too (no direction).
  # The raw output is a MIXED sfc_GEOMETRY (LINESTRING + single-path
  # MULTILINESTRING); st_line_merge.sfc requires sfc_MULTILINESTRING, so cast to a
  # uniform MULTILINESTRING first, then merge to a LINESTRING copy for endpoint ops.
  geom  <- sf::st_cast(sf::st_geometry(rec_sf), "MULTILINESTRING", warn = FALSE)
  mgeom <- sf::st_line_merge(geom)
  ids   <- rec_sf$ID
  toids <- rec_sf$toID
  total_reversed <- 0L

  is_line <- as.character(sf::st_geometry_type(mgeom)) == "LINESTRING" &
    as.numeric(sf::st_length(mgeom)) > 0

  # Iterate orientation + toid re-derivation to a FIXED POINT. The re-derivation
  # step (below) can reverse reaches and re-point toids AFTER the orientation
  # pass, leaving a reach mis-oriented relative to its NEW toid (e.g. one whose
  # reference toid was set to a sibling it meets at the start, not the end). A
  # single orientation->re-derivation pass therefore does not converge; loop until
  # a round makes no flips and no toid changes (bounded by max_passes).
  for (outer in seq_len(max_passes)) {
    round_flips <- 0L
    round_chg   <- 0L

    for (pass in seq_len(max_passes)) {
      has_ds  <- is_line & (toids %in% ids)
      chain   <- which(has_ds)
      j_idx   <- match(toids[chain], ids)

      ep    <- sf::st_coordinates(lwgeom::st_endpoint(mgeom[chain]))[, 1:2, drop = FALSE]
      sp    <- sf::st_coordinates(lwgeom::st_startpoint(mgeom[chain]))[, 1:2, drop = FALSE]
      sp_ds <- sf::st_coordinates(lwgeom::st_startpoint(mgeom[j_idx]))[, 1:2, drop = FALSE]

      d_ep <- sqrt(rowSums((ep - sp_ds)^2))
      d_sp <- sqrt(rowSums((sp - sp_ds)^2))

      # Flip rows where startpoint is closer to downstream startpoint than endpoint is
      wrong <- chain[d_sp + tol_m < d_ep]
      if (length(wrong) == 0L) break

      geom[wrong]  <- sf::st_reverse(geom[wrong])
      mgeom[wrong] <- sf::st_reverse(mgeom[wrong])
      total_reversed <- total_reversed + length(wrong)
      round_flips    <- round_flips + length(wrong)
      if (warn) message(sprintf(
        "  .fix_reconciled_flowdir outer %d pass %d: reversed %d line(s)", outer, pass, length(wrong)))
    }

    # -- Toid re-derivation ----------------------------------------------------
    # Orientation alone cannot fix a WRONG toid: a reach wired to a confluence
    # sibling (which it never touches) or a reversed headwater wired to a sibling
    # has no orientation that makes its end meet the toid's start. For each reach
    # whose stored toid still fails to connect (end !~ toid.start, > conn_tol),
    # re-derive from node geometry: the unique reach whose START sits on this
    # reach's connected end IS the true downstream (sibling inflows END at the
    # node, so they are never candidates). Reverse the reach if the connection is
    # at its stored start. Guarded by conn_tol so the (already-valid) ~98% of
    # edges are never touched; ambiguous (divergent) nodes are left as-is.
    conn_tol <- 100
    # Operate on the LINESTRING subset directly (reuse the orientation pass's mask):
    # the reconcile output is a mixed sfc_GEOMETRY, which both st_line_merge and
    # lwgeom's endpoint funcs reject, but lwgeom works on the LINESTRING subset, as
    # the orientation pass proves. Single-part-MULTILINESTRING reaches are left as-is.
    li <- which(is_line)
    if (length(li) > 1L) {
      li_ids <- as.character(ids[li])
      sp <- sf::st_coordinates(lwgeom::st_startpoint(mgeom[li]))[, 1:2, drop = FALSE]
      ep <- sf::st_coordinates(lwgeom::st_endpoint(mgeom[li]))[, 1:2, drop = FALSE]
      rownames(sp) <- rownames(ep) <- li_ids
      nkey <- function(x, y) paste(round(x / tol_m), round(y / tol_m))
      starts_at <- split(li_ids, nkey(sp[, 1], sp[, 2]))   # node -> reach ids STARTING there
      lp_of <- stats::setNames(as.character(rec_sf$levelpathid)[li], li_ids)
      toid_chr <- as.character(toids)
      orig_toids <- toids
      pick <- function(cands, rid) {
        cands <- cands[cands != rid]
        if (length(cands) == 0L) return(NA_character_)
        if (length(cands) == 1L) return(cands)
        same <- cands[lp_of[cands] == lp_of[rid]]            # prefer same levelpath
        if (length(same) == 1L) same else NA_character_      # else ambiguous -> skip
      }
      chg_idx <- integer(0)
      chg_rev <- logical(0)
      for (k in seq_along(li)) {
        rid <- li_ids[k]
        t <- toid_chr[li[k]]
        if (is.na(t) || !(t %in% li_ids)) next               # terminal / non-line toid
        d_now <- sqrt(sum((ep[k, ] - sp[match(t, li_ids), ])^2))
        if (!is.finite(d_now) || d_now <= conn_tol) next      # already connects -> leave
        cand_end   <- pick(starts_at[[nkey(ep[k, 1], ep[k, 2])]], rid)
        cand_start <- pick(starts_at[[nkey(sp[k, 1], sp[k, 2])]], rid)
        if (!is.na(cand_end) && is.na(cand_start)) {
          toids[li[k]] <- as.numeric(cand_end)
          chg_idx <- c(chg_idx, li[k])
          chg_rev <- c(chg_rev, FALSE)
        } else if (is.na(cand_end) && !is.na(cand_start)) {
          geom[li[k]]  <- sf::st_reverse(geom[li[k]])
          mgeom[li[k]] <- sf::st_reverse(mgeom[li[k]])
          toids[li[k]] <- as.numeric(cand_start)
          chg_idx <- c(chg_idx, li[k])
          chg_rev <- c(chg_rev, TRUE)
        }
        # both / neither -> ambiguous or genuine terminal: leave unchanged
      }
      # -- Cycle guard ----------------------------------------------------------
      # The original reconciled topology is a DAG. A re-derived edge can only break
      # that for pairs that connect the same two nodes (A.end=B.start & B.end=A.start),
      # where re-deriving both points them at each other (a 2-cycle). Revert any
      # re-derived reach that lands in a strongly-connected component (toid + the
      # reversal) until the graph is acyclic again.
      if (length(chg_idx) > 0L) {
        repeat {
          e <- data.frame(f = as.character(ids),
            t = as.character(as.integer(toids)), stringsAsFactors = FALSE)
          e <- e[!is.na(e$t) & e$t %in% as.character(ids), ]
          sc <- igraph::components(igraph::graph_from_data_frame(e, directed = TRUE), "strong")
          cyc_nodes <- names(sc$membership)[sc$membership %in% which(sc$csize > 1L)]
          rv <- which(as.character(ids[chg_idx]) %in% cyc_nodes)
          if (length(rv) == 0L) break
          j <- rv[1L]                                  # revert ONE per cycle; re-check
          toids[chg_idx[j]] <- orig_toids[chg_idx[j]]  # frees the other end of a 2-cycle
          if (chg_rev[j]) {
            geom[chg_idx[j]]  <- sf::st_reverse(geom[chg_idx[j]])
            mgeom[chg_idx[j]] <- sf::st_reverse(mgeom[chg_idx[j]])
          }
          chg_idx <- chg_idx[-j]
          chg_rev <- chg_rev[-j]
        }
      }
      if (length(chg_idx) > 0L) {
        message(sprintf(
          "refactor_flowlines: re-derived toid for %d disconnected reach(es) (%d via reversal)",
          length(chg_idx), sum(chg_rev)))
      }
      round_chg <- length(chg_idx)
    }

    if (round_flips == 0L && round_chg == 0L) break
  }

  # Rebuild as a clean, uniform sfc_MULTILINESTRING: in-place element flips above
  # leave the sfc class stale/mixed, and downstream (generate_events) asserts
  # sfc_MULTILINESTRING. st_sfc recomputes the class from the elements; st_cast
  # re-wraps any element st_reverse returned as LINESTRING back to MULTILINESTRING.
  geom <- sf::st_sfc(geom, crs = sf::st_crs(rec_sf))
  if (!all(as.character(sf::st_geometry_type(geom)) == "MULTILINESTRING"))
    geom <- sf::st_cast(geom, "MULTILINESTRING")
  sf::st_geometry(rec_sf) <- geom
  rec_sf$toID <- toids
  if (total_reversed > 0L) message(sprintf(
    "refactor_flowlines:reversed %d reconciled flowline(s) for upstream-first orientation",
    total_reversed))
  rec_sf
}

# Reference-faithful reconciled-toid derivation. Reconcile builds toids from a
# collapse "representative external toid" heuristic that, at a confluence, can
# pick a co-located sibling instead of the NHD mainstem continuation: severing a
# levelpath that is contiguous in the reference. This
# pass sets each reconciled reach's toid from the reference COMID->toCOMID chain:
# the downstream reach is the one carrying the reference toCOMID of the reach's
# most-downstream (min NHD hydroseq) member COMID. Authoritative over geometry,
# which cannot distinguish mainstem from tributary at a shared node. Guarded by a
# minimal cycle-revert so the (DAG) topology is preserved.
.fix_reconciled_toid_from_reference <- function(rec_sf, ref_toCOMID, ref_hydseq) {
  if (!all(c("ID", "toID", "member_flowpath_id") %in% names(rec_sf))) return(rec_sf)
  rid  <- as.character(rec_sf$ID)
  memb <- rec_sf$member_flowpath_id                       # list of refactored ids
  base <- function(x) format(as.numeric(sub("\\..*$", "", as.character(x))),
    scientific = FALSE, trim = TRUE)   # -> base COMID
  # Everything resolves on the FULL refactored id (split suffix intact). A base
  # COMID is not an address after splitting: it names a reach that may now be
  # several fragments in several groups, so a base-keyed map has to pick one and
  # is wrong for the others.
  # Built vectorised: growing a named vector one member at a time is O(n^2) and
  # costs minutes at VPU scale (~143k members). `fromLast` keeps the
  # last-writer-wins behaviour of the element-by-element assignment it replaces.
  .mem_all <- unlist(lapply(memb, as.character), use.names = FALSE)
  f2r <- stats::setNames(rep(rid, lengths(memb)), .mem_all)  # full refactored id -> rec_id
  f2r <- f2r[!duplicated(names(f2r), fromLast = TRUE)]
  # Parse the split suffix once, with the package's own id parser (it is the one
  # place that knows the encoding, and it handles multi-digit parts like ".10").
  # Unsplit ids come back as part 1, so `is_split` is tracked separately.
  all_ids  <- names(f2r)
  parsed   <- parse_member_part(all_ids)
  base_of  <- stats::setNames(base(all_ids), all_ids)
  part_of  <- stats::setNames(parsed$hf_id_part, all_ids)
  split_of <- stats::setNames(grepl(".", all_ids, fixed = TRUE), all_ids)
  # ENTRY fragment of each reference reach: where flow enters it. That is ".1"
  # for a split reach and the bare id for an unsplit one, the smallest
  # surviving part per base, not whichever member happened to be seen last.
  ent_ord  <- order(base_of, part_of)
  entry_of <- stats::setNames(all_ids[ent_ord], base_of[ent_ord])
  entry_of <- entry_of[!duplicated(names(entry_of))]   # first = smallest part

  orig_toid <- as.numeric(rec_sf$toID)
  new_toid  <- orig_toid
  changed   <- integer(0)
  for (i in seq_along(rid)) {
    fm  <- as.character(memb[[i]])
    cms <- unname(base_of[fm])
    pt  <- unname(part_of[fm])
    ok  <- which(cms %in% names(ref_hydseq))
    if (!length(ok)) next
    # Most-downstream member: min NHD hydroseq, ties broken by the LARGEST split
    # part. Every fragment of one split parent inherits the parent's identical
    # hydroseq, so without the part tie-break the outlet of a group holding two
    # fragments is picked arbitrarily.
    ord    <- ok[order(ref_hydseq[cms[ok]], -pt[ok])]
    j      <- ord[1L]
    outlet <- cms[j]

    # (1) SPLIT-SIBLING CHAIN. split_flowlines cuts a reference reach into
    # ".1" (upstream) ... ".N" (downstream) and chains .k -> .(k+1), only .N
    # inheriting the parent's toid. Resolving a fragment through its BASE COMID
    # throws that away: every fragment then looks up the parent's single
    # reference toCOMID and the upstream fragment is pointed past its own
    # sibling, so the two straddle whichever group holds the reference
    # downstream. Honour the sibling first; it is the reference's own topology
    # at finer grain, not a competing heuristic.
    dnr <- NA_character_
    if (isTRUE(unname(split_of[fm[j]]))) {
      sib <- paste0(sub("\\.[0-9]+$", "", fm[j]), ".", pt[j] + 1L)
      if (sib %in% names(f2r) && f2r[[sib]] != rid[i]) dnr <- f2r[[sib]]
    }

    # (2) Otherwise the reference COMID->toCOMID chain. The reference is
    # base-grained (it predates the split), so its toCOMID names a reach, not a
    # fragment: resolve it to that reach's ENTRY fragment, which is where flow
    # arrives. For an unsplit reach the entry is the reach itself, so this is
    # the ordinary path and is unchanged for the vast majority of reaches.
    if (is.na(dnr)) {
      tc <- ref_toCOMID[outlet]                           # its reference toCOMID
      if (is.na(tc) || tc == "0" || !(tc %in% names(entry_of))) next # terminal / outside VPU -> keep
      dnr <- f2r[[entry_of[[tc]]]]
    }
    if (dnr == rid[i]) next                               # intra-group -> internal, keep
    if (!is.na(new_toid[i]) && as.character(as.integer(new_toid[i])) == dnr) next  # already correct
    new_toid[i] <- as.numeric(dnr)
    changed <- c(changed, i)
  }
  if (!length(changed)) return(rec_sf)

  # Cycle guard: reference is a DAG, but a collapse mapping could close a loop,
  # revert the minimum changed edges until acyclic.
  repeat {
    e <- data.frame(f = rid, t = as.character(as.integer(new_toid)), stringsAsFactors = FALSE)
    e <- e[!is.na(e$t) & e$t %in% rid, ]
    sc <- igraph::components(igraph::graph_from_data_frame(e, directed = TRUE), "strong")
    cyc <- names(sc$membership)[sc$membership %in% which(sc$csize > 1L)]
    rv <- changed[as.character(rid[changed]) %in% cyc]
    if (length(rv) == 0L) break
    j <- rv[1L]
    new_toid[j] <- orig_toid[j]
    changed <- changed[changed != j]
  }
  rec_sf$toID <- new_toid
  message(sprintf(
    "refactor_flowlines: reference-toCOMID derivation corrected %d reconciled toid(s)",
    length(changed)))
  rec_sf
}
