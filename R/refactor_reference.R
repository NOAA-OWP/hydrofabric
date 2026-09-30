#' Refactor a reference hydrofabric (split, collapse, reconcile; refactor phase)
#'
#' Consumes a reference fabric + FDR/FAC rasters; produces refactored
#' flowpaths, reconciled flowpaths + divides (raster-partitioned), and
#' a COMID->reconciled_id lineage table. This is the raster-heavy stage;
#' aggregate phase has no raster dependency.
#'
#' Writes to `{out_dir}/vpu_{vpuid}/work/`:
#' - `refactored.gpkg` (single layer: flowpaths)
#' - `reconciled.gpkg` (layers: flowpaths, divides, events)
#' - `pois.gpkg` (event-bearing POIs, deduplicated)
#' - `hydrolocations.gpkg` (copy of the reference layer for the nextgen phase)
#' - `lineage.parquet` (COMID -> reconciled_id, list cols coerced to CSV)
#' - `timings.json`, `invariants.json`
#' - `provenance.json` (written by the caller: [refactor_fabric()] only
#'   returns the record; writing happens once the caller knows the
#'   image digest + upstream ids).
#'
#' @param vpuid VPU identifier (e.g. `"01"`, `"21"`).
#' @param reference_gpkg Path to the reference fabric. May be multi-VPU
#'   (auto-sliced) or a single-VPU fabric.
#' @param exclude_hl_types Hydrolocation `hl_type` values dropped at ingest;
#'   see [run_pipeline()].
#' @param hydrolocations_gpkg Optional path to a separate hydrolocations
#'   GeoPackage (POIs to embed). When `NULL`, hydrolocations are read from
#'   `reference_gpkg` if present.
#' @param fdr_path Path (or `/vsis3/...`) to the D8 flow-direction raster.
#' @param fac_path Path (or `/vsis3/...`) to the flow-accumulation raster.
#' @param fallback_fdr_path,fallback_fac_path Optional paths to fallback
#'   FDR / FAC rasters used when the primary rasters return NoData. `NULL`
#'   disables the fallback (default); ignored when equal to the primary paths.
#' @param out_dir Parent output directory; `vpu_{vpuid}/` is created.
#' @param tolerance_m POI snap tolerance (m), used in `generate_events`.
#' @param snap_threshold_m Warn when a POI snap distance exceeds this (m).
#' @param split_flines_meters Split flowlines longer than this (m); the default
#'   (1e9) effectively disables Stage-1 splitting (routing-unit splitting is done
#'   later by [add_flowlines()]).
#' @param collapse_flines_meters Minimum inter-confluence segment length (m).
#' @param force_lr Force left/right event assignment.
#' @param n_cores Worker count for `reconcile_divides`.
#' @param use_checkpoints Logical; if `TRUE`, reuse existing stage-1 checkpoints
#'   instead of recomputing.
#' @return Invisibly, a list with `files` (named character, output paths),
#'   `timings` (named list of elapsed seconds), `invariants` (named list).
#' @examples
#' \dontrun{
#' refactor_fabric(vpuid          = "01",
#'                 reference_gpkg = "reference_01.gpkg",
#'                 fdr_path       = "fdr.tif",
#'                 fac_path       = "fac.tif",
#'                 out_dir        = "out",
#'                 collapse_flines_meters = 300)
#' }
#' @export
refactor_fabric <- function(vpuid, reference_gpkg, fdr_path, fac_path,
                            fallback_fdr_path = NULL,
                            fallback_fac_path = NULL,
                            hydrolocations_gpkg = NULL,
                            exclude_hl_types = NULL,
                            out_dir,
                            tolerance_m            = 500,
                            snap_threshold_m       = 250,
                            split_flines_meters    = 1e9,
                            collapse_flines_meters = 250,
                            force_lr               = TRUE,
                            n_cores = NULL,
                            use_checkpoints = FALSE) {
  params <- list(tolerance_m            = tolerance_m,
                 snap_threshold_m       = snap_threshold_m,
                 split_flines_meters    = split_flines_meters,
                 collapse_flines_meters = collapse_flines_meters,
                 force_lr               = force_lr)
  stopifnot(
    is.character(vpuid), length(vpuid) == 1L, nzchar(vpuid),
    file.exists(reference_gpkg),
    .raster_accessible(fdr_path),
    .raster_accessible(fac_path))

  if (is.null(n_cores))
    n_cores <- max(1L, parallel::detectCores() - 2L)

  vpu_dir <- file.path(out_dir, paste0("vpu_", vpuid))
  dir.create(vpu_dir, recursive = TRUE, showWarnings = FALSE)
  work_dir <- file.path(vpu_dir, "work")
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  timings <- list()
  .t <- function(label, expr) {
    t0 <- proc.time()[["elapsed"]]
    result <- force(expr)
    timings[[label]] <<- proc.time()[["elapsed"]] - t0
    message(sprintf("[stage1][%s][%5.1fs] %s",
      vpuid, timings[[label]], label))
    result
  }

  # -- 0. Extract VPU slice if the reference is multi-VPU --------------
  ref_slice <- file.path(work_dir, sprintf("ref_vpu%s.gpkg", vpuid))
  gpkg_in <- .extract_vpu_slice_if_needed(reference_gpkg, ref_slice,
    vpuid, .t_fun = .t, hydrolocations_gpkg = hydrolocations_gpkg,
    exclude_hl_types = exclude_hl_types)

  n_ref_fps <- NA_integer_  # captured in step 1b; used in summary line

  # Checkpoint paths: written immediately after each step so a restart skips
  # completed work. Deleted only when the caller passes overwrite = TRUE.
  ckpt_events     <- file.path(work_dir, "ckpt_events.gpkg")
  ckpt_refactored <- file.path(work_dir, "refactored.gpkg")
  ckpt_reconciled <- file.path(work_dir, "reconciled.gpkg")
  ckpt_rec_div    <- file.path(work_dir, "rec_div.gpkg")

  # When NOT reusing checkpoints, clear any stale stage-1 checkpoint files up front
  # so the writes below don't collide with a prior run's layers ("Layer ... already
  # exists" on a re-run into a non-empty work dir). These are regenerated
  # intermediates, never outputs.
  if (!use_checkpoints) {
    unlink(file.path(work_dir, c("ckpt_events.gpkg", "refactored.gpkg",
      "reconciled.gpkg", "rec_div.gpkg", "lineage.parquet", "rec_fp_agg.gpkg",
      "pois.gpkg")))
  }

  # -- 1. Events ------------------------------------------------------
  if (use_checkpoints && file.exists(ckpt_events)) {
    message(sprintf("[stage1][%s] Skipping Step 1: events checkpoint exists", vpuid))
    events_in <- sf::read_sf(ckpt_events, "events", quiet = TRUE)
  } else {
    events_in <- .t("Step 1: generate_events",
      rf_generate_events_from_hydrolocations_fast(
        hydrolocations_gpkg = gpkg_in, flowpaths_gpkg = gpkg_in,
        hydroloc_layer   = "hydrolocations", flowpath_layer = "flowpaths",
        tolerance_m      = params$tolerance_m,
        snap_threshold_m = params$snap_threshold_m,
        force_lr         = params$force_lr))
    sf::st_write(.strip_list_cols(events_in), ckpt_events, "events", quiet = TRUE)
  }

  # -- 1b. Event integrity gate ---------------------------------------
  # Every event must reference exactly one flowpath that exists in the
  # reference network. Events with NA or unresolved flowpath_id are
  # dropped here rather than silently propagating into refactor_flowlines.
  # The poi/junction designation (event_kind) is reported so any
  # collapse of that distinction is visible before the split step.
  .t("Step 1b: validate_events", {
    fp_ids_ref <- as.character(sf::st_drop_geometry(
      sf::read_sf(gpkg_in, "flowpaths", quiet = TRUE))$flowpath_id)

    n_total <- nrow(events_in)
    ev_fp   <- as.character(events_in$flowpath_id)

    # Mask 1: NA or empty flowpath_id
    na_mask <- is.na(ev_fp) | !nzchar(ev_fp)
    if (any(na_mask)) {
      message(sprintf("[events][%s] dropping %d event(s) with NA/empty flowpath_id",
        vpuid, sum(na_mask)))
      events_in <- events_in[!na_mask, , drop = FALSE]
      ev_fp     <- ev_fp[!na_mask]
    }

    # Mask 2: flowpath_id not in reference network
    bad_mask <- !(ev_fp %in% fp_ids_ref)
    if (any(bad_mask)) {
      message(sprintf("[events][%s] dropping %d event(s) with flowpath_id absent from reference",
        vpuid, sum(bad_mask)))
      events_in <- events_in[!bad_mask, , drop = FALSE]
      ev_fp     <- ev_fp[!bad_mask]
    }

    # Mask 3: remove NBI bridge hydrolocations.
    # Bridges are too dense (1000+ per VPU) and not relevant to NOAA routing.
    # Retaining them creates large numbers of sub-threshold POI-split sub-segments
    # that refactor_flowlines cannot collapse (excl protection), resulting in thousands
    # of connectors and low-coverage divide failures in flat terrain.
    if ("hl_type" %in% names(events_in)) {
      bridge_mask <- !is.na(events_in$hl_type) &
        grepl("^nbi$", as.character(events_in$hl_type),
          ignore.case = TRUE, perl = TRUE)
      if (any(bridge_mask)) {
        message(sprintf("[events][%s] dropping %d NBI bridge event(s) (hl_type=nbi)",
          vpuid, sum(bridge_mask)))
        events_in <- events_in[!bridge_mask, , drop = FALSE]
        ev_fp     <- ev_fp[!bridge_mask]
      }
    } else if ("hl_class" %in% names(events_in)) {
      bridge_mask <- !is.na(events_in$hl_class) &
        grepl("^bridge$", as.character(events_in$hl_class),
          ignore.case = TRUE, perl = TRUE)
      if (any(bridge_mask)) {
        message(sprintf("[events][%s] dropping %d bridge event(s) (hl_class=bridge)",
          vpuid, sum(bridge_mask)))
        events_in <- events_in[!bridge_mask, , drop = FALSE]
        ev_fp     <- ev_fp[!bridge_mask]
      }
    }

    # Check: no hydrolocation (hl_id) assigned to >1 flowpath
    if ("hl_id" %in% names(events_in) && "event_kind" %in% names(events_in)) {
      poi_df  <- sf::st_drop_geometry(
        events_in[!is.na(events_in$event_kind) & events_in$event_kind == "poi",
          c("hl_id", "flowpath_id"), drop = FALSE])
      uniq_df <- unique(poi_df)
      dup_hl  <- unique(uniq_df$hl_id[duplicated(uniq_df$hl_id)])
      if (length(dup_hl) > 0L)
        .hf_log_warn(sprintf("[events][%s] %d hl_id(s) assigned to >1 flowpath",
          vpuid, length(dup_hl)))
    }

    # Report breakdown by event_kind (poi / junction): confirms designation intact
    kind_tbl <- table(events_in$event_kind)
    message(sprintf("[events][%s] %d/%d events valid: %s",
      vpuid, nrow(events_in), n_total,
      paste(names(kind_tbl), as.integer(kind_tbl), sep = "=", collapse = ", ")))

    # Capture reference count here so summary line can use it without re-reading
    n_ref_fps <- length(fp_ids_ref)

    # Overwrite checkpoint with validated events
    sf::st_write(.strip_list_cols(events_in), ckpt_events, "events",
      delete_layer = TRUE, quiet = TRUE)
  })

  # -- 2. Refactor ----------------------------------------------------
  collapse_stuck_warns <- character()
  if (use_checkpoints && file.exists(ckpt_refactored) && file.exists(ckpt_reconciled) &&
    hfutils::layer_exists(ckpt_reconciled, "events")) {
    message(sprintf("[stage1][%s] Skipping Step 2: refactored/reconciled checkpoints exist", vpuid))
    # Layer name: current fast_reactor.R uses filename-derived names ("refactored",
    # "reconciled"); older checkpoints used "flowpaths". Try both.
    .read_ckpt_layer <- function(gpkg, preferred, fallback = "flowpaths") {
      lyr <- if (hfutils::layer_exists(gpkg, preferred)) preferred else fallback
      sf::read_sf(gpkg, lyr, quiet = TRUE)
    }
    out <- list(
      refactored = .read_ckpt_layer(ckpt_refactored, "refactored"),
      reconciled = .read_ckpt_layer(ckpt_reconciled, "reconciled"),
      events     = sf::read_sf(ckpt_reconciled, "events",      quiet = TRUE))
  } else {
    out <- .t("Step 2: refactor_flowlines",
      withCallingHandlers(
        refactor_flowlines(
          flines                 = sf::read_sf(gpkg_in, "flowpaths", quiet = TRUE),
          split_flines_meters    = params$split_flines_meters,
          collapse_flines_meters = params$collapse_flines_meters,
          out_refactored         = ckpt_refactored,
          out_reconciled         = ckpt_reconciled,
          exclude_cats           = NULL,
          events                 = events_in[events_in$event_kind == "poi", ],
          cores                  = 1L),
        warning = function(w) {
          if (grepl("Phase [0-9]+ stuck", conditionMessage(w))) {
            collapse_stuck_warns <<- c(collapse_stuck_warns, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        }
    ))
    if (length(collapse_stuck_warns) > 0L)
      .hf_log_warn(sprintf("[stage1][%s] %d collapse convergence failure(s): %s",
        vpuid, length(collapse_stuck_warns),
        paste(collapse_stuck_warns, collapse = "; ")))
  }

  invariants <- list()
  if (length(collapse_stuck_warns) > 0L)
    invariants$collapse_stuck <- collapse_stuck_warns

  invariants$refactored <- tryCatch(
    check_invariants("refactored",
      refactored = out$refactored, reconciled = out$reconciled),
    error = function(e) list(error = conditionMessage(e)))

  # -- 2b. DAG assertion + hydroseq recompute on reconciled network -----
  # Cycles here indicate a bug in fast_collapse/reconcile_divides toid logic.
  # The break is an emergency fallback only; if it fires, investigate upstream.
  .t("Step 2b: dag_assert_reconciled", {
    dag_v <- hfutils::hf_assert_network_dag(out$reconciled,
      id_col = "ID", toid_col = "toID")
    if (!dag_v$is_dag) {
      .hf_log_warn(sprintf(
        "[reconcile][%s] BUG: cycles in reconciled network (%s) -- collapse toid logic is incorrect; breaking as fallback",
        vpuid, dag_v$message))
      out$reconciled <- hfutils::hf_break_cycles(out$reconciled,
        id_col = "ID", toid_col = "toID")
      invariants$reconciled_dag_fixed <- dag_v$cycle_ids
    } else {
      message(sprintf("[reconcile][%s] DAG OK", vpuid))
    }
    out$reconciled <- hfutils::hf_recompute_hydroseq(out$reconciled,
      id_col = "ID", toid_col = "toID")
    message(sprintf("[reconcile][%s] hydroseq recomputed (%d fps)", vpuid, nrow(out$reconciled)))
  })

  ckpt_lineage    <- file.path(work_dir, "lineage.parquet")
  ckpt_rec_fp_agg <- file.path(work_dir, "rec_fp_agg.gpkg")
  ckpt_pois       <- file.path(work_dir, "pois.gpkg")

  if (use_checkpoints && file.exists(ckpt_lineage)) {
    lineage <- NULL  # not needed downstream; file will be used directly
  } else {
    lineage <- lineage_build(out$refactored, out$reconciled)
  }

  # Prep the aggregation-facing flowpath table (columns renamed for
  # downstream aggregation): same transforms as run_pipeline.R:104-122.
  if (use_checkpoints && file.exists(ckpt_rec_fp_agg)) {
    rec_fp_agg <- sf::read_sf(ckpt_rec_fp_agg, "flowpaths", quiet = TRUE)
  } else {
    rec_fp_agg <- out$reconciled
    names(rec_fp_agg)[names(rec_fp_agg) == "ID"]   <- "flowpath_id"
    names(rec_fp_agg)[names(rec_fp_agg) == "toID"] <- "flowpath_toid"
    rec_fp_agg$flowpath_id   <- as.character(rec_fp_agg$flowpath_id)
    rec_fp_agg$flowpath_toid <- as.character(rec_fp_agg$flowpath_toid)
    if (is.list(rec_fp_agg$member_flowpath_id)) {
      rec_fp_agg$member_comid <- vapply(rec_fp_agg$member_flowpath_id,
        paste, collapse = ",", character(1))
    } else if (is.character(rec_fp_agg$member_flowpath_id)) {
      rec_fp_agg$member_comid <- rec_fp_agg$member_flowpath_id
    }
    rec_fp_agg$member_flowpath_id <- NULL
    rec_fp_agg$flowline_id <- rec_fp_agg$flowpath_id
  }

  if (use_checkpoints && file.exists(ckpt_pois)) {
    poi_ev2 <- sf::read_sf(ckpt_pois, "pois", quiet = TRUE)
    poi_ev2_dedup <- poi_ev2
  } else {
    poi_ev2 <- events_in[events_in$event_kind == "poi", ]
    # Ensure poi_id exists: required by aggregate_to_distribution. Reference
    # fabrics without it (e.g. AK TDX) use hl_id as a stable surrogate. A
    # present-but-all-NA poi_id column (AK reference_fabric_NWMv4 ships it empty)
    # is treated identically to a missing one.
    if (!"poi_id" %in% names(poi_ev2) || all(is.na(poi_ev2$poi_id))) {
      if ("hl_id" %in% names(poi_ev2) && !all(is.na(poi_ev2$hl_id))) {
        poi_ev2$poi_id <- as.character(poi_ev2$hl_id)
      } else {
        poi_ev2$poi_id <- as.character(seq_len(nrow(poi_ev2)))
      }
    }
    # Match reconciled_id back to poi_ev2. Prefer poi_id join (CONUS-style events
    # that carry a poi_id); fall back to orig_flowpath_id / flowpath_id join for
    # reference fabrics (e.g. AK TDX) that don't generate poi_id in events.
    if (!is.null(out$events) && "poi_id" %in% names(out$events) &&
      !all(is.na(out$events$poi_id)) && "poi_id" %in% names(poi_ev2)) {
      poi_ev2$flowpath_id <- as.character(
        out$events$reconciled_id[match(poi_ev2$poi_id, out$events$poi_id)])
    } else if (!is.null(out$events) && "orig_flowpath_id" %in% names(out$events)) {
      poi_ev2$flowpath_id <- as.character(
        out$events$reconciled_id[match(as.character(poi_ev2$flowpath_id),
          out$events$orig_flowpath_id)])
    } else {
      # length-safe: a domain with no POI events (e.g. Guam) has a 0-row poi_ev2,
      # where a scalar assignment would error ("replacement has 1 row, data has 0").
      poi_ev2$flowpath_id <- rep(NA_character_, nrow(poi_ev2))
    }

    # Put the POI exactly on the node it is attributed to.
    #
    # Up to here the geometry is still the raw hydrolocation coordinate, so a POI
    # that split a reach sits near the resulting boundary rather than on it: the
    # measure and the split are computed on different vertex sets for the same
    # reach, and the same along-line distance lands in two slightly different
    # places.
    #
    # Taking the location from the produced network instead of from a measure
    # avoids that, the POI is placed at the outlet vertex of the reconciled
    # reach it is attributed to, so coincidence holds by construction. The
    # original coordinate is kept in the hydrolocations layer as hl_x / hl_y.
    if (nrow(poi_ev2) > 0L && !is.null(out$reconciled) &&
      inherits(poi_ev2, "sf") && inherits(out$reconciled, "sf")) {
      rid <- as.character(if ("ID" %in% names(out$reconciled))
        out$reconciled$ID else out$reconciled$rec_id)
      ix  <- match(as.character(poi_ev2$flowpath_id), rid)
      hit <- which(!is.na(ix))
      if (length(hit)) {
        rg <- sf::st_geometry(out$reconciled)
        pg <- sf::st_geometry(poi_ev2)
        if (sf::st_crs(pg) != sf::st_crs(rg))
          pg <- sf::st_transform(pg, sf::st_crs(rg))
        # Bounded by the placement rule's own 200 m threshold.
        #
        # Unbounded, this relocated 56 POIs by more than 200 m and 6 by more than
        # 5 km (max 8847 m). Those are not precision corrections: a POI that far
        # from the outlet it is attributed to is MIS-ATTRIBUTED, and moving it
        # there converts a visible indexing error into an invisible relocation,
        # and makes poi_at_flowline_outlet report a perfect 100% for the wrong
        # reason, since placing POIs on outlets guarantees the check passes.
        #
        # The split is diagnostic: of the 187 moves over 30 m, 86 are nwis gages
        # and only 3 exceed 200 m, so the sub-200 m population is the genuine
        # geometry-basis offset this fix exists for. Of the 56 over 200 m, 41 are
        # wbout: waterbody outlets, the class that snaps poorly because it sits
        # on zero-area connectors. Those are left exactly where they are so they
        # keep failing the outlet checks and stay visible as indexing defects.
        max_move_m <- 200
        moved <- numeric(0); skipped <- numeric(0)
        newg <- as.list(pg)
        for (k in hit) {
          m <- sf::st_coordinates(rg[[ix[k]]])
          if (is.null(dim(m)) || nrow(m) < 1L) next
          tgt <- unname(m[nrow(m), 1:2])
          d_k <- sqrt(sum((unlist(newg[[k]])[1:2] - tgt)^2))
          if (!is.finite(d_k) || d_k > max_move_m) {
            skipped <- c(skipped, d_k)
            next
          }
          moved <- c(moved, d_k)
          newg[[k]] <- sf::st_point(tgt)
        }
        # Rebuild the column rather than assigning elements in place: in-place sfc
        # element assignment leaves the class stale and breaks a later transform.
        sf::st_geometry(poi_ev2) <- sf::st_sfc(newg, crs = sf::st_crs(rg))
        message(sprintf(paste0(
          "[stage1][%s] POIs placed on their reconciled outlet: %d moved ",
          "(median %.1fm, max %.1fm); %d left in place as > %dm from their ",
          "attributed outlet (max %.0fm) -- mis-attribution, not precision"),
          vpuid, sum(moved > 1),
          if (length(moved)) stats::median(moved) else 0,
          if (length(moved)) max(moved) else 0,
          length(skipped), max_move_m,
          if (length(skipped)) max(skipped, na.rm = TRUE) else 0))
      }
    }

    # Merge hydrolocations sharing the same poi_id (e.g. a coastal outlet
    # and a co-located nwis gauge both assigned the same poi_id in the
    # reference fabric).  Prefer the nwis row as the anchor (gauge snap
    # position); combine hl_type / hl_reference / hl_link across the group.
    dup_ids <- unique(poi_ev2$poi_id[duplicated(poi_ev2$poi_id)])
    dup_ids <- dup_ids[!is.na(dup_ids)]   # NA poi_id is not a real shared-id group
    if (length(dup_ids)) {
      .combine_pipe <- function(vals) {
        parts <- unique(trimws(unlist(strsplit(as.character(vals[!is.na(vals)]), "|", fixed = TRUE))))
        paste(sort(parts[nzchar(parts)]), collapse = "|")
      }
      for (pid in dup_ids) {
        rows <- which(poi_ev2$poi_id == pid)
        nwis_row <- rows[grepl("nwis", poi_ev2$hl_type[rows], ignore.case = TRUE)]
        keep_row  <- if (length(nwis_row) > 0L) nwis_row[1L] else rows[1L]
        # hl_id is combined alongside the descriptive columns. Without it the
        # merged row kept both hl_references but only one hl_id, so the surviving
        # POI could not be joined back to the second source hydrolocation by id,
        # the lineage this merge exists to preserve. Measured before the fix: 801
        # minted hl_ids reached the product as 790, while only 2 hl_references
        # were missing, i.e. 9 hydrolocations survived by reference but lost their
        # id. Coerced to character first so the column type change is explicit.
        if ("hl_id" %in% names(poi_ev2))
          poi_ev2$hl_id <- as.character(poi_ev2$hl_id)
        for (col in c("hl_id", "hl_type", "hl_reference", "hl_link")) {
          if (col %in% names(poi_ev2))
            poi_ev2[keep_row, col] <- .combine_pipe(poi_ev2[[col]][rows])
        }
        poi_ev2 <- poi_ev2[-setdiff(rows, keep_row), ]
      }
      message(sprintf("stage1: merged %d shared-poi_id hydrolocation pair(s) -> hl_type combined",
        length(dup_ids)))
    }
    # Collocated hydrolocations COLLAPSE into a single POI: one poi_id, but both
    # hydrolocations preserved in the merged row's lineage. The shared-poi_id merge
    # above already does this; this applies the same rule when two hydrolocations
    # resolve to one reconciled reach carrying DIFFERENT poi_ids.
    #
    # Merge rather than drop by row order, which would let a coastal marker outrank
    # a gage purely by sorting first. Anchor on the gage row, as the shared-poi_id
    # merge does, a gage is what the fabric is reported at, so it must not lose
    # an id contest to a marker.
    fp_dups <- unique(poi_ev2$flowpath_id[
      duplicated(poi_ev2$flowpath_id) & !is.na(poi_ev2$flowpath_id)])
    if (length(fp_dups)) {
      .combine_pipe2 <- function(vals) {
        parts <- unique(trimws(unlist(strsplit(
          as.character(vals[!is.na(vals)]), "|", fixed = TRUE))))
        paste(sort(parts[nzchar(parts)]), collapse = "|")
      }
      drop_rows <- integer(0)
      for (fpid in fp_dups) {
        rows <- which(!is.na(poi_ev2$flowpath_id) & poi_ev2$flowpath_id == fpid)
        if (length(rows) < 2L) next
        gage <- rows[grepl("nwis|rfc", poi_ev2$hl_type[rows], ignore.case = TRUE)]
        keep <- if (length(gage)) gage[1L] else rows[1L]
        if ("hl_id" %in% names(poi_ev2))
          poi_ev2$hl_id <- as.character(poi_ev2$hl_id)
        for (col in c("hl_id", "hl_type", "hl_reference", "hl_link", "hl_uri")) {
          if (col %in% names(poi_ev2))
            poi_ev2[keep, col] <- .combine_pipe2(poi_ev2[[col]][rows])
        }
        drop_rows <- c(drop_rows, setdiff(rows, keep))
      }
      message(sprintf(paste0(
        "stage1: collapsed %d hydrolocation(s) into %d shared POI(s) on one reach ",
        "(identities combined, none dropped)"),
        length(drop_rows), length(fp_dups)))
      if (length(drop_rows)) poi_ev2 <- poi_ev2[-drop_rows, ]
    }
    poi_ev2_dedup <- poi_ev2

    # POI coverage invariant.
    # n_poi_raw:      total POI events in ckpt_events (includes shared-poi_id duplicates)
    # n_truly_lost:   POIs with NA flowpath_id after reconciled_id lookup, the only
    #                 actionable case; warn on these
    # n_poi_merged:   rows removed by shared-poi_id consolidation (legitimate)
    # n_poi_deduped:  rows removed because >1 POI mapped to the same reconciled group
    # n_poi_out:      unique reconciled groups with a POI (what goes in pois.gpkg)
    n_poi_raw    <- sum(events_in$event_kind == "poi", na.rm = TRUE)
    n_truly_lost <- sum(is.na(poi_ev2_dedup$flowpath_id))
    n_poi_out    <- sum(!is.na(poi_ev2_dedup$flowpath_id))
    n_poi_merged <- n_poi_raw - nrow(poi_ev2)
    n_poi_deduped <- nrow(poi_ev2) - nrow(poi_ev2_dedup)
    if (n_truly_lost > 0L) {
      .hf_log_warn(sprintf(paste0(
        "[stage1][%s] %d of %d POI event(s) have no reconciled_id ",
        "(no member match and >%gm from any endpoint)"),
        vpuid, n_truly_lost, n_poi_raw, 250))
      invariants$unsnapped_pois <- n_truly_lost
    } else {
      message(sprintf(paste0(
        "[stage1][%s] all POI events resolved: %d -> %d unique reconciled ",
        "groups (%d shared-poi merged, %d fp-deduped)"),
        vpuid, n_poi_raw, n_poi_out, n_poi_merged, n_poi_deduped))
    }
  }  # end else (pois checkpoint not found)

  # -- 3. Reconcile (raster partitioning) -----------------------------
  if (use_checkpoints && file.exists(ckpt_rec_div)) {
    message(sprintf("[stage1][%s] Skipping Step 3: rec_div checkpoint exists", vpuid))
    rec_div <- sf::read_sf(ckpt_rec_div, "divides", quiet = TRUE)
  } else {
    rec_div <- .t("Step 3: reconcile_divides",
      reconcile_divides(
        flowline          = out$refactored,
        divides           = sf::read_sf(gpkg_in, "divides", quiet = TRUE),
        fdr_path          = fdr_path,
        fac_path          = fac_path,
        fallback_fdr_path = fallback_fdr_path,
        fallback_fac_path = fallback_fac_path,
        lookup_table      = out$reconciled,
        lr                = FALSE,
        parallel          = TRUE,
        num_workers       = n_cores))
    sf::st_write(.strip_list_cols(rec_div), ckpt_rec_div, "divides", quiet = TRUE)
  }

  rec_div$divide_id   <- as.character(rec_div$divide_id)
  rec_div$flowpath_id <- as.character(rec_div$flowpath_id)

  # -- 3b. Post-reconcile integrity gate ------------------------------
  .t("Step 3b: validate_reconciled", {
    rec_fp <- out$reconciled

    # 1. Dotted IDs are acceptable if both sides carry them.
    #    Every fp with a dotted ID must have a matching divide with the same ID.
    fp_ids  <- as.character(rec_fp$ID)
    div_ids <- as.character(rec_div$divide_id)
    dotted  <- grepl(".", fp_ids, fixed = TRUE)
    if (any(dotted)) {
      dotted_fp  <- fp_ids[dotted]
      unmatched  <- dotted_fp[!dotted_fp %in% div_ids]
      if (length(unmatched) > 0L)
        .hf_log_warn(sprintf(
          "[reconcile][%s] %d dotted-ID flowpath(s) have no matching divide (IDs: %s ...)",
          vpuid, length(unmatched),
          paste(utils::head(unmatched, 5), collapse = ", ")))
      else
        message(sprintf("[reconcile][%s] %d dotted-ID flowpath(s) -- all matched by divide",
          vpuid, sum(dotted)))
    } else {
      message(sprintf("[reconcile][%s] all %d flowpath IDs are clean integers", vpuid, length(fp_ids)))
    }

    # 2. Geometry types.
    #    Reconciled flowpaths are MULTILINESTRING (collapsed groups): correct.
    #    Refactored (split) flowpaths must be LINESTRING.
    #    Reconciled divides must be POLYGON (no MULTIPOLYGON).
    ref_fac_lyr <- if (hfutils::layer_exists(ckpt_refactored, "refactored")) "refactored" else "flowpaths"
    ref_fac <- sf::read_sf(ckpt_refactored, ref_fac_lyr, quiet = TRUE)
    rfac_types  <- as.character(sf::st_geometry_type(ref_fac))
    bad_rfac_gt <- rfac_types != "LINESTRING"
    if (any(bad_rfac_gt))
      .hf_log_warn(sprintf("[reconcile][%s] %d refactored flowpath(s) are not LINESTRING: %s",
        vpuid, sum(bad_rfac_gt),
        paste(unique(rfac_types[bad_rfac_gt]), collapse = ", ")))
    else
      message(sprintf("[reconcile][%s] refactored geometry OK: %d LINESTRING", vpuid, nrow(ref_fac)))

    div_types  <- as.character(sf::st_geometry_type(rec_div))
    bad_div_gt <- div_types != "POLYGON"
    n_multi    <- sum(div_types == "MULTIPOLYGON")
    if (any(bad_div_gt)) {
      .hf_log_warn(sprintf("[reconcile][%s] %d divide(s) are not POLYGON: %s",
        vpuid, sum(bad_div_gt),
        paste(unique(div_types[bad_div_gt]), collapse = ", ")))
      if (n_multi > 0L)
        invariants$rec_div_multipolygon <- n_multi
    } else {
      message(sprintf("[reconcile][%s] divide geometry OK: %d POLYGON", vpuid, nrow(rec_div)))
    }

    # 3. 1:1 flowpath/divide: count orphans on each side.
    #    Match by divide_id == fp ID (shared integer / dotted convention).
    fp_id_set  <- unique(fp_ids)
    div_id_set <- unique(div_ids)
    fp_no_div  <- setdiff(fp_id_set,  div_id_set)   # flowpaths without a divide
    div_no_fp  <- setdiff(div_id_set, fp_id_set)    # divides without a flowpath

    # Reference orphan count (pre-reconcile) for comparison.
    # "divides" may be absent from a per-VPU presplit slice (geometry column
    # named "geom" confuses sf::st_read on some GPKG writers). Fall back to 0
    # so the comparison degrades gracefully rather than aborting.
    ref_fp_ids  <- as.character(sf::st_drop_geometry(
      sf::read_sf(gpkg_in, "flowpaths", quiet = TRUE))$flowpath_id)
    ref_div_raw <- tryCatch(
      sf::st_drop_geometry(sf::read_sf(gpkg_in, "divides", quiet = TRUE)),
      error = function(e) NULL)
    ref_div_ids <- if (!is.null(ref_div_raw)) as.character(ref_div_raw$flowpath_id) else character(0L)
    ref_no_div  <- length(setdiff(ref_fp_ids, ref_div_ids))
    ref_no_fp   <- length(setdiff(ref_div_ids, ref_fp_ids))

    message(sprintf("[reconcile][%s] fp-without-divide: %d (ref had %d)",
      vpuid, length(fp_no_div), ref_no_div))

    # 3a. Lineage check, every fp-without-divide must trace back exclusively
    #     to reference flowlines that themselves had no divide (has_divide = FALSE
    #     or divide_id = NA). Any fp-without-divide whose members had divides in
    #     the reference signals a reconcile error, not a legitimate connector.
    #     Only warn after the lineage check: if all orphans trace to no-divide
    #     reference features the count increase is expected (new connectors created
    #     by reconciliation) and should not surface as a warning.
    if (length(fp_no_div) > 0L && "member_flowpath_id" %in% names(rec_fp)) {
      ref_fp_raw <- sf::st_drop_geometry(sf::read_sf(gpkg_in, "flowpaths", quiet = TRUE))
      ref_fp_tbl <- ref_fp_raw[, intersect(c("flowpath_id", "has_divide", "divide_id"), names(ref_fp_raw)), drop = FALSE]
      if (!"has_divide" %in% names(ref_fp_tbl)) ref_fp_tbl$has_divide <- FALSE
      if (!"divide_id"  %in% names(ref_fp_tbl)) ref_fp_tbl$divide_id  <- NA_character_
      ref_fp_tbl$flowpath_id <- as.character(ref_fp_tbl$flowpath_id)

      # Cross-check: only treat a reference fp as having a divide when the polygon
      # actually exists in the divides layer (mirrors the predictive block's cross-check).
      ref_div_fps_cc <- tryCatch(
        as.character(sf::st_drop_geometry(
          sf::read_sf(gpkg_in, "divides", quiet = TRUE))$flowpath_id),
        error = function(e) NULL)
      if (!is.null(ref_div_fps_cc)) {
        ref_fp_tbl$has_divide <- ref_fp_tbl$has_divide &
          ref_fp_tbl$flowpath_id %in% ref_div_fps_cc
        ref_fp_tbl$divide_id[!ref_fp_tbl$flowpath_id %in% ref_div_fps_cc] <- NA_character_
      }

      orphan_rows <- rec_fp[as.character(rec_fp$ID) %in% fp_no_div, , drop = FALSE]
      members_csv <- sf::st_drop_geometry(orphan_rows)$member_flowpath_id

      bad_lineage <- vapply(members_csv, function(csv) {
        raw <- trimws(unlist(strsplit(as.character(csv), ",")))
        # Strip ".N" split-child suffixes so "12345.1" resolves to parent "12345"
        ids <- unique(sub("\\.[0-9]+$", "", raw))
        m_ref <- ref_fp_tbl[ref_fp_tbl$flowpath_id %in% ids, , drop = FALSE]
        if (nrow(m_ref) == 0L) return(FALSE)
        any(m_ref$has_divide, na.rm = TRUE) || any(!is.na(m_ref$divide_id))
      }, logical(1L))

      n_bad <- sum(bad_lineage)
      if (n_bad > 0L)
        .hf_log_warn(sprintf(
          "[reconcile][%s] %d fp-without-divide have member(s) that HAD divides in reference -- unexpected connector(s)",
          vpuid, n_bad))
      else
        message(sprintf(
          "[reconcile][%s] lineage OK: all %d fp-without-divide trace to reference features with no divide",
          vpuid, length(fp_no_div)))
    } else if (length(fp_no_div) > ref_no_div) {
      # No member_flowpath_id: cannot do lineage check; fall back to count comparison
      .hf_log_warn(sprintf(
        "[reconcile][%s] fp-without-divide: %d (was %d in reference -- %d new orphans; no lineage column for check)",
        vpuid, length(fp_no_div), ref_no_div, length(fp_no_div) - ref_no_div))
    }

    if (length(div_no_fp) > ref_no_fp)
      .hf_log_warn(sprintf(
        "[reconcile][%s] divide-without-fp: %d (was %d in reference -- %d new orphans)",
        vpuid, length(div_no_fp), ref_no_fp, length(div_no_fp) - ref_no_fp))
    else
      message(sprintf("[reconcile][%s] divide-without-fp: %d (ref had %d -- OK)",
        vpuid, length(div_no_fp), ref_no_fp))

    # 4. Events at outlets: use reconciled_id from out$events (post-split IDs).
    #    events_in carries pre-split reference flowpath_ids; after refactor_flowlines
    #    those map to new reconciled IDs stored in out$events$reconciled_id.
    rec_ev <- if (!is.null(out$events) && "reconciled_id" %in% names(out$events))
      sf::st_drop_geometry(out$events)
    else
      NULL

    if (!is.null(rec_ev) && nrow(rec_ev) > 0L) {
      poi_rec_ids <- as.character(rec_ev$reconciled_id)
      in_rec      <- poi_rec_ids %in% fp_id_set
      if (!all(in_rec, na.rm = TRUE))
        .hf_log_warn(sprintf("[reconcile][%s] %d POI event(s) have reconciled_id not found in reconciled flowpaths",
          vpuid, sum(!in_rec, na.rm = TRUE)))

      snap_tol <- params$collapse_flines_meters
      rec_fp_t <- sf::st_transform(rec_fp, sf::st_crs(out$events))
      ev_geom  <- sf::st_geometry(out$events)[in_rec]
      fp_row   <- match(poi_rec_ids[in_rec], as.character(rec_fp_t$ID))
      valid    <- !is.na(fp_row)
      if (any(valid)) {
        outlets <- sf::st_sfc(lapply(fp_row[valid], function(i) {
          cc <- sf::st_coordinates(rec_fp_t[i, ])
          sf::st_point(cc[nrow(cc), c("X", "Y")])
        }), crs = sf::st_crs(out$events))
        d <- suppressWarnings(as.numeric(
          sf::st_distance(ev_geom[valid], outlets, by_element = TRUE)))
        far        <- d > snap_tol
        snap_types <- out$events$snap_type[in_rec][valid]
        far_member <- far & !is.na(snap_types) & snap_types == "member"
        far_other  <- far & (is.na(snap_types) | snap_types != "member")

        if (any(far_member, na.rm = TRUE))
          message(sprintf(
            "[reconcile][%s] %d POI event(s) >%gm from group outlet (snap_type=member -- collapsed reach, expected)",
            vpuid, sum(far_member, na.rm = TRUE), snap_tol))
        if (any(far_other, na.rm = TRUE))
          .hf_log_warn(sprintf(
            "[reconcile][%s] %d non-member POI event(s) are >%gm from their reconciled flowpath outlet",
            vpuid, sum(far_other, na.rm = TRUE), snap_tol))
        if (!any(far_other, na.rm = TRUE))
          message(sprintf(paste0(
            "[reconcile][%s] all %d member-resolved POI events are at or ",
            "near their reconciled flowpath outlet"),
            vpuid, sum(valid) - sum(far_member, na.rm = TRUE)))
      }
    }
  })

  # Tag each reconciled flowpath with has_divide = TRUE if any reference member
  # had a divide, but only if that divide polygon actually exists in the reference
  # data. The reference flowpaths layer may carry has_divide=TRUE for IDs whose
  # polygon was never written to the divides layer (e.g. TDX-Hydro phantom divides),
  # which would produce spurious data-loss reports in the block below.
  if (!is.null(out$reconciled) && "member_flowpath_id" %in% names(out$reconciled)) {
    ref_hd <- tryCatch(
      {
        tmp      <- sf::st_drop_geometry(sf::read_sf(gpkg_in, "flowpaths", quiet = TRUE))
        hd_flag  <- as.logical(tmp$has_divide)
        # If reference has no has_divide column (e.g. NHD-HR for HI/PRVI), assume
        # complete divide coverage and derive purely from the reference divides layer.
        if (length(hd_flag) == 0L) hd_flag <- rep(TRUE, nrow(tmp))
        # Cross-check: only TRUE where the divide polygon actually exists in reference.
        ref_div_fps <- tryCatch(
          as.character(sf::st_drop_geometry(
            sf::read_sf(gpkg_in, "divides", quiet = TRUE))$flowpath_id),
          error = function(e) NULL)
        if (!is.null(ref_div_fps))
          hd_flag <- hd_flag & as.character(tmp$flowpath_id) %in% ref_div_fps
        stats::setNames(hd_flag, as.character(tmp$flowpath_id))
      },
      error = function(e) NULL)
    if (!is.null(ref_hd)) {
      mems <- as.character(sf::st_drop_geometry(out$reconciled)$member_flowpath_id)
      out$reconciled$has_divide <- vapply(mems, function(csv) {
        ids <- unique(sub("[.][0-9]+$", "", trimws(unlist(strsplit(csv, ",", fixed = TRUE)))))
        any(ref_hd[ids], na.rm = TRUE)
      }, logical(1L))
    }
  }

  # Reconcile predictive has_divide against actual Step 3 output.
  # Lost paths (predicted TRUE by reference, no polygon produced) represent
  # genuine algorithmic data loss (the reference divide existed but raster
  # partitioning gave all the area to sibling paths). Log with lineage and trace
  # where the area went via rec_div$orig_divide_id.
  if (!is.null(rec_div) && !is.null(out$reconciled) && "has_divide" %in% names(out$reconciled)) {
    actual_ids <- as.character(rec_div$flowpath_id)
    rec_ids    <- as.character(out$reconciled$ID)
    rfc_ids    <- as.character(out$refactored$flowpath_id)
    lost_mask  <- out$reconciled$has_divide & !rec_ids %in% actual_ids

    if (any(lost_mask)) {
      lost_df <- sf::st_drop_geometry(out$reconciled)[lost_mask,
        intersect(c("ID", "member_flowpath_id"), names(out$reconciled)),
        drop = FALSE]
      message(sprintf(
        "[stage1][%s] DATA LOSS: %d flowpath(s) predicted has_divide=TRUE (reference) but received no divide polygon from Step 3:",
        vpuid, nrow(lost_df)))
      for (i in seq_len(nrow(lost_df))) {
        mems  <- trimws(unlist(strsplit(as.character(lost_df$member_flowpath_id[i]), ",", fixed = TRUE)))
        refs  <- unique(sub("[.][0-9]+$", "", mems))
        n_rfc <- sum(mems %in% rfc_ids)
        message(sprintf("  reconciled_id=%-20s  members_in_refactored=%d/%d  reference_comids=%s",
          lost_df$ID[i], n_rfc, length(mems), paste(refs, collapse = ", ")))
      }

      # Track where each lost divide area ended up using orig_divide_id lineage.
      # orig_divide_id records the majority reference COMID for each rec_div polygon;
      # when a lost path's reference COMID appears there, it identifies which sibling
      # reconciled path absorbed the area during raster partitioning.
      if ("orig_divide_id" %in% names(rec_div)) {
        oid_by_ref <- split(as.character(rec_div$flowpath_id),
          as.character(rec_div$orig_divide_id))
        for (i in seq_len(nrow(lost_df))) {
          refs     <- unique(sub("[.][0-9]+$", "",
            trimws(unlist(strsplit(as.character(lost_df$member_flowpath_id[i]),
              ",", fixed = TRUE)))))
          absorbed <- sort(unique(unlist(oid_by_ref[refs[refs %in% names(oid_by_ref)]])))
          if (length(absorbed) > 0L) {
            message(sprintf("    [trace] reconciled_id=%-20s -> area absorbed by: %s",
              lost_df$ID[i], paste(absorbed, collapse = ", ")))
          } else {
            message(sprintf("    [trace] reconciled_id=%-20s -> area not traceable via orig_divide_id",
              lost_df$ID[i]))
          }
        }
      }

      # Repair: merge each divideless path into its immediate downstream path.
      # The downstream path's divide already covers this path's watershed (the
      # raster gave it all the area), so the fix is topological only: absorb
      # the flowline geometry and members into toID, rewire upstream, drop the
      # divideless path.
      rec_df    <- sf::st_drop_geometry(out$reconciled)
      geom_list <- as.list(sf::st_geometry(out$reconciled))
      id_chr    <- as.character(rec_df$ID)
      to_chr    <- as.character(rec_df$toID)
      drop_rows <- integer(0)
      # Where each dropped path's identity goes, so anything keyed on it can follow.
      remap     <- character(0)

      for (lid in as.character(lost_df$ID)) {
        src <- match(lid, id_chr)
        if (is.na(src)) next
        tgt <- match(to_chr[src], id_chr)
        if (is.na(tgt)) {
          message(sprintf("    [repair] reconciled_id=%-20s  downstream not in reconciled -- skipped", lid))
          next
        }

        # Union flowline geometry
        geom_list[[tgt]] <- sf::st_union(
          sf::st_sfc(geom_list[[src]], crs = sf::st_crs(out$reconciled)),
          sf::st_sfc(geom_list[[tgt]], crs = sf::st_crs(out$reconciled)))[[1L]]

        # Absorb member_flowpath_id
        src_m <- trimws(unlist(strsplit(as.character(rec_df$member_flowpath_id[src]), ",", fixed = TRUE)))
        tgt_m <- trimws(unlist(strsplit(as.character(rec_df$member_flowpath_id[tgt]), ",", fixed = TRUE)))
        rec_df$member_flowpath_id[tgt] <- paste(unique(c(tgt_m, src_m)), collapse = ",")

        # Rewire: anything upstream of the lost path now points to the target
        rewire <- !is.na(to_chr) & to_chr == id_chr[src]
        to_chr[rewire] <- id_chr[tgt]

        message(sprintf("    [repair] merged reconciled_id=%-20s into %s  (rewired %d upstream path(s))",
          lid, id_chr[tgt], sum(rewire)))
        remap[id_chr[src]] <- id_chr[tgt]
        drop_rows <- c(drop_rows, src)
      }

      # A POI attributed to a dropped path has to move with it.
      #
      # The repair rewires upstream toIDs onto the absorbing path. The POI layer is
      # keyed on the same reconciled ids, so without this it points at ids that no
      # longer exist and stage 2 drops those POIs when it joins by flowpath_id.
      #
      # The absorbing path is the correct destination, not a fallback: the raster
      # partition already gave this path's watershed to it, so that is where the
      # POI's drainage now is.
      if (length(remap)) {
        # Resolve chains (A absorbed into B, B later absorbed into C).
        for (. in seq_len(10L)) {
          nxt <- remap[unname(remap)]
          upd <- !is.na(nxt)
          if (!any(upd)) break
          remap[upd] <- unname(nxt[upd])
        }
        .repoint <- function(x) {
          if (is.null(x) || !"flowpath_id" %in% names(x)) return(x)
          v   <- as.character(x$flowpath_id)
          hit <- !is.na(v) & v %in% names(remap)
          if (any(hit)) x$flowpath_id[hit] <- unname(remap[v[hit]])
          x
        }
        n_moved <- sum(as.character(poi_ev2_dedup$flowpath_id) %in% names(remap),
          na.rm = TRUE)
        poi_ev2       <- .repoint(poi_ev2)
        poi_ev2_dedup <- .repoint(poi_ev2_dedup)
        if (n_moved > 0L)
          message(sprintf(paste0(
            "[stage1][%s] re-pointed %d POI(s) from merged divide-less path(s) ",
            "onto the path that absorbed their drainage"), vpuid, n_moved))
      }

      rec_df$toID <- to_chr
      keep <- setdiff(seq_len(nrow(rec_df)), drop_rows)
      out$reconciled <- sf::st_sf(
        rec_df[keep, , drop = FALSE],
        geometry = sf::st_sfc(geom_list[keep], crs = sf::st_crs(out$reconciled)))

      # Re-derive rec_fp_agg from the repaired reconciled and overwrite the
      # checkpoint so stage2 reads the corrected network.
      rec_fp_agg <- out$reconciled
      names(rec_fp_agg)[names(rec_fp_agg) == "ID"]   <- "flowpath_id"
      names(rec_fp_agg)[names(rec_fp_agg) == "toID"] <- "flowpath_toid"
      rec_fp_agg$flowpath_id   <- as.character(rec_fp_agg$flowpath_id)
      rec_fp_agg$flowpath_toid <- as.character(rec_fp_agg$flowpath_toid)
      rec_fp_agg$member_comid  <-
        if (is.list(rec_fp_agg$member_flowpath_id))
          vapply(rec_fp_agg$member_flowpath_id, paste, collapse = ",", character(1))
        else
          as.character(rec_fp_agg$member_flowpath_id)
      rec_fp_agg$member_flowpath_id <- NULL
      rec_fp_agg$flowline_id <- rec_fp_agg$flowpath_id
      sf::st_write(.strip_list_cols(rec_fp_agg), ckpt_rec_fp_agg, "flowpaths",
        delete_layer = TRUE, quiet = TRUE)
      message(sprintf("[stage1][%s] repaired %d divideless path(s); rec_fp_agg updated",
        vpuid, length(drop_rows)))
    }

    # Update has_divide to reflect what Step 3 actually produced (factual, not predictive).
    rec_ids <- as.character(out$reconciled$ID)
    out$reconciled$has_divide <- rec_ids %in% actual_ids
  }

  # Remove degenerate (<1 m2) divide polygons that would fail network_divides_have_area.
  # These are artifacts from partition_reconcile_divides at complex boundaries.
  if (!is.null(rec_div) && nrow(rec_div) > 0L) {
    area_m2  <- as.numeric(sf::st_area(sf::st_transform(rec_div, 5070)))
    tiny     <- !is.na(area_m2) & area_m2 < 1
    if (any(tiny)) {
      n_tiny <- sum(tiny)
      message(sprintf("[stage1][%s] removing %d degenerate divide(s) (< 1 m2) from rec_div",
        vpuid, n_tiny))
      rec_div <- rec_div[!tiny, , drop = FALSE]
      sf::st_write(.strip_list_cols(rec_div), ckpt_rec_div, "divides",
        delete_layer = TRUE, quiet = TRUE)
      # Re-apply factual has_divide after removing tiny divides
      if (!is.null(out$reconciled) && "has_divide" %in% names(out$reconciled)) {
        actual_ids_clean <- as.character(rec_div$flowpath_id)
        rec_ids_clean    <- as.character(out$reconciled$ID)
        out$reconciled$has_divide <- rec_ids_clean %in% actual_ids_clean
      }
    }
  }

  # Fill sub-pixel interior rings (holes) in reconciled divides. Raster->vector
  # reconciliation on fine DEMs (e.g. 30m GEDTM30) leaves degenerate ~single-cell
  # pinhole rings; a single divide can accumulate dozens, enough to trip the
  # no_excessive_rings invariant.
  # Stage 2/3 already run fill_small_holes, but the Stage-1 reconciled invariant
  # gate runs first: clean here so artifacts don't halt the build while genuine
  # ring problems (holes >= threshold) still fail the gate. Divides are projected
  # (m), so the m2 threshold applies directly; 1000 m2 > one 30m cell (900 m2).
  if (!is.null(rec_div) && nrow(rec_div) > 0L) {
    rec_div <- fill_small_holes(rec_div, threshold_m2 = 1000)
    sf::st_write(.strip_list_cols(rec_div), ckpt_rec_div, "divides",
      delete_layer = TRUE, quiet = TRUE)
  }

  invariants$reconciled <- tryCatch(
    check_invariants("reconciled",
      reconciled = out$reconciled, divides = rec_div,
      pois = if (exists("poi_ev2", inherits = FALSE)) poi_ev2 else NULL,
      strict = FALSE),
    error = function(e) list(error = conditionMessage(e)))

  # -- Stage 1 summary invariant --------------------------------------
  n_split <- if (!is.null(out$refactored)) nrow(out$refactored) else NA_integer_
  n_rec   <- if (!is.null(out$reconciled)) nrow(out$reconciled) else NA_integer_
  n_div   <- nrow(rec_div)
  message(sprintf("[stage1][%s] summary: %d ref -> %d split -> %d reconciled -> %d divides",
    vpuid, n_ref_fps, n_split, n_rec, n_div))

  # -- Persist artifacts ----------------------------------------------
  # refactored.gpkg / reconciled.gpkg / rec_div.gpkg were already written
  # as step checkpoints above. Only derive + write the remaining artifacts.
  files <- list(
    refactored       = ckpt_refactored,
    reconciled       = ckpt_reconciled,
    pois             = file.path(work_dir, "pois.gpkg"),
    hydrolocations   = file.path(work_dir, "hydrolocations.gpkg"),
    lineage          = file.path(work_dir, "lineage.parquet"),
    rec_fp_agg       = file.path(work_dir, "rec_fp_agg.gpkg"),
    rec_div          = ckpt_rec_div,
    timings          = file.path(work_dir, "timings.json"),
    invariants       = file.path(work_dir, "invariants.json"))

  if (!file.exists(files$pois))
    .t("persist: pois.gpkg",
      sf::st_write(.strip_list_cols(poi_ev2_dedup), files$pois,
        "pois", quiet = TRUE))

  if (!file.exists(files$hydrolocations))
    .t("persist: hydrolocations.gpkg", {
      hl <- sf::read_sf(gpkg_in, "hydrolocations", quiet = TRUE)
      sf::st_write(.strip_list_cols(hl), files$hydrolocations,
        "hydrolocations", quiet = TRUE)
    })

  if (!file.exists(files$lineage))
    .t("persist: lineage.parquet", {
      lin <- as.data.frame(lineage)
      for (nm in names(lin))
        if (is.list(lin[[nm]]))
          lin[[nm]] <- vapply(lin[[nm]], paste, collapse = ",", character(1))
      if (requireNamespace("arrow", quietly = TRUE))
        arrow::write_parquet(lin, files$lineage)
      else
        utils::write.csv(lin, sub("\\.parquet$", ".csv", files$lineage),
          row.names = FALSE)
    })

  if (!file.exists(files$rec_fp_agg))
    .t("persist: rec_fp_agg.gpkg",
      sf::st_write(.strip_list_cols(rec_fp_agg), files$rec_fp_agg,
        "flowpaths", quiet = TRUE))

  jsonlite::write_json(timings,    files$timings,    auto_unbox = TRUE,
    pretty = TRUE)
  jsonlite::write_json(invariants, files$invariants, auto_unbox = TRUE,
    pretty = TRUE, force = TRUE)

  invisible(list(files = unlist(files), timings = timings,
    invariants = invariants))
}

# Accept plain filesystem paths OR GDAL virtual URIs (/vsis3/, /vsicurl/, etc).
# file.exists() returns FALSE for virtual URIs even when GDAL can read them
# fine: we only reject plain paths that genuinely don't resolve.
.raster_accessible <- function(path) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path))
    return(FALSE)
  if (grepl("^/vsi[a-z0-9_]+/|^s3://|^https?://", path))
    return(TRUE)
  file.exists(path)
}

# Coerce list columns to comma-separated strings so sf::st_write doesn't choke.
.strip_list_cols <- function(x) {
  if (!inherits(x, "sf")) return(x)
  gcol <- attr(x, "sf_column")
  for (nm in setdiff(names(x), gcol)) {
    if (is.list(x[[nm]]))
      x[[nm]] <- vapply(x[[nm]], paste, collapse = ",", character(1))
  }
  x
}
