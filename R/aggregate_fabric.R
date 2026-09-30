#' Equal-area CRS for divide-geometry validity checks / heals
#'
#' Pick the CRS in which to test `st_is_valid()` and run geometry heals so the
#' check is meaningful and does not manufacture invalidity through projection
#' distortion. CONUS data is stored in EPSG:5070 (Conus Albers) and is valid
#' there by construction, so we return 5070 unchanged (no reprojection drift or
#' false self-intersections). oCONUS domains (AK/HI/PRVI/PACIS) fall far outside
#' the Conus Albers envelope, so we build a local Albers Equal Area conic centred
#' on the data with standard parallels at 1/6 and 5/6 of the latitude span:
#' metrically faithful to that domain and area-preserving.
#' @param x An `sf`/`sfc` object (any CRS).
#' @return An `sf::st_crs` object.
#' @keywords internal
#' @noRd
.heal_crs <- function(x) {
  bb   <- sf::st_bbox(sf::st_transform(x, 4326))
  lon0 <- mean(c(bb[["xmin"]], bb[["xmax"]]))
  lat0 <- mean(c(bb[["ymin"]], bb[["ymax"]]))
  # CONUS Albers (5070) envelope: lon -125..-66, lat 24..50
  if (lon0 > -125 && lon0 < -66 && lat0 > 24 && lat0 < 50)
    return(sf::st_crs(5070))
  lat_span <- bb[["ymax"]] - bb[["ymin"]]
  sf::st_crs(sprintf(
    paste0("+proj=aea +lat_1=%f +lat_2=%f +lat_0=%f +lon_0=%f ",
      "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"),
    bb[["ymin"]] + lat_span / 6, bb[["ymax"]] - lat_span / 6, lat0, lon0))
}

#' Aggregate a refactored fabric to the distribution network (aggregate phase)
#'
#' Consumes stage-1 artifacts; produces the aggregated network + divides
#' + remapped POIs in a single `out_vpu{id}.gpkg`. No raster access.
#'
#' Reads stage-1 outputs from disk, splits the `rec_div` layer into landscape
#' vs network on read, and exposes the stage's tunable thresholds as documented
#' arguments.
#'
#' Writes to `{out_dir}/vpu_{vpuid}/work/`:
#' - `out_vpu{id}.gpkg`: layers: `flowpaths` (aggregated), `divides`,
#'   `network` (the aspatial topology table), `pois`, `small_segments`, and
#'   optional `landscape_divides`
#' - `timings.json`
#'
#' @param vpuid VPU identifier.
#' @param fabric_dir Path to the vpu working directory holding the refactored fabric.
#' @param out_dir Parent output directory; `vpu_{vpuid}/` is created.
#' @param min_fragment_m Sub-threshold flowpath collapse length (m), Step 3c.
#' @param threshold_m2 `fill_small_holes` area threshold on divides (m^2).
#' @param gap_fill_area Gap-fill area for `hfutils::clean_gpkg_layer`.
#' @param overlap_rule Divide overlap resolution rule.
#' @param snap_interval mapshaper snap interval (degrees); collapses degenerate edges.
#' @param min_length_km Aggregate-to-distribution flowpath length threshold (km).
#' @param min_area_sqkm Aggregate-to-distribution divide area threshold (km^2).
#'   The minimum a unit must reach, applied ahead of `max_area_sqkm`.
#' @param ideal_size_sqkm Divide area the collapse aims a merged unit at (km^2).
#' @param max_area_sqkm Divide area a merged unit is not carried past (km^2),
#'   except where no receiver has room for a unit below `min_area_sqkm`.
#' @param verbose Propagated to internal helpers.
#' @return Invisibly, a list with `files`, `timings`.
#' @examples
#' \dontrun{
#' aggregate_fabric(vpuid      = "01",
#'                  fabric_dir = "out/vpu_01",
#'                  out_dir    = "out",
#'                  min_area_sqkm = 5)
#' }
#' @export
aggregate_fabric <- function(vpuid, fabric_dir, out_dir,
                             min_fragment_m = 5,
                             threshold_m2   = 1000,
                             gap_fill_area  = 0,
                             overlap_rule   = "max-area",
                             snap_interval  = 1e-8,
                             min_length_km  = 0.25,
                             min_area_sqkm  = 3,
                             ideal_size_sqkm = 10,
                             max_area_sqkm  = 15,
                             verbose = FALSE) {
  params <- list(min_fragment_m = min_fragment_m,
                 threshold_m2   = threshold_m2,
                 gap_fill_area  = gap_fill_area,
                 overlap_rule   = overlap_rule,
                 snap_interval  = snap_interval,
                 min_length_km  = min_length_km,
                 min_area_sqkm  = min_area_sqkm,
                 ideal_size_sqkm = ideal_size_sqkm,
                 max_area_sqkm  = max_area_sqkm)
  stopifnot(
    is.character(vpuid), length(vpuid) == 1L, nzchar(vpuid),
    dir.exists(fabric_dir))

  vpu_dir <- file.path(out_dir, paste0("vpu_", vpuid))
  dir.create(vpu_dir, recursive = TRUE, showWarnings = FALSE)
  work_dir <- file.path(vpu_dir, "work")
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  timings <- list()
  .t <- function(label, expr) {
    t0 <- proc.time()[["elapsed"]]
    result <- force(expr)
    timings[[label]] <<- proc.time()[["elapsed"]] - t0
    message(sprintf("[stage2][%s][%5.1fs] %s",
      vpuid, timings[[label]], label))
    result
  }

  # Stage-2 divide coverage probe (HF_PROBE_AREA), mirroring the stage-3 one.
  # Overlap between aggregate divides is the seed of the whole stage-3 cleanup
  # cascade, and it is introduced somewhere between construct_agg_divides and the
  # closing gate. Sum alone cannot locate it; sum-vs-union can.
  .probe_area_s2 <- function(tag, gpkg = NULL, layer = "divides") {
    if (!nzchar(Sys.getenv("HF_PROBE_AREA"))) return(invisible())
    tryCatch({
      g <- gpkg %||% out_gpkg
      if (!hfutils::layer_exists(g, layer)) return(invisible())
      d <- sf::st_read(g, layer, quiet = TRUE)
      s <- sum(suppressWarnings(as.numeric(sf::st_area(d))), na.rm = TRUE) / 1e6
      u <- as.numeric(sf::st_area(sf::st_union(sf::st_geometry(d)))) / 1e6
      message(sprintf(
        "[AREA-S2 %-16s][%s] divs=%d | sum=%.4f | union=%.4f | overlap=%.4f",
        tag, vpuid, nrow(d), s, u, s - u))
    }, error = function(e)
      message(sprintf("[AREA-S2 %s] error: %s", tag, conditionMessage(e))))
    invisible()
  }

  out_gpkg <- file.path(work_dir, sprintf("out_vpu%s.gpkg", vpuid))
  if (file.exists(out_gpkg)) file.remove(out_gpkg)

  # -- Read stage-1 artifacts -----------------------------------------
  s1_fp_gpkg <- file.path(fabric_dir, "work", "rec_fp_agg.gpkg")
  s1_rec_div    <- file.path(fabric_dir, "work", "rec_div.gpkg")
  s1_pois       <- file.path(fabric_dir, "work", "pois.gpkg")
  stopifnot(
    file.exists(s1_fp_gpkg), file.exists(s1_rec_div),
    file.exists(s1_pois))

  # -- Load Stage-1 reconciled flowpaths ---------------------------------
  # Naming convention used through Step 4 (aggregate_to_distribution):
  #   * `s1_fp`  : Stage-1 reconciled flowpaths (one row per rec_id).
  #                 `s1_fp$flowpath_id` IS the Stage-1 rec_id (we rename
  #                 `rec_id` -> `flowpath_id` here only because the aggregation
  #                 step expects that column name as input).
  #   * `fp_agg` : Stage-2 aggregated flowpaths produced by
  #                 `aggregate_to_distribution`, read back from `out_gpkg`.
  #                 `fp_agg$flowpath_id` is a FRESH aggregate id assigned by
  #                 the aggregation step, NOT comparable to `s1_fp$flowpath_id`.
  # The bridge between the two spaces is `fp_agg$rec_id_list` (or
  # `network$rec_id_list`), a comma-separated list of the Stage-1 rec_ids
  # absorbed into each Stage-2 aggregate.
  s1_fp <- sf::read_sf(s1_fp_gpkg, "flowpaths", quiet = TRUE)
  if ("rec_id"   %in% names(s1_fp)) names(s1_fp)[names(s1_fp) == "rec_id"]   <- "flowpath_id"
  if ("rec_toid" %in% names(s1_fp)) names(s1_fp)[names(s1_fp) == "rec_toid"] <- "flowpath_toid"
  # Sci-notation-safe id normalization. Stage-1 can persist round-number toids
  # (e.g. AK's x00000 TDX region ids) as "1e+05" while the matching flowpath_id
  # is "100000", so every `toid %in% flowpath_id` test mis-fires -> spurious
  # "dangling" toids that trip the collapse redirect and check_network_validity.
  # hfutils::hf_fmt_id re-parses "1e+05" -> "100000", making ids/toids comparable.
  s1_fp$flowpath_id   <- hfutils::hf_fmt_id(s1_fp$flowpath_id)
  s1_fp$flowpath_toid <- hfutils::hf_fmt_id(s1_fp$flowpath_toid)
  rec_div    <- sf::read_sf(s1_rec_div,    "divides",   quiet = TRUE)
  # A reference with no POIs writes no `pois` layer. Downstream already treats
  # NULL as "no POIs", so read it only when it is there.
  poi_ev2_dedup <- if (hfutils::layer_exists(s1_pois, "pois"))
    sf::read_sf(s1_pois, "pois", quiet = TRUE) else NULL
  # Flowpath ids carrying a POI, resolved once at function scope. Steps below run
  # inside .t() wrappers, which force their expression in .t()'s own frame, so a
  # value looked up there with exists(inherits = FALSE) is invisible, that is
  # exactly how the first version of the 3c.5 protection silently did nothing.
  poi_fp_ids <- if (!is.null(poi_ev2_dedup) &&
    "flowpath_id" %in% names(poi_ev2_dedup)) {
    v <- hfutils::hf_fmt_id(poi_ev2_dedup$flowpath_id)
    unique(v[!is.na(v) & nzchar(v)])
  } else character(0)
  message(sprintf("[stage2][%s] %d flowpath(s) carry a POI and are protected from collapse",
    vpuid, length(poi_fp_ids)))
  .ridc <- if ("flowpath_id" %in% names(rec_div)) "flowpath_id" else "divide_id"
  s1_fp$is_connector <- !(as.character(s1_fp$flowpath_id) %in% as.character(rec_div[[.ridc]]))

  # -- 2a. Heal invalid divide geometries AT LOAD TIME ------------------
  # Must happen before ANY spatial operation (st_union, clip, st_intersects).
  # Invalid inputs to st_union produce GEOMETRYCOLLECTION/degenerate slivers
  # that cascade through all subsequent steps. The validity check / heal runs
  # in an equal-area CRS matched to the domain (.heal_crs): EPSG:5070 for CONUS
  # (where the data already lives, so no reprojection drift), a local Albers
  # Equal Area for oCONUS.
  {
    orig_crs <- sf::st_crs(rec_div)
    heal_crs <- .heal_crs(rec_div)   # 5070 for CONUS (in-place), local AEA oCONUS
    div_proj <- sf::st_transform(rec_div, heal_crs)
    inv      <- !sf::st_is_valid(div_proj)
    if (any(inv, na.rm = TRUE)) {
      message(sprintf(
        "[stage2][%s] Step 2a: healing %d invalid divide(s) at load (equal-area check)",
        vpuid, sum(inv, na.rm = TRUE)))
      # Area-preserving heal: st_make_valid, not buffer(0). buffer(0) sheds area
      # and fails on the self-touching geometries that most need fixing. Keep all
      # polygonal parts and union into one (MULTI)POLYGON. Never drop a divide
      # here: dropping silently loses its catchment area. Any geometry still
      # flagged after re-heal is a transform-precision self-intersection on an
      # otherwise-valid polygon, resolved downstream by Step 4b clean_gpkg_layer;
      # keep it so the area is conserved.
      heal_div <- function(g) {
        mv <- tryCatch(sf::st_make_valid(g), error = function(e) NULL)
        if (is.null(mv)) return(NULL)
        p  <- tryCatch(suppressWarnings(sf::st_collection_extract(mv, "POLYGON")),
          error = function(e) mv)
        if (length(p) == 0L) return(NULL)
        u  <- tryCatch(sf::st_union(sf::st_geometry(p)), error = function(e) NULL)
        if (is.null(u) || length(u) == 0L || sf::st_is_empty(u) ||
          !isTRUE(sf::st_is_valid(u))) return(NULL)
        u
      }
      gp <- sf::st_geometry(div_proj)
      for (i in which(inv)) {
        h <- heal_div(sf::st_sfc(gp[[i]], crs = heal_crs))
        if (!is.null(h)) gp[[i]] <- h[[1L]]
      }
      div_proj <- sf::st_set_geometry(div_proj, sf::st_sfc(as.list(gp), crs = heal_crs))
      rec_div  <- sf::st_transform(div_proj, orig_crs)

      still_inv <- !sf::st_is_valid(sf::st_transform(rec_div, heal_crs))
      if (any(still_inv, na.rm = TRUE)) {
        # Round-trip re-heal (transform precision re-invalidated a few). Rebuild
        # the sfc via st_sfc so in-place element edits don't leave the class stale.
        ri        <- which(still_inv)
        rec_proj2 <- sf::st_transform(rec_div[ri, ], heal_crs)
        gvec      <- sf::st_geometry(rec_div)
        for (j in seq_along(ri)) {
          h <- heal_div(sf::st_sfc(sf::st_geometry(rec_proj2)[[j]], crs = heal_crs))
          if (!is.null(h))
            gvec[[ri[j]]] <- sf::st_geometry(sf::st_transform(h, orig_crs))[[1L]]
        }
        rec_div   <- sf::st_set_geometry(rec_div,
          sf::st_sfc(as.list(gvec), crs = orig_crs))
        still_inv <- !sf::st_is_valid(sf::st_transform(rec_div, heal_crs))
        if (any(still_inv, na.rm = TRUE))
          message(sprintf(
            "[stage2][%s] Step 2a: kept %d residual-invalid divide(s) for downstream cleanup (Step 4b) -- area conserved, not dropped",
            vpuid, sum(still_inv, na.rm = TRUE)))
      }
    }
    rm(div_proj, inv, orig_crs)
  }

  # -- Landscape vs network split (was run_pipeline.R:144-147) ---------
  lnd_mask    <- is.na(rec_div$flowpath_id) |
    grepl("^lnd-", as.character(rec_div$divide_id))
  rec_div_lnd <- rec_div[lnd_mask, ]
  rec_div_net <- rec_div[!lnd_mask, ]
  collapsed_div_ids <- character(0)   # set by Step 3c.5 below
  dedup_dropped_ids <- character(0)  # set by Step 3b.5 below

  # Baseline drainage area to conserve through aggregation: total area of the
  # full Stage-1 network divide set as written on disk, NOT the in-memory
  # rec_div_net (which the Step 2a invalid-geometry purge has pruned and the
  # Step 4_div supplement later restores from this same on-disk set). Measuring
  # the pruned set would count restored-and-aggregated divides as spurious area
  # drift; reading the full Stage-1 set makes divides_area_conserved compare
  # like-for-like against the aggregation input. Falls back to the in-memory set
  # if the on-disk read fails.
  # Summed AND union baselines. The sum guards area appearing or vanishing; the
  # union guards ground being lost, which a sum cannot see because removing an
  # overlap reduces it just as a hole does.
  .baseline_div0 <- local({
    rd0 <- tryCatch(sf::read_sf(s1_rec_div, "divides", quiet = TRUE),
      error = function(e) NULL)
    if (is.null(rd0)) rec_div_net else {
      lnd0 <- is.na(rd0$flowpath_id) | grepl("^lnd-", as.character(rd0$divide_id))
      rd0[!lnd0, ]
    }
  })
  baseline_div_area_m2 <- sum(
    suppressWarnings(as.numeric(sf::st_area(.baseline_div0))), na.rm = TRUE)
  baseline_div_cov_m2 <- tryCatch(
    as.numeric(sf::st_area(sf::st_union(sf::st_geometry(.baseline_div0)))),
    error = function(e) NA_real_)

  # -- 3a. DAG assertion + hydroseq recompute before aggregation --------
  # If the network is correct, this should always be a DAG. A cycle here
  # indicates a bug in stage-1 topology propagation; the break is an
  # emergency fallback so the batch does not silently corrupt output.
  .t("Step 3a: dag_assert_s1_fp", {
    dag_v <- hfutils::hf_assert_network_dag(s1_fp,
      id_col   = "flowpath_id",
      toid_col = "flowpath_toid")
    if (!dag_v$is_dag) {
      .hf_log_warn(sprintf(
        "[stage2][%s] BUG: cycles in s1_fp (%s) -- stage-1 toid propagation is incorrect; breaking as fallback",
        vpuid, dag_v$message))
      s1_fp <- hfutils::hf_break_cycles(s1_fp,
        id_col   = "flowpath_id",
        toid_col = "flowpath_toid")
      dag_after <- hfutils::hf_assert_network_dag(s1_fp,
        id_col   = "flowpath_id",
        toid_col = "flowpath_toid")
      if (!dag_after$is_dag)
        .hf_log_warn(sprintf("[stage2][%s] s1_fp still not a DAG after break: %s",
          vpuid, dag_after$message))
    } else {
      message(sprintf("[stage2][%s] DAG OK", vpuid))
    }
    # Always recompute hydroseq so aggregation sees topo-consistent order
    s1_fp <- hfutils::hf_recompute_hydroseq(s1_fp,
      id_col   = "flowpath_id",
      toid_col = "flowpath_toid")
    message(sprintf("[stage2][%s] hydroseq recomputed (%d fps)", vpuid, nrow(s1_fp)))
  })

  # Steps 3a.5 / 3b / 3d / 4d3 removed: the disjoint-repair cascade. It sat behind
  # HF_STAGE2_DISJOINT_REPAIR, which defaulted to "off", and all 24 VPUs of the
  # current build logged it disabled, it has not executed in any delivered run.
  # Its own note recorded why it was switched off: against a clean stage-1 partition
  # the disjoint detection fires on false positives and corrupts correct assignments
  # into shared or overlapping divides. The pure-partition work made a clean partition
  # the guaranteed input, and divides_non_overlapping plus divide_flowpath_consistent
  # now assert what the cascade used to chase.

  # -- 3b.5. Enforce unique-geometry invariant on divides -------------
  # Invariant: no two divide rows may share identical geometry.
  # When Stage-1 lineage allows a single reference polygon to be carried
  # by two flowpaths (e.g. an upstream fp whose member_comid is also
  # contained in its downstream aggregate fp's member_comid), the polygon
  # cannot be partitioned between them because only one fp's line actually
  # lies inside it. Left alone the duplicate-geometry pair produces
  # sum > union for the s2 divides layer.
  # Resolution: for each duplicate-geometry group, keep the divide
  # assigned to the fp whose line maximally intersects the polygon;
  # drop the other divide row(s). Step 3c's no_divide_idx collapse
  # then folds the now-divide-less fp into its downstream toid.
  .t("Step 3b.5: resolve duplicate-geometry divides", {
    if (nrow(rec_div_net) > 1L) {
      div_crs_b5  <- sf::st_crs(rec_div_net)
      div_geoms_b5 <- sf::st_geometry(rec_div_net)
      # Cheap key: (rounded area, rounded centroid x/y) to bin candidates;
      # confirm true equality only within each bin to avoid O(n2) cost.
      areas_km2 <- as.numeric(sf::st_area(
        sf::st_transform(rec_div_net, 5070))) / 1e6
      ctr_xy <- suppressWarnings(sf::st_coordinates(
        sf::st_centroid(rec_div_net)))
      bin_key <- paste(round(areas_km2, 6),
        round(ctr_xy[, 1], 0),
        round(ctr_xy[, 2], 0), sep = "|")
      bin_tab <- table(bin_key)
      cand_keys <- names(bin_tab)[bin_tab > 1L]
      loser_rows <- integer(0L)
      n_groups   <- 0L
      fp_ids_b5  <- as.character(s1_fp$flowpath_id)
      fp_crs_b5  <- sf::st_crs(s1_fp)
      for (k in cand_keys) {
        rows <- which(bin_key == k)
        if (length(rows) < 2L) next
        # Confirm true geometry equality (numeric tolerance via st_equals).
        ref_g <- sf::st_sfc(div_geoms_b5[[rows[1L]]], crs = div_crs_b5)
        eq_mat <- suppressWarnings(sf::st_equals(
          ref_g,
          sf::st_sfc(as.list(div_geoms_b5[rows]), crs = div_crs_b5),
          sparse = FALSE))
        eq_vec <- as.logical(eq_mat[1L, ])
        rows   <- rows[eq_vec]
        if (length(rows) < 2L) next
        n_groups <- n_groups + 1L
        # Pick the fp whose line maximally intersects the shared polygon
        # as the winner. Losers' divides are dropped.
        div_fps  <- as.character(rec_div_net$flowpath_id[rows])
        fp_ix    <- match(div_fps, fp_ids_b5)
        poly_sfc <- sf::st_sfc(div_geoms_b5[[rows[1L]]], crs = div_crs_b5)
        overlap_m <- vapply(seq_along(rows), function(i) {
          fi <- fp_ix[i]
          if (is.na(fi)) return(0)
          fp_sfc <- sf::st_sfc(sf::st_geometry(s1_fp)[[fi]],
            crs = fp_crs_b5)
          if (!identical(fp_crs_b5, div_crs_b5))
            fp_sfc <- sf::st_transform(fp_sfc, div_crs_b5)
          inter <- tryCatch(
            suppressWarnings(sf::st_intersection(fp_sfc, poly_sfc)),
            error = function(e) NULL)
          if (is.null(inter) || length(inter) == 0L) return(0)
          sum(as.numeric(sf::st_length(inter)), na.rm = TRUE)
        }, numeric(1L))
        winner <- which.max(overlap_m)
        loser_rows <- c(loser_rows, rows[-winner])
      }
      if (length(loser_rows) > 0L) {
        # Record what was dropped on purpose. These reaches carry geometry
        # IDENTICAL to the winner's, so their divide must not come back: the
        # Step 4_div supplementation reads divides off disk and would re-introduce
        # the duplicate this step just removed.
        dedup_dropped_ids <- unique(c(dedup_dropped_ids,
          as.character(rec_div_net$flowpath_id)[loser_rows]))
        rec_div_net <- rec_div_net[-loser_rows, ]
        message(sprintf(
          "[stage2][%s] Step 3b.5: dropped %d duplicate-geometry divide row(s) across %d group(s)",
          vpuid, length(loser_rows), n_groups))
      }
    }
  })

  # -- 3c. lengthkm refresh + sub-threshold fp collapse -------------------
  # Flowpath geometry comes straight from the reconciled member lines (the
  # line is intentionally NOT clipped to its divide: divides are built
  # independently by construct_agg_divides, and a stream centerline is not
  # contained by its catchment by definition: clipping only fragments
  # contiguous reaches, costing real stream length for no divide-coverage
  # gain). Lengths must therefore be (re)computed from the actual geometry,
  # since aggregate_to_distribution uses lengthkm to decide what to collapse
  # and stale values cause under-collapsing. Additionally, any fp below
  # min_fragment_m (noise-level, < 5 m by default) with a valid downstream is
  # merged into that downstream now: before the fp-divide disjoint check and
  # aggregation ever see it, to prevent near-zero divide slivers from
  # reaching mapshaper.
  .t("Step 3c: lengthkm refresh", {
    # Recompute lengths in native CRS (projected metres).
    new_len_m  <- as.numeric(sf::st_length(s1_fp))
    new_len_km <- new_len_m / 1000
    n_changed  <- sum(abs(new_len_km - s1_fp$lengthkm) > 1e-6, na.rm = TRUE)
    s1_fp$lengthkm <- new_len_km

    # Collapse fps shorter than min_fragment_m OR fps that are already
    # sub-threshold (< min_length_km) from stage1 with no upstream: isolated
    # single-levelpath heads (levelpathid == flowpath_id) that fast_collapse
    # protects from merging across levelpath boundaries but are too short to
    # stand alone. Both groups are merged into their downstream.
    noise_thresh_km <- params$min_fragment_m / 1000
    toid_col  <- "flowpath_toid"
    fp_ids    <- as.character(s1_fp$flowpath_id)
    fp_toids  <- as.character(s1_fp[[toid_col]])
    idx_by_id <- stats::setNames(seq_len(nrow(s1_fp)), fp_ids)
    # Stage1-isolated heads: levelpathid == flowpath_id AND below min_length_km
    # AND true headwaters (no fp points to them) AND have a valid downstream.
    lp_ids      <- as.character(s1_fp$levelpathid)
    up_count    <- tabulate(match(fp_toids, fp_ids), nbins = length(fp_ids))
    isolated_hw <- which(
      new_len_km < params$min_length_km &
        lp_ids == fp_ids &
        up_count == 0L
    )

    noise_idx <- union(
      which(new_len_km < noise_thresh_km),
      isolated_hw
    )

    # Also collapse fps that lost their divide in step 3b (disjoint clear).
    # These fps have no watershed territory of their own; allowing them to
    # survive into aggregate_to_distribution can create duplicate stage-2
    # divide records that represent the same parent watershed polygon.
    # Collapsing them into downstream keeps the fp/divide mapping one-to-one
    # for active watershed-bearing segments.
    div_fp_ids_after3b <- as.character(rec_div_net$flowpath_id)
    no_divide_idx <- which(!fp_ids %in% div_fp_ids_after3b)
    # EXCEPT reference connectors (reaches with areasqkm=0 in the source fabric,
    # waterbody / artificial paths that NEVER had a divide). Collapsing those into
    # a downstream catchment makes a flowpath span the waterbody (flowline far
    # from its divide). Keep them as standalone connector flowpaths (has_divide=0
    # in the final model); only collapse reaches that LOST their divide in 3b.
    if ("is_connector" %in% names(s1_fp))
      no_divide_idx <- setdiff(no_divide_idx, which(s1_fp$is_connector))
    noise_idx <- union(noise_idx, no_divide_idx)

    # A POI-bearing flowpath is never collapsed. Collapsing it removes its row
    # from the aggregation input, leaving the POI nothing to join to; and even
    # where the join survived, the merged reach's outlet is the recipient's, which
    # would leave the POI interior.
    #
    # This is the same rule the stage-1 collapse enforces by cutting its graph at
    # event nodes: a merge may not destroy a node a POI is indexed to. Keeping a
    # short reach is the accepted cost of that, and it is bounded, the reaches in
    # question are already below the size the collapse was trying to remove.
    n_before_poi <- length(noise_idx)
    noise_idx <- noise_idx[!fp_ids[noise_idx] %in% poi_fp_ids]
    # Reported unconditionally: a protection step that matches nothing looks
    # identical to a working one in the log.
    message(sprintf(
      "[stage2][%s] Step 3c.5: kept %d of %d collapse candidate(s) for carrying a POI",
      vpuid, n_before_poi - length(noise_idx), n_before_poi))

    # Keep only those with a resolvable downstream recipient.
    noise_idx <- noise_idx[
      !is.na(fp_toids[noise_idx]) &
        nzchar(fp_toids[noise_idx]) &
        fp_toids[noise_idx] != "0" &
        fp_toids[noise_idx] %in% fp_ids
    ]
    if (length(noise_idx) > 0L) {
      drop_ids <- fp_ids[noise_idx]
      # -- Vectorized merge: resolve recipient rows once, then group-union
      # A per-row loop would make one sf::st_union + one geometry assignment per
      # noise fp, each triggering CRS validation and GEOS overhead. Instead:
      # 1. Look up recipient row index for every noise fp in one shot.
      # 2. Group source indices by recipient.
      # 3. Union all sources for a given recipient in a single sf::st_union
      #    call, writing the result back once per recipient group.
      recv_rows <- idx_by_id[fp_toids[noise_idx]]
      valid_p   <- !is.na(recv_rows)
      src_rows  <- noise_idx[valid_p]
      rec_rows  <- as.integer(recv_rows[valid_p])

      if (length(src_rows) > 0L) {
        recv_uniq_r <- unique(rec_rows)
        # -- Flowpath geometry: deliberately NOT unioned ----------------------
        # The collapsed sources are sub-threshold tributary heads / no-divide
        # stubs (their own levelpath) that JOIN the recipient, often mid-reach or
        # geometrically disconnected, rather than continuing its mainstem
        # colinearly. Unioning + st_line_merge would either split the recipient
        # at the tributary junction or graft on a disconnected stub, producing a
        # multi-part geometry that Keep the
        # recipient's own single, contiguous line; the sources' area and lineage
        # are folded below (divide union + member_comid) so drainage area is
        # conserved. lengthkm already equals the recipient's geometry length, so
        # no length adjustment is needed.

        # -- Divide geometry
        div_fp_ids_v <- as.character(rec_div_net$flowpath_id)
        div_geoms    <- sf::st_geometry(rec_div_net)
        div_src_rows <- match(fp_ids[src_rows],  div_fp_ids_v)
        div_rec_rows <- match(fp_toids[src_rows], div_fp_ids_v)
        valid_d      <- !is.na(div_src_rows) & !is.na(div_rec_rows)
        if (any(valid_d)) {
          dsr <- div_src_rows[valid_d]
          drr <- div_rec_rows[valid_d]
          recv_uniq_d <- unique(drr)
          for (ri_d in recv_uniq_d) {
            sdrs <- dsr[drr == ri_d]
            unioned <- tryCatch(sf::st_union(div_geoms[c(ri_d, sdrs)]),
              error = function(e) NULL)
            if (!is.null(unioned) && !sf::st_is_empty(unioned)) {
              u1 <- unioned[[1L]]
              # If union produced a GEOMETRYCOLLECTION, extract only polygons
              u_type <- as.character(sf::st_geometry_type(sf::st_sfc(u1)))
              if (u_type == "GEOMETRYCOLLECTION") {
                u_parts <- tryCatch(
                  sf::st_collection_extract(sf::st_sfc(u1, crs = sf::st_crs(div_geoms)), "POLYGON"),
                  error = function(e) NULL)
                if (!is.null(u_parts) && length(u_parts) > 0L)
                  u1 <- u_parts[[1L]]
                else
                  next  # keep original div_geom[ri_d]
              } else if (!grepl("POLYGON", u_type)) {
                next  # non-polygon result -- keep original
              }
              div_geoms[[ri_d]] <- u1
            }
          }
          # Rebuild the sfc rather than assigning the element-mutated one. In-place
          # `[[<-` on an sfc leaves its cached per-element class attribute stale, and
          # a later st_transform then mis-dispatches with "Not compatible with
          # requested type: [type=list; target=double]". This was latent: the Step
          # 3c.5 subset immediately below happened to rebuild the column as a side
          # effect, so the staleness never surfaced until that subset was removed.
          sf::st_geometry(rec_div_net) <- sf::st_sfc(as.list(div_geoms),
            crs = sf::st_crs(div_geoms))
        }

        # -- absorbed_recs: carry the collapsed reach's own id downstream ------
        # member_comid below carries the reference COMIDs a collapsed fp represented,
        # but not the fp's own reconciled id. That id is what aggregation pastes into
        # rec_id_list, and construct_agg_divides looks up member divides by it, so
        # dropping the row dropped the id from the survivor's membership entirely.
        #
        # Step 4_lineage_repair existed to put it back afterwards, reading the
        # original stage-1 file to find rec_ids no aggregate covered and walking
        # rec_toid to a covered ancestor. Measured across 24 VPUs it reattached
        # exactly the 38741 ids dropped here, on every VPU. Recording them at the
        # point of absorption removes the need to reconstruct them later.
        #
        # Transitive by construction: a source's own absorbed list is folded in too,
        # so a chain X -> W -> Y lands X on Y.
        {
          ab_col <- if ("absorbed_recs" %in% names(s1_fp))
            as.character(s1_fp$absorbed_recs) else rep(NA_character_, nrow(s1_fp))
          fp_id_col <- as.character(s1_fp$flowpath_id)
          for (ri in recv_uniq_r) {
            srcs  <- src_rows[rec_rows == ri]
            parts <- expand_member_csv(
              c(ab_col[ri], ab_col[srcs], fp_id_col[srcs]), unique = TRUE)
            if (length(parts) > 0L) ab_col[ri] <- compose_member_csv(parts)
          }
          s1_fp$absorbed_recs <- ab_col
        }

        # -- member_comid: carry collapsed COMIDs into the downstream fp ------
        # When a sub-threshold fp is removed, its member_comid (the reference COMIDs
        # it represents) must be appended to the receiving downstream fp so that
        # aggregate_to_distribution can include them in the stage-2 lineage.
        if ("member_comid" %in% names(s1_fp)) {
          mc_col <- as.character(s1_fp$member_comid)
          for (ri in recv_uniq_r) {
            srcs   <- src_rows[rec_rows == ri]
            all_mc <- c(mc_col[ri], mc_col[srcs])
            parts  <- expand_member_csv(all_mc, unique = TRUE)
            if (length(parts) > 0L) mc_col[ri] <- compose_member_csv(parts)
          }
          s1_fp$member_comid <- mc_col
        }
      }
      # Redirect surviving fps that point AT a dropped fp to that dropped fp's
      # downstream recipient. A collapsed fp X is folded into its toid (its
      # recipient); any survivor Z with toid==X must follow to X's recipient,
      # else Z is left with a dangling toid that fails check_network_validity
      # ("All toIDs are not present in network") during aggregation. Collapses
      # chain (X->W->Y with X,W both dropped), so resolve transitively to the
      # first surviving id or terminal ("0"). DAG => no cycle, but guard anyway.
      drop_set <- drop_ids
      toid_of  <- stats::setNames(fp_toids, fp_ids)
      resolve_target <- function(t) {
        seen <- character(0L)
        while (!is.na(t) && nzchar(t) && t %in% drop_set) {
          if (t %in% seen) return("0")
          seen <- c(seen, t)
          t <- toid_of[[t]]
        }
        if (is.na(t) || !nzchar(t)) "0" else t
      }
      cur_toids <- as.character(s1_fp[[toid_col]])
      need_fix  <- which(cur_toids %in% drop_set)
      if (length(need_fix) > 0L) {
        remap <- vapply(unique(cur_toids[need_fix]), resolve_target, character(1L))
        s1_fp[[toid_col]][need_fix] <- unname(remap[cur_toids[need_fix]])
        message(sprintf(
          "[stage2][%s] Step 3c.5: redirected %d surviving toid(s) through collapsed fp(s) to their recipient",
          vpuid, length(need_fix)))
      }

      # Remove the collapsed flowpath rows. Their DIVIDES stay: a collapse folds a
      # reach's catchment into its downstream recipient, so that catchment is still
      # a member divide of the survivor and construct_agg_divides needs it to
      # dissolve the survivor's geometry.
      #
      # Dropping the divides here would only mean putting them back: the lineage
      # repair reattaches every orphaned rec_id, then the divide build re-reads
      # them off disk. Keeping them routes the geometry through the Step 3f
      # healing below instead of raw stage-1 geometry.
      #
      # aggregate_network still sees only surviving flowpaths' divides, it uses
      # them for area-based collapse decisions, which must not change here, so the
      # collapsed ids are recorded and excluded at that call site.
      s1_fp  <- s1_fp[!fp_ids %in% drop_ids, ]
      collapsed_div_ids <- unique(as.character(
        rec_div_net$flowpath_id)[as.character(rec_div_net$flowpath_id) %in% drop_ids])
      n_isolated <- sum(noise_idx %in% isolated_hw)
      n_noise    <- length(drop_ids) - n_isolated
      message(sprintf(
        "[stage2][%s] Step 3c.5: collapsed %d sub-threshold fp(s) into downstream (%d noise <%.0fm, %d isolated hw <%.2fkm)",
        vpuid, length(drop_ids), n_noise, params$min_fragment_m,
        n_isolated, params$min_length_km))
    }
    if (n_changed > 0L)
      message(sprintf("[stage2][%s] Step 3c.5: refreshed lengthkm for %d fp(s) after clip",
        vpuid, n_changed))
  })

  # -- 3c.6 TRIPWIRE: collapse must not leave a TRUE dangling toid -------------
  # 3c.5 removes noise/headwater/divide-less fps but does not rewrite the toid of
  # a surviving fp that pointed at a dropped fp. A subgraph of a DAG is itself a
  # DAG, so the collapse cannot create a cycle, but a survivor pointing to a
  # *removed regular fp* is a real defect. Terminal toids (>=1e9 nexus ids, or
  # "0") legitimately are not flowpath ids, so exclude them; flag only toids that
  # reference a non-terminal id absent from the surviving fp set.
  {
    .ids   <- as.character(s1_fp$flowpath_id)
    .toids <- as.character(s1_fp$flowpath_toid)
    .tnum  <- suppressWarnings(as.numeric(.toids))
    .dang  <- which(!is.na(.toids) & nzchar(.toids) & .toids != "0" &
      !(.toids %in% .ids) &
      !(is.finite(.tnum) & .tnum >= 1e9))   # exclude terminals
    if (length(.dang) > 0L)
      .hf_log_warn(sprintf(
        "[stage2][%s] Step 3c.6: %d surviving fp(s) have a dangling toid to a removed non-terminal fp (e.g. %s -> %s)",
        vpuid, length(.dang),
        paste(utils::head(.ids[.dang], 3), collapse = ","),
        paste(utils::head(.toids[.dang], 3), collapse = ",")))
    rm(.ids, .toids, .tnum, .dang)
  }

  # -- 3d (was 3c). Post-clip disjoint check, any residual disjoint pairs ----
  # A second pass after clip catches any fp-divide assignments that clip's
  # min_inside_frac guard left in place (e.g., exactly at the 10% boundary).
  # Step 3d removed with the rest of the disjoint-repair cascade; see Step 3b.

  # -- 3e. Repair invalid divide geometries before aggregation ----------
  # GEOS TopologyException ("unable to assign free hole to a shell") and other
  # crashes in st_intersects inside aggregate_to_distribution come from
  # floating-point precision artifacts that survive st_make_valid. Check / heal
  # in an equal-area CRS matched to the domain (.heal_crs): 5070 for CONUS (in
  # place, since the data is already 5070 and valid there by construction), a
  # local Albers Equal Area for oCONUS. Reprojecting CONUS data to a distant CRS
  # would manufacture false self-intersections that the heal then collapses to
  # empty geometries.
  orig_crs    <- sf::st_crs(rec_div_net)
  heal_crs    <- .heal_crs(rec_div_net)
  div_proj    <- sf::st_transform(rec_div_net, heal_crs)
  inv_mask    <- !sf::st_is_valid(div_proj)
  n_invalid   <- sum(inv_mask, na.rm = TRUE)
  if (n_invalid > 0L) {
    message(sprintf("[stage2][%s] Step 3e: healing %d invalid divide geometry(ies) (equal-area check)",
      vpuid, n_invalid))
    for (i in which(inv_mask)) {
      g <- sf::st_sfc(sf::st_geometry(div_proj)[[i]], crs = heal_crs)
      g_fix <- tryCatch(sf::st_buffer(g, 0), error = function(e) NULL)
      if (!is.null(g_fix) && sf::st_is_valid(g_fix) &&
        grepl("POLYGON", as.character(sf::st_geometry_type(g_fix)))) {
        sf::st_geometry(div_proj)[i] <- g_fix
      } else {
        g_mv <- tryCatch(
          sf::st_collection_extract(sf::st_make_valid(g), "POLYGON"),
          error = function(e) NULL)
        if (!is.null(g_mv) && nrow(g_mv) > 0 && sf::st_is_valid(g_mv[1L, ]))
          sf::st_geometry(div_proj)[i] <- sf::st_geometry(g_mv)[1L]
      }
    }
    rec_div_net <- sf::st_transform(div_proj, orig_crs)
    # second pass in projected CRS to catch any round-trip residuals
    div_proj2  <- sf::st_transform(rec_div_net, heal_crs)
    inv_mask2  <- !sf::st_is_valid(div_proj2)
    if (any(inv_mask2, na.rm = TRUE)) {
      div_proj2[inv_mask2, ] <- sf::st_make_valid(div_proj2[inv_mask2, ])
      rec_div_net <- sf::st_transform(div_proj2, orig_crs)
    }
    rec_div_net <- rec_div_net[!sf::st_is_empty(rec_div_net), ]
  }

  # -- 3f. Force all divide geometries to MULTIPOLYGON -----------------
  # Per-row handler: covers GEOMETRYCOLLECTION, degenerate rings (<4 pts),
  # and any other non-standard type produced by st_union / clip / buffer.
  # Never calls sf::st_geometry(x)[i] <- directly (that assignment triggers
  # sf geometry validation and can itself throw "polygons require at least
  # 4 points"). Instead, rebuild the entire sfc from a list of fixed sfg.
  #
  # NOTE: aggregate_to_distribution casts MULTIPOLYGON -> POLYGON
  # internally, which is where degenerate rings (< 4 pts) are detected even when
  # MULTIPOLYGON cast succeeds. After all geometries are MULTIPOLYGON, a batch
  # st_buffer(0) strips any remaining degenerate rings.
  {
    div_crs   <- sf::st_crs(rec_div_net)
    # Phase 1: per-row: unwrap GEOMETRYCOLLECTIONs and cast to MULTIPOLYGON
    geom_list <- lapply(seq_len(nrow(rec_div_net)), function(i) {
      g     <- sf::st_sfc(sf::st_geometry(rec_div_net)[[i]], crs = div_crs)
      gtype <- as.character(sf::st_geometry_type(g))
      # Unwrap GEOMETRYCOLLECTION -> take first polygon part
      if (grepl("COLLECTION", gtype)) {
        gp <- tryCatch(sf::st_collection_extract(g, "POLYGON"),
          error = function(e) NULL)
        if (is.null(gp) || length(sf::st_geometry(gp)) == 0L) return(NULL)
        g     <- gp[1L, ]
        gtype <- as.character(sf::st_geometry_type(g))
      }
      if (!grepl("POLYGON", gtype)) return(NULL)
      # Cast to MULTIPOLYGON; on error try st_buffer(0) heal first
      tryCatch(
        sf::st_cast(g, "MULTIPOLYGON", warn = FALSE)[[1L]],
        error = function(e) {
          g2 <- tryCatch(sf::st_buffer(g, 0), error = function(e2) NULL)
          if (!is.null(g2))
            tryCatch(sf::st_cast(g2, "MULTIPOLYGON", warn = FALSE)[[1L]],
              error = function(e3) NULL)
          else
            NULL
        })
    })
    castable  <- !vapply(geom_list, is.null, logical(1L))
    n_dropped <- sum(!castable)
    n_non_poly <- sum(!grepl("POLYGON",
      as.character(sf::st_geometry_type(rec_div_net))))
    if (n_non_poly > 0L || n_dropped > 0L)
      message(sprintf(
        "[stage2][%s] Step 3f: %d non-polygon type(s) normalised; %d dropped (uncastable)",
        vpuid, n_non_poly, n_dropped))
    if (n_dropped > 0L) {
      rec_div_net <- rec_div_net[castable, ]
      geom_list   <- geom_list[castable]
    }
    if (nrow(rec_div_net) > 0L)
      sf::st_geometry(rec_div_net) <- sf::st_sfc(geom_list, crs = div_crs)
    # Phase 2: batch st_buffer(0), all geometries are now MULTIPOLYGON.
    # This strips any degenerate interior rings (< 4 pts) that survived the
    # per-row cast and would crash the MULTIPOLYGON -> POLYGON enforcer.
    if (nrow(rec_div_net) > 0L) {
      healed <- tryCatch(sf::st_buffer(rec_div_net, 0), error = function(e) NULL)
      if (!is.null(healed) && inherits(healed, "sf") &&
        nrow(healed) == nrow(rec_div_net))
        rec_div_net <- healed
    }
    # st_buffer(0) (and the 3e heal round-trip) can collapse a sliver/degenerate
    # divide to an EMPTY geometry while keeping its attribute row. Step 3e drops
    # empties, but this Phase-2 buffer runs after it, so drop again here: an empty
    # divide that survives to Step 4_div feeds union_polygons -> terra::aggregate,
    # which drops the empty geometry but keeps the row, raising "nrow dataframe
    # does not match nrow geometry". Empty geoms carry no area, so this is
    # area-safe.
    if (nrow(rec_div_net) > 0L) {
      emp <- sf::st_is_empty(rec_div_net)
      if (any(emp)) {
        message(sprintf(
          "[stage2][%s] Step 3f: dropped %d divide(s) emptied by heal/buffer",
          vpuid, sum(emp)))
        rec_div_net <- rec_div_net[!emp, ]
      }
    }
  }

  # -- 4. Aggregate network (flowpaths + topology, no divide geometry)
  # Replaces aggregate_to_distribution: runs the same mainstem + headwater
  # collapse logic but skips the slow divide union / st_snap path entirely.
  # Divides are constructed from rec_div in Step 4_div below.
  agg_net <- .t("Step 4: aggregate_network", {
    aggregate_network(
      flowpath        = s1_fp,
      # Surviving flowpaths' divides only. The divides of flowpaths collapsed at
      # Step 3c.5 are retained in rec_div_net as member divides for
      # construct_agg_divides, but they must not reach set formation, whose
      # area-based collapse decisions are scoped to flowpaths still in the network.
      divide          = if (length(collapsed_div_ids))
        rec_div_net[!as.character(rec_div_net$flowpath_id) %in% collapsed_div_ids, ]
        else rec_div_net,
      pois            = poi_ev2_dedup,
      vpu             = vpuid,
      min_length_km   = params$min_length_km,
      min_area_sqkm   = params$min_area_sqkm,
      ideal_size_sqkm = params$ideal_size_sqkm,
      max_area_sqkm   = params$max_area_sqkm,
      # connector reaches survive Step 3c; their CURRENT flowpath_ids gate
      # the cs_group break + collapse/absorb guards (no catchment<->connector merge)
      connector_ids   = if ("is_connector" %in% names(s1_fp))
        as.character(s1_fp$flowpath_id[which(s1_fp$is_connector)]) else NULL)
  })

  # Write flowpaths + network + pois to out_gpkg (divides written in Step 4_div)
  .t("Step 4_write: write aggregated network to out_gpkg", {
    if (file.exists(out_gpkg)) file.remove(out_gpkg)
    sf::st_write(agg_net$flowpaths, out_gpkg, layer = "flowpaths",
      delete_layer = TRUE, quiet = TRUE)
    {
      con <- DBI::dbConnect(RSQLite::SQLite(), out_gpkg)
      DBI::dbWriteTable(con, "network", agg_net$network, overwrite = TRUE)
      DBI::dbDisconnect(con)
    }
    if (!is.null(agg_net$pois) && nrow(agg_net$pois) > 0L)
      sf::st_write(agg_net$pois, out_gpkg, layer = "pois",
        delete_layer = TRUE, quiet = TRUE)
    message(sprintf("[stage2][%s] Step 4_write: %d flowpaths, %d network rows",
      vpuid, nrow(agg_net$flowpaths), nrow(agg_net$network)))
  })

  rm(agg_net)

  # -- 4_rename. Rename `flowline_id` -> `rec_id_list` on Stage-2 layers
  # The column emitted by the aggregation step as `flowline_id` is semantically a
  # comma-separated list of Stage-1 reconciled IDs (rec_id), NOT a single
  # per-row flowline id (which is what Stage-3 emits in its own `flowline_id`
  # column on the `flowlines` layer). Rename at this boundary so downstream
  # code is no longer overloaded.
  .t("Step 4_rename: flowline_id -> rec_id_list", {
    for (lyr in c("network", "flowpaths")) {
      df <- sf::st_read(out_gpkg, lyr, quiet = TRUE)
      if ("flowline_id" %in% names(df)) {
        names(df)[names(df) == "flowline_id"] <- "rec_id_list"
        sf::st_write(df, out_gpkg, layer = lyr,
          delete_layer = TRUE, quiet = TRUE)
      }
    }
  })

  # Step 4_toid removed: reported "toids OK post-aggregation (no repair needed)" on
  # all 24 VPUs. Asserted by aggregated_toid_resolves.

  # Step 4_selfloop removed: reported "no self-loop toids (skipping)" on all 24 VPUs.
  # Asserted by aggregated_toid_resolves.

  # Step 4_dedup_recid removed: it enforced that no reconciled reach is claimed by
  # more than one aggregate. It reported "rec_id partition OK (no duplicates)" on all
  # 24 VPUs, it never had anything to enforce. The property is load-bearing, since
  # construct_agg_divides dissolves rec_div by member, so it is now asserted by the
  # rec_id_partition invariant instead of repaired.

  # -- 4_lineage_repair. Attach orphan rec_id lineage to downstream survivor
  # Invariant: every Stage-1 reconciled rec_id (s1_fp$flowpath_id) must
  # be represented somewhere in the Stage-2 network (either as a member of a
  # network row's rec_id_list, or, for sub-threshold tribs that
  # the aggregation step dropped: absorbed into its downstream survivor's lineage).
  # The aggregate-to-distribution step can drop a rec_id (typically a no-divide
  # headwater or a split residual) without propagating its member_comid into any
  # downstream aggregate. The area is folded into a downstream divide, so this is
  # a lineage loss rather than a spatial one, and the repair is: walk rec_toid
  # downstream to the first covered rec_id, then
  # attach the orphan's member_comid to that survivor fp's lineage columns
  # and append one network row per orphan COMID with hf_id = COMID and
  # flowpath_id = the survivor aggregate.
  .t("Step 4_lineage_repair: attach orphan rec_id lineage downstream", {
    tryCatch(
      {
        net <- sf::st_read(out_gpkg, "network", quiet = TRUE)
        fp_agg <- sf::st_read(out_gpkg, "flowpaths", quiet = TRUE)
        if (!all(c("flowpath_id", "rec_id_list") %in% names(net)) ||
          !"flowpath_id" %in% names(fp_agg)) {
          message(sprintf("[stage2][%s] Step 4_lineage_repair: missing required cols (skipping)",
            vpuid))
        } else {
          net$flowpath_id <- as.character(net$flowpath_id)
          net$rec_id_list <- as.character(net$rec_id_list)
          fp_agg$flowpath_id <- as.character(fp_agg$flowpath_id)

          # (fp, rec_id) edge set: expand each net row's comma-list once.
          edges <- do.call(rbind, lapply(seq_len(nrow(net)), function(i) {
            rs <- expand_member_csv(net$rec_id_list[i])
            if (length(rs) == 0L) return(NULL)
            data.frame(fp = net$flowpath_id[i], rec = rs, stringsAsFactors = FALSE)
          }))
          edges    <- unique(edges)
          covered  <- unique(edges$rec)
          # Read the ORIGINAL Stage-1 file: runtime s1_fp has been filtered
          # by intermediate steps (Step 3b.5 dedup, Step 3c clip, no-divide
          # collapse, etc.), so it already excludes the rec_ids the aggregation step
          # dropped.  We need the full Stage-1 rec set to detect attrition.
          s1_full <- sf::st_drop_geometry(sf::read_sf(s1_fp_gpkg, "flowpaths",
            quiet = TRUE))
          if ("rec_id"   %in% names(s1_full)) names(s1_full)[names(s1_full) == "rec_id"]   <- "flowpath_id"
          if ("rec_toid" %in% names(s1_full)) names(s1_full)[names(s1_full) == "rec_toid"] <- "flowpath_toid"
          s1_recs  <- as.character(s1_full$flowpath_id)
          lost     <- setdiff(s1_recs, covered)

          if (length(lost) == 0L) {
            message(sprintf("[stage2][%s] Step 4_lineage_repair: no orphan rec_ids (lineage complete)",
              vpuid))
          } else {
            toid_by_rec <- stats::setNames(as.character(s1_full$flowpath_toid),
              s1_recs)
            mc_by_rec   <- stats::setNames(as.character(s1_full$member_comid),
              s1_recs)
            lp_by_rec   <- stats::setNames(as.character(s1_full$levelpathid),
              s1_recs)

            # Walk rec_toid until we hit a covered ancestor (bounded hops).
            covered_set <- new.env(hash = TRUE, parent = emptyenv())
            for (c in covered) assign(c, TRUE, envir = covered_set)
            walk <- function(r) {
              cur <- r
              hops <- 0L
              while (hops < 200L) {
                nx <- toid_by_rec[cur]
                if (is.na(nx) || !nzchar(nx) || nx == "0" || nx == cur)
                  return(NA_character_)
                if (exists(nx, envir = covered_set, inherits = FALSE)) return(nx)
                cur <- nx
                hops <- hops + 1L
              }
              NA_character_
            }
            landing  <- vapply(lost, walk, character(1L))
            resolved <- lost[!is.na(landing)]

            if (length(resolved) == 0L) {
              message(sprintf(
                "[stage2][%s] Step 4_lineage_repair: %d orphan rec_id(s) with no covered downstream ancestor -- lineage gap unrecoverable",
                vpuid, length(lost)))
            } else {
              # rec_id -> owner fp: pick fp whose mainstem matches the landing
              # rec_id's levelpathid (its same-levelpath aggregate); fall back to
              # the first fp that lists the landing rec_id.
              ms_by_fp  <- stats::setNames(as.character(fp_agg$mainstem),
                fp_agg$flowpath_id)
              fp_by_rec <- split(edges$fp, edges$rec)
              owner_of  <- function(rid) {
                fps <- unique(fp_by_rec[[rid]])
                if (length(fps) == 0L) return(NA_character_)
                lp <- lp_by_rec[rid]
                if (!is.na(lp)) {
                  hit <- which(ms_by_fp[fps] == lp)
                  if (length(hit) >= 1L) return(fps[hit[1L]])
                }
                fps[1L]
              }

              adds <- do.call(rbind, lapply(resolved, function(r) {
                owner <- owner_of(landing[r])
                if (is.na(owner)) return(NULL)
                cmids <- expand_member_csv(mc_by_rec[r])
                if (length(cmids) == 0L) return(NULL)
                data.frame(fp = owner, rec = r, comid = cmids,
                  stringsAsFactors = FALSE)
              }))

              if (is.null(adds) || nrow(adds) == 0L) {
                message(sprintf(
                  "[stage2][%s] Step 4_lineage_repair: %d orphan resolved but no comids to attach",
                  vpuid, length(resolved)))
              } else {
                n_attached_recs   <- length(unique(adds$rec))
                n_attached_comids <- nrow(adds)

                # 1. Update flowpaths: append COMIDs to member_comid and append
                #    the orphan rec_ids to rec_id_list (comma-list).
                for (f in unique(adds$fp)) {
                  row_ix <- which(fp_agg$flowpath_id == f)
                  if (length(row_ix) != 1L) next
                  new_cmids <- unique(adds$comid[adds$fp == f])
                  new_recs  <- unique(adds$rec[adds$fp == f])
                  if ("member_comid" %in% names(fp_agg)) {
                    cur <- expand_member_csv(fp_agg$member_comid[row_ix])
                    fp_agg$member_comid[row_ix] <-
                      compose_member_csv(c(cur, new_cmids))
                  }
                  if ("rec_id_list" %in% names(fp_agg)) {
                    cur_r <- expand_member_csv(fp_agg$rec_id_list[row_ix])
                    fp_agg$rec_id_list[row_ix] <-
                      compose_member_csv(c(cur_r, new_recs))
                  }
                }

                # 2. Update network: append one row per (owner_fp, orphan COMID),
                #    using an existing row of that fp as a template, and also
                #    rewrite the rec_id_list comma-list on existing rows to
                #    include the absorbed rec_ids.  Parse .N suffix on the
                #    orphan COMID so split-children get hf_id_part > 1 instead
                #    of being silently flattened to part 1.
                new_rows <- vector("list", nrow(adds))
                parts    <- parse_member_part(adds$comid)
                for (i in seq_len(nrow(adds))) {
                  fp_i <- adds$fp[i]
                  tmpl_ix <- which(net$flowpath_id == fp_i)
                  if (length(tmpl_ix) == 0L) next
                  tmpl <- net[tmpl_ix[1L], , drop = FALSE]
                  if ("hf_id"      %in% names(tmpl)) tmpl$hf_id      <- parts$hf_id[i]
                  if ("hf_id_part" %in% names(tmpl)) tmpl$hf_id_part <- parts$hf_id_part[i]
                  new_rows[[i]] <- tmpl
                }
                new_rows <- new_rows[!vapply(new_rows, is.null, logical(1L))]
                net2 <- if (length(new_rows) > 0L) rbind(net, do.call(rbind, new_rows))
                else net
                for (f in unique(adds$fp)) {
                  new_recs <- unique(adds$rec[adds$fp == f])
                  ix <- which(net2$flowpath_id == f)
                  if (length(ix) == 0L) next
                  cur_r  <- expand_member_csv(net2$rec_id_list[ix[1L]])
                  merged <- compose_member_csv(c(cur_r, new_recs))
                  net2$rec_id_list[ix] <- merged
                }

                sf::st_write(net2,   out_gpkg, layer = "network",
                  delete_layer = TRUE, quiet = TRUE)
                sf::st_write(fp_agg, out_gpkg, layer = "flowpaths",
                  delete_layer = TRUE, quiet = TRUE)

                n_unresolved <- length(lost) - length(resolved)
                message(sprintf(
                  "[stage2][%s] Step 4_lineage_repair: attached %d orphan rec_id(s) (%d COMID(s)) to downstream survivor fps%s",
                  vpuid, n_attached_recs, n_attached_comids,
                  if (n_unresolved > 0L)
                    sprintf("; %d unresolvable (no covered ancestor)", n_unresolved)
                  else ""))
              }
            }
          }
        }
      },
      error = function(e) {
        message(sprintf("[stage2][%s] Step 4_lineage_repair ERROR: %s",
          vpuid, conditionMessage(e)))
        stop(e)
      })
  })

  # -- 4_div. Construct divides from rec_div ----------------------------
  # Run AFTER lineage_repair so rec_id_list has had orphan rec_ids appended.
  # Pass Stage-1 hydroseq so construct_agg_divides uses the downstream-most
  # rec_id per fp, the one whose divide polygon contains the fp outlet.
  # -- 4_junction. Close the routing set on coincident flowlines ---------------
  # Runs BEFORE construct_agg_divides so the divide split follows from the corrected
  # member partition and is exact by construction. Gated while it is proven out.
  if (identical(Sys.getenv("HF_JUNCTION_PROMOTE", ""), "1")) {
    .t("Step 4_junction: split absorbed chains that receive routed flow", {
      tryCatch(.split_junction_absorbed(out_gpkg, s1_fp, vpuid = vpuid),
        error = function(e) {
          .hf_log_warn(sprintf("[stage2][%s] Step 4_junction non-fatal: %s", vpuid,
            conditionMessage(e)))
          invisible(0L)
        })
    })
  }

  .t("Step 4_div: construct_agg_divides", {
    fp_for_div  <- sf::st_read(out_gpkg, "flowpaths", quiet = TRUE)
    net_for_div <- as.data.frame(sf::st_read(out_gpkg, "network", quiet = TRUE))

    # Supplement rec_div_net with any Stage-1 rec_ids that Stage-2
    # preprocessing dropped (Step 2a invalid-geometry purge, Step 3b/3c.5
    # disjoint cascade).  construct_agg_divides silently skips rec_ids it
    # can't find, leaving those fps without a divide.  Re-load the missing
    # entries from the on-disk Stage-1 rec_div so every fp in the network
    # gets a divide attempt. A genuinely misassigned pairing is reported by the
    # divide_flowpath_consistent and divide_has_member_flowline invariants.
    {
      all_recs <- unique(trimws(unlist(strsplit(
        as.character(net_for_div$rec_id_list), ","))))
      all_recs <- all_recs[nzchar(all_recs)]
      present  <- as.character(rec_div_net$flowpath_id)
      # Reaches whose divide Step 3b.5 removed as a duplicate are excluded: their
      # geometry is already present under the winner's id, so re-adding it would
      # re-create the duplicate and double-count that catchment.
      missing  <- setdiff(setdiff(all_recs, present), dedup_dropped_ids)
      if (length(missing) > 0L) {
        rd_full <- tryCatch(
          sf::read_sf(s1_rec_div, "divides", quiet = TRUE),
          error = function(e) NULL)
        if (!is.null(rd_full)) {
          lnd_m  <- is.na(rd_full$flowpath_id) |
            grepl("^lnd-", as.character(rd_full$divide_id))
          rd_full_net <- rd_full[!lnd_m &
            as.character(rd_full$flowpath_id) %in% missing, ]
          if (nrow(rd_full_net) > 0L) {
            rd_full_net <- tryCatch(
              sf::st_transform(rd_full_net, sf::st_crs(rec_div_net)),
              error = function(e) NULL)
            if (!is.null(rd_full_net)) {
              # Align to rec_div_net schema geometry-safely. A column-reorder +
              # rbind can leave the active sf_column pointing at an NA-filled data
              # column when rd_full's geometry column name differs, yielding an sf
              # whose data nrow != geometry length ("nrow dataframe does not match
              # nrow geometry"). Build the supplement explicitly: data columns
              # aligned, geometry carried via st_geometry(), reassembled with
              # st_sf, deduped on rec_id.
              gcol      <- attr(rec_div_net, "sf_column")
              data_cols <- setdiff(names(rec_div_net), gcol)
              sup_g     <- sf::st_geometry(rd_full_net)
              sup_d     <- sf::st_drop_geometry(rd_full_net)
              for (nm in setdiff(data_cols, names(sup_d))) sup_d[[nm]] <- NA
              sup_d     <- sup_d[, data_cols, drop = FALSE]
              sup_d[[gcol]] <- sup_g
              sup_sf    <- sf::st_sf(sup_d, sf_column_name = gcol,
                crs = sf::st_crs(rec_div_net))
              # never re-add a rec_id already present (avoids duplicate divides
              # that would desync the construct_agg_divides union)
              sup_sf    <- sup_sf[!(as.character(sup_sf$flowpath_id) %in%
                as.character(rec_div_net$flowpath_id)), , drop = FALSE]
              if (nrow(sup_sf) > 0L) {
                rec_div_net <<- rbind(rec_div_net, sup_sf)
                message(sprintf(
                  "[stage2][%s] Step 4_div: supplemented %d missing rec_div(s) from Stage-1",
                  vpuid, nrow(sup_sf)))
              }
            }
          }
        }
      }
    }

    agg_divides <- construct_agg_divides(
      network   = net_for_div,
      rec_div   = rec_div_net,
      flowpaths = fp_for_div)
    sf::st_write(agg_divides, out_gpkg, layer = "divides",
      delete_layer = TRUE, quiet = TRUE)
    message(sprintf("[stage2][%s] Step 4_div: %d divide(s) written",
      vpuid, nrow(agg_divides)))
    rm(agg_divides, fp_for_div, net_for_div)
  })

  .probe_area_s2("after-4_div")
  # A multi-part aggregate line is left as it is. The merge is already attempted
  # upstream in union_linestrings, so anything still multi-part has genuinely
  # non-contiguous members; picking one component would discard the rest. Stage-3
  # extend_flowpaths_to_mainstem rebuilds these lines from member reaches anyway.
  # Non-contiguous membership surfaces through the
  # aggregated_flowpath_single_part measurement.

  # -- 4a. Landscape divides layer ------------------------------------
  if (nrow(rec_div_lnd) > 0L) {
    .t("Step 4a: write landscape_divides layer", {
      existing <- sf::st_read(out_gpkg, "divides", quiet = TRUE)
      lnd_out  <- sf::st_transform(rec_div_lnd, sf::st_crs(existing))
      keep_nm  <- names(lnd_out)[!vapply(sf::st_drop_geometry(lnd_out),
        is.list, logical(1))]
      lnd_out  <- lnd_out[, keep_nm]
      lnd_out$areasqkm     <- as.numeric(sf::st_area(lnd_out)) / 1e6
      lnd_out$has_flowline <- FALSE
      lnd_out$type         <- "landscape"
      sf::st_write(lnd_out, out_gpkg, "landscape_divides",
        delete_layer = TRUE, quiet = TRUE)
    })
  }

  # -- 4a2-4d. Post-aggregation cleanup -------------------------------
  .t("Step 4a2: restore_aggregated_divides", {
    # Pass reference divides so restore_aggregated_divides can detect and fill
    # Stage-1 rec_div dropout (large reference catchments where reconcile_divides
    # only delineated a fraction of the true area, leaving spatial gaps).
    s1_ref_gpkg <- file.path(fabric_dir, "work", sprintf("ref_vpu%s.gpkg", vpuid))
    ref_div_for_restore <- if (file.exists(s1_ref_gpkg) &&
      hfutils::layer_exists(s1_ref_gpkg, "divides")) {
      tryCatch(sf::read_sf(s1_ref_gpkg, "divides", quiet = TRUE),
        error = function(e) NULL)
    } else NULL
    restore_aggregated_divides(out_gpkg, rec_div_net,
      reference_divides = ref_div_for_restore,
      verbose = verbose)
  })
  .t("Step 4a3: remap_poi_flowpath_ids", remap_poi_flowpath_ids(out_gpkg))

  # -- 4b-pre. Pre-clean geometry repair --------------------------------
  # Mapshaper's -clean snap-to-interval collapses geometries with floating-
  # point precision artifacts into zero-area slivers, which it then drops.
  # This step normalises all divide coordinates to a 1e4 grid (0.1 mm in
  # projected metres) and repairs validity before clean runs, so lineage divides
  # that would otherwise be dropped survive the clean pass intact. Using 1e4
  # (slightly coarser than the 1e5 snap used in restore_aggregated_divides) gives
  # mapshaper a cleaner input without visibly shifting divide boundaries.
  .probe_area_s2("after-4a2")
  .t("Step 4b-pre: pre-clean divide geometry repair", {
    d_pre <- sf::st_read(out_gpkg, "divides", quiet = TRUE)
    n_invalid_pre <- sum(!sf::st_is_valid(d_pre), na.rm = TRUE)

    # Precision normalisation + validity repair in one pass.
    d_pre_r <- tryCatch(
      sf::st_make_valid(sf::st_set_precision(d_pre, 1e4)),
      error = function(e) d_pre)

    # Guard: if make_valid empties any geometry, revert that row to original.
    empty_mask <- sf::st_is_empty(d_pre_r)
    if (any(empty_mask, na.rm = TRUE)) {
      sf::st_geometry(d_pre_r)[empty_mask] <- sf::st_geometry(d_pre)[empty_mask]
      message(sprintf("[stage2][%s] Step 4b-pre: reverted %d empty geometry(ies) to original",
        vpuid, sum(empty_mask, na.rm = TRUE)))
    }

    n_repaired <- sum(!sf::st_is_valid(d_pre) & sf::st_is_valid(d_pre_r), na.rm = TRUE)

    # Near-miss buffer: divides that are within 250m of their flowpath but don't
    # intersect it is reported by the divide_flowpath_consistent invariant.
    # Buffer these divides slightly to bridge the raster-snap gap before clean runs.
    tryCatch(
      {
        fp_pre <- sf::st_read(out_gpkg, "flowpaths", quiet = TRUE)
        fp_df  <- sf::st_drop_geometry(fp_pre)
        div_df <- sf::st_drop_geometry(d_pre_r)
        # match divide -> fp by divide_id == flowpath_id
        fp_idx <- match(as.character(div_df$divide_id), as.character(fp_df$flowpath_id))
        valid  <- !is.na(fp_idx)
        if (any(valid)) {
          div_t  <- sf::st_transform(d_pre_r[valid, ], 5070)
          fp_t   <- sf::st_transform(fp_pre[fp_idx[valid], ], 5070)
          inter  <- suppressWarnings(as.logical(
            sf::st_intersects(div_t, fp_t, sparse = FALSE)[
              cbind(seq_len(sum(valid)), seq_len(sum(valid)))]))
          inter[is.na(inter)] <- FALSE
          # for non-intersecting, check distance
          cands <- which(!inter)
          if (length(cands) > 0L) {
            dists <- suppressWarnings(as.numeric(
              sf::st_distance(div_t[cands, ], fp_t[cands, ], by_element = TRUE)))
            near  <- cands[is.finite(dists) & dists > 0 & dists <= 250]
            if (length(near) > 0L) {
              # Union the divide with a thin buffer around its flowpath to
              # guarantee intersection regardless of raster-snap offset.
              fp_near <- fp_t[near, ]
              fp_buf  <- suppressWarnings(sf::st_buffer(fp_near, 50))
              div_near <- div_t[near, ]
              buffered <- tryCatch(
                sf::st_sf(
                  sf::st_drop_geometry(div_near)[, intersect(c("divide_id", "flowpath_toid"), names(div_near)), drop = FALSE],
                  geometry = sf::st_sfc(mapply(function(d, f) {
                    tryCatch(sf::st_geometry(sf::st_union(d, f))[[1L]],
                      error = function(e) sf::st_geometry(d)[[1L]])
                  }, sf::st_geometry(div_near), sf::st_geometry(fp_buf), SIMPLIFY = FALSE),
                  crs = 5070)),
                error = function(e) div_near)
              buffered <- tryCatch(sf::st_make_valid(buffered), error = function(e) buffered)
              # cast back to POLYGON (buffer of POLYGON is POLYGON)
              geom_full <- sf::st_geometry(d_pre_r)
              orig_idx  <- which(valid)[near]
              for (ki in seq_along(near)) {
                bg <- sf::st_transform(buffered[ki, ], sf::st_crs(d_pre_r))
                geom_full[[orig_idx[ki]]] <- sf::st_geometry(bg)[[1L]]
              }
              sf::st_geometry(d_pre_r) <- sf::st_sfc(as.list(geom_full),
                crs = sf::st_crs(d_pre_r))
              message(sprintf("[stage2][%s] Step 4b-pre: buffered %d near-miss divide(s) (within 250m of fp)",
                vpuid, length(near)))
            }
          }
        }
      },
      error = function(e) NULL)

    # Polygonal coercion. st_make_valid (precision snap above) and the buffer
    # union can turn a self-touching divide into a GEOMETRYCOLLECTION of a
    # POLYGON plus dangling LINESTRING spikes. ogr2ogr writes those collections
    # to GeoJSON and `mapshaper -clean` (Step 4b) then aborts with "Mismatch
    # between number of properties and number of shapes" because the mixed
    # collection desyncs its shape array from the property table. Extract the
    # polygonal parts and union them per row so each divide keeps exactly one
    # (MULTI)POLYGON geometry and a 1:1 divide_id mapping.
    gtypes <- as.character(sf::st_geometry_type(d_pre_r))
    nonpoly <- which(!gtypes %in% c("POLYGON", "MULTIPOLYGON"))
    if (length(nonpoly) > 0L) {
      gvec <- sf::st_geometry(d_pre_r)
      for (i in nonpoly) {
        gvec[[i]] <- tryCatch(
          {
            poly <- suppressWarnings(
              sf::st_collection_extract(sf::st_sfc(gvec[[i]]), "POLYGON"))
            if (length(poly) == 0L) gvec[[i]]
            else sf::st_geometry(sf::st_union(poly))[[1L]]
          },
          error = function(e) gvec[[i]])
      }
      sf::st_geometry(d_pre_r) <- sf::st_sfc(as.list(gvec),
        crs = sf::st_crs(d_pre_r))
      message(sprintf("[stage2][%s] Step 4b-pre: coerced %d non-polygonal divide(s) to polygon (mapshaper-safe)",
        vpuid, length(nonpoly)))
    }

    d_pre_r <- .drop_list_cols(d_pre_r)
    sf::st_write(d_pre_r, out_gpkg, "divides", delete_layer = TRUE, quiet = TRUE)
    if (n_invalid_pre > 0L || n_repaired > 0L)
      message(sprintf("[stage2][%s] Step 4b-pre: normalised %d divide(s) (%d were invalid)",
        vpuid, nrow(d_pre_r), n_invalid_pre))
  })

  .probe_area_s2("after-4b-pre")
  .t("Step 4b: clean_gpkg_layer(divides)", {
    # Snapshot ids before clean so we can restore any divide that mapshaper
    # drops: those divides carry lineage-derived ids that flowpaths reference.
    div_before_clean <- sf::st_read(out_gpkg, "divides", quiet = TRUE)
    old_ogr_ct_opts <- Sys.getenv("OGR_CT_OPTS", unset = "")
    on.exit(Sys.setenv(OGR_CT_OPTS = old_ogr_ct_opts), add = TRUE)
    # Suppress the "several coordinate operations" informational warning from
    # GDAL reprojection steps inside clean_gpkg_layer. ALLOW_BALLPARK=NO /
    # ONLY_BEST=YES are intentionally omitted: they cause ogr2ogr to exit with
    # status 1 when PROJ lacks NADCON grid files for NAD83->WGS84, which is
    # common in conda/brew PROJ installs.
    Sys.setenv(OGR_CT_OPTS = "WARN_ABOUT_DIFFERENT_COORD_OP=NO")
    hfutils::clean_gpkg_layer(out_gpkg, layer = "divides",
      snap_interval = params$snap_interval,
      gap_fill_area = params$gap_fill_area,
      overlap_rule  = params$overlap_rule, verbose = FALSE)
    div_after_clean <- sf::st_read(out_gpkg, "divides", quiet = TRUE)
    dropped_ids <- setdiff(as.character(div_before_clean$divide_id),
      as.character(div_after_clean$divide_id))
    if (length(dropped_ids) > 0L) {
      # Only restore divides that are referenced by a flowpath's divide_id.
      # Orphan divides (fp_id=NA for every flowpath) come from all-disjoint
      # partition groups where multiple aggregated fp units share an identical
      # geometry. Restoring them re-introduces duplicate polygons that inflate
      # Stage 2 overlap.
      fp_out <- sf::st_drop_geometry(sf::st_read(out_gpkg, "flowpaths", quiet = TRUE))
      fp_div_ids <- as.character(fp_out$divide_id)
      fp_div_ids <- fp_div_ids[!is.na(fp_div_ids) & nzchar(fp_div_ids)]
      n_orphan <- sum(!dropped_ids %in% fp_div_ids)
      dropped_ids <- dropped_ids[dropped_ids %in% fp_div_ids]
      if (n_orphan > 0L)
        message(sprintf("[stage2][%s] Step 4b: skipped %d orphan divide(s) (no fp reference) -- not restored",
          vpuid, n_orphan))
    }
    if (length(dropped_ids) > 0L) {
      restored <- div_before_clean[
        as.character(div_before_clean$divide_id) %in% dropped_ids, , drop = FALSE]
      restored <- sf::st_transform(restored, sf::st_crs(div_after_clean))
      restored <- .drop_list_cols(restored)
      sf::st_write(restored, out_gpkg, "divides", append = TRUE, quiet = TRUE)
      message(sprintf("[stage2][%s] Step 4b: restored %d lineage divide(s) dropped by clean",
        vpuid, length(dropped_ids)))
    }
    # Net effect of the whole step, clean plus restore. This step reported only
    # its run time, so what it did to the fabric was unmeasured. Overlap removal
    # is its purpose and overlap is only visible as sum minus union.
    # Re-read only when a restore actually happened; otherwise div_after_clean is
    # already the final state. The extra read of a 36k-divide layer is not free.
    .report_divide_effect("Step 4b clean_gpkg_layer",
      .divide_metrics(div_before_clean),
      .divide_metrics(if (length(dropped_ids) > 0L)
        sf::st_read(out_gpkg, "divides", quiet = TRUE) else div_after_clean),
      vpuid = vpuid)
  })
  .probe_area_s2("after-4b")
  .t("Step 4c: fill_small_holes(divides)",
    fill_small_holes(out_gpkg, threshold_m2 = params$threshold_m2,
      layer = "divides"))

  # -- 4d. Normalize divide CRS + enforce POLYGON geometry ---------------
  .probe_area_s2("after-4c")
  .t("Step 4d: normalize divide CRS", {
    fp_crs <- sf::st_crs(sf::st_read(out_gpkg, "flowpaths", quiet = TRUE))
    for (lyr in c("divides", "landscape_divides")) {
      if (!hfutils::layer_exists(out_gpkg, lyr)) next
      d <- sf::st_read(out_gpkg, lyr, quiet = TRUE)
      changed <- FALSE
      if (sf::st_crs(d) != fp_crs) {
        d <- sf::st_transform(d, fp_crs)
        changed <- TRUE
      }
      # Cleaning can bisect a narrow divide into disjoint parts. Rejoin them into
      # one multipart divide rather than keeping a single part, which discarded
      # the rest and destroyed catchment area.
      .mb <- .merge_bisected_divides(d)
      if (.mb$n_merged > 0L)
        message(sprintf(
          "[stage2][%s] Step 4d: re-merged %d bisected divide(s), area-conserving",
          vpuid, .mb$n_merged))
      if (isTRUE(.mb$changed)) { d <- .mb$divides; changed <- TRUE }
      if (changed) {
        d <- .drop_list_cols(d)
        sf::st_write(d, out_gpkg, lyr, delete_layer = TRUE, quiet = TRUE)
      }
    }
  })

  # -- 4d2. Refresh derived geometry metrics ----------------------------
  # Geometry changed through clip/union/clean/hole-fill/cast; lengthkm and
  # areasqkm must be recomputed from final geometry before validation.
  fp_current <- .t("Step 4d2: refresh lengthkm/areasqkm", {
    fp <- sf::st_read(out_gpkg, "flowpaths", quiet = TRUE)
    if ("lengthkm" %in% names(fp)) {
      new_len <- as.numeric(sf::st_length(fp)) / 1000
      n_diff <- sum(abs(new_len - suppressWarnings(as.numeric(fp$lengthkm))) > 1e-6,
        na.rm = TRUE)
      fp$lengthkm <- new_len
      if (n_diff > 0L)
        message(sprintf("[stage2][%s] Step 4d2: refreshed lengthkm for %d flowpath(s)",
          vpuid, n_diff))
    }
    fp <- .drop_list_cols(fp)
    sf::st_write(fp, out_gpkg, "flowpaths", delete_layer = TRUE, quiet = TRUE)

    div <- sf::st_read(out_gpkg, "divides", quiet = TRUE)
    if ("areasqkm" %in% names(div)) {
      new_area <- as.numeric(sf::st_area(div)) / 1e6
      n_diff_a <- sum(abs(new_area - suppressWarnings(as.numeric(div$areasqkm))) > 1e-6,
        na.rm = TRUE)
      div$areasqkm <- new_area
      if (n_diff_a > 0L)
        message(sprintf("[stage2][%s] Step 4d2: refreshed areasqkm for %d divide(s)",
          vpuid, n_diff_a))
    }
    div <- .drop_list_cols(div)
    sf::st_write(div, out_gpkg, "divides", delete_layer = TRUE, quiet = TRUE)
    fp   # return the refreshed flowpaths for reuse by Step 4e
  })

  # -- 4d2b. Retract a divide claim the layer cannot honour -------------
  # has_divide is set before any divide exists, from the reference connector
  # flag, so it is a prediction. construct_agg_divides then silently skips
  # rec_ids it cannot find, and an aggregate whose members have no Stage-1
  # catchment gets no divide while the flag and divide_id still say it has one.
  # On VPU 04 that leaves 69 flowpaths naming a divide_id with no divide row,
  # eight of them also carrying an areasqkm taken from a different network row.
  #
  # Stage 3 nulls the dangling id much later (Step 6d5), so the product is
  # self-consistent while every intermediate is not, and the aggregate-grain
  # invariant cannot state the property. Reconcile it here, where both layers
  # are final and in hand.
  #
  # Claims are only ever retracted, never granted: a divide that exists but is
  # not claimed may belong to another flowpath, and promoting it would create
  # the shared-divide state flowpaths_no_shared_divide forbids.
  fp_current <- .t("Step 4d2b: retract unhonourable divide claims", {
    fp <- fp_current
    if (all(c("divide_id", "has_divide") %in% names(fp))) {
      dids <- as.character(sf::st_read(out_gpkg, "divides", quiet = TRUE)$divide_id)
      fdid <- as.character(fp$divide_id)
      claim <- fp$has_divide %in% c(TRUE, 1L, "1", "TRUE") &
        !is.na(fdid) & nzchar(fdid)
      dang  <- which(claim & !fdid %in% dids)
      if (length(dang)) {
        fp$has_divide[dang] <- FALSE
        fp$divide_id[dang]  <- NA_character_
        if ("areasqkm" %in% names(fp)) fp$areasqkm[dang] <- 0
        fp <- .drop_list_cols(fp)
        sf::st_write(fp, out_gpkg, "flowpaths", delete_layer = TRUE, quiet = TRUE)
        message(sprintf(paste0("[stage2][%s] Step 4d2b: %d flowpath(s) claimed a ",
          "divide with no divide row; claim retracted (has_divide=FALSE, ",
          "areasqkm=0)"), vpuid, length(dang)))
      }
    }
    fp
  })

  # Step 4d3 removed with the rest of the disjoint-repair cascade; see Step 3b.
  # fp_current stays as Step 4d2b left it.

  # -- 4e. Small-segment diagnostic layer --------------------------------
  .t("Step 4e: write small_segments layer", {
    # Reuse the in-memory flowpaths threaded from Step 4d2 instead of
    # re-reading the whole layer from disk (the Step 4e hotspot). Fall back to
    # a read if it is somehow unavailable.
    fp_all <- if (!is.null(fp_current)) fp_current else
      sf::read_sf(out_gpkg, "flowpaths", quiet = TRUE)
    small  <- .small_segments_layer(fp_all)
    sf::st_write(small, out_gpkg, "small_segments",
      delete_layer = TRUE, quiet = TRUE)
    message(sprintf("[stage2][%s] small_segments: %d total (%d poi_protected, %d terminal, %d other)",
      vpuid, nrow(small),
      sum(small$small_type == "poi_protected"),
      sum(small$small_type == "terminal"),
      sum(small$small_type == "other")))
  })

  # -- Stage 2 integrity gate ---------------------------------------------
  .t("Step 5: integrity_gate", {
    fp  <- sf::st_drop_geometry(sf::read_sf(out_gpkg, "flowpaths", quiet = TRUE))
    div <- sf::st_drop_geometry(sf::read_sf(out_gpkg, "divides",   quiet = TRUE))
    poi <- tryCatch(sf::st_drop_geometry(sf::read_sf(out_gpkg, "pois", quiet = TRUE)),
      error = function(e) NULL)

    fp_ids  <- as.character(fp$flowpath_id)
    div_ids <- as.character(div$divide_id)

    # 1. DAG: all flowpath_toids must resolve or be terminal.
    # Terminal convention: "0", NA, or numeric value >= 1e9 (outlet sentinel used by the aggregation step).
    toids     <- as.character(fp$flowpath_toid)
    toids_num <- suppressWarnings(as.numeric(toids))
    is_terminal <- is.na(toids) | toids == "0" | (!is.na(toids_num) & toids_num >= 1e9)
    bad_toids <- toids[!is_terminal & !toids %in% fp_ids]
    if (length(bad_toids) > 0L)
      .hf_log_warn(sprintf("[stage2][%s] %d flowpath(s) have toid not in flowpath_ids (dangling)",
        vpuid, length(bad_toids)))
    else
      message(sprintf("[stage2][%s] toid integrity OK (%d fps, %d terminal)",
        vpuid, nrow(fp), sum(is_terminal)))

    # 2. Every flowpath must have a matching divide_id
    fp_div_id <- as.character(fp$divide_id)
    orphan_fp <- sum(!is.na(fp_div_id) & !fp_div_id %in% div_ids)
    if (orphan_fp > 0L)
      .hf_log_warn(sprintf("[stage2][%s] %d flowpath(s) reference a non-existent divide_id",
        vpuid, orphan_fp))
    else
      message(sprintf("[stage2][%s] fp->div join OK (%d fps have valid divide_id)", vpuid,
        sum(!is.na(fp_div_id))))

    # 3. No orphan divides
    orphan_div <- sum(!div_ids %in% fp_div_id)
    if (orphan_div > 0L)
      .hf_log_warn(sprintf("[stage2][%s] %d divide(s) not referenced by any flowpath",
        vpuid, orphan_div))
    else
      message(sprintf("[stage2][%s] div->fp join OK (no orphan divides)", vpuid))

    # 4. POI flowpath_id remapping completeness
    if (!is.null(poi) && nrow(poi) > 0L) {
      poi_fp <- as.character(poi$flowpath_id)
      bad_poi <- sum(!is.na(poi_fp) & !poi_fp %in% fp_ids)
      na_poi  <- sum(is.na(poi_fp))
      if (bad_poi > 0L || na_poi > 0L)
        .hf_log_warn(sprintf("[stage2][%s] POI remap: %d with invalid fp_id, %d with NA fp_id",
          vpuid, bad_poi, na_poi))
      else
        message(sprintf("[stage2][%s] POI remap OK (%d pois all reference valid flowpaths)",
          vpuid, nrow(poi)))
    }

    # 5. No MULTILINESTRING/MULTIPOLYGON remnants
    fp_geom  <- sf::read_sf(out_gpkg, "flowpaths", quiet = TRUE)
    div_geom <- sf::read_sf(out_gpkg, "divides",   quiet = TRUE)
    n_multi_fp  <- sum(sf::st_geometry_type(fp_geom)  == "MULTILINESTRING")
    n_multi_div <- sum(sf::st_geometry_type(div_geom) == "MULTIPOLYGON")
    if (n_multi_fp > 0L)
      .hf_log_warn(sprintf("[stage2][%s] %d MULTILINESTRING flowpath(s) remain",
        vpuid, n_multi_fp))
    if (n_multi_div > 0L)
      .hf_log_warn(sprintf("[stage2][%s] %d MULTIPOLYGON divide(s) remain",
        vpuid, n_multi_div))
    if (n_multi_fp == 0L && n_multi_div == 0L)
      message(sprintf("[stage2][%s] geometry types OK (all LINESTRING / POLYGON)", vpuid))

    # 6. areasqkm recompute: update any divide whose stored area deviates by >1%
    div_geom_area <- as.numeric(sf::st_area(div_geom)) / 1e6
    stored_area   <- div_geom$areasqkm
    if (!is.null(stored_area)) {
      stale <- abs(div_geom_area - stored_area) / pmax(stored_area, 1e-9) > 0.01
      stale[is.na(stale)] <- FALSE
      if (any(stale)) {
        stale_ids  <- as.character(div_geom$divide_id[stale])
        stale_area <- div_geom_area[stale]
        hfutils::gpkg_update_col(out_gpkg, "divides", "divide_id",
          stale_ids, "areasqkm", stale_area)
        message(sprintf("[stage2][%s] areasqkm refreshed for %d divide(s) (>1%% stale)",
          vpuid, sum(stale)))
      } else {
        message(sprintf("[stage2][%s] areasqkm OK (all divides within 1%% of geometry)",
          vpuid))
      }
    }

    # 7. Minimum length / area floors. A network flowpath should not be a
    # degenerate stub: anything shorter than min_fragment_m that is NOT a
    # legitimately-protected short segment (POI-pinned or terminal, which are
    # tracked separately in small_segments) indicates geometry loss during
    # aggregation. Likewise a network divide with near-zero area is a sliver.
    # These are loud diagnostics, not hard stops, because the upstream
    # fragmentation that produces them is a known spectrum; the gate surfaces
    # them per-VPU for triage.
    fp_geom_len_m <- as.numeric(sf::st_length(fp_geom))
    poi_col       <- if ("poi_id" %in% names(fp)) as.character(fp$poi_id) else rep(NA_character_, nrow(fp))
    protected     <- is_terminal | (!is.na(poi_col) & nzchar(poi_col))
    tiny_fp <- which(fp_geom_len_m < params$min_fragment_m & !protected)
    if (length(tiny_fp) > 0L)
      .hf_log_warn(sprintf(
        "[stage2][%s] min-length: %d unprotected flowpath(s) below %.0fm (degenerate stub) -- e.g. %s",
        vpuid, length(tiny_fp), params$min_fragment_m,
        paste(utils::head(sprintf("%s(%.1fm)", fp_ids[tiny_fp],
          fp_geom_len_m[tiny_fp]), 5), collapse = ", ")))
    else
      message(sprintf("[stage2][%s] min-length OK (all unprotected flowpaths >= %.0fm)",
        vpuid, params$min_fragment_m))

    min_div_area_km2 <- 1e-4  # 100 m2 -- below this is a degenerate sliver
    tiny_div <- which(div_geom_area < min_div_area_km2)
    if (length(tiny_div) > 0L)
      .hf_log_warn(sprintf(
        "[stage2][%s] min-area: %d divide(s) below %.0f m2 (sliver) -- e.g. %s",
        vpuid, length(tiny_div), min_div_area_km2 * 1e6,
        paste(utils::head(sprintf("%s(%.1fm2)", as.character(div_geom$divide_id[tiny_div]),
          div_geom_area[tiny_div] * 1e6), 5), collapse = ", ")))
    else
      message(sprintf("[stage2][%s] min-area OK (all divides >= %.0f m2)",
        vpuid, min_div_area_km2 * 1e6))
  })

  # -- Stage 2 invariant check --------------------------------------------
  agg_inv <- .t("Step 6: check_invariants(aggregated)", {
    fp_inv  <- sf::read_sf(out_gpkg, "flowpaths", quiet = TRUE)
    div_inv <- sf::read_sf(out_gpkg, "divides",   quiet = TRUE)
    net_inv <- tryCatch(as.data.frame(sf::st_read(out_gpkg, "network", quiet = TRUE)),
      error = function(e) NULL)
    check_invariants("aggregated",
      flowpaths = fp_inv, divides = div_inv, network = net_inv,
      strict = FALSE, coverage = FALSE,
      area = TRUE, area_baseline_m2 = baseline_div_area_m2,
      coverage_baseline_m2 = baseline_div_cov_m2)
  })

  files <- list(
    out_gpkg = out_gpkg,
    timings  = file.path(work_dir, "timings.json"))

  jsonlite::write_json(timings, files$timings, auto_unbox = TRUE,
    pretty = TRUE)

  invisible(list(files = unlist(files), timings = timings,
    invariants = agg_inv))
}


# Build the small-segment diagnostic subset (Step 4e helper).
# Selects flowpaths below the size thresholds and tags each by why it is small.
# @noRd
.small_segments_layer <- function(fp) {
  small <- fp[fp$areasqkm < 3 | fp$lengthkm < 0.25, ]
  small$small_type <- ifelse(
    !is.na(small$poi_id), "poi_protected",
    ifelse(as.numeric(small$flowpath_toid) >= 1e9, "terminal", "other"))
  small
}
