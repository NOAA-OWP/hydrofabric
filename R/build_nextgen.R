#' Build the NextGen data model from refactored + aggregated layers (nextgen phase)
#'
#' Consumes stage-1 (refactored + reconciled + hydrolocations) AND
#' stage-2 (`out_vpu{id}.gpkg`) artifacts; produces the ngen-ready
#' GeoPackage with model attributes and aspatial network.
#'
#' `ngen_vpu{id}.gpkg` writes to `{out_dir}/vpu_{vpuid}/`; all other
#' files write to `{out_dir}/vpu_{vpuid}/work/`:
#' - `ngen_vpu{id}.gpkg`: layers: `flowpaths`, `divides`, `nexus`,
#'   plus aspatial network and model-attribute layers
#' - `timings.json`, `invariants.json`
#'
#' @param vpuid VPU identifier.
#' @param fabric_dir Path to the vpu working directory holding the refactored +
#'   aggregated fabric.
#' @param fdr_path,fac_path Raster paths (local or `/vsis3/...`).
#' @param fallback_fdr_path,fallback_fac_path Optional paths to a fallback
#'   FDR/FAC pair used when the primary returns NA.  `NULL` disables the
#'   fallback (default).  Ignored when equal to the primary paths.
#' @param out_dir Parent output directory; `vpu_{vpuid}/` is created.
#' @param coverage_min Divide-coverage invariant + coalesce/clip floor (0-1).
#' @param threshold_m2 `fill_small_holes` area threshold on ngen divides (m^2).
#' @param threshold_sqkm `prune_ghost_divides` area threshold (km^2).
#' @param gap_fill_area Gap-fill area for `hfutils::clean_gpkg_layer` (ngen phase).
#' @param overlap_rule Divide overlap resolution rule.
#' @param n_cores Worker count for `repair_split_divides_at_pois`.
#' @param hydraulics_path Optional path to a RiverML channel table parquet,
#'   either its v2 output (e.g. `superconus_ml_attributes_v2.parquet`) or the
#'   earlier channel-geometry table.
#'   When supplied, Step 6g computes flow attributes for both `flowpath_id` and
#'   `flowline_id` and writes `flowpath-attributes` and `flowline-attributes`
#'   layers to the ngen GeoPackage. Reaches with no match in the table fall
#'   back through a geomorphic proxy. `NULL` skips Step 6g (default).
#' @param lakes_path Optional path to the NWM lakes GeoPackage
#'   (e.g. `nwm_lakes.gpkg`).  When supplied, Step 6h identifies the
#'   outlet flowline of each lake (minimum `flowline_hydroseq` among
#'   COMIDs belonging to the lake), joins NWM lake attributes, and
#'   writes a `lakes` layer to the ngen GeoPackage.  The reach -> `lake_id`
#'   membership used for outlet fallback and member stamping comes from
#'   `lake_mask_path`.  `NULL` skips Step 6h (default).
#' @param lake_mask_path Optional path to a lake-mask parquet: reach -> waterbody
#'   membership with integer columns `comid` (reference flowline id) and
#'   `wbareacomi` (waterbody / `lake_id`). Required when `lakes_path` is
#'   supplied. The CONUS mask is materialized from NHDPlusV2 wbareacomi;
#'   oCONUS / superconus domains supply an equivalent with the same schema.
#' @param slope_parquet Optional path to a per-reference-id slope parquet,
#'   passed to [finalize_data_model()]. `NULL` (default) takes slope from
#'   `hydraulics_path` when that table carries a `slope_m_m` column, and
#'   otherwise leaves `slope` `NA`.
#' @param reference_gpkg Optional path to the reference fabric GeoPackage (source
#'   of member slope / divide geometry). `NULL` disables features needing it.
#' @param verbose Propagated to internal helpers.
#' @return Invisibly, a list with `files`, `timings`, `invariants`.
#' @examples
#' \dontrun{
#' nextgen_fabric(vpuid      = "01",
#'                fabric_dir = "out/vpu_01",
#'                fdr_path   = "fdr.tif",
#'                fac_path   = "fac.tif",
#'                out_dir    = "out",
#'                coverage_min = 0.95)
#' }
#' @export
nextgen_fabric <- function(vpuid, fabric_dir,
                           fdr_path, fac_path,
                           fallback_fdr_path = NULL,
                           fallback_fac_path = NULL,
                           out_dir,
                           coverage_min   = 0.90,
                           threshold_m2   = 1000,
                           threshold_sqkm = 0.001,
                           gap_fill_area  = 1000,
                           overlap_rule   = "max-area",
                           n_cores = NULL,
                           hydraulics_path = NULL,
                           lakes_path = NULL,
                           lake_mask_path = NULL,
                           slope_parquet   = NULL,
                           reference_gpkg  = NULL,
                           verbose = FALSE) {
  params <- list(coverage_min   = coverage_min,
                 threshold_m2   = threshold_m2,
                 threshold_sqkm = threshold_sqkm,
                 gap_fill_area  = gap_fill_area,
                 overlap_rule   = overlap_rule)
  stopifnot(
    is.character(vpuid), length(vpuid) == 1L, nzchar(vpuid),
    dir.exists(fabric_dir),
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
    message(sprintf("[stage3][%s][%5.1fs] %s",
      vpuid, timings[[label]], label))
    result
  }

  # -- Stage-1 inputs -------------------------------------------------
  s1_refactored <- file.path(fabric_dir, "work", "refactored.gpkg")
  s1_reconciled <- file.path(fabric_dir, "work", "reconciled.gpkg")
  s1_hl         <- file.path(fabric_dir, "work", "hydrolocations.gpkg")
  stopifnot(file.exists(s1_refactored), file.exists(s1_reconciled),
    file.exists(s1_hl))

  # Layer name: current fast_reactor.R uses filename-derived names ("refactored",
  # "reconciled"); older checkpoints used "flowpaths". Try both.
  .s3_read_ckpt <- function(gpkg, preferred, fallback = "flowpaths") {
    lyr <- if (hfutils::layer_exists(gpkg, preferred)) preferred else fallback
    sf::read_sf(gpkg, lyr, quiet = TRUE)
  }

  refactored_sf <- .s3_read_ckpt(s1_refactored, "refactored")
  reconciled_sf <- .s3_read_ckpt(s1_reconciled, "reconciled")
  # Normalize: rec_id/rec_toid (new schema) -> ID/toID (add_flowlines internal convention).
  names(reconciled_sf)[names(reconciled_sf) == "rec_id"]   <- "ID"
  names(reconciled_sf)[names(reconciled_sf) == "rec_toid"] <- "toID"
  events_sf     <- if (hfutils::layer_exists(s1_reconciled, "events"))
    sf::read_sf(s1_reconciled, "events", quiet = TRUE) else NULL
  hl_sf         <- if (hfutils::layer_exists(s1_hl, "hydrolocations"))
    sf::read_sf(s1_hl, "hydrolocations", quiet = TRUE) else NULL

  # Assert reconciled network is a DAG before add_flowlines builds sub-flowline
  # topology. Cycles here signal a stage-1 bug; break is emergency fallback only.
  {
    dag_v <- hfutils::hf_assert_network_dag(reconciled_sf, id_col = "ID", toid_col = "toID")
    if (!dag_v$is_dag) {
      .hf_log_warn(sprintf(
        "[stage3][%s] BUG: cycles in reconciled_sf (%s) -- stage-1 toid propagation is incorrect; breaking as fallback",
        vpuid, dag_v$message))
      reconciled_sf <- hfutils::hf_break_cycles(reconciled_sf, id_col = "ID", toid_col = "toID")
    } else {
      message(sprintf("[stage3][%s] reconciled_sf DAG OK", vpuid))
    }
    reconciled_sf <- hfutils::hf_recompute_hydroseq(reconciled_sf, id_col = "ID", toid_col = "toID")
    message(sprintf("[stage3][%s] reconciled_sf hydroseq recomputed (%d fps)",
      vpuid, nrow(reconciled_sf)))
  }

  # -- Stage-2 inputs -------------------------------------------------
  out_gpkg <- file.path(fabric_dir, "work", sprintf("out_vpu%s.gpkg", vpuid))
  if (!file.exists(out_gpkg))
    stop("stage2 output missing: ", out_gpkg)

  # Copy the stage-2 gpkg into our working dir so we don't mutate the
  # upstream artifact. Stage-3 modifies out_gpkg in place during cleanup.
  # Skip when fabric_dir == vpu_dir (local run_pipeline where all three
  # stages share a directory: mutating is fine there).
  out_gpkg_local <- file.path(work_dir, basename(out_gpkg))
  if (normalizePath(out_gpkg, mustWork = FALSE) !=
    normalizePath(out_gpkg_local, mustWork = FALSE)) {
    file.copy(out_gpkg, out_gpkg_local, overwrite = TRUE)
    out_gpkg <- out_gpkg_local
  }

  # Baseline drainage area to conserve through the ngen build: total area of
  # the Stage-2 aggregated network divides. The closing invariant check asserts
  # the ngen divides stay within tolerance of it.
  # Both a summed and a UNION baseline: the sum guards against area appearing or
  # vanishing, the union guards against ground being lost. Only the union can see
  # a hole, since removing an overlap also reduces the sum.
  .baseline_div <- tryCatch(
    .hf_network_divides(sf::read_sf(out_gpkg, "divides", quiet = TRUE)),
    error = function(e) NULL)
  baseline_div_area_m2 <- tryCatch(
    sum(suppressWarnings(as.numeric(sf::st_area(.baseline_div))), na.rm = TRUE),
    error = function(e) NA_real_)
  baseline_div_cov_m2 <- tryCatch(
    as.numeric(sf::st_area(sf::st_union(sf::st_geometry(.baseline_div)))),
    error = function(e) NA_real_)
  rm(.baseline_div)

  ngen_gpkg_final <- file.path(vpu_dir, sprintf("ngen_vpu%s.gpkg", vpuid))
  if (file.exists(ngen_gpkg_final)) file.remove(ngen_gpkg_final)

  # Use a local (APFS) temp path for all intermediate ngen_gpkg operations.
  # This avoids NTFS SQLite locking issues where GDAL can't open the file
  # for UPDATE after a cross-device file.copy from write_hydrofabric's tmpfile.
  # write_hydrofabric will also benefit: tmpfile -> same-device rename succeeds.
  # The PID makes the path process-unique so forked parallel workers (e.g.
  # concurrent gage runs of the same VPU via mclapply) never clobber each
  # other's working ngen gpkg in the shared tempdir().
  ngen_gpkg <- file.path(tempdir(), sprintf("ngen_vpu%s_p%d_work.gpkg", vpuid, Sys.getpid()))
  if (file.exists(ngen_gpkg)) file.remove(ngen_gpkg)

  # DIAGNOSTIC probe (gated on HF_PROBE env): logs how many fps have a
  # divide_id whose polygon covers <90% of the fp line, after a given step.
  # Used to localize which stage-3 divide step strips coverage.
  .probe_cov <- function(tag) {
    if (!nzchar(Sys.getenv("HF_PROBE"))) return(invisible())
    tryCatch(
      {
        if (!file.exists(ngen_gpkg) || !hfutils::layer_exists(ngen_gpkg, "divides")) return(invisible())
        fpp <- sf::st_read(ngen_gpkg, "flowpaths", quiet = TRUE)
        dvp <- sf::st_read(ngen_gpkg, "divides",   quiet = TRUE)
        has <- !is.na(fpp$divide_id) & nchar(as.character(fpp$divide_id)) > 0L
        fps <- fpp[has, ]
        di  <- match(as.character(fps$divide_id), as.character(dvp$divide_id))
        m   <- !is.na(di)
        lg  <- geos::as_geos_geometry(sf::st_transform(sf::st_geometry(fps)[m], sf::st_crs(dvp)))
        dg  <- geos::as_geos_geometry(sf::st_geometry(dvp)[di[m]])
        ln  <- as.numeric(geos::geos_length(lg))
        ins <- suppressWarnings(as.numeric(geos::geos_length(geos::geos_intersection(lg, dg))))
        ins[is.na(ins)] <- 0
        fr  <- ifelse(ln > 0, ins / ln, NA_real_)
        # ORPHAN divides: rows in the divides layer that no flowpath references.
        # A sever-while-keep bug (a step clears an fp's divide_id but leaves the
        # divide row) creates orphans; tracking the orphan count per tag pins the
        # step that creates them.
        used    <- unique(as.character(fpp$divide_id[has]))
        orphans <- sum(!(as.character(dvp$divide_id) %in% used))
        # Coverage accounting (HF_PROBE_AREA, separate because the union is slow).
        # Sum-of-areas cannot separate overlap removal from a hole; the UNION can,
        # so tracking union per tag pins the step that actually loses coverage.
        if (nzchar(Sys.getenv("HF_PROBE_AREA"))) {
          .s <- sum(suppressWarnings(as.numeric(sf::st_area(dvp))), na.rm = TRUE) / 1e6
          .u <- tryCatch(as.numeric(sf::st_area(sf::st_union(sf::st_geometry(dvp)))) / 1e6,
            error = function(e) NA_real_)
          message(sprintf(
            "[AREA  %-18s][%s] divs=%d | sum=%.4f | union=%.4f | overlap=%.4f",
            tag, vpuid, nrow(dvp), .s, .u, .s - .u))
        }
        # Divide-less REAL-catchment fps (area>0) = the invariant's "unexplained gap"
        # population. Tracking this per substep pins which step destroys real divides.
        ar      <- if ("areasqkm" %in% names(fpp)) suppressWarnings(as.numeric(fpp$areasqkm)) else rep(NA_real_, nrow(fpp))
        gap_rc  <- sum(!has & !is.na(ar) & ar > 0)
        message(sprintf(
          "[PROBE %-18s][%s] fps=%d divs=%d w/div=%d matched=%d | cov<0.90=%d | dangling=%d | ORPHANS=%d | GAP_realcat=%d | med_cov=%.3f",
          tag, vpuid, nrow(fpp), nrow(dvp), sum(has), sum(m), sum(fr < 0.90, na.rm = TRUE),
          sum(!m), orphans, gap_rc, stats::median(fr, na.rm = TRUE)))
      },
      error = function(e) message(sprintf("[PROBE %s] error: %s", tag, conditionMessage(e))))
    invisible()
  }

  # -- 5. Flowlines + cleanup -----------------------------------------
  # Coverage probe on out_gpkg (HF_PROBE_AREA). The .probe_cov() probes read
  # ngen_gpkg, which does not exist until Step 6, so the whole Steps 5-5d window
  # is invisible to them, and that window is where the coverage loss actually
  # happens. This one reads out_gpkg, which those steps mutate in place.
  .probe_area_out <- function(tag) {
    if (!nzchar(Sys.getenv("HF_PROBE_AREA"))) return(invisible())
    tryCatch({
      d <- .hf_network_divides(sf::st_read(out_gpkg, "divides", quiet = TRUE))
      s <- sum(suppressWarnings(as.numeric(sf::st_area(d))), na.rm = TRUE) / 1e6
      u <- as.numeric(sf::st_area(sf::st_union(sf::st_geometry(d)))) / 1e6
      message(sprintf(
        "[AREAOUT %-14s][%s] divs=%d | sum=%.4f | union=%.4f | overlap=%.4f",
        tag, vpuid, nrow(d), s, u, s - u))
    }, error = function(e)
      message(sprintf("[AREAOUT %s] error: %s", tag, conditionMessage(e))))
    invisible()
  }
  .probe_area_out("stage3-entry")

  # DIAGNOSTIC fast-replay (gated on HF_SKIP_ADDFL): steps 5-5d only mutate
  # out_gpkg; when out_gpkg already carries their output ("flowlines" layer) a
  # rerun can skip the expensive add_flowlines/FAC-sampling phase and go straight
  # to the divide cascade (Step 6+). Safe no-op unless the env var is set.
  .skip_addfl <- nzchar(Sys.getenv("HF_SKIP_ADDFL")) &&
    hfutils::layer_exists(out_gpkg, "flowlines")
  if (.skip_addfl) {
    message(sprintf("[stage3][%s] HF_SKIP_ADDFL: reusing existing out_gpkg flowlines; skipping Steps 5-5d", vpuid))
  } else {
    .t("Step 5: add_flowlines",
      add_flowlines(flowlines = reconciled_sf, input_gpkg = out_gpkg,
        fac          = terra::rast(fac_path),
        fallback_fac = fallback_fac_path,
        poi_events   = events_sf))
    .probe_area_out("after-5")
    # Step 5a2 removed: the junction rule moved to stage 2 (Step 4_junction), which
    # runs before construct_agg_divides. Here it had to re-dissolve divides that
    # 4a2 and 4b had already corrected, so it could not conserve area exactly
    # (+14.09 km2 transient on VPU 01). In stage 2 the divide follows from the
    # corrected member partition on first construction.
    .t("Step 5b: extend_flowpaths_to_mainstem",
      extend_flowpaths_to_mainstem(out_gpkg))
    .probe_cov("after-5-add_flowlines-5b")
    .probe_area_out("after-5b")
    .t("Step 5b2: resolve_disjoint_flowpath_divides",
      tryCatch(
        resolve_disjoint_flowpath_divides(out_gpkg),
        error = function(e) {
          .hf_log_warn(sprintf("[stage3][%s] Step 5b2 non-fatal error: %s", vpuid,
            conditionMessage(e)))
          invisible(NULL)
        }))
    .probe_area_out("after-5b2")
    .t("Step 5c: coalesce_fragmented_flowpaths",
      coalesce_fragmented_flowpaths(out_gpkg,
        coverage_threshold = params$coverage_min))
    .probe_area_out("after-5c")
    .t("Step 5c-dag: dag_check + hydroseq after coalesce",
      .hf_enforce_dag_gpkg(out_gpkg, layer = "flowpaths", vpuid = vpuid))
    # A flowpath's line is NOT clipped to its own divide. An aggregated line is a
    # mainstem walk through its members and is not required to lie inside the
    # catchment: where it runs outside, so do the member reaches themselves, so
    # clipping discards real routing geometry rather than trimming an overshoot.
    # Containment is asserted by divide_has_member_flowline and the routeable
    # flag, not by the line.
  }

  # Step 5e removed: it enforced one divide per flowpath by demoting the
  # non-owner when two flowpaths referenced the same divide_id. Its input came
  # from Step 5b2's Group-B swap, which handed a disjoint flowpath a divide that
  # another flowpath already owned; the ownership guard added to that swap removed
  # the only source. Measured on all 10 VPUs where 5b2 and 5e ever fired, every
  # one now reports "no shared flowpath divides", and the other 14 never fired it.
  # The condition is asserted by flowpaths_no_shared_divide and
  # divide_flowpath_consistent, both gates, so a shared divide fails the build
  # instead of being silently repaired.

  # -- 6. NextGen -----------------------------------------------------
  .t("Step 6: nextgen_from_gpkg", nextgen_from_gpkg(out_gpkg, outfile = ngen_gpkg))
  local({
    con <- DBI::dbConnect(RSQLite::SQLite(), ngen_gpkg)
    on.exit(DBI::dbDisconnect(con))
    DBI::dbExecute(con, "REINDEX")
  })

  if (hfutils::layer_exists(out_gpkg, "landscape_divides")) {
    .t("Step 6a: append landscape divides", {
      lnd    <- sf::st_read(out_gpkg, "landscape_divides", quiet = TRUE)
      ng_div <- sf::st_read(ngen_gpkg, "divides", quiet = TRUE)
      lnd <- sf::st_transform(lnd, sf::st_crs(ng_div))
      # NA-safe area filter; may leave zero rows (no qualifying landscape divide).
      lnd <- lnd[!is.na(lnd$areasqkm) & lnd$areasqkm >= 0.001, ]
      # Only append when there is something to add: assigning columns to a
      # zero-row sf ("lnd[[nm]] <- NA") errors with "replacement has 1 row, data
      # has 0". An empty landscape set is a no-op; leave ngen divides unchanged.
      if (nrow(lnd) > 0L) {
        for (nm in setdiff(names(ng_div), names(lnd))) lnd[[nm]] <- NA
        for (nm in setdiff(names(lnd), names(ng_div))) ng_div[[nm]] <- NA
        lnd <- lnd[, names(ng_div)]
        combined <- rbind(ng_div, lnd)
        combined <- .drop_list_cols(combined)
        sf::st_write(combined, ngen_gpkg, "divides",
          delete_layer = TRUE, quiet = TRUE)
      } else {
        message(sprintf(
          "[stage3][%s] Step 6a: no landscape divides >= 0.001 km2 to append", vpuid))
      }
    })
  }

  # Step 6b0: Remove duplicate-geometry divides before mapshaper sees them.
  # Stage 2 can assign one catchment polygon to several consecutive mainstem
  # flowpaths. mapshaper's -clean would drop the duplicates and leave those
  # flowpaths divide-less, so the groups are collapsed here first: the flowpath
  # whose line best intersects the shared polygon keeps it, the rest are cleared.
  .probe_cov("after-5b2")
  .t("Step 6b0: dedup identical-geometry divides", local({
    div <- sf::st_read(ngen_gpkg, "divides", quiet = TRUE)
    fp  <- sf::st_read(ngen_gpkg, "flowpaths", quiet = TRUE)
    if (nrow(div) == 0) return(invisible(NULL))

    div_t  <- sf::st_transform(div, 5070)
    area_r <- round(as.numeric(sf::st_area(div_t)) / 100)  # 100 m2 bins
    cents  <- suppressWarnings(sf::st_centroid(div_t))
    cx     <- round(sf::st_coordinates(cents)[, 1L] / 100)  # 100 m grid
    cy     <- round(sf::st_coordinates(cents)[, 2L] / 100)
    fp_key <- paste(area_r, cx, cy, sep = "_")

    key_tab  <- table(fp_key)
    dup_keys <- names(key_tab[key_tab > 1L])
    if (length(dup_keys) == 0L) return(invisible(NULL))

    fp_did  <- as.character(fp$divide_id)
    div_did <- as.character(div$divide_id)
    fp_t    <- sf::st_transform(fp, 5070)

    remove_cats <- character(0)
    for (dk in dup_keys) {
      grp_idx  <- which(fp_key == dk)
      grp_div  <- div_t[grp_idx, , drop = FALSE]
      grp_cats <- as.character(grp_div$divide_id)

      # Representative polygon (group members are ~identical)
      rep_div_geom <- sf::st_geometry(grp_div)[[1L]]

      # Intersection length of each fp through the shared polygon
      inter_len <- vapply(grp_cats, function(cat) {
        m <- which(fp_did == cat)
        if (length(m) == 0L) return(0)
        fp_row <- fp_t[m[1L], , drop = FALSE]
        inter  <- suppressWarnings(sf::st_intersection(
          sf::st_geometry(fp_row)[[1L]], rep_div_geom))
        if (length(inter) == 0L || inherits(inter, c("POINT", "MULTIPOINT")))
          return(0)
        as.numeric(sf::st_length(inter))
      }, numeric(1L))

      winner   <- which.max(inter_len)
      remove_cats <- c(remove_cats, grp_cats[-winner])
    }

    if (length(remove_cats) == 0L) return(invisible(NULL))

    # Clear divide_id for fps that referenced removed cats
    clear_fp_ids <- as.character(fp$flowpath_id[fp_did %in% remove_cats])
    if (length(clear_fp_ids) > 0L) {
      hfutils::gpkg_exec(ngen_gpkg,
        list(drop_triggers_for = "flowpaths"),
        list(table = "._dedup_clr",
          df    = data.frame(`._id` = clear_fp_ids,
            stringsAsFactors = FALSE)),
        'UPDATE "flowpaths" SET "divide_id" = NULL, "has_divide" = 0
         WHERE CAST("flowpath_id" AS TEXT) IN (SELECT "._id" FROM "._dedup_clr")')
    }

    # Remove duplicate divides
    keep <- !(div_did %in% remove_cats)
    div_out <- div[keep, , drop = FALSE]
    div_out <- .drop_list_cols(div_out)
    sf::st_write(div_out, ngen_gpkg, "divides",
      delete_layer = TRUE, quiet = TRUE)

    message(sprintf(
      "[stage3][%s] Step 6b0: removed %d duplicate-geometry divide(s) across %d group(s)",
      vpuid, length(remove_cats), length(dup_keys)))
  }))
  .probe_cov("after-6b0")

  .t("Step 6b: fill_small_holes(ngen divides)",
    fill_small_holes(ngen_gpkg, threshold_m2 = params$threshold_m2,
      layer = "divides"))
  .probe_cov("after-6b")
  .t("Step 6c: prune_ghost_divides",
    prune_ghost_divides(ngen_gpkg,
      threshold_sqkm = params$threshold_sqkm))
  .probe_cov("after-6c")
  # The pre-clean area snapshot that fed Step 6d1 is gone with it: it existed only
  # so 6d1 could spot divides the clean pass had collapsed, and nothing collapses
  # now. It was a full read of the divides layer per VPU.
  .probe_cov("after-6b0-6b-6c")
  # Gated geometry snapshot around the clean, so the area a clean pass removes can
  # be attributed: absorbed by a neighbour (coverage intact) or left as a hole
  # (coverage lost). Sum-of-areas cannot tell those apart. Off unless HF_SNAP_6D.
  .snap6d <- function(tag) {
    if (!nzchar(Sys.getenv("HF_SNAP_6D"))) return(invisible())
    tryCatch({
      d <- sf::st_read(ngen_gpkg, "divides", quiet = TRUE)
      out <- file.path(work_dir, sprintf("snap6d_%s.gpkg", tag))
      sf::st_write(.drop_list_cols(d), out, layer = "divides",
        delete_dsn = TRUE, quiet = TRUE)
      message(sprintf("[stage3][%s] snap6d %-6s %d divides | sum %.4f km2 | -> %s",
        vpuid, tag, nrow(d),
        sum(suppressWarnings(as.numeric(sf::st_area(d))), na.rm = TRUE) / 1e6,
        basename(out)))
    }, error = function(e)
      message(sprintf("[stage3][%s] snap6d %s failed: %s", vpuid, tag,
        conditionMessage(e))))
    invisible()
  }
  .snap6d("pre")
  .t("Step 6d: clean_gpkg_layer(ngen divides)", {
    old_ogr_ct_opts <- Sys.getenv("OGR_CT_OPTS", unset = "")
    on.exit(Sys.setenv(OGR_CT_OPTS = old_ogr_ct_opts), add = TRUE)
    # Force a single high-quality coordinate operation in GDAL reprojection
    # steps used internally by clean_gpkg_layer.
    Sys.setenv(OGR_CT_OPTS = paste(
      "ALLOW_BALLPARK=NO",
      "ONLY_BEST=YES",
      "WARN_ABOUT_DIFFERENT_COORD_OP=NO"
    ))
    # This step reported only its run time, so its effect on the fabric was
    # unmeasured. Overlap removal is its stated purpose, and overlap is only
    # visible as sum minus union, so the union is worth the cost here.
    .eff_before <- .divide_metrics(
      sf::st_read(ngen_gpkg, "divides", quiet = TRUE))
    hfutils::clean_gpkg_layer(ngen_gpkg, layer = "divides",
      gap_fill_area = params$gap_fill_area,
      overlap_rule  = params$overlap_rule, verbose = FALSE)
    .report_divide_effect("Step 6d clean_gpkg_layer", .eff_before,
      .divide_metrics(sf::st_read(ngen_gpkg, "divides", quiet = TRUE)),
      vpuid = vpuid)
  })
  .snap6d("post")
  local({
    d <- sf::st_read(ngen_gpkg, "divides", quiet = TRUE)
    changed <- FALSE
    # Repair invalid geometries
    invalid <- !sf::st_is_valid(d)
    if (any(invalid, na.rm = TRUE)) {
      message(sprintf("[stage3][%s] repairing %d invalid divide(s) after clean_gpkg_layer",
        vpuid, sum(invalid, na.rm = TRUE)))
      d[invalid, ] <- sf::st_make_valid(d[invalid, ])
      changed <- TRUE
    }
    # Cleaning can bisect a narrow divide into disjoint parts. Rejoin them into
    # one multipart divide instead of keeping a single part, which discarded the
    # rest and destroyed catchment area.
    .mb <- .merge_bisected_divides(d)
    if (.mb$n_merged > 0L)
      message(sprintf(
        "[stage3][%s] Step 6d: re-merged %d bisected divide(s), area-conserving",
        vpuid, .mb$n_merged))
    if (isTRUE(.mb$changed)) { d <- .mb$divides; changed <- TRUE }
    if (changed) {
      d <- .drop_list_cols(d)
      sf::st_write(d, ngen_gpkg, "divides", delete_layer = TRUE, quiet = TRUE)
    }
  })
  # Step 6d1: clean_gpkg_layer resolves a fully-overlapping (duplicate / nested)
  # divide by assigning its area to the container, but leaves the loser as a
  # near-zero boundary sliver: above the 6d6 1 m2 floor, so it would survive as
  # a 0-area-divide "connector with a ghost divide". Detect divides whose area
  # collapsed from meaningful (>1000 m2) to negligible (<10 m2) across the clean
  # step, drop them, and mark their flowpaths as proper connectors (the catchment
  # now lives in the container divide: area is conserved, just re-owned).
  # Step 6d1 removed: it dropped divides that clean_gpkg_layer had collapsed.
  # Across all 24 VPUs on current code it dropped none. The collapse it cleaned
  # up after came from divides arriving overlapped; stage 2 now hands over a
  # clean tiling, so nothing collapses and there is nothing to drop.
  .probe_cov("after-6d-6d1")
  # Step 6d2 removed: a second clip_flowpath_lines_to_divides pass, here because
  # the cleanup steps between 5d and this point could re-break line coverage of a
  # divide. On all 24 VPUs it reported "no failing fps to clip": 5d clips 9 to 262
  # flowpaths per VPU and nothing after it reintroduces a failing one. The condition
  # is reported by the flowpath_coverage measurement and the pairing itself is gated
  # by divide_has_member_flowline. The 5d pass is unchanged.
  .t("Step 6d3: repair_split_divides_at_pois (post-cleanup)",
    repair_split_divides_at_pois(ngen_gpkg,
      fdr_path          = fdr_path,
      fac_path          = fac_path,
      fallback_fdr_path = fallback_fdr_path,
      fallback_fac_path = fallback_fac_path,
      num_workers       = n_cores))
  # Conform divide ownership to the specification, once the geometry is final.
  #
  # `divide_has_member_flowline` requires a divide to carry one of its flowpath's
  # flowlines. Everything upstream of here derives divides from a grid, and where
  # that grid and the vector network come from unrelated sources: Alaska pairs
  # TDX-Hydro with HydroSHEDS, the derivation has an accuracy floor set by the
  # inputs, not by the algorithm. This enforces the rule the gate states rather
  # than hoping the split achieved it, either by moving the divide to the
  # flowpath whose line it actually holds, or by absorbing it into a neighbour
  # and demoting its owner to a connector.
  #
  # It acts ONLY on divides that fail that test, so a domain already conforming
  # is untouched, which is why it cannot regress the 23 VPUs that pass today.
  # Placed after the last step that rewrites geometry and before Step 6d8, which
  # derives drainage area from these divides.
  .t("Step 6d4: normalize_divide_ownership",
    normalize_divide_ownership(ngen_gpkg))

  # No second prune pass here. A divide collapsed to zero area by overlap
  # resolution is a symptom of overlap that originates in stage 1, which the
  # splitter's reach-containment test addresses at source. Add one back only with
  # evidence that a collapsed divide survives a clean run.
  .t("Step 6d5: sync fp divide_id after cleanup", {
    fp      <- sf::st_drop_geometry(sf::st_read(ngen_gpkg, "flowpaths", quiet = TRUE))
    div_ids <- sf::st_read(ngen_gpkg, "divides",   quiet = TRUE)$divide_id
    bad     <- !is.na(fp$divide_id) & !(fp$divide_id %in% div_ids)
    if (any(bad)) {
      bad_ids <- as.character(fp$flowpath_id[bad])
      hfutils::gpkg_exec(ngen_gpkg,
        list(drop_triggers_for = "flowpaths"),
        list(table = "._bad", df = data.frame(._id = bad_ids,
          stringsAsFactors = FALSE)),
        'UPDATE "flowpaths" SET "divide_id" = NULL, "has_divide" = 0 WHERE CAST("flowpath_id" AS TEXT) IN (SELECT "._id" FROM "._bad")')
      message(sprintf("[stage3][%s] Step 6d5: cleared %d dangling fp divide_id(s)",
        vpuid, sum(bad)))
    }
  })

  .t("Step 6d5b: absorb multi-part divide slivers", {
    d0 <- sf::st_read(ngen_gpkg, "divides", quiet = TRUE)
    fl0 <- tryCatch(sf::st_read(ngen_gpkg, "flowlines", quiet = TRUE),
      error = function(e) NULL)
    r  <- .absorb_divide_slivers(d0, flowlines = fl0)
    if (r$n_absorbed > 0L) {
      sf::st_write(.drop_list_cols(r$div), ngen_gpkg, "divides",
        delete_layer = TRUE, quiet = TRUE)
      message(sprintf(paste0("[stage3][%s] Step 6d5b: absorbed %d sliver part(s) ",
        "(%.4f km2) into the neighbour sharing the most boundary; multi-part ",
        "divides %d -> %d%s"), vpuid, r$n_absorbed, r$area_moved_km2,
        r$n_multi_before, r$n_multi_after,
        if (r$n_kept) sprintf("; %d kept (no adjacent divide)", r$n_kept) else ""))
    } else if (r$n_multi_before > 0L) {
      message(sprintf(paste0("[stage3][%s] Step 6d5b: %d multi-part divide(s), ",
        "none eligible for absorption (parts are real catchment lobes)"),
        vpuid, r$n_multi_before))
    }
  })

  # Step 6d6 removed: it pruned sub-1m2 divide slivers. None exist on any of the
  # 24 VPUs now. The slivers were residue of the overlap-removal cascade; with
  # the seed fixed they are not produced. (It had already been converted from
  # deleting them to absorbing them, so this removes an absorber with no input.)

  # Step 6d6-sync removed: it cleared flowpath divide_id values left dangling by
  # 6c / 6d / 6d6 deleting divides. It never fired on any of the 24 VPUs, and
  # those three steps no longer delete, they absorb into a neighbour and clear
  # the owning flowpath themselves. Any residual dangling reference is still
  # caught downstream by Step 6j, which clears exactly that case, and asserted by
  # the divide_flowpath_exists invariant.

  .probe_cov("after-6d2-thru-6d6")
  # Step 6d6b (recover_missing_divides) removed: it recovered nothing on any of
  # the 24 VPUs, before or after the divide-conservation fixes, while costing
  # 132s of build time across the fabric. It existed to rebuild divides that
  # earlier steps had deleted, and the flowpaths it reports as "missing divides"
  # are connectors, which are divide-less by definition and not recoverable.
  # Steps that used to delete divides now absorb them into a neighbour, so there
  # is nothing left for it to put back.

  .t("Step 6d7: refresh lengthkm/areasqkm", {
    fp <- sf::st_read(ngen_gpkg, "flowpaths", quiet = TRUE)
    if ("lengthkm" %in% names(fp)) {
      fp$lengthkm <- as.numeric(sf::st_length(fp)) / 1000
    }
    fp <- .drop_list_cols(fp)
    sf::st_write(fp, ngen_gpkg, "flowpaths", delete_layer = TRUE, quiet = TRUE)

    d <- sf::st_read(ngen_gpkg, "divides", quiet = TRUE)
    # Guarantee divide validity once, here, the stage-3 mutation steps (6d
    # keep-fp-covering-part, 6d6b recover-via-rec_div) overlay/clip geometry and
    # can leave a self-intersection that crashes any downstream st_union (e.g.
    # the Step 6h basin mask: "side location conflict"). st_make_valid is
    # area-preserving for the fixable cases; keep polygonal parts only so the
    # layer stays (MULTI)POLYGON.
    inval <- !sf::st_is_valid(d)
    if (any(inval, na.rm = TRUE)) {
      n_inv <- sum(inval, na.rm = TRUE)
      d <- sf::st_make_valid(d)
      gc_rows <- as.character(sf::st_geometry_type(d)) == "GEOMETRYCOLLECTION"
      if (any(gc_rows)) {
        d[gc_rows, ] <- sf::st_collection_extract(d[gc_rows, ], "POLYGON", warn = FALSE)
      }
      message(sprintf("[stage3][%s] Step 6d7: repaired %d invalid divide geometry(ies)",
        vpuid, n_inv))
    }
    if ("areasqkm" %in% names(d)) {
      d$areasqkm <- as.numeric(sf::st_area(d)) / 1e6
    }
    d <- .drop_list_cols(d)
    sf::st_write(d, ngen_gpkg, "divides", delete_layer = TRUE, quiet = TRUE)
  })

  # Divide-based area model (FAC-free). Both flowline area fields are derived
  # from the reconciled-grain divides (rec_div) instead of the FAC grid, which
  # is unreliable at the sub-flowline grain (gaps / off-network cells leave many
  # segments near-zero).
  #   incremental: each reconciled line's rec_div catchment area, distributed
  #     across its flowline segments by length: reconciled/refactored segments
  #     (one per line) get the full area, length-split children subdivide it,
  #     connectors (no rec_div) get 0. Conserves basin area by construction.
  #   total_dasqkm: accumulated DOWN the (acyclic) flowline_id -> flowline_toid
  #     graph from the incremental, so total = incremental + sum(upstream totals);
  #     total >= incremental by construction and the outlet total = basin area.
  # Runs after 6d7 (lengthkm final) and before finalize (so the network inherits
  # both values).
  .t("Step 6d8: divide-based incremental + total_dasqkm (FAC-free)", {
    rec_div_gpkg <- file.path(work_dir, "rec_div.gpkg")
    if (!file.exists(rec_div_gpkg)) {
      message(sprintf("[stage3][%s] Step 6d8: rec_div.gpkg not found -- keeping FAC areas", vpuid))
    } else {
      fl <- sf::st_read(ngen_gpkg, "flowlines", quiet = TRUE)
      rd <- sf::st_read(rec_div_gpkg, quiet = TRUE)
      if (!all(c("reconciled_id", "flowline_id", "flowline_toid") %in% names(fl))) {
        message(sprintf("[stage3][%s] Step 6d8: flowlines lack id/toid -- keeping FAC areas", vpuid))
      } else {
        # --- incremental: length-weighted rec_div area ---
        rd_key  <- if ("divide_id" %in% names(rd)) "divide_id" else
          if ("flowpath_id" %in% names(rd)) "flowpath_id" else names(rd)[1]
        rd_area <- tapply(as.numeric(sf::st_area(rd)) / 1e6,
          as.character(rd[[rd_key]]), sum)
        rid <- as.character(fl$reconciled_id)
        len <- as.numeric(fl$lengthkm)
        len[is.na(len) | len <= 0] <- 0
        denom <- unname(tapply(len, rid, sum)[rid])      # total member length per line
        share <- ifelse(!is.na(denom) & denom > 0, len / denom, 0)
        inc   <- unname(rd_area[rid]) * share            # rec_div area * length share
        inc[is.na(inc)] <- 0                             # connectors / unmatched -> 0

        # --- total_dasqkm: accumulate incremental down the flowline DAG ---
        id <- as.character(fl$flowline_id)
        to <- as.character(fl$flowline_toid)
        n  <- length(id)
        idx    <- stats::setNames(seq_len(n), id)
        to_idx <- idx[to]
        valid  <- !is.na(to_idx) & to_idx != seq_len(n)  # edge to an in-set, non-self target
        tot <- inc
        if (any(valid)) {
          g <- igraph::graph_from_data_frame(
            data.frame(from = id[valid], to = to[valid], stringsAsFactors = FALSE),
            directed = TRUE, vertices = data.frame(name = id, stringsAsFactors = FALSE))
          ord <- as.integer(igraph::topo_sort(g, mode = "out"))  # upstream-first
          up_by <- vector("list", n)                              # recipient row -> donor rows
          sp <- split(which(valid), to_idx[valid])
          up_by[as.integer(names(sp))] <- sp
          for (v in ord) {
            d <- up_by[[v]]
            if (length(d)) tot[v] <- inc[v] + sum(tot[d])
          }
        }

        fl$flowline_incr <- inc
        fl$total_dasqkm  <- tot
        fl <- .drop_list_cols(fl)
        sf::st_write(fl, ngen_gpkg, "flowlines", delete_layer = TRUE, quiet = TRUE)
        message(sprintf(paste0(
          "[stage3][%s] Step 6d8: divide-based areas -- incremental sum=%.1f ",
          "km2, max total=%.1f km2 (outlet~basin), %d segs, %d incr>total"),
          vpuid, sum(inc), max(tot), n, sum(inc > tot + 1e-6)))
      }
    }
  })

  .probe_cov("after-6d6b-6d7")
  .t("Step 6e: build_aspatial_network",
    build_aspatial_network(ngen_gpkg,
      refactored     = refactored_sf,
      reconciled     = reconciled_sf,
      hydrolocations = hl_sf,
      verbose        = verbose))
  # Slope (length-weighted over member COMIDs) needs per-COMID lengthkm from the
  # reference. Prefer the per-VPU slice the batch presplit writes; otherwise fall
  # back to the reference_gpkg passed through from run_pipeline (e.g. a single-VPU
  # subset, where no ref_vpu<id> slice is ever created).
  ref_slice  <- file.path(fabric_dir, "work", paste0("ref_vpu", vpuid, ".gpkg"))
  ref_for_fdm <- if (file.exists(ref_slice)) ref_slice else reference_gpkg
  .t("Step 6f: finalize_data_model",
    finalize_data_model(ngen_gpkg,
      reference_gpkg = ref_for_fdm,
      vpuid = vpuid, events_gpkg = out_gpkg,
      # A channel table carrying slope serves as the slope source when no
      # separate slope table is given.
      slope_parquet = slope_parquet %||% hydraulics_path))

  # -- 6g. Routing attributes -----------------------------------------
  if (!is.null(hydraulics_path) && file.exists(hydraulics_path) &&
    !requireNamespace("arrow", quietly = TRUE))
    message(sprintf("[stage3][%s] Step 6g: hydraulics_path set but 'arrow' not installed; skipping routing attributes",
      vpuid))
  if (!is.null(hydraulics_path) && file.exists(hydraulics_path) &&
    requireNamespace("arrow", quietly = TRUE)) {
    .t("Step 6g: routing attributes (flowpath + flowline)", {
      hydraulics <- arrow::read_parquet(hydraulics_path, as_data_frame = FALSE)
      network    <- sf::st_read(ngen_gpkg, "network",   quiet = TRUE)
      flowpaths  <- sf::st_read(ngen_gpkg, "flowpaths", quiet = TRUE)

      # Per-reach weights for a table that carries no length of its own: the
      # network holds the length each reference reach contributes.
      reach_lengths <- dplyr::summarise(
        dplyr::group_by(dplyr::filter(sf::st_drop_geometry(network),
          !is.na(reference_id)), reference_id),
        lengthkm = sum(flowline_lengthkm, na.rm = TRUE), .groups = "drop")

      fp_atts <- .calculate_flow_atts(
        network   = network,
        hydraulics = hydraulics,
        hf_id      = "flowpath_id",
        flowpaths  = flowpaths,
        reach_lengths = reach_lengths)

      fl_atts <- .calculate_flow_atts(
        network   = network,
        hydraulics = hydraulics,
        hf_id      = "flowline_id",
        flowpaths  = flowpaths,
        reach_lengths = reach_lengths)

      sf::st_write(fp_atts, ngen_gpkg, "flowpath-attributes",
        delete_layer = TRUE, quiet = TRUE)
      sf::st_write(fl_atts, ngen_gpkg, "flowline-attributes",
        delete_layer = TRUE, quiet = TRUE)

      n_proxy_fp <- sum(fp_atts$attr_source == "geomorphic_proxy", na.rm = TRUE)
      n_proxy_fl <- sum(fl_atts$attr_source == "geomorphic_proxy", na.rm = TRUE)
      if (n_proxy_fp > 0 || n_proxy_fl > 0)
        message(sprintf(
          "[stage3][%s] routing proxy applied: %d flowpath(s), %d flowline(s)",
          vpuid, n_proxy_fp, n_proxy_fl))
    })
  } else if (!is.null(hydraulics_path)) {
    .hf_log_warn(sprintf("[stage3][%s] hydraulics_path not found, skipping Step 6g: %s",
      vpuid, hydraulics_path))
  }

  # -- 6h. Lakes ------------------------------------------------------
  if (!is.null(lakes_path) && file.exists(lakes_path)) {
    .t("Step 6h: lakes (outlet identification + network mask)", {
      s       <- sf::st_read(lakes_path, "lakes", quiet = TRUE)
      network <- sf::st_read(ngen_gpkg,  "network", quiet = TRUE)
      divs    <- sf::st_read(ngen_gpkg,  "divides", quiet = TRUE)

      s_local <- sf::st_transform(s, sf::st_crs(divs))

      # Primary outlet identification: use s$comid, the NWM-specified outlet
      # reach for each lake.  This avoids relying on NHDPlus wbareacomi which
      # can mis-index headwater reaches to distant lakes (known VAA error class).
      # Because the network is a long table where every flowpath on the same
      # levelpathid carries ALL of that group's member COMIDs, the outlet COMID
      # can appear in dozens of flowpaths.  min-hydroseq selects the one true
      # outlet flowpath (lowest hydroseq = most downstream).
      net_distinct <- network |>
        sf::st_drop_geometry() |>
        dplyr::select(flowline_id, flowpath_id, reference_id, flowline_hydroseq) |>
        dplyr::distinct()

      outlet_join <- dplyr::inner_join(
        net_distinct,
        dplyr::transmute(sf::st_drop_geometry(s_local),
          lake_id      = as.character(lake_id),
          outlet_comid = as.character(comid)),
        by = c("reference_id" = "outlet_comid")) |>
        dplyr::group_by(lake_id) |>
        dplyr::slice_min(flowline_hydroseq, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(-flowline_hydroseq)

      # Reach -> waterbody membership (versioned lake mask) used for the outlet
      # fallback below and for member stamping; replaces a live VAA fetch.
      lake_mask <- .read_lake_mask(lake_mask_path)

      # Fall back: if the NWM outlet COMID is not in this basin's network (e.g.
      # VPU-boundary lake whose outlet is in the adjacent VPU), use the lake
      # mask with min-hydroseq for any lakes that didn't match above.
      matched_lake_ids <- unique(outlet_join$lake_id)
      unmatched <- s_local[!as.character(s_local$lake_id) %in% matched_lake_ids, ]
      if (nrow(unmatched) > 0L) {
        comid_to_lake_fb <- lake_mask |>
          dplyr::rename(hf_id = comid, lake_id = wbareacomi) |>
          dplyr::filter(lake_id %in% as.numeric(unmatched$lake_id)) |>
          dplyr::mutate(hf_id = as.character(hf_id), lake_id = as.character(lake_id))
        fallback_join <- net_distinct |>
          dplyr::inner_join(comid_to_lake_fb, by = c("reference_id" = "hf_id")) |>
          dplyr::group_by(lake_id) |>
          dplyr::slice_min(flowline_hydroseq, with_ties = FALSE) |>
          dplyr::ungroup() |>
          dplyr::select(-flowline_hydroseq)
        outlet_join <- dplyr::bind_rows(outlet_join, fallback_join)
      }

      # Build lake-mask member lookup for lake_id stamping on network rows.
      all_matched_lake_ids <- as.integer(unique(outlet_join$lake_id))
      comid_members <- lake_mask |>
        dplyr::rename(hf_id = comid, lake_id = wbareacomi) |>
        dplyr::filter(lake_id %in% all_matched_lake_ids) |>
        dplyr::mutate(hf_id = as.character(hf_id), lake_id = as.character(lake_id))

      lakes <- outlet_join |>
        dplyr::inner_join(
          sf::st_drop_geometry(s_local) |>
            dplyr::mutate(lake_id = as.character(lake_id)),
          by = "lake_id") |>
        sf::st_sf(geometry = sf::st_geometry(
          s_local[match(outlet_join$lake_id, as.character(s_local$lake_id)), ]))

      lakes <- sf::st_as_sf(lakes)

      # Drop lakes that don't intersect the basin footprint: these are VAA index
      # errors (a COMID mapped to a lake that is spatially outside the watershed).
      basin_mask <- suppressWarnings(sf::st_union(
        sf::st_make_valid(sf::st_geometry(divs))))
      in_basin <- lengths(sf::st_intersects(lakes, sf::st_sf(geometry = basin_mask))) > 0L
      if (any(!in_basin, na.rm = TRUE)) {
        bad_ids <- as.character(unique(lakes$lake_id[!in_basin]))
        .hf_log_warn(sprintf(
          "[stage3][%s] lakes: dropped %d lake(s) that do not intersect basin footprint (bad VAA index): %s",
          vpuid, length(bad_ids), paste(utils::head(bad_ids, 10), collapse = ", ")))
        lakes <- lakes[in_basin, , drop = FALSE]
      }

      sf::st_write(lakes, ngen_gpkg, "lakes", delete_layer = TRUE, quiet = TRUE)

      # -- Waterbody crosswalk: stamp lake_id on EVERY in-lake reach ---------
      # t-route replaces a waterbody's underlying reaches with reservoir routing,
      # so it needs lake_id on all member flowpaths/flowlines, not just the outlet.
      # `member_comid` on the flowpaths layer is a precise per-flowpath COMID list
      # (unlike the network's levelpath-broad reference_id), so we assign each
      # flowpath to a lake from its own members and propagate. A flowpath is in a
      # lake when the MAJORITY of its member-COMID length belongs to that lake's
      # VAA members (member-count majority when per-COMID lengths are unavailable).
      final_lake_ids <- as.character(unique(lakes$lake_id))
      cm         <- comid_members[comid_members$lake_id %in% final_lake_ids, , drop = FALSE]
      comid2lake <- stats::setNames(cm$lake_id, cm$hf_id)

      # Per-COMID length for the majority rule (COMID == reference flowpath_id).
      ref_for_lakes <- if (file.exists(ref_slice)) ref_slice else reference_gpkg
      comid_len <- if (!is.null(ref_for_lakes) && file.exists(ref_for_lakes)) tryCatch(
        {
          cc <- DBI::dbConnect(RSQLite::SQLite(), ref_for_lakes)
          q  <- DBI::dbGetQuery(cc, sprintf(
            "SELECT flowpath_id, lengthkm FROM flowpaths WHERE vpuid = '%s'", vpuid))
          DBI::dbDisconnect(cc)
          stats::setNames(as.numeric(q$lengthkm), as.character(q$flowpath_id))
        },
        error = function(e) NULL) else NULL

      .assign_lake <- function(mc_str) {
        if (is.na(mc_str) || !nzchar(mc_str)) return(NA_character_)
        cids <- sub("[.][0-9]+$", "", trimws(strsplit(mc_str, ",", fixed = TRUE)[[1]]))
        lk   <- unname(comid2lake[cids])
        if (all(is.na(lk))) return(NA_character_)
        w <- if (!is.null(comid_len)) {
          ww <- unname(comid_len[cids])
          ww[is.na(ww) | ww <= 0] <- 0
          ww
        } else rep(1, length(cids))
        if (all(w == 0)) w <- rep(1, length(cids))          # no usable lengths -> count
        inlake <- !is.na(lk)
        if (sum(w[inlake]) <= sum(w[!inlake])) return(NA_character_)   # majority must be in-lake
        agg <- tapply(w[inlake], lk[inlake], sum)           # dominant lake among members
        names(agg)[which.max(agg)]
      }

      fp_for_lake <- sf::st_read(ngen_gpkg, "flowpaths", quiet = TRUE)
      fp_lake_vec <- vapply(as.character(fp_for_lake$member_comid), .assign_lake,
        character(1L), USE.NAMES = FALSE)

      # Guarantee each lake's outlet flowpath carries its lake_id even if the
      # majority rule didn't (e.g. a lake that is only a minor part of a large
      # flowpath): t-route needs at least the outlet reach of every waterbody.
      # The outlet flowpath is identified upstream via NWM comid (reliable).
      out_lake   <- stats::setNames(as.character(lakes$lake_id), as.character(lakes$flowpath_id))
      fp_ids_chr <- as.character(fp_for_lake$flowpath_id)
      need_out   <- is.na(fp_lake_vec) & fp_ids_chr %in% names(out_lake)
      fp_lake_vec[need_out] <- unname(out_lake[fp_ids_chr[need_out]])

      # Membership is the VAA `comid -> wbareacomi` mapping, and that is all of it.
      #
      # There is no distance bound on a stamp. The lakes are POINTS, so there is
      # no extent to scale one against, and a flat bound deletes reaches that are
      # legitimately upstream of their own outlet. The test that replaces it is
      # applied below, once the outlet is corrected: a member must be reachable
      # from its lake's outlet.

      # Close lake-id gaps, so a lake's flowpaths form a contiguous subgraph.
      #
      # Routing collapses every flowpath sharing a lake_id into one waterbody
      # node, which is only valid if they are connected: a non-member between two
      # members makes node_L -> gap -> node_L, a cycle.
      #
      # The majority rule above decides each flowpath alone, so a reach mostly
      # outside the lake is rejected even when both its neighbours are members. A
      # bridging reach belongs to the waterbody because both its ends do. Bounded
      # to a short run, so a genuinely long river between two arms of one
      # waterbody is reported rather than absorbed.
      if (any(!is.na(fp_lake_vec))) {
        nx_lk <- tryCatch(sf::st_read(ngen_gpkg, "nexus", quiet = TRUE),
                          error = function(e) NULL)
        if (!is.null(nx_lk) && nrow(nx_lk) > 0L) {
          .n2f  <- stats::setNames(as.character(nx_lk$nexus_toid),
                                   as.character(nx_lk$nexus_id))
          .fid  <- as.character(fp_for_lake$flowpath_id)
          .nxt  <- match(unname(.n2f[as.character(fp_for_lake$flowpath_toid)]), .fid)
          .len  <- suppressWarnings(as.numeric(fp_for_lake$lengthkm))
          .MAXSTEP <- 5L
          .n_fill <- 0L; .km_fill <- 0
          for (.i in which(!is.na(fp_lake_vec))) {
            .L <- fp_lake_vec[.i]; .j <- .nxt[.i]; .path <- integer(0)
            for (.s in seq_len(.MAXSTEP)) {
              if (is.na(.j)) break
              if (!is.na(fp_lake_vec[.j])) {
                if (identical(fp_lake_vec[.j], .L) && length(.path)) {
                  fp_lake_vec[.path] <- .L
                  .n_fill  <- .n_fill + length(.path)
                  .km_fill <- .km_fill + sum(.len[.path], na.rm = TRUE)
                }
                break
              }
              .path <- c(.path, .j); .j <- .nxt[.j]
            }
          }
          if (.n_fill > 0L)
            message(sprintf(paste0(
              "[stage3][%s] lakes: closed %d lake_id gap(s) totalling %.3f km -- ",
              "flowpaths bridging two members of one lake, which t-route would ",
              "otherwise collapse into a cycle"), vpuid, .n_fill, .km_fill))
        }
      }

      # Carry the outlet down to the lake's most downstream member.
      #
      # The outlet is chosen from the single flowpath containing the RouteLink
      # `comid`, so it cannot see that the waterbody continues past it. On VPU 01,
      # lake 120052046's outlet resolved to fp-32570 while fp-32552 and fp-32550
      # sit 1 and 2 hops below it with every one of their 8 and 4 member reaches
      # inside the lake. The lake really does extend downstream; the outlet was
      # simply picked from the wrong flowpath.
      #
      # Walk down only while the next flowpath is still a member of the SAME lake,
      # so the outlet moves only where the membership demonstrates that it should.
      if (any(!is.na(fp_lake_vec)) && nrow(lakes) > 0L && exists(".nxt")) {
        .fidv  <- as.character(fp_for_lake$flowpath_id)
        .moved <- 0L
        .lk_out <- as.character(lakes$flowpath_id)
        .lk_id  <- as.character(lakes$lake_id)
        for (.q in seq_along(.lk_out)) {
          .cur <- match(.lk_out[.q], .fidv)
          if (is.na(.cur)) next
          repeat {
            .d <- .nxt[.cur]
            if (is.na(.d)) break
            if (is.na(fp_lake_vec[.d]) || !identical(fp_lake_vec[.d], .lk_id[.q])) break
            .cur <- .d
          }
          if (!identical(.fidv[.cur], .lk_out[.q])) {
            .lk_out[.q] <- .fidv[.cur]; .moved <- .moved + 1L
          }
        }
        if (.moved > 0L) {
          lakes$flowpath_id <- .lk_out
          # The lakes layer was written further up, before membership was known,
          # so the corrected outlet has to be persisted or it exists only here:
          # the drop below would use it while every downstream reader, the
          # invariants included: still saw the reach it replaced.
          try(sf::st_write(.drop_list_cols(lakes), ngen_gpkg, "lakes",
            delete_layer = TRUE, quiet = TRUE), silent = TRUE)
          message(sprintf(paste0(
            "[stage3][%s] lakes: moved %d outlet(s) downstream to the lake's ",
            "most downstream member (RouteLink comid names a reach the waterbody ",
            "continues past)"), vpuid, .moved))
        }
      }

      # Drop members the network cannot reach from their lake's outlet.
      #
      # `wbareacomi` is NHDPlus, and its known mis-index class puts a reach in a
      # waterbody it has no connection to. t-route collapses a lake's members
      # into one node, so a member on an unrelated branch invents a routing
      # connection, the same corruption as a lake-id gap, in a different shape,
      # and one the contiguity check cannot see because the member is not BETWEEN
      # anything.
      #
      # Judged after the outlet correction above, so a member that merely sits
      # below a mis-picked outlet is not mistaken for a mis-index.
      #
      # Reachability, not distance: the previous bound measured metres to a lake
      # POINT and was 99% false positives. Dropped rows are kept in
      # `lake_id_flags` and are worth reporting upstream, this is a workaround
      # for a base-data error, not a repair of it.
      if (any(!is.na(fp_lake_vec)) && nrow(lakes) > 0L &&
          all(c("upstream_id", "num_upstreams") %in% names(fp_for_lake))) {
        .u   <- suppressWarnings(as.numeric(fp_for_lake$upstream_id))
        .nu  <- suppressWarnings(as.numeric(fp_for_lake$num_upstreams))
        .fpc <- as.character(fp_for_lake$flowpath_id)
        .out <- stats::setNames(as.character(lakes$flowpath_id),
                                as.character(lakes$lake_id))
        .k   <- which(!is.na(fp_lake_vec))
        .oi  <- match(unname(.out[fp_lake_vec[.k]]), .fpc)
        .ok  <- !is.na(.oi) & is.finite(.u[.k]) & is.finite(.u[.oi]) &
                .u[.k] >= .u[.oi] & .u[.k] <= .u[.oi] + .nu[.oi]
        .bad <- .k[!is.na(.oi) & !.ok]
        if (length(.bad) > 0L) {
          flagged <- data.frame(
            flowpath_id        = .fpc[.bad],
            lake_id            = fp_lake_vec[.bad],
            outlet_flowpath_id = unname(.out[fp_lake_vec[.bad]]),
            lengthkm           = suppressWarnings(as.numeric(fp_for_lake$lengthkm))[.bad],
            stringsAsFactors   = FALSE)
          .hf_log_warn(sprintf(paste0(
            "[stage3][%s] lakes: dropped %d of %d member(s) unreachable from ",
            "their lake's outlet (NHDPlus wbareacomi mis-index) -- FLAGGED to ",
            "`lake_id_flags`"), vpuid, length(.bad), length(.k)))
          try(sf::st_write(flagged, ngen_gpkg, "lake_id_flags",
            delete_layer = TRUE, quiet = TRUE), silent = TRUE)
          fp_lake_vec[.bad] <- NA_character_
        }
      }

      fp_for_lake$lake_id <- fp_lake_vec
      fp_for_lake <- .drop_list_cols(fp_for_lake)
      sf::st_write(fp_for_lake, ngen_gpkg, "flowpaths", delete_layer = TRUE, quiet = TRUE)

      fp2lake <- stats::setNames(fp_lake_vec, as.character(fp_for_lake$flowpath_id))

      # Propagate to flowlines (inherit parent flowpath's lake).
      if (hfutils::layer_exists(ngen_gpkg, "flowlines")) {
        fl_for_lake <- sf::st_read(ngen_gpkg, "flowlines", quiet = TRUE)
        fl_for_lake$lake_id <- unname(fp2lake[as.character(fl_for_lake$flowpath_id)])
        fl_for_lake <- .drop_list_cols(fl_for_lake)
        sf::st_write(fl_for_lake, ngen_gpkg, "flowlines", delete_layer = TRUE, quiet = TRUE)
      }

      # Stamp the network (every row of an in-lake flowpath): keyed on the
      # flowpath's own lake, so no over-stamping from shared reference_id.
      net_df <- sf::st_drop_geometry(network)
      net_df$lake_id <- unname(fp2lake[as.character(net_df$flowpath_id)])
      sf::st_write(net_df, ngen_gpkg, "network", delete_layer = TRUE, quiet = TRUE)

      message(sprintf("[stage3][%s] lakes written: %d | flowpaths in lakes: %d",
        vpuid, nrow(lakes), sum(!is.na(fp_lake_vec))))
    })
  } else if (!is.null(lakes_path)) {
    .hf_log_warn(sprintf("[stage3][%s] lakes_path not found, skipping Step 6h: %s",
      vpuid, lakes_path))
  }

  # -- 6i. Global identification validation ------------------------------
  .probe_cov("after-6e-6f-6g-6h")
  .t("Step 6i: validate_global_ids", {
    fp  <- sf::st_read(ngen_gpkg, "flowpaths",      quiet = TRUE)
    div <- sf::st_read(ngen_gpkg, "divides",         quiet = TRUE)
    nex <- sf::st_read(ngen_gpkg, "nexus",           quiet = TRUE)
    net <- sf::st_read(ngen_gpkg, "network",         quiet = TRUE)
    hl  <- if (hfutils::layer_exists(ngen_gpkg, "hydrolocations"))
      sf::st_read(ngen_gpkg, "hydrolocations", quiet = TRUE) else NULL

    .chk <- function(layer_name, x, col, prefix = NULL, allow_na = FALSE) {
      if (!col %in% names(x)) {
        .hf_log_warn(sprintf("[ids][%s] %s.%s MISSING", vpuid, layer_name, col))
        return(invisible(NULL))
      }
      vals <- x[[col]]
      na_n <- sum(is.na(vals))
      if (!allow_na && na_n > 0)
        .hf_log_warn(sprintf("[ids][%s] %s.%s: %d NA value(s)", vpuid, layer_name, col, na_n))
      if (!is.null(prefix)) {
        bad <- !is.na(vals) & !grepl(paste0("^", prefix), as.character(vals))
        if (any(bad))
          .hf_log_warn(sprintf("[ids][%s] %s.%s: %d value(s) missing '%s' prefix",
            vpuid, layer_name, col, sum(bad), prefix))
      }
    }

    # Flowpaths
    .chk("flowpaths", fp, "flowpath_id",   prefix = "fp-")
    .chk("flowpaths", fp, "flowpath_toid")
    .chk("flowpaths", fp, "divide_id",     prefix = "cat-", allow_na = TRUE)
    .chk("flowpaths", fp, "mainstem_id",   allow_na = TRUE)
    .chk("flowpaths", fp, "vpuid")
    .chk("flowpaths", fp, "hydroseq")
    .chk("flowpaths", fp, "streamorder",   allow_na = TRUE)
    .chk("flowpaths", fp, "has_divide")
    .chk("flowpaths", fp, "member_comid")

    # Divides
    .chk("divides", div, "divide_id")
    .chk("divides", div, "flowpath_toid",  allow_na = TRUE)
    .chk("divides", div, "vpuid")
    .chk("divides", div, "areasqkm")
    .chk("divides", div, "type")

    # Nexus
    .chk("nexus", nex, "nexus_id")
    .chk("nexus", nex, "nexus_toid",  allow_na = TRUE)
    .chk("nexus", nex, "vpuid")

    # Network: key crosswalk columns
    .chk("network", net, "flowpath_id",   allow_na = TRUE)
    .chk("network", net, "flowline_id")
    .chk("network", net, "reference_id",  allow_na = TRUE)
    .chk("network", net, "hf_source",     allow_na = TRUE)
    .chk("network", net, "hf_toid",       allow_na = TRUE)
    .chk("network", net, "vpuid")
    .chk("network", net, "divide_id",     allow_na = TRUE)
    .chk("network", net, "mainstem_id",   allow_na = TRUE)

    # Hydrolocations
    if (!is.null(hl)) {
      .chk("hydrolocations", hl, "poi_id")
      .chk("hydrolocations", hl, "flowpath_id", prefix = "fp-")
      .chk("hydrolocations", hl, "hl_uri",      allow_na = TRUE)
      .chk("hydrolocations", hl, "hl_type",     allow_na = TRUE)
      .chk("hydrolocations", hl, "hl_reference", allow_na = TRUE)
      .chk("hydrolocations", hl, "vpuid")
    }

    # Referential integrity: fp divide_id resolves to a divide
    fp_div <- stats::na.omit(fp$divide_id)
    div_ids <- div$divide_id
    unresolved_div <- sum(!fp_div %in% div_ids)
    if (unresolved_div > 0)
      .hf_log_warn(sprintf("[ids][%s] %d flowpath divide_id(s) do not resolve to a divide",
        vpuid, unresolved_div))

    # Referential integrity: network flowpath_id resolves to fp
    net_fp <- stats::na.omit(net$flowpath_id)
    unresolved_fp <- sum(!net_fp %in% fp$flowpath_id)
    if (unresolved_fp > 0)
      .hf_log_warn(sprintf("[ids][%s] %d network flowpath_id(s) do not resolve to a flowpath",
        vpuid, unresolved_fp))

    # Geometry types
    fp_gt  <- as.character(sf::st_geometry_type(fp))
    div_gt <- as.character(sf::st_geometry_type(div))
    bad_fp  <- fp_gt  != "LINESTRING"
    bad_div <- div_gt != "POLYGON"
    if (any(bad_fp))
      .hf_log_warn(sprintf("[ids][%s] %d flowpath(s) not LINESTRING: %s",
        vpuid, sum(bad_fp), paste(unique(fp_gt[bad_fp]), collapse = ", ")))
    if (any(bad_div))
      .hf_log_warn(sprintf("[ids][%s] %d divide(s) not POLYGON: %s",
        vpuid, sum(bad_div), paste(unique(div_gt[bad_div]), collapse = ", ")))

    message(sprintf(
      "[ids][%s] fp=%d div=%d nex=%d net=%d hl=%d | connectors=%d (%.1f%%)",
      vpuid, nrow(fp), nrow(div), nrow(nex), nrow(net),
      if (is.null(hl)) 0L else nrow(hl),
      sum(!fp$has_divide, na.rm = TRUE),
      100 * mean(!fp$has_divide, na.rm = TRUE)))
  })

  # Steps 6j and 6k2 removed. 6j cleared flowpaths whose divide_id no longer
  # overlapped them; once narrowed to genuine dangling references only, it cleared
  # none on any of the 24 VPUs, no dangling refs, no low-coverage flowpaths, no
  # connector areas to zero. 6k2 re-homed the divides 6j and 6k orphaned, and
  # re-homed nothing. Both were compensating for the stage-2 overlap seed; with
  # that fixed neither has input. A dangling reference is now asserted by the
  # divide_flowpath_exists invariant rather than silently repaired.
  .probe_cov("after-6i")

  # -- 6k. Demote divide-less headwater-leaf flowpaths to flowline-only -----
  # Reaches whose reference flowlines were divide-less in the source fabric
  # surface as catchment-less flowpaths (has_divide=FALSE, areasqkm=0). The
  # connector set is only finalized at Step 6j (zero-coverage divide clearing),
  # so this runs immediately after it. Where such a reach is a headwater leaf,
  # lift it out of the flowpath grain: keep its flowlines (re-pointed to the
  # catchment that contains them), drop its flowpath/attributes/nexus rows.
  # Iterate: demoting a leaf can expose a new headwater leaf upstream. Each
  # pass peels one layer off a chain, so the passes needed scale with the
  # longest chain and the bound is generous. Reaching it warns rather than
  # exiting quietly, since survivors would trip
  # no_divideless_leaf_flowpaths.
  .t("Step 6k: demote divide-less leaf flowpaths", {
    .converged <- FALSE
    for (.k in 1:50) {
      n_dem <- .demote_divideless_leaf_flowpaths(ngen_gpkg, vpuid = vpuid,
        lakes_path = lakes_path)
      if (!isTRUE(n_dem > 0L)) { .converged <- TRUE; break }
    }
    if (!.converged)
      .hf_log_warn(sprintf(paste0("[stage3][%s] Step 6k did not converge in 50 passes; ",
        "divide-less leaf flowpath(s) may survive"), vpuid))
  })

  # Step 6k2 removed: it re-homed the divides that 6j (clearing a divide_id) and
  # 6k (deleting a demoted flowpath) left owned by a flowpath that no longer
  # claimed them. It re-homed nothing on any of the 24 VPUs. It was added earlier
  # in this same effort as a safety net for those two steps; fixing the stage-2
  # overlap seed removed the orphans it existed to repair. The condition it
  # repaired is asserted by the divide_flowpath_exists and
  # divide_flowpath_consistent gates.

  # Step 6l removed: it re-oriented aggregated flowpath geometry to outlet-last,
  # because geos_line_merge does not preserve orientation. It reversed 0 of every
  # VPU's flowpaths on all 24, the geometry already arrives outlet-last, so there
  # was nothing to re-orient. The condition is asserted by the
  # flowpath_flow_direction gate, which is green on all 24.
  #
  # One caveat, recorded rather than glossed: that gate and the removed step derive
  # each outlet from the same flowline signal, so they are not fully independent.
  # Removing a step that reversed nothing does not widen that blind spot, but the
  # gate should not be read as proof the signal itself is right. A truly
  # independent orientation check would come from hydroseq or elevation.
  #
  # enforce_flowpath_flowdir stays exported. It is a public entry point with no
  # remaining internal caller, usable for manual repair; removing it would be a
  # breaking API change for no benefit.

  # -- 6m. Flag flowpath constituents ---------------------------------------
  # Mark each flowline TRUE/FALSE for whether it is directly part of its
  # flowpath's traced line (a stem constituent) vs a tributary that drains into
  # the flowpath. Runs LAST, after every flowpath-geometry (incl. 6l) and
  # flowline re-pointing step, so the flag matches the final line + membership.
  if (hfutils::layer_exists(ngen_gpkg, "flowlines"))
    .t("Step 6m: flag flowpath constituents",
      tryCatch(flag_flowpath_constituents(ngen_gpkg),
        error = function(e) {
          .hf_log_warn(sprintf("[stage3][%s] Step 6m non-fatal: %s", vpuid,
            conditionMessage(e)))
          invisible(NULL)
        }))

  # -- 6n. Flag the routeable network ---------------------------------------
  # Mark the flowlines OWP routes on: the coincident (flowpath_constituent)
  # reaches plus the incremental reaches that bridge each coincident reach to
  # the next coincident reach downstream, so the routeable subgraph is closed
  # under flowline_toid. Runs after 6m (needs flowpath_constituent).
  if (hfutils::layer_exists(ngen_gpkg, "flowlines"))
    .t("Step 6n: flag routeable flowlines",
      tryCatch(flag_routeable_flowlines(ngen_gpkg),
        error = function(e) {
          .hf_log_warn(sprintf("[stage3][%s] Step 6n non-fatal: %s", vpuid,
            conditionMessage(e)))
          invisible(NULL)
        }))

  # Routing channel-slope `So` lives in the aspatial *-attributes tables, not the
  # geometry layers; join it on so the So_valid invariant actually validates it.
  .join_So <- function(sfx, attr_tbl, key) {
    if (!is.null(sfx) && hfutils::layer_exists(ngen_gpkg, attr_tbl)) {
      at <- tryCatch(sf::st_read(ngen_gpkg, attr_tbl, quiet = TRUE),
        error = function(e) NULL)
      if (!is.null(at) && all(c("So", key) %in% names(at)) && key %in% names(sfx))
        sfx$So <- at$So[match(as.character(sfx[[key]]), as.character(at[[key]]))]
    }
    sfx
  }
  # -- 6o. Re-apply id typing + join indexes --------------------------------
  # Steps 6g/6m/6n rewrite flowpaths and flowlines with delete_layer = TRUE,
  # which recreates the table and drops the INTEGER types and join indexes
  # finalize_data_model established, leaving R's doubles as REAL. Identifiers
  # stored as doubles are how ids get lost to overflow and scientific notation.
  # Runs after the last layer rewrite, so it is the state that ships.
  # Re-close lake gaps now that the topology is final. Step 6h closes them, but
  # 6k deletes divide-less leaf flowpaths and rewires `flowpath_toid` around
  # them, which can place a non-member between two members that were not adjacent
  # before. Runs ahead of id enforcement because it rewrites
  # flowpaths/flowlines/network.
  .t("Step 6n2: close lake id gaps (final topology)", {
    .lk_fix <- close_lake_id_gaps(ngen_gpkg)
    if (length(.lk_fix))
      message(sprintf(paste0(
        "[stage3][%s] lakes: closed %d further lake_id gap(s) after the graph ",
        "settled -- flowpaths bridging two members of one lake (e.g. %s)"),
        vpuid, length(.lk_fix), paste(utils::head(.lk_fix, 3), collapse = ", ")))
  })

  # Step 6d5b runs before the steps that reshape divides: 6k demotes
  # divide-less leaf flowpaths, and finalize rewrites the layer, so slivers
  # created after it survived to the product. Absorbing again once nothing else
  # will touch the geometry catches those, and a degenerate part (under a square
  # metre, an artefact of union arithmetic rather than land) is dropped outright
  # because there is no area to give anyone.
  .t("Step 6n3: final multi-part divide cleanup", {
    d <- sf::st_read(ngen_gpkg, "divides", quiet = TRUE)
    a0 <- sum(suppressWarnings(as.numeric(sf::st_area(d))), na.rm = TRUE)
    g  <- sf::st_geometry(d)
    multi <- which(lengths(g) > 1L)
    n_deg <- 0L
    for (i in multi) {
      parts <- g[[i]]
      ar <- vapply(parts, function(pp)
        as.numeric(sf::st_area(sf::st_sfc(sf::st_polygon(pp), crs = sf::st_crs(d)))),
        numeric(1))
      keep <- ar >= 1                      # square metres
      if (all(keep) || !any(keep)) next
      n_deg <- n_deg + sum(!keep)
      g[[i]] <- sf::st_multipolygon(parts[keep])
    }
    if (n_deg) {
      sf::st_geometry(d) <- sf::st_cast(sf::st_sfc(g, crs = sf::st_crs(d)), "MULTIPOLYGON")
      message(sprintf("[stage3][%s] Step 6n3: dropped %d degenerate part(s) (<1 m2)",
        vpuid, n_deg))
    }
    fl <- tryCatch(sf::st_read(ngen_gpkg, "flowlines", quiet = TRUE),
      error = function(e) NULL)
    r <- .absorb_divide_slivers(d, flowlines = fl)
    d <- r$div
    if (r$n_absorbed || n_deg) {
      sf::st_write(.drop_list_cols(d), ngen_gpkg, "divides",
        delete_layer = TRUE, quiet = TRUE)
      a1 <- sum(suppressWarnings(as.numeric(sf::st_area(d))), na.rm = TRUE)
      message(sprintf(paste0("[stage3][%s] Step 6n3: absorbed %d late sliver(s) ",
        "(%.6f km2); multi-part divides %d -> %d; area delta %.6f km2"),
        vpuid, r$n_absorbed, r$area_moved_km2, r$n_multi_before, r$n_multi_after,
        (a1 - a0) / 1e6))
    }
  })

  .t("Step 6o: enforce id types", enforce_integer_ids(ngen_gpkg))

  fp_chk <- .join_So(sf::st_read(ngen_gpkg, "flowpaths", quiet = TRUE),
    "flowpath-attributes", "flowpath_id")
  fl_chk <- if (hfutils::layer_exists(ngen_gpkg, "flowlines"))
    .join_So(sf::st_read(ngen_gpkg, "flowlines", quiet = TRUE),
      "flowline-attributes", "flowline_id") else NULL

  invariants <- tryCatch(
    check_invariants("ngen",
      gpkg       = ngen_gpkg,
      reference_gpkg = ref_for_fdm,
      flowpaths  = fp_chk,
      divides    = sf::st_read(ngen_gpkg, "divides",    quiet = TRUE),
      nexus      = sf::st_read(ngen_gpkg, "nexus",      quiet = TRUE),
      network    = if (hfutils::layer_exists(ngen_gpkg, "network"))
        sf::st_read(ngen_gpkg, "network",   quiet = TRUE) else NULL,
      flowlines  = fl_chk,
      lakes      = if (hfutils::layer_exists(ngen_gpkg, "lakes"))
        sf::st_read(ngen_gpkg, "lakes", quiet = TRUE) else NULL,
      # Stage-1 rec_div gates divides_no_enclosed_gaps: distinguishes real
      # dropped-catchment gaps from open-water/outside-domain voids (Great Lakes,
      # closed basins) that were never catchments.
      rec_div    = local({
        p <- file.path(work_dir, "rec_div.gpkg")
        if (file.exists(p)) tryCatch(sf::st_read(p, "divides", quiet = TRUE),
          error = function(e) NULL) else NULL
      }),
      # POI contract inputs: the delivered POI and hydrolocation layers, plus the
      # stage-1 event set they must all still be present in.
      pois       = if (hfutils::layer_exists(ngen_gpkg, "pois"))
        tryCatch(sf::st_read(ngen_gpkg, "pois", quiet = TRUE),
          error = function(e) NULL) else NULL,
      hydrolocations = if (hfutils::layer_exists(ngen_gpkg, "hydrolocations"))
        tryCatch(sf::st_read(ngen_gpkg, "hydrolocations", quiet = TRUE),
          error = function(e) NULL) else NULL,
      pois_input = local({
        p <- file.path(work_dir, "ckpt_events.gpkg")
        if (file.exists(p)) tryCatch(sf::st_read(p, "events", quiet = TRUE),
          error = function(e) NULL) else NULL
      }),
      refactored = refactored_sf,
      coverage  = TRUE, coverage_min = params$coverage_min, strict = FALSE,
      area = TRUE, area_baseline_m2 = baseline_div_area_m2,
      coverage_baseline_m2 = baseline_div_cov_m2),
    error = function(e) list(error = conditionMessage(e)))

  # Copy the completed ngen_gpkg from local temp to the NTFS destination.
  if (!file.copy(ngen_gpkg, ngen_gpkg_final, overwrite = TRUE))
    .hf_log_warn(sprintf("[stage3][%s] Failed to copy ngen_gpkg to %s", vpuid, ngen_gpkg_final))
  else
    unlink(ngen_gpkg)
  ngen_gpkg <- ngen_gpkg_final

  files <- list(
    ngen_gpkg  = ngen_gpkg,
    timings    = file.path(work_dir, "timings.json"),
    invariants = file.path(work_dir, "invariants.json"))

  jsonlite::write_json(timings,    files$timings,    auto_unbox = TRUE,
    pretty = TRUE)
  jsonlite::write_json(invariants, files$invariants, auto_unbox = TRUE,
    pretty = TRUE, force = TRUE)

  invisible(list(files = unlist(files), timings = timings,
    invariants = invariants))
}


# ---------------------------------------------------------------------------
# Demote divide-less headwater-leaf flowpaths to flowline-only
# ---------------------------------------------------------------------------
# A flowpath with no divide (`has_divide = FALSE`, `areasqkm = 0`) carries no
# incremental drainage; its reference flowlines bounded no catchment in the
# source fabric. Where such a reach is a headwater leaf, nothing drains into
# it and it can leave the flowpath grain without affecting connectivity:
#   - its flowlines stay in `flowlines`/`network`, re-pointed to the flowpath
#     whose divide contains them, which is the one its nexus drains to
#   - its `flowpaths` row is dropped
#   - a nexus left with no inflows is dropped, and its downstream flowpath
#     becomes a headwater
#
# Guarded: only LEAF (no upstream inflow), non-POI, non-lake divide-less
# flowpaths are demoted. Interior divide-less reaches (with real upstream
# flowpaths) are left as connectors; splicing their upstream across the removed
# node is a separate, riskier operation. The fold is area-conserving: every
# demoted stub has incremental area 0, so no divide area moves and no downstream
# total drainage area changes.
#
# Runs as Step 6k: AFTER the connector set is finalized (Step 6j clears
# zero-coverage divides) and after finalize/lakes, so `has_divide`, `lake_id`,
# the *-attributes layers, and the fully-denormalized `network` are all present.
# The lake guard uses the stamped `lake_id` column (with a `lakes_path`
# geometric fallback).
.demote_divideless_leaf_flowpaths <- function(gpkg, vpuid = NA_character_,
                                              lakes_path = NULL) {
  if (!hfutils::layer_exists(gpkg, "flowpaths") ||
    !hfutils::layer_exists(gpkg, "nexus")) return(invisible(NULL))

  fp_sf <- sf::st_read(gpkg, "flowpaths", quiet = TRUE)
  fpd   <- sf::st_drop_geometry(fp_sf)
  fpd$flowpath_id   <- as.character(fpd$flowpath_id)
  fpd$flowpath_toid <- as.character(fpd$flowpath_toid)
  nex_sf <- sf::st_read(gpkg, "nexus", quiet = TRUE)
  nexd   <- sf::st_drop_geometry(nex_sf)
  nexd$nexus_id   <- as.character(nexd$nexus_id)
  nexd$nexus_toid <- as.character(nexd$nexus_toid)

  # candidate = divide-less, non-POI, non-lake
  has_div <- if ("has_divide" %in% names(fpd)) {
    v <- suppressWarnings(as.logical(fpd$has_divide)); v[is.na(v)] <- FALSE; v
  } else !is.na(fpd$divide_id) & nzchar(as.character(fpd$divide_id))
  is_poi  <- "poi_id"  %in% names(fpd) &
    !is.na(fpd$poi_id)  & nzchar(as.character(fpd$poi_id))
  is_lake <- "lake_id" %in% names(fpd) &
    !is.na(fpd$lake_id) & nzchar(as.character(fpd$lake_id))
  cand <- !has_div & !is_poi & !is_lake

  # leaf = no flowpath drains into this fp's inflow nexus. fp X's inflow nexus is
  # the nexus whose nexus_toid == X; a nexus's donors are the flowpaths whose
  # flowpath_toid == that nexus_id.
  inflow_nex_of  <- stats::setNames(nexd$nexus_id, nexd$nexus_toid)
  donors_per_nex <- table(fpd$flowpath_toid)
  fp_inflow_nex  <- unname(inflow_nex_of[fpd$flowpath_id])
  n_up <- as.integer(donors_per_nex[fp_inflow_nex]); n_up[is.na(n_up)] <- 0L
  is_leaf <- n_up == 0L

  # downstream flowpath each stub drains to (topology proxy for the catchment
  # containing it).
  nex_to_fp <- stats::setNames(nexd$nexus_toid, nexd$nexus_id)
  target_fp <- unname(nex_to_fp[fpd$flowpath_toid])

  demote <- cand & is_leaf & !is.na(target_fp) & nzchar(target_fp)
  # never re-point into another demoted stub (targets carry divides, so this is a
  # belt-and-braces guard against chains)
  demote <- demote & !(target_fp %in% fpd$flowpath_id[which(demote)])

  # Terminal divide-less leaves (no downstream catchment via topology) are
  # resolved by SPATIAL containment, the demotion's actual intent ("the
  # catchment containing it"); topology is only a proxy that breaks at terminals.
  # Fold the stub into the divide its flowline mostly lies within. A stub outside
  # every divide (over-water coastal connector: no catchment, no upstream,
  # terminal) is an isolated artifact -> drop it (area- & connectivity-neutral).
  # demote_drop is retained for the action-section plumbing but is intentionally
  # never populated: dropping a divide-less leaf would discard its flowline AND
  # its reference_id from the (flowline-grain) network -> network_covers failure.
  # A reference-bearing line is never dropped; it always folds into a catchment.
  demote_drop <- character(0)
  term_i <- which(cand & is_leaf & (is.na(target_fp) | !nzchar(target_fp)))
  if (length(term_i)) {
    dv_sp <- tryCatch(sf::st_read(gpkg, "divides", quiet = TRUE),
      error = function(e) NULL)
    if (!is.null(dv_sp) && nrow(dv_sp) > 0L && "flowpath_id" %in% names(dv_sp)) {
      sg       <- sf::st_transform(fp_sf[term_i, ], sf::st_crs(dv_sp))
      ixl      <- sf::st_intersects(sg, dv_sp)
      need_near <- integer(0)
      for (k in seq_along(term_i)) {
        cont <- ixl[[k]]
        if (length(cont) == 0L) { need_near <- c(need_near, k); next }
        len <- suppressWarnings(as.numeric(sf::st_length(
          sf::st_intersection(sg[k, ], dv_sp[cont, ]))))
        target_fp[term_i[k]] <- as.character(dv_sp$flowpath_id[cont[which.max(len)]])
        demote[term_i[k]]    <- TRUE
      }
      # Stubs lying outside EVERY divide (coastal / over-water connectors) fold
      # into the NEAREST divide: preserves the flowline + reference_id in the
      # network (never drop a reference-bearing line) and gives the connector a
      # sensible containing catchment for routing.
      if (length(need_near)) {
        nf <- sf::st_nearest_feature(sg[need_near, ], dv_sp)
        target_fp[term_i[need_near]] <- as.character(dv_sp$flowpath_id[nf])
        demote[term_i[need_near]]    <- TRUE
      }
    }
    # No divides layer at all (should never happen): leave stubs as-is rather
    # than drop them; no_divideless would flag it, surfacing the real problem.
  }

  # hard lake guard: a stub touching a lake is a lake connector: never demote
  # (fold) or drop it.
  guard_i <- union(which(demote),
    which(as.character(fpd$flowpath_id) %in% demote_drop))
  if (length(guard_i) && !is.null(lakes_path) && file.exists(lakes_path) &&
    hfutils::layer_exists(lakes_path, "lakes")) {
    lk <- tryCatch(sf::st_read(lakes_path, "lakes", quiet = TRUE),
      error = function(e) NULL)
    if (!is.null(lk) && nrow(lk) > 0L) {
      cg  <- sf::st_transform(fp_sf[guard_i, ], sf::st_crs(lk))
      hit <- lengths(sf::st_intersects(cg, lk)) > 0L
      demote[guard_i[hit]] <- FALSE
      demote_drop <- setdiff(demote_drop,
        as.character(fpd$flowpath_id[guard_i[hit]]))
    }
  }

  if (!any(demote) && length(demote_drop) == 0L) {
    message(sprintf(
      "[stage3][%s] Step 6k: no divide-less leaf flowpaths to demote", vpuid))
    return(invisible(0L))
  }

  di        <- which(demote)
  demote_fp <- fpd$flowpath_id[di]
  src_nex   <- fpd$flowpath_toid[di]            # nexus each stub drains to
  tgt_fp    <- target_fp[di]
  fp_toid   <- stats::setNames(fpd$flowpath_toid, fpd$flowpath_id)
  map       <- data.frame(src = demote_fp, tgt = tgt_fp,
    tgt_toid = unname(fp_toid[tgt_fp]), stringsAsFactors = FALSE)
  mp_fp   <- stats::setNames(map$tgt,      map$src)
  mp_toid <- stats::setNames(map$tgt_toid, map$src)

  # nexuses orphaned by the demotion: a stub's target nexus when no SURVIVING
  # flowpath drains to it, plus the stub's own (empty) inflow nexus. Accounts for
  # both folded (di) and dropped (drop_di) stubs.
  drop_di    <- which(as.character(fpd$flowpath_id) %in% demote_drop)
  gone_i     <- c(di, drop_di)
  remaining_toid <- if (length(gone_i)) fpd$flowpath_toid[-gone_i] else fpd$flowpath_toid
  stub_nex   <- c(src_nex, fpd$flowpath_toid[drop_di])   # nexuses stubs drained to
  drop_nex <- c(stub_nex[!(stub_nex %in% remaining_toid)],
    unname(inflow_nex_of[c(demote_fp, demote_drop)]))
  drop_nex <- unique(drop_nex[!is.na(drop_nex) & nzchar(drop_nex)])

  # generic re-point helper for a geometry/attribute layer carrying flowpath_id
  # (+ optional flowpath_toid): rewrite the demoted stubs' parent pointers.
  .repoint_layer <- function(layer) {
    if (!hfutils::layer_exists(gpkg, layer)) return(invisible(NULL))
    x <- sf::st_read(gpkg, layer, quiet = TRUE)
    if (!("flowpath_id" %in% names(x))) return(invisible(NULL))
    fid  <- as.character(x$flowpath_id)
    hit  <- fid %in% demote_fp      # folded stubs -> re-point parent pointers
    drop <- fid %in% demote_drop    # over-water artifacts -> remove rows
    if (!any(hit) && !any(drop)) return(invisible(NULL))
    x$flowpath_id      <- fid
    x$flowpath_id[hit] <- unname(mp_fp[fid[hit]])
    if ("flowpath_toid" %in% names(x)) {
      x$flowpath_toid      <- as.character(x$flowpath_toid)
      x$flowpath_toid[hit] <- unname(mp_toid[fid[hit]])
    }
    if (any(drop)) x <- x[!drop, ]
    sf::st_write(.drop_list_cols(x), gpkg, layer, delete_layer = TRUE,
      quiet = TRUE)
  }

  # --- 1. flowpaths + flowpath-attributes: drop folded AND artifact stubs ----
  gone_fp <- c(demote_fp, demote_drop)
  keep <- !(as.character(fp_sf$flowpath_id) %in% gone_fp)
  sf::st_write(.drop_list_cols(fp_sf[keep, ]), gpkg, "flowpaths",
    delete_layer = TRUE, quiet = TRUE)
  if (hfutils::layer_exists(gpkg, "flowpath-attributes")) {
    fa  <- sf::st_read(gpkg, "flowpath-attributes", quiet = TRUE)
    if ("flowpath_id" %in% names(fa))
      sf::st_write(.drop_list_cols(
        fa[!(as.character(fa$flowpath_id) %in% gone_fp), ]),
        gpkg, "flowpath-attributes", delete_layer = TRUE, quiet = TRUE)
  }

  # --- 2. nexus: drop orphaned / emptied nexuses -----------------------------
  if (length(drop_nex)) {
    nkeep <- !(as.character(nex_sf$nexus_id) %in% drop_nex)
    sf::st_write(.drop_list_cols(nex_sf[nkeep, ]), gpkg, "nexus",
      delete_layer = TRUE, quiet = TRUE)
  }

  # --- 3. flowlines + flowline-attributes: re-point member parent pointers,
  .repoint_layer("flowlines")
  .repoint_layer("flowline-attributes")

  # --- 4. network: re-point + copy the target's flowpath-scale columns -------
  # `network` is fully denormalized by this stage, so the stub's rows carry the
  # OLD flowpath's scale attributes. Copy every flowpath-grain column from a
  # prototype row of the target (its member rows all share these), keep the
  # flowline-grain columns, and mark the folded members off-mainstem.
  if (hfutils::layer_exists(gpkg, "network")) {
    con <- DBI::dbConnect(RSQLite::SQLite(), gpkg)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    net <- DBI::dbReadTable(con, "network")
    # dbReadTable returns the GPKG integer PK as a regular `fid` column;
    # dbWriteTable(overwrite) would then persist it as data, colliding with the
    # auto-PK on a later ogr2ogr -append (stage-4 merge). Drop it so the rewritten
    # table keeps a clean implicit rowid PK.
    net$fid <- NULL
    nid  <- as.character(net$flowpath_id)
    hit  <- nid %in% demote_fp      # folded -> re-point to target catchment
    drop <- nid %in% demote_drop    # over-water artifact -> remove rows
    if (any(hit) || any(drop)) {
      net$flowpath_id <- nid
      proto   <- net[!duplicated(nid), , drop = FALSE]   # one row per flowpath
      tgt_row <- unname(mp_fp[nid[hit]])
      pidx    <- match(tgt_row, as.character(proto$flowpath_id))
      fg_cols <- intersect(c("mainstem_id", "hydroseq", "areasqkm", "lengthkm",
        "total_dasqkm", "slope", "divide_id", "divide_areasqkm", "poi_id",
        "lake_id", "hl_reference", "hydrolocation_id"), names(net))
      for (col in fg_cols) net[[col]][hit] <- proto[[col]][pidx]
      net$flowpath_id[hit] <- tgt_row
      if ("flowpath_toid" %in% names(net))
        net$flowpath_toid[hit] <- unname(mp_toid[nid[hit]])
      if ("is_mainstem" %in% names(net)) net$is_mainstem[hit] <- 0L
      if (any(drop)) net <- net[!drop, ]
      DBI::dbWriteTable(con, "network", net, overwrite = TRUE)
    }
  }

  message(sprintf(paste0(
    "[stage3][%s] Step 6k: demoted %d divide-less leaf flowpath(s) to ",
    "flowline-only (folded into %d catchment(s)); dropped %d over-water ",
    "artifact(s); removed %d nexus(es)"),
    vpuid, length(demote_fp), length(unique(tgt_fp)), length(demote_drop),
    length(drop_nex)))
  invisible(length(demote_fp) + length(demote_drop))
}
