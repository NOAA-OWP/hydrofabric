# reconcile_divides.R
# Best-of-both divide reconciliation:
#   - Respects FAC/FDR raster CRS (no hardcoded 5070)
#   - fast_* column schema throughout (flowpath_id / member_flowpath_id)
#   - lr = FALSE flag for optional left/right bank split via lwgeom::st_split
#
# Entry point: reconcile_divides()

# -- reconcile_divides ------------------------------------------------------------

#' Reconcile divides to refactored flowlines (CRS-aware, lr-capable)
#'
#' @description
#' Raster-aware divide reconciliation with two key properties:
#'
#' (1) CRS from rasters: all metric spatial operations (gap fill, area
#'     calculations, MULTIPOLYGON cleanup) use the native CRS of the FDR
#'     raster instead of a hardcoded EPSG:5070.  Required for datasets
#'     outside CONUS (e.g. PRVI, Alaska).
#'
#' (2) lr flag: when lr = TRUE, each reconciled divide polygon is further
#'     split left/right along its associated flowline using lwgeom::st_split,
#'     matching the reference split_catchment_divide(lr = TRUE) approach.
#'     Output flowpath_id values are suffixed "_L" and "_R".
#'
#' @param flowline sf LINESTRING. Refactored flowlines with flowpath_id and
#'   flowpath_toid columns. Typically out$refactored from refactor_flowlines().
#' @param divides sf POLYGON. Input divide polygons with flowpath_id.
#' @param fdr_path Character. Path to D8 flow-direction raster (required).
#' @param fac_path Character. Path to flow-accumulation raster (required).
#' @param lookup_table The out$reconciled sf from refactor_flowlines(), or NULL.
#'   When provided (must have ID and member_flowpath_id columns) it is used
#'   for reconciled-aware deduplication before splitting.
#' @param lr Logical. If TRUE, apply left/right bank split along the flowline
#'   after D8 reconciliation. Requires lwgeom. Default FALSE.
#' @param lr_densify_m Numeric. dfMaxLength for st_segmentize when densifying
#'   divide boundary before lr split. Default 100. Units match raster CRS.
#' @param lr_snap_thresh_m Numeric. Distance threshold (raster CRS units) for
#'   snapping flowline endpoints to divide boundary. Default 100.
#' @param fallback_fdr_path,fallback_fac_path Optional paths to fallback
#'   FDR / FAC rasters, used when the primary rasters return NoData at a split
#'   point. `NULL` disables the fallback (default); ignored when equal to the
#'   primary paths.
#' @param parallel Logical. Enable parallel multi-parent divide splitting. Default TRUE.
#' @param num_workers Integer or NULL. Parallel workers. NULL uses 4.
#'
#' @return sf POLYGON with columns divide_id, flowpath_id, and geometry.
#'   When lr = TRUE, flowpath_id values are suffixed "_L" / "_R".
#'
#' @seealso refactor_flowlines()
#' @importFrom sf st_crs st_transform st_as_sf st_geometry st_union
#' @importFrom lwgeom lwgeom_make_valid
#'   st_difference st_cast st_sf st_sfc st_is_empty st_area st_length
#'   st_intersects st_nearest_feature st_set_geometry st_geometry_type
#'   st_intersection st_boundary st_coordinates st_linestring st_segmentize
#'   st_as_sfc st_make_valid
#' @importFrom terra rast crs crop buffer vect writeRaster
#' @importFrom hfutils union_linestrings
#' @importFrom logger log_info log_warn
#' @importFrom data.table as.data.table
#' @importFrom geos as_geos_geometry geos_make_collection geos_unary_union
#' @examples
#' \dontrun{
#' out <- refactor_flowlines(flines, split_flines_meters = 2000,
#'                           collapse_flines_meters = 1000,
#'                           out_refactored = "refactored.gpkg",
#'                           out_reconciled = "reconciled.gpkg")
#' rec_div <- reconcile_divides(flowline      = out$refactored,
#'                              divides       = divides_sf,
#'                              fdr_path      = "fdr.tif",
#'                              fac_path      = "fac.tif",
#'                              lookup_table  = out$reconciled)
#' }
#' @export
reconcile_divides <- function(flowline,
                              divides,
                              fdr_path,
                              fac_path,
                              fallback_fdr_path  = NULL,
                              fallback_fac_path  = NULL,
                              lookup_table       = NULL,
                              lr                 = FALSE,
                              lr_densify_m       = 100,
                              lr_snap_thresh_m   = 100,
                              parallel           = TRUE,
                              num_workers        = NULL) {

  if (!nzchar(fdr_path) || !.raster_accessible(fdr_path))
    stop("fdr_path is required and must exist: ", fdr_path)
  if (!nzchar(fac_path) || !.raster_accessible(fac_path))
    stop("fac_path is required and must exist: ", fac_path)
  if (is.null(flowline) || nrow(flowline) == 0) stop("flowline is empty")
  if (is.null(divides)  || nrow(divides)  == 0) stop("divides is empty")

  if (lr && !requireNamespace("lwgeom", quietly = TRUE))
    stop("lr = TRUE requires the lwgeom package: install.packages('lwgeom')")

  # -- 1. Detect raster CRS, used for ALL metric spatial ops -----------------
  fdr_hdr <- terra::rast(fdr_path)
  rast_crs <- sf::st_crs(terra::crs(fdr_hdr))
  rm(fdr_hdr)
  if (is.na(rast_crs) || is.null(rast_crs$wkt))
    stop("Unable to determine CRS from FDR raster: ", fdr_path)
  logger::log_info("reconcile_divides: raster CRS = '{rast_crs$input}'")

  in_crs <- sf::st_crs(divides)

  # Keep the original input divides for post-hoc gap fill.
  orig_input_divides <- divides

  # -- 2. Reconciled-aware deduplication (same as fast_reconcile_divides) -----
  divide_id_set   <- unique(as.character(divides$flowpath_id))
  divide_id_set   <- divide_id_set[!is.na(divide_id_set) & nzchar(divide_id_set)]
  rfc_ids_in_fl   <- as.character(flowline$flowpath_id)

  use_reconciled_dedup <- !is.null(lookup_table) &&
    inherits(lookup_table, "data.frame") &&
    all(c("ID", "member_flowpath_id") %in% names(lookup_table))

  orig_reconciled_lt <- if (use_reconciled_dedup) {
    if (inherits(lookup_table, "sf")) sf::st_drop_geometry(lookup_table)
    else lookup_table
  } else NULL

  if (use_reconciled_dedup) {
    rec_df <- orig_reconciled_lt

    member_map <- do.call(rbind, lapply(seq_len(nrow(rec_df)), function(i) {
      members <- trimws(strsplit(as.character(rec_df$member_flowpath_id[i]), ",")[[1]])
      members <- members[nzchar(members) & members %in% rfc_ids_in_fl]
      if (!length(members)) return(NULL)
      data.frame(reconciled_id         = rec_df$ID[i],
        refactored_flowpath_id = members,
        stringsAsFactors       = FALSE)
    }))

    if (!is.null(member_map) && nrow(member_map) > 0) {
      base_ids <- sub("[.][0-9]+$", "", member_map$refactored_flowpath_id)
      member_map$parent_divide_id <- ifelse(base_ids %in% divide_id_set,
        base_ids,
        member_map$refactored_flowpath_id)

      lkm_idx             <- match(member_map$refactored_flowpath_id, rfc_ids_in_fl)
      member_map$lkm_sort <- sf::st_drop_geometry(flowline)$lengthkm[lkm_idx]
      member_map$lkm_sort[is.na(member_map$lkm_sort)] <- 0
      member_map          <- member_map[order(-member_map$lkm_sort), , drop = FALSE]
      dedup_key           <- paste(member_map$refactored_flowpath_id,
        member_map$parent_divide_id, sep = "|||")
      member_map          <- member_map[!duplicated(dedup_key), , drop = FALSE]
      member_map$lkm_sort <- NULL

      flowline     <- flowline[rfc_ids_in_fl %in% member_map$refactored_flowpath_id, ,
        drop = FALSE]
      lookup_table <- data.frame(
        refactored_flowpath_id = member_map$refactored_flowpath_id,
        parent_divide_id       = member_map$parent_divide_id,
        stringsAsFactors       = FALSE
      )
      logger::log_info(
        "reconcile_divides: dedup -> {nrow(flowline)} representative flowlines (was {length(rfc_ids_in_fl)})")
    } else {
      use_reconciled_dedup <- FALSE
    }
  }

  if (!use_reconciled_dedup) {
    rfc_ids    <- unique(rfc_ids_in_fl[!is.na(rfc_ids_in_fl) & nzchar(rfc_ids_in_fl)])
    base_ids   <- sub("[.][0-9]+$", "", rfc_ids)
    parent_ids <- ifelse(base_ids %in% divide_id_set, base_ids, rfc_ids)
    lookup_table <- data.frame(
      refactored_flowpath_id = rfc_ids,
      parent_divide_id       = parent_ids,
      stringsAsFactors       = FALSE
    )
  }

  # -- 3. Pre-crop: one small raster tile for all multi-parent divides ---------
  fdr_work      <- fdr_path
  fac_work      <- fac_path
  precrop_files <- character(0)
  t_precrop     <- proc.time()[["elapsed"]]

  tryCatch(
    {
      child_map <- .prd_build_child_parent_map(
        flowline     = flowline,
        divides      = divides,
        lookup_table = lookup_table
      )
      cm_dt     <- data.table::as.data.table(child_map)
      multi_pid <- cm_dt[, .N, by = "parent_divide_id"][N > 1L, as.character(parent_divide_id)]

      if (length(multi_pid) > 0) {
        multi_div <- divides[as.character(divides$flowpath_id) %in% multi_pid, , drop = FALSE]
        if (nrow(multi_div) > 0) {
          fdr_r    <- terra::rast(fdr_path)
          fac_r    <- terra::rast(fac_path)
          bb       <- sf::st_bbox(sf::st_transform(multi_div, rast_crs))
          bbox_ext <- terra::extend(
            terra::ext(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]),
            500)

          # Guard: when multi-parent divides are scattered across the VPU their
          # combined bbox spans the whole domain (billions of cells). Materializing
          # that as a single tile OOMs, and most of it is empty. Skip the pre-crop
          # tile above a cell-count cap and let process_parent crop small per-parent
          # windows straight from the VRT instead.
          res_xy    <- terra::res(fdr_r)
          est_cells <- (as.numeric(terra::xmax(bbox_ext) - terra::xmin(bbox_ext)) / res_xy[1]) *
            (as.numeric(terra::ymax(bbox_ext) - terra::ymin(bbox_ext)) / res_xy[2])
          precrop_cell_cap <- 2.5e8

          if (is.finite(est_cells) && est_cells > precrop_cell_cap) {
            logger::log_info(
              "reconcile_divides: pre-crop skipped -- combined multi-parent bbox ~{format(round(est_cells), big.mark=',')} cells > {format(precrop_cell_cap, scientific=FALSE, big.mark=',')}; using full rasters with per-parent crops")
          } else {
            fdr_crop <- terra::crop(fdr_r, bbox_ext, snap = "out")
            fac_crop <- terra::crop(fac_r, bbox_ext, snap = "out")

            fdr_work <- tempfile(pattern = "fdr_pre_", fileext = ".tif")
            fac_work <- tempfile(pattern = "fac_pre_", fileext = ".tif")
            terra::writeRaster(fdr_crop, fdr_work, overwrite = TRUE,
              gdal = c("COMPRESS=NONE", "BIGTIFF=IF_NEEDED"))
            terra::writeRaster(fac_crop, fac_work, overwrite = TRUE,
              gdal = c("COMPRESS=NONE", "BIGTIFF=IF_NEEDED"))
            precrop_files <- c(fdr_work, fac_work)
            on.exit(unlink(precrop_files[file.exists(precrop_files)]), add = TRUE)

            logger::log_info(paste0(
              "reconcile_divides: pre-cropped FDR/FAC to {length(multi_pid)} ",
              "multi-parent divide(s) ",
              "[{sprintf('%.2f', proc.time()[['elapsed']] - t_precrop)}s]"))
          }
        }
      } else {
        logger::log_info("reconcile_divides: no multi-parent divides; skipping pre-crop")
      }
    },
    error = function(e) {
      logger::log_warn("reconcile_divides: pre-crop failed ({conditionMessage(e)}); using full rasters")
      # Reset work paths so partition step uses the originals, not a stale temp path
      fdr_work <<- fdr_path
      fac_work <<- fac_path
    })

  # -- 4. Core partitioning ----------------------------------------------------
  result <- partition_reconcile_divides_from_objects(
    flowline          = flowline,
    divides           = divides,
    lookup_table      = lookup_table,
    fdr_path          = fdr_work,
    fac_path          = fac_work,
    fallback_fdr_path = fallback_fdr_path,
    fallback_fac_path = fallback_fac_path,
    write_outputs     = FALSE,
    num_workers       = num_workers
  )
  divides_out <- result$divides

  # -- 5. dissolve_multipart_divides ------------------------------------------
  # Replaces hfutils::clean_geometry for intra-divide cleanup. Inter-divide
  # topology (overlaps, gaps) is handled downstream by clean_gpkg_layer.
  if (!is.null(divides_out) && nrow(divides_out) > 0) {
    tryCatch(
      {
        na_fp <- is.na(divides_out$flowpath_id)
        # Remember which divides arrived with no flowpath before the id is
        # filled in below. A landscape divide is defined by having no flowpath,
        # and once flowpath_id has been back-filled from divide_id there is no
        # way to tell it apart from a network divide again.
        divides_out$.hf_landscape <- na_fp
        if (any(na_fp) && "divide_id" %in% names(divides_out))
          divides_out$flowpath_id[na_fp] <- as.character(divides_out$divide_id[na_fp])
        # Validity repair + multi-part reduction.
        sf::st_geometry(divides_out) <- lwgeom::lwgeom_make_valid(sf::st_geometry(divides_out))
        divides_out <- dissolve_multipart_divides(divides_out, id_col = "flowpath_id")
        divides_out <- sf::st_as_sf(divides_out)
        if (!"divide_id" %in% names(divides_out))
          divides_out$divide_id <- as.character(divides_out$flowpath_id)
        logger::log_info("reconcile_divides: dissolve_multipart_divides complete")
      },
      error = function(e) {
        logger::log_warn("reconcile_divides: dissolve_multipart_divides failed ({conditionMessage(e)})")
      })
  }

  # -- 6. Reconciled-group union (1:1 divide per reconciled flowpath) ----------
  # Build a deterministic mapping: each divide's flowpath_id -> exactly ONE
  # reconciled group ID. For each group, register only the REFACTORED flowpath
  # IDs (those that actually appear in out$refactored) as members. Parent-expansion
  # IDs are excluded: they correspond to unsplit originals already consumed by the
  # FAC/FDR partition in Step 4, and a divide matching multiple groups via a
  # parent-expansion id would union into the wrong one.
  if (!is.null(orig_reconciled_lt) && !is.null(divides_out) && nrow(divides_out) > 0) {
    div_fp  <- as.character(divides_out$flowpath_id)
    out_crs <- sf::st_crs(divides_out)

    # Build fp_id -> group_id mapping (first-match wins for duplicates)
    # data.table string-split avoids two nested R for-loops over ~11k rows.
    rfc_fp_set <- hfutils::hf_fmt_id(sf::st_drop_geometry(flowline)$flowpath_id)
    rec_lt_dt  <- data.table::as.data.table(
      sf::st_drop_geometry(orig_reconciled_lt)[, c("ID", "member_flowpath_id")])
    rec_lt_dt[, ID                := hfutils::hf_fmt_id(ID)]
    rec_lt_dt[, member_flowpath_id := as.character(member_flowpath_id)]
    members_dt <- rec_lt_dt[,
      {
        mems <- trimws(strsplit(member_flowpath_id, ",", fixed = TRUE)[[1]])
        mems <- mems[nzchar(mems)]
        .(member = mems)
      },
      by = ID]
    # Only refactored IDs; first row wins for any duplicate member
    members_dt <- members_dt[member %in% rfc_fp_set][!duplicated(member)]
    fp_to_group <- stats::setNames(members_dt$ID, members_dt$member)

    # Map each divide to its group via its flowpath_id
    div_group <- fp_to_group[div_fp]

    # Spatial recovery for divides whose flowpath_id names no reconciled group.
    # Match them to a group by finding the flowline most-contained in the divide's
    # polygon.
    #
    # A landscape divide, flowpath_id NA, has_flowpath 0, is excluded. It has
    # no flowpath by definition, so there is no lost id link to recover, and
    # "most-contained" answers the wrong question for one: a coastal divide is a
    # ring around the catchments it encloses and contains dozens of flowlines, so
    # handing the whole ring to one of them claims land that reach does not
    # drain. On VPU 01 that put 71 interior rings into a single network divide,
    # against a maximum of 1 without it, and failed no_excessive_rings; on VPU 04
    # it moved 3,169 of 7,935 km2 of landscape area into network divides.
    # Landscape divides are carried through as `lnd-` rows instead, which is what
    # they are and what the delivered non-dendritic products already hold.
    all_gids     <- as.character(orig_reconciled_lt$ID)
    group_ids    <- unique(stats::na.omit(div_group))
    missing_gids <- setdiff(all_gids, group_ids)
    is_landscape <- if (".hf_landscape" %in% names(divides_out))
      !is.na(divides_out$.hf_landscape) & divides_out$.hf_landscape else
      is.na(div_fp) | !nzchar(as.character(div_fp))
    na_idx <- which(is.na(div_group) & !is_landscape)
    if (length(na_idx) > 0 && length(missing_gids) > 0 && !is.null(flowline)) {
      mems_per_group <- split(names(fp_to_group), unname(fp_to_group))
      missing_mems   <- unlist(mems_per_group[missing_gids], use.names = FALSE)
      fl_fp          <- as.character(sf::st_drop_geometry(flowline)$flowpath_id)
      fl_sub         <- flowline[fl_fp %in% missing_mems, ]
      if (nrow(fl_sub) > 0) {
        fl_geom_sub <- sf::st_geometry(sf::st_transform(fl_sub, out_crs))
        fl_to_gid   <- unname(fp_to_group[
          as.character(sf::st_drop_geometry(fl_sub)$flowpath_id)])
        na_geoms    <- sf::st_geometry(divides_out[na_idx, ])
        ints        <- sf::st_intersects(na_geoms, fl_geom_sub, sparse = TRUE)

        # Flatten (div, fl) candidate pairs filtered by missing_gids, compute
        # all intersections in a single vectorized geos call, pick per-div
        # the fl whose intersection length is longest. Replaces the prior
        # doubly-nested loop (~O(n_divs x candidates) serial sf calls).
        miss_set <- missing_gids
        pair_di  <- integer(0)
        pair_fi  <- integer(0)
        for (di in seq_along(ints)) {
          cand <- ints[[di]]
          if (length(cand) == 0L) next
          cand <- cand[!is.na(fl_to_gid[cand]) & fl_to_gid[cand] %in% miss_set]
          if (length(cand) == 0L) next
          pair_di <- c(pair_di, rep.int(di, length(cand)))
          pair_fi <- c(pair_fi, cand)
        }

        if (length(pair_di) > 0L) {
          na_geos <- geos::as_geos_geometry(na_geoms)
          fl_geos <- geos::as_geos_geometry(fl_geom_sub)
          inter_g <- geos::geos_intersection(na_geos[pair_di], fl_geos[pair_fi])
          lens    <- suppressWarnings(as.numeric(geos::geos_length(inter_g)))
          lens[is.na(lens)] <- -1

          # Per-div argmax via data.table (tapply + match is O(n log n) and
          # keeps one call). Use base R to avoid a hard dep if data.table
          # isn't present in this code path's env.
          order_ix  <- order(pair_di, -lens)
          pair_di_s <- pair_di[order_ix]
          pair_fi_s <- pair_fi[order_ix]
          first_ix  <- !duplicated(pair_di_s)
          best_di   <- pair_di_s[first_ix]
          best_fi   <- pair_fi_s[first_ix]

          div_group[na_idx[best_di]] <- fl_to_gid[best_fi]
          if (length(best_di) > 0L) logger::log_info(
            "reconcile_divides: spatially recovered {length(best_di)} NA-id reference divide(s)")
        }
      }
    }
    group_ids    <- unique(stats::na.omit(div_group))
    missing_gids <- setdiff(all_gids, group_ids)
    if (length(missing_gids) > 0) {
      logger::log_warn(
        "reconcile_divides: {length(missing_gids)} reconciled group(s) had no matching divide")
    }

    # Union divides per group: split geometry indices by group, union each
    # subset, then build ONE sf from parallel vectors. Avoids the O(n^2)
    # schema-check overhead of iterative single-row rbind.

    has_group  <- !is.na(div_group)
    unified_sf <- NULL
    if (any(has_group)) {
      grp_vals  <- div_group[has_group]
      geom_list <- sf::st_geometry(divides_out)[has_group]
      # Exclude NA/empty geometries (PRD may return NA-polygon rows for
      # flowlines with no matching divide).
      non_empty <- !sf::st_is_empty(geom_list) & !is.na(sf::st_dimension(geom_list))
      idx_by_grp  <- split(seq_along(grp_vals)[non_empty], grp_vals[non_empty])
      grp_names   <- names(idx_by_grp)
      geom_geos   <- geos::as_geos_geometry(geom_list)
      union_geoms <- lapply(idx_by_grp, function(idx) {
        g <- geom_geos[idx]
        r <- if (length(g) == 1L) g
        else geos::geos_unary_union(geos::geos_make_collection(g))
        sf::st_as_sfc(r)[[1L]]
      })
      unified_sf  <- sf::st_sf(
        divide_id   = grp_names,
        flowpath_id = grp_names,
        geometry    = sf::st_sfc(union_geoms, crs = out_crs))
    }

    # Post-union MULTIPOLYGON fix: keep largest component, release rest to
    # gap fill. Vectorized: flatten every MULTIPOLYGON into its component
    # POLYGONs once, compute all areas in one st_area call, then pick the
    # largest component per source row. Skip the transform when unified_sf
    # is already in a projected metres CRS.
    if (!is.null(unified_sf) && nrow(unified_sf) > 0) {
      tryCatch(
        {
          multi_mask <- sf::st_geometry_type(unified_sf) == "MULTIPOLYGON"
          if (any(multi_mask)) {
            u_crs      <- sf::st_crs(unified_sf)
            mi_idx     <- which(multi_mask)
            multi_geoms <- sf::st_geometry(unified_sf)[multi_mask]
            # Each MULTIPOLYGON sfg is a list of polygon ring-lists.
            comps_per_multi <- lapply(multi_geoms,
              function(g) lapply(g, sf::st_polygon))
            group_vec  <- rep(seq_along(mi_idx), lengths(comps_per_multi))
            flat_comps <- unlist(comps_per_multi, recursive = FALSE,
              use.names = FALSE)
            flat_sfc   <- sf::st_sfc(flat_comps, crs = u_crs)
            area_sfc   <- if (!sf::st_is_longlat(unified_sf)) flat_sfc
            else sf::st_transform(flat_sfc, rast_crs)
            areas      <- suppressWarnings(as.numeric(sf::st_area(area_sfc)))
            ord            <- order(group_vec, -areas)
            first_per_grp  <- !duplicated(group_vec[ord])
            best_flat_ix   <- ord[first_per_grp]
            ugeom_list     <- as.list(sf::st_geometry(unified_sf))
            for (k in seq_along(mi_idx))
              ugeom_list[[mi_idx[k]]] <- flat_comps[[best_flat_ix[k]]]
            unified_sf <- sf::st_set_geometry(
              unified_sf,
              sf::st_sfc(ugeom_list, crs = u_crs))
            logger::log_info(
              "reconcile_divides: post-union MULTIPOLYGON fix: {sum(multi_mask)} group(s)")
          }
        },
        error = function(e) {
          logger::log_warn("reconcile_divides: post-union MULTIPOLYGON fix failed ({conditionMessage(e)})")
        })
    }

    # -- Post-union degenerate interior-ring cleanup ----------------------------
    # geos_unary_union of a reconciled group's member catchments can emit
    # zero-area (degenerate) interior rings at shared edges, a floating-point
    # topology artifact, NOT a real enclosed catchment. The source reference
    # divides are valid with only genuine, positive-area holes, so drop interior
    # rings below a sub-epsilon area; every real hole (>> 1 m2) is preserved.
    # Reference fabrics are in projected metres, so the shoelace area is exact;
    # skip (can't threshold in m2) if somehow long/lat.
    if (!is.null(unified_sf) && nrow(unified_sf) > 0 &&
        !sf::st_is_longlat(unified_sf)) {
      ring_eps_m2 <- 1
      .ring_area <- function(m) {                 # shoelace, projected metres
        n <- nrow(m); if (n < 4L) return(0)
        abs(sum(m[-n, 1L] * m[-1L, 2L] - m[-1L, 1L] * m[-n, 2L])) / 2
      }
      .clean_rings <- function(p) {               # p: list of ring matrices
        if (length(p) <= 1L) return(p)
        keep <- c(TRUE, vapply(p[-1L],
          function(r) .ring_area(r) >= ring_eps_m2, logical(1)))
        p[keep]
      }
      glist     <- as.list(sf::st_geometry(unified_sf))
      n_dropped <- 0L
      for (i in seq_along(glist)) {
        gi <- glist[[i]]
        if (inherits(gi, "POLYGON") && length(gi) > 1L) {
          cleaned <- .clean_rings(unclass(gi))
          if (length(cleaned) < length(gi)) {
            n_dropped  <- n_dropped + (length(gi) - length(cleaned))
            glist[[i]] <- sf::st_polygon(cleaned)
          }
        } else if (inherits(gi, "MULTIPOLYGON")) {
          parts  <- lapply(unclass(gi), .clean_rings)
          before <- sum(vapply(unclass(gi), length, 0L))
          after  <- sum(vapply(parts, length, 0L))
          if (after < before) {
            n_dropped  <- n_dropped + (before - after)
            glist[[i]] <- sf::st_multipolygon(parts)
          }
        }
      }
      if (n_dropped > 0L) {
        unified_sf <- sf::st_set_geometry(unified_sf,
          sf::st_sfc(glist, crs = sf::st_crs(unified_sf)))
        logger::log_info(
          "reconcile_divides: dropped {n_dropped} degenerate interior ring(s) (<{ring_eps_m2} m2)")
      }
    }

    # Add orig_divide_id: traces each reconciled divide back to the original
    # NHD catchment COMID, following the same lineage pattern as orig_levelpathID
    # on reconciled flowpaths.  Built by inverting fp_to_group through lookup_table.
    if (!is.null(unified_sf) && nrow(unified_sf) > 0 &&
      !is.null(lookup_table) &&
      all(c("refactored_flowpath_id", "parent_divide_id") %in% names(lookup_table))) {
      # data.table join: invert fp_to_group, join with lookup_table, pick most
      # common parent per group: avoids O(n*m) vapply scan.
      fp_grp_dt <- data.table::data.table(
        member   = names(fp_to_group),
        group_id = as.character(fp_to_group))
      lt_dt <- data.table::data.table(
        member     = as.character(lookup_table$refactored_flowpath_id),
        parent_div = as.character(lookup_table$parent_divide_id))
      merged_dt  <- fp_grp_dt[lt_dt, on = "member", nomatch = 0L]
      merged_dt  <- merged_dt[!is.na(parent_div) & nzchar(parent_div)]
      orig_div_dt <- merged_dt[
        , .(orig_divide_id = names(which.max(table(parent_div)))), by = group_id]
      unified_sf$orig_divide_id <- orig_div_dt$orig_divide_id[
        match(as.character(unified_sf$divide_id), orig_div_dt$group_id)]
    } else if (!is.null(unified_sf) && nrow(unified_sf) > 0) {
      unified_sf$orig_divide_id <- NA_character_
    }

    # -- 6b. Re-home displaced partitions of a split parent catchment ----------
    # When a reference COMID is split at POI events (refactor_flowlines creates
    # ".1", ".2", ".3" sub-segments), the partition divides the parent catchment
    # among them using FAC/FDR delineation. Raster snapping and D8 connectivity
    # breaks at masked polygon boundaries can return sub-polygons that sit
    # nowhere near the fragment they were assigned to. The land is right and no
    # area is missing; the assignment of land to fragment is wrong.
    #
    # Resolve it the way the partition is itself trying to: each polygon part of
    # the displaced land goes to the fragment whose own flowline is nearest it.
    # Fragments whose partition already meets their flowline are left alone, so
    # only displaced land moves and the parent's area is conserved exactly.
    #
    # Rate is strongly domain-dependent: Alaska displaces 222 of 631 split
    # groups, VPU 04 displaces 45 of 659, VPU 18 none of 488. Domains that
    # displace carry reconciled reaches whose catchment sits kilometres from
    # their own line; domains that do not stay within a few DEM cells.
    if (!is.null(unified_sf) && nrow(unified_sf) > 0 &&
      "orig_divide_id" %in% names(unified_sf)) {
      oid_col  <- as.character(unified_sf$orig_divide_id)
      grp_cnt  <- table(oid_col[!is.na(oid_col)])
      split_oids <- names(grp_cnt[grp_cnt > 1L])
      if (length(split_oids) > 0L) {
        div_geoms <- sf::st_geometry(unified_sf)
        # Bring flowlines into the same CRS as the divide geometries so that
        # intersection tests are geometrically correct.
        fl_t       <- sf::st_transform(flowline, sf::st_crs(unified_sf))
        fl_fp_col  <- as.character(sf::st_drop_geometry(fl_t)$flowpath_id)
        n_fixed    <- 0L

        for (oid in split_oids) {
          rows <- which(!is.na(oid_col) & oid_col == oid)
          if (length(rows) < 2L) next

          # Member flowlines per row, computed once: they decide both which
          # partitions are displaced and, when they are, where the land belongs.
          # A row whose members are absent from the deduplicated flowline set
          # cannot be located, so it neither gives up land nor receives any.
          row_lines <- lapply(rows, function(r) {
            grp_id  <- as.character(unified_sf$flowpath_id[r])
            members <- names(fp_to_group)[fp_to_group == grp_id]
            fl_sub  <- fl_t[fl_fp_col %in% members, , drop = FALSE]
            if (nrow(fl_sub) == 0L) NULL else sf::st_geometry(fl_sub)
          })

          disjoint_rows <- vapply(seq_along(rows), function(k) {
            if (is.null(row_lines[[k]])) return(FALSE)
            poly_sfc <- sf::st_sfc(div_geoms[[rows[k]]], crs = sf::st_crs(fl_t))
            !any(suppressWarnings(
              lengths(sf::st_intersects(poly_sfc, row_lines[[k]]))) > 0L)
          }, logical(1L))

          if (!any(disjoint_rows)) next   # all partitions valid -- nothing to do

          # Recompute the whole parent's partition as a nearest-reach cut. The
          # trace that produced these pieces has demonstrably failed for this
          # parent, so repairing its output piece by piece inherits its errors;
          # cutting the parent afresh gives every fragment the land nearest its
          # own line, conserves area exactly and cannot starve a fragment.
          parent <- tryCatch(sf::st_union(div_geoms[rows])[[1L]],
                             error = function(e) NULL)
          alloc <- if (is.null(parent)) NULL else
            .cut_parent_by_nearest_line(parent, row_lines, sf::st_crs(fl_t))

          # Fall back to moving displaced pieces when the cut cannot be built.
          if (is.null(alloc)) alloc <- .allocate_parent_parts(
            geoms    = lapply(rows, function(r) div_geoms[[r]]),
            lines    = row_lines,
            disjoint = disjoint_rows,
            crs      = sf::st_crs(fl_t))
          if (is.null(alloc)) next

          for (k in seq_along(rows)) div_geoms[[rows[k]]] <- alloc[[k]]
          n_fixed <- n_fixed + sum(disjoint_rows)
        }

        # Rebuild sfc from scratch so the class attribute reflects actual geometry
        # types.  sfc[[i]] <- sfg does NOT update the class, so if st_union()
        # returned a MULTIPOLYGON for a group of disjoint polygons, the sfc
        # would still advertise sfc_POLYGON while containing MULTIPOLYGON
        # elements: causing sf::st_transform() to fail with a vctrs type error.
        sf::st_geometry(unified_sf) <- sf::st_sfc(
          as.list(div_geoms), crs = sf::st_crs(unified_sf))
        logger::log_info(
          "reconcile_divides: reallocated {n_fixed} displaced partition(s) to their nearest reach across {length(split_oids)} split-COMID group(s)")
      }
    }

    if (!is.null(unified_sf) && nrow(unified_sf) > 0) {
      divides_out <- sf::st_as_sf(unified_sf)
      logger::log_info(
        "reconcile_divides: unified divides -> {nrow(divides_out)} reconciled group(s)")
    }

    # Carry-through: NA-flowpath reference divides that never got matched to a
    # reconciled group. These are "landscape" catchments (coastal drainage,
    # closed basins, etc.) with no associated NHD flowline. Preserve them in
    # the output with synthetic divide_ids so downstream steps can emit them
    # as has_flowline=FALSE entries.
    # The landscape marker is internal bookkeeping, not part of the data model.
    if (".hf_landscape" %in% names(divides_out))
      divides_out$.hf_landscape <- NULL

    still_na <- is.na(div_group)
    if (any(still_na) && !is.null(divides_out)) {
      lnd_geoms <- sf::st_geometry(sf::st_as_sf(result$divides))[still_na]
      non_empty <- !sf::st_is_empty(lnd_geoms) &
        !is.na(sf::st_dimension(lnd_geoms))
      if (any(non_empty)) {
        lnd_ids <- sprintf("lnd-%d", seq_len(sum(non_empty)))
        lnd_sf  <- sf::st_sf(
          divide_id      = lnd_ids,
          flowpath_id    = NA_character_,
          orig_divide_id = NA_character_,
          geometry       = lnd_geoms[non_empty])
        for (nm in setdiff(names(divides_out), names(lnd_sf))) lnd_sf[[nm]] <- NA
        for (nm in setdiff(names(lnd_sf), names(divides_out))) divides_out[[nm]] <- NA
        lnd_sf <- lnd_sf[, names(divides_out)]
        divides_out <- rbind(divides_out, lnd_sf)
        logger::log_info(
          "reconcile_divides: carried {nrow(lnd_sf)} landscape (NA-flowpath) divide(s)")
      }
    }
  }

  # -- 7. Left/right bank split (lr = TRUE) ------------------------------------
  # Gap-fill and MULTIPOLYGON cleanup are deferred to hfutils::clean_gpkg_layer
  # (mapshaper -clean) after the divides layer is written to GeoPackage.
  # That path is faster and more robust than the equivalent pure-R loops.
  # For each reconciled divide, split the polygon L/R along its flowline using
  # lwgeom::st_split. Mirrors the reference split_lr():
  #   1. Densify the polygon boundary (st_segmentize)
  #   2. Snap the upstream and downstream flowline endpoints to the nearest
  #      boundary vertex
  #   3. Force those coordinates onto the flowline endpoints
  #   4. lwgeom::st_split(polygon, snapped_flowline)
  #   5. Output two rows per input row, flowpath_id suffixed "_L" and "_R"
  # All work is done in rast_crs then transformed back.
  if (lr && !is.null(divides_out) && nrow(divides_out) > 0) {
    logger::log_info("reconcile_divides: applying lr split to {nrow(divides_out)} divide(s)")

    # Work in rast_crs for metric snapping/densification
    div_r  <- sf::st_transform(divides_out, rast_crs)
    fl_r   <- sf::st_transform(flowline,    rast_crs)

    lr_rows <- lapply(seq_len(nrow(div_r)), function(i) {
      fp_id <- as.character(div_r$flowpath_id[i])
      fl_i  <- fl_r[as.character(fl_r$flowpath_id) == fp_id, , drop = FALSE]

      # No matching flowline: return unchanged
      if (nrow(fl_i) == 0L) return(div_r[i, , drop = FALSE])

      tryCatch(
        {
          cat_poly <- sf::st_geometry(div_r)[[i]]
          cat_sfc  <- sf::st_sfc(cat_poly, crs = rast_crs)

          # Use the first (longest) flowline piece if there are multiple
          line_raw <- suppressWarnings(
            sf::st_cast(sf::st_geometry(fl_i)[[1L]], "LINESTRING"))

          # Densify polygon boundary
          cat_dense <- sf::st_sfc(
            sf::st_segmentize(cat_poly, dfMaxLength = lr_densify_m),
            crs = rast_crs)
          bnd_pts <- suppressWarnings(sf::st_cast(cat_dense, "POINT"))

          # Extract endpoint coordinates
          lm <- sf::st_coordinates(line_raw)[, 1:2, drop = FALSE]
          un <- sf::st_sfc(sf::st_point(lm[1L, ]), crs = rast_crs)
          dn <- sf::st_sfc(sf::st_point(lm[nrow(lm), ]), crs = rast_crs)

          snap_un <- bnd_pts[[sf::st_nearest_feature(un, bnd_pts)]]
          snap_dn <- bnd_pts[[sf::st_nearest_feature(dn, bnd_pts)]]

          # If upstream snap is far away, extend line upstream via trace (skip
          # the full FDR trace for compactness; fall back on a simpler re-snap)
          dist_un <- as.numeric(sf::st_distance(un,
            sf::st_sfc(snap_un, crs = rast_crs)))
          if (dist_un > lr_snap_thresh_m) {
            # Try to snap to any boundary point within threshold
            dists <- as.numeric(sf::st_distance(un, bnd_pts))
            best  <- which.min(dists)
            snap_un <- bnd_pts[[best]]
          }

          lm[1L, ] <- sf::st_coordinates(
            sf::st_sfc(snap_un, crs = rast_crs))[1L, 1:2]
          lm[nrow(lm), ] <- sf::st_coordinates(
            sf::st_sfc(snap_dn, crs = rast_crs))[1L, 1:2]

          blade <- sf::st_sfc(sf::st_linestring(lm), crs = rast_crs)

          # Split
          split_result <- suppressWarnings(lwgeom::st_split(cat_sfc, blade))
          pieces       <- suppressWarnings(
            sf::st_collection_extract(split_result, "POLYGON"))
          geoms        <- sf::st_geometry(pieces)

          if (length(geoms) < 2L) {
            logger::log_warn(
              "reconcile_divides: lr split produced <2 pieces for {fp_id}; keeping unsplit")
            return(div_r[i, , drop = FALSE])
          }

          # Build two output rows
          row_L <- div_r[i, , drop = FALSE]
          row_R <- div_r[i, , drop = FALSE]
          sf::st_geometry(row_L) <- geoms[1L]
          sf::st_geometry(row_R) <- geoms[2L]
          row_L$flowpath_id <- paste0(fp_id, "_L")
          row_R$flowpath_id <- paste0(fp_id, "_R")
          row_L$divide_id   <- row_L$flowpath_id
          row_R$divide_id   <- row_R$flowpath_id

          rbind(row_L, row_R)
        },
        error = function(e) {
          logger::log_warn(
            "reconcile_divides: lr split failed for {fp_id} ({conditionMessage(e)}); keeping unsplit")
          div_r[i, , drop = FALSE]
        })
    })

    lr_sf       <- do.call(rbind, lr_rows)
    divides_out <- suppressMessages(sf::st_transform(sf::st_as_sf(lr_sf), in_crs))
    logger::log_info("reconcile_divides: lr split complete -> {nrow(divides_out)} divide(s)")
  }

  # -- 10. Ensure output CRS matches input divides -----------------------------
  if (!is.null(divides_out) && !is.na(sf::st_crs(divides_out)) &&
    sf::st_crs(divides_out) != in_crs) {
    divides_out <- suppressMessages(sf::st_transform(divides_out, in_crs))
  }

  sf::st_as_sf(divides_out)
}

#' Allocate a split parent's displaced land to the child reach nearest each part
#'
#' @description
#' When a reference catchment is split among the fragments of a split reach, the
#' D8 partition can return sub-polygons that sit nowhere near the fragment they
#' were assigned to. The land itself is correct, it is the parent catchment and
#' no area is missing, but the assignment of land to fragment is not.
#'
#' This resolves the assignment the way the partition is itself trying to: every
#' polygon part of the displaced land goes to the fragment whose own flowline is
#' nearest it. Fragments whose partition already meets their flowline keep what
#' they hold and may gain parts, so only displaced land moves, the parent's area
#' is conserved exactly, and no two fragments can claim the same ground.
#'
#' The alternatives it replaces were both spatially arbitrary: handing the whole
#' parent to whichever fragment had the longest flowline, and handing every
#' displaced fragment the same leftover polygon.
#'
#' A fragment with no flowline cannot receive land, since there is nothing to
#' measure against; a group in which no fragment has one is left untouched.
#'
#' @param geoms List of `sfg`, the current partition polygon of each fragment.
#' @param lines List of `sfc` or `NULL`, each fragment's member flowlines.
#' @param disjoint Logical marking fragments whose polygon does not meet their
#'   own flowline.
#' @param crs CRS shared by `geoms` and `lines`.
#' @returns A list of `sfg` the same length as `geoms`, or `NULL` when there is
#'   nothing to reallocate.
#' @keywords internal
.allocate_parent_parts <- function(geoms, lines, disjoint, crs) {
  n <- length(geoms)
  if (n < 2L || !any(disjoint)) return(NULL)

  # A fragment can only receive land if it has a line to measure against.
  recip <- which(vapply(lines,
    function(x) !is.null(x) && length(x) > 0L, logical(1L)))
  if (!length(recip)) return(NULL)

  # Explode the displaced land into polygon parts. Parts are the unit of
  # reassignment: displacement moves whole sub-polygons, not fractions of them.
  # Normalise to MULTIPOLYGON first: st_cast on a mixed POLYGON/MULTIPOLYGON
  # sfc leaves the multi-part rows unexploded and warns rather than erroring.
  # Empty geometries are dropped before the cast, not after: st_cast leaves a
  # mixed sfc containing one unexploded and warns instead of erroring.
  held <- sf::st_sfc(as.list(geoms[disjoint]), crs = crs)
  held <- held[!sf::st_is_empty(held)]
  if (!length(held)) return(NULL)
  pool <- tryCatch(
    sf::st_cast(sf::st_cast(held, "MULTIPOLYGON"), "POLYGON"),
    error = function(e) NULL)
  if (is.null(pool) || !length(pool)) return(NULL)

  # One representative line per candidate recipient, so the distances between
  # every part and every candidate come from a single call.
  rep_lines <- sf::st_sfc(
    lapply(recip, function(i) sf::st_union(lines[[i]])[[1L]]), crs = crs)

  d <- suppressWarnings(sf::st_distance(pool, rep_lines))
  d <- matrix(as.numeric(d), nrow = length(pool))
  # Every part must be placeable. If any distance is unmeasurable the allocation
  # would have to drop or guess at land, so decline the group entirely and leave
  # the caller's existing geometry alone.
  if (anyNA(d) || !all(is.finite(d))) return(NULL)

  owner <- recip[max.col(-d, ties.method = "first")]

  # Nearest-part-wins on its own can strip a fragment of its catchment entirely:
  # if a sibling's line is closer to every part of the parent, a short fragment
  # receives nothing and ships with no divide. Give each fragment that gave up
  # land its own nearest part back, taking only from a fragment that can spare
  # one. Fragments only ever donate while holding more than one part, so no
  # transfer creates a new starved fragment and the loop terminates.
  gave <- intersect(which(disjoint), recip)
  col  <- match(seq_len(n), recip)
  repeat {
    starved <- setdiff(gave, unique(owner))
    if (!length(starved)) break
    moved <- FALSE
    for (s in starved) {
      cnt  <- tabulate(owner, nbins = n)
      cand <- which(cnt[owner] > 1L)
      if (!length(cand)) next
      owner[cand[which.min(d[cand, col[s]])]] <- s
      moved <- TRUE
    }
    if (!moved) break
  }

  # Displaced fragments surrender their land; every fragment keeps what it held.
  kept <- geoms
  kept[disjoint] <- list(sf::st_polygon())

  lapply(seq_len(n), function(i) {
    add <- pool[which(owner == i)]
    if (!length(add)) return(kept[[i]])
    g <- if (sf::st_is_empty(sf::st_sfc(kept[[i]], crs = crs))) add
         else c(sf::st_sfc(kept[[i]], crs = crs), add)
    tryCatch(sf::st_union(g)[[1L]], error = function(e) kept[[i]])
  })
}

#' Cut a split parent's catchment among its fragments by nearest reach
#'
#' @description
#' The partition of a split reference catchment among the fragments of the reach
#' is a nearest-channel assignment: land drains to the piece of channel closest
#' to it. The D8 trace computes that by following flow directions, which needs
#' the vector reach and the flow grid to agree on where the channel is. Where
#' they do not, Alaska routes TDX-Hydro reaches over HydroSHEDS grids, the
#' trace lands elsewhere and returns sub-polygons displaced from the fragment
#' they were assigned to, or none at all.
#'
#' This computes the same assignment geometrically. Each fragment's line is
#' sampled at `step` metres, the samples generate a nearest-point tessellation,
#' cells are dissolved by the fragment that generated them, and the result is
#' clipped to the parent. Three properties hold by construction rather than by
#' repair: the pieces tile the parent exactly, so area is conserved; no two
#' fragments can claim the same ground; and every fragment whose line meets the
#' parent receives a share, so none is left without a catchment.
#'
#' This is a partition of one real catchment among its own children, clipped to
#' that catchment. It creates no land and fills no gap between divides.
#'
#' @param parent `sfg`, the parent catchment to divide.
#' @param lines List of `sfc` or `NULL`, each fragment's member flowlines.
#' @param crs CRS shared by `parent` and `lines`.
#' @param step Spacing in CRS units at which each line is sampled. 25 m sits
#'   below the 30 m DEM cell, so the cut resolves detail the grid cannot.
#' @returns A list of `sfg` the same length as `lines`, or `NULL` if the cut
#'   cannot be computed and the caller should leave the group alone.
#' @keywords internal
.cut_parent_by_nearest_line <- function(parent, lines, crs, step = 25) {
  n <- length(lines)
  keep <- which(vapply(lines,
    function(x) !is.null(x) && length(x) > 0L, logical(1L)))
  if (!length(keep)) return(NULL)

  psfc <- sf::st_sfc(parent, crs = crs)
  if (as.numeric(sf::st_area(psfc)) <= 0) return(NULL)

  out <- rep(list(sf::st_polygon()), n)
  if (length(keep) == 1L) { out[[keep]] <- parent; return(out) }

  # Sample every fragment's line; the tessellation of the samples is the
  # nearest-line partition, to within the sample spacing.
  pts <- lapply(keep, function(i) {
    l <- tryCatch(sf::st_segmentize(sf::st_union(lines[[i]]), dfMaxLength = step),
                  error = function(e) sf::st_union(lines[[i]]))
    suppressWarnings(sf::st_cast(l, "POINT"))
  })
  owner_of_pt <- rep(keep, vapply(pts, length, integer(1L)))
  allpts <- do.call(c, pts)
  if (length(allpts) < 2L) return(NULL)

  env <- sf::st_geometry(sf::st_buffer(psfc, max(step, 100)))[[1L]]
  vor <- tryCatch(
    sf::st_collection_extract(
      sf::st_voronoi(sf::st_union(allpts), envelope = env), "POLYGON"),
    error = function(e) NULL)
  if (is.null(vor) || !length(vor)) return(NULL)
  vor <- sf::st_sfc(vor, crs = crs)

  # Each cell belongs to the fragment whose sample generated it.
  hit <- suppressWarnings(sf::st_intersects(vor, allpts))
  own <- vapply(hit, function(h) if (length(h)) owner_of_pt[h[1L]] else NA_integer_,
                integer(1L))
  if (all(is.na(own))) return(NULL)

  for (i in keep) {
    cells <- vor[which(!is.na(own) & own == i)]
    if (!length(cells)) next
    g <- tryCatch(
      sf::st_intersection(sf::st_union(cells), psfc),
      error = function(e) NULL)
    if (is.null(g) || !length(g) || sf::st_is_empty(g)[1L]) next
    g <- suppressWarnings(sf::st_collection_extract(g, "POLYGON"))
    if (!length(g)) next
    out[[i]] <- sf::st_union(g)[[1L]]
  }
  out
}
