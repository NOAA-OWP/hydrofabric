# =====================================================================
# Why divides are partitioned from the grids rather than clipped:
# the package documentation.
# =============================================================================
# Point-based divide splitting via in-file FAC/FDR D8 delineation
# -----------------------------------------------------------------------------
# Lightweight replacement for reconcile_divides. Given a set of
# divide polygons and split points (one per desired sub-catchment outlet), it
# re-delineates each multi-member parent divide directly from the D8 flow-
# direction (FDR) and flow-accumulation (FAC) rasters, with no dependency on the
# heavier raster stack. Called by partition_reconcile_divides (Stage 1).
#
# Pipeline, per parent divide (see the helpers below):
#   1. .snap_outlet_cell      : snap each split point to a high-FAC cell
#   2. .label_parent_d8       : D8-walk upstream from each outlet, labelling
#                                every raster cell with the outlet it drains to
#   3. .labels_to_pieces      : polygonise the labelled cells into sub-divides
#   4. .assign_residual_to_neighbor: fold unlabelled residual area into the
#                                best-touching neighbour piece (area-conserving)
# Falls back to a geometric (non-raster) split when delineation fails and
# `fallback_geometric = TRUE`.
#
# NOTE: ids on `divides`/`split_points` (divide_id, flowpath_id) are Lynker
# Spatial reference-fabric ids, not raw NHD COMIDs.
# =============================================================================


# Does a split piece carry the reach it is being cut for? This is the same
# question the rest of the fabric asks of any flowpath/divide pairing,
# intersects, else within 30 m, as used by `.clear_disjoint_fp_divides` and the
# `divide_has_member_flowline` gate. Asking it here rather than downstream is the
# point: a piece that fails is not a divide for this child, and stage 3 deletes
# the pairing anyway.
#
# TRUE when no reach geometry is available, so callers that supply no
# `reach_lines` keep the previous behaviour exactly.
.piece_holds_reach <- function(piece_geom, reach_line) {
  if (is.null(reach_line) || is.null(piece_geom)) return(TRUE)
  rl <- tryCatch(sf::st_transform(reach_line, sf::st_crs(piece_geom)),
                 error = function(e) NULL)
  if (is.null(rl)) return(TRUE)
  d <- tryCatch(suppressWarnings(as.numeric(sf::st_distance(
         sf::st_union(piece_geom), sf::st_union(rl)))),
       error = function(e) NA_real_)
  length(d) == 1L && is.finite(d) && d <= 30
}

.piece_reach_coverage <- function(piece_geom, reach_line) {
  if (is.null(reach_line) || is.null(piece_geom)) return(NA_real_)
  rl <- tryCatch(sf::st_transform(reach_line, sf::st_crs(piece_geom)),
                 error = function(e) NULL)
  if (is.null(rl)) return(NA_real_)
  total <- tryCatch(sum(as.numeric(sf::st_length(rl)), na.rm = TRUE),
                    error = function(e) NA_real_)
  if (!is.finite(total) || total <= 0) return(NA_real_)
  inside <- tryCatch(
    {
      inter <- suppressWarnings(sf::st_intersection(rl, sf::st_union(piece_geom)))
      if (length(inter) == 0L) 0 else sum(as.numeric(sf::st_length(inter)), na.rm = TRUE)
    },
    error = function(e) NA_real_)
  if (!is.finite(inside)) return(NA_real_)
  inside / total
}

.line_tangent_near_point <- function(reach_line, point_geom) {
  if (is.null(reach_line) || is.null(point_geom)) return(NULL)
  line <- tryCatch(sf::st_transform(reach_line, sf::st_crs(point_geom)),
                   error = function(e) NULL)
  if (is.null(line)) return(NULL)
  coords <- tryCatch(sf::st_coordinates(line), error = function(e) NULL)
  pxy <- tryCatch(as.numeric(sf::st_coordinates(point_geom)[1, c("X", "Y")]),
                  error = function(e) NULL)
  if (is.null(coords) || is.null(pxy) || nrow(coords) < 2L) return(NULL)
  d <- sqrt((coords[, "X"] - pxy[1])^2 + (coords[, "Y"] - pxy[2])^2)
  j <- which.min(d)
  if (!length(j) || !is.finite(d[j])) return(NULL)
  if (j == 1L) {
    a <- coords[j, c("X", "Y")]
    b <- coords[j + 1L, c("X", "Y")]
  } else if (j == nrow(coords)) {
    a <- coords[j - 1L, c("X", "Y")]
    b <- coords[j, c("X", "Y")]
  } else {
    a <- coords[j - 1L, c("X", "Y")]
    b <- coords[j + 1L, c("X", "Y")]
  }
  tangent <- as.numeric(b - a)
  n <- sqrt(sum(tangent^2))
  if (!is.finite(n) || n <= 0) return(NULL)
  tangent / n
}

# Assign two split parts to the two children by the reach each one carries. The
# D8 path builds (intersection, difference) in child order, so row position is
# meaningful there; `lwgeom::st_split` returns parts in an order it does not
# define, and the selection downstream reads row position. Ordering by reach
# makes the geometric path correct on its own terms.
#
# Returns c(upstream, downstream) row indices, or NULL when the parts cannot be
# told apart or either child's reach is unheld, which is itself a failed split,
# not something to paper over with an arbitrary order.
# Order a parent's members from most upstream to most downstream by walking the
# within-group chain: each member names its successor in `dsids`, and the head is
# the one no other member points at.
#
# The splitter used to sort these by id and describe the result as "upstream-most
# by point order". The ids are character, so that sort was lexicographic,
# "10" before "2", and carried no topological meaning at all. It matters
# because when two outlets snap to the same raster cell the FIRST is kept, so the
# sort decides which sibling keeps the ground.
#
# Returns a permutation of `seq_along(ids)`, or NULL when the members do not form
# a simple chain (a fork, a cycle, or a break), which is the caller's signal that
# no topological order exists and it should not claim one.
.topo_order_members <- function(ids, dsids) {
  ids <- as.character(ids); dsids <- as.character(dsids)
  n <- length(ids)
  if (n == 0L) return(integer(0))
  if (anyDuplicated(ids)) return(NULL)
  succ <- match(dsids, ids)                    # NA = leaves the group (the sink)
  # a fork (two members with the same successor) has no single chain
  inner <- succ[!is.na(succ)]
  if (anyDuplicated(inner)) return(NULL)
  head <- setdiff(seq_len(n), inner)           # pointed at by nobody
  if (length(head) != 1L) return(NULL)
  out <- integer(n)
  cur <- head
  for (k in seq_len(n)) {
    out[k] <- cur
    nxt <- succ[cur]
    if (is.na(nxt)) {
      if (k != n) return(NULL)                 # chain ended early: a break
      break
    }
    if (nxt %in% out[seq_len(k)]) return(NULL) # cycle
    cur <- nxt
  }
  if (anyDuplicated(out) || any(out == 0L)) return(NULL)
  out
}

.order_parts_by_reach <- function(parts, up_line, dn_line) {
  if (is.null(up_line) || is.null(dn_line) || nrow(parts) < 2L) return(NULL)
  g <- sf::st_geometry(parts)
  cov_up <- vapply(seq_along(g), function(i) {
    .piece_reach_coverage(g[i], up_line)
  }, numeric(1))
  cov_dn <- vapply(seq_along(g), function(i) {
    .piece_reach_coverage(g[i], dn_line)
  }, numeric(1))
  if (any(is.finite(cov_up)) && any(is.finite(cov_dn))) {
    best <- NULL
    best_score <- -Inf
    for (iu in seq_along(g)) {
      for (id in seq_along(g)) {
        if (iu == id) next
        if (!is.finite(cov_up[iu]) || !is.finite(cov_dn[id])) next
        score <- cov_up[iu] + cov_dn[id]
        if (is.finite(score) && score > best_score) {
          best_score <- score
          best <- c(iu, id)
        }
      }
    }
    if (!is.null(best) && best_score > 0) {
      return(best)
    }
  }

  dist_to <- function(line) {
    l <- tryCatch(sf::st_union(sf::st_transform(line, sf::st_crs(parts))),
                  error = function(e) NULL)
    if (is.null(l)) return(rep(NA_real_, length(g)))
    tryCatch(suppressWarnings(as.numeric(sf::st_distance(g, l))),
             error = function(e) rep(NA_real_, length(g)))
  }
  du <- dist_to(up_line); dd <- dist_to(dn_line)
  if (all(is.na(du)) || all(is.na(dd))) return(NULL)
  iu <- which.min(du); id <- which.min(dd)
  if (!length(iu) || !length(id) || iu == id) return(NULL)
  if (!(is.finite(du[iu]) && du[iu] <= 30)) return(NULL)
  if (!(is.finite(dd[id]) && dd[id] <= 30)) return(NULL)
  c(iu, id)
}

#' Split multi-member divides via point-based FAC/FDR D8 delineation
#'
#' Re-delineates each parent divide that carries more than one member flowpath
#' into per-outlet sub-divides, working directly from D8 flow-direction (FDR)
#' and flow-accumulation (FAC) rasters. This is the lightweight in-file splitter
#' used by the refactor phase's `partition_reconcile_divides`; it intentionally does not
#' route through `reconcile_divides`.
#'
#' @param divides `sf` polygons to split. Must contain `divide_id` and
#'   `flowpath_id` (Lynker Spatial reference-fabric ids, not NHD COMIDs).
#' @param split_points `sf`/data.frame of outlet points, one per desired
#'   sub-divide. Must contain `divide_id`, `flowpath_id`, `downstream_id`.
#' @param fdr_path,fac_path Paths to the D8 flow-direction and flow-accumulation
#'   rasters (any CRS; the splitter reads the CRS from the FDR raster).
#' @param fallback_fdr_path,fallback_fac_path Optional secondary rasters used
#'   when a point falls outside the primary raster footprint.
#' @param reach_lines Optional `sf` of reach linestrings carrying a `line_id`
#'   column. When supplied, and when `split_points` carries `reach_line_id` /
#'   `sib_line_id`, each split piece keeps the component that holds its own
#'   reach instead of every component above an area threshold. Omit to retain
#'   the area-based behaviour exactly.
#' @param min_area_m2 Minimum parent-divide area to attempt a raster split
#'   (smaller parents are returned unsplit).
#' @param snap_distance_m Max distance to snap a split point onto a high-FAC cell.
#' @param simplify_tolerance_m Douglas-Peucker tolerance applied to output pieces.
#' @param parallel,num_workers,parallel_backend Parallelisation controls for the
#'   per-parent delineation (furrr/future).
#' @param min_piece_area_m2 Pieces smaller than this are merged into a neighbour
#'   rather than emitted as standalone sub-divides.
#' @param fallback_geometric If `TRUE`, fall back to a geometric split when
#'   raster delineation fails for a parent.
#' @param use_legacy_reconcile Retained for signature compatibility; the in-file
#'   splitter is always used regardless of this flag.
#' @return An `sf` of divides with the same schema as `divides`, with split
#'   parents replaced by their per-outlet sub-divides (total area conserved).
#' @keywords internal
#' @noRd
#' @importFrom dplyr mutate filter bind_rows arrange select ungroup group_by desc all_of
#' @importFrom sf st_drop_geometry
split_divides_d8 <- function(divides,
                             split_points,
                             fdr_path,
                             fac_path,
                             fallback_fdr_path = NULL,
                             fallback_fac_path = NULL,
                             reach_lines = NULL,
                             min_area_m2 = 800000,
                             snap_distance_m = 100,
                             simplify_tolerance_m = 50,
                             parallel = TRUE,
                             num_workers = 4,
                             parallel_backend = c("auto", "multisession", "multicore", "sequential"),
                             min_piece_area_m2 = 50,
                             fallback_geometric = TRUE,
                             use_legacy_reconcile = TRUE) {

  if (is.null(divides) || nrow(divides) == 0) {
    .hf_log_warn("No divides provided")
    return(divides)
  }

  if (is.null(split_points) || nrow(split_points) == 0) {
    message("No split points - returning original divides")
    return(divides)
  }

  # Check required columns
  required_cols <- c("divide_id", "flowpath_id")
  if (any(!required_cols %in% names(divides))) {
    stop("divides missing required columns: divide_id, flowpath_id")
  }

  required_point_cols <- c("divide_id", "flowpath_id", "downstream_id")
  if (any(!required_point_cols %in% names(split_points))) {
    stop("split_points missing required columns: divide_id, flowpath_id, downstream_id")
  }

  # Check rasters exist
  if (!.raster_accessible(fdr_path)) stop("FDR not found: ", fdr_path)
  if (!.raster_accessible(fac_path)) stop("FAC not found: ", fac_path)

  # This implementation uses the lightweight in-file FAC/FDR approach.
  # Keep compatibility flags in signature, but do not route through reconcile_divides.
  if (isTRUE(use_legacy_reconcile)) {
    message("Using lightweight point-based FAC/FDR splitter (not reconcile_divides).")
  }

  # Determine raster CRS without retaining large raster objects in the enclosing
  # environment (they should be opened inside worker processes to avoid heavy exports).
  tmp_rast <- tryCatch(terra::rast(fdr_path), error = function(e) NULL)
  if (is.null(tmp_rast)) stop("Unable to read FDR raster for CRS: ", fdr_path)
  raster_crs <- tryCatch(sf::st_crs(terra::crs(tmp_rast)), error = function(e) NA)
  rm(tmp_rast)
  gc()
  if (is.na(raster_crs)) {
    stop("Unable to determine CRS from FDR raster")
  }

  # Normalize IDs
  divides <- divides |>
    dplyr::mutate(
      divide_id = as.character(divide_id),
      flowpath_id = as.character(flowpath_id)
    )

  split_points <- split_points |>
    dplyr::mutate(
      divide_id = as.character(divide_id),
      flowpath_id = as.character(flowpath_id),
      downstream_id = as.character(downstream_id),
      downstream_x = if ("downstream_x" %in% names(split_points)) as.numeric(downstream_x) else NA_real_,
      downstream_y = if ("downstream_y" %in% names(split_points)) as.numeric(downstream_y) else NA_real_
    )
  split_points_raster <- sf::st_transform(split_points, raster_crs)

  # Get unique parent divides to split (only those present in divides)
  divide_ids_available <- unique(as.character(divides$divide_id))
  parents_to_split <- split_points |>
    sf::st_drop_geometry() |>
    dplyr::transmute(divide_id = as.character(divide_id)) |>
    dplyr::distinct() |>
    dplyr::filter(divide_id %in% divide_ids_available) |>
    dplyr::pull(divide_id)

  if (length(parents_to_split) == 0) {
    message("No matching parent divides found for split points - returning original divides")
    return(divides)
  }

  message(sprintf("Splitting %d divides at %d points...",
    length(parents_to_split), nrow(split_points)))

  # Partition engine selector (transition flag). "legacy" = the iterative
  # per-outlet delineation + fallback cascade; "pure" = single masked
  # multi-source D8 partition per parent (non-overlap + area-conservation by
  # construction). Lives entirely inside the splitter so the call site and the
  # output contract are unchanged.
  engine <- getOption("hydrofabric.partition_engine",
    Sys.getenv("HF_PARTITION_ENGINE", "pure"))
  if (!engine %in% c("legacy", "pure")) engine <- "pure"
  message(sprintf("Divide partition engine: %s", engine))

  parallel_backend <- match.arg(parallel_backend)
  effective_workers <- max(1L, min(as.integer(num_workers), as.integer(length(parents_to_split))))
  selected_backend <- "sequential"

  # Setup parallelization if requested.
  # Uses parallel::mclapply (fork-based) on Unix/macOS for near-zero
  # serialization overhead. Workers inherit raster handles via fork CoW;
  # the lazy-init design (rasters opened inside the worker on first call)
  # ensures terra/GDAL objects are created after fork, which is safe.
  # Falls back to sequential on Windows or on any fork error.
  use_parallel <- isTRUE(parallel) && effective_workers > 1L &&
                  .Platform$OS.type != "windows" && parallel_backend != "sequential"
  selected_backend <- if (use_parallel) "multicore" else "sequential"
  if (use_parallel) {
    message(sprintf(
      "Parallel divide splitting enabled: backend=%s workers=%d parents=%d",
      selected_backend, effective_workers, length(parents_to_split)
    ))
  } else {
    message("Divide splitting executing sequentially.")
  }

  geom_col <- attr(divides, "sf_column")
  non_geom_cols <- setdiff(names(divides), geom_col)

  coerce_to_divide_schema <- function(x) {
    x <- sf::st_as_sf(x)
    missing_cols <- setdiff(non_geom_cols, names(x))
    if (length(missing_cols) > 0) {
      for (nm in missing_cols) x[[nm]] <- NA
    }

    x <- x |>
      dplyr::mutate(
        divide_id = as.character(divide_id),
        flowpath_id = as.character(flowpath_id)
      )

    x <- x[, c(non_geom_cols, geom_col), drop = FALSE]
    sf::st_as_sf(x)
  }

  select_piece_for_point <- function(piece_list, point_row) {
    pt_geom <- sf::st_geometry(point_row)

    # Build simple sfc for pieces (one geometry per piece)
    piece_geoms <- lapply(piece_list, function(psf) sf::st_geometry(psf)[[1]])
    pieces_sfc <- sf::st_sfc(piece_geoms, crs = sf::st_crs(piece_list[[1]]))

    # Prefer geos STRtree intersection + nearest queries for speed when available
    if (requireNamespace("geos", quietly = TRUE)) {
      # create tree and query intersects
      tr <- tryCatch(geos::as_geos_strtree(pieces_sfc), error = function(e) NULL)
      if (!is.null(tr)) {
        hits <- tryCatch(geos::geos_strtree_query(tr, pt_geom), error = function(e) integer(0))
        if (length(hits) > 0) {
          # geos_strtree_query returns a list (one element per query point)
          hit_vals <- if (is.list(hits)) unlist(hits) else hits
          if (length(hit_vals) > 0) return(as.integer(hit_vals[1]))
        }
        # fallback to nearest
        nearest_idx <- tryCatch(geos::geos_nearest(points = pt_geom, targets = pieces_sfc), error = function(e) NULL)
        if (!is.null(nearest_idx) && length(nearest_idx) >= 1) return(as.integer(nearest_idx[1]))
      }
    }

    # Fallback to sf-based intersects/nearest
    intersects <- vapply(piece_list, function(piece_sf) {
      hit <- tryCatch(
        sf::st_intersects(sf::st_geometry(piece_sf), pt_geom, sparse = FALSE),
        error = function(e) matrix(FALSE, nrow = 1, ncol = 1)
      )
      isTRUE(any(hit))
    }, logical(1))

    if (any(intersects)) {
      return(which(intersects)[1])
    }

    dists <- vapply(piece_list, function(piece_sf) {
      as.numeric(sf::st_distance(sf::st_geometry(piece_sf), pt_geom, by_element = FALSE))[1]
    }, numeric(1))

    which.min(dists)
  }

  delineate_contributing_area <- function(fdr_rast, fac_rast, pour_x, pour_y) {
    pour_cell <- terra::cellFromXY(fdr_rast, cbind(pour_x, pour_y))
    if (is.na(pour_cell)) return(NULL)

    watershed <- terra::rast(fdr_rast)
    fdr_matrix <- terra::as.matrix(fdr_rast, wide = TRUE)
    if (any(dim(fdr_matrix) == 0)) return(NULL)

    start_rc <- terra::rowColFromCell(fdr_rast, pour_cell)
    if (any(!is.finite(start_rc))) return(NULL)

    # D8 offsets around the current cell and the corresponding flow code
    # that would indicate the neighbor drains into the current cell.
    nbr_row_off <- c(0L, 1L, 1L, 1L, 0L, -1L, -1L, -1L)
    nbr_col_off <- c(1L, 1L, 0L, -1L, -1L, -1L, 0L, 1L)
    nbr_flow_to_current <- c(16, 32, 64, 128, 1, 2, 4, 8)
    nrow_mat <- nrow(fdr_matrix)
    ncol_mat <- ncol(fdr_matrix)

    visited <- matrix(FALSE, nrow = nrow(fdr_matrix), ncol = ncol(fdr_matrix))
    queued <- matrix(FALSE, nrow = nrow(fdr_matrix), ncol = ncol(fdr_matrix))

    q_cap <- max(1024L, min(500000L, nrow_mat * ncol_mat))
    q_row <- integer(q_cap)
    q_col <- integer(q_cap)
    q_head <- 1L
    q_tail <- 1L

    sr <- as.integer(start_rc[1])
    sc <- as.integer(start_rc[2])
    if (is.na(sr) || is.na(sc) || sr < 1L || sc < 1L || sr > nrow_mat || sc > ncol_mat) return(NULL)
    q_row[q_tail] <- sr
    q_col[q_tail] <- sc
    queued[sr, sc] <- TRUE

    enqueue_cell <- function(rr, cc) {
      q_tail <<- q_tail + 1L
      if (q_tail > length(q_row)) {
        new_len <- as.integer(length(q_row) * 2L)
        length(q_row) <<- new_len
        length(q_col) <<- new_len
      }
      q_row[q_tail] <<- rr
      q_col[q_tail] <<- cc
      queued[rr, cc] <<- TRUE
    }

    # Traverse only connected upstream cells instead of scanning full raster each iteration.
    while (q_head <= q_tail) {
      r <- q_row[q_head]
      c <- q_col[q_head]
      q_head <- q_head + 1L
      if (r < 1L || c < 1L || r > nrow_mat || c > ncol_mat || visited[r, c]) next

      visited[r, c] <- TRUE
      for (k in seq_along(nbr_row_off)) {
        rr <- r + nbr_row_off[k]
        cc <- c + nbr_col_off[k]
        if (rr < 1L || cc < 1L || rr > nrow_mat || cc > ncol_mat) next
        if (visited[rr, cc] || queued[rr, cc]) next

        v <- fdr_matrix[rr, cc]
        if (!is.na(v) && v == nbr_flow_to_current[k]) {
          enqueue_cell(rr, cc)
        }
      }
    }

    visited_rc <- which(visited, arr.ind = TRUE)
    if (nrow(visited_rc) == 0) return(NULL)

    ws_vals <- rep(NA_real_, terra::ncell(watershed))
    ws_cells <- terra::cellFromRowCol(watershed, visited_rc[, "row"], visited_rc[, "col"])
    ws_vals[ws_cells] <- 1
    terra::values(watershed) <- ws_vals
    watershed
  }

  delineate_watershed_simple <- function(point_row_raster, output_crs, fdr_rast, fac_rast,
                                         divide_poly = NULL) {
    coords <- sf::st_coordinates(point_row_raster)
    x <- as.numeric(coords[1, "X"])
    y <- as.numeric(coords[1, "Y"])
    if (!is.finite(x) || !is.finite(y)) return(NULL)

    res_x <- terra::res(fac_rast)[1]
    if (!is.finite(res_x) || res_x <= 0) return(NULL)

    # Progressive focal search: start at +/-3 cells, doubling each step until
    # valid FAC data is found. Keeps the snap near the split point so a distant
    # high-accumulation mainstem cell in the same parent divide cannot win the
    # max-FAC selection. Working in multiples of res_x keeps the window in native
    # raster units regardless of CRS.
    fac_crop <- NULL
    n_cells  <- 3L
    for (.snap_iter in seq_len(12L)) {       # 3 -> 6 -> 12 -> ... -> 3x2^11 = 6144 cells
      r_dist <- n_cells * res_x
      fe <- terra::ext(x - r_dist, x + r_dist, y - r_dist, y + r_dist)
      fc <- tryCatch(terra::crop(fac_rast, fe), error = function(e) NULL)
      if (!is.null(fc) && !all(is.na(terra::values(fc, mat = FALSE)))) {
        fac_crop <- fc
        break
      }
      n_cells <- n_cells * 2L
    }
    if (is.null(fac_crop)) return(NULL)

    # Secondary mask: restrict to cells inside divide_poly when provided. Helps
    # when piece_sf is already a small post-split polygon, so masking excludes a
    # neighbouring piece's channel cells. For a large unsplit parent the
    # progressive search above is the primary constraint and the mask is a no-op.
    fac_for_snap <- fac_crop
    if (!is.null(divide_poly) && inherits(divide_poly, c("sf", "sfc"))) {
      tryCatch(
        {
          dp_rast_crs <- sf::st_transform(
            if (inherits(divide_poly, "sf")) sf::st_geometry(divide_poly) else divide_poly,
            sf::st_crs(terra::crs(fac_crop)))
          dp_vect <- terra::vect(sf::st_as_sf(dp_rast_crs))
          masked  <- terra::mask(fac_crop, dp_vect)
          m_vals  <- terra::values(masked, mat = FALSE)
          if (!all(is.na(m_vals))) fac_for_snap <- masked
        },
        error = function(e) NULL)  # on any error keep fac_for_snap = fac_crop
    }

    fac_vals <- terra::values(fac_for_snap, mat = FALSE)
    # Some accumulation products carry a large positive NoData sentinel
    # (HydroSHEDS v2 uses 1.79e308) that survives as a value rather than NA when
    # the sentinel is not declared on every member of a mosaic. Left in place it
    # wins which.max outright and snaps every pour point onto a void cell.
    fac_vals[!is.finite(fac_vals) | fac_vals > 1e300] <- NA
    if (length(fac_vals) == 0 || all(is.na(fac_vals))) return(NULL)

    max_coords <- terra::xyFromCell(fac_for_snap, which.max(fac_vals))
    snap_xy <- c(as.numeric(max_coords[1, "x"]), as.numeric(max_coords[1, "y"]))
    snap_x <- snap_xy[1]
    snap_y <- snap_xy[2]
    if (!is.finite(snap_x) || !is.finite(snap_y)) return(NULL)

    watershed_raster <- delineate_contributing_area(fdr_rast, fac_rast, snap_x, snap_y)
    if (is.null(watershed_raster)) return(NULL)

    watershed_poly <- tryCatch(
      terra::as.polygons(watershed_raster, dissolve = TRUE) |> sf::st_as_sf(),
      error = function(e) NULL
    )
    if (is.null(watershed_poly) || nrow(watershed_poly) == 0) return(NULL)

    sf::st_transform(watershed_poly, output_crs)
  }

  split_piece_once <- function(piece_sf, point_row, point_row_raster, fdr_parent, fac_parent, local_window_cache = NULL,
                               use_fdr_fallback = FALSE, fdr_parent_fb = NULL, fac_parent_fb = NULL,
                               is_two_way_split = TRUE) {
    # Fast geometric set ops via geos when available; fall back to sf.
    .sfc_intersection <- function(a, b) {
      if (requireNamespace("geos", quietly = TRUE)) {
        tryCatch(
          {
            r <- geos::geos_intersection(geos::as_geos_geometry(a), geos::as_geos_geometry(b))
            sf::st_set_crs(sf::st_as_sfc(r), sf::st_crs(a))
          },
          error = function(e) suppressWarnings(sf::st_intersection(a, b)))
      } else {
        suppressWarnings(sf::st_intersection(a, b))
      }
    }
    .sfc_difference <- function(a, b) {
      if (requireNamespace("geos", quietly = TRUE)) {
        tryCatch(
          {
            r <- geos::geos_difference(geos::as_geos_geometry(a), geos::as_geos_geometry(b))
            sf::st_set_crs(sf::st_as_sfc(r), sf::st_crs(a))
          },
          error = function(e) suppressWarnings(sf::st_difference(a, b)))
      } else {
        suppressWarnings(sf::st_difference(a, b))
      }
    }

    geometric_split_once <- function(piece_sf_local, point_row_local) {
      if (!requireNamespace("lwgeom", quietly = TRUE)) return(NULL)

      poly <- sf::st_geometry(piece_sf_local)[[1]]
      up_pt <- sf::st_geometry(point_row_local)[[1]]
      down_x <- suppressWarnings(as.numeric(point_row_local$downstream_x[1]))
      down_y <- suppressWarnings(as.numeric(point_row_local$downstream_y[1]))
      down_pt <- if (is.finite(down_x) && is.finite(down_y)) sf::st_point(c(down_x, down_y)) else up_pt

      up_xy <- as.numeric(up_pt)
      down_xy <- as.numeric(down_pt)
      v <- .line_tangent_near_point(.up_line, sf::st_sfc(up_pt, crs = sf::st_crs(point_row_local)))
      if (is.null(v)) v <- down_xy - up_xy
      if (!all(is.finite(v)) || sqrt(sum(v^2)) < 1e-6) v <- c(1, 0)
      n <- c(-v[2], v[1])
      n_norm <- sqrt(sum(n^2))
      if (!is.finite(n_norm) || n_norm < 1e-9) n <- c(0, 1) else n <- n / n_norm

      bb <- sf::st_bbox(sf::st_sfc(poly, crs = sf::st_crs(piece_sf_local)))
      span <- max(bb$xmax - bb$xmin, bb$ymax - bb$ymin)
      if (!is.finite(span) || span <= 0) span <- 1000
      L <- span * 3

      p1 <- up_xy - n * L
      p2 <- up_xy + n * L
      cut_line <- sf::st_sfc(sf::st_linestring(rbind(p1, p2)), crs = sf::st_crs(piece_sf_local))

      split_geom <- tryCatch(
        lwgeom::st_split(sf::st_sfc(poly, crs = sf::st_crs(piece_sf_local)), cut_line),
        error = function(e) NULL
      )
      if (is.null(split_geom)) return(NULL)

      parts <- tryCatch(
        sf::st_collection_extract(split_geom, "POLYGON"),
        error = function(e) NULL
      )
      if (is.null(parts) || length(parts) < 2) return(NULL)

      out <- sf::st_sf(
        piece_sf_local[rep(1, length(parts)), non_geom_cols, drop = FALSE],
        geometry = parts
      )
      out <- coerce_to_divide_schema(out)
      out
    }

    split_result <- NULL
    piece_raster <- tryCatch(sf::st_transform(piece_sf, raster_crs), error = function(e) NULL)
    if (is.null(piece_raster) || nrow(piece_raster) == 0) return(NULL)

    piece_vect <- tryCatch(terra::vect(piece_raster), error = function(e) NULL)
    if (is.null(piece_vect)) return(NULL)
    piece_buf <- tryCatch(terra::buffer(piece_vect, width = 120), error = function(e) NULL)
    if (is.null(piece_buf)) return(NULL)

    # Cache local raster windows by extent to avoid repeated crop I/O when
    # successive split points target the same (or equivalent) piece window.
    cache_key <- NULL
    if (!is.null(local_window_cache)) {
      ex <- tryCatch(terra::ext(piece_buf), error = function(e) NULL)
      if (!is.null(ex)) {
        ex_num <- c(ex$xmin, ex$xmax, ex$ymin, ex$ymax)
        # Round key to mm precision in projected CRS to improve cache-hit stability.
        cache_key <- paste(sprintf("%.3f", ex_num), collapse = "|")
      }
    }

    if (!is.null(cache_key) && exists(cache_key, envir = local_window_cache, inherits = FALSE)) {
      cached <- get(cache_key, envir = local_window_cache, inherits = FALSE)
      fdr_local <- cached$fdr
      fac_local <- cached$fac
    } else {
      fdr_local <- tryCatch(
        terra::crop(fdr_parent, piece_buf, snap = "out", mask = TRUE),
        error = function(e) NULL
      )
      fac_local <- tryCatch(
        terra::crop(fac_parent, piece_buf, snap = "out", mask = TRUE),
        error = function(e) NULL
      )
      if (!is.null(cache_key) && !is.null(fdr_local) && !is.null(fac_local)) {
        assign(cache_key, list(fdr = fdr_local, fac = fac_local), envir = local_window_cache)
      }
    }
    if (is.null(fdr_local) || is.null(fac_local)) return(NULL)

    # Own and sibling reach geometry, resolved by the ids the caller attached to
    # the split point. `reach_lines` comes from the enclosing split_divides_d8
    # frame. Absent for callers that supply none, in which case collapse_piece
    # keeps its area-based component rules unchanged.
    .lookup_line <- function(col) {
      if (is.null(reach_lines) || !(col %in% names(point_row))) return(NULL)
      if (!("line_id" %in% names(reach_lines))) return(NULL)
      id <- as.character(point_row[[col]][1L])
      if (is.na(id) || !nzchar(id)) return(NULL)
      k <- which(as.character(reach_lines$line_id) == id)
      if (!length(k)) return(NULL)
      sf::st_geometry(reach_lines)[k]
    }
    .up_line <- .lookup_line("reach_line_id")
    .dn_line <- .lookup_line("sib_line_id")

    watershed_poly <- tryCatch(
      delineate_watershed_simple(point_row_raster, sf::st_crs(piece_sf), fdr_local, fac_local,
        divide_poly = piece_sf),
      error = function(e) NULL
    )


    # Fallback: retry delineation when the primary FDR returns no valid
    # cells OR captures too little of the piece polygon. The null/empty case
    # covers out-of-domain points (Canadian border, Great Lakes). The
    # low-coverage case covers flat terrain where D8 routing traces only a thin
    # ribbon along the channel instead of the full catchment. If less than 10%
    # of the piece area is captured the primary result is unreliable; keep
    # whichever of primary vs fallback delineates the larger area.
    primary_area_m2 <- if (!is.null(watershed_poly) && nrow(watershed_poly) > 0L)
      tryCatch(as.numeric(sf::st_area(watershed_poly)), error = function(e) 0)
    else 0
    piece_area_m2 <- tryCatch(as.numeric(sf::st_area(piece_sf)), error = function(e) 0)
    low_coverage   <- piece_area_m2 > 0 && (primary_area_m2 / piece_area_m2) < 0.10

    if (low_coverage && !use_fdr_fallback)
      message(sprintf("fallback skipped (not configured): fp %s coverage %.1f%%",
        as.character(point_row$flowpath_id[1L]),
        100 * primary_area_m2 / piece_area_m2))

    if ((is.null(watershed_poly) || nrow(watershed_poly) == 0L || low_coverage) &&
      use_fdr_fallback) {
      pt_fb <- tryCatch(
        sf::st_transform(sf::st_sfc(sf::st_geometry(point_row_raster)[[1L]],
          crs = sf::st_crs(point_row_raster)),
        sf::st_crs(fdr_parent_fb)),
        error = function(e) NULL)
      if (!is.null(pt_fb)) {
        pt_row_fb <- point_row_raster
        sf::st_geometry(pt_row_fb) <- pt_fb
        fb_poly <- tryCatch(
          delineate_watershed_simple(pt_row_fb, sf::st_crs(piece_sf),
            fdr_parent_fb, fac_parent_fb,
            divide_poly = piece_sf),
          error = function(e) NULL)
        if (!is.null(fb_poly) && nrow(fb_poly) > 0L) {
          fb_area_m2 <- tryCatch(as.numeric(sf::st_area(fb_poly)),
            error = function(e) 0)
          if (fb_area_m2 > primary_area_m2) {
            trigger <- if (low_coverage) "low-coverage" else "null-primary"
            message(sprintf(
              "fallback [%s]: %.3f -> %.3f sqkm (fp %s)",
              trigger, primary_area_m2 / 1e6, fb_area_m2 / 1e6,
              as.character(point_row$flowpath_id[1L])))
            watershed_poly <- fb_poly
          } else {
            message(sprintf(
              "fallback tried but no improvement: primary %.3f sqkm, fallback %.3f sqkm (fp %s)",
              primary_area_m2 / 1e6, fb_area_m2 / 1e6,
              as.character(point_row$flowpath_id[1L])))
          }
        } else {
          message(sprintf(
            "fallback delineation returned null/empty (fp %s)",
            as.character(point_row$flowpath_id[1L])))
        }
      }
    }

    if (!is.null(watershed_poly) && nrow(watershed_poly) > 0) {
      collapse_piece <- function(x, template_row, reach_line = NULL) {
        if (is.null(x) || nrow(x) == 0) return(NULL)

        x_valid <- tryCatch(sf::st_is_valid(x), error = function(e) rep(TRUE, nrow(x)))
        if (any(!x_valid, na.rm = TRUE)) {
          x <- tryCatch(sf::st_make_valid(x), error = function(e) x)
        }

        # Always explode to polygon components before filtering by area so
        # tiny MULTIPOLYGON fragments don't survive inside a large parent row.
        x_geom <- sf::st_geometry(x)
        x_types <- unique(as.character(sf::st_geometry_type(x_geom, by_geometry = TRUE)))
        poly_parts <- if (all(x_types %in% c("POLYGON", "MULTIPOLYGON"))) {
          x_geom
        } else {
          tryCatch(
            suppressWarnings(sf::st_collection_extract(x_geom, "POLYGON")),
            error = function(e) NULL
          )
        }
        if (is.null(poly_parts) || length(poly_parts) == 0) return(NULL)
        poly_parts <- tryCatch(
          suppressWarnings(sf::st_cast(poly_parts, "POLYGON")),
          error = function(e) poly_parts
        )
        if (length(poly_parts) == 0) return(NULL)

        x <- sf::st_sf(geometry = poly_parts)
        x$.__area <- suppressWarnings(as.numeric(sf::st_area(x)))
        x <- x[is.finite(x$.__area) & x$.__area > min_piece_area_m2, , drop = FALSE]
        if (nrow(x) == 0) return(NULL)

        # Reach-aware component selection. The area rules below keep every
        # component above a relative threshold and union them, which can leave one
        # divide in several detached pieces, the shape construct_agg_divides is
        # known to drop. Where the owning reach's geometry is known the right
        # component is unambiguous: the one that carries that reach. Intersecting
        # a raster watershed with a vector parent routinely leaves dozens of
        # stray cell clusters (median 38 per case, max 200 when measured on AK),
        # and none of them belong to this child. Orphans are not lost: they fall
        # to the sibling via reconcile_split_coverage, which enforces exact
        # coverage of the parent, so area is still conserved.
        #
        # 30 m is the fabric's universal tolerance, the DEM resolution and the
        # snap distance, and is how proximity is judged everywhere else.
        if (!is.null(reach_line) && nrow(x) > 1L) {
          rl <- tryCatch(sf::st_transform(reach_line, sf::st_crs(x)),
                         error = function(e) NULL)
          if (!is.null(rl)) {
            dd <- suppressWarnings(as.numeric(
              sf::st_distance(sf::st_geometry(x), sf::st_union(rl))))
            on_reach <- which(is.finite(dd) & dd <= 30)
            if (length(on_reach))
              x <- x[on_reach[which.max(x$.__area[on_reach])], , drop = FALSE]
          }
        }

        # DEM-derived delineations can include tiny disconnected artifact components.
        # Keep only meaningful components before union.
        component_min_m2 <- max(
          as.numeric(min_piece_area_m2),
          as.numeric(min_area_m2) * 0.005
        )
        largest_area <- suppressWarnings(max(x$.__area, na.rm = TRUE))
        rel_min_m2 <- if (is.finite(largest_area) && largest_area > 0) largest_area * 0.02 else 0
        keep_min_m2 <- max(component_min_m2, rel_min_m2)
        keep_idx <- which(is.finite(x$.__area) & x$.__area >= keep_min_m2)
        if (length(keep_idx) == 0) keep_idx <- which.max(x$.__area)
        x <- x[keep_idx, , drop = FALSE]

        geom <- sf::st_union(sf::st_geometry(x))
        geom_sfc <- sf::st_sfc(geom, crs = sf::st_crs(template_row))

        # Re-check components after union/cleanup to drop tiny slivers.
        geom_types <- unique(as.character(sf::st_geometry_type(geom_sfc, by_geometry = TRUE)))
        post_parts <- if (all(geom_types %in% c("POLYGON", "MULTIPOLYGON"))) {
          geom_sfc
        } else {
          tryCatch(
            suppressWarnings(sf::st_collection_extract(geom_sfc, "POLYGON")),
            error = function(e) NULL
          )
        }
        if (!is.null(post_parts) && length(post_parts) > 0) {
          post_parts <- tryCatch(
            suppressWarnings(sf::st_cast(post_parts, "POLYGON")),
            error = function(e) post_parts
          )
          post_areas <- suppressWarnings(as.numeric(sf::st_area(post_parts)))
          keep_post <- which(is.finite(post_areas) & post_areas >= keep_min_m2)
          if (length(keep_post) == 0 && length(post_areas) > 0) {
            keep_post <- which.max(post_areas)
          }
          if (length(keep_post) > 0) {
            geom_sfc <- sf::st_sfc(sf::st_union(post_parts[keep_post]), crs = sf::st_crs(template_row))
          }
        }

        out <- template_row
        sf::st_geometry(out) <- geom_sfc
        out
      }

      upstream_raw <- tryCatch(
        sf::st_sf(geometry = .sfc_intersection(sf::st_geometry(piece_sf), sf::st_geometry(watershed_poly))),
        error = function(e) NULL
      )

      downstream_raw <- tryCatch(
        sf::st_sf(geometry = .sfc_difference(sf::st_geometry(piece_sf), sf::st_geometry(watershed_poly))),
        error = function(e) NULL
      )

      upstream_divide <- collapse_piece(upstream_raw, piece_sf[1, , drop = FALSE],
                                        reach_line = .up_line)
      downstream_divide <- collapse_piece(downstream_raw, piece_sf[1, , drop = FALSE],
                                          reach_line = .dn_line)

      if (!is.null(upstream_divide) && !is.null(downstream_divide)) {
        split_result <- dplyr::bind_rows(upstream_divide, downstream_divide)
      }
    }

    if (!is.null(split_result) && nrow(split_result) >= 2L &&
      (!is.null(.up_line) || !is.null(.dn_line))) {
      cov_min <- as.numeric(getOption("hydrofabric.split_min_reach_coverage", 0.9))
      if (!is.finite(cov_min) || cov_min < 0) cov_min <- 0.9
      cov_up <- .piece_reach_coverage(sf::st_geometry(split_result)[1], .up_line)
      cov_dn <- .piece_reach_coverage(sf::st_geometry(split_result)[2], .dn_line)
      bad_up <- is.finite(cov_up) && cov_up < cov_min
      bad_dn <- is.finite(cov_dn) && cov_dn < cov_min
      if (bad_up || bad_dn) {
        message(sprintf(
          "split_divides_d8: raster split rejected by reach coverage (up=%.3f, down=%.3f); trying geometric fallback",
          cov_up, cov_dn))
        split_result <- NULL
      }
    }

    used_geometric_fallback <- FALSE
    if ((is.null(split_result) || nrow(split_result) < 2) && isTRUE(fallback_geometric)) {
      split_result <- geometric_split_once(piece_sf, point_row)
      used_geometric_fallback <- !is.null(split_result) && nrow(split_result) >= 2L
      # `geometric_split_once` returns parts in `lwgeom::st_split` order, which it
      # does not define, while the selection below reads row position. Order by
      # the reach each part carries when that is knowable.
      if (!is.null(split_result) && nrow(split_result) >= 2L) {
        ord <- .order_parts_by_reach(split_result, .up_line, .dn_line)
        if (!is.null(ord)) split_result <- split_result[ord, , drop = FALSE]
      }
    }

    # The masked multi-source D8 partition is a total function over the parent
    # (every cell labeled by its first downstream sub-outlet), so the raster
    # split is correct and oriented by construction. No geometric-bisector
    # orientation correction is needed here; the genuine <2-piece geometric
    # fallback above (fallback_geometric) is retained.

    if (is.null(split_result) || nrow(split_result) < 2) {
      return(NULL)
    }

    split_result <- coerce_to_divide_schema(split_result)
    areas <- suppressWarnings(as.numeric(sf::st_area(split_result)))

    finite_areas <- areas[is.finite(areas)]
    if (length(finite_areas) < 2 || sum(finite_areas > min_piece_area_m2) < 2) {
      return(NULL)
    }

    total_area <- sum(finite_areas)
    if (!is.finite(total_area) || total_area <= 0) {
      return(NULL)
    }
    if ((max(finite_areas) / total_area) > 0.9999) {
      return(NULL)
    }

    if (!any(is.finite(areas))) {
      upstream_idx <- 1L
      downstream_idx <- 2L
    } else {
      # Row 1 = watershed intersection = upstream contributing area (LARGE).
      # Row 2 = residual difference    = downstream sliver area (SMALL).
      # Order is deterministic from the bind_rows(upstream_divide, downstream_divide)
      # assembly above; select by row position, not by area (area-based selection
      # can invert and assign a parent divide to the wrong child flowpath).
      upstream_idx <- 1L
      downstream_idx <- 2L
      if (upstream_idx == downstream_idx && nrow(split_result) >= 2) {
        downstream_idx <- setdiff(seq_len(nrow(split_result)), upstream_idx)[1]
      }
    }

    # Enforce complete topological coverage: the two retained split pieces must
    # be disjoint and exactly cover the parent divide footprint.
    reconcile_split_coverage <- function(parent_piece, split_sf, up_idx, down_idx, min_area_keep) {
      crs <- sf::st_crs(parent_piece)

      to_sfc <- function(g) {
        if (inherits(g, "sfc")) return(g)
        sf::st_sfc(g, crs = crs)
      }

      normalize_poly <- function(g) {
        if (is.null(g)) return(NULL)
        g <- to_sfc(g)
        if (length(g) == 0) return(NULL)

        g_types <- unique(as.character(sf::st_geometry_type(g, by_geometry = TRUE)))
        g <- if (all(g_types %in% c("POLYGON", "MULTIPOLYGON"))) {
          g
        } else {
          tryCatch(
            suppressWarnings(sf::st_collection_extract(g, "POLYGON")),
            error = function(e) NULL
          )
        }
        if (is.null(g) || length(g) == 0) return(NULL)

        g_valid <- tryCatch(sf::st_is_valid(g), error = function(e) rep(TRUE, length(g)))
        if (any(!g_valid, na.rm = TRUE)) {
          g <- tryCatch(sf::st_make_valid(g), error = function(e) g)
        }

        g <- tryCatch(suppressWarnings(sf::st_collection_extract(g, "POLYGON")), error = function(e) g)
        g <- tryCatch(suppressWarnings(sf::st_cast(g, "POLYGON")), error = function(e) g)
        if (length(g) == 0) return(NULL)

        ga <- suppressWarnings(as.numeric(sf::st_area(g)))
        keep <- which(is.finite(ga) & ga > 0)
        if (length(keep) == 0) return(NULL)
        if (length(keep) == 1L) return(g[keep])
        sf::st_union(g[keep])
      }

      clip_to_parent <- function(g, parent_geom) {
        g <- normalize_poly(g)
        if (is.null(g)) return(NULL)
        out <- tryCatch(
          .sfc_intersection(g, parent_geom),
          error = function(e) g
        )
        normalize_poly(out)
      }

      parent_geom <- normalize_poly(sf::st_geometry(parent_piece))
      if (is.null(parent_geom)) return(NULL)

      up_geom <- clip_to_parent(sf::st_geometry(split_sf)[[up_idx]], parent_geom)
      down_geom <- clip_to_parent(sf::st_geometry(split_sf)[[down_idx]], parent_geom)
      if (is.null(up_geom) || is.null(down_geom)) return(NULL)

      # Remove overlap from downstream side first to keep a clean partition.
      down_geom <- normalize_poly(tryCatch(
        .sfc_difference(down_geom, up_geom),
        error = function(e) down_geom
      ))
      if (is.null(down_geom)) return(NULL)

      coverage <- normalize_poly(tryCatch(
        sf::st_union(up_geom, down_geom),
        error = function(e) NULL
      ))
      if (is.null(coverage)) return(NULL)

      gap_geom <- normalize_poly(tryCatch(
        .sfc_difference(parent_geom, coverage),
        error = function(e) NULL
      ))
      if (!is.null(gap_geom)) {
        gap_area <- suppressWarnings(as.numeric(sf::st_area(gap_geom))[1])
        if (is.finite(gap_area) && gap_area > 0) {
          # Assign uncovered residual to downstream piece to guarantee full coverage.
          down_geom <- normalize_poly(tryCatch(
            sf::st_union(down_geom, gap_geom),
            error = function(e) down_geom
          ))
        }
      }

      # Final clip + disjoint pass after residual assignment.
      up_geom <- clip_to_parent(up_geom, parent_geom)
      down_geom <- clip_to_parent(down_geom, parent_geom)
      down_geom <- normalize_poly(tryCatch(
        .sfc_difference(down_geom, up_geom),
        error = function(e) down_geom
      ))

      if (is.null(up_geom) || is.null(down_geom)) return(NULL)
      up_area <- suppressWarnings(as.numeric(sf::st_area(up_geom))[1])
      down_area <- suppressWarnings(as.numeric(sf::st_area(down_geom))[1])
      if (!is.finite(up_area) || !is.finite(down_area) ||
        up_area <= min_area_keep || down_area <= min_area_keep) {
        return(NULL)
      }

      out <- split_sf[c(up_idx, down_idx), , drop = FALSE]
      sf::st_geometry(out) <- sf::st_sfc(
        sf::st_geometry(up_geom)[[1]],
        sf::st_geometry(down_geom)[[1]],
        crs = crs
      )
      list(split_sf = out, upstream_idx = 1L, downstream_idx = 2L)
    }

    needs_reconcile <- TRUE
    if (nrow(split_result) == 2) {
      parent_area <- suppressWarnings(as.numeric(sf::st_area(piece_sf))[1])
      split_areas <- suppressWarnings(as.numeric(sf::st_area(split_result)))
      if (is.finite(parent_area) &&
        parent_area > 0 &&
        length(split_areas) == 2 &&
        all(is.finite(split_areas)) &&
        all(split_areas > min_piece_area_m2)) {
        rel_area_diff <- abs(sum(split_areas) - parent_area) / parent_area
        has_overlap <- tryCatch(
          {
            ov <- sf::st_overlaps(
              sf::st_geometry(split_result)[upstream_idx],
              sf::st_geometry(split_result)[downstream_idx],
              sparse = FALSE
            )
            isTRUE(ov[1, 1])
          },
          error = function(e) TRUE)
        if (is.finite(rel_area_diff) && rel_area_diff <= 0.002 && !isTRUE(has_overlap)) {
          needs_reconcile <- FALSE
        }
      }
    }

    if (isTRUE(needs_reconcile)) {
      reconciled <- reconcile_split_coverage(
        parent_piece = piece_sf,
        split_sf = split_result,
        up_idx = upstream_idx,
        down_idx = downstream_idx,
        min_area_keep = min_piece_area_m2
      )
      if (!is.null(reconciled)) {
        split_result <- reconciled$split_sf
        upstream_idx <- reconciled$upstream_idx
        downstream_idx <- reconciled$downstream_idx
      }
    }

    # Normalize multipart outputs:
    # 1) explode all multipolygons into polygon parts
    # 2) retain largest part of each feature as its primary geometry
    # 3) dissolve remaining parts into the primary feature sharing the longest border
    redistribute_multipart_fragments <- function(split_sf) {
      if (is.null(split_sf) || nrow(split_sf) == 0) return(split_sf)
      split_types <- unique(as.character(sf::st_geometry_type(split_sf, by_geometry = TRUE)))
      # Fast path: no multipart fragments to redistribute.
      if (all(split_types == "POLYGON")) return(split_sf)

      crs <- sf::st_crs(split_sf)
      base <- split_sf |>
        dplyr::mutate(.feature_id = seq_len(nrow(split_sf)))

      parts <- tryCatch(
        suppressWarnings(sf::st_cast(sf::st_cast(base, "MULTIPOLYGON"), "POLYGON")),
        error = function(e) base
      )
      if (is.null(parts) || nrow(parts) == 0) return(split_sf)

      parts$.part_id <- seq_len(nrow(parts))
      parts$.area    <- suppressWarnings(as.numeric(sf::st_area(parts)))
      parts <- parts[is.finite(parts$.area) & parts$.area > 0, , drop = FALSE]
      if (nrow(parts) == 0) return(split_sf)

      primary_parts <- parts |>
        dplyr::group_by(.feature_id) |>
        dplyr::arrange(dplyr::desc(.area), .by_group = TRUE) |>
        dplyr::slice_head(n = 1) |>
        dplyr::ungroup() |>
        dplyr::mutate(.primary_id = .feature_id)

      if (nrow(primary_parts) != nrow(split_sf)) {
        return(split_sf)
      }

      fragments <- parts |>
        dplyr::filter(!.part_id %in% primary_parts$.part_id)
      if (nrow(fragments) == 0) {
        out <- primary_parts |>
          dplyr::arrange(.feature_id)
        sf::st_geometry(out) <- sf::st_sfc(sf::st_geometry(out), crs = crs)
        out <- out |> dplyr::select(dplyr::all_of(names(split_sf)))
        return(sf::st_as_sf(out))
      }

      primary_geom <- sf::st_geometry(primary_parts)
      fragment_geom <- sf::st_geometry(fragments)

      assign_idx <- as.integer(fragments$.feature_id)
      frag_empty <- tryCatch(sf::st_is_empty(fragment_geom), error = function(e) rep(FALSE, length(fragment_geom)))
      non_empty_idx <- which(!frag_empty)
      if (length(non_empty_idx) > 0) {
        if (requireNamespace("geos", quietly = TRUE)) {
          nearest <- tryCatch(geos::geos_nearest(fragment_geom[non_empty_idx], primary_geom), error = function(e) rep(NA_integer_, length(non_empty_idx)))
          good <- is.finite(nearest) & nearest >= 1L & nearest <= nrow(primary_parts)
          if (any(good)) {
            assign_idx[non_empty_idx[good]] <- as.integer(primary_parts$.feature_id[nearest[good]])
          }
        } else {
          nearest <- tryCatch(
            sf::st_nearest_feature(fragment_geom[non_empty_idx], primary_geom),
            error = function(e) rep(NA_integer_, length(non_empty_idx))
          )
          good <- is.finite(nearest) & nearest >= 1L & nearest <= nrow(primary_parts)
          if (any(good)) {
            assign_idx[non_empty_idx[good]] <- as.integer(primary_parts$.feature_id[nearest[good]])
          }
        }
      }

      fragments$.target_feature_id <- assign_idx

      out_geoms <- vector("list", nrow(primary_parts))
      for (k in seq_len(nrow(primary_parts))) {
        fid <- primary_parts$.feature_id[k]
        geom_list <- list(primary_geom[k][[1]])
        add_idx <- which(fragments$.target_feature_id == fid)
        if (length(add_idx) == 0) {
          out_geoms[[k]] <- primary_geom[k][[1]]
          next
        }
        geom_list <- c(geom_list, lapply(add_idx, function(ii) fragment_geom[ii][[1]]))
        geom_union <- tryCatch(
          sf::st_union(sf::st_sfc(geom_list, crs = crs)),
          error = function(e) sf::st_sfc(primary_geom[k][[1]], crs = crs)
        )
        out_geoms[[k]] <- sf::st_geometry(geom_union)[[1]]
      }

      out <- primary_parts |>
        dplyr::arrange(.feature_id)
      sf::st_geometry(out) <- sf::st_sfc(out_geoms, crs = crs)
      out <- out |> dplyr::select(dplyr::all_of(names(split_sf)))
      sf::st_as_sf(out)
    }
    split_result <- redistribute_multipart_fragments(split_result)

    up_id <- as.character(point_row$flowpath_id[1])
    down_id <- as.character(point_row$downstream_id[1])

    # Validate the upstream/downstream assignment against known geometry. The
    # FDR-traced watershed can come out inverted (Whitebox-vs-Esri FDR encoding
    # mismatch, or the max-FAC pour-point snap biasing downstream), leaving
    # split_result rows with the wrong orientation. Use the downstream
    # flowline's outlet coordinate, which must lie inside the downstream-side
    # catchment, as arbiter: if it is closer to the piece currently tagged
    # upstream, flip the indices.
    if (!isTRUE(used_geometric_fallback) && !is.null(point_row) && nrow(point_row) >= 1L) {
      try(
        {
          dn_coord <- c(as.numeric(point_row$downstream_x[1]),
            as.numeric(point_row$downstream_y[1]))
          if (all(is.finite(dn_coord))) {
            dn_pt <- sf::st_transform(
              sf::st_sfc(sf::st_point(dn_coord), crs = sf::st_crs(point_row)),
              sf::st_crs(split_result))
            d_up <- as.numeric(sf::st_distance(
              dn_pt, sf::st_geometry(split_result)[upstream_idx]))
            d_dn <- as.numeric(sf::st_distance(
              dn_pt, sf::st_geometry(split_result)[downstream_idx]))
            if (is.finite(d_up) && is.finite(d_dn) && d_up < d_dn) {
              tmp <- upstream_idx
              upstream_idx <- downstream_idx
              downstream_idx <- tmp
            }
          }
        },
        silent = TRUE)
    }

    split_result$divide_id <- as.character(split_result$divide_id)
    split_result$flowpath_id <- as.character(split_result$flowpath_id)
    split_result$divide_id[upstream_idx] <- up_id
    split_result$flowpath_id[upstream_idx] <- up_id
    split_result$divide_id[downstream_idx] <- down_id
    split_result$flowpath_id[downstream_idx] <- down_id

    other_idx <- setdiff(seq_len(nrow(split_result)), c(upstream_idx, downstream_idx))
    if (length(other_idx) > 0) {
      split_result$divide_id[other_idx] <- down_id
      split_result$flowpath_id[other_idx] <- down_id
    }

    split_result
  }

  process_parent <- function(parent_divide, parent_points, parent_points_raster, parent_id, fdr_base, fac_base,
                             fdr_fallback_base = NULL, fac_fallback_base = NULL) {
    parent_divide <- coerce_to_divide_schema(parent_divide)
    key <- paste(parent_points$flowpath_id, parent_points$downstream_id, sep = "|")
    keep <- !duplicated(key)
    parent_points <- parent_points[keep, ]
    parent_points_raster <- parent_points_raster[keep, ]
    ord <- order(parent_points$flowpath_id, parent_points$downstream_id)
    parent_points <- parent_points[ord, ]
    parent_points_raster <- parent_points_raster[ord, ]

    attempted <- nrow(parent_points)
    successful <- 0L
    pieces <- list(parent_divide)

    parent_raster <- tryCatch(sf::st_transform(parent_divide, raster_crs), error = function(e) NULL)
    if (is.null(parent_raster) || nrow(parent_raster) == 0) {
      return(list(parent_id = parent_id, attempted = attempted, successful = 0L, pieces = parent_divide))
    }
    parent_vect <- tryCatch(terra::vect(parent_raster), error = function(e) NULL)
    if (is.null(parent_vect)) {
      return(list(parent_id = parent_id, attempted = attempted, successful = 0L, pieces = parent_divide))
    }
    parent_buf <- tryCatch(terra::buffer(parent_vect, width = 250), error = function(e) NULL)
    if (is.null(parent_buf)) {
      return(list(parent_id = parent_id, attempted = attempted, successful = 0L, pieces = parent_divide))
    }
    fdr_parent <- tryCatch(
      terra::crop(fdr_base, parent_buf, snap = "out", mask = TRUE),
      error = function(e) NULL
    )
    fac_parent <- tryCatch(
      terra::crop(fac_base, parent_buf, snap = "out", mask = TRUE),
      error = function(e) NULL
    )
    if (is.null(fdr_parent) || is.null(fac_parent)) {
      return(list(parent_id = parent_id, attempted = attempted, successful = 0L, pieces = parent_divide))
    }

    # Crop fallback rasters to parent extent if provided (used when primary yields no watershed).
    use_fdr_fallback <- !is.null(fdr_fallback_base) && !is.null(fac_fallback_base)
    fdr_parent_fb <- if (use_fdr_fallback) tryCatch(
      terra::crop(fdr_fallback_base, terra::project(parent_buf, terra::crs(fdr_fallback_base)),
        snap = "out", mask = TRUE),
      error = function(e) NULL) else NULL
    fac_parent_fb <- if (use_fdr_fallback) tryCatch(
      terra::crop(fac_fallback_base, terra::project(parent_buf, terra::crs(fac_fallback_base)),
        snap = "out", mask = TRUE),
      error = function(e) NULL) else NULL
    use_fdr_fallback <- use_fdr_fallback && !is.null(fdr_parent_fb) && !is.null(fac_parent_fb)

    # parent-local cache for cropped raster windows used by split_piece_once
    local_window_cache <- new.env(parent = emptyenv())

    for (i in seq_len(nrow(parent_points))) {
      point_row <- parent_points[i, , drop = FALSE]
      point_row_raster <- parent_points_raster[i, , drop = FALSE]
      piece_idx <- select_piece_for_point(pieces, point_row)

      split_result <- split_piece_once(
        pieces[[piece_idx]],
        point_row,
        point_row_raster,
        fdr_parent,
        fac_parent,
        local_window_cache  = local_window_cache,
        use_fdr_fallback    = use_fdr_fallback,
        fdr_parent_fb       = fdr_parent_fb,
        fac_parent_fb       = fac_parent_fb,
        is_two_way_split    = nrow(parent_points) == 1L
      )
      if (is.null(split_result)) next

      replacement <- lapply(seq_len(nrow(split_result)), function(k) split_result[k, , drop = FALSE])
      pieces <- c(pieces[-piece_idx], replacement)
      successful <- successful + 1L
    }

    # Reference-residual fallback: if fewer splits succeeded than attempted,
    # some sibling flowpaths have no piece assigned. For POI-split siblings of
    # the same parent COMID, the upstream-most sibling's watershed IS the
    # residual of the parent polygon after removing all successfully-delineated
    # sibling pieces. Assign the residual to any unrepresented sibling.
    if (successful < attempted) {
      assigned_ids <- unique(as.character(
        do.call(rbind, pieces)$flowpath_id))
      expected_ids <- unique(as.character(parent_points$flowpath_id))
      missing_ids  <- setdiff(expected_ids, assigned_ids)

      if (length(missing_ids) > 0L && length(pieces) > 0L) {
        tryCatch(
          {
            all_pieces_sf <- do.call(rbind, pieces)
            all_pieces_sf <- sf::st_as_sf(all_pieces_sf)
            pieces_union  <- sf::st_union(all_pieces_sf)
            parent_geom   <- sf::st_geometry(parent_divide)[[1L]]
            parent_sfc    <- sf::st_sfc(parent_geom, crs = sf::st_crs(all_pieces_sf))
            residual      <- suppressWarnings(
              sf::st_difference(parent_sfc, pieces_union))
            residual_area <- as.numeric(sf::st_area(residual))
            if (is.finite(residual_area) && residual_area > 0) {
              # Assign residual to the first missing sibling (most upstream)
              stub <- all_pieces_sf[1L, , drop = FALSE]
              sf::st_geometry(stub) <- residual
              stub$divide_id    <- missing_ids[1L]
              stub$flowpath_id  <- missing_ids[1L]
              pieces <- c(pieces, list(stub))
              message(sprintf(
                "split fallback: assigned reference residual (%.3f sqkm) to %s",
                residual_area / 1e6, missing_ids[1L]))
            }
          },
          error = function(e) NULL)
      }
    }

    list(
      parent_id = parent_id,
      attempted = attempted,
      successful = successful,
      pieces = do.call(rbind, pieces)
    )
  }

  # -- PURE ENGINE: single masked multi-source D8 partition per parent --------
  process_parent_pure <- function(parent_divide, parent_points, parent_points_raster,
                                  parent_id, fdr_base, fac_base,
                                  fdr_fallback_base = NULL, fac_fallback_base = NULL) {
    parent_divide <- coerce_to_divide_schema(parent_divide)
    .tr <- Sys.getenv("HF_SPLIT_TRACE") == "1"
    .trace <- function(msg) if (.tr) message(sprintf("[pp %s] %s", parent_id, msg))
    key  <- paste(parent_points$flowpath_id, parent_points$downstream_id, sep = "|")
    keep <- !duplicated(key)
    parent_points <- parent_points[keep, ]
    parent_points_raster <- parent_points_raster[keep, ]
    # Most upstream first, by the members' own chain. The collision rule below
    # keeps the FIRST outlet, so this order decides which sibling keeps contested
    # ground, it has to mean something. Falls back to the previous id sort only
    # when the members do not form a simple chain, and says so rather than
    # letting an arbitrary order pass for a topological one.
    ord <- .topo_order_members(parent_points$flowpath_id,
                               parent_points$downstream_id)
    if (is.null(ord)) {
      ord <- order(parent_points$flowpath_id, parent_points$downstream_id)
      .trace("members do not form a simple chain; falling back to id order")
    }
    parent_points <- parent_points[ord, ]
    parent_points_raster <- parent_points_raster[ord, ]
    attempted <- nrow(parent_points)
    # Every abandonment names its reason. A parent that returns unsplit has its
    # whole catchment handed to the downstream-most child by PRD, so the upstream
    # fragment silently ships owning no land, and with no reason recorded, that
    # took a night to attribute across seven reaches on VPU 04. The reason is
    # always reported, not gated behind a trace flag, because it is the only
    # record that a catchment was surrendered.
    fail <- function(why = NA_character_) {
      if (!is.na(why))
        message(sprintf("split_divides_d8: parent %s kept unsplit -- %s",
          parent_id, why))
      list(parent_id = parent_id, attempted = attempted,
        successful = 0L, pieces = parent_divide, reason = why)
    }

    parent_raster <- tryCatch(sf::st_transform(parent_divide, raster_crs), error = function(e) NULL)
    if (is.null(parent_raster) || nrow(parent_raster) == 0) return(fail("reproject parent"))
    parent_vect <- tryCatch(terra::vect(parent_raster), error = function(e) NULL)
    if (is.null(parent_vect)) return(fail("vectorize parent"))
    parent_buf  <- tryCatch(terra::buffer(parent_vect, width = 250), error = function(e) NULL)
    if (is.null(parent_buf)) return(fail("buffer parent"))
    fdr_parent  <- tryCatch(terra::crop(fdr_base, parent_buf, snap = "out", mask = TRUE),
      error = function(e) NULL)
    fac_parent  <- tryCatch(terra::crop(fac_base, parent_buf, snap = "out", mask = TRUE),
      error = function(e) NULL)

    # Whole-parent emptiness is the cheap half of the fallback rule: there is
    # nothing to snap to anywhere, so the doomed snap is skipped entirely. The
    # other half, a grid that is sound over the parent but holed exactly at an
    # outlet, cannot be seen from here and is handled after the snap has tried.
    fb_in_use <- FALSE
    primary_empty <- is.null(fdr_parent) ||
      all(is.na(terra::values(fdr_parent, mat = FALSE)))
    if (primary_empty && !is.null(fdr_fallback_base) && !is.null(fac_fallback_base)) {
      fdr_parent <- tryCatch(terra::crop(
        fdr_fallback_base, terra::project(parent_buf, terra::crs(fdr_fallback_base)),
        snap = "out", mask = TRUE), error = function(e) NULL)
      fac_parent <- tryCatch(terra::crop(
        fac_fallback_base, terra::project(parent_buf, terra::crs(fac_fallback_base)),
        snap = "out", mask = TRUE), error = function(e) NULL)
      primary_empty <- is.null(fdr_parent) ||
        all(is.na(terra::values(fdr_parent, mat = FALSE)))
      # The outlets were projected to the PRIMARY grid before this point, so they
      # have to follow the raster onto the fallback's CRS; otherwise the snap
      # looks up coordinates that do not exist in the grid it is searching.
      if (!primary_empty) {
        moved <- tryCatch(
          sf::st_transform(parent_points_raster, sf::st_crs(terra::crs(fdr_parent))),
          error = function(e) NULL)
        if (is.null(moved)) return(fail("reproject outlets onto the fallback grid"))
        parent_points_raster <- moved
        fb_in_use <- TRUE
      }
    }
    # Three different failures used to share one message, which made the log
    # unable to say whether the grid was missing, unreadable, or merely empty.
    if (is.null(fdr_parent)) return(fail("FDR could not be cropped to the parent"))
    if (is.null(fac_parent)) return(fail("FAC could not be cropped to the parent"))
    if (primary_empty) return(fail("no flow direction anywhere over the parent"))
    .trace(sprintf("cropped FDR ncell=%d (%dx%d)", terra::ncell(fdr_parent),
      terra::nrow(fdr_parent), terra::ncol(fdr_parent)))

    res_x <- terra::res(fac_parent)[1]
    .trace(sprintf("res_x=%.6g; snapping %d outlet(s); fdr_inmem=%s",
      res_x, nrow(parent_points), terra::inMemory(fdr_parent)))
    pts_xy <- sf::st_coordinates(parent_points_raster)

    # The outlet snap is deliberately reach-blind. Constraining candidates to
    # cells within 30 m of the reach, with or without a drainage-area preference
    # among them, was tried and reverted: on Alaska it left the containment rate
    # unchanged at 1945 of 2113 pieces and starved children at 253 against 260,
    # while stage-1 divide overlap went from 34.81 to 195.42 km2 because the
    # pieces stopped tiling their parent cleanly. Narrowing the search to the same
    # two-cell window the blind rule uses changed nothing (byte-identical
    # containment, 101.73 km2 of overlap at Step 6d), so the widening was not the
    # harm: picking a different near cell is. See the tally below for how the
    # accuracy of this snap is measured.
    outlet_cells <- integer(0)
    outlet_labels <- integer(0)
    for (i in seq_len(nrow(parent_points))) {
      .trace(sprintf("snap outlet %d/%d", i, nrow(parent_points)))
      cell <- .snap_outlet_cell(pts_xy[i, "X"], pts_xy[i, "Y"],
        fac_parent, fdr_parent, snap_distance_m, res_x)
      if (is.na(cell)) next
      outlet_cells  <- c(outlet_cells, cell)
      outlet_labels <- c(outlet_labels, i)
    }
    # Second half of the fallback rule. Asking only whether the primary is empty
    # over the whole parent misses the case that actually surrenders a
    # catchment: a hole in the grid exactly over an outlet, with sound data all
    # around it. Measured on VPU 04, three parents fail that way: 613 to 3002
    # valid cells in the window, NoData at the split point, and their outlets
    # carry a valid direction in the fallback grid.
    #
    # The retry runs only where the primary has already failed to snap anything,
    # so it can take nothing away: without it the parent is abandoned, PRD hands
    # the whole catchment to the downstream-most child, and every other member
    # ships owning no land.
    if (length(outlet_cells) == 0L && !fb_in_use &&
        !is.null(fdr_fallback_base) && !is.null(fac_fallback_base)) {
      fdr_fb <- tryCatch(terra::crop(fdr_fallback_base,
        terra::project(parent_buf, terra::crs(fdr_fallback_base)),
        snap = "out", mask = TRUE), error = function(e) NULL)
      fac_fb <- tryCatch(terra::crop(fac_fallback_base,
        terra::project(parent_buf, terra::crs(fac_fallback_base)),
        snap = "out", mask = TRUE), error = function(e) NULL)
      moved <- if (is.null(fdr_fb) || is.null(fac_fb)) NULL else tryCatch(
        sf::st_transform(parent_points_raster, sf::st_crs(terra::crs(fdr_fb))),
        error = function(e) NULL)
      if (!is.null(moved)) {
        xy_fb  <- sf::st_coordinates(moved)
        res_fb <- terra::res(fac_fb)[1]
        for (i in seq_len(nrow(parent_points))) {
          cell <- .snap_outlet_cell(xy_fb[i, "X"], xy_fb[i, "Y"],
            fac_fb, fdr_fb, snap_distance_m, res_fb)
          if (is.na(cell)) next
          outlet_cells  <- c(outlet_cells, cell)
          outlet_labels <- c(outlet_labels, i)
        }
        # Adopt the grid only if it answered. When both are dry the primary is
        # kept, so nothing is traded for a grid that is no better.
        if (length(outlet_cells)) {
          fdr_parent <- fdr_fb
          fac_parent <- fac_fb
          parent_points_raster <- moved
          pts_xy <- xy_fb
          res_x  <- res_fb
          fb_in_use <- TRUE
          message(sprintf(paste0("split_divides_d8: parent %s split on the ",
            "fallback grid -- the primary has no flow direction at its outlet(s)"),
            parent_id))
        }
      }
    }
    if (length(outlet_cells) == 0L) return(fail("no outlet snapped to a grid cell"))
    .trace(sprintf("snapped %d outlet cell(s)", length(outlet_cells)))
    # Colliding snaps: keep the first, which `.topo_order_members` has made the
    # upstream-most member of the chain; the other sub-piece's cells resolve to
    # the next outlet / terminal, no piece lost.
    dup <- duplicated(outlet_cells)
    if (any(dup)) {
      outlet_cells <- outlet_cells[!dup]
      outlet_labels <- outlet_labels[!dup]
    }

    fps <- as.character(parent_points$flowpath_id)
    dss <- as.character(parent_points$downstream_id)
    term_candidates <- setdiff(dss, fps)               # the sink (chain tail)
    terminal_id    <- if (length(term_candidates) >= 1L) term_candidates[1L] else dss[length(dss)]
    terminal_label <- nrow(parent_points) + 1L
    label_to_id    <- stats::setNames(fps, as.character(seq_along(fps)))

    fdr_mat   <- terra::as.matrix(fdr_parent, wide = TRUE)
    rc        <- terra::rowColFromCell(fdr_parent, outlet_cells)
    if (is.null(dim(rc))) rc <- matrix(rc, ncol = 2L)
    .trace("as.matrix done; starting D8 label")
    label_mat <- .label_parent_d8(fdr_mat, rc, outlet_labels, terminal_label)
    if (is.null(label_mat)) return(fail("D8 labelling returned nothing"))
    .trace("D8 label done; starting labels_to_pieces (polygonize+clip)")

    lab_rast <- terra::rast(fdr_parent)
    lv <- as.vector(t(label_mat))
    lv[lv == 0L] <- NA_integer_
    terra::values(lab_rast) <- lv

    # Grow labels into unlabelled cells that lie inside the parent.
    #
    # A cell with no flow direction (FDR no-data inside the parent) never gets
    # a label, so polygonisation skips it and its area is left unclaimed, to be
    # picked up later by the residual pass. Measured on VPU 04: one parent
    # (15568145) held 11 such cells, 9936 m2, the whole of that VPU's unclaimed
    # area. Filling them moves coverage toward being a property of the labelling
    # instead of something a later pass repairs; it took the unclaimed total to
    # 3064 m2. It does not reach zero, because a strip of the parent narrower than
    # a cell has no cell centre to claim it, a rasterisation limit no amount of
    # label growing removes, so the residual pass still has work to do.
    #
    # Restricted to cells inside the parent, so the ordinary NA ring outside the
    # crop buffer is untouched and parents with no interior gap do no work.
    # Bounded passes, and pieces are clipped to the parent afterwards, so no area
    # outside it can be invented.
    pmask <- tryCatch(terra::rasterize(terra::vect(parent_raster), lab_rast,
      field = 1), error = function(e) NULL)
    if (!is.null(pmask)) {
      inside <- !is.na(terra::values(pmask, mat = FALSE))
      need   <- is.na(terra::values(lab_rast, mat = FALSE)) & inside
      if (any(need, na.rm = TRUE)) {
        if (Sys.getenv("HF_PARTITION_LEDGER") == "1")
          message(sprintf("[ledger-gap] parent=%s unlabelled_cells_in_parent=%d",
            as.character(parent_divide$divide_id[1L]), sum(need, na.rm = TRUE)))
        for (.pass in seq_len(8L)) {
          fl <- tryCatch(terra::focal(lab_rast, w = 3, fun = "modal",
            na.policy = "only", na.rm = TRUE), error = function(e) NULL)
          if (is.null(fl)) break
          nv <- terra::values(lab_rast, mat = FALSE)
          fv <- terra::values(fl, mat = FALSE)
          take <- is.na(nv) & inside & !is.na(fv)
          if (!any(take, na.rm = TRUE)) break
          nv[take] <- fv[take]
          terra::values(lab_rast) <- nv
          if (!any(is.na(nv) & inside, na.rm = TRUE)) break
        }
      }
    }

    pieces <- tryCatch(
      .labels_to_pieces(lab_rast, parent_divide, label_to_id, terminal_label, terminal_id),
      error = function(e) {
        message(sprintf("pure partition %s: %s", parent_id, conditionMessage(e)))
        NULL
      })
    if (is.null(pieces) || nrow(pieces) == 0L) return(fail("polygonize produced no pieces"))

    # Attach the parent's other attributes; set ids per piece; coerce to schema.
    # Replicate the parent row (preserving its geometry-column NAME) then swap in
    # the piece geometry, so coerce_to_divide_schema's geom_col lookup matches
    # regardless of whether the source used "geom" or "geometry".
    out <- parent_divide[rep(1L, nrow(pieces)), , drop = FALSE]
    sf::st_geometry(out) <- sf::st_geometry(pieces)
    out$divide_id   <- as.character(pieces$divide_id)
    out$flowpath_id <- as.character(pieces$flowpath_id)
    out <- coerce_to_divide_schema(out)

    # Does each piece hold the reach it was cut for? Reported, not repaired: the
    # partition is a total function over the parent, so a piece in the wrong place
    # means its OUTLET was snapped to the wrong cell, and re-cutting the parent
    # here would only move the error. Counting it is what makes the snap's
    # accuracy visible: stage 3 otherwise clears these pairings silently.
    n_checked <- 0L; n_on_reach <- 0L
    if (!is.null(reach_lines) && "line_id" %in% names(reach_lines)) {
      lid <- as.character(reach_lines$line_id)
      for (r in seq_len(nrow(out))) {
        kk <- match(as.character(out$flowpath_id[r]), lid)
        if (is.na(kk)) next
        n_checked  <- n_checked + 1L
        n_on_reach <- n_on_reach +
          isTRUE(.piece_holds_reach(sf::st_geometry(out)[r], sf::st_geometry(reach_lines)[kk]))
      }
    }

    list(parent_id = parent_id, attempted = attempted,
      successful = attempted, pieces = out,
      checked = n_checked, on_reach = n_on_reach)
  }

  # Engine dispatch: both runners call this.
  .pp <- if (engine == "pure") process_parent_pure else process_parent

  # Process each parent divide
  divides_to_keep <- divides |> dplyr::filter(!divide_id %in% parents_to_split)
  divides_subset <- divides |> dplyr::filter(divide_id %in% parents_to_split)
  points_subset <- split_points |> dplyr::filter(divide_id %in% parents_to_split)
  points_subset_raster <- split_points_raster |> dplyr::filter(divide_id %in% parents_to_split)
  divides_by_id <- split(divides_subset, divides_subset$divide_id)
  split_points_by_id <- split(points_subset, points_subset$divide_id)
  split_points_raster_by_id <- split(points_subset_raster, points_subset_raster$divide_id)
  parent_tasks <- lapply(parents_to_split, function(pid) {
    list(
      parent_id = pid,
      parent_divide = divides_by_id[[pid]],
      parent_points = split_points_by_id[[pid]],
      parent_points_raster = split_points_raster_by_id[[pid]]
    )
  })

  use_rast_fallback <- !is.null(fallback_fdr_path) && nzchar(fallback_fdr_path) &&
    !is.null(fallback_fac_path) && nzchar(fallback_fac_path) &&
    !identical(normalizePath(fdr_path, mustWork = FALSE),
      normalizePath(fallback_fdr_path, mustWork = FALSE))

  if (!use_parallel) {
    # In sequential mode open rasters once here (no worker export overhead).
    fdr_main <- tryCatch(terra::rast(fdr_path), error = function(e) NULL)
    fac_main <- tryCatch(terra::rast(fac_path), error = function(e) NULL)
    fdr_fb_main <- if (use_rast_fallback) tryCatch(terra::rast(fallback_fdr_path), error = function(e) NULL) else NULL
    fac_fb_main <- if (use_rast_fallback) tryCatch(terra::rast(fallback_fac_path), error = function(e) NULL) else NULL
    runner_seq <- function(task) {
      # Diagnostic trace (gated, off by default): names the parent BEFORE .pp
      # runs, so a non-terminating parent, which never returns and so never
      # hits the SLOW/ABORTED logs below: is identifiable as the last START
      # line in the log. Enable with HF_SPLIT_TRACE=1.
      if (Sys.getenv("HF_SPLIT_TRACE") == "1")
        message(sprintf("[simple_split] START parent %s (%d pts)",
          task$parent_id, nrow(task$parent_points)))
      .t0 <- Sys.time()
      # Per-parent hard timeout. A masked-D8 partition that never returns
      # (pathological FDR / degenerate geometry) would otherwise hang the whole
      # VPU at 100% CPU; the >20s timer below can't fire because .pp never
      # returns. On timeout keep the parent UNSPLIT (whole divide, no area loss,
      # matching .pp's own fail() fallback) and log the offender.
      setTimeLimit(elapsed = 300, transient = TRUE)
      .res <- tryCatch(
        .pp(
          parent_divide = task$parent_divide,
          parent_points = task$parent_points,
          parent_points_raster = task$parent_points_raster,
          parent_id = task$parent_id,
          fdr_base = fdr_main,
          fac_base = fac_main,
          fdr_fallback_base = fdr_fb_main,
          fac_fallback_base = fac_fb_main
        ),
        error = function(e) {
          message(sprintf(
            "[simple_split] parent %s ABORTED (300s cap): %s -- kept unsplit",
            task$parent_id, conditionMessage(e)))
          list(parent_id = task$parent_id,
            attempted = nrow(task$parent_points),
            successful = 0L, pieces = task$parent_divide)
        }
      )
      setTimeLimit()  # clear the limit
      # Also flag merely-slow (but returning) parents for profiling.
      .el <- as.numeric(difftime(Sys.time(), .t0, units = "secs"))
      if (.el > 20) message(sprintf(
        "[simple_split] SLOW parent %s: %.1fs (%d split point(s))",
        task$parent_id, .el, nrow(task$parent_points)))
      .res
    }
  } else {
    # In parallel mode, do not create large raster objects here. Workers will
    # lazily initialize rasters from disk inside runner_parallel.
    runner_seq <- function(task) {
      .pp(
        parent_divide = task$parent_divide,
        parent_points = task$parent_points,
        parent_points_raster = task$parent_points_raster,
        parent_id = task$parent_id,
        fdr_base = NULL,
        fac_base = NULL,
        fdr_fallback_base = NULL,
        fac_fallback_base = NULL
      )
    }
  }

  runner_parallel <- local({
    fdr_worker    <- NULL
    fac_worker    <- NULL
    fdr_fb_worker <- NULL
    fac_fb_worker <- NULL
    function(task) {
      if (is.null(fdr_worker)) fdr_worker <<- terra::rast(fdr_path)
      if (is.null(fac_worker)) fac_worker <<- terra::rast(fac_path)
      if (use_rast_fallback && is.null(fdr_fb_worker))
        fdr_fb_worker <<- tryCatch(terra::rast(fallback_fdr_path), error = function(e) NULL)
      if (use_rast_fallback && is.null(fac_fb_worker))
        fac_fb_worker <<- tryCatch(terra::rast(fallback_fac_path), error = function(e) NULL)
      .pp(
        parent_divide = task$parent_divide,
        parent_points = task$parent_points,
        parent_points_raster = task$parent_points_raster,
        parent_id = task$parent_id,
        fdr_base = fdr_worker,
        fac_base = fac_worker,
        fdr_fallback_base = fdr_fb_worker,
        fac_fallback_base = fac_fb_worker
      )
    }
  })
  runner <- if (use_parallel) runner_parallel else runner_seq

  results <- if (use_parallel) {
    tryCatch(
      parallel::mclapply(parent_tasks, runner, mc.cores = effective_workers,
        mc.silent = FALSE),
      error = function(e) {
        .hf_log_warn(sprintf("mclapply failed (%s); retrying sequentially", conditionMessage(e)))
        lapply(parent_tasks, runner_seq)
      }
    )
  } else {
    lapply(parent_tasks, runner)
  }

  collect_results <- function(results_list) {
    new_divides_list <- lapply(results_list, `[[`, "pieces")
    n_attempted <- sum(vapply(results_list, `[[`, integer(1), "attempted"))
    n_successful <- sum(vapply(results_list, `[[`, integer(1), "successful"))
    n_parents_success <- sum(vapply(results_list, function(x) x$successful > 0, logical(1)))

    # Reported here rather than from inside the workers: parents are processed
    # under parallel::mclapply, whose message() output never reaches the build
    # log, so any instrumentation down there is invisible. The counts ride back
    # with the results instead.
    .tally <- function(f) sum(vapply(results_list,
      function(x) { v <- x[[f]]; if (is.null(v)) 0L else as.integer(v) }, integer(1)))
    n_checked <- .tally("checked")
    if (n_checked > 0L) {
      n_on <- .tally("on_reach")
      message(sprintf(
        "split_divides_d8: %d of %d split piece(s) hold their own reach (%.1f%%); %d do not",
        n_on, n_checked, 100 * n_on / n_checked, n_checked - n_on))
    }

    if (length(new_divides_list) > 0) {
      new_divides <- do.call(rbind, new_divides_list) |> coerce_to_divide_schema()
      result <- dplyr::bind_rows(divides_to_keep, new_divides) |> coerce_to_divide_schema()
    } else {
      result <- divides_to_keep |> coerce_to_divide_schema()
    }

    # Whole-function area ledger: everything handed in against everything handed
    # back. Measured exact (0.000000 km2) on VPU 04 and VPU 13, which is what
    # located the remaining coverage loss downstream of this function.
    if (Sys.getenv("HF_PARTITION_LEDGER") == "1") {
      a_in  <- sum(as.numeric(sf::st_area(divides)), na.rm = TRUE)
      a_out <- sum(as.numeric(sf::st_area(result)), na.rm = TRUE)
      message(sprintf(paste0(
        "[ledger-sd] in_divides=%d in_area_km2=%.6f | out_divides=%d ",
        "out_area_km2=%.6f | delta_km2=%.6f"),
        nrow(divides), a_in / 1e6, nrow(result), a_out / 1e6, (a_out - a_in) / 1e6))
    }

    list(
      result = result,
      n_attempted = n_attempted,
      n_successful = n_successful,
      n_parents_success = n_parents_success
    )
  }

  stats <- collect_results(results)

  if (use_parallel && stats$n_attempted > 0 && stats$n_successful == 0) {
    .hf_log_warn("Parallel split returned 0 successful splits; retrying sequentially inside splitter.")
    results <- lapply(parent_tasks, runner_seq)
    stats <- collect_results(results)
  }

  message(sprintf(
    "Split complete: %d/%d successful point splits across %d/%d parent divides; %d total divides (was %d)",
    stats$n_successful,
    stats$n_attempted,
    stats$n_parents_success,
    length(parents_to_split),
    nrow(stats$result),
    nrow(divides)
  ))

  return(stats$result)
}

# -- PURE ENGINE helpers (masked multi-source D8 partition) -------------------
# These are standalone (no closure deps) so they are unit-testable in isolation.

# Esri D8 downstream offset: for a cell's own flow code, the (drow, dcol) to the
# cell it flows INTO. Inverse of the inbound-neighbour table used above.
.d8_ds_offset <- function() {
  list(`1` = c(0L, 1L),  `2` = c(1L, 1L),  `4` = c(1L, 0L),  `8` = c(1L, -1L),
    `16` = c(0L, -1L), `32` = c(-1L, -1L), `64` = c(-1L, 0L), `128` = c(-1L, 1L))
}

#' Label every parent cell by the first sub-outlet downstream of it (D8).
#'
#' In a D8 grid every cell has exactly one downstream cell, so "first outlet
#' reached along the downstream chain" is a total function, each cell gets
#' exactly one label. Non-overlap and full coverage therefore hold by
#' construction. Cells whose chain exits the parent (NA / edge / invalid code /
#' cycle from a flat) get `terminal_label`. Memoized stack unwind -> O(cells).
#'
#' @param fdr_mat integer matrix (terra::as.matrix wide=TRUE), NA outside parent.
#' @param outlet_rc 2-col matrix (row, col) of snapped sub-outlet cells.
#' @param outlet_labels integer label for each outlet row (1..K).
#' @param terminal_label integer label for the terminal (chain-tail) piece.
#' @return integer matrix of labels (same dims as fdr_mat); 0 = NA/unlabeled.
#' @keywords internal
#' @noRd
.label_parent_d8 <- function(fdr_mat, outlet_rc, outlet_labels, terminal_label) {
  nr <- nrow(fdr_mat)
  nc <- ncol(fdr_mat)
  ncell <- nr * nc
  if (ncell == 0L) return(NULL)
  ds <- .d8_ds_offset()

  nonNA <- which(!is.na(fdr_mat))                 # column-major linear indices
  if (length(nonNA) == 0L) return(NULL)
  r0   <- ((nonNA - 1L) %% nr) + 1L
  c0   <- ((nonNA - 1L) %/% nr) + 1L
  code <- fdr_mat[nonNA]

  drow <- rep(NA_integer_, length(code))
  dcol <- drow
  for (k in names(ds)) {
    sel <- !is.na(code) & code == as.integer(k)
    if (any(sel)) {
      drow[sel] <- ds[[k]][1]
      dcol[sel] <- ds[[k]][2]
    }
  }
  r1  <- r0 + drow
  c1 <- c0 + dcol
  inb <- !is.na(drow) & r1 >= 1L & r1 <= nr & c1 >= 1L & c1 <= nc
  nxt <- integer(ncell)                            # 0 = exit the parent
  if (any(inb)) {
    wl   <- which(inb)
    dlin <- (c1[wl] - 1L) * nr + r1[wl]            # column-major linear of downstream
    ok   <- !is.na(fdr_mat[dlin])
    nxt[nonNA[wl[ok]]] <- dlin[ok]
  }

  label <- integer(ncell)                          # 0 = unassigned
  out_lin <- (outlet_rc[, 2] - 1L) * nr + outlet_rc[, 1]
  label[out_lin] <- as.integer(outlet_labels)      # seed outlets (inclusive)

  onpath <- logical(ncell)
  path   <- integer(ncell)
  for (s in nonNA) {
    if (label[s] != 0L) next
    np <- 0L
    cur <- s
    resolved <- terminal_label
    repeat {
      if (cur == 0L) {
        resolved <- terminal_label
        break
      }   # exited parent
      if (label[cur] != 0L) {
        resolved <- label[cur]
        break
      }   # hit seeded/resolved
      if (onpath[cur]) {
        resolved <- terminal_label
        break
      }   # cycle (flat)
      onpath[cur] <- TRUE
      np <- np + 1L
      path[np] <- cur
      cur <- nxt[cur]
    }
    if (np > 0L) {
      idx <- path[1:np]
      label[idx] <- resolved
      onpath[idx] <- FALSE
    }
  }
  matrix(label, nr, nc)
}

#' Snap a sub-outlet (x,y in raster CRS) to a channel cell inside the parent.
#' Bounded progressive focal max-FAC search; falls back to the nearest in-parent
#' FDR cell to the raw coordinate (never a bisector). Returns a terra cell index
#' in `fdr_parent`, or NA.
#' @keywords internal
#' @noRd
.snap_outlet_cell <- function(x, y, fac_parent, fdr_parent, snap_distance_m, res_x) {
  if (!is.finite(x) || !is.finite(y) || !is.finite(res_x) || res_x <= 0) return(NA_integer_)
  cell_at <- function(cx, cy) {
    cell <- tryCatch(terra::cellFromXY(fdr_parent, cbind(cx, cy)), error = function(e) NA_integer_)
    if (length(cell) != 1L || is.na(cell)) return(NA_integer_)
    v <- tryCatch(terra::extract(fdr_parent, cell)[1, 1], error = function(e) NA)
    if (is.na(v)) NA_integer_ else as.integer(cell)
  }

  # 1. Local-thalweg snap: max-FAC within a TIGHT radius (<= 2 cells). A small
  #    radius nudges the outlet onto its own reach's channel cell while making
  #    it impossible for a low-order sub-outlet to jump to a higher-accumulation
  #    river that merely passes within snap_distance_m (which would misroute a
  #    tiny headwater reach onto an adjacent mainstem).
  #
  #    hyRefactor's `get_row_col()` instead keeps the outlet's own cell unless an
  #    8-neighbour carries more than twice its accumulation. That was tried here
  #    and backed out on measurement: holding everything else constant it moved
  #    starved split children only 200 -> 194, while overlap into the clean step
  #    went 19.67 -> 35.60 km2 and divides collapsed to nothing went 4 -> 36. The
  #    same run recovered 60 flowpaths from a change to the POI split repair, so
  #    the improvement first attributed to this rule was almost entirely not it.
  tight <- 2 * res_x
  fc <- tryCatch(terra::crop(fac_parent,
    terra::ext(x - tight, x + tight, y - tight, y + tight)),
  error = function(e) NULL)
  if (!is.null(fc)) {
    fv <- terra::values(fc, mat = FALSE)
    fv[!is.finite(fv) | fv > 1e300] <- NA
    if (length(fv) && !all(is.na(fv))) {
      xy   <- terra::xyFromCell(fc, which.max(fv))
      cell <- cell_at(xy[1, 1], xy[1, 2])
      if (!is.na(cell)) return(cell)
    }
  }

  # 2. Tight window had no FAC: use the outlet's own cell if it carries FDR.
  c0 <- cell_at(x, y)
  if (!is.na(c0)) return(c0)

  # 3. Else snap to the NEAREST valid-FDR cell within snap_distance_m (minimal
  #    movement, not max-FAC) so we never reach for a distant channel. Crop ONCE
  #    to the full search window (bounded by the parent raster) and take the
  #    global nearest: equivalent to a radius-expanding search but a SINGLE
  #    terra::crop rather than up to floor(snap_distance_m/res_x) crops. A
  #    per-radius loop would grind through thousands of crops for degenerate
  #    parents whose outlets sit in NA-FDR territory far from any channel (e.g.
  #    synthetic terminal/coastal divides).
  r  <- snap_distance_m
  fd <- tryCatch(terra::crop(fdr_parent, terra::ext(x - r, x + r, y - r, y + r)),
    error = function(e) NULL)
  if (!is.null(fd)) {
    v   <- terra::values(fd, mat = FALSE)
    idx <- which(!is.na(v))
    if (length(idx)) {
      xy <- terra::xyFromCell(fd, idx)
      nn <- which.min((xy[, 1] - x)^2 + (xy[, 2] - y)^2)
      cell <- cell_at(xy[nn, 1], xy[nn, 2])
      if (!is.na(cell)) return(cell)
    }
  }
  NA_integer_
}

#' Assign each residual polygon part to the piece with the longest shared edge.
#' (parent -  union pieces, the thin raster-vs-vector boundary band). Area-conserving.
#' @keywords internal
#' @noRd
.assign_residual_to_neighbor <- function(pieces, resid_geom) {
  crs   <- sf::st_crs(pieces)
  parts <- suppressWarnings(sf::st_cast(sf::st_sfc(resid_geom, crs = crs), "POLYGON"))
  geo   <- sf::st_geometry(pieces)
  pg_bnd <- sf::st_boundary(geo)
  for (p in seq_along(parts)) {
    pb <- sf::st_boundary(parts[p])
    ln <- vapply(seq_along(geo), function(t) {
      ii <- suppressWarnings(sf::st_intersection(pg_bnd[t], pb))
      sum(suppressWarnings(as.numeric(sf::st_length(ii))), na.rm = TRUE)
    }, numeric(1))
    if (length(ln) == 0L) next
    # A part sharing no boundary with any piece was skipped, which discarded it.
    # Fall back to the nearest piece: it belongs somewhere, and the partition has
    # to stay exhaustive.
    t <- if (max(ln) > 0) which.max(ln) else {
      dd <- suppressWarnings(as.numeric(sf::st_distance(parts[p], geo)))
      if (all(!is.finite(dd))) next else which.min(dd)
    }
    geo[t] <- sf::st_union(sf::st_sfc(c(list(geo[[t]]), list(parts[[p]])),
      crs = crs))[[1L]]
  }
  sf::st_set_geometry(pieces, geo)
}

#' Polygonize a label raster into one clean POLYGON per owning id.
#' Clips to the vector parent, assigns the boundary residual to longest-edge
#' neighbours, and enforces single-POLYGON + sliver-merge via `.normalize_divides`
#' (which also asserts area conservation). Returns sf with divide_id, flowpath_id.
#' @keywords internal
#' @noRd
.labels_to_pieces <- function(lab_rast, parent_divide, label_to_id,
                              terminal_label, terminal_id) {
  polys <- tryCatch(sf::st_as_sf(terra::as.polygons(lab_rast, dissolve = TRUE)),
    error = function(e) NULL)
  if (is.null(polys) || nrow(polys) == 0L) return(NULL)
  polys$.lab <- as.integer(polys[[names(polys)[1]]])
  polys$owning_id <- ifelse(polys$.lab == terminal_label, terminal_id,
    unname(label_to_id[as.character(polys$.lab)]))
  polys <- polys[!is.na(polys$owning_id), , drop = FALSE]
  if (nrow(polys) == 0L) return(NULL)

  out_crs <- sf::st_crs(parent_divide)
  polys <- sf::st_transform(polys, out_crs)
  pg    <- sf::st_geometry(parent_divide)
  polys <- suppressWarnings(sf::st_intersection(polys, pg))
  polys <- suppressWarnings(sf::st_collection_extract(polys, "POLYGON"))
  if (nrow(polys) == 0L) return(NULL)

  ids <- unique(polys$owning_id)
  g   <- lapply(ids, function(id)
    sf::st_union(sf::st_geometry(polys[polys$owning_id == id, ]))[[1L]])
  pieces <- sf::st_sf(divide_id = ids, flowpath_id = ids,
    geometry = sf::st_sfc(g, crs = out_crs))

  # Close the raster-vs-vector boundary band so union pieces == parent exactly.
  # "Exactly" is the point: these pieces are the partition every downstream area
  # figure rests on, so anything outside them is coverage that silently vanishes.
  #
  # Every residual part is assigned, at any size. Thresholding the band or taking
  # only its first part drops area that is individually invisible and additive.
  #
  # The band is inherent to comparing a rasterised label field against a vector
  # parent: cells outside the parent are masked to NA, so a thin strip can have no
  # cell to claim it. Assigning it decides where the band goes rather than
  # repairing a computed value. Labelling a window that fully covers the parent
  # would make this a no-op; left as follow-up.
  ua    <- sf::st_union(sf::st_geometry(pieces))
  resid <- suppressWarnings(sf::st_difference(pg, ua))

  # Per-parent exhaustiveness ledger. The partition's defining property is that
  # its pieces tile the parent, so the quantity to watch is parent area minus the
  # union of the pieces, per parent, before anything is assigned. Summed over a VPU
  # this is exactly the coverage a build can lose, and per-parent detail is what
  # separates "the labelling did not cover the parent" from "a later step dropped a
  # piece", a distinction that took several wrong guesses to establish without it.
  if (Sys.getenv("HF_PARTITION_LEDGER") == "1") {
    p_area <- sum(as.numeric(sf::st_area(pg)), na.rm = TRUE)
    r_area <- if (length(resid) > 0L)
      sum(as.numeric(sf::st_area(resid)), na.rm = TRUE) else 0
    message(sprintf(
      "[ledger] parent=%s pieces=%d parent_m2=%.1f unclaimed_m2=%.1f frac=%.3e",
      as.character(parent_divide$divide_id[1L]), nrow(pieces), p_area, r_area,
      if (p_area > 0) r_area / p_area else NA_real_))
  }

  if (length(resid) > 0L &&
    sum(as.numeric(sf::st_area(resid)), na.rm = TRUE) > 0) {
    resid_all <- suppressWarnings(sf::st_union(resid))
    if (length(resid_all) > 0L)
      pieces <- .assign_residual_to_neighbor(pieces, resid_all[[1L]])
  }

  # Single POLYGON per piece + merge internal slivers (+ area assertion).
  .normalize_divides(pieces, id_col = "divide_id")
}
