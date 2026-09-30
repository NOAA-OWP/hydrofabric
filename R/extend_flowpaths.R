#' Extend aggregated flowpath geometry along its main-stem sub-flowlines
#'
#' @description
#' `aggregate_to_distribution` often returns a flowpath whose geometry only
#' covers the outlet portion of its catchment. Sub-flowlines attributed to the
#' flowpath (via `add_flowlines`) may extend upstream along eligible branches
#' to the divide headwater but are not part of the aggregated flowpath line.
#' This function rebuilds each flowpath's geometry as the main-stem path
#' through its sub-flowlines: from the outlet, walking upstream at each
#' confluence along the branch with the largest `tot_drainage_areasqkm` until
#' reaching a headwater.
#'
#' @param input_gpkg Path to the intermediate GeoPackage (must contain
#'   `flowpaths` and `flowlines` layers).
#' @param flowpath_layer Flowpath layer name. Default `"flowpaths"`.
#' @param flowlines_layer Sub-flowline layer name. Default `"flowlines"`.
#' @param divides_layer Retained for backward compatibility; no longer read.
#'
#' @return Invisibly, the updated flowpaths sf object. Writes back to
#'   `input_gpkg` (overwrites flowpath_layer).
#'
#' @examples
#' \dontrun{
#' extend_flowpaths_to_mainstem("out_vpu01.gpkg")
#' }
#'
#' @export
extend_flowpaths_to_mainstem <- function(input_gpkg,
                                         flowpath_layer            = "flowpaths",
                                         flowlines_layer           = "flowlines",
                                         divides_layer             = "divides") {
  fp <- sf::st_read(input_gpkg, flowpath_layer,  quiet = TRUE)
  fl <- sf::st_read(input_gpkg, flowlines_layer, quiet = TRUE)

  # Column-name compatibility: intermediate schema uses `id`/`toid`; NextGen
  # schema uses `flowline_id`/`flowline_toid`. Normalize.
  if (!"flowline_id"   %in% names(fl) && "id"   %in% names(fl)) fl$flowline_id   <- fl$id
  if (!"flowline_toid" %in% names(fl) && "toid" %in% names(fl)) fl$flowline_toid <- fl$toid
  fl$flowline_id   <- as.character(fl$flowline_id)
  fl$flowline_toid <- as.character(fl$flowline_toid)
  fp$flowpath_id   <- as.character(fp$flowpath_id)

  target_crs <- sf::st_crs(fp)
  if (sf::st_crs(fl) != target_crs) fl <- sf::st_transform(fl, target_crs)

  fl_by_fp <- split(seq_len(nrow(fl)), as.character(fl$flowpath_id))

  # Lengths for the single-sub decision below, computed once in bulk rather than
  # per flowpath. Both layers are already in target_crs here.
  fl_len_m <- suppressWarnings(as.numeric(sf::st_length(sf::st_geometry(fl))))
  fp_len_m <- suppressWarnings(as.numeric(sf::st_length(sf::st_geometry(fp))))
  fl_len_m[is.na(fl_len_m)] <- 0
  fp_len_m[is.na(fp_len_m)] <- 0

  # Parent flowpath's mainstem levelpathid (canonical NHD main-stem key).
  fp_mainstem <- if ("mainstem" %in% names(fp))
    stats::setNames(as.character(fp$mainstem), fp$flowpath_id)
  else stats::setNames(rep(NA_character_, nrow(fp)), fp$flowpath_id)
  sub_lp <- if ("levelpathid" %in% names(fl))
    as.character(fl$levelpathid) else rep(NA_character_, nrow(fl))

  # Phase 1. Fast walk (no geometry ops): compute path_rows per flowpath
  # as a list of integer indices into `fl`. The per-flowpath work is cheap
  # pure R (sub id comparisons + drainage-area tie-breaks); heavy geometry
  # operations are deferred to Phase 2 where they run in bulk via GEOS.
  total       <- 0L
  path_by_fp  <- vector("list", nrow(fp))
  # Branch counters. This function decides the geometry of every aggregated
  # flowpath and used to report only how many it touched, so which decision
  # produced a given line was not recoverable from a run.
  n_single_skp <- 0L   # single-sub walk left unwritten: it would shorten the line
  n_truncated  <- 0L   # walk stopped early: children existed but none on-mainstem
  for (i in seq_len(nrow(fp))) {
    fp_id <- fp$flowpath_id[i]
    rows  <- fl_by_fp[[fp_id]]
    if (is.null(rows) || length(rows) == 0L) next
    total <- total + 1L

    sub_ids <- fl$flowline_id[rows]
    sub_to  <- fl$flowline_toid[rows]
    sub_tda <- as.numeric(fl$tot_drainage_areasqkm[rows])
    sub_tda[is.na(sub_tda)] <- 0
    sub_lp_local <- sub_lp[rows]
    fp_lp <- fp_mainstem[fp_id]

    outlet_local <- which(!(sub_to %in% sub_ids))
    if (length(outlet_local) == 0L) next

    # Multi-terminal all-subs branch removed. When a flowpath had several subs that
    # were all true network terminals, it used every sub instead of one walked arm,
    # which by construction yields multi-part geometry. It was gated on divide
    # coverage below 0.9 and fired 0 times on all 24 VPUs: including Hawaii,
    # Puerto Rico and Alaska, the coastal domains it was written for. Such
    # flowpaths already took the single-arm path below, so this is a no-op on
    # current data, and reference flowlines are LINESTRING so a derived line
    # should not be multi-part by design.

    if (length(outlet_local) > 1L)
      outlet_local <- outlet_local[which.max(sub_tda[outlet_local])]

    up_by_target <- split(seq_along(sub_ids), sub_to)

    # Mainstem-strict walk: stop the walk at any sub whose levelpathid does
    # not match the parent's mainstem. Relaxing this to accept single-child
    # transitions into lp-mismatched subs does more harm than good
    # because the "single-child" pattern does not reliably distinguish
    # reconciliation renumbering (legit continuation) from cross-mainstem
    # aggregation (walking into the wrong stream).
    pick_next <- function(children) {
      if (length(children) == 1L) {
        c_lp <- sub_lp_local[children]
        lp_ok <- is.na(c_lp) || is.na(fp_lp) || c_lp == fp_lp
        return(if (lp_ok) children else integer(0))
      }
      lp_match <- is.na(sub_lp_local[children]) |
        is.na(fp_lp) |
        sub_lp_local[children] == fp_lp
      if (!any(lp_match)) return(integer(0))
      cand <- children[lp_match]
      cand[which.max(sub_tda[cand])]
    }

    path_local <- outlet_local
    current    <- outlet_local
    repeat {
      cid      <- sub_ids[current]
      children <- up_by_target[[cid]]
      if (is.null(children) || length(children) == 0L) break
      nxt <- pick_next(children)
      if (length(nxt) == 0L) { n_truncated <- n_truncated + 1L; break }
      path_local <- c(path_local, nxt)
      current    <- nxt
      if (length(path_local) > length(sub_ids)) break
    }
    # Single-sub rewrite, decided on length rather than divide coverage.
    #
    # When the walk ends with only one sub, that sub is usually shorter than the
    # aggregated original, so rewriting loses line. The old rule kept the original
    # only when its divide coverage was >= 0.9, which meant a flowpath whose line
    # legitimately left its catchment got rewritten to the shorter sub and lost
    # geometry. That is the same retired property Step 5d was removed for
    # enforcing: an aggregated flowpath's line is a mainstem walk through its
    # members and is not required to lie inside its own catchment.
    #
    # The rule that matters is simply never to shorten a line. Rewrite only when
    # the walked sub is at least as long as what is already there, so this can add
    # geometry but never remove it.
    if (length(path_local) == 0L) next
    if (length(path_local) == 1L) {
      if (is.finite(fp_len_m[i]) && fp_len_m[i] > 0 &&
        fl_len_m[rows[path_local]] < fp_len_m[i]) {
        n_single_skp <- n_single_skp + 1L
        next
      }
    }
    path_by_fp[[i]] <- rows[path_local]
  }

  # -- Phase 2. Bulk union + line_merge via GEOS ---------------------------
  # For each fp with a path, build a GEOS MultiLineString from the path subs
  # then run geos_unary_union + geos_line_merge in a single vectorized pass.
  has_path <- !vapply(path_by_fp, is.null, logical(1L))
  new_geoms <- as.list(sf::st_geometry(fp))
  extended  <- 0L

  if (any(has_path)) {
    fl_geos <- geos::as_geos_geometry(sf::st_geometry(fl))
    idx_with_path <- which(has_path)

    # Build a per-fp MultiLineString (as GEOS geometry collection) holding
    # that fp's path subs. geos_make_collection accepts a GEOS vector and a
    # grouping vector: assigning one group per fp lets it build N
    # collections in one call.
    flat_subs     <- unlist(path_by_fp[idx_with_path], use.names = FALSE)
    group_vec     <- rep(seq_along(idx_with_path),
      lengths(path_by_fp[idx_with_path]))
    sub_geos      <- fl_geos[flat_subs]
    collections   <- geos::geos_make_collection(
      sub_geos, type_id = "multilinestring", feature_id = group_vec)
    merged_all    <- tryCatch(
      geos::geos_line_merge(geos::geos_unary_union(collections)),
      error = function(e) NULL)

    if (!is.null(merged_all) && length(merged_all) == length(idx_with_path)) {
      merged_sfc <- sf::st_as_sfc(merged_all, crs = target_crs)
      for (k in seq_along(idx_with_path)) {
        g <- merged_sfc[[k]]
        if (!is.null(g)) {
          new_geoms[[idx_with_path[k]]] <- g
          extended <- extended + 1L
        }
      }
    }
  }

  new_sfc <- sf::st_sfc(new_geoms, crs = target_crs)
  fp_out  <- sf::st_sf(sf::st_drop_geometry(fp), geometry = new_sfc)
  if ("lengthkm" %in% names(fp_out)) {
    fp_out$lengthkm <- as.numeric(sf::st_length(
      sf::st_transform(fp_out, 5070))) / 1000
  }
  # Strip list-columns (gpkg can't serialize them)
  fp_out <- .drop_list_cols(fp_out)

  # Parts per feature, counted explicitly: on an sfc_LINESTRING lengths() returns
  # the coordinate count, so it cannot be used to detect multi-part geometry.
  .np <- vapply(new_geoms, function(g)
    if (inherits(g, "MULTILINESTRING")) length(g) else 1L, integer(1L))
  message(sprintf("extend_flowpaths_to_mainstem: extended %d of %d flowpath(s)",
    extended, total))
  message(sprintf(paste0(
    "extend_flowpaths_to_mainstem: %d multi-part result(s) (max %d parts) | ",
    "%d single-sub rewrite(s) skipped to avoid shortening, ",
    "%d walk(s) truncated off-mainstem"),
    sum(.np > 1L), max(.np), n_single_skp, n_truncated))

  sf::st_write(fp_out, input_gpkg, flowpath_layer,
    delete_layer = TRUE, quiet = TRUE)
  invisible(fp_out)
}
