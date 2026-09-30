# Internal pipeline repair + cleanup steps.
#
# These functions fix or tidy geometry/topology *within* the pipeline stages
# (Stage 1 reconcile, Stage 2 aggregate, Stage 3 nextgen). They are not public
# entry points, so they are documented `@keywords internal` (no longer exported).
# Merged here from per-function files; see git history for original provenance.

# recover_missing_divides was removed: it recovered nothing on any of the 24
# VPUs, before or after the divide-conservation fixes. It existed to rebuild
# divides that earlier steps had deleted; those steps now absorb into a
# neighbour instead, so there is nothing to put back. The flowpaths it reported
# as missing divides are connectors, which are divide-less by definition.

# =====================================================================
# prune_ghost_divides  (was R/prune_ghost_divides.R)
# =====================================================================
#' Make every network divide carry one of its flowpath's flowlines
#'
#' @description
#' NextGen requires a divide to be the catchment of the flowpath that owns it,
#' and the `divide_has_member_flowline` gate states that as: at least one of the
#' owning flowpath's flowlines must intersect the divide, or lie within 30 m of
#' it. This enforces exactly that rule, in the same terms, once the geometry is
#' final.
#'
#' It exists because a domain whose grid and vector network come from unrelated
#' sources cannot be split exactly. Alaska pairs TDX-Hydro reaches with
#' HydroSHEDS grids, and where the two disagree the partition can hand a flowpath
#' ground its own line never touches. Chasing that to zero in the splitter has a
#' floor set by the input data; conforming the finished fabric to the
#' specification does not.
#'
#' Two remedies, in order of preference:
#'
#' * **Re-own.** If the divide carries some *other* flowpath's flowlines and that
#'   flowpath has no divide, the pairing was simply swapped: move the divide to
#'   the flowpath whose line it actually carries. Nothing moves geometrically.
#' * **Absorb.** Otherwise the divide is nobody's catchment. Merge it into the
#'   neighbour sharing the longest boundary and mark its flowpath a connector.
#'   Area is conserved, the ground stays in the layer, it just changes owner.
#'
#' Only divides that FAIL the gate are touched, so a domain that already
#' conforms is left byte-identical. That is the safety property: this cannot
#' regress a passing VPU, because it does nothing on one.
#'
#' @param input_gpkg Path to the NextGen GeoPackage.
#' @param tol_m Distance below which a flowline counts as touching its divide.
#'   Defaults to the fabric's 30 m tolerance, matching the gate.
#' @param divide_layer,flowpath_layer,flowline_layer Layer names.
#' @return Invisibly, a data.frame of the divides acted on and what was done.
#' @keywords internal
normalize_divide_ownership <- function(input_gpkg,
                                       tol_m          = 30,
                                       divide_layer   = "divides",
                                       flowpath_layer = "flowpaths",
                                       flowline_layer = "flowlines") {
  if (!hfutils::layer_exists(input_gpkg, flowline_layer)) {
    message("normalize_divide_ownership: no flowlines layer; skipped")
    return(invisible(data.frame()))
  }
  d  <- sf::st_read(input_gpkg, divide_layer,   quiet = TRUE)
  fp <- sf::st_read(input_gpkg, flowpath_layer, quiet = TRUE)
  fl <- sf::st_read(input_gpkg, flowline_layer, quiet = TRUE)
  if (!all(c("flowpath_id", "divide_id") %in% names(d)) ||
      !("flowpath_id" %in% names(fl))) {
    message("normalize_divide_ownership: layers lack flowpath_id; skipped")
    return(invisible(data.frame()))
  }

  did   <- as.character(d$divide_id)
  owner <- as.character(d$flowpath_id)
  is_net <- !is.na(owner) & nzchar(owner) &
    !(!is.na(did) & startsWith(did, "lnd-"))
  if (!any(is_net)) return(invisible(data.frame()))

  fl <- sf::st_transform(fl, sf::st_crs(d))
  fl_owner <- as.character(fl$flowpath_id)
  fl_by <- split(seq_len(nrow(fl)), fl_owner)

  # Same test the gate applies: intersects, else within tol_m.
  holds <- function(i) {
    ix <- fl_by[[owner[i]]]
    if (is.null(ix) || !length(ix)) return(NA)      # untestable, leave alone
    g <- sf::st_geometry(d)[i]
    if (any(lengths(sf::st_intersects(sf::st_geometry(fl)[ix], g)) > 0L)) return(TRUE)
    suppressWarnings(any(as.numeric(
      sf::st_distance(sf::st_geometry(fl)[ix], g)) <= tol_m, na.rm = TRUE))
  }
  ok <- vapply(which(is_net), holds, logical(1))
  bad <- which(is_net)[!is.na(ok) & !ok]
  if (!length(bad)) {
    message("normalize_divide_ownership: all network divides carry a member flowline")
    return(invisible(data.frame()))
  }

  fp_id <- as.character(fp$flowpath_id)

  # Who does each failing divide actually belong to? The flowpath with the most
  # of its line inside the polygon. Computed for ALL failures first and applied
  # together, because the common case is a mutual swap: two divides each holding
  # the other's line. Resolving them one at a time cannot work, since neither
  # candidate is free at the moment it is considered.
  best <- rep(NA_character_, nrow(d)); best_len <- rep(0, nrow(d))
  for (i in bad) {
    g <- sf::st_geometry(d)[i]
    near <- which(lengths(sf::st_intersects(sf::st_geometry(fl), g)) > 0L)
    if (!length(near)) next
    inter <- suppressWarnings(sf::st_intersection(sf::st_geometry(fl)[near], g))
    len   <- suppressWarnings(as.numeric(sf::st_length(inter)))
    if (!length(len)) next
    by_fp <- tapply(len, fl_owner[near][seq_along(len)], sum, na.rm = TRUE)
    by_fp <- by_fp[is.finite(by_fp) & by_fp > 0]
    if (!length(by_fp)) next
    best[i]     <- names(by_fp)[which.max(by_fp)]
    best_len[i] <- max(by_fp)
  }

  # A flowpath owns at most one divide. Where two divides name the same flowpath,
  # the one carrying more of its line wins and the other is absorbed.
  proposed <- owner
  proposed[bad] <- best[bad]
  keep <- rep(TRUE, nrow(d))
  for (f in unique(stats::na.omit(proposed[bad]))) {
    ix <- which(proposed == f & !is.na(proposed))
    if (length(ix) > 1L) {
      win <- ix[which.max(ifelse(ix %in% bad, best_len[ix], Inf))]
      keep[setdiff(ix, win)] <- FALSE
    }
  }
  reowned  <- bad[!is.na(best[bad]) & keep[bad]]
  absorbed <- setdiff(bad, reowned)

  report <- rbind(
    if (length(reowned)) data.frame(divide_id = did[reowned], action = "reowned",
      from = owner[reowned], to = best[reowned], stringsAsFactors = FALSE),
    if (length(absorbed)) data.frame(divide_id = did[absorbed], action = "absorbed",
      from = owner[absorbed], to = NA_character_, stringsAsFactors = FALSE))

  if (length(reowned)) d$flowpath_id[reowned] <- best[reowned]
  if (length(absorbed)) {
    mask <- rep(FALSE, nrow(d)); mask[absorbed] <- TRUE
    d <- .absorb_orphan_divides(d, mask)$div
  }

  # Rebuild flowpath ownership from the divides layer rather than editing it in
  # step with the loop above: the divides are the authority on who owns what, and
  # deriving it once cannot leave the two disagreeing.
  own_now <- as.character(d$flowpath_id)
  own_now <- own_now[!is.na(own_now) & nzchar(own_now)]
  div_of  <- stats::setNames(as.character(d$divide_id)[
    !is.na(d$flowpath_id) & nzchar(as.character(d$flowpath_id))], own_now)
  newdiv <- unname(div_of[fp_id])
  fp$divide_id <- newdiv
  if ("has_divide" %in% names(fp)) fp$has_divide <- !is.na(newdiv)
  if ("areasqkm"  %in% names(fp)) fp$areasqkm[is.na(newdiv)] <- 0

  sf::st_write(.drop_list_cols(d),  input_gpkg, divide_layer,
    delete_layer = TRUE, quiet = TRUE)
  sf::st_write(.drop_list_cols(fp), input_gpkg, flowpath_layer,
    delete_layer = TRUE, quiet = TRUE)

  rep_df <- report
  message(sprintf(paste0(
    "normalize_divide_ownership: %d divide(s) did not carry a member flowline; ",
    "%d re-owned to the flowpath whose line they hold, %d absorbed into a ",
    "neighbour (owner demoted to connector)"),
    nrow(rep_df), sum(rep_df$action == "reowned"), sum(rep_df$action == "absorbed")))
  invisible(rep_df)
}

#' Prune ghost divides (areas below a threshold)
#'
#' @description
#' Removes network divides whose `areasqkm` is below a threshold. These are
#' usually produced by imbalanced PRD watershed splits where two sibling
#' reconciled flowlines share a parent reference COMID and one side receives
#' essentially no area (often <500 m2). The corresponding aggregated
#' flowpath is a through-reach with no meaningful local drainage; keeping a
#' microscopic "ghost" polygon makes the visualization misleading and the
#' schema inconsistent (has_divide = TRUE with a 100 m2 polygon).
#'
#' For each ghost divide: the divide row is dropped and its flowpath is
#' updated with `has_divide = FALSE`, `divide_id = NA`, `areasqkm = 0`. A
#' subsequent `hfutils::clean_gpkg_layer(..., gap_fill_area = > threshold)`
#' call will absorb the resulting tiny gaps into neighbor divides.
#'
#' @param input_gpkg Path to the GeoPackage (must contain `divides` and
#'   `flowpaths` layers).
#' @param threshold_sqkm Divides with `areasqkm` strictly less than this are
#'   pruned. Default `0.001` (1000 m2).
#' @param divide_layer Name of the divides layer. Default `"divides"`.
#' @param flowpath_layer Name of the flowpaths layer. Default `"flowpaths"`.
#'
#' @return Invisibly, a list of pruned `divide_id` values.
#'
#' @keywords internal
prune_ghost_divides <- function(input_gpkg,
                                threshold_sqkm  = 0.001,
                                divide_layer    = "divides",
                                flowpath_layer  = "flowpaths") {
  d  <- sf::st_read(input_gpkg, divide_layer,   quiet = TRUE)
  fp <- sf::st_read(input_gpkg, flowpath_layer, quiet = TRUE)

  is_network <- if ("type" %in% names(d))
    (is.na(d$type) | d$type == "network") else rep(TRUE, nrow(d))
  # Measure the geometry rather than trusting the stored `areasqkm`. This prune
  # runs wherever geometry has just been rewritten, while the column is refreshed
  # at one fixed point in the build, so the two are only in step by coincidence of
  # ordering. Deciding "too small to be a catchment" from the polygon itself makes
  # the call correct wherever it is made, and costs one st_area over a layer this
  # function already reads in full.
  area_km2 <- suppressWarnings(as.numeric(sf::st_area(d))) / 1e6
  if (all(is.na(area_km2))) area_km2 <- suppressWarnings(as.numeric(d$areasqkm))
  ghost_mask <- is_network & !is.na(area_km2) & area_km2 < threshold_sqkm
  if (!any(ghost_mask)) {
    message(sprintf("prune_ghost_divides: no divides below %.4f sqkm",
      threshold_sqkm))
    return(invisible(character(0)))
  }

  ghost_ids   <- as.character(d$divide_id[ghost_mask])
  ghost_fpids <- as.character(d$flowpath_id[ghost_mask])

  # A ghost divide is too small to be a catchment, but it is still land. Merge
  # it into a neighbour instead of deleting it: the flowpath still becomes a
  # connector (below), which is the point of this step, while the ground stays in
  # the layer. Deleting was bounded at threshold_sqkm per divide, but bounded is
  # not zero, and the coverage gate is meant to reach zero.
  ab    <- .absorb_orphan_divides(d, ghost_mask)
  d_out <- ab$div
  d_out <- .drop_list_cols(d_out)
  sf::st_write(d_out, input_gpkg, divide_layer,
    delete_layer = TRUE, quiet = TRUE)

  # Update flowpaths: mark has_divide=FALSE, clear divide_id/areasqkm
  fp_mask <- as.character(fp$flowpath_id) %in% ghost_fpids
  if (any(fp_mask)) {
    if ("has_divide" %in% names(fp)) fp$has_divide[fp_mask] <- FALSE
    if ("divide_id"  %in% names(fp)) fp$divide_id[fp_mask]  <- NA_character_
    if ("areasqkm"   %in% names(fp)) fp$areasqkm[fp_mask]   <- 0
    fp <- .drop_list_cols(fp)
    sf::st_write(fp, input_gpkg, flowpath_layer,
      delete_layer = TRUE, quiet = TRUE)
  }

  message(sprintf(paste0("prune_ghost_divides: %d divide(s) below %.4f sqkm; ",
    "absorbed %d into neighbours, kept %d (area preserved)"),
    length(ghost_ids), threshold_sqkm, ab$n_absorbed, ab$n_kept))
  invisible(ghost_ids)
}

# =====================================================================
# resolve_disjoint_flowpath_divides  (was R/resolve_disjoint_flowpath_divides.R)
# =====================================================================
#' Clear disjoint flowpath-divide assignments (in-memory)
#'
#' @description
#' Core detection and clearing logic for flowpath-divide misassignments.
#' Operates on in-memory sf objects; used by \code{resolve_disjoint_flowpath_divides}
#' and by the aggregate phase before aggregation.
#'
#' @param fp  sf with flowpath_id, divide_id, and line geometry.
#' @param div sf with divide_id and polygon geometry.
#' @param drop_orphan_cats If TRUE, remove divides no longer referenced.
#' @return list(fp = sf, div = sf, report = data.frame of cleared pairs)
#' @keywords internal
.clear_disjoint_fp_divides <- function(fp, div, drop_orphan_cats = TRUE,
                                       flowlines = NULL) {

  target_crs <- sf::st_crs(fp)
  if (sf::st_crs(div) != target_crs) div <- sf::st_transform(div, target_crs)

  fp  <- sf::st_make_valid(fp)
  div <- sf::st_make_valid(div)

  fp$flowpath_id  <- as.character(fp$flowpath_id)
  if ("flowpath_id" %in% names(div)) div$flowpath_id <- as.character(div$flowpath_id)
  if ("divide_id"   %in% names(div)) div$divide_id   <- as.character(div$divide_id)

  # Join fp to div via flowpath_id (the common key at stage 2).
  # fp may also carry divide_id directly (stage 3 gpkg path); prefer that when
  # present to stay consistent with the original function semantics.
  if ("divide_id" %in% names(fp)) {
    fp$divide_id <- as.character(fp$divide_id)
    d_ix  <- match(fp$divide_id, div$divide_id)
    valid <- !is.na(d_ix) & !is.na(fp$divide_id) & nchar(fp$divide_id) > 0L
    get_div_id <- function(i) fp$divide_id[i]
  } else {
    d_ix  <- match(fp$flowpath_id, div$flowpath_id)
    valid <- !is.na(d_ix)
    get_div_id <- function(i) div$divide_id[d_ix[i]]
  }

  if (!any(valid)) return(list(fp = fp, div = div, report = data.frame()))

  fp_t    <- sf::st_transform(fp, sf::st_crs(div))
  lines_g <- geos::as_geos_geometry(sf::st_geometry(fp_t)[valid])
  divs_g  <- geos::as_geos_geometry(sf::st_geometry(div)[d_ix[valid]])
  inter_g <- geos::geos_intersection(lines_g, divs_g)
  lt  <- suppressWarnings(as.numeric(geos::geos_length(lines_g)))
  li  <- suppressWarnings(as.numeric(geos::geos_length(inter_g)))
  cov <- ifelse(lt > 0, li / lt, NA_real_)

  # For zero-coverage candidates, compute distance from the fp line to its divide.
  # A fp that terminates exactly ON the divide boundary has cov=0 but dist=0
  # (the endpoint touches the boundary). This is expected for POI-split upstream
  # sub-segments whose outlet is at the POI location on the divide edge.
  # Only flag as disjoint when the fp is genuinely AWAY from the divide (dist > 1m).
  dists_all <- tryCatch(
    suppressWarnings(as.numeric(geos::geos_distance(lines_g, divs_g))),
    error = function(e) rep(Inf, length(cov)))
  # Ask the question the gate asks before taking a divide away. This pass judged
  # the AGGREGATED line: zero overlap with its own polygon means misassigned.
  # But an aggregated line is a mainstem walk through its members and is not
  # required to lie inside the catchment: build_nextgen states exactly that,
  # and declines to clip for the same reason. The property that actually has to
  # hold is divide_has_member_flowline: at least one of the flowpath's MEMBER
  # flowlines within 30 m. Measured on VPU 09, all 5 divides this pass cleared
  # satisfied that, so every one was absorbed away from a flowpath the shipped
  # contract considers correctly paired.
  #
  # Consulting the members costs nothing here and makes the repair and the gate
  # agree by construction. Without a flowlines layer the behaviour is unchanged.
  member_holds <- rep(FALSE, length(cov))
  if (!is.null(flowlines) && inherits(flowlines, "sf") &&
      "flowpath_id" %in% names(flowlines) && nrow(flowlines)) {
    flt   <- sf::st_transform(flowlines, sf::st_crs(div))
    fl_by <- split(seq_len(nrow(flt)), as.character(flt$flowpath_id))
    vix   <- which(valid)
    member_holds <- vapply(seq_along(vix), function(k) {
      ix <- fl_by[[as.character(fp$flowpath_id[vix[k]])]]
      if (is.null(ix) || !length(ix)) return(FALSE)
      gd <- sf::st_geometry(div)[d_ix[vix[k]]]
      gl <- sf::st_geometry(flt)[ix]
      if (any(lengths(sf::st_intersects(gl, gd)) > 0L)) return(TRUE)
      suppressWarnings(any(as.numeric(sf::st_distance(gl, gd)) <= 30, na.rm = TRUE))
    }, logical(1L))
  }

  clear <- which(!is.na(cov) & cov == 0 & dists_all > 1.0 & !member_holds)
  if (length(clear) == 0L)
    return(list(fp = fp, div = div, report = data.frame()))

  dists <- dists_all[clear]
  fp_row_ix <- which(valid)[clear]
  cleared_div_ids <- vapply(fp_row_ix, get_div_id, character(1L))

  # -- Phase 1 metrics: classify Group A (split-flowline) vs Group B -----
  is_groupA   <- grepl("\\.", fp$flowpath_id[fp_row_ix])
  n_groupA    <- sum(is_groupA)
  n_groupB    <- sum(!is_groupA)

  report    <- data.frame(
    flowpath_id   = fp$flowpath_id[fp_row_ix],
    old_divide_id = cleared_div_ids,
    distance_m    = round(dists),
    group         = ifelse(is_groupA, "A", "B"),
    stringsAsFactors = FALSE)

  # -- Phase 2: overlap-driven divide reassignment (Group B swap) ---------
  # For each disjoint fp, compute intersection length against ALL divides to
  # find the polygon that actually spatially corresponds to this flowpath.
  # Group B (FDR raster headwater misattributions): the fp line is correct but
  # its delineated divide landed in the wrong location. The correct divide is
  # the one whose polygon contains or overlaps the fp. We find it via geos
  # intersection, swap the assignment, and let the old (wrong) divide become
  # an orphan to be cleaned up below.
  #
  # The target divide must be unowned, on both counts: not claimed by an earlier
  # swap in this pass (`swapped_to`), and not already held by a flowpath outside
  # this repair (`held_elsewhere`). The polygon a misplaced flowpath overlaps is
  # usually its neighbour's correctly assigned divide, so without the second
  # check the swap leaves two flowpaths pointing at one divide. Divides held by
  # the flowpaths under repair are fair game: those assignments are the ones
  # being vacated.
  swapped_to   <- character(0L)     # divide_ids claimed in this pass
  swapped_rows <- integer(0L)       # fp_row_ix entries that got a swap
  n_swapped    <- 0L

  if ("divide_id" %in% names(fp)) {
    .under_repair <- rep(FALSE, nrow(fp))
    .under_repair[fp_row_ix] <- TRUE
    held_elsewhere <- unique(as.character(fp$divide_id[!.under_repair]))
    held_elsewhere <- held_elsewhere[!is.na(held_elsewhere) & nzchar(held_elsewhere)]

    all_divs_g <- geos::as_geos_geometry(sf::st_geometry(div))
    fp_t_swap  <- sf::st_transform(fp, sf::st_crs(div))

    for (k in seq_along(fp_row_ix)) {
      ri <- fp_row_ix[k]

      fp_g <- tryCatch(
        geos::as_geos_geometry(sf::st_geometry(fp_t_swap)[ri]),
        error = function(e) NULL)
      if (is.null(fp_g)) {
        # Can't intersect: clear assignment
        fp$divide_id[ri] <- NA_character_
        if ("has_divide" %in% names(fp)) fp$has_divide[ri] <- FALSE
        next
      }

      inter_lens <- tryCatch(
        suppressWarnings(as.numeric(
          geos::geos_length(geos::geos_intersection(fp_g, all_divs_g)))),
        error = function(e) rep(0, nrow(div)))
      inter_lens[is.na(inter_lens)] <- 0

      # Exclude every divide that already has an owner: claimed by a prior swap
      # in this pass, or held by a flowpath outside the repair set. Without the
      # second exclusion the swap manufactures a shared divide.
      .taken <- c(swapped_to, held_elsewhere)
      if (length(.taken) > 0L)
        inter_lens[as.character(div$divide_id) %in% .taken] <- 0

      best_idx <- which.max(inter_lens)
      best_len <- if (length(best_idx) > 0L) inter_lens[best_idx] else 0

      # Threshold: 0.5 m: eliminates floating-point boundary-touch noise
      # while accepting tiny headwater fps (Group B are raster-cell scale).
      if (best_len >= 0.5) {
        new_did <- as.character(div$divide_id[best_idx])
        old_did <- cleared_div_ids[k]
        fp$divide_id[ri] <- new_did
        if ("has_divide" %in% names(fp)) fp$has_divide[ri] <- TRUE
        swapped_to   <- c(swapped_to, new_did)
        swapped_rows <- c(swapped_rows, k)
        n_swapped    <- n_swapped + 1L
        report$group[k] <- paste0(report$group[k], "_swap")
      } else {
        # No overlapping divide found: clear the stale assignment.
        # Log the best overlap seen (diagnostic: helps tune threshold).
        if (best_len > 0 && best_len < 0.5)
          message(sprintf(
            "resolve_disjoint_flowpath_divides: fp %s best overlap=%.4f m (below 0.5m threshold; cleared)",
            fp$flowpath_id[ri], best_len))
        fp$divide_id[ri] <- NA_character_
        if ("has_divide" %in% names(fp)) fp$has_divide[ri] <- FALSE
      }
    }
  } else {
    # -- Phase 2 (else-branch): overlap-driven swap for Group B -----------
    # At stage2, fp has no divide_id column; join is via
    # div$flowpath_id == fp$flowpath_id. The correct divide for a disjoint
    # fp is the one whose polygon contains (or best overlaps) the fp line.
    # We redirect that divide's flowpath_id to the disjoint fp, and track
    # both the old wrong div AND unclaimed disjoint fp divs as orphans by
    # row-index (not fp_id) to handle duplicate-id ambiguity after swap.
    all_divs_g   <- geos::as_geos_geometry(sf::st_geometry(div))
    fp_t_swap    <- sf::st_transform(fp, sf::st_crs(div))
    claimed_rows <- integer(0L)  # div row indices already claimed by a swap
    swap_targets <- rep(NA_integer_, length(fp_row_ix))

    # Protect divides already correctly assigned to non-disjoint fps.
    # A non-disjoint fp has demonstrated coverage > 0, so its divide should not
    # be "stolen" by a disjoint fp swap. Stealing creates a new disjoint pair for
    # the previously-correct fp (e.g. giving a split-COMID's watershed polygon to
    # a sibling transit fp while leaving the correct fp without a divide).
    non_disjoint_fp_ix <- setdiff(which(valid), fp_row_ix)
    protected_rows     <- d_ix[non_disjoint_fp_ix]
    protected_rows     <- protected_rows[!is.na(protected_rows)]

    for (k in seq_along(fp_row_ix)) {
      ri <- fp_row_ix[k]
      fp_g <- tryCatch(
        geos::as_geos_geometry(sf::st_geometry(fp_t_swap)[ri]),
        error = function(e) NULL)
      if (is.null(fp_g)) next

      inter_lens <- tryCatch(
        suppressWarnings(as.numeric(
          geos::geos_length(geos::geos_intersection(fp_g, all_divs_g)))),
        error = function(e) rep(0, nrow(div)))
      inter_lens[is.na(inter_lens)] <- 0

      # Exclude divides already claimed by a prior swap.
      if (length(claimed_rows) > 0L) inter_lens[claimed_rows] <- 0
      # Exclude the fp's own (wrong) div to prevent self-swap.
      own_row <- d_ix[fp_row_ix[k]]
      if (!is.na(own_row)) inter_lens[own_row] <- 0
      # Exclude divides already correctly assigned to non-disjoint fps.
      if (length(protected_rows) > 0L) inter_lens[protected_rows] <- 0

      best_idx <- which.max(inter_lens)
      best_len <- if (length(best_idx) > 0L) inter_lens[best_idx] else 0

      if (best_len >= 0.5) {
        swap_targets[k] <- best_idx
        claimed_rows    <- c(claimed_rows, best_idx)
      } else if (best_len > 0) {
        message(sprintf(
          "resolve_disjoint_flowpath_divides: fp %s best overlap=%.4f m (below 0.5 m threshold; cleared)",
          fp$flowpath_id[ri], best_len))
      }
    }

    # Apply swaps and collect orphan row indices.
    orphan_rows <- integer(0L)
    for (k in seq_along(fp_row_ix)) {
      old_row <- d_ix[fp_row_ix[k]]   # row of the (wrong) div currently assigned to this fp
      if (!is.na(swap_targets[k])) {
        div$flowpath_id[swap_targets[k]] <- fp$flowpath_id[fp_row_ix[k]]
        if ("divide_id" %in% names(div))
          div$divide_id[swap_targets[k]] <- fp$flowpath_id[fp_row_ix[k]]
        if (!is.na(old_row)) orphan_rows <- c(orphan_rows, old_row)
        n_swapped    <- n_swapped + 1L
        swapped_rows <- c(swapped_rows, k)
        report$group[k] <- paste0(report$group[k], "_swap")
      } else {
        if (!is.na(old_row)) orphan_rows <- c(orphan_rows, old_row)
      }
    }
    orphan_rows <- unique(orphan_rows)
  }

  if (n_swapped > 0L)
    message(sprintf(
      "resolve_disjoint_flowpath_divides: swapped %d divide(s) to overlapping polygon (Group B fix)",
      n_swapped))
  if (n_groupA > 0L)
    message(sprintf(
      "resolve_disjoint_flowpath_divides: %d Group A (split-flowline) pair(s) cleared -- fix via stage1 progressive-focal snap",
      n_groupA))

  if (isTRUE(drop_orphan_cats)) {
    # Old (wrong) divides whose fp has been swapped away or cleared are now
    # orphaned. Drop them, they are genuinely unowned at this point.
    if ("divide_id" %in% names(fp)) {
      still_ref    <- unique(fp$divide_id[!is.na(fp$divide_id)])
      orphans_mask <- !div$divide_id %in% still_ref & div$divide_id %in% cleared_div_ids
    } else {
      # Use row-index based orphan tracking (set by Phase 2 swap above)
      # to avoid ambiguity when two divides share a flowpath_id after swap.
      if (exists("orphan_rows") && length(orphan_rows) > 0L) {
        orphans_mask            <- logical(nrow(div))
        orphans_mask[orphan_rows] <- TRUE
      } else {
        cleared_fp_ids <- fp$flowpath_id[fp_row_ix]
        orphans_mask   <- div$flowpath_id %in% cleared_fp_ids
      }
    }
    if (any(orphans_mask)) {
      # These divides are unowned, not empty. Deleting them deleted their land:
      # 97 km2 on one CONUS VPU, hidden because divides_area_conserved compares
      # sums and cannot tell a hole from an overlap removal. Merge each into the
      # neighbour it shares the most boundary with, which conserves coverage and
      # still leaves one divide per flowpath.
      n_orphan <- sum(orphans_mask)
      ab  <- .absorb_orphan_divides(div, orphans_mask)
      div <- ab$div
      message(sprintf(paste0("resolve_disjoint_flowpath_divides: absorbed %d of %d ",
        "orphan divide(s) into their neighbours%s"),
        ab$n_absorbed, n_orphan,
        if (ab$n_kept > 0L)
          sprintf("; kept %d with no neighbour (area preserved)", ab$n_kept) else ""))
    }
  }

  list(fp = fp, div = div, report = report)
}


#' Detect and clear flowpath-divide misassignments
#'
#' @description
#' Finds flowpaths whose line geometry is spatially disjoint from the divide
#' they're assigned to (coverage = 0, distance above `distance_threshold_m`)
#' and clears the bad assignment. The flowpath remains in the network
#' (topology is preserved via attribute IDs) but loses its catchment.
#'
#' Root cause is an aggregation join artifact that appears at VPU scale:
#' reconciled flowlines with missing/default levelpathid (=1) can be grouped
#' into aggregate flowpaths whose member subs are geographically unrelated,
#' leading to a bogus cat pairing. The proper fix is to call this check
#' before aggregation (aggregate step 3c) so the bad pair never enters
#' aggregate_to_distribution.
#'
#' @param input_gpkg Path to the GeoPackage (mutated in place).
#' @param flowpath_layer,divides_layer Layer names.
#' @param distance_threshold_m Kept for back-compat; no longer used as a gate
#'   (coverage = 0 is the sole criterion). Reported in output for triage.
#' @param drop_orphan_cats If TRUE, remove the now-orphan divide polygons
#'   that no flowpath references. Default TRUE.
#'
#' @return Invisibly, a data.frame of cleared `{flowpath_id, old_divide_id,
#'   distance_m}` rows.
#' @keywords internal
resolve_disjoint_flowpath_divides <- function(input_gpkg,
                                              flowpath_layer = "flowpaths",
                                              divides_layer  = "divides",
                                              distance_threshold_m = 100,
                                              drop_orphan_cats = TRUE) {
  fp  <- sf::st_read(input_gpkg, flowpath_layer, quiet = TRUE)
  div <- sf::st_read(input_gpkg, divides_layer,  quiet = TRUE)
  div_pre <- div                      # rows this pass may delete; kept for the probe

  fl_layer <- if ("flowlines" %in% sf::st_layers(input_gpkg)$name) {
    tryCatch(sf::st_read(input_gpkg, "flowlines", quiet = TRUE), error = function(e) NULL)
  } else NULL

  result <- .clear_disjoint_fp_divides(fp, div, drop_orphan_cats, flowlines = fl_layer)

  if (nrow(result$report) == 0L) {
    message("resolve_disjoint_flowpath_divides: no disjoint fps found")
    return(invisible(data.frame()))
  }

  fp  <- result$fp
  div <- result$div

  fp  <- .drop_list_cols(fp)
  div <- .drop_list_cols(div)

  sf::st_write(fp,  input_gpkg, flowpath_layer,
    delete_layer = TRUE, quiet = TRUE)
  sf::st_write(div, input_gpkg, divides_layer,
    delete_layer = TRUE, quiet = TRUE)

  # Diagnostic (HF_PROBE_5B2=1): this pass judges the AGGREGATED line at 0 m
  # overlap, while divide_has_member_flowline judges the flowpath's MEMBER
  # flowlines within 30 m, and build_nextgen states that an aggregated line is
  # not required to lie inside its own catchment. Where the two disagree, a
  # divide is absorbed that the gate would have accepted. Count that here so the
  # criterion can be changed on evidence rather than on the reading of it.
  if (nzchar(Sys.getenv("HF_PROBE_5B2")) &&
      "flowlines" %in% sf::st_layers(input_gpkg)$name) {
    try({
      fl <- sf::st_read(input_gpkg, "flowlines", quiet = TRUE)
      if ("flowpath_id" %in% names(fl)) {
        d0  <- sf::st_transform(div_pre, sf::st_crs(fl))   # BEFORE any row was deleted
        cleared <- as.character(result$report$flowpath_id)
        oldid   <- as.character(result$report$old_divide_id)
        flsp <- split(seq_len(nrow(fl)), as.character(fl$flowpath_id))
        dix  <- match(oldid, as.character(d0$divide_id))
        no_fl  <- vapply(cleared, function(x) is.null(flsp[[x]]), logical(1L))
        no_div <- is.na(dix)
        keepable <- vapply(seq_along(cleared), function(i) {
          if (no_fl[i] || no_div[i]) return(NA)
          any(lengths(sf::st_intersects(
            sf::st_geometry(fl)[flsp[[cleared[i]]]],
            sf::st_buffer(sf::st_geometry(d0)[dix[i]], 30))) > 0L)
        }, logical(1L))
        message(sprintf(paste0("resolve_disjoint_flowpath_divides [probe]: of %d cleared, ",
          "%d had member flowlines within 30 m of the divide taken away ",
          "(the gate would have accepted them); %d genuinely away; ",
          "%d had no flowlines; %d had no pre-repair divide row"),
          length(cleared), sum(keepable, na.rm = TRUE), sum(!keepable, na.rm = TRUE),
          sum(no_fl), sum(no_div & !no_fl)))
      }
    }, silent = TRUE)
  }

  message(sprintf(
    "resolve_disjoint_flowpath_divides: cleared %d bogus fp-cat assignment(s)",
    nrow(result$report)))

  invisible(result$report)
}

# =====================================================================
# restore_aggregated_divides  (was R/restore_aggregated_divides.R)
# =====================================================================
#' Restore aggregated divide areas dropped by PRD / aggregation
#'
#' @description
#' Two kinds of area-loss can affect an aggregated divide between the
#' reference NHD input and the final output:
#'
#' 1. **Aggregation-stage loss**: the aggregate-to-distribution step
#'    internally runs `hfutils::clean_geometry` on the aggregated divides,
#'    whose small-part glue logic can crush a reconciled-split-sibling
#'    polygon (e.g. the `.2` half of an reference COMID) to a ~100 m2 residual.
#'    Detected by comparing each aggregated cat's area to the **sum of its
#'    member reconciled-divide areas** (from `rec_div`).
#'
#' 2. **PRD-stage loss**: the refactor-split may drop one of the split
#'    children entirely (typically a narrow sliver that falls below a
#'    cleanup threshold inside `reconcile_divides`), so its area never
#'    appears in `rec_div` at all. The sibling's divide survives; the lost
#'    half vanishes into no catchment. Detected only when `reference_divides`
#'    is supplied, by comparing the aggregated cat's area to the **sum of
#'    reference NHD divide areas for its member COMIDs**.
#'
#' Whichever reference is larger (rec_div-based or reference-based) wins.
#' Shortfalls beyond `threshold` trigger a rebuild from the union of:
#' - rec_div geometries for members whose reconciled divide exists, plus
#' - reference-divide geometries for members missing from rec_div (the PRD
#'   dropouts). Divide polygons are not otherwise modified.
#'
#' @param input_gpkg Path to the intermediate GeoPackage (contains `divides`
#'   and `network` layers).
#' @param rec_div sf object of reconciled divides returned by
#'   `reconcile_divides()` (before aggregation). Must contain `flowpath_id`
#'   (= reconciled integer ID) and geometry.
#' @param reference_divides Optional sf of reference NHD divides (the pipeline's
#'   pre-refactor input). Must contain `divide_id` (the reference COMID as
#'   character) and geometry. When supplied, enables detection of PRD-stage
#'   drops where a split child was never produced.
#' @param threshold Fraction. Aggregated divides with area < `threshold * max(
#'   expected_rec, expected_ref)` are restored. Default `0.95`.
#' @param abs_shortfall_m2 numeric (m^2). A divide is also restored when its
#'   absolute area shortfall vs expected exceeds this, catching small-percentage
#'   gaps on large divides. Default `1e4`.
#' @param split_gap_m2 numeric (m^2). Minimum reference-divide gap area that
#'   qualifies as a genuine PRD split dropout worth filling. Default `5e5`.
#' @param divide_layer Name of the divides layer. Default `"divides"`.
#' @param network_layer Name of the network layer. Default `"network"`.
#' @param verbose If TRUE, emit a per-divide breakdown of area shortfalls and
#'   recovery deltas, plus a source-breakdown (rec_div / ref / both) per fix.
#'
#' @return Invisibly, character vector of restored aggregated divide ids.
#'
#' @keywords internal
restore_aggregated_divides <- function(input_gpkg, rec_div,
                                       reference_divides = NULL,
                                       threshold        = 0.9,
                                       abs_shortfall_m2 = 1e4,
                                       split_gap_m2     = 5e5,
                                       divide_layer     = "divides",
                                       network_layer    = "network",
                                       verbose          = FALSE) {
  d  <- sf::st_read(input_gpkg, divide_layer,  quiet = TRUE)
  nw <- sf::st_drop_geometry(
    sf::st_read(input_gpkg, network_layer, quiet = TRUE))

  rec_areas <- stats::setNames(as.numeric(sf::st_area(rec_div)),
    as.character(rec_div$flowpath_id))
  target_crs <- sf::st_crs(d)
  rec_crs    <- sf::st_crs(rec_div)

  # (A) rec-based expected: divide_id -> reconciled members via network.rec_id_list
  nw_rec <- unique(data.frame(
    divide_id    = as.character(nw$divide_id),
    rec_id_list  = as.character(nw$rec_id_list),
    stringsAsFactors = FALSE))
  nw_rec <- nw_rec[!is.na(nw_rec$divide_id) & !is.na(nw_rec$rec_id_list), ]
  splits <- strsplit(nw_rec$rec_id_list, ",", fixed = TRUE)
  long_rec <- unique(data.frame(
    divide_id = rep(nw_rec$divide_id, lengths(splits)),
    rec_id    = trimws(unlist(splits)),
    stringsAsFactors = FALSE))
  long_rec <- long_rec[long_rec$rec_id %in% names(rec_areas), , drop = FALSE]
  long_rec$rec_area <- unname(rec_areas[long_rec$rec_id])
  expected_rec_by_did <- tapply(long_rec$rec_area, long_rec$divide_id,
    sum, na.rm = TRUE)
  recs_by_did <- split(long_rec$rec_id, long_rec$divide_id)
  rec_idx <- stats::setNames(seq_len(nrow(rec_div)),
    as.character(rec_div$flowpath_id))

  # (B) reference-based expected: divide_id -> reference COMIDs via network.hf_id
  expected_ref_by_did <- NULL
  nhds_by_did         <- list()
  ref_areas           <- NULL
  ref_idx             <- NULL
  if (!is.null(reference_divides) && nrow(reference_divides) > 0L &&
    (any(c("reference_id", "hf_id") %in% names(nw)))) {
    .rid_col <- if ("reference_id" %in% names(nw)) "reference_id" else "hf_id"
    nw$hf_id_chr <- format(as.numeric(nw[[.rid_col]]),
      scientific = FALSE, trim = TRUE)
    nw_ref <- unique(data.frame(
      divide_id = as.character(nw$divide_id),
      nhd_comid = nw$hf_id_chr,
      stringsAsFactors = FALSE))
    nw_ref <- nw_ref[!is.na(nw_ref$divide_id) & !is.na(nw_ref$nhd_comid) &
      nzchar(nw_ref$nhd_comid), ]

    ref_areas <- stats::setNames(as.numeric(sf::st_area(reference_divides)),
      as.character(reference_divides$divide_id))
    ref_idx   <- stats::setNames(seq_len(nrow(reference_divides)),
      as.character(reference_divides$divide_id))
    nw_ref <- nw_ref[nw_ref$nhd_comid %in% names(ref_areas), , drop = FALSE]
    nw_ref$ref_area <- unname(ref_areas[nw_ref$nhd_comid])
    expected_ref_by_did <- tapply(nw_ref$ref_area, nw_ref$divide_id,
      sum, na.rm = TRUE)
    nhds_by_did <- split(nw_ref$nhd_comid, nw_ref$divide_id)
  }

  did_all <- as.character(d$divide_id)
  d_area_m2 <- as.numeric(sf::st_area(d))
  d_geoms   <- as.list(sf::st_geometry(d))

  expected_rec <- unname(expected_rec_by_did[did_all])
  expected_rec[is.na(expected_rec)] <- 0
  expected_ref <- if (!is.null(expected_ref_by_did))
    unname(expected_ref_by_did[did_all]) else rep(0, length(did_all))
  expected_ref[is.na(expected_ref)] <- 0

  # Aggregate-level detection against rec_div expected area only (NOT ref,
  # reference NHD polygons have NHD-digitization boundaries that often
  # exceed their PRD-clipped rec_div counterparts, creating systematic
  # pmax-based false positives that break cat shape alignment).
  # Two triggers: (a) cat below `threshold` fraction of rec expected, or
  # (b) absolute shortfall exceeds `abs_shortfall_m2` (catches small-percentage
  # drops that are hydrologically meaningful).
  agg_needs <- is.finite(expected_rec) & expected_rec > 0 &
    (d_area_m2 < threshold * expected_rec |
      expected_rec - d_area_m2 > abs_shortfall_m2)

  # (B) Split-member / PRD-dropout detection: cats whose flowpath member_comid
  # lists a `.N` split sibling AND whose reference-based total meaningfully
  # exceeds their rec-based total (> 1000 m2 gap) indicate PRD lost one of
  # the split halves before it reached rec_div. The aggregate-level check
  # cannot see this because the cat's current area matches `expected_rec`
  # exactly: both are missing the dropped half. Reference divides are the
  # only ground truth that detects this.
  split_needs <- rep(FALSE, length(did_all))
  if (!is.null(reference_divides) && hfutils::layer_exists(input_gpkg, "flowpaths")) {
    fp <- sf::st_drop_geometry(sf::st_read(input_gpkg, "flowpaths",
      quiet = TRUE))
    fp_did <- if ("divide_id" %in% names(fp)) as.character(fp$divide_id)
    else as.character(fp$flowpath_id)
    has_split <- grepl("\\.[0-9]", as.character(fp$member_comid))
    split_did <- fp_did[has_split]
    # Combined relative + absolute threshold: a PRD-stage drop worth the
    # shape-shift risk of a ref-based rebuild must be both (a) a meaningful
    # fraction of the cat's reference total AND (b) large in absolute terms.
    # This isolates genuine "missing half a sqkm" cases from boundary noise.
    # A split COMID is a TRUE dropout only when the PARENT reference catchment
    # is not fully covered by the union of ALL its split children's rec_div
    # pieces. Detect dropout at the parent level (comparing a single child's rec
    # area to the full parent reference would misfire on a correctly partitioned
    # split, rebuilding the child to the whole parent polygon and engulfing its
    # siblings): ref_area(P) - sum(rec_div pieces with orig_divide_id == P).
    if ("orig_divide_id" %in% names(rec_div)) {
      rec_cover_by_parent <- tapply(as.numeric(sf::st_area(rec_div)),
        as.character(rec_div$orig_divide_id),
        sum, na.rm = TRUE)
      genuine_ref_gap <- vapply(did_all, function(did) {
        cs <- nhds_by_did[[did]]
        if (is.null(cs) || !length(cs)) return(0)
        ra <- unname(ref_areas[cs])
        ra[is.na(ra)] <- 0
        rc <- unname(rec_cover_by_parent[cs])
        rc[is.na(rc)] <- 0
        sum(pmax(0, ra - rc), na.rm = TRUE)
      }, numeric(1))
    } else {
      genuine_ref_gap <- pmax(0, expected_ref - expected_rec)
    }
    ref_gap_frac <- ifelse(expected_ref > 0, genuine_ref_gap / expected_ref, 0)
    split_needs <- did_all %in% split_did &
      ref_gap_frac > 0.30 &
      genuine_ref_gap > split_gap_m2
  }

  # Only use reference-based rebuild when split_needs fires (cat has .N
  # member AND ref > rec meaningfully). Broadening to all ref>rec cases
  # causes shape-shift failures because reference NHD polygons have
  # different boundaries than PRD-clipped rec_div polygons, and the
  # flowpath geometry is aligned to the latter.
  use_ref_for <- split_needs
  needs_restore <- which(agg_needs | split_needs)

  ref_crs <- if (!is.null(reference_divides)) sf::st_crs(reference_divides) else NA

  restored    <- character(0)
  restore_log <- vector("list", 0L)
  # Live geometry vector for the occupancy test below, kept in step with d_geoms
  # so a divide restored earlier in this loop is treated as occupied by a later
  # one. Built once: rebuilding it per candidate would be O(candidates x divides).
  d_sfc     <- sf::st_sfc(d_geoms, crs = target_crs)
  # Bounding boxes for the occupancy test below. Intersecting one new geometry
  # against the whole divide set indexes the single geometry and then walks
  # every divide, so the cost is O(divides) per restore and the step scales with
  # the product of the two. Box overlap is a strict superset of intersection, so
  # screening on it first cannot change which divides are found, only how many
  # reach GEOS. The boxes track d_geoms, since a divide restored earlier in the
  # loop counts as occupied for a later one.
  .bx <- function(g) {
    b <- sf::st_bbox(sf::st_sfc(g, crs = target_crs))
    c(b[["xmin"]], b[["ymin"]], b[["xmax"]], b[["ymax"]])
  }
  bb <- vapply(d_geoms, .bx, numeric(4))          # 4 x n: xmin ymin xmax ymax
  n_trimmed <- 0L
  for (i in needs_restore) {
    did  <- did_all[i]
    recs <- recs_by_did[[did]]

    # Collect rec_div polygons for this divide's reconciled members.
    rec_geoms <- if (length(recs) > 0L) {
      sf::st_geometry(rec_div)[rec_idx[recs]]
    } else NULL

    # Reference-divide augmentation: only for cats flagged as split-dropout.
    # (For agg-level shortfalls, rec_div is the authoritative source; mixing
    # ref_div there would over-reclaim from neighbors.)
    ref_geoms <- NULL
    ref_used  <- character(0)
    if (use_ref_for[i] && did %in% names(nhds_by_did)) {
      ref_used  <- nhds_by_did[[did]]
      if (length(ref_used) > 0L) {
        ref_geoms <- sf::st_geometry(reference_divides)[ref_idx[ref_used]]
        if (!is.na(ref_crs) && ref_crs != target_crs)
          ref_geoms <- sf::st_transform(ref_geoms, target_crs)
      }
    }

    # Union: preserve the existing cat polygon and add the reconciled-member
    # polygons (which recover aggregate-stage losses) plus reference-divide
    # polygons (which recover PRD-stage losses). Using the current cat as the
    # base keeps existing boundaries intact: we only ADD area, never drop it.
    parts_list <- list(sf::st_sfc(d_geoms[[i]], crs = target_crs))
    if (!is.null(rec_geoms) && length(rec_geoms) > 0L) {
      if (rec_crs != target_crs)
        rec_geoms <- sf::st_transform(rec_geoms, target_crs)
      parts_list[[length(parts_list) + 1L]] <- rec_geoms
    }
    if (!is.null(ref_geoms) && length(ref_geoms) > 0L) {
      parts_list[[length(parts_list) + 1L]] <- ref_geoms
    }
    parts <- do.call(c, parts_list)
    if (length(parts) == 0L) next

    new_geom <- tryCatch(sf::st_union(parts), error = function(e) NULL)
    if (is.null(new_geom) || length(new_geom) == 0L) next

    # Normalize to POLYGON: st_union can return GEOMETRYCOLLECTION when
    # unioning polygons with non-polygon parts (e.g. reference divide edges).
    new_geom <- tryCatch(
      {
        g <- sf::st_make_valid(new_geom)
        g <- tryCatch(sf::st_collection_extract(g, "POLYGON"), error = function(e) g)
        g <- tryCatch(suppressWarnings(sf::st_cast(g, "MULTIPOLYGON")), error = function(e) g)
        if (length(g) > 0L && !sf::st_is_empty(g[[1L]])) g else new_geom
      },
      error = function(e) new_geom)
    if (is.null(new_geom) || length(new_geom) == 0L) next

    # Fill only genuine gaps. The restore unions the current divide with its
    # reconciled members and reference divides, which recovers real losses but
    # also reaches into ground that neighbouring divides already hold, the
    # comment above says it "only ADD area, never drop it", and that is precisely
    # the problem. Measured on two CONUS VPUs it added 263 and 501 km2 of pure
    # overlap while adding 0.0005 km2 of actual coverage, and a later clean pass
    # existed to scrub that off. Subtracting the occupied ground keeps the
    # gap-filling and drops the trespass.
    occ <- tryCatch({
      nb   <- .bx(new_geom[[1L]])
      near <- which(bb[1L, ] <= nb[3L] & bb[3L, ] >= nb[1L] &
                    bb[2L, ] <= nb[4L] & bb[4L, ] >= nb[2L])
      near <- setdiff(near, i)
      cand <- if (length(near))
        near[lengths(sf::st_intersects(d_sfc[near], new_geom)) > 0L] else integer(0)
      if (length(cand))
        sf::st_union(sf::st_sfc(d_geoms[cand], crs = target_crs)) else NULL
    }, error = function(e) NULL)
    if (!is.null(occ)) {
      trimmed <- tryCatch(sf::st_difference(new_geom, occ), error = function(e) NULL)
      if (!is.null(trimmed) && length(trimmed) > 0L &&
          !sf::st_is_empty(trimmed[[1L]])) {
        if (!isTRUE(all.equal(as.numeric(sf::st_area(trimmed)),
                              as.numeric(sf::st_area(new_geom)))))
          n_trimmed <- n_trimmed + 1L
        new_geom <- trimmed
      }
    }

    d_geoms[[i]] <- new_geom[[1L]]
    d_sfc[i]     <- new_geom[[1L]]
    bb[, i]      <- .bx(new_geom[[1L]])   # keep the boxes in step
    restored <- c(restored, did)

    if (isTRUE(verbose)) {
      source_tag <- if (length(ref_used) > 0L) "rec+ref" else "rec"
      restore_log[[length(restore_log) + 1L]] <- data.frame(
        divide_id      = did,
        before_sqkm    = d_area_m2[i] / 1e6,
        expected_sqkm  = max(expected_rec[i], expected_ref[i]) / 1e6,
        expected_rec   = expected_rec[i] / 1e6,
        expected_ref   = expected_ref[i] / 1e6,
        n_rec_members  = length(recs),
        n_ref_comids   = length(ref_used),
        source         = source_tag,
        stringsAsFactors = FALSE)
    }
  }

  if (length(restored) == 0L) {
    message("restore_aggregated_divides: no divides needed restoring")
    return(invisible(character(0)))
  }

  # Rebuild sfc from list so mixed POLYGON/MULTIPOLYGON rows get the correct
  # sfc_GEOMETRY class. In-place `[[i]] <-` leaves the sfc class stale, which
  # makes sf::st_write fail with "type=list; target=double".
  sf::st_geometry(d) <- sf::st_sfc(d_geoms, crs = target_crs)

  # Normalize precision on the restored geometries so that mapshaper's
  # clean step (clean_gpkg_layer / Step 4b) does not drop them due to
  # floating-point topology artifacts introduced by the multi-part union.
  # Snap to 0.01 mm grid (1e5 in projected units), then make_valid to
  # handle any residual ring-order errors from the precision rounding.
  restored_idx <- needs_restore[did_all[needs_restore] %in% restored]
  if (length(restored_idx) > 0L) {
    geom_r <- sf::st_geometry(d)[restored_idx]
    geom_r <- tryCatch(
      sf::st_make_valid(sf::st_set_precision(
        sf::st_sfc(geom_r, crs = target_crs), 1e5)),
      error = function(e) geom_r)
    # Only apply if make_valid didn't drop anything.
    if (length(geom_r) == length(restored_idx) && !any(sf::st_is_empty(geom_r))) {
      geom_list <- sf::st_geometry(d)
      for (k in seq_along(restored_idx)) geom_list[[restored_idx[k]]] <- geom_r[[k]]
      sf::st_geometry(d) <- sf::st_sfc(geom_list, crs = target_crs)
    }
  }

  if ("areasqkm" %in% names(d)) {
    # Recompute in bulk: only the restored rows changed geometry, so leave
    # unchanged rows at their existing value. Vectorized st_area (single sf call)
    # is much faster than per-row vapply at VPU scale. restored_idx is unchanged
    # since the geometry-restore block above; reuse it.
    if (length(restored_idx) > 0L) {
      new_areas_m2 <- suppressWarnings(as.numeric(
        sf::st_area(sf::st_geometry(d)[restored_idx])))
      d$areasqkm[restored_idx] <- new_areas_m2 / 1e6
    }
  }
  d <- .drop_list_cols(d)
  sf::st_write(d, input_gpkg, divide_layer, delete_layer = TRUE, quiet = TRUE)
  message(sprintf(paste0("restore_aggregated_divides: restored %d divide(s) from ",
    "reconciled members; %d trimmed to avoid overlapping a neighbour"),
    length(restored), n_trimmed))

  if (isTRUE(verbose) && length(restore_log) > 0L) {
    rl <- do.call(rbind, restore_log)
    rl$shortfall_pct <- 100 * (1 - rl$before_sqkm / rl$expected_sqkm)
    rl <- rl[order(-rl$shortfall_pct), , drop = FALSE]
    message("-- restore_aggregated_divides detail ----------------------")
    message(sprintf("  by source: %s",
      paste(sprintf("%s=%d", names(table(rl$source)),
        as.integer(table(rl$source))),
      collapse = ", ")))
    message(sprintf("  median shortfall: %.1f%%   worst: %.1f%% (%s)",
      stats::median(rl$shortfall_pct),
      rl$shortfall_pct[1], rl$divide_id[1]))
    message(sprintf("  total area recovered: %.3f sqkm",
      sum(rl$expected_sqkm - rl$before_sqkm)))
    for (r in seq_len(min(5, nrow(rl)))) {
      message(sprintf(
        "  %-12s  %.3f -> %.3f sqkm (%s, %.0f%% shortfall, %d rec + %d ref)",
        rl$divide_id[r], rl$before_sqkm[r], rl$expected_sqkm[r],
        rl$source[r], rl$shortfall_pct[r],
        rl$n_rec_members[r], rl$n_ref_comids[r]))
    }
    message("------------------------------------------------------------")
  }
  invisible(restored)
}

# =====================================================================
# repair_split_divides_at_pois  (was R/repair_split_divides_at_pois.R)
# =====================================================================
#' Repair fp<->cat joins where PRD missed a POI-triggered divide split
#'
#' @description
#' After `resolve_disjoint_flowpath_divides` clears bogus fp<->cat joins, a
#' small class of legitimate flowpaths is left with `has_divide = FALSE` even
#' though their LINE sits entirely inside the downstream neighbour's cat on
#' the same NHD base. Pattern: upstream half (`.1`) of a split NHD, with a
#' POI at its outlet: PRD should have carved a new divide boundary at that
#' POI but the full divide ended up with the downstream half.
#'
#' This pass re-runs the raster-based splitter (`split_divides_d8`)
#' for just those cases, using the POI location (= fp_up outlet = fp_down
#' inlet) as the split point. Upstream piece becomes fp_up's new cat;
#' downstream piece replaces fp_down's existing cat.
#'
#' @param input_gpkg Path to the pre-ngen GeoPackage (layers: flowpaths,
#'   divides, pois; optional: network).
#' @param flowpath_layer Name of the flowpaths layer in `input_gpkg`. Default
#'   `"flowpaths"`.
#' @param divides_layer Name of the divides layer in `input_gpkg`. Default
#'   `"divides"`.
#' @param pois_layer Name of the POIs layer in `input_gpkg`. Default `"pois"`.
#' @param network_layer Name of the network layer in `input_gpkg`. Default
#'   `"network"`.
#' @param fdr_path,fac_path Paths to FDR / FAC rasters (same files PRD used
#'   earlier in the pipeline).
#' @param fallback_fdr_path,fallback_fac_path Optional paths to a fallback
#'   FDR/FAC pair.  When a split attempt returns NULL with the primary
#'   rasters (split point lands on NoData), the splitter is retried with the
#'   fallback.  `NULL` disables the retry (default).  Ignored when equal to
#'   the primary paths (path equality check, same as `add_flowlines`).
#' @param coverage_threshold Minimum fraction of fp_up's line that must sit
#'   inside fp_down's cat for the pair to qualify. Default 0.9.
#' @param num_workers Parallel workers for the splitter. Default 4.
#'
#' @return Invisibly, a data.frame of repaired `{fp_up, fp_down, new_cat,
#'   area_km2}` rows.
#' @keywords internal
repair_split_divides_at_pois <- function(input_gpkg,
                                         flowpath_layer    = "flowpaths",
                                         divides_layer     = "divides",
                                         pois_layer        = "pois",
                                         network_layer     = "network",
                                         fdr_path,
                                         fac_path,
                                         fallback_fdr_path = NULL,
                                         fallback_fac_path = NULL,
                                         coverage_threshold = 0.9,
                                         num_workers = 4L) {
  stopifnot(.raster_accessible(fdr_path), .raster_accessible(fac_path))

  if (!hfutils::layer_exists(input_gpkg, pois_layer)) {
    message("repair_split_divides_at_pois: no '", pois_layer, "' layer -> skipped (0 pois)")
    return(invisible(NULL))
  }
  fp   <- sf::st_read(input_gpkg, flowpath_layer, quiet = TRUE)
  div  <- sf::st_read(input_gpkg, divides_layer,  quiet = TRUE)
  pois <- sf::st_read(input_gpkg, pois_layer,     quiet = TRUE)
  net  <- if (hfutils::layer_exists(input_gpkg, network_layer))
    sf::st_read(input_gpkg, network_layer, quiet = TRUE) else NULL

  fp$flowpath_id   <- as.character(fp$flowpath_id)
  fp$flowpath_toid <- as.character(fp$flowpath_toid)
  fp$divide_id     <- as.character(fp$divide_id)
  div$divide_id    <- as.character(div$divide_id)

  target_crs <- sf::st_crs(div)
  if (sf::st_crs(fp)   != target_crs) fp   <- sf::st_transform(fp,   target_crs)
  if (sf::st_crs(pois) != target_crs) pois <- sf::st_transform(pois, target_crs)

  # Candidate fp_ups: no cat, has POI at outlet
  need_cat <- is.na(fp$divide_id) | !isTRUE(fp$has_divide) |
    !(fp$divide_id %in% div$divide_id)
  if ("has_divide" %in% names(fp))
    need_cat <- is.na(fp$divide_id) | !fp$has_divide
  has_poi <- fp$flowpath_id %in% as.character(pois$flowpath_id)
  cand <- which(need_cat & has_poi)
  if (length(cand) == 0L) {
    message("repair_split_divides_at_pois: no candidate fps")
    return(invisible(data.frame()))
  }

  # Downstream fp lookup: pre-ngen: toid is another fp_id. Post-ngen: toid
  # is a nexus_id which maps to another fp via the nexus layer.
  fp_ix <- stats::setNames(seq_len(nrow(fp)), fp$flowpath_id)
  nex <- if (hfutils::layer_exists(input_gpkg, "nexus"))
    sf::st_drop_geometry(sf::st_read(input_gpkg, "nexus", quiet = TRUE)) else NULL
  toid_to_fp_ix <- if (!is.null(nex) && "nexus_id" %in% names(nex)) {
    # nex$nexus_toid is the downstream flowpath_id
    nex_lu <- stats::setNames(as.character(nex$nexus_toid), as.character(nex$nexus_id))
    function(to) {
      ix <- fp_ix[to]
      if (!is.na(ix)) return(unname(ix))
      fp2 <- nex_lu[to]
      if (is.na(fp2)) return(NA_integer_)
      unname(fp_ix[fp2])
    }
  } else {
    function(to) unname(fp_ix[to])
  }

  # Build per-candidate parent/split-point records
  parent_list <- list()
  split_list  <- list()
  report      <- list()
  for (i in cand) {
    fp_up   <- fp$flowpath_id[i]
    to      <- fp$flowpath_toid[i]
    dj      <- toid_to_fp_ix(to)
    if (is.na(dj)) next                           # terminal
    fp_down <- fp$flowpath_id[dj]
    cat_dn  <- fp$divide_id[dj]
    if (is.na(cat_dn)) next                       # downstream also no-cat
    cat_row <- match(cat_dn, div$divide_id)
    if (is.na(cat_row)) next

    up_line <- sf::st_geometry(fp)[[i]]
    cat_g   <- sf::st_geometry(div)[[cat_row]]

    # Containment: fp_up line must be (almost) inside cat_down
    total_m  <- as.numeric(sf::st_length(sf::st_sfc(up_line, crs = target_crs)))
    if (total_m <= 0) next
    inter_g  <- sf::st_intersection(up_line, cat_g)
    inside_m <- as.numeric(sf::st_length(sf::st_sfc(inter_g, crs = target_crs)))
    if (inside_m / total_m < coverage_threshold) next

    # Split point = POI location. In nominal outlet-semantics this is the same
    # as fp_up's outlet, but real snapped POIs can sit along the segment by
    # enough that a perpendicular fallback through the endpoint cuts the wrong
    # cross-section of the parent divide.
    up_c <- sf::st_coordinates(sf::st_sfc(up_line, crs = target_crs))
    jp   <- up_c[nrow(up_c), 1:2]
    poi_rows <- which(as.character(pois$flowpath_id) == fp_up)
    split_pt <- if (length(poi_rows) > 0L) {
      poi_geoms <- sf::st_geometry(pois[poi_rows, ])
      outlet_pt <- sf::st_sfc(sf::st_point(jp), crs = target_crs)
      if (sf::st_crs(poi_geoms) != target_crs) {
        poi_geoms <- sf::st_transform(poi_geoms, target_crs)
      }
      nearest <- which.min(as.numeric(sf::st_distance(poi_geoms, outlet_pt)))
      poi_geoms[nearest]
    } else {
      sf::st_sfc(sf::st_point(jp), crs = target_crs)
    }
    if (sf::st_crs(split_pt) != target_crs) {
      split_pt <- sf::st_transform(split_pt, target_crs)
    }
    down_c <- sf::st_coordinates(sf::st_sfc(sf::st_geometry(fp)[[dj]], crs = target_crs))
    down_start <- down_c[1, 1:2]

    # Parent divide entry: splitter wants divide_id + flowpath_id + geom
    parent_list[[length(parent_list) + 1L]] <- sf::st_sf(
      divide_id   = cat_dn,
      flowpath_id = cat_dn,
      geometry    = sf::st_sfc(cat_g, crs = target_crs))

    split_list[[length(split_list) + 1L]] <- sf::st_sf(
      divide_id     = cat_dn,
      flowpath_id   = fp_up,     # upstream piece gets fp_up's id
      downstream_id = fp_down,   # downstream piece stays fp_down
      reach_line_id = fp_up,
      sib_line_id   = fp_down,
      downstream_x  = as.numeric(down_start[1]),
      downstream_y  = as.numeric(down_start[2]),
      geometry      = split_pt)

    report[[length(report) + 1L]] <- data.frame(
      fp_up = fp_up, fp_down = fp_down, parent_cat = cat_dn,
      coverage_in_parent = round(inside_m / total_m, 4),
      stringsAsFactors = FALSE)
  }

  if (length(parent_list) == 0L) {
    message("repair_split_divides_at_pois: no valid candidates")
    return(invisible(data.frame()))
  }

  parent_divides <- do.call(rbind, parent_list)
  split_points   <- do.call(rbind, split_list)
  # Dedup parent divides (same cat possibly shared by >1 fp_up, unusual)
  parent_divides <- parent_divides[!duplicated(parent_divides$divide_id), ]
  reach_lines <- fp[, intersect(names(fp), c("flowpath_id")), drop = FALSE]
  reach_lines$line_id <- as.character(reach_lines$flowpath_id)
  reach_lines <- reach_lines[, c("line_id", attr(reach_lines, "sf_column")), drop = FALSE]

  message(sprintf(
    "repair_split_divides_at_pois: attempting %d split(s) across %d parent cat(s)",
    nrow(split_points), nrow(parent_divides)))

  # Parallel setup (raster serialization to each worker) dominates cost at
  # small N. For <=16 parents, run sequentially: worker setup overhead exceeds
  # the sequential runtime.
  use_parallel <- num_workers > 1L && nrow(parent_divides) > 16L
  eff_workers  <- if (use_parallel)
    min(as.integer(num_workers), nrow(parent_divides)) else 1L

  .run_splitter <- function(fdr, fac) {
    old_engine <- getOption("hydrofabric.partition_engine", NULL)
    old_cov <- getOption("hydrofabric.split_min_reach_coverage", NULL)
    on.exit({
      if (is.null(old_engine)) options(hydrofabric.partition_engine = NULL)
      else options(hydrofabric.partition_engine = old_engine)
      if (is.null(old_cov)) options(hydrofabric.split_min_reach_coverage = NULL)
      else options(hydrofabric.split_min_reach_coverage = old_cov)
    }, add = TRUE)
    options(
      hydrofabric.partition_engine = "legacy",
      hydrofabric.split_min_reach_coverage = coverage_threshold
    )
    tryCatch(
      split_divides_d8(
        divides              = parent_divides,
        split_points         = split_points,
        fdr_path             = fdr,
        fac_path             = fac,
        reach_lines          = reach_lines,
        snap_distance_m      = 500,
        simplify_tolerance_m = 25,
        min_piece_area_m2    = 1,
        min_area_m2          = 10000,
        fallback_geometric   = TRUE,
        parallel             = use_parallel,
        num_workers          = eff_workers,
        use_legacy_reconcile = FALSE),
      error = function(e) {
        message("repair_split_divides_at_pois: splitter error: ",
          conditionMessage(e))
        NULL
      })
  }

  split_raw <- .run_splitter(fdr_path, fac_path)

  # If the primary rasters yielded no result (split point on NoData domain),
  # retry with the GEDTM30 fallback when available and distinct from primary.
  if ((is.null(split_raw) || nrow(split_raw) == 0L) &&
    !is.null(fallback_fdr_path) && nzchar(fallback_fdr_path) &&
    !identical(normalizePath(fdr_path,          mustWork = FALSE),
      normalizePath(fallback_fdr_path,  mustWork = FALSE))) {
    message("repair_split_divides_at_pois: retrying with GEDTM30 fallback rasters")
    split_raw <- .run_splitter(fallback_fdr_path, fallback_fac_path)
  }

  if (is.null(split_raw) || nrow(split_raw) == 0L) {
    return(invisible(data.frame()))
  }
  split_raw <- sf::st_as_sf(split_raw)

  # The splitter labels each output piece with the flowpath_id from the
  # split_points table: upstream = fp_up, downstream = fp_down. Use that to
  # decide which pieces become new fp_up cats vs replace the fp_down cat.
  split_raw$flowpath_id <- as.character(split_raw$flowpath_id)
  target_up_ids <- unique(split_points$flowpath_id)

  # Ensure output CRS matches div CRS
  if (sf::st_crs(split_raw) != target_crs)
    split_raw <- sf::st_transform(split_raw, target_crs)

  # Partition: upstream pieces, downstream pieces
  up_pieces <- split_raw[split_raw$flowpath_id %in% target_up_ids, , drop = FALSE]
  dn_pieces <- split_raw[!(split_raw$flowpath_id %in% target_up_ids), , drop = FALSE]

  # Build new rows for div layer: one per fp_up (= a new cat-fp_up) and
  # replace the geom of each fp_down's cat with the downstream piece.
  repaired <- data.frame(
    fp_up      = character(0),
    fp_down    = character(0),
    new_cat    = character(0),
    area_km2   = numeric(0),
    stringsAsFactors = FALSE)

  # fp_up -> fp_down mapping from split_points
  up_to_down <- stats::setNames(as.character(split_points$downstream_id),
    as.character(split_points$flowpath_id))
  up_to_parent <- stats::setNames(as.character(split_points$divide_id),
    as.character(split_points$flowpath_id))

  # Build upstream-cat rows (new divide entries). Validate that the resulting
  # upstream piece actually contains fp_up's line (coverage >= threshold);
  # reject tiny slivers that don't represent a real catchment.
  new_div_rows <- list()
  fp_up_geoms <- stats::setNames(
    lapply(match(as.character(split_points$flowpath_id), fp$flowpath_id),
      function(i) sf::st_geometry(fp)[[i]]),
    as.character(split_points$flowpath_id))
  rejected <- character(0)

  # Id convention: pre-ngen uses plain integers, ngen uses "fp-N"/"cat-N"
  # prefixes. Derive the cat_id from fp_up_id so new divides fit the schema.
  fp_to_cat_id <- function(fp_id) {
    if (grepl("^fp-", fp_id)) sub("^fp-", "cat-", fp_id) else fp_id
  }

  for (k in seq_len(nrow(up_pieces))) {
    fp_up_id <- up_pieces$flowpath_id[k]
    g_sfc    <- sf::st_geometry(up_pieces)[k]  # already in target_crs
    # Ensure valid polygon: raster-to-vector splits sometimes produce rings
    # that self-intersect or have near-duplicate vertices along the cut line.
    g_sfc <- tryCatch(sf::st_make_valid(g_sfc), error = function(e) g_sfc)
    area_m2  <- as.numeric(sf::st_area(g_sfc))
    if (is.na(area_m2) || area_m2 <= 0) next

    # Coverage check: fp_up's line inside the upstream piece
    up_line  <- fp_up_geoms[[fp_up_id]]
    line_sfc <- sf::st_sfc(up_line, crs = target_crs)
    total_m  <- as.numeric(sf::st_length(line_sfc))
    inside_m <- tryCatch(
      {
        inter <- sf::st_intersection(line_sfc, g_sfc)
        if (length(inter) == 0L) 0
        else sum(as.numeric(sf::st_length(inter)), na.rm = TRUE)
      },
      error = function(e) 0)
    new_cov <- if (total_m > 0) inside_m / total_m else 0
    if (is.na(new_cov) || new_cov < coverage_threshold) {
      rejected <- c(rejected, fp_up_id)
      next
    }
    g <- g_sfc[[1L]]

    new_cat_id <- fp_to_cat_id(fp_up_id)
    # Structural collision: fp_up lost cat-N to its downstream neighbor (which
    # became the parent divide we're carving). fp_to_cat_id(fp_up) == cat_dn is
    # the expected case here, the new upstream piece can't reuse the parent id
    # since the parent still exists as fp_down's divide until the split completes.
    # Also guard against same-call id reuse if multiple fp_ups carve the same parent.
    existing_ids <- c(as.character(div$divide_id),
      vapply(new_div_rows, function(x)
        as.character(sf::st_drop_geometry(x)$divide_id),
      character(1L)))
    if (new_cat_id %in% existing_ids) {
      # When the collision is exactly the parent being carved, use "-up" suffix
      # to make the structural reason explicit; otherwise increment numerically.
      candidate <- if (new_cat_id == cat_dn) {
        paste0(new_cat_id, "-up")
      } else {
        suffix <- 1L
        repeat {
          s <- sprintf("%s.%d", new_cat_id, suffix)
          if (!s %in% existing_ids) break
          suffix <- suffix + 1L
        }
        s
      }
      if (candidate %in% existing_ids) {
        suffix <- 1L
        repeat {
          s <- sprintf("%s.%d", candidate, suffix)
          if (!s %in% existing_ids) break
          suffix <- suffix + 1L
        }
        candidate <- s
      }
      message(sprintf(
        "repair_split_divides_at_pois: id collision %s -> %s for fp_up=%s",
        new_cat_id, candidate, fp_up_id))
      new_cat_id <- candidate
    }
    # Stamp the owner. A cat- divide with no flowpath_id is classified as
    # landscape by .hf_network_divides regardless of its prefix, so an unstamped
    # divide silently leaves every area, coverage and overlap figure, its land
    # off the books with no geometry lost and no check failing. That has cost
    # real area before (0.2249 km2 on VPU 03N, 0.4362 km2 on VPU 12).
    row <- data.frame(
      divide_id    = new_cat_id,
      areasqkm     = area_m2 / 1e6,
      has_flowline = TRUE,
      ds_id        = NA_character_,
      type         = "network",
      flowpath_id  = fp_up_id,
      flowpath_toid = NA_character_,
      stringsAsFactors = FALSE)
    new_div_rows[[length(new_div_rows) + 1L]] <-
      sf::st_sf(row, geometry = g_sfc)
    repaired <- rbind(repaired, data.frame(
      fp_up      = fp_up_id,
      fp_down    = up_to_down[fp_up_id],
      new_cat    = new_cat_id,
      area_km2   = area_m2 / 1e6,
      stringsAsFactors = FALSE))
  }

  if (length(rejected) > 0L) {
    message(sprintf(
      "repair_split_divides_at_pois: %d split(s) rejected (upstream piece failed coverage): %s",
      length(rejected), paste(rejected, collapse = ",")))
  }

  # Validate the downstream piece ALSO contains fp_down's line with
  # adequate coverage. If not, the split got the geometry wrong (e.g. raster
  # snap placed the POI in a pixel whose upstream contributing area eats
  # into fp_down's reach): revert the repair for that pair.
  fp_down_geoms <- stats::setNames(
    lapply(match(as.character(split_points$downstream_id), fp$flowpath_id),
      function(i) sf::st_geometry(fp)[[i]]),
    as.character(split_points$downstream_id))

  accepted_parents <- up_to_parent[repaired$fp_up]
  revert_fp_ups <- character(0)
  for (k in seq_len(nrow(dn_pieces))) {
    fp_down_id <- dn_pieces$flowpath_id[k]
    # parent divide uses "cat-" prefix in ngen schema; fp_down_id may use "fp-"
    parent_cat <- fp_to_cat_id(fp_down_id)
    if (!(parent_cat %in% accepted_parents)) next

    # fp_down coverage in the new downstream piece
    down_line <- fp_down_geoms[[fp_down_id]]
    if (is.null(down_line)) next
    dn_g_sfc  <- sf::st_geometry(dn_pieces)[k]
    line_sfc  <- sf::st_sfc(down_line, crs = target_crs)
    total_m   <- as.numeric(sf::st_length(line_sfc))
    inside_m  <- tryCatch(
      {
        inter <- sf::st_intersection(line_sfc, dn_g_sfc)
        if (length(inter) == 0L) 0
        else sum(as.numeric(sf::st_length(inter)), na.rm = TRUE)
      },
      error = function(e) 0)
    new_cov <- if (total_m > 0) inside_m / total_m else 0

    if (is.na(new_cov) || new_cov < coverage_threshold) {
      # Identify the fp_up(s) whose parent this is, mark for revert
      revert_fp_ups <- c(revert_fp_ups,
        repaired$fp_up[repaired$fp_down == fp_down_id])
      next
    }

    # Enforce strict non-overlap: downstream piece = original_cat - up_piece.
    # Post-cleanup repair bypasses mapshaper's corruption but we still want
    # a clean carve: use parent polygon minus the union of upstream pieces.
    row <- match(parent_cat, div$divide_id)
    if (is.na(row)) next
    up_for_parent <- repaired[repaired$fp_down == fp_down_id, , drop = FALSE]
    if (nrow(up_for_parent) == 0L) {
      dn_clean <- dn_g_sfc
    } else {
      up_pieces_for_parent <- up_pieces[
        up_pieces$flowpath_id %in% up_for_parent$fp_up, , drop = FALSE]
      up_union <- sf::st_union(sf::st_geometry(up_pieces_for_parent))
      parent_g <- sf::st_geometry(div)[row]
      dn_clean <- tryCatch(
        sf::st_make_valid(sf::st_difference(parent_g, up_union)),
        error = function(e) sf::st_make_valid(dn_g_sfc))
      if (length(dn_clean) == 0L ||
        sf::st_is_empty(dn_clean[[1L]]) ||
        as.numeric(sf::st_area(dn_clean)) < 1e3) {
        dn_clean <- sf::st_make_valid(dn_g_sfc)
      }
    }
    sf::st_geometry(div)[[row]] <- dn_clean[[1L]]
    if ("areasqkm" %in% names(div)) {
      div$areasqkm[row] <- as.numeric(sf::st_area(dn_clean)) / 1e6
    }
  }

  # Revert failed repairs: drop the new upstream cats, un-set fp_up assignment
  if (length(revert_fp_ups) > 0L) {
    message(sprintf(
      "repair_split_divides_at_pois: %d pair(s) reverted (downstream piece failed coverage): %s",
      length(revert_fp_ups), paste(revert_fp_ups, collapse = ",")))
    # Translate the reverted FLOWPATH ids into the DIVIDE ids they produced
    # before `repaired` is filtered, that mapping is the only link between the
    # two id namespaces.
    #
    # Testing a new row's `divide_id` (cat-N) against the reverted flowpath ids
    # (fp-N) can never match, so a reverted split's divide is appended anyway: an
    # unowned duplicate nested inside an owned divide, which
    # divides_non_overlapping cannot see because it skips divides with no
    # flowpath_id.
    revert_cats   <- as.character(repaired$new_cat[repaired$fp_up %in% revert_fp_ups])
    keep_repaired <- !(repaired$fp_up %in% revert_fp_ups)
    repaired <- repaired[keep_repaired, , drop = FALSE]
    new_div_rows <- new_div_rows[vapply(new_div_rows, function(x)
      !(as.character(sf::st_drop_geometry(x)$divide_id) %in% revert_cats),
    logical(1L))]
  }

  # Append new upstream cat rows: build a combined sf matching div's schema.
  if (length(new_div_rows) > 0L) {
    div_cols    <- names(sf::st_drop_geometry(div))
    div_geom_nm <- attr(div, "sf_column")

    # Build a conforming attribute data.frame + geometry list
    attr_rows <- lapply(new_div_rows, function(x) {
      xd <- sf::st_drop_geometry(x)
      miss <- setdiff(div_cols, names(xd))
      for (nm in miss) xd[[nm]] <- NA
      xd[, div_cols, drop = FALSE]
    })
    new_attrs <- do.call(rbind, attr_rows)
    new_geoms <- sf::st_sfc(
      lapply(new_div_rows, function(x) sf::st_geometry(x)[[1L]]),
      crs = target_crs)

    new_sf <- sf::st_sf(new_attrs, geometry = new_geoms)
    names(new_sf)[names(new_sf) == "geometry"] <- div_geom_nm
    sf::st_geometry(new_sf) <- div_geom_nm

    # Stamp the owner. Conforming to div's schema above fills every absent column
    # with NA, and flowpath_id is absent from the splitter's output, so a new
    # divide arrived unowned while the flowpath side was set below. An unowned
    # cat- divide is excluded from every area and coverage figure by
    # .hf_network_divides, so its land silently left the books: 0.2249 km2 on
    # VPU 03N and 0.4362 km2 on VPU 12, the only two VPUs that run this repair.
    # The divide is created for a specific upstream flowpath, so record it.
    if ("flowpath_id" %in% names(new_sf) && nrow(repaired) > 0L) {
      own <- stats::setNames(as.character(repaired$fp_up),
        as.character(repaired$new_cat))
      hit <- as.character(new_sf$divide_id) %in% names(own)
      new_sf$flowpath_id[hit] <-
        unname(own[as.character(new_sf$divide_id)[hit]])
    }

    div <- rbind(div, new_sf)
  }

  # Update fp rows for the newly-catted upstream fps. new_cat was stored in
  # `repaired$new_cat` using fp_to_cat_id mapping so it matches div schema.
  if (nrow(repaired) > 0L) {
    ix <- match(repaired$fp_up, fp$flowpath_id)
    fp$divide_id[ix] <- repaired$new_cat
    if ("has_divide" %in% names(fp)) fp$has_divide[ix] <- TRUE
  }

  # Strip list-cols and write back
  for (obj_name in c("fp", "div")) {
    x <- get(obj_name)
    x <- .drop_list_cols(x)
    assign(obj_name, x)
  }

  sf::st_write(fp,  input_gpkg, flowpath_layer, delete_layer = TRUE, quiet = TRUE)
  sf::st_write(div, input_gpkg, divides_layer,  delete_layer = TRUE, quiet = TRUE)

  message(sprintf(
    "repair_split_divides_at_pois: repaired %d fp<->cat pair(s)",
    nrow(repaired)))

  invisible(repaired)
}

# =====================================================================
# coalesce_fragmented_flowpaths  (was R/coalesce_fragmented_flowpaths.R)
# =====================================================================
#' Coalesce fragmented flowpaths with their downstream same-mainstem neighbor
#'
#' @description
#' After aggregation + PRD + restore, some flowpaths end up with severely
#' fragmented divides (sliver-dominated MULTIPOLYGONs with tiny total area)
#' while their line geometry is valid. The missing territory was absorbed by
#' the same-mainstem downstream neighbor, so the line overshoots its divide.
#' This pass merges such flowpaths into their immediate downstream neighbor
#' when the merge is safe.
#'
#' Merge rule: fp_up is absorbed into fp_down when ALL hold:
#'   1. fp_up.flowpath_toid == fp_down.flowpath_id (directly adjacent)
#'   2. fp_up.mainstem == fp_down.mainstem (same stream)
#'   3. No POI attached to fp_up (preserves POI=outlet semantics)
#'   4. fp_up's current line coverage by its own divide is below
#'      `coverage_threshold`
#'
#' Chains cascade: if fp_up -> fp_mid -> fp_down and both fp_up and fp_mid are
#' absorbable, all three merge into fp_down.
#'
#' @param input_gpkg Path to the intermediate GeoPackage.
#' @param flowpath_layer,divides_layer,flowlines_layer,pois_layer,network_layer
#'   Layer names. Defaults match the pre-ngen pipeline.
#' @param coverage_threshold fp_up coverage below this triggers merge candidacy
#'   (default 0.9).
#'
#' @return Invisibly, a data.frame of `{up, down}` merge pairs.
#' @keywords internal
coalesce_fragmented_flowpaths <- function(input_gpkg,
                                          flowpath_layer   = "flowpaths",
                                          divides_layer    = "divides",
                                          flowlines_layer  = "flowlines",
                                          pois_layer       = "pois",
                                          network_layer    = "network",
                                          coverage_threshold = 0.9) {
  fp  <- sf::st_read(input_gpkg, flowpath_layer,  quiet = TRUE)
  div <- sf::st_read(input_gpkg, divides_layer,   quiet = TRUE)
  fl  <- sf::st_read(input_gpkg, flowlines_layer, quiet = TRUE)
  pois <- if (hfutils::layer_exists(input_gpkg, pois_layer))
    sf::st_read(input_gpkg, pois_layer, quiet = TRUE) else NULL
  net  <- if (hfutils::layer_exists(input_gpkg, network_layer))
    sf::st_read(input_gpkg, network_layer, quiet = TRUE) else NULL

  target_crs <- sf::st_crs(fp)
  if (sf::st_crs(div) != target_crs) div <- sf::st_transform(div, target_crs)
  if (sf::st_crs(fl)  != target_crs) fl  <- sf::st_transform(fl,  target_crs)

  fp$flowpath_id    <- as.character(fp$flowpath_id)
  fp$flowpath_toid  <- as.character(fp$flowpath_toid)
  fp$divide_id      <- as.character(fp$divide_id)
  div$divide_id     <- as.character(div$divide_id)
  fl$flowpath_id    <- as.character(fl$flowpath_id)
  if ("flowpath_toid" %in% names(fl))
    fl$flowpath_toid <- as.character(fl$flowpath_toid)

  # -- Coverage per fp (vectorized GEOS on equal-area) -----------------------
  fp5070  <- sf::st_transform(fp,  5070)
  div5070 <- sf::st_transform(div, 5070)
  d_ix    <- match(fp$divide_id, div$divide_id)
  covers  <- rep(NA_real_, nrow(fp))
  valid   <- !is.na(d_ix)
  if (any(valid)) {
    lines_g <- geos::as_geos_geometry(sf::st_geometry(fp5070)[valid])
    divs_g  <- geos::as_geos_geometry(sf::st_geometry(div5070)[d_ix[valid]])
    inter_g <- tryCatch(geos::geos_intersection(lines_g, divs_g),
      error = function(e) NULL)
    if (!is.null(inter_g)) {
      lt <- suppressWarnings(as.numeric(geos::geos_length(lines_g)))
      li <- suppressWarnings(as.numeric(geos::geos_length(inter_g)))
      covers[valid] <- ifelse(lt > 0, li / lt, NA_real_)
    }
  }

  has_poi <- stats::setNames(logical(nrow(fp)), fp$flowpath_id)
  if (!is.null(pois) && "flowpath_id" %in% names(pois)) {
    fp_ids_with_poi <- unique(as.character(pois$flowpath_id[!is.na(pois$flowpath_id)]))
    has_poi[intersect(fp_ids_with_poi, names(has_poi))] <- TRUE
  }

  is_absorbable <- !is.na(covers) & covers < coverage_threshold & !has_poi

  # -- Resolve leader per fp by walking downstream through absorbable chain
  fp_ix   <- stats::setNames(seq_len(nrow(fp)), fp$flowpath_id)
  leader  <- seq_len(nrow(fp))
  for (k in seq_len(nrow(fp))) {
    if (!is_absorbable[k]) next
    cur <- k
    seen <- integer(0)
    repeat {
      to_id <- fp$flowpath_toid[cur]
      if (is.na(to_id) || !(to_id %in% names(fp_ix))) {
        leader[k] <- k
        break
      }
      j <- fp_ix[[to_id]]
      if (!identical(as.character(fp$mainstem[j]),
        as.character(fp$mainstem[k]))) {
        leader[k] <- k
        break
      }
      if (j %in% seen) {
        leader[k] <- k
        break
      }
      seen <- c(seen, j)
      if (!is_absorbable[j]) {
        leader[k] <- j
        break
      }
      cur <- j
    }
  }

  absorbed_mask    <- leader != seq_len(nrow(fp))
  absorbed_fp_ids  <- fp$flowpath_id[absorbed_mask]
  absorbed_div_ids <- fp$divide_id[absorbed_mask]

  merges_df <- data.frame(
    up   = absorbed_fp_ids,
    down = fp$flowpath_id[leader[absorbed_mask]],
    stringsAsFactors = FALSE)

  if (nrow(merges_df) == 0L) {
    message("coalesce_fragmented_flowpaths: no merges needed")
    return(invisible(merges_df))
  }

  # -- Group fps by leader, union geoms + attrs per group --------------------
  groups       <- split(seq_len(nrow(fp)), leader)
  fp_geoms_new  <- sf::st_geometry(fp)
  div_geoms_new <- sf::st_geometry(div)
  folded_div_ids <- character(0)   # divides whose geometry reached a leader

  for (key in names(groups)) {
    grp <- groups[[key]]
    if (length(grp) < 2L) next
    lead_i <- as.integer(key)
    if (is.na(lead_i) || lead_i < 1L || lead_i > nrow(fp)) next
    others <- setdiff(grp, lead_i)

    # Union lines + line_merge (fp_up's downstream endpoint matches fp_down's
    # upstream endpoint by construction, so this yields a single LINESTRING).
    all_line <- sf::st_union(fp_geoms_new[grp])
    merged_g <- tryCatch(
      geos::geos_line_merge(geos::as_geos_geometry(all_line)),
      error = function(e) NULL)
    if (!is.null(merged_g)) {
      merged_sfc <- sf::st_as_sfc(merged_g, crs = target_crs)
      fp_geoms_new[[lead_i]] <- merged_sfc[[1]]
    }

    # Union divides into the leader's row
    d_grp <- match(fp$divide_id[grp], div$divide_id)
    d_grp <- d_grp[!is.na(d_grp)]
    # Record which divides actually had their geometry folded into the leader.
    # The removal below must drop ONLY those: this union is conditional (the
    # leader may have no resolvable divide, or fewer than two may resolve) while
    # the removal was not, so any group that skipped the union still lost its
    # absorbed divides, and their land with them. That cost VPU 10U 60.97 km2.
    if (length(d_grp) >= 2L) {
      lead_d <- match(fp$divide_id[lead_i], div$divide_id)
      if (!is.na(lead_d)) {
        all_div <- sf::st_union(div_geoms_new[d_grp])
        div_geoms_new[[lead_d]] <- all_div[[1]]
        folded_div_ids <- c(folded_div_ids,
          setdiff(div$divide_id[d_grp], div$divide_id[lead_d]))
      }
    }

    # Attribute merges on the leader row
    fp$member_comid[lead_i] <- paste(
      fp$member_comid[c(lead_i, others)], collapse = ",")
    if ("areasqkm" %in% names(fp))
      fp$areasqkm[lead_i] <- sum(fp$areasqkm[grp], na.rm = TRUE)

    # Retarget: fps whose toid is an absorbed id now points to leader
    absorbed_here <- fp$flowpath_id[others]
    fp$flowpath_toid[fp$flowpath_toid %in% absorbed_here] <- fp$flowpath_id[lead_i]
    fl$flowpath_id[fl$flowpath_id %in% absorbed_here] <- fp$flowpath_id[lead_i]
    if ("flowpath_toid" %in% names(fl))
      fl$flowpath_toid[fl$flowpath_toid %in% absorbed_here] <- fp$flowpath_id[lead_i]
    if (!is.null(pois) && "flowpath_id" %in% names(pois))
      pois$flowpath_id[pois$flowpath_id %in% absorbed_here] <- fp$flowpath_id[lead_i]
    if (!is.null(net)) {
      if ("flowpath_id"   %in% names(net))
        net$flowpath_id[net$flowpath_id %in% absorbed_here] <- fp$flowpath_id[lead_i]
      if ("flowpath_toid" %in% names(net))
        net$flowpath_toid[net$flowpath_toid %in% absorbed_here] <- fp$flowpath_id[lead_i]
      if ("divide_id" %in% names(net)) {
        absorbed_d_here <- fp$divide_id[others]
        net$divide_id[net$divide_id %in% absorbed_d_here] <- fp$divide_id[lead_i]
      }
    }
  }

  # -- Rebuild fp / div with new geoms, drop absorbed rows -------------------
  fp_out <- sf::st_sf(sf::st_drop_geometry(fp),
    geometry = sf::st_sfc(fp_geoms_new, crs = target_crs))
  fp_out <- fp_out[!(fp_out$flowpath_id %in% absorbed_fp_ids), ]
  if ("lengthkm" %in% names(fp_out)) {
    fp_out$lengthkm <- as.numeric(sf::st_length(
      sf::st_transform(fp_out, 5070))) / 1000
  }

  div_out <- sf::st_sf(sf::st_drop_geometry(div),
    geometry = sf::st_sfc(div_geoms_new, crs = target_crs))
  # Drop only the divides whose geometry actually reached a leader. The rest
  # keep their land: their flowpath is gone, so they are unowned, and unowned
  # divides get merged into a neighbour rather than deleted.
  div_out <- div_out[!(div_out$divide_id %in% folded_div_ids), ]
  stranded <- setdiff(absorbed_div_ids, folded_div_ids)
  stranded <- stranded[!is.na(stranded) & nzchar(stranded)]
  if (length(stranded)) {
    smask <- div_out$divide_id %in% stranded
    if (any(smask)) {
      # Pass the mask, do NOT stamp flowpath_id: this layer has no such column
      # at this point in the pipeline, and creating one silently changes the
      # schema. An all-NULL flowpath_id also makes .hf_network_divides treat
      # every divide as unowned, which removes the whole layer from the area
      # and coverage accounting.
      ab <- .absorb_orphan_divides(div_out, smask)
      div_out <- ab$div
      message(sprintf(paste0("coalesce_fragmented_flowpaths: %d divide(s) were ",
        "not folded into a leader; absorbed %d into neighbours, kept %d"),
        sum(smask), ab$n_absorbed, ab$n_kept))
    }
  }
  if ("areasqkm" %in% names(div_out)) {
    div_out$areasqkm <- as.numeric(sf::st_area(
      sf::st_transform(div_out, 5070))) / 1e6
  }

  # Strip list-cols (gpkg can't serialize them)
  for (x_name in c("fp_out", "div_out", "fl", "pois", "net")) {
    x <- get(x_name)
    if (is.null(x)) next
    x <- .drop_list_cols(x)
    assign(x_name, x)
  }

  sf::st_write(fp_out,  input_gpkg, flowpath_layer,
    delete_layer = TRUE, quiet = TRUE)
  sf::st_write(div_out, input_gpkg, divides_layer,
    delete_layer = TRUE, quiet = TRUE)
  sf::st_write(fl,      input_gpkg, flowlines_layer,
    delete_layer = TRUE, quiet = TRUE)
  if (!is.null(pois))
    sf::st_write(pois, input_gpkg, pois_layer,
      delete_layer = TRUE, quiet = TRUE)
  if (!is.null(net))
    sf::st_write(net, input_gpkg, network_layer,
      delete_layer = TRUE, quiet = TRUE)

  message(sprintf(
    "coalesce_fragmented_flowpaths: merged %d flowpath(s) into %d leader(s)",
    nrow(merges_df), length(unique(merges_df$down))))

  invisible(merges_df)
}

# =====================================================================
# dissolve_multipart_divides  (was R/dissolve_multipart_divides.R)
# =====================================================================
#' Dissolve multi-part divide polygons to their largest component
#'
#' @description
#' Lightweight replacement for `hfutils::clean_geometry`'s intra-divide
#' fragmentation handling. For each divide whose geometry is a MULTIPOLYGON
#' with disjoint components, keeps the LARGEST component as the sole
#' geometry for that ID and drops the smaller components.
#'
#' Rationale: `clean_geometry`'s "glue small orphan parts to nearest
#' main_part by shared-boundary length" logic is buggy in the presence of
#' single-part sibling divides, it can crush a single-part polygon to a
#' ~100 m2 residual. For our pipeline, inter-divide topology (overlaps,
#' gaps between IDs) is handled separately by `hfutils::clean_gpkg_layer`
#' (Step 4b + 6d), so we only need the intra-divide "pick largest piece"
#' behavior here.
#'
#' Orphan fragments that are dropped become geometry gaps, which the
#' downstream mapshaper pass absorbs into neighboring divides by
#' `gap_fill_area`.
#'
#' @param divides An `sf` object of divides with a `flowpath_id` (or
#'   equivalent ID) column and polygon geometry.
#' @param id_col Name of the ID column. Default `"flowpath_id"`.
#'
#' @return Input `sf` with MULTIPOLYGON divides reduced to their largest
#'   POLYGON component. Single-part polygons are untouched.
#'
#' @keywords internal
dissolve_multipart_divides <- function(divides, id_col = "flowpath_id") {
  if (!inherits(divides, "sf") || nrow(divides) == 0L) return(divides)

  geom <- sf::st_geometry(divides)
  types <- as.character(sf::st_geometry_type(geom))
  multi_idx <- which(types == "MULTIPOLYGON")
  if (length(multi_idx) == 0L) return(divides)

  crs <- sf::st_crs(geom)
  n_reduced <- 0L
  for (i in multi_idx) {
    mp <- geom[[i]]
    if (length(mp) <= 1L) next
    # Each mp[[j]] is a POLYGON (list of coordinate rings). Pick the one
    # with the largest exterior-ring area via the shoelace formula: avoids
    # an sf::st_area round-trip per component.
    areas <- vapply(mp, function(poly) {
      ring <- poly[[1L]]
      if (nrow(ring) < 4L) return(0)
      x <- ring[, 1L]
      y <- ring[, 2L]
      n <- nrow(ring)
      abs(sum(x[-n] * y[-1L] - x[-1L] * y[-n])) / 2
    }, numeric(1L))
    best <- which.max(areas)
    geom[[i]] <- sf::st_polygon(mp[[best]])
    n_reduced <- n_reduced + 1L
  }
  sf::st_geometry(divides) <- sf::st_sfc(geom, crs = crs)

  if (n_reduced > 0L) message(sprintf(
    "dissolve_multipart_divides: reduced %d multi-part divide(s) to largest component",
    n_reduced))
  divides
}

# =====================================================================
# fill_small_holes  (was R/fill_small_holes.R)
# =====================================================================
#' Fill tiny interior holes in polygons
#'
#' @description
#' Removes interior rings (holes) whose area is below a threshold. Common
#' fix for raster-to-vector artifacts where polygon unioning leaves sub-pixel
#' slivers as holes.
#'
#' @param x An `sf` object with POLYGON or MULTIPOLYGON geometry, or a file
#'   path to a GeoPackage layer.
#' @param threshold_m2 Numeric. Interior rings with area below this threshold
#'   (in square meters, using the object's projected CRS) are removed.
#'   Default `1000` (0.001 sqkm).
#' @param layer If `x` is a file path, the layer name to read/write.
#' @param inplace If `x` is a file path, whether to write the filled result
#'   back to the same layer. Default `TRUE`.
#'
#' @return An `sf` object with tiny holes removed. If `inplace = TRUE` and
#'   `x` is a file path, writes to disk and returns invisibly.
#'
#' @keywords internal
fill_small_holes <- function(x, threshold_m2 = 1000, layer = NULL, inplace = TRUE) {
  read_path <- is.character(x) && length(x) == 1L
  if (read_path) {
    if (is.null(layer)) stop("layer must be specified when x is a file path")
    out_path <- x
    x <- sf::st_read(x, layer, quiet = TRUE)
  }

  ring_area <- function(ring) {
    if (nrow(ring) < 4L) return(0)
    n <- nrow(ring)
    abs(sum(ring[-n, 1L] * ring[-1L, 2L] -
      ring[-1L, 1L] * ring[-n, 2L])) / 2
  }
  filter_poly <- function(p, thr) {
    if (length(p) <= 1L) return(p)
    ha <- vapply(p[-1L], ring_area, numeric(1L))
    c(p[1L], p[-1L][ha >= thr])
  }

  geom <- sf::st_geometry(x)
  n_filled <- 0L
  new_geom <- lapply(geom, function(g) {
    if (length(g) == 0L || sf::st_is_empty(g)) return(g)
    if (inherits(g, "POLYGON")) {
      if (length(g) > 1L) {
        new_p <- filter_poly(unclass(g), threshold_m2)
        if (length(new_p) < length(g)) n_filled <<- n_filled + 1L
        return(sf::st_polygon(new_p))
      }
      g
    } else if (inherits(g, "MULTIPOLYGON")) {
      parts <- unclass(g)
      had_hole <- any(vapply(parts, length, integer(1L)) > 1L)
      new_mp <- lapply(parts, filter_poly, thr = threshold_m2)
      if (had_hole) n_filled <<- n_filled + 1L
      sf::st_multipolygon(new_mp)
    } else g
  })
  sf::st_geometry(x) <- sf::st_sfc(new_geom, crs = sf::st_crs(geom))

  if (n_filled > 0L)
    message(sprintf("fill_small_holes: cleaned %d polygon(s)", n_filled))

  if (read_path && inplace) {
    sf::st_write(x, out_path, layer, delete_layer = TRUE, quiet = TRUE)
    return(invisible(x))
  }
  x
}

# =====================================================================
# remap_poi_flowpath_ids  (was R/remap_poi_flowpath_ids.R)
# =====================================================================
#' Remap POI flowpath_ids to aggregated flowpath_ids
#'
#' @description
#' The aggregate-to-distribution step aggregates reconciled flowlines
#' into aggregated flowpaths but does not update each POI's `flowpath_id`
#' from the pre-aggregation reconciled flowline ID to the post-aggregation
#' aggregated flowpath ID. The result: the `pois` layer's `flowpath_id`
#' points to a reconciled ID (e.g. `6603`) that doesn't exist in the output
#' flowpaths layer; the POI should point to the aggregated flowpath (e.g.
#' `4325`) whose `poi_id` column correctly identifies it.
#'
#' This helper reads the `pois` and `flowpaths` layers from `input_gpkg`
#' and rewrites `pois$flowpath_id` using the aggregated flowpath's
#' `poi_id` column as the authoritative mapping.
#'
#' @param input_gpkg Path to the GeoPackage (must have `pois` and
#'   `flowpaths` layers).
#' @param pois_layer Name of the pois layer. Default `"pois"`.
#' @param flowpath_layer Name of the flowpaths layer. Default `"flowpaths"`.
#'
#' @return Invisibly, integer count of remapped rows.
#'
#' @keywords internal
remap_poi_flowpath_ids <- function(input_gpkg,
                                   pois_layer     = "pois",
                                   flowpath_layer = "flowpaths") {
  if (!hfutils::layer_exists(input_gpkg, pois_layer)) {
    message("remap_poi_flowpath_ids: no '", pois_layer, "' layer -> skipped (0 pois)")
    return(invisible(0L))
  }
  pois <- sf::st_read(input_gpkg, pois_layer,     quiet = TRUE)
  fp   <- sf::st_read(input_gpkg, flowpath_layer, quiet = TRUE)

  if (!"poi_id" %in% names(fp)) {
    .hf_log_warn("flowpaths layer has no poi_id column; skip remap")
    return(invisible(0L))
  }

  poi_to_fp <- stats::setNames(as.character(fp$flowpath_id),
    as.character(fp$poi_id))
  poi_ids   <- as.character(pois$poi_id)
  new_fp    <- poi_to_fp[poi_ids]

  # Only update where we have a match; leave NA where not found
  changed <- !is.na(new_fp) &
    (as.character(pois$flowpath_id) != new_fp |
      is.na(pois$flowpath_id))
  if (!any(changed)) {
    message("remap_poi_flowpath_ids: no changes needed")
    return(invisible(0L))
  }

  pois$flowpath_id <- ifelse(is.na(new_fp), as.character(pois$flowpath_id),
    new_fp)

  pois <- .drop_list_cols(pois)
  sf::st_write(pois, input_gpkg, pois_layer,
    delete_layer = TRUE, quiet = TRUE)
  message(sprintf("remap_poi_flowpath_ids: updated %d POI row(s)",
    sum(changed)))
  invisible(sum(changed))
}

#' Make each lake's flowpaths a contiguous subgraph
#'
#' @description
#' t-route collapses every flowpath sharing a `lake_id` into one waterbody node.
#' That is only valid if those flowpaths are connected: a non-member sitting
#' between two members turns `lake -> gap -> lake` into `node_L -> gap -> node_L`,
#' a cycle, and routing with waterbodies enabled fails on it.
#'
#' Membership is a length-weighted majority of a flowpath's member reaches, a
#' rule that decides each flowpath alone and cannot see the network: a reach
#' mostly outside the lake is rejected even when both its neighbours are members.
#' An earlier closure runs before divide-less leaf flowpaths are demoted, which
#' rewires `flowpath_toid` and can open a gap that was not there before, so this
#' runs once the topology is final.
#'
#' It only ever *adds* a `lake_id` to a flowpath already bracketed by one lake,
#' so it cannot change a fabric with no gaps and cannot move a flowpath between
#' lakes.
#'
#' The stamp is propagated to `flowlines` and `network`, which inherit
#' `lake_id` from their flowpath and would otherwise disagree with it.
#'
#' @param input_gpkg Path to the NextGen GeoPackage.
#' @param max_step Longest run of non-members to bridge. Beyond this the gap is
#'   a genuine river between two arms of one waterbody and is left to be
#'   reported rather than silently absorbed.
#' @param flowpath_layer,nexus_layer,flowline_layer,network_layer Layer names.
#' @return Invisibly, a character vector of the `flowpath_id`s stamped.
#' @keywords internal
close_lake_id_gaps <- function(input_gpkg,
                               max_step       = 5L,
                               flowpath_layer = "flowpaths",
                               nexus_layer    = "nexus",
                               flowline_layer = "flowlines",
                               network_layer  = "network") {
  if (!hfutils::layer_exists(input_gpkg, flowpath_layer) ||
      !hfutils::layer_exists(input_gpkg, nexus_layer)) return(invisible(character(0)))

  fp <- sf::st_read(input_gpkg, flowpath_layer, quiet = TRUE)
  nx <- sf::st_read(input_gpkg, nexus_layer, quiet = TRUE)
  if (!all(c("flowpath_id", "flowpath_toid", "lake_id") %in% names(fp)) ||
      !all(c("nexus_id", "nexus_toid") %in% names(nx))) return(invisible(character(0)))

  lake <- as.character(fp$lake_id)
  lake[is.na(lake) | !nzchar(lake) | lake == "0"] <- NA_character_
  if (!any(!is.na(lake))) return(invisible(character(0)))

  fid <- as.character(fp$flowpath_id)
  n2f <- stats::setNames(as.character(nx$nexus_toid), as.character(nx$nexus_id))
  nxt <- match(unname(n2f[as.character(fp$flowpath_toid)]), fid)

  stamped <- character(0)
  for (i in which(!is.na(lake))) {
    L <- lake[i]; j <- nxt[i]; path <- integer(0)
    for (s in seq_len(max_step)) {
      if (is.na(j)) break
      if (!is.na(lake[j])) {
        if (identical(lake[j], L) && length(path)) {
          lake[path] <- L
          stamped <- c(stamped, fid[path])
        }
        break
      }
      path <- c(path, j); j <- nxt[j]
    }
  }
  if (!length(stamped)) return(invisible(character(0)))
  stamped <- unique(stamped)

  fp$lake_id <- lake
  sf::st_write(.drop_list_cols(fp), input_gpkg, flowpath_layer,
    delete_layer = TRUE, quiet = TRUE)

  fp2lake <- stats::setNames(lake, fid)
  if (hfutils::layer_exists(input_gpkg, flowline_layer)) {
    fl <- sf::st_read(input_gpkg, flowline_layer, quiet = TRUE)
    if ("flowpath_id" %in% names(fl)) {
      fl$lake_id <- unname(fp2lake[as.character(fl$flowpath_id)])
      sf::st_write(.drop_list_cols(fl), input_gpkg, flowline_layer,
        delete_layer = TRUE, quiet = TRUE)
    }
  }
  if (hfutils::layer_exists(input_gpkg, network_layer)) {
    net <- sf::st_drop_geometry(sf::st_read(input_gpkg, network_layer, quiet = TRUE))
    if ("flowpath_id" %in% names(net)) {
      net$lake_id <- unname(fp2lake[as.character(net$flowpath_id)])
      sf::st_write(net, input_gpkg, network_layer, delete_layer = TRUE, quiet = TRUE)
    }
  }
  invisible(stamped)
}
