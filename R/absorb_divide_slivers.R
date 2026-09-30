#' Absorb a divide's sliver parts into the neighbour they share the most boundary with
#'
#' @description
#' `construct_agg_divides` unions each aggregate's member catchments, so an
#' aggregate whose members are not spatially contiguous gets a genuinely
#' multi-part divide. Most of that multi-part-ness is not a real disjoint
#' catchment: on CONUS, 1,192 of 1,300 multi-part divides have a second part
#' under 1% of the divide (median 900 m2). Those are slivers left by clipping,
#' snapping and union arithmetic.
#'
#' Nothing downstream reduces them. `dissolve_multipart_divides` runs in
#' `reconcile_divides`, which is BEFORE the aggregation that creates the
#' fragmentation, and `hfutils::clean_gpkg_layer(gap_fill_area = )` fills holes
#' BETWEEN divides, a divide's own second part is an owned polygon, not a gap,
#' so it is never a gap-fill candidate.
#'
#' This pass moves each qualifying sliver into the adjacent divide it shares the
#' most boundary with, which is the divide that land most plausibly drains to.
#' It is area-conserving by construction: the polygon changes owner, it is never
#' dropped. That distinction matters, the dissolve-and-drop approach used at
#' the refactor phase would forfeit 285 km2 nationally, and for a genuinely disjoint divide
#' it would discard a real catchment.
#'
#' Deliberately conservative. A part must be BOTH absolutely small
#' (`min_part_km2`) and small relative to its own divide (`max_part_frac`)
#' before it is eligible, so the ~53 divides whose second part is a substantial
#' fraction of a real catchment keep their multi-part geometry. A sliver with no
#' adjacent divide is kept as well: carrying a multi-part divide is better than
#' losing the land.
#'
#' @param div Divides `sf`.
#' @param min_part_km2 Absolute ceiling on an absorbable part, km2. Default 0.01
#'   (1 hectare), which covers 1,194 of the 1,300 CONUS cases.
#' @param max_part_frac Ceiling on an absorbable part as a fraction of its own
#'   divide's total area. Default 0.01.
#' @param min_shared_m Minimum shared boundary with the recipient, metres. A
#'   part touching only at a point would union into a still-disjoint polygon.
#' @returns A list with `div` (repaired layer), `n_absorbed`, `n_kept` (eligible
#'   slivers with no usable neighbour), `area_moved_km2`, and `n_multi_before` /
#'   `n_multi_after`.
#' @keywords internal
.absorb_divide_slivers <- function(div, min_part_km2 = 0.01,
                                   max_part_frac = 0.01, min_shared_m = 1,
                                   flowlines = NULL) {
  nil <- list(div = div, n_absorbed = 0L, n_kept = 0L, area_moved_km2 = 0,
              n_multi_before = 0L, n_multi_after = 0L, n_seam_dissolved = 0L)
  if (!inherits(div, "sf") || !nrow(div)) return(nil)

  # Thresholds are metric and the recipient test buffers by metres, so do the
  # geometry in an equal-area projection and hand the result back in the layer's
  # own CRS. The ngen divides ship in EPSG:4326, where a bare `1` would be a
  # degree.
  in_crs <- sf::st_crs(div)
  geographic <- isTRUE(sf::st_is_longlat(div))
  if (geographic) div <- sf::st_transform(div, 5070)

  g  <- sf::st_geometry(div)
  np <- vapply(g, function(x) if (is.null(x)) 0L else length(x), integer(1L))
  multi <- which(np > 1L)
  nil$n_multi_before <- length(multi)
  nil$n_multi_after  <- length(multi)
  if (!length(multi)) return(nil)

  # HAIRLINE SEAMS FIRST. Half of the multi-part divides are not disjoint at all:
  # their parts touch at 0.000 m but meet at a point or along a zero-width seam
  # rather than sharing an edge, so st_union leaves them as separate parts. They
  # are one polygon and should be stored as one. A 1 cm close (buffer out, union,
  # buffer back) merges exactly those and nothing else: 1 cm is three orders of
  # magnitude below the 30 m DEM cell, so it cannot bridge a real separation. The
  # area guard makes that a fact rather than an assumption: the result is kept
  # only if it has fewer parts AND the same area.
  n_seam <- 0L
  for (i in multi) {
    one <- sf::st_sfc(g[[i]], crs = sf::st_crs(div))
    a0  <- as.numeric(sf::st_area(one))
    fix <- tryCatch(
      sf::st_buffer(sf::st_union(sf::st_buffer(sf::st_make_valid(one), 0.01)), -0.01),
      error = function(e) NULL)
    if (is.null(fix) || !length(fix)) next
    np_new <- length(sf::st_geometry(fix)[[1]])
    a1     <- as.numeric(sf::st_area(fix))
    if (np_new < length(g[[i]]) && abs(a1 - a0) <= 1e-6 * max(a0, 1)) {
      g[i] <- sf::st_cast(fix, "MULTIPOLYGON")
      n_seam <- n_seam + 1L
    }
  }
  if (n_seam) {
    sf::st_geometry(div) <- g
    np <- vapply(g, function(x) if (is.null(x)) 0L else length(x), integer(1L))
    multi <- which(np > 1L)
    nil$n_multi_after <- length(multi)
    nil$div <- if (geographic) sf::st_transform(div, in_crs) else div
    nil$n_seam_dissolved <- n_seam
    if (!length(multi)) return(nil)
  }

  crs <- sf::st_crs(div)
  .km2 <- function(x) as.numeric(sf::st_area(x)) / 1e6

  # Explode the multi-part rows into single polygons, tagged with their owner row.
  parts <- list(); owner <- integer(0)
  for (i in multi) {
    for (k in seq_along(g[[i]])) {
      parts[[length(parts) + 1L]] <- sf::st_polygon(g[[i]][[k]])
      owner <- c(owner, i)
    }
  }
  psfc <- sf::st_sfc(parts, crs = crs)
  parea <- .km2(psfc)

  # Largest part of each row is the keeper; the rest are absorption candidates.
  keep <- rep(FALSE, length(owner))
  for (i in unique(owner)) {
    ix <- which(owner == i)
    keep[ix[which.max(parea[ix])]] <- TRUE
  }
  row_tot <- vapply(unique(owner), function(i) sum(parea[owner == i]), numeric(1L))
  names(row_tot) <- as.character(unique(owner))
  frac <- parea / unname(row_tot[as.character(owner)])

  cand <- which(!keep & parea < min_part_km2 & frac < max_part_frac)

  # A part carrying its own flowpath's flowline is not a sliver: it is the piece
  # of the catchment the reach actually runs through, and moving it leaves the
  # divide with no member flowline at all. On Alaska that turned five conforming
  # divides into gate failures, every one of them because the absorbed part was
  # the one holding the line. Judged within the fabric's 30 m tolerance, the same
  # distance divide_has_member_flowline uses, so the two agree by construction.
  if (length(cand) && !is.null(flowlines) && inherits(flowlines, "sf") &&
      all(c("flowpath_id") %in% names(flowlines)) &&
      "flowpath_id" %in% names(div)) {
    fl <- sf::st_transform(flowlines, crs)
    fl_fp <- as.character(sf::st_drop_geometry(fl)$flowpath_id)
    holds_line <- vapply(cand, function(j) {
      fp <- as.character(div$flowpath_id[owner[j]])
      if (is.na(fp) || !nzchar(fp)) return(FALSE)
      own <- sf::st_geometry(fl)[fl_fp == fp]
      if (!length(own)) return(FALSE)
      near <- function(x) any(suppressWarnings(lengths(
        sf::st_intersects(own, sf::st_buffer(x, 30)))) > 0L)
      if (!near(psfc[j])) return(FALSE)
      # The part is only worth protecting if moving it would leave the divide
      # with NO member flowline, which is what the Alaska failures were. The
      # keeper part stays by construction, so when that one also holds a line
      # the divide keeps one either way and the sliver is free to move. Asking
      # only whether the sliver is near a line protected far more than the
      # invariant needs: a 30 m buffer reaches lines lying in the main part, and
      # it left three VPU 01 divides multi-part whose keepers held three, three
      # and eight flowlines between them.
      keeper <- which(owner == owner[j] & keep)
      if (!length(keeper)) return(TRUE)
      !near(psfc[keeper[1L]])
    }, logical(1L))
    cand <- cand[!holds_line]
  }

  if (!length(cand)) return(nil)

  # Pick a recipient per sliver: an adjacent DIFFERENT divide, preferring the one
  # sharing the most boundary, the same rule .absorb_orphan_divides uses.
  recipient <- rep(NA_integer_, length(cand))
  for (n in seq_along(cand)) {
    j  <- cand[n]
    pg <- psfc[j]
    nb <- tryCatch(which(lengths(sf::st_intersects(g, sf::st_buffer(pg, 1))) > 0L),
                   error = function(e) integer(0))
    nb <- setdiff(nb, owner[j])
    if (!length(nb)) next
    shared <- vapply(nb, function(q) tryCatch(suppressWarnings(sum(as.numeric(
      sf::st_length(sf::st_intersection(sf::st_boundary(g[q]), sf::st_boundary(pg)))))),
      error = function(e) 0), numeric(1L))
    shared[is.na(shared)] <- 0
    if (!any(shared >= min_shared_m)) next
    recipient[n] <- nb[[which.max(shared)]]
  }

  moved <- which(!is.na(recipient))
  if (!length(moved)) {
    nil$n_kept <- length(cand)
    return(nil)
  }

  # Rebuild only the rows that changed. Donors keep the parts that stayed;
  # recipients gain the slivers assigned to them.
  gone <- cand[moved]
  newg <- g
  for (i in unique(owner[gone])) {                       # donors
    stay <- which(owner == i & !(seq_along(owner) %in% gone))
    if (!length(stay)) next                              # cannot happen: keeper stays
    newg[i] <- sf::st_union(psfc[stay])
  }
  for (r in unique(recipient[moved])) {                  # recipients
    add <- gone[recipient[moved] == r]
    newg[r] <- sf::st_union(c(sf::st_sfc(newg[[r]], crs = crs), psfc[add]))
  }
  newg <- suppressWarnings(sf::st_cast(newg, "MULTIPOLYGON"))
  sf::st_geometry(div) <- newg
  if (geographic) div <- sf::st_transform(div, in_crs)

  np2 <- vapply(sf::st_geometry(div),
                function(x) if (is.null(x)) 0L else length(x), integer(1L))
  list(div = div, n_absorbed = length(moved),
       n_kept = length(cand) - length(moved),
       area_moved_km2 = sum(parea[gone]),
       n_multi_before = nil$n_multi_before, n_multi_after = sum(np2 > 1L),
       n_seam_dissolved = n_seam)
}
