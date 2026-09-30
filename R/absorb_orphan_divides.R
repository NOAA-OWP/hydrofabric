#' Absorb an unowned divide into its neighbour instead of deleting it
#'
#' When a divide's flowpath is swapped away or cleared, that divide is left
#' unowned. Deleting it removes its polygon from the layer, and with it that
#' land's coverage: on a single CONUS VPU this cost 97 km2, invisible because
#' `divides_area_conserved` compares sums at a 0.5% tolerance and a hole and an
#' overlap removal look alike in a sum.
#'
#' Unowned is not the same as empty. The land still drains somewhere, so the
#' polygon is merged into the neighbouring divide it shares the most boundary
#' with, which conserves coverage while keeping one divide per flowpath. A divide
#' with no neighbour at all is kept rather than dropped, since losing area is
#' worse than carrying an unowned polygon forward.
#'
#' @param div Divides `sf`.
#' @param orphans_mask Logical vector over `div` rows marking unowned divides.
#' @returns A list with `div` (the repaired layer), `n_absorbed` and `n_kept`
#'   (orphans that found no neighbour and were retained).
#' @keywords internal
.absorb_orphan_divides <- function(div, orphans_mask) {
  if (!inherits(div, "sf") || !any(orphans_mask))
    return(list(div = div, n_absorbed = 0L, n_kept = 0L))

  orphan_ix <- which(orphans_mask)
  host_ix   <- which(!orphans_mask)
  if (!length(host_ix))            # nothing to absorb into; keep everything
    return(list(div = div, n_absorbed = 0L, n_kept = length(orphan_ix)))

  g        <- sf::st_geometry(div)
  absorbed <- integer(0)

  for (i in orphan_ix) {
    cand <- tryCatch(
      host_ix[lengths(sf::st_intersects(g[host_ix], g[i])) > 0L],
      error = function(e) integer(0))
    if (!length(cand)) {
      # Fall back to the nearest host: still area-conserving, just not adjacent.
      cand <- tryCatch({
        d <- suppressWarnings(as.numeric(sf::st_distance(g[host_ix], g[i])))
        if (all(is.na(d))) integer(0) else host_ix[which.min(d)]
      }, error = function(e) integer(0))
    }
    if (!length(cand)) next        # keep the orphan; never delete its area

    # Prefer the neighbour sharing the most boundary; that is the divide this
    # land most plausibly belongs to.
    pick <- cand[[1L]]
    if (length(cand) > 1L) {
      shared <- vapply(cand, function(j) tryCatch(
        suppressWarnings(sum(as.numeric(sf::st_length(
          sf::st_intersection(sf::st_boundary(g[j]), sf::st_boundary(g[i])))))),
        error = function(e) 0), numeric(1L))
      shared[is.na(shared)] <- 0
      if (any(shared > 0)) {
        pick <- cand[[which.max(shared)]]
      } else {
        ov <- vapply(cand, function(j) tryCatch(
          suppressWarnings(sum(as.numeric(sf::st_area(
            sf::st_intersection(g[j], g[i]))))),
          error = function(e) 0), numeric(1L))
        ov[is.na(ov)] <- 0
        pick <- cand[[which.max(ov)]]
      }
    }

    merged <- tryCatch(
      suppressWarnings(sf::st_cast(sf::st_union(c(g[pick], g[i])), "MULTIPOLYGON")),
      error = function(e) NULL)
    if (is.null(merged) || !length(merged)) next
    g[pick]  <- merged[[1L]]
    absorbed <- c(absorbed, i)
  }

  if (!length(absorbed))
    return(list(div = div, n_absorbed = 0L, n_kept = length(orphan_ix)))

  # Element assignment leaves the sfc class stale, which breaks the later write.
  g <- sf::st_sfc(g, crs = sf::st_crs(div))
  sf::st_geometry(div) <- g
  if ("areasqkm" %in% names(div))
    div$areasqkm <- suppressWarnings(as.numeric(sf::st_area(div))) / 1e6
  div <- div[-absorbed, , drop = FALSE]

  list(div = div, n_absorbed = length(absorbed),
       n_kept = length(orphan_ix) - length(absorbed))
}
