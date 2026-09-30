#' Re-merge divide parts split by geometry cleaning
#'
#' `mapshaper -clean`, reached through [hfutils::clean_gpkg_layer()], can bisect
#' a narrow divide into disjoint parts. Casting that result to `POLYGON` turns
#' the parts into separate rows sharing one `divide_id`, and keeping only one of
#' them discards the rest. That silently destroys catchment area: across a
#' national build it cost 886 km2, with every VPU still inside the 0.5%
#' conservation tolerance, so nothing flagged it. It also strands the owning
#' flowpath outside what remains of its catchment, which later reads as low
#' coverage and gets the flowpath cleared and then deleted.
#'
#' A bisected divide is a cleaning artifact rather than two catchments: the parts
#' were one divide before cleaning and belong to one flowpath after it. So
#' normalize to `MULTIPOLYGON`, which is lossless and already the layer's type
#' elsewhere, and union any rows sharing a `divide_id` back into one multipart
#' feature.
#'
#' @param d Divides `sf`.
#' @returns A list with `divides` (the repaired layer), `n_merged` (divides whose
#'   parts were rejoined) and `changed` (whether anything was altered).
#' @keywords internal
.merge_bisected_divides <- function(d) {
  none <- list(divides = d, n_merged = 0L, changed = FALSE)
  if (!inherits(d, "sf") || nrow(d) == 0L || !"divide_id" %in% names(d))
    return(none)

  # Lossless type normalization. Casting to POLYGON would be the lossy step.
  recast <- any(sf::st_geometry_type(d) != "MULTIPOLYGON")
  if (recast)
    d <- suppressWarnings(sf::st_cast(d, "MULTIPOLYGON", warn = FALSE))

  ids <- as.character(d$divide_id)
  if (!anyDuplicated(ids))
    return(list(divides = d, n_merged = 0L, changed = recast))

  dup_ids <- unique(ids[duplicated(ids)])
  g       <- sf::st_geometry(d)
  drop    <- integer(0)
  for (did in dup_ids) {
    ix <- which(ids == did)
    u  <- suppressWarnings(sf::st_cast(sf::st_union(g[ix]), "MULTIPOLYGON"))
    g[ix[1L]] <- u[[1L]]
    drop <- c(drop, ix[-1L])
  }
  # Assigning into an sfc element leaves its class stale, which breaks the write.
  g <- sf::st_sfc(g, crs = sf::st_crs(d))
  sf::st_geometry(d) <- g
  if (length(drop)) d <- d[-drop, , drop = FALSE]

  list(divides = d, n_merged = length(dup_ids), changed = TRUE)
}
