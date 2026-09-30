# Does a supplied grid actually cover the domain being built?
#
# Grid paths are commonly set as defaults and overridden per domain, so a domain
# whose block omits the override inherits a grid for somewhere else. Nothing
# downstream objects: the reference clips to its own extent, the partition finds
# no valid cells, and the build produces divides from whatever the fallback or
# the floor supplies. Comparing extents up front turns that into a stop.
#
# The reference bbox comes from gpkg_contents rather than the geometry, so this
# costs one small query rather than a read of the flowpaths layer.
# A national reference spans the whole country, so its layer extent is the wrong
# thing to measure a single VPU's grids against: a legitimate per-VPU grid tile
# would cover a few percent of it and be rejected. When a vpuid is supplied and
# the layer carries that column, the extent of just that VPU is taken from the
# rtree index, an index scan, not a geometry read.
.reference_bbox <- function(reference_gpkg, layer = "flowpaths",
                            vpuid = NA_character_) {
  con <- DBI::dbConnect(RSQLite::SQLite(), reference_gpkg)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  q <- tryCatch(DBI::dbGetQuery(con, sprintf(
    "SELECT min_x, min_y, max_x, max_y, srs_id FROM gpkg_contents
      WHERE table_name = '%s'", layer)), error = function(e) NULL)
  if (is.null(q) || !nrow(q) || anyNA(q[1, c("min_x", "min_y", "max_x", "max_y")]))
    return(NULL)
  srs <- tryCatch(DBI::dbGetQuery(con, sprintf(
    "SELECT organization, organization_coordsys_id FROM gpkg_spatial_ref_sys
      WHERE srs_id = %d", q$srs_id[1])), error = function(e) NULL)
  crs <- if (!is.null(srs) && nrow(srs))
    sprintf("%s:%s", srs$organization[1], srs$organization_coordsys_id[1]) else NA_character_
  bbox <- c(q$min_x[1], q$min_y[1], q$max_x[1], q$max_y[1])
  scope <- "reference"

  if (!is.na(vpuid) && nzchar(vpuid)) {
    v <- tryCatch(DBI::dbGetQuery(con, sprintf(
      "SELECT MIN(minx) x1, MIN(miny) y1, MAX(maxx) x2, MAX(maxy) y2
         FROM rtree_%s_geom
        WHERE id IN (SELECT fid FROM \"%s\" WHERE vpuid = '%s')",
      layer, layer, vpuid)), error = function(e) NULL)
    if (!is.null(v) && nrow(v) && !anyNA(unlist(v[1, ]))) {
      bbox  <- c(v$x1[1], v$y1[1], v$x2[1], v$y2[1])
      scope <- sprintf("VPU %s", vpuid)
    }
  }
  list(bbox = bbox, crs = crs, scope = scope)
}

# Fraction of the reference bbox that falls inside the raster's extent.
.grid_coverage <- function(raster_path, ref) {
  if (is.null(ref) || is.na(ref$crs)) return(NA_real_)
  r <- tryCatch(terra::rast(raster_path), error = function(e) NULL)
  if (is.null(r)) return(NA_real_)
  b <- tryCatch({
    p <- sf::st_as_sfc(sf::st_bbox(c(xmin = ref$bbox[1], ymin = ref$bbox[2],
      xmax = ref$bbox[3], ymax = ref$bbox[4]), crs = sf::st_crs(ref$crs)))
    sf::st_bbox(sf::st_transform(p, sf::st_crs(terra::crs(r))))
  }, error = function(e) NULL)
  if (is.null(b)) return(NA_real_)
  e <- as.vector(terra::ext(r))          # xmin, xmax, ymin, ymax
  ow <- max(0, min(b["xmax"], e[2]) - max(b["xmin"], e[1]))
  oh <- max(0, min(b["ymax"], e[4]) - max(b["ymin"], e[3]))
  area <- (b["xmax"] - b["xmin"]) * (b["ymax"] - b["ymin"])
  if (!is.finite(area) || area <= 0) return(NA_real_)
  as.numeric((ow * oh) / area)
}

#' Check that the supplied grids cover the domain
#'
#' @description
#' Compares each raster's extent against the reference fabric's extent. A
#' primary grid that does not cover the domain is an error: the most common
#' cause is a per-domain configuration inheriting a grid meant for somewhere
#' else, which otherwise builds without complaint. A fallback that does not
#' reach the domain is reported and left alone, since it simply never fires.
#'
#' @param reference_gpkg Path to the reference fabric.
#' @param fdr_path,fac_path Primary grid paths.
#' @param fallback_fdr_path,fallback_fac_path Optional fallback grid paths.
#' @param vpuid Optional VPU being built. When supplied and the reference carries
#'   a `vpuid` column, coverage is measured against that VPU's extent rather than
#'   the whole reference, so a per-VPU grid tile is judged on the ground it is
#'   meant to cover.
#' @param min_coverage Minimum share of the reference extent a primary grid must
#'   cover. Default `0.99`.
#' @param strict If `TRUE` (default) an uncovered primary grid raises an error.
#' @returns Invisibly, a named numeric of coverage fractions.
#' @export
check_grid_coverage <- function(reference_gpkg, fdr_path, fac_path,
                               fallback_fdr_path = NULL, fallback_fac_path = NULL,
                               vpuid = NA_character_,
                               min_coverage = 0.99, strict = TRUE) {
  ref <- .reference_bbox(reference_gpkg, vpuid = vpuid)
  if (is.null(ref)) {
    .hf_log_warn("check_grid_coverage: reference extent unavailable; skipped")
    return(invisible(c()))
  }
  prim <- c(fdr = fdr_path, fac = fac_path)
  fall <- c(fallback_fdr = fallback_fdr_path, fallback_fac = fallback_fac_path)
  fall <- fall[!vapply(fall, is.null, logical(1L))]
  cov  <- vapply(c(prim, unlist(fall)), .grid_coverage, numeric(1L), ref = ref)

  for (nm in names(cov)) {
    v <- cov[[nm]]
    lab <- sprintf("%s covers %.1f%% of the %s extent", nm,
                   100 * ifelse(is.na(v), NA_real_, v), ref$scope)
    if (is.na(v)) { message("[grids] ", nm, ": coverage not determined"); next }
    if (nm %in% names(prim) && v < min_coverage) {
      msg <- sprintf(paste(
        "%s. A primary grid must cover the domain (floor %.0f%%).",
        "Check that this domain sets its own fdr_path/fac_path rather than",
        "inheriting another domain's."), lab, 100 * min_coverage)
      if (strict) stop(msg, call. = FALSE) else .hf_log_warn(msg)
    } else {
      message("[grids] ", lab)
    }
  }
  invisible(cov)
}
