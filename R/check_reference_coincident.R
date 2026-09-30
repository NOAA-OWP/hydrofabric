#' Find reference reaches drawn on top of one another
#'
#' @description
#' Two reaches that occupy the same ground are indistinguishable to anyone
#' reading the fabric: one click returns both, and they appear to be one channel
#' counted twice. The pipeline carries such geometry through faithfully, so the
#' only place to catch it is the input.
#'
#' The failure is not exact duplication, which is rare and easy to test for. It
#' is near-coincidence: reaches running within a metre of each other for
#' kilometres without sharing a vertex. A strict topological overlap test sees
#' almost nothing there, while a mouse click cannot tell the difference. So this
#' measures the fraction of each reach that lies within `tolerance_m` of a
#' different reach.
#'
#' Observed on the CONUS reference at a confluence in VPU 04: four reaches
#' terminating at one point with the same bearing, sharing 2.3 to 5.5 km with
#' each other, three of them from a preliminary source and the fourth a
#' 522,054 km2 river drawn as a single straight 9.5 km segment. Coincident
#' geometry also misleads aggregation, which reads spatial adjacency: those four
#' tributaries were contracted into one flowpath carrying four outlets.
#'
#' Report-only, and deliberately not wired into `run_pipeline()`. The geometry
#' belongs to the reference and is not ours to change, and this asks a question
#' about a reference release rather than about a VPU, so it runs once per input
#' rather than once per build.
#'
#' @param reference_gpkg Path to the reference GeoPackage.
#' @param layer Flowline layer name.
#' @param id Identifier column.
#' @param vpuid Optional VPU to restrict to, when the layer carries `vpuid`.
#' @param tolerance_m Distance within which two reaches count as coincident.
#'   1 m is three orders of magnitude below the 30 m DEM cell, so it cannot
#'   catch channels that are merely close.
#' @param min_frac Fraction of a reach that must lie within `tolerance_m` of
#'   another before it is reported.
#' @param sample Number of reaches to test, or `NULL` for all. The default keeps
#'   the cost bounded on a national reference; the reported rate is an estimate
#'   and the message says so.
#' @param strict If `TRUE`, warn when any are found.
#' @returns Invisibly, a data frame of `id`, `frac`, `length_m` for the reaches
#'   found, most coincident first.
#' @export
check_reference_coincident_geometry <- function(reference_gpkg,
                                                layer  = "flowpaths",
                                                id     = "flowpath_id",
                                                vpuid  = NA_character_,
                                                tolerance_m = 1,
                                                min_frac = 0.5,
                                                sample = 5000L,
                                                strict = FALSE) {
  stopifnot(is.character(reference_gpkg), length(reference_gpkg) == 1L,
    file.exists(reference_gpkg))

  q <- if (!is.na(vpuid))
    sprintf('SELECT "%s", geom FROM "%s" WHERE vpuid = \'%s\'', id, layer, vpuid)
  else
    sprintf('SELECT "%s", geom FROM "%s"', id, layer)

  fl <- tryCatch(sf::st_read(reference_gpkg, query = q, quiet = TRUE),
                 error = function(e) NULL)
  if (is.null(fl) || !nrow(fl)) {
    message("check_reference_coincident_geometry: no features read; skipped")
    return(invisible(data.frame()))
  }

  # Distances are metres, so work in an equal-area projection when the layer is
  # geographic. Conus Albers is the reference's own CRS for CONUS.
  if (isTRUE(sf::st_is_longlat(fl))) fl <- sf::st_transform(fl, 5070)
  g <- sf::st_geometry(fl)
  n <- length(g)

  ix <- if (is.null(sample) || sample >= n) seq_len(n) else
    sort(sample.int(n, as.integer(sample)))

  own  <- as.numeric(sf::st_length(g[ix]))
  cand <- sf::st_is_within_distance(g[ix], g, dist = tolerance_m, sparse = TRUE)

  frac <- vapply(seq_along(ix), function(k) {
    js <- setdiff(cand[[k]], ix[k])
    if (!length(js) || !is.finite(own[k]) || own[k] <= 0) return(0)
    inb <- tryCatch(suppressWarnings(sf::st_intersection(
      g[ix[k]], sf::st_buffer(sf::st_union(g[js]), tolerance_m))),
      error = function(e) NULL)
    if (is.null(inb) || !length(inb)) return(0)
    sum(as.numeric(sf::st_length(inb))) / own[k]
  }, numeric(1L))

  hit <- which(frac >= min_frac)
  out <- data.frame(
    id       = as.character(sf::st_drop_geometry(fl)[[id]])[ix][hit],
    frac     = round(frac[hit], 4),
    length_m = round(own[hit]),
    stringsAsFactors = FALSE)
  out <- out[order(-out$frac), , drop = FALSE]

  est <- if (length(ix) < n) sprintf(", est. %d of %d in the layer",
    round(n * length(hit) / length(ix)), n) else ""
  msg <- sprintf(paste0("check_reference_coincident_geometry: %d of %d reach(es) ",
    "tested have at least %.0f%% of their length within %g m of a different ",
    "reach%s"), length(hit), length(ix), 100 * min_frac, tolerance_m, est)
  message(msg)
  if (nrow(out))
    message(sprintf("  most coincident: %s",
      paste(sprintf("%s (%.0f%% of %.0f m)", utils::head(out$id, 3),
        100 * utils::head(out$frac, 3), utils::head(out$length_m, 3)),
        collapse = ", ")))
  if (strict && nrow(out))
    warning(msg, call. = FALSE)

  invisible(out)
}
