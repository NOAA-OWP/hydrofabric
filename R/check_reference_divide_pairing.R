#' Find reference reaches that end outside their own divide
#'
#' @description
#' A reach and its divide are a pair: the reach drains the divide, and its
#' downstream end is the divide's outlet. Where the pair is wrong the two
#' describe different places, and nothing downstream can reconcile them.
#'
#' The splitter is where it shows. A parent divide is partitioned by seeding a
#' masked D8 flood at each member's downstream endpoint, and the mask is the
#' parent buffered by 250 m. An endpoint beyond that buffer has no grid cell to
#' seed from, the partition is abandoned, and the whole parent passes to the
#' downstream-most child while its siblings ship owning no land. Widening the
#' mask is not the cure: it would let every parent's flood escape its own
#' catchment.
#'
#' Most endpoints sit marginally outside their divide because an outlet lies on
#' the boundary, so the measure is distance rather than containment. Observed on
#' the CONUS reference for VPU 04: of 150,744 paired reaches, 104,830 endpoints
#' are within 30 m of their divide, and 147 (0.098%) are beyond 250 m, the worst
#' 10.4 km away.
#'
#' Report-only, and deliberately not wired into `run_pipeline()`. The pairing
#' belongs to the reference and is not ours to change, and this asks a question
#' about a reference release rather than about a VPU.
#'
#' @param reference_gpkg Path to the reference GeoPackage.
#' @param flowpath_layer,divide_layer Layer names.
#' @param id Identifier column, present in both layers.
#' @param vpuid Optional VPU to restrict to, when the layers carry `vpuid`.
#' @param tolerance_m Distance beyond which a pairing is reported. The default
#'   is the 250 m buffer the splitter masks with, so it is the distance at which
#'   the pairing stops being usable rather than merely imprecise.
#' @param sample Number of pairs to test, or `NULL` for all.
#' @param strict If `TRUE`, warn when any are found.
#' @returns Invisibly, a data frame of `id`, `dist_m` for the reaches found,
#'   furthest first.
#' @export
check_reference_divide_pairing <- function(reference_gpkg,
                                           flowpath_layer = "flowpaths",
                                           divide_layer   = "divides",
                                           id             = "flowpath_id",
                                           vpuid          = NA_character_,
                                           tolerance_m    = 250,
                                           sample         = NULL,
                                           strict         = FALSE) {
  stopifnot(is.character(reference_gpkg), length(reference_gpkg) == 1L,
    file.exists(reference_gpkg))

  rd <- function(layer) {
    q <- if (!is.na(vpuid))
      sprintf('SELECT "%s", geom FROM "%s" WHERE vpuid = \'%s\'', id, layer, vpuid)
    else
      sprintf('SELECT "%s", geom FROM "%s"', id, layer)
    tryCatch(sf::st_read(reference_gpkg, query = q, quiet = TRUE),
             error = function(e) NULL)
  }
  fl <- rd(flowpath_layer)
  dv <- rd(divide_layer)
  if (is.null(fl) || is.null(dv) || !nrow(fl) || !nrow(dv)) {
    message("check_reference_divide_pairing: no features read; skipped")
    return(invisible(data.frame()))
  }

  # Distances are metres, so work in an equal-area projection when the layers
  # are geographic. Conus Albers is the reference's own CRS for CONUS.
  if (isTRUE(sf::st_is_longlat(fl))) fl <- sf::st_transform(fl, 5070)
  if (!identical(sf::st_crs(dv), sf::st_crs(fl))) dv <- sf::st_transform(dv, sf::st_crs(fl))

  fid <- as.character(sf::st_drop_geometry(fl)[[id]])
  did <- as.character(sf::st_drop_geometry(dv)[[id]])
  j   <- match(fid, did)
  ix  <- which(!is.na(j))
  if (!length(ix)) {
    message("check_reference_divide_pairing: no reach shares an id with a divide; skipped")
    return(invisible(data.frame()))
  }
  if (!is.null(sample) && sample < length(ix))
    ix <- ix[sort(sample.int(length(ix), as.integer(sample)))]

  fg <- sf::st_geometry(fl)[ix]
  dg <- sf::st_geometry(dv)[j[ix]]

  # The downstream end is the last vertex of the last part, which is where the
  # splitter seeds from.
  ep <- sf::st_sfc(lapply(fg, function(g) {
    m <- unclass(g)
    while (is.list(m)) m <- m[[length(m)]]
    sf::st_point(m[nrow(m), 1:2])
  }), crs = sf::st_crs(fl))

  # Paired, not pairwise: st_distance without by_element compares every reach to
  # every divide, which answers a different and useless question.
  d <- suppressWarnings(as.numeric(sf::st_distance(ep, dg, by_element = TRUE)))

  hit <- which(is.finite(d) & d > tolerance_m)
  out <- data.frame(id = fid[ix][hit], dist_m = round(d[hit]),
    stringsAsFactors = FALSE)
  out <- out[order(-out$dist_m), , drop = FALSE]

  message(sprintf(paste0("check_reference_divide_pairing: %d of %d paired ",
    "reach(es) end more than %g m outside their own divide (%.3f%%); %d end ",
    "within 30 m of it"), nrow(out), length(ix), tolerance_m,
    100 * nrow(out) / length(ix), sum(d <= 30, na.rm = TRUE)))
  if (nrow(out))
    message(sprintf("  furthest: %s",
      paste(sprintf("%s (%.0f m)", utils::head(out$id, 3),
        utils::head(out$dist_m, 3)), collapse = ", ")))
  if (strict && nrow(out))
    warning("reference reaches end outside their own divide", call. = FALSE)

  invisible(out)
}
