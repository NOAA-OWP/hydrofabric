#' Find reference reaches that terminate on a channel that carries on
#'
#' @description
#' A reach with no downstream becomes a terminal nexus, and flow arriving there
#' is not routed on. That is right for an outlet: a coast, a closed basin. It is
#' wrong where the channel plainly continues, and a routing consumer sees the
#' flow simply stop.
#'
#' The test is geometric because the attributes cannot answer it. `terminalfl`
#' is unpopulated for whole sources, and `dnhydroseq` is stored too narrowly to
#' hold the values it points at. What does distinguish the two is what lies at
#' the reach's downstream end: nothing, or another reach carrying more water
#' than the one that just stopped.
#'
#' Observed on the CONUS reference: of 15,877 terminal reaches, 822 end within a
#' DEM cell of another reach, and 446 of those end on one that drains more,
#' including a 324,094 km2 river in VPU 13 that stops 30 m short of its own
#' continuation and the Saint John at Fort Kent, where 12,334 km2 stops 13 m from
#' a reach carrying 12,363 km2. They are spread nationwide rather than gathered
#' at the borders: VPU 17 has 76, VPU 18 69, VPU 16 68.
#'
#' Report-only, and deliberately not wired into `run_pipeline()`. The topology
#' belongs to the reference, and repairing it means writing edges into an input
#' rather than working around it in a build. This reports what a staged
#' reference still carries.
#'
#' @param reference_gpkg Path to the reference GeoPackage.
#' @param layer Flowline layer name.
#' @param id,toid Identifier and downstream-identifier columns.
#' @param area Column holding total upstream drainage area, used to tell a
#'   receiver from a neighbouring headwater.
#' @param vpuid Optional VPU to restrict to, when the layer carries `vpuid`.
#' @param tolerance_m Distance within which a reach counts as ending on another.
#'   The default is the 30 m DEM cell. A wider radius reports unrelated channel
#'   heads: at 100 m VPU 01 scores 17 and at 5 m it scores none.
#' @param strict If `TRUE`, warn when any are found.
#' @returns Invisibly, a data frame of `id`, `area`, `continues_into` and
#'   `into_area`, largest first.
#' @export
check_reference_severed_terminals <- function(reference_gpkg,
                                              layer = "flowpaths",
                                              id    = "flowpath_id",
                                              toid  = "flowpath_toid",
                                              area  = "totdasqkm",
                                              vpuid = NA_character_,
                                              tolerance_m = 30,
                                              strict = FALSE) {
  stopifnot(is.character(reference_gpkg), length(reference_gpkg) == 1L,
    file.exists(reference_gpkg))

  q <- if (!is.na(vpuid))
    sprintf('SELECT "%s", "%s", "%s", geom FROM "%s" WHERE vpuid = \'%s\'',
      id, toid, area, layer, vpuid)
  else
    sprintf('SELECT "%s", "%s", "%s", geom FROM "%s"', id, toid, area, layer)
  f <- tryCatch(sf::st_read(reference_gpkg, query = q, quiet = TRUE),
                error = function(e) NULL)
  if (is.null(f) || !nrow(f)) {
    message("check_reference_severed_terminals: no features read; skipped")
    return(invisible(data.frame()))
  }
  if (isTRUE(sf::st_is_longlat(f))) f <- sf::st_transform(f, 5070)

  ids <- as.character(f[[id]])
  tos <- as.character(f[[toid]])
  da  <- suppressWarnings(as.numeric(f[[area]]))
  term <- which(is.na(tos) | tos == "0" | !nzchar(tos))
  if (!length(term)) {
    message("check_reference_severed_terminals: no terminal reaches; skipped")
    return(invisible(data.frame()))
  }

  # The downstream end is where the channel is cut.
  ends <- sf::st_sfc(lapply(sf::st_geometry(f)[term], function(x) {
    m <- unclass(x); while (is.list(m)) m <- m[[length(m)]]
    sf::st_point(m[nrow(m), 1:2])
  }), crs = sf::st_crs(f))
  near <- suppressWarnings(
    sf::st_is_within_distance(ends, sf::st_geometry(f), dist = tolerance_m))

  hit <- lapply(seq_along(term), function(k) {
    i  <- term[k]
    js <- setdiff(near[[k]], i)
    # A receiver drains more than the reach handed to it; a neighbour that
    # drains less is a nearby headwater, not somewhere the water was going.
    js <- js[is.finite(da[js]) & is.finite(da[i]) & da[js] > da[i]]
    if (!length(js)) return(NULL)
    j <- js[which.max(da[js])]
    data.frame(id = ids[i], area = da[i], continues_into = ids[j],
      into_area = da[j], stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, hit)
  if (is.null(out)) out <- data.frame(id = character(0), area = numeric(0),
    continues_into = character(0), into_area = numeric(0),
    stringsAsFactors = FALSE)
  out <- out[order(-out$area), , drop = FALSE]

  message(sprintf(paste0("check_reference_severed_terminals: %d of %d terminal ",
    "reach(es) end within %g m of a channel that drains more"),
    nrow(out), length(term), tolerance_m))
  if (nrow(out))
    message(sprintf("  largest: %s",
      paste(sprintf("%s (%.0f km2) -> %s (%.0f km2)", utils::head(out$id, 3),
        utils::head(out$area, 3), utils::head(out$continues_into, 3),
        utils::head(out$into_area, 3)), collapse = ", ")))
  if (strict && nrow(out))
    warning("reference reaches terminate on a channel that continues",
      call. = FALSE)

  invisible(out)
}
