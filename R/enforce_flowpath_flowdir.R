#' Enforce upstream->downstream orientation on aggregated flowpath geometry
#'
#' @description
#' [extend_flowpaths_to_mainstem()] rebuilds each aggregate flowpath line with
#' `geos_line_merge`, which does not preserve flow direction, and (unlike the
#' reconciled `flowlines`, guarded by the `reconciled_flow_direction` invariant)
#' the aggregated `flowpaths` layer has no orientation check, so a fraction of
#' flowpaths ship reversed (their geometry runs downstream->upstream). This pass
#' re-orients every flowpath so its OUTLET is the LAST vertex.
#'
#' The true outlet is located from the reliably-oriented `flowlines` layer: each
#' flowpath's outlet sub-flowline is the member whose `flowline_toid` leaves the
#' flowpath, and its downstream (end) vertex is the flowpath outlet. A flowpath
#' whose START is nearer that outlet than its END is reversed and is flipped with
#' [sf::st_reverse()]. Genuinely branching multi-part geometry (no single
#' meaningful direction) is left untouched.
#'
#' @param gpkg Path to a GeoPackage containing `flowpath_layer` and
#'   `flowlines_layer`.
#' @param flowpath_layer,flowlines_layer Layer names. Defaults `"flowpaths"` /
#'   `"flowlines"`.
#'
#' @return Invisibly, the number of flowpaths reversed. Overwrites
#'   `flowpath_layer` in place.
#' @examples
#' \dontrun{
#' enforce_flowpath_flowdir("ngen_vpu01.gpkg")
#' }
#' @export
enforce_flowpath_flowdir <- function(gpkg,
                                     flowpath_layer  = "flowpaths",
                                     flowlines_layer = "flowlines") {
  fp <- sf::st_read(gpkg, flowpath_layer, quiet = TRUE)
  if (!nrow(fp)) return(invisible(0L))
  if (!hfutils::layer_exists(gpkg, flowlines_layer)) {
    .hf_log_warn("enforce_flowpath_flowdir: no flowlines layer; skipped")
    return(invisible(0L))
  }
  fl <- sf::st_read(gpkg, flowlines_layer, quiet = TRUE)
  fl$flowline_id   <- as.character(fl$flowline_id)
  fl$flowline_toid <- as.character(fl$flowline_toid)
  fl$flowpath_id   <- as.character(fl$flowpath_id)
  if (sf::st_crs(fl) != sf::st_crs(fp)) fl <- sf::st_transform(fl, sf::st_crs(fp))

  fpids <- as.character(fp$flowpath_id)
  fl_by <- split(seq_len(nrow(fl)), fl$flowpath_id)
  tda   <- suppressWarnings(as.numeric(fl$tot_drainage_areasqkm)); tda[is.na(tda)] <- 0

  # Downstream (end) vertex of every flowline: flowlines are single LINESTRINGs
  # with enforced up->down orientation, so the end vertex is the true outlet node.
  fl_end <- sf::st_coordinates(lwgeom::st_endpoint(sf::st_geometry(fl)))[, 1:2, drop = FALSE]

  # Flowpath endpoints from a merged LINESTRING copy (lwgeom segfaults on
  # multipart). Genuinely branching geometry stays MULTILINESTRING and is skipped.
  gm    <- sf::st_line_merge(
    sf::st_cast(sf::st_geometry(fp), "MULTILINESTRING", warn = FALSE))
  is_ls <- as.character(sf::st_geometry_type(gm)) == "LINESTRING" &
    as.numeric(sf::st_length(gm)) > 0
  sp <- matrix(NA_real_, nrow(fp), 2L); ep <- sp
  if (any(is_ls)) {
    sp[is_ls, ] <- sf::st_coordinates(lwgeom::st_startpoint(gm[is_ls]))[, 1:2, drop = FALSE]
    ep[is_ls, ] <- sf::st_coordinates(lwgeom::st_endpoint(gm[is_ls]))[, 1:2, drop = FALSE]
  }

  geom   <- sf::st_geometry(fp)
  rev_ct <- 0L
  for (k in which(is_ls)) {
    r <- fl_by[[fpids[k]]]; if (is.null(r)) next
    ids  <- fl$flowline_id[r]
    o    <- r[!(fl$flowline_toid[r] %in% ids)]        # sub whose toid leaves the fp
    orow <- if (length(o)) o[which.max(tda[o])] else r[which.max(tda[r])]
    ov   <- fl_end[orow, ]
    d_start <- sum((sp[k, ] - ov)^2)
    d_end   <- sum((ep[k, ] - ov)^2)
    if (is.finite(d_start) && is.finite(d_end) && d_start < d_end) {
      geom[k] <- sf::st_reverse(geom[k])
      rev_ct  <- rev_ct + 1L
    }
  }
  sf::st_geometry(fp) <- geom
  message(sprintf(
    "enforce_flowpath_flowdir: reversed %d of %d flowpath(s) to outlet-last",
    rev_ct, nrow(fp)))
  sf::st_write(fp, gpkg, flowpath_layer, delete_layer = TRUE, quiet = TRUE)
  invisible(rev_ct)
}
