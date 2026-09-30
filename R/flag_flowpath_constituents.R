#' Flag whether each flowline is a direct constituent of its flowpath's line
#'
#' @description
#' Adds a logical `flowpath_constituent` attribute to the `flowlines` layer. A
#' flowline is a *constituent* (`TRUE`) when its geometry lies on the line that
#' its aggregate flowpath traces from outlet to head: i.e. it is directly part
#' of `flowpath_id`. Member flowlines that enter the flowpath as tributaries (a
#' branch that drains into the stem but whose channel the flowpath line does not
#' trace) are `FALSE`.
#'
#' The test is geometric and matches the flowpath line by construction:
#' [extend_flowpaths_to_mainstem()] builds each flowpath line by merging its
#' stem member reaches, so a constituent's geometry is a sub-segment of the line
#' (coverage ~1) while a tributary branches away (coverage low). A flowline is a
#' constituent when at least `cov_threshold` of its length lies within
#' `tol_m` of its flowpath's line. Every flowpath is guaranteed at least one
#' constituent: if no member clears the threshold (a flowpath still carrying
#' the aggregated original, non-member line), its single best-covered member is
#' flagged.
#'
#' @section Routing contract:
#' OWP routes on the flowlines coincident with flowpaths
#' (`flowpath_constituent == TRUE`), that subset IS the routing network, not a
#' descriptive overlay. It must be closed under `flowline_toid`: a coincident
#' flowline's downstream must resolve to another coincident flowline. Do not fix a
#' routing gap by filtering the mask; dropping a non-coincident connector reach
#' that carries a coincident reach's downstream link would leave the routing set
#' disconnected. A coincident flowline whose `flowline_toid` points to a
#' non-coincident reach is a routing-set gap (detected by the
#' `coincident_flowline_toid_closed` invariant); the remedy is to re-point the
#' routing toid across the intervening non-coincident reaches to the next
#' coincident flowline. `flowpath_constituent == FALSE` reaches are off-network
#' tributary / connector geometry; a flowline's `flowline_toid` may legitimately
#' name a reach in a different flowpath (a tributary confluence).
#'
#' @param input_gpkg Path to a GeoPackage containing `flowpath_layer` and
#'   `flowlines_layer`.
#' @param flowpath_layer,flowlines_layer Layer names. Defaults `"flowpaths"` /
#'   `"flowlines"`.
#' @param tol_m Buffer (m) around the flowpath line within which a flowline
#'   segment counts as "on the line". Default 5.
#' @param cov_threshold Fraction of a flowline's length that must lie on the
#'   line for it to be a constituent. Default 0.5.
#'
#' @return Invisibly, the number of constituent flowlines. Adds
#'   `flowpath_constituent` (logical) to `flowlines_layer`, overwriting it in
#'   place.
#' @examples
#' \dontrun{
#' flag_flowpath_constituents("ngen_vpu01.gpkg")
#' }
#' @export
flag_flowpath_constituents <- function(input_gpkg,
                                       flowpath_layer  = "flowpaths",
                                       flowlines_layer = "flowlines",
                                       tol_m           = 5,
                                       cov_threshold   = 0.5) {
  fp <- sf::st_read(input_gpkg, flowpath_layer,  quiet = TRUE)
  fl <- sf::st_read(input_gpkg, flowlines_layer, quiet = TRUE)
  if (!nrow(fl)) return(invisible(0L))

  fpid <- as.character(fl$flowpath_id)
  if (sf::st_crs(fl) != sf::st_crs(fp)) fl <- sf::st_transform(fl, sf::st_crs(fp))

  fp_ix <- match(fpid, as.character(fp$flowpath_id))
  frac  <- rep(0, nrow(fl))
  v     <- !is.na(fp_ix)
  if (any(v)) {
    # Fraction of each flowline lying within tol_m of its OWN flowpath's line.
    # Vectorised GEOS: paired flowline vs its flowpath-line buffer.
    flg   <- geos::as_geos_geometry(sf::st_geometry(fl)[v])
    fpb   <- geos::geos_buffer(
      geos::as_geos_geometry(sf::st_geometry(fp)[fp_ix[v]]), tol_m)
    inter <- tryCatch(geos::geos_intersection(flg, fpb), error = function(e) NULL)
    if (!is.null(inter)) {
      ll <- suppressWarnings(as.numeric(geos::geos_length(flg)))
      il <- suppressWarnings(as.numeric(geos::geos_length(inter)))
      frac[v] <- ifelse(ll > 0, il / ll, 0)
    }
  }

  constit <- frac >= cov_threshold

  # Guarantee >= 1 constituent per flowpath: a flowpath still carrying
  # the aggregated original (non-member) line can have no member clear the
  # threshold: flag its single best-covered member so the stem is never empty.
  for (rows in split(seq_len(nrow(fl)), fpid)) {
    if (!any(constit[rows])) constit[rows[which.max(frac[rows])]] <- TRUE
  }

  fl$flowpath_constituent <- constit
  sf::st_write(fl, input_gpkg, flowlines_layer, delete_layer = TRUE, quiet = TRUE)
  message(sprintf(
    "flag_flowpath_constituents: %d of %d flowline(s) flagged flowpath_constituent = TRUE",
    sum(constit), nrow(fl)))
  invisible(sum(constit))
}
