#' @title Find a Refernece Feature
#' @param network table from network file default is NULL
#' @param id hydrofabric id (relevant only to nextgen fabrics)
#' @param comid NHDPlusV2 COMID
#' @param hl_id hydrolocation URI (relevant only to nextgen fabrics)
#' @param poi_id POI identifier
#' @param nldi_feature list with names 'featureSource' and 'featureID' where 'featureSource' is derived from the "source" column of the response of dataRetrieval::get_nldi_sources() and the 'featureID' is a known identifier from the specified 'featureSource'.
#' @param xy Location given as vector of XY in CRS 4326 (longitude, latitude)

input_to_reference_feature = function(net, 
                                      id, 
                                      comid, 
                                      hl_id, 
                                      poi_id,
                                      nldi_feature, 
                                      xy){
  
  toid <- divide_id <- hf_hydroseq <- hf_id <- hydroseq <- vpu <- NULL
  
  if (!is.null(id) & !is.null(net)) {
    comid = dplyr::filter(net, id == !!id | toid == !!id | divide_id == !!id) |>
      dplyr::slice_max(hf_hydroseq) |>
      dplyr::pull(hf_id) |>
      unique()
  }
  
  if (!is.null(nldi_feature)) {
    if (length(nldi_feature) == 1) {
      if (is.url(nldi_feature)) {
        comid = sf::read_sf(nldi_feature)$nhdpv2_comid
      }
    } else {
      comid = nhdplusTools::get_nldi_feature(nldi_feature)$comid
    }
  }
  
  if (!is.null(xy)) {
    comid = nhdplusTools::discover_nhdplus_id(point = sf::st_sfc(sf::st_point(c(xy[1], xy[2])), crs = 4326))
  }
  
  if (!is.null(hl_uri) & !is.null(net)) {
    if (is.url(hl_uri)) {
      comid = sf::read_sf(hl_uri)$nhdpv2_comid
    } else {
      origin = dplyr::filter(net, hl_uri == !!hl_uri) |>
        dplyr::slice_max(hf_hydroseq) |>
        dplyr::pull(toid) |>
        unique()
    }
  }
  
  if (!is.null(comid) & !is.null(net)) {
    origin = dplyr::filter(net, hf_id == comid) |>
      dplyr::arrange(hf_hydroseq, hydroseq) |>
      dplyr::slice(1) |>
      dplyr::pull(id) |>
      unique()
  } else if (is.null(net)) {
    origin = comid
  }
  
  if (is.null(origin) | length(origin) > 1) {
    stop("Single origin not found")
    print(origin)
  }
  
  if (is.null(net)) {
    xx = suppressMessages({
      nhdplusTools::get_nhdplus(comid = comid)
    })
    
    v = nhdplusTools::vpu_boundaries
    
    vpuid = v$VPUID[which(lengths(sf::st_intersects(
      sf::st_transform(v, sf::st_crs(xx)), xx
    )) > 0)]
  } else {
    vpuid = unique(dplyr::pull(dplyr::filter(net, id == origin |
                                               toid == origin), vpu))
  }
}
