mask_hydrofabric = function(gpkg, mask) {
  
  sf::st_delete(gpkg, "layer_styles")
  lyrs = sf::st_layers(gpkg)
  
  s_lyrs  = lyrs$name[!is.na(lyrs$geomtype)]
  as_lyrs = lyrs$name[is.na(lyrs$geomtype)]
  
  hydrofabric = list()
  ids = list()
  
  # Spatial Data Subset
  for (j in 1:length(s_lyrs)) {

    ind = which(lyrs$name %in% s_lyrs[j])
    crs = lyrs$crs[[ind]]
    
    tmp = sf::read_sf(gpkg, s_lyrs[j]) |>
      sf::st_filter(sf::st_transform(mask, crs))
    
    ids[[length(ids) + 1]] = dplyr::select(sf::st_drop_geometry(tmp), dplyr::any_of(c(
      'COMID',  'FEATUREID', 'divide_id', 'id', 'ds_id', "ID"
    ))) |>
      unlist() |>
      unique()
    
    if (!is.null(outfile)) {
      sf::write_sf(tmp, outfile, s_lyrs[j])
    } else {
      hydrofabric[[s_lyrs[j]]] = tmp
    }
  }
  
  # A-spatial Data Subset
  all_ids = unique(do.call("c", ids))
  all_ids = all_ids[!is.na(all_ids)]
  
  for (j in 1:length(as_lyrs)) {
    tmp = sf::read_sf(gpkg, as_lyrs[j]) %>%
      dplyr::filter(dplyr::if_any(dplyr::any_of(c('COMID',  'FEATUREID', 'divide_id', 'id', 'ds_id', "ID"
      )), ~ . %in% all_ids))
    
    if (!is.null(outfile)) {
      sf::write_sf(tmp, outfile, s_lyrs[j])
    } else {
      hydrofabric[[s_lyrs[j]]] = tmp
    }
  }
  
  if (!is.null(outfile)) {
    outfile = append_style(outfile, qml_dir = qml_dir, layer_names = lyrs)
    return(outfile)
  } else {
    return(hydrofabric)
  }
}