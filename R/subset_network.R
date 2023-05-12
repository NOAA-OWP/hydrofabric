#' Subset Hydrofabric Network
#' @param id hydrofabric id
#' @param comid NHDPlusV2 COMID
#' @param hl_id Hydrolocation URI
#' @param network Network Parquet file
#' @param pattern Pattern for distributed VPU based GPKGS
#' @param lyrs layers to extract
#' @param export_gpkg file path to write to
#'
#' @return
#' @export

subset_network = function(id = NULL, 
                          comid = NULL, 
                          hl_id = NULL,
                          network = 'data/conus_net.parquet',
                          pattern = '/Volumes/Transcend/ngen/CONUS-hydrofabric/05_nextgen/nextgen_{vpu}.gpkg',
                          lyrs  = c("divides", "nexus", "flowpaths", "network", "hydrolocations"),
                          export_gpkg = NULL
                         ){
  
  origin = NULL
  
  net = open_dataset(network) %>% 
      select(id, toid, hf_id, hl_uri, hydroseq, vpu) %>% 
      collect()
  
  if(!is.null(comid)){
    origin = filter(net, hf_id == comid) %>% 
      slice_min(hydroseq) %>% 
      pull(id) %>% 
      unique()
  }
  
  if(!is.null(hl_id)){
    origin = filter(net, hl_uri == !!hl_id) %>% 
      slice_max(hydroseq) %>% 
      pull(id) %>% 
      unique()
  }
  
  if(!is.null(id)){ origin = id }
  
  message("Starting from: `",  origin, "`")
  
  #toid = unique(pull(filter(net, id == !!origin), toid))
  
  vpuid = unique(pull(filter(net, id == origin | toid == origin), vpu))
  
  sub_net = select(filter(net, vpu == vpuid), -vpu)

  ids = unique(c(unlist(get_sorted(
    distinct(select(sub_net, id, toid)), outlets = origin
  ))))
  
  gpkg = glue(pattern, vpu = vpuid)
  
  db <- dbConnect(SQLite(), gpkg)
  on.exit(dbDisconnect(db))
  
  hydrofabric = list()
  
  for(j in 1:length(lyrs)){
    hyaggregate_log("INFO", glue("Subsetting: {lyrs[j]} ({j}/{length(lyrs)})"))
    
    crs = st_layers(gpkg)$crs
    
    t = tbl(db, lyrs[j]) %>%
      filter(if_any(contains('id'), ~ . %in% ids)) %>%
      collect()
    
    if(all(!any(is.na(as.character(crs[[j]]))), nrow(t) > 0 )){
      if(any(c("geometry", "geom") %in% names(t))){
        t = st_as_sf(t, crs = crs[[j]])
      } else {
        t = t
      }
    }
    
    if(!is.null(export_gpkg)){
      write_sf(t, export_gpkg, lyrs[j])
    } else {
      hydrofabric[[lyrs[j]]] = t
    }
  }
  
  if(!is.null(export_gpkg)){ 
    export_gpkg = append_style(export_gpkg, layer_names = lyrs)
    return(export_gpkg)
  } else {
    hydrofabric
  }
}

