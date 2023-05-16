#' Subset Hydrofabric Network
#' @param id hydrofabric id
#' @param comid NHDPlusV2 COMID
#' @param hl_id Hydrolocation URI
#' @param network Network Parquet file
#' @param pattern Pattern for distributed VPU based GPKGS
#' @param lyrs layers to extract
#' @param export_gpkg file path to write to
#'
#' @return file path or list 
#' @export

subset_reference = function(nldi_feature,
                            gpkg = NULL,
                            pattern = '/Volumes/Transcend/ngen/CONUS-hydrofabric/01_reference/reference_{vpu}.gpkg',
                            export_gpkg = NULL) {
  
  feat =  get_nldi_feature(nldi_feature)
  vpu =  st_join(feat, st_transform(vpu_boundaries, st_crs(feat))) %>% 
    pull(VPUID)
  
  if(is.null(gpkg)){
    gpkg = glue(pattern, vpu = vpu)
  }
  

  flowpaths = grep("flowpath|flowline", st_layers(gpkg)$name, value = TRUE)
  flowpaths = flowpaths[!grepl("attributes|edge_list", flowpaths)]
  
  if (length(flowpaths) > 1) {
    stop("Multiple flowpath names found.")
  }
  
  catchments = grep("divide|catchment", st_layers(gpkg)$name, value = TRUE)
  catchments = catchments[!grepl("network", catchments)]
  if (length(catchments) > 1) {
    stop("Multiple catchment names found.")
  }
  
  UT_COMIDs <-  navigate_nldi(nldi_feature = nldi_feature, mode = "UT", distance_km = 1000)$UT$nhdplus_comid
  
  id_list      <- paste0("'", UT_COMIDs, "'", collapse = ",")
  
  query_suffix <- paste(glue("COMID IN ({id_list})"), collapse = " OR ")
  query        <- glue("SELECT * FROM {flowpaths} WHERE {query_suffix}")
  flowlines            <- st_read(gpkg, query = query, quiet = TRUE)
  
  query_suffix <- paste(glue("FEATUREID IN ({id_list})"), collapse = " OR ")
  query        <- glue("SELECT * FROM {catchments} WHERE {query_suffix}")
  catchments           <- st_read(gpkg, query = query, quiet = TRUE)
  
  
  if(!is.null(export_gpkg)){
    write_hydrofabric(list(flowlines = flowlines, catchments = catchments), export_gpkg, enforce_dm = FALSE)
  } else {
    return(list(flowpaths = flowlines, catchments = catchments))
  }
}


#' Subset Hydrofabric Network
#' @param id hydrofabric id
#' @param comid NHDPlusV2 COMID
#' @param hl_id Hydrolocation URI
#' @param network Network Parquet file
#' @param pattern Pattern for distributed VPU based GPKGS
#' @param lyrs layers to extract
#' @param export_gpkg file path to write to
#'
#' @return file path or list 
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
      slice_max(hydroseq) %>% 
      # We want everthing that flows to this ID
      pull(id) %>% 
      unique()
  }
  
  if(!is.null(hl_id)){
    origin = filter(net, hl_uri == !!hl_id) %>% 
      slice_max(hydroseq) %>% 
      # We want everything that flows to this HL, HL are nexus.
      pull(toid) %>% 
      unique()
  }
  
  if(!is.null(id)){ origin = id }
  
  message("Starting from: `",  origin, "`")

  vpuid = unique(pull(filter(net, id == origin | toid == origin), vpu))
  
  sub_net = select(filter(net, vpu == vpuid), -vpu)

  tmap = get_sorted(
    distinct(select(sub_net, id, toid)), 
    outlets = origin
  )
  
  if(!grepl("nex-", tmap$toid[nrow(tmap)])){
    tmap = tmap[-nrow(tmap),]
  }
  
  ids = unique(c(unlist(tmap)))
  
  gpkg = glue(pattern, vpu = vpuid)
  
  db <- dbConnect(SQLite(), gpkg)
  on.exit(dbDisconnect(db))
  
  hydrofabric = list()
  
  for(j in 1:length(lyrs)){
    hyaggregate_log("INFO", glue("Subsetting: {lyrs[j]} ({j}/{length(lyrs)})"))
    
    crs = st_layers(gpkg)$crs
    
    t = tbl(db, lyrs[j]) %>%
      filter(if_any(any_of(c('divide_id', 'id', 'ds_id')), ~ . %in% ids)) %>%
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

