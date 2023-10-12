#' Subset Hydrofabric Network
#' @param id hydrofabric id (relevant only to nextgen fabrics)
#' @param comid NHDPlusV2 COMID 
#' @param hl_id hydrolocation URI (relevant only to nextgen fabrics)
#' @param nldi_feature list with names 'featureSource' and 'featureID' where 'featureSource' is derived from the "source" column of the response of dataRetrieval::get_nldi_sources() and the 'featureID' is a known identifier from the specified 'featureSource'.
#' @param loc Location given as vector of XY in CRS 4326 (long, lat)
#' @param base_s3 the base hydrofabric directory to access in Lynker's s3
#' @param lyrs layers to extract. Default is all possible in the hydrofabric GPKG data model
#' @param outfile file path to write to. Must have ".gpkg" extension
#' @param cache_dir should data be cached to a locat directory? Will speed up multiple subsets in the same region
#' @param cache_overwrite description. Should a cached file be overwritten
#' @return file path (outfile) or list of features
#' @export

subset_network = function(id = NULL,
                          comid = NULL,
                          hl_id = NULL,
                          nldi_feature = NULL,
                          loc = NULL,
                          base_s3 = 's3://nextgen-hydrofabric/pre-release/',
                          lyrs  = c(
                            "divides",
                            "nexus",
                            "flowpaths",
                            "network",
                            "hydrolocations",
                            "reference_flowline",
                            "reference_catchment"
                          ),
                          outfile = NULL,
                          cache_dir = NULL,
                          cache_overwrite = FALSE) {
  
  net = tryCatch({
    open_dataset(glue(base_s3, "conus_net.parquet")) %>%
      select(id, toid, hf_id, hl_uri, hf_hydroseq, hydroseq, vpu) %>%
      collect() %>%
      distinct()
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(id) & !is.null(net)) {
    comid = filter(net, id == !!id | toid == !!id) %>%
      slice_max(hf_hydroseq) %>%
      pull(hf_id)
  }
  
  if (!is.null(nldi_feature)) {
    comid = get_nldi_feature(nldi_feature)$comid
  }
  
  if (!is.null(loc)) {
    comid = discover_nhdplus_id(point = st_sfc(st_point(c(loc[1], loc[2])), crs = 4326))
  }
  
  if (!is.null(hl_id) & !is.null(net)) {
    tmp = filter(net, hl_uri == !!hl_id) %>%
      slice_max(hf_hydroseq) %>%
      pull(toid) %>%
      unique()
    
    origin = filter(net, id == tmp) %>%
      pull(toid) %>%
      unique()
  }
  
  if (!is.null(comid) & !is.null(net)) {
    origin = filter(net, hf_id == comid) %>%
      slice_max(hf_hydroseq) %>%
      pull(id) %>%
      unique()
  } else if (is.null(net)) {
    origin = comid
  }
  
  if (is.null(origin)) {
    stop("origin not found")
  }
  
  message("Starting from: `",  origin, "`")
  
  if (is.null(net)) {
    xx = suppressMessages({
      get_nhdplus(comid = comid)
    })

    vpuid = vpu_boundaries$VPUID[which(lengths(st_intersects(st_transform(vpu_boundaries, st_crs(xx)), xx)) > 0)]
  } else {
    vpuid = unique(pull(filter(net, id == origin |
                                 toid == origin), vpu))
  }
  
  if (!is.null(base_s3)) {
    xx = get_bucket_df(bucket = dirname(base_s3), prefix = basename(base_s3)) %>%
      filter(grepl(basename(base_s3), Key) &
               grepl(paste0(vpuid, ".gpkg$"), Key)) %>%
      filter(!grepl("[.]_", Key)) %>%
      filter(!grepl("/", dirname(Key)))
    
    if (!is.null(cache_dir)) {
      dir.create(cache_dir,
                 recursive = TRUE,
                 showWarnings = FALSE)
      gpkg = glue("{cache_dir}/{basename(xx$Key)}")
      if (cache_overwrite) {
        unlink(gpkg)
      }
      temp = FALSE
    } else {
      gpkg = tempfile(fileext  = ".gpkg")
      temp = TRUE
    }
    
    if (!file.exists(gpkg)) {
      save_object(bucket = xx$Bucket,
                  object = xx$Key,
                  file = gpkg)
    }
    
    lyrs = lyrs[lyrs %in% st_layers(gpkg)$name]
    flowpaths = lyrs[grepl("flowline|flowpath", lyrs)]
    catchments = lyrs[grepl("divides|catchment", lyrs)]
  }
  
  db <- dbConnect(SQLite(), gpkg)
  on.exit(dbDisconnect(db))
  
  if (!is.null(net)) {
    sub_net = distinct(select(filter(net, vpu == vpuid), id, toid))
  } else {
    sub_net = tbl(db, flowpaths) %>%
      select(any_of(c("id", "toid", 'COMID', 'toCOMID'))) %>%
      collect() %>%
      rename(id = COMID, toid = toCOMID)
  }
  
  tmap = suppressWarnings({ get_sorted(distinct(sub_net), outlets = origin) })
  
  ids = unique(c(unlist(tmap)))
  
  hydrofabric = list()
  
  for (j in 1:length(lyrs)) {
    hyaggregate_log("INFO", glue("Subsetting: {lyrs[j]} ({j}/{length(lyrs)})"))
    
    crs = st_layers(gpkg)$crs
    
    t = tbl(db, lyrs[j]) %>%
      filter(if_any(any_of(
        c('COMID',  'FEATUREID', 'divide_id', 'id', 'ds_id')
      ), ~ . %in% ids)) %>%
      collect()
    
    if (all(!any(is.na(as.character(crs[[j]]))), nrow(t) > 0)) {
      if (any(c("geometry", "geom") %in% names(t))) {
        t = st_as_sf(t, crs = crs[[j]])
      } else {
        t = t
      }
    }
    
    if (!is.null(outfile)) {
      write_sf(t, outfile, lyrs[j])
    } else {
      hydrofabric[[lyrs[j]]] = t
    }
  }
  
  if(temp){
    unlink(gpkg)
  }
  
  if (!is.null(outfile)) {
    outfile = append_style(outfile, layer_names = lyrs)
    return(outfile)
  } else {
    hydrofabric
  }
  
}
