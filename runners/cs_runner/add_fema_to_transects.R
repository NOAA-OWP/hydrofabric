library(dplyr) 

# Generate the flowlines layer for the final cross_sections_<VPU>.gpkg for each VPU
source("runners/cs_runner/config.R")

# transect bucket prefix
transects_prefix <- paste0(s3_bucket, version_prefix, "/3D/transects/")

# paths to nextgen datasets and model attribute parquet files
nextgen_files   <- list.files(nextgen_dir, full.names = FALSE)
# FEMA_files      <- list.files(FEMA_FGB_PATH, full.names = FALSE)
# FEMA_BB_files   <- list.files(FEMA_FGB_BB_PATH, full.names = FALSE)
transects_files <- list.files(transects_dir, full.names = FALSE)
transects_files <- transects_files[!grepl("updated", transects_files)]

FEMA_VPU_SUBFOLDERS
# string to fill in "cs_source" column in output datasets
FEMA_files
net_source <- "hydrofabric3D"

# ensure the files are in the same order and matched up by VPU
path_df <- align_files_by_vpu(
  x    = nextgen_files,
  y    = transects_files,
  base = base_dir
)

path_df

us_states <- 
  USAboundaries::us_states() %>% 
  sf::st_transform(5070)

# loop over each VPU and generate cross sections, then save locally and upload to S3 bucket
# for(i in 1:nrow(path_df)) {
  
  i = 8
  
  # nextgen file and full path
  nextgen_file <- path_df$x[i]
  nextgen_path <- paste0(nextgen_dir, nextgen_file)
  
  transect_file <- path_df$y[i]
  transect_path <- paste0(transects_dir, transect_file)
  
  VPU <- path_df$vpu[i]
  transect_path
  
  # # model attributes file and full path
  # model_attr_file <- path_df$y[i]
  # model_attr_path <- paste0(model_attr_dir, model_attr_file)
  
  message("Creating VPU ", path_df$vpu[i], "\n - transects: ", transect_file, "\n - flowpaths: '", nextgen_file, "'")
  # message("Creating VPU ", path_df$vpu[i], " transects:\n - flowpaths: '", nextgen_file, "'\n - model attributes: '", model_attr_file, "'")
  
  # read in nextgen data
  flines <- sf::read_sf(nextgen_path, layer = "flowpaths")

  
  FEMA_vpu_dir <- FEMA_VPU_SUBFOLDERS[grepl(paste0("VPU_", VPU), basename(FEMA_VPU_SUBFOLDERS))]
  
  vpu_fema_file <- list.files(FEMA_vpu_dir, full.names = TRUE)
  
  FEMA_vpu_dir
  VPU
  
  fema <- sf::read_sf(vpu_fema_file)
    
   
  transects <- sf::read_sf(transect_path)
  
  transects_geos <- geos::as_geos_geometry(transects)
  
  fema_geos <- geos::as_geos_geometry(fema)     

  
  fema_geos  
  fema_transects_matrix <- geos::geos_intersects_matrix(transects_geos, fema_geos) 
  transects_fema_matrix <- geos::geos_intersects_matrix(fema_geos, transects_geos) 
  
  fema_transects_matrix
  # get the polygons that have atleast 1 intersection with the 'lines'
  fema_transects_matrix[[557]]
  fema_tmp <- fema[fema_transects_matrix[[557]], ]
  trans_tmp <- transects[557, ]
  
 fema_dissolve <-  rmapshaper::ms_dissolve(fema_tmp)
 fema_simple <- rmapshaper::ms_simplify(fema_tmp, keep = 1) %>% 
   # rmapshaper::ms_dissolve(field = "state")
 sf::st_union()
  mapview::mapview(trans_tmp, color = "green") + 
    mapview::mapview(fema_tmp[1, ], col.region = "red") + 
    mapview::mapview(fema_tmp[2, ],  col.region = "dodgerblue") + 
    mapview::mapview(fema_simple,  col.region = "yellow")
  transects
  
  
  fema_transects_matrix
  
  transects_with_fema <- transects[lengths(fema_transects_matrix) != 0, ]
  fema_with_transects <- fema[lengths(transects_fema_matrix) != 0, ]
  
  lengths(transects_fema_matrix)
  mapview::mapview(transects_with_fema, color = "green") + fema_with_transects
  unique(hydrofabric3D::add_tmp_id(transects)$tmp_id)[1:30]
  transects  %>% 
    hydrofabric3D::add_tmp_id(transects) %>% 
    .$tmp_id %>% 
    unique() %>% .[1:30]
  trans_subset <- 
    transects %>% 
    hydrofabric3D::add_tmp_id() %>% 
    dplyr::filter(tmp_id %in%  unique(hydrofabric3D::add_tmp_id(transects)$tmp_id)[1:30])
  
  fema_subset <- 
    fema %>% 
    dplyr::filter(fema_id == "1268")

  extended <- hydrofabric3D:::extend_by_length(trans_subset,    rep(500, nrow(trans_subset)))
  extended
  clipped_trans <- rmapshaper::ms_clip(extended, fema)
  
  rmapshaper::ms_clip(extended, fema_subset)
  mapview::mapview(trans_subset, color = "red") + 
    mapview::mapview(extended, color = "yellow") +
    # mapview::mapview(  sf::st_difference(extended, fema_subset), color = "green") +
    mapview::mapview(clipped_trans,  color = "green") +
    fema 
  
   rep(50, nrow(trans_subset))
  extended <- hydrofabric3D:::extend_by_length(trans_subset,    rep(50, nrow(trans_subset)))
  
  
  hydrofabric3D::geos_extend_line(trans_subset, 50) %>% 
    sf::st_as_sf() %>% mapview::mapview()
  
  
  
  
  
  
  
  
  
  
  
  
  
 
  
  
       
  
     
  
     
  
     
  
    
  
    