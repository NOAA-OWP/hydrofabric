library(dplyr) 

# Generate the flowlines layer for the final cross_sections_<VPU>.gpkg for each VPU
source("runners/cs_runner/config.R")

# # # # load libraries
# library(hydrofabric3D)
# library(dplyr)
# library(sf)
# install.packages("devtools")

# transect bucket prefix
transects_prefix <- paste0(s3_bucket, version_prefix, "/3D/transects/")

# paths to nextgen datasets and model attribute parquet files
nextgen_files   <- list.files(nextgen_dir, full.names = FALSE)
FEMA_files      <- list.files(FEMA_FGB_PATH, full.names = FALSE)
transects_files <- list.files(transects_dir, full.names = FALSE)

transects_files <- transects_files[!grepl("updated", transects_files)]

# string to fill in "cs_source" column in output datasets
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
for(i in 1:nrow(path_df)) {
  
  i = 8
  
  # nextgen file and full path
  nextgen_file <- path_df$x[i]
  nextgen_path <- paste0(nextgen_dir, nextgen_file)
  
  transect_file <- path_df$y[i]
  transect_path <- paste0(transects_dir, transect_file)
  
  transect_path
  
  # # model attributes file and full path
  # model_attr_file <- path_df$y[i]
  # model_attr_path <- paste0(model_attr_dir, model_attr_file)
  
  message("Creating VPU ", path_df$vpu[i], "\n - transects: ", transect_file, "\n - flowpaths: '", nextgen_file, "'")
  # message("Creating VPU ", path_df$vpu[i], " transects:\n - flowpaths: '", nextgen_file, "'\n - model attributes: '", model_attr_file, "'")
  
  # read in nextgen data
  flines <- sf::read_sf(nextgen_path, layer = "flowpaths")
  flines_bb <- 
    flines %>% 
    sf::st_bbox() %>% 
    sf::st_as_sfc() %>% 
    sf::st_as_sf()
  
  transects <- sf::read_sf(transect_path)
  
  
  # find the states intersecting with the given VPU flowlines
  intersecting_states <- 
    sf::st_intersection(us_states, flines_bb) %>% 
    sf::st_drop_geometry() %>% 
    .$name %>% 
    gsub(" ", "-", .)
  
  # get the matching FEMA floodplain FGB file names
  matching_fema_files <- unlist(lapply(intersecting_states, function(state_name) {
    FEMA_files[grepl(state_name, FEMA_files)]
    }))
  
  # full paths
  files_of_interest <-  paste0(FEMA_FGB_PATH, "/", matching_fema_files)
  
  # Iterate over each FEMA file and determine optimal widths for cross sections.....
  # for (file in rev(files_of_interest)) {
  
    fema_fgb <- 
      file %>% 
      sf::read_sf() %>% 
      sf::st_transform(5070)
    
    fema_bb <- 
      fema_fgb %>% 
      sf::st_bbox() %>% 
      sf::st_as_sfc() %>% 
      sf::st_as_sf()
    
  
  }