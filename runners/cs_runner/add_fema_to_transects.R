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
  
  list.files(FEMA_vpu_dir)
  
  VPU
  
  fema <- sf::read_sf()
  
  
  
    