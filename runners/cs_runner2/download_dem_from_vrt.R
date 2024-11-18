library(terra)
library(dplyr)

# get the main config variables / paths 
# source("runners/cs_runner/config_vars.R")
source("runners/cs_runner2/base_variables.R")
source("runners/cs_runner2/utils.R")

# specify where VRT should go
DEM_VRT_DIR <- BASE_DIRS_LIST$dem_vrt_dir
# # DEM_TIF_DIR <- base_dirs$dem_tif_dir
# DEM_PATH        <- "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13/TIFF/USGS_Seamless_DEM_13.vrt"
# # /vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13/TIFF/current/n70w148/USGS_13_n70w148.tif
# message("Syncing 3DEP elevation VRTs...")
# s3_cp_cmd <- paste0("aws s3 cp ", "s3://prd-tnm/StagedProducts/Elevation/13/TIFF/USGS_Seamless_DEM_13.vrt", " ", "/Users/anguswatters/Desktop/USGS_Seamless_DEM_13.vrt", " --profile ", AWS_PROFILE, " --no-sign-request")
# s3_cp_cmd
# system(s3_cp_cmd)
# "s3://prd-tnm/StagedProducts/Elevation/13/TIFF/USGS_Seamless_DEM_13.vrt"
# 
# VRT_file <- "/Users/anguswatters/Desktop/USGS_Seamless_DEM_13.vrt"
# vrt_tiles <- terra::vrt_tiles(VRT_file)
# 
# dem <- terra::rast(vrt_tiles[999])
# dem
# terra::writeRaster(dem, "/Users/anguswatters/Desktop/tile1_USGS_Seamless_DEM_13.tif", overwrite = T)
# length(vrt_tiles)
# s3://lynker-nypa-inputs/ryan_test/
  

for (i in seq_along(nrow(NED_META))) {
  
  NED_META
  i = 1
  ned_uri         <- NED_META$s3_uri[i]  
  ned_current_uri <- paste0(ned_uri, "current/")
  ned_current_uri
  
  output_dir <- file.path(DEM_VRT_DIR,  NED_META$id[i])
  output_dir
  current_output_dir <- file.path(DEM_VRT_DIR,  NED_META$id[i], "current")
  
  # create specific VRT directory if it doesnt already exist
  create_if_not_exists(output_dir)
  create_if_not_exists(current_output_dir)
  
  message("[", i, "] - ", "Syncing\n '", ned_current_uri, "' > '", current_output_dir, "'")
  
  s3_sync_cmd <- paste0("aws s3 sync ", ned_current_uri, " ", current_output_dir, " --profile ", AWS_PROFILE, " --no-sign-request")
  # s3_sync_cmd <- paste0("aws s3 sync s3://prd-tnm/StagedProducts/Elevation/1/TIFF/ ", DEM_VRT_DIR, " --profile ", AWS_PROFILE, " --no-sign-request")
  
  
  tryCatch({
    
    # system(s3_sync_cmd)
    
  }, error = function(e) {
    
    message("Error syncing ", ned_uri, " to local directory...")
    message(e)
    
  })
  
}

# vrt_file <- file.path(DEM_VRT_DIR, "USGS_Seamless_DEM_1.vrt")
# 
# # TODO: use AWS S3 sync command like below
# s3_sync_cmd <- paste0("aws s3 sync s3://prd-tnm/StagedProducts/Elevation/1/TIFF/ ", DEM_VRT_DIR, " --profile ", AWS_PROFILE, " --no-sign-request")
# # s3_sync_cmd <- paste0("aws s3 sync s3://prd-tnm/StagedProducts/Elevation/1/TIFF/ ", DEM_VRT_DIR, " --profile ", AWS_PROFILE, " --no-sign-request --only-show-errors")
# 
# system(s3_sync_cmd)

# ------------------------------------
# ---- Build VRTs ----
# ------------------------------------

message("Building 3DEP elevation '.vrt' files...")

for (i in seq_along(nrow(NED_META))) {
  
  i = 1
  NED_META$id[i]
  
  output_dir <- file.path(DEM_VRT_DIR,  NED_META$id[i])
  
  # output_dir 
  
  
  # Define output text file path
  txt_file  <- file.path(output_dir, paste0("ned_list_", NED_META$id[i], ".txt"))
  
  # Define output VRT path
  vrt_file  <- file.path(output_dir, paste0("ned_", NED_META$id[i], ".vrt"))
  # vrt_file <- paste0(paste0("data-raw/ned_", NED_META$id[i], ".vrt"))
  
  # If VRT does NOT exist, build VRT
  if (!file.exists(vrt_file)) {
    
    # txt_file
  
    current_dir <- file.path(output_dir, "current")
    
    # all the current/ directories that have the TIFs
    tile_dirs   <- list.files(current_dir, full.names = TRUE)
    # list.files(tile_dirs[1:2],  full.names = T)
    # index_gpkg <- sf::read_sf(list.files(tile_dirs[999], full.names = T, pattern = ".gpkg"))
    
    # get list of all the .tif files in the "current/" dir
    tif_paths   <- list.files(tile_dirs,  pattern = ".tif", full.names = T)
    
    # write list of files to text file
    write.table(tif_paths, txt_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # build VRT from text file input using GDAL system call ...
    system(paste("gdalbuildvrt -input_file_list", txt_file, vrt_file))
    
  }
  
}


# file.path(DEM_TIF_DIR, "USGS_Seamless_DEM_1.vrt")
# length(list.files(file.path(DEM_TIF_DIR, "current"), full.names = T))
# vrt_tiles <- terra::vrt_tiles(file.path(DEM_TIF_DIR, "USGS_Seamless_DEM_1.vrt"))
# vrt_tiles[1]
# # "aws s3 sync s3://prd-tnm/StagedProducts/Elevation/1/TIFF/ /Volumes/T7SSD/lynker-spatial/dem/tif/ --no-sign-request --only-show-errors"
# 
# # Parse the VRT file
# vrt_file <- list.files(DEM_VRT_DIR, full.names = TRUE)
# vrt_file <- list.files(DEM_VRT_DIR, full.names = TRUE)
# file.path(DEM_TIF_DIR, "USGS_Seamless_DEM_1.vrt")
# vrt_tiles <- terra::vrt_tiles(vrt_file)
# 
# # tile_path <- vrt_tiles[500]
# tile_paths <- gsub("/vsicurl/", "",  vrt_tiles)
# 
# # TODO: use AWS S3 sync command like below
# s3_sync_cmd <- paste0("aws s3 sync s3://prd-tnm/StagedProducts/Elevation/1/TIFF/ ", DEM_VRT_DIR, " --profile ", AWS_PROFILE, " --no-sign-request --only-show-errors")
# # aws s3 sync s3://prd-tnm/StagedProducts/Elevation/1/TIFF/ /Volumes/T7SSD/lynker-spatial/dem/tif/ --no-sign-request --only-show-errors
# # "aws s3 sync s3://prd-tnm/StagedProducts/Elevation/1/TIFF/ /Volumes/T7SSD/lynker-spatial/dem/tif/ --no-sign-request --only-show-errors"
# 
# # TODO: Old method curl request each file individually.... super slow
# # error_tiles <- download_tiles(tile_paths, DEM_TIF_DIR)