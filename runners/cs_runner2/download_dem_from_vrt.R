library(raster)
library(httr)
library(terra)
library(dplyr)

# get the main config variables / paths 
# source("runners/cs_runner/config_vars.R")
source("runners/cs_runner2/base_variables.R")
source("runners/cs_runner2/utils.R")

base_dirs   <- get_base_dir_paths(BASE_DIR)

DEM_VRT_DIR <- base_dirs$dem_vrt_dir
DEM_TIF_DIR <- base_dirs$dem_tif_dir

# "aws s3 sync s3://prd-tnm/StagedProducts/Elevation/1/TIFF/ /Volumes/T7SSD/lynker-spatial/dem/tif/ --no-sign-request --only-show-errors"

# Parse the VRT file
vrt_file <- list.files(DEM_VRT_DIR, full.names = TRUE)

vrt_tiles <- terra::vrt_tiles(vrt_file)

# tile_path <- vrt_tiles[500]
tile_paths <- gsub("/vsicurl/", "",  vrt_tiles)

# TODO: use AWS S3 sync command like below
s3_sync_cmd <- paste0("aws s3 sync s3://prd-tnm/StagedProducts/Elevation/1/TIFF/ ", DEM_TIF_DIR, " --profile ", AWS_PROFILE, " --no-sign-request --only-show-errors")

# "aws s3 sync s3://prd-tnm/StagedProducts/Elevation/1/TIFF/ /Volumes/T7SSD/lynker-spatial/dem/tif/ --no-sign-request --only-show-errors"

# TODO: Old method curl request each file individually.... super slow
# error_tiles <- download_tiles(tile_paths, DEM_TIF_DIR)