### EDIT BASE_DIR, AWS_PROFILE, and DEM_PATH ###
# load utils 
source("runners/cs_runner2/utils.R")

# ---------------------------------------------------------------------------------
# ---- General paths and constants variables ----
# - edit to match your local environment
# - BASE_DIR: base directory for local file storage
# - AWS_PROFILE: AWS profile to run CLI commands
# - VERSION: S3 prefix/folder of version to run / generate hydrofabric data for
# ---------------------------------------------------------------------------------
# Base directory for local file storage
BASE_DIR    <- '/Volumes/T7SSD/lynker-spatial'
# BASE_DIR           <- '/Users/anguswatters/Desktop/lynker-spatial'

BASE_DIRS_LIST   <- get_base_dir_paths(BASE_DIR)

# AWS profile to run CLI commands 
AWS_PROFILE        <- "angus-lynker"

# S3 prefix/folder of version run
VERSION            <- "v3.0"
CONUS_VERSION      <- "v2.2"

VERSION_DIRS_LIST            <- get_version_base_dir_paths(BASE_DIR, VERSION)

# string to fill in "CS_SOURCE" column in output datasets
CS_SOURCE          <- "hydrofabric3D"

# -------------------------------------------------------------------------------------
# ---- S3 BUCKET NAMES ----
# -------------------------------------------------------------------------------------

# name of bucket with nextgen data
LYNKER_SPATIAL_S3_BUCKET_NAME   <- "lynker-spatial"
LYNKER_SPATIAL_HF_S3_PREFIX     <- "hydrofabric"

# AWS S3 bucket URI 
LYNKER_SPATIAL_BASE_S3_URI  <- paste0("s3://", LYNKER_SPATIAL_S3_BUCKET_NAME, "/")
LYNKER_SPATIAL_HF_S3_URI    <- paste0(LYNKER_SPATIAL_BASE_S3_URI, LYNKER_SPATIAL_HF_S3_PREFIX, "/")

# name of bucket with nextgen data
LYNKER_HF_S3_BUCKET_NAME        <- "lynker-hydrofabric"

# -------------------------------------------------------------------------------------
# ---- VPU IDs  ----
# -------------------------------------------------------------------------------------
VPU_IDS              <- c('01', '02', '03N', '03S', '03W', 
                          '04', '05', '06', '07', '08', '09', 
                          '10L', '10U', '11', '12', '13', '14', 
                          '15', '16', '17', '18', '20', '21')
# VPU_IDS              <- sf::st_drop_geometry(nhdplusTools::get_boundaries())$VPUID

# -------------------------------------------------------------------------------------
# ---- CONUS NEXTGEN  ----
# -------------------------------------------------------------------------------------

CONUS_NEXTGEN_S3_URI    <- paste0("s3://", LYNKER_SPATIAL_S3_BUCKET_NAME, "/", LYNKER_SPATIAL_HF_S3_PREFIX, "/", CONUS_VERSION, "/conus/conus_nextgen.gpkg")
CONUS_NEXTGEN_GPKG_PATH <- file.path(VERSION_DIRS_LIST$network_dir, "conus_nextgen.gpkg")

# -------------------------------------------------------------------------------------
# ---- CONUS REFEREMCE FEATURES  ----
# -------------------------------------------------------------------------------------

CONUS_REF_FEATURES_S3_URI    <- paste0("s3://", LYNKER_SPATIAL_S3_BUCKET_NAME, "/", LYNKER_SPATIAL_HF_S3_PREFIX, "/", CONUS_VERSION, "/conus/conus_reference.gpkg")
CONUS_REF_FEATURES_GPKG_PATH <- file.path(VERSION_DIRS_LIST$ref_features_dir, "conus_reference.gpkg")

# -------------------------------------------------------------------------------------
# ---- CONUS ML DATA  ----
# -------------------------------------------------------------------------------------

CONUS_ML_S3_URI          <- paste0("s3://", LYNKER_SPATIAL_S3_BUCKET_NAME, "/", LYNKER_SPATIAL_HF_S3_PREFIX, "/", CONUS_VERSION, "/conus/bathymetry/ml_auxiliary_data.parquet")
CONUS_ML_PARQUET_PATH    <- file.path(VERSION_DIRS_LIST$ml_dir, "ml_auxiliary_data.parquet")

# -------------------------------------------------------------------------------------
# ----  FEMA 100 year flood plain data  ----
# -------------------------------------------------------------------------------------

LYNKER_HF_FEMA_S3_PREFIX        <- "FEMA100"
LYNKER_HF_FEMA_S3_URI           <- paste0("s3://", LYNKER_HF_S3_BUCKET_NAME, "/", LYNKER_HF_FEMA_S3_PREFIX, "/")

FEMA_FGB_PATH       = BASE_DIRS_LIST$fema_fgb_dir
FEMA_GEOJSON_PATH   = BASE_DIRS_LIST$fema_geojson_dir
FEMA_CLEAN_PATH     = BASE_DIRS_LIST$fema_clean_dir
FEMA_GPKG_PATH      = BASE_DIRS_LIST$fema_gpkg_dir
FEMA_BY_VPU_PATH    = BASE_DIRS_LIST$fema_by_vpu_dir
FEMA_VPU_SUBFOLDERS = BASE_DIRS_LIST$fema_by_vpu_subdirs

# -------------------------------------------------------------------------------------
# ---- CONUS CS extension polygons ----
# Data derived from FEMA 100 year flood plain data, puts all FEMA VPU datasets into a 
# single geopackage for CONUS
# -------------------------------------------------------------------------------------

CS_EXTENSION_POLYGONS_DIR  <- BASE_DIRS_LIST$cs_extension_polygons_dir
CONUS_FEMA_GPKG_PATH       <- file.path(CS_EXTENSION_POLYGONS_DIR, 'conus_fema.gpkg')

# -------------------------------------------------------------------------------------
# ---- CONUS CS extension polygons ----
# -------------------------------------------------------------------------------------

DEM_DIR        <- BASE_DIRS_LIST$dem_dir
DEM_VRT_DIR    <- BASE_DIRS_LIST$dem_vrt_dir
DEM_TIF_DIR    <- BASE_DIRS_LIST$dem_tif_dir

# DEM URL
DEM_PATH        <- "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/TIFF/USGS_Seamless_DEM_1.vrt"

# ----------------------------------------------------------------------------
# ---- Cross section point extraction constant variables ----
# ----------------------------------------------------------------------------

# scale argument for cross_section_pts() function. 
# The percentage of the length of the transect line to try and extend a transect to see if viable Z values can be found by extending transect line
# Default setting is 50% of the original transect lines length (0.5)
EXTENSION_PCT  <- 0.5

# percentage of the length each cross section that should be used as a threshold for classifying a cross section as having relief or not
# 1% of the cross sections length is the default value we are using 
# (i.e. a 100m long cross section needs a minimum of 1 meter (1%) of relief in its cross section points to be classified as "having relief")
PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF <- 0.01

# # -------------------------------------------------------------------------------------
# # ---- Create local directory / path variables (FEMA data) ----
# # -------------------------------------------------------------------------------------

# LYNKER_HF_FEMA_S3_PREFIX        <- "FEMA100"
# LYNKER_HF_FEMA_S3_URI           <- paste0("s3://", LYNKER_HF_S3_BUCKET_NAME, "/", LYNKER_HF_FEMA_S3_PREFIX, "/")

# # location of FEMA 100 year flood plain FGB files
# FEMA_S3_BUCKET         <- "s3://lynker-hydrofabric/"
# LYNKER_HF_FEMA_S3_PREFIX  <- "FEMA100/"
# LYNKER_HF_FEMA_S3_URI            <- paste0(FEMA_S3_BUCKET, LYNKER_HF_FEMA_S3_PREFIX)

# # FEMA100 year flood map FGB save location (temporary, will be deleted after processing)
# FEMA_FGB_PATH        <- file.path(BASE_DIRS_LIST$fema_dir, "fema_fgb")
# FEMA_GEOJSON_PATH    <- file.path(BASE_DIRS_LIST$fema_dir, "fema_geojson")
# FEMA_CLEAN_PATH      <- file.path(BASE_DIRS_LIST$fema_dir, "fema_clean")
# FEMA_GPKG_PATH       <- file.path(BASE_DIRS_LIST$fema_dir, "fema_gpkg")
# FEMA_BY_VPU_PATH     <- file.path(BASE_DIRS_LIST$fema_dir, "FEMA_BY_VPU")

# VPU_IDS              <- sf::st_drop_geometry(nhdplusTools::get_boundaries())$VPUID

# FEMA_VPU_SUBFOLDERS  <- paste0(FEMA_BY_VPU_PATH, "/VPU_", VPU_IDS)

# CS_EXTENSION_POLYGONS_DIR  <- paste0(BASE_DIR, "/cs-extension-polygons")
# CONUS_FEMA_GPKG_PATH       <- file.path(CS_EXTENSION_POLYGONS_DIR, 'conus_fema.gpkg')

# DEM_DIR      <- BASE_DIRS_LIST$dem_dir
# DEM_VRT_DIR  <- BASE_DIRS_LIST$dem_vrt_dir
# DEM_TIF_DIR  <- BASE_DIRS_LIST$dem_tif_dir

# # ----------------------------------------------------------------------------
# # ---- Cross section point extraction constant variables ----
# # ----------------------------------------------------------------------------

# # DEM URL
# DEM_PATH        <- "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/TIFF/USGS_Seamless_DEM_1.vrt"

# # scale argument for cross_section_pts() function. 
# # The percentage of the length of the transect line to try and extend a transect to see if viable Z values can be found by extending transect line
# # Default setting is 50% of the original transect lines length (0.5)
# EXTENSION_PCT  <- 0.5

# # percentage of the length each cross section that should be used as a threshold for classifying a cross section as having relief or not
# # 1% of the cross sections length is the default value we are using 
# # (i.e. a 100m long cross section needs a minimum of 1 meter (1%) of relief in its cross section points to be classified as "having relief")
# PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF <- 0.01


# # -------------------------------------------------------------------------------------
# # ---- S3 output directories -----
# # - transects 
# # - cross section points
# # - ML cross section points
# # -------------------------------------------------------------------------------------
# 
# # transect bucket prefix
# S3_TRANSECTS_DIR   <- paste0(LYNKER_SPATIAL_HF_S3_URI, VERSION, "/3D/transects/")
# 
# # cross section bucket prefix
# S3_CS_PTS_DIR      <- paste0(LYNKER_SPATIAL_HF_S3_URI, VERSION, "/3D/dem-cross-sections/")
# 
# # cross section bucket prefix
# S3_CS_ML_PTS_DIR   <- paste0(LYNKER_SPATIAL_HF_S3_URI, VERSION, "/3D/cross-sections/")
# 
# # -------------------------------------------------------------------------------------
# # ---- S3 nextgen data paths / directories -----
# # -------------------------------------------------------------------------------------
# 
# # the name of the folder in the S3 bucket with the nextgen data
# S3_BUCKET_NEXTGEN_DIR       <- paste0(VERSION, "/gpkg/")
# # S3_BUCKET_NEXTGEN_DIR <- "v20.1/gpkg/"
# 
# # full URI to the S3 bucket folder with the nextgen data 
# S3_BUCKET_NEXTGEN_DIR_URI   <- paste0(LYNKER_SPATIAL_HF_S3_URI, S3_BUCKET_NEXTGEN_DIR)
# 
# # reference features S3 bucket prefix
# S3_BUCKET_REF_FEATURES_URI  <- paste0("s3://", LYNKER_SPATIAL_S3_BUCKET_NAME, "/00_reference_features/gpkg/")
# # S3_BUCKET_REF_FEATURES_URI  <- "s3://lynker-spatial/00_reference_features/gpkg/"
# 
# # ----------------------------------------------------------------------------
# # ---- Machine learning data path variables ----
# # ----------------------------------------------------------------------------
# 
# ML_OUTPUTS_S3_FILE   <- "channel_ml_outputs.parquet"
# 
# # ML_OUTPUTS_S3_DIR    <- paste0(VERSION, "/3D/ml-outputs/")
# # ML_OUTPUTS_S3_DIR  <- "v20.1/3D/ml-outputs/"
# 
# ML_OUTPUTS_S3_URI      <- paste0(LYNKER_SPATIAL_HF_S3_URI, VERSION, "/3D/ml-outputs/", ML_OUTPUTS_S3_FILE)
# # ML_OUTPUTS_S3_URI    <- paste0(LYNKER_SPATIAL_HF_S3_URI, ML_OUTPUTS_S3_DIR, ML_OUTPUTS_S3_FILE)
# 
# ML_OUTPUTS_PATH      <- paste0(BASE_DIR, "/ml-outputs/", ML_OUTPUTS_S3_FILE)
# 
# # path to the remote CONUS net parquet file
# CONUS_NETWORK_FILE   <- "conus_net.parquet"
# CONUS_NETWORK_URI    <- paste0(LYNKER_SPATIAL_HF_S3_URI, VERSION, "/", CONUS_NETWORK_FILE)
# 
# # ----------------------------------------------------------------------------
# 
# 
# # -------------------------------------------------------------------------------------
# # ---- Local directory / path variables ----
# # -------------------------------------------------------------------------------------
# 
# # directory to copy nextgen bucket data too
# NEXTGEN_DIR      <- paste0(BASE_DIR, "/", S3_BUCKET_NEXTGEN_DIR)
# # NEXTGEN_DIR <- paste0(BASE_DIR, "/pre-release/")
# 
# # # model attributes directory
# # MODEL_ATTR_DIR   <- paste0(BASE_DIR, "/model_attributes/")
# 
# # cross-section data model data directories
# TRANSECTS_DIR    <- paste0(BASE_DIR, "/01_transects/")
# CS_PTS_DIR       <- paste0(BASE_DIR, "/02_cs_pts/")
# 
# # final output directory with geopackages per VPU
# CS_OUTPUT_DIR    <- paste0(BASE_DIR, "/cross_sections/")
# 
# # directory to copy nextgen bucket data too
# REF_FEATURES_DIR      <- paste0(BASE_DIR, "/00_reference_features/")
# REF_FEATURES_GPKG_DIR <- paste0(REF_FEATURES_DIR, "gpkg/")
# 
# # make a directory for the ML outputs data
# ML_OUTPUTS_DIR   <- paste0(BASE_DIR, "/ml-outputs/")
# 
# DEM_DIR      <- paste0(BASE_DIR, "/dem")
# DEM_VRT_DIR  <- paste0(DEM_DIR, "/vrt")
# DEM_TIF_DIR  <- paste0(DEM_DIR, "/tif")
# 
# # -------------------------------------------------------------------------------------
# # ---- Create local directory / path variables (FEMA data) ----
# # -------------------------------------------------------------------------------------
# 
# # location of FEMA 100 year flood plain FGB files
# FEMA_S3_BUCKET         <- "s3://lynker-hydrofabric/"
# LYNKER_HF_FEMA_S3_PREFIX  <- "FEMA100/"
# LYNKER_HF_FEMA_S3_URI            <- paste0(FEMA_S3_BUCKET, LYNKER_HF_FEMA_S3_PREFIX)
# 
# # FEMA100 year flood map FGB save location (temporary, will be deleted after processing)
# FEMA_FGB_PATH        <- paste0(BASE_DIR, "/FEMA100")
# FEMA_GEOJSON_PATH    <- paste0(BASE_DIR, "/FEMA100_geojson")
# FEMA_CLEAN_PATH      <- paste0(BASE_DIR, "/FEMA100_clean")
# FEMA_GPKG_PATH       <- paste0(BASE_DIR, "/FEMA100_gpkg")
# FEMA_GPKG_BB_PATH    <- paste0(BASE_DIR, "/FEMA100_bounding_box") # TODO: Probably can be deleted too, not sure yet
# 
# FEMA_BY_VPU_PATH     <- paste0(BASE_DIR, "/FEMA_BY_VPU")
# VPU_IDS              <- sf::st_drop_geometry(nhdplusTools::get_boundaries())$VPUID
# 
# FEMA_VPU_SUBFOLDERS  <- paste0(FEMA_BY_VPU_PATH, "/VPU_", VPU_IDS)
# # FEMA_VPU_SUBFOLDERS <- paste0(
# #                         FEMA_BY_VPU_PATH, "/VPU_",
# #                         unlist(
# #                           lapply(list.files(NEXTGEN_DIR, full.names = FALSE), function(vpu_file_names) {
# #                             unlist(regmatches(vpu_file_names,  gregexpr("\\d+[A-Za-z]*", vpu_file_names)))})
# #                           )
# #                         )
# 
# # -------------------------------------------------------------------------------------
# # ---- OVERWRITE_FEMA_FILES constant logicals----
# # ---- > if TRUE, processing steps will be run again 
# #          and overwrite existing previously processed files
# # TODO: Describe these variables
# # -------------------------------------------------------------------------------------
# 
# # Default is TRUE (i.e. a fresh processing run is done from start to finish)
# OVERWRITE_FEMA_FILES  <- TRUE
# DELETE_STAGING_GPKGS  <- TRUE # remove intermediary files from the main output folder
# 

# # -------------------------------------------------------------------------------------
# # ---- (New single domain) Local directory / path variables ----
# # -------------------------------------------------------------------------------------
# 
# # directory for new domain data 
# NEW_DOMAIN_DIRNAME  <- "new_domain"
# NEW_DOMAIN_DIR      <- paste0(BASE_DIR, "/", NEW_DOMAIN_DIRNAME)
# 
# NEW_DOMAIN_FLOWLINES_DIRNAME <- "flowlines"
# NEW_DOMAIN_FLOWLINES_DIR     <- paste0(NEW_DOMAIN_DIR, "/", NEW_DOMAIN_FLOWLINES_DIRNAME)
# 
# NEW_DOMAIN_DEM_DIRNAME <- "dem"
# NEW_DOMAIN_DEM_DIR     <- paste0(NEW_DOMAIN_DIR, "/", NEW_DOMAIN_DEM_DIRNAME)
# 
# NEW_DOMAIN_FLOWLINES_FILE  <- "AllDiffusiveCombined.gpkg"
# NEW_DOMAIN_FLOWLINES_PATH  <- paste0(NEW_DOMAIN_FLOWLINES_DIR, "/", NEW_DOMAIN_FLOWLINES_FILE)
# 
# # # Local DEM file
# # NEW_DOMAIN_DEM_FILE        <- "hi_dem.tif"
# # NEW_DOMAIN_DEM_PATH        <- paste0(NEW_DOMAIN_DEM_DIR, "/", NEW_DOMAIN_DEM_FILE)
# 
# # Remote DEM file
# NEW_DOMAIN_DEM_PATH        <- "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/TIFF/USGS_Seamless_DEM_1.vrt"
# 
# NEW_DOMAIN_TRANSECTS_DIRNAME         <- "transects"
# NEW_DOMAIN_CS_PTS_DIRNAME            <- "cs_pts"
# NEW_DOMAIN_CROSS_SECTIONS_DIRNAME    <- "cross_sections"
# 
# NEW_DOMAIN_TRANSECTS_DIR          <- paste0(NEW_DOMAIN_DIR, "/", NEW_DOMAIN_TRANSECTS_DIRNAME)
# NEW_DOMAIN_CS_PTS_DIR             <- paste0(NEW_DOMAIN_DIR, "/", NEW_DOMAIN_CS_PTS_DIRNAME)
# NEW_DOMAIN_CROSS_SECTIONS_DIR     <- paste0(NEW_DOMAIN_DIR, "/", NEW_DOMAIN_CROSS_SECTIONS_DIRNAME)
# 
# # -------------------------------------------------------------------------------------
# # ---- (New single domain) Local directory / path variables ----
# # -------------------------------------------------------------------------------------
# 
# # directory for new domain data 
# DOMAIN_WITH_FEMA_DIRNAME  <- "domain_with_fema"
# DOMAIN_WITH_FEMA_DIR      <- paste0(BASE_DIR, "/", DOMAIN_WITH_FEMA_DIRNAME)
# 
# DOMAIN_WITH_FEMA_FLOWLINES_DIRNAME <- "flowlines"
# DOMAIN_WITH_FEMA_FLOWLINES_DIR     <- paste0(DOMAIN_WITH_FEMA_DIR, "/", DOMAIN_WITH_FEMA_FLOWLINES_DIRNAME)
# 
# DOMAIN_WITH_FEMA_SUBSET_DIRNAME  <- "domain_subset"
# DOMAIN_WITH_FEMA_SUBSET_DIR      <- paste0(DOMAIN_WITH_FEMA_DIR, "/", DOMAIN_WITH_FEMA_SUBSET_DIRNAME)
# 
# DOMAIN_WITH_FEMA_DEM_DIRNAME     <- "dem"
# DOMAIN_WITH_FEMA_DEM_DIR         <- paste0(DOMAIN_WITH_FEMA_DIR, "/", DOMAIN_WITH_FEMA_DEM_DIRNAME)
# 
# DOMAIN_WITH_FEMA_FLOWLINES_FILE  <- "ls_conus.gpkg"
# DOMAIN_WITH_FEMA_FLOWLINES_PATH  <- paste0(DOMAIN_WITH_FEMA_FLOWLINES_DIR, "/", DOMAIN_WITH_FEMA_FLOWLINES_FILE)
# 
# # Geopackage containing area to subset flowlines to before processing
# DOMAIN_WITH_FEMA_SUBSET_FILE     <- "AllDiffusiveCombined.gpkg"
# DOMAIN_WITH_FEMA_SUBSET_PATH     <- paste0(DOMAIN_WITH_FEMA_SUBSET_DIR, "/", DOMAIN_WITH_FEMA_SUBSET_FILE)
# 
# # # Local DEM file
# # DOMAIN_WITH_FEMA_DEM_FILE        <- "hi_dem.tif"
# # DOMAIN_WITH_FEMA_DEM_PATH        <- paste0(DOMAIN_WITH_FEMA_DEM_DIR, "/", DOMAIN_WITH_FEMA_DEM_FILE)
# 
# # Remote DEM file
# DOMAIN_WITH_FEMA_DEM_PATH        <- "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/TIFF/USGS_Seamless_DEM_1.vrt"
# 
# DOMAIN_WITH_FEMA_TRANSECTS_DIRNAME         <- "transects"
# DOMAIN_WITH_FEMA_CS_PTS_DIRNAME            <- "cs_pts"
# DOMAIN_WITH_FEMA_CROSS_SECTIONS_DIRNAME    <- "cross_sections"
# DOMAIN_WITH_FEMA_OUTPUT_DIRNAME            <- "outputs"
# DOMAIN_WITH_FEMA_VPU_SUBSETS_DIRNAME       <- "vpu-subsets"
# DOMAIN_WITH_FEMA_ML_DIRNAME                <- "ml"
# # aws s3 cp s3://prd-tnm/StagedProducts/Elevation/1/TIFF/USGS_Seamless_DEM_1.vrt /Users/anguswatters/Desktop/3DEP/3DEP.vrt
# DOMAIN_WITH_FEMA_TRANSECTS_DIR          <- paste0(DOMAIN_WITH_FEMA_DIR, "/", DOMAIN_WITH_FEMA_TRANSECTS_DIRNAME)
# DOMAIN_WITH_FEMA_CS_PTS_DIR             <- paste0(DOMAIN_WITH_FEMA_DIR, "/", DOMAIN_WITH_FEMA_CS_PTS_DIRNAME)
# DOMAIN_WITH_FEMA_CROSS_SECTIONS_DIR     <- paste0(DOMAIN_WITH_FEMA_DIR, "/", DOMAIN_WITH_FEMA_CROSS_SECTIONS_DIRNAME)
# DOMAIN_WITH_FEMA_OUTPUT_DIR             <- paste0(DOMAIN_WITH_FEMA_DIR, "/", DOMAIN_WITH_FEMA_OUTPUT_DIRNAME)
# DOMAIN_WITH_FEMA_VPU_SUBSETS_DIR        <- paste0(DOMAIN_WITH_FEMA_DIR, "/", DOMAIN_WITH_FEMA_VPU_SUBSETS_DIRNAME)
# DOMAIN_WITH_FEMA_ML_DIR                 <- paste0(DOMAIN_WITH_FEMA_DIR, "/", DOMAIN_WITH_FEMA_ML_DIRNAME)
# 
# ML_AUXILIARY_DATA_S3_URI <- paste0(LYNKER_SPATIAL_BASE_S3_URI, "bathymetry/ml_auxiliary_data")
# # LYNKER_SPATIAL_BASE_S3_URI
# LYNKER_HYDROFABRIC_S3_BUCKET_BASE_URI   <- "s3://lynker-hydrofabric/"
# ML_BATHYMETRY_S3_DATA_DIR               <- "hydrofabric/nextgen/bathymetry/multisource_river_attributes/"
# ML_BATHYMETRY_S3_URI                    <- paste0(LYNKER_HYDROFABRIC_S3_BUCKET_BASE_URI, ML_BATHYMETRY_S3_DATA_DIR)
# VPU_ML_BATHYMETRY_S3_DIRS               <- paste0(ML_BATHYMETRY_S3_URI, "vpuid=", nhdplusTools::vpu_boundaries$VPUID, "/")
# 
# COASTAL_BATHY_DEM_S3_DIR     <- "coastal_bathy/diffusive_domain/"
# COASTAL_BATHY_DEM_S3_DIR_URI <- paste0(LYNKER_HYDROFABRIC_S3_BUCKET_BASE_URI, COASTAL_BATHY_DEM_S3_DIR)
# # COASTAL_BATHY_DEM_S3_URIS    <- paste0(COASTAL_BATHY_DEM_S3_DIR_URI, list_s3_objects(COASTAL_BATHY_DEM_S3_DIR_URI, ".tif$", AWS_PROFILE))
# 
# 
# # # # --------------------------------------------------------------------------
# # # # ---- Get locations of diffusive domain DEM files in S3 ----
# # # # --------------------------------------------------------------------------
# # COASTAL_BATHY_DEM_S3_URIS    <- paste0(COASTAL_BATHY_DEM_S3_DIR_URI, list_s3_objects(COASTAL_BATHY_DEM_S3_DIR_URI, ".tif$", AWS_PROFILE))
# # 
# # # COASTAL_BATHY_DEM_S3_URIS    <- paste0(COASTAL_BATHY_DEM_S3_DIR_URI, 
# # #                                        list_s3_objects(COASTAL_BATHY_DEM_S3_DIR_URI, ".tif$", AWS_PROFILE)
# # #                                        )
# 
# 






