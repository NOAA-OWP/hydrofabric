### EDIT base_dir, aws_profile, and DEM_URL ###

# ----------------------------------------------------------------------------
# ---- General paths and constants variables ----
# ----------------------------------------------------------------------------

base_dir     <- '/Users/anguswatters/Desktop/lynker-spatial'

# AWS profile to run CLI commands 
aws_profile  <- "angus-lynker"

# name of S3 bucket
s3_bucket    <- "s3://lynker-spatial/"

# S3 prefix/folder of version run
version_prefix <- "v20.1"

# location of FEMA 100 year flood plain FGB files
FEMA_S3_BUCKET <- "s3://lynker-hydrofabric/"
FEMA_S3_BUCKET_PREFIX <- "FEMA100/"
FEMA_S3_DIR <- paste0(FEMA_S3_BUCKET, FEMA_S3_BUCKET_PREFIX)

# FEMA100 year flood map FGB save location (temporary, will be deleted after processing)
FEMA_FGB_PATH <- paste0(base_dir, "/FEMA100")
FEMA_GEOJSON_PATH    <- paste0(base_dir, "/FEMA100_geojson")
FEMA_CLEAN_PATH      <- paste0(base_dir, "/FEMA100_clean")
FEMA_GPKG_PATH       <- paste0(base_dir, "/FEMA100_gpkg")
FEMA_GPKG_BB_PATH    <- paste0(base_dir, "/FEMA100_bounding_box") # TODO: Probably can be deleted too, not sure yet

# TODO: these can be deleted
# FEMA_SIMPLIFIED_PATH <- paste0(base_dir, "/FEMA100_simplified")
# FEMA_DISSOLVED_PATH  <- paste0(base_dir, "/FEMA100_dissolved")
# FEMA_EXPLODED_PATH   <- paste0(base_dir, "/FEMA100_exploded")


# ----------------------------------------------------------------------------
# ---- Cross section point extraction constant variables ----
# ----------------------------------------------------------------------------

# DEM URL
DEM_URL      <- "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/TIFF/USGS_Seamless_DEM_1.vrt"

# scale argument for cross_section_pts() function. 
# The percentage of the length of the transect line to try and extend a transect to see if viable Z values can be found by extending transect line
# Default setting is 50% of the original transect lines length (0.5)
EXTENSION_PCT <- 0.5

# percentage of the length each cross section that should be used as a threshold for classifying a cross section as having relief or not
# 1% of the cross sections length is the default value we are using 
# (i.e. a 100m long cross section needs a minimum of 1 meter (1%) of relief in its cross section points to be classified as "having relief")
PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF <- 0.01

# Whether to collect meta data from runs to generate an output CSV (currently only being created in 02_cs_pts.R)
# TODO: Probably delete this 
COLLECT_META <- TRUE

# Where should meta data CSVs be saved to? 
# Local path to save CSVs of cross section meta data during each iteration
# TODO: Probably delete this 
META_PATH <- '/Users/anguswatters/Desktop/cs_meta/'
# META_PATH <-  "/local/path/to/save/cross_section_meta_data/"

# ----------------------------------------------------------------------------
# ---- Machine learning data path variables ----
# ----------------------------------------------------------------------------

ML_OUTPUTS_FILE   = "channel_ml_outputs.parquet"
ML_OUTPUTS_PREFIX = "v20.1/3D/ml-outputs/"
ML_OUTPUTS_URI    = paste0(s3_bucket, ML_OUTPUTS_PREFIX, ML_OUTPUTS_FILE)
# ML_OUTPUTS_URI    = "s3://lynker-spatial/v20.1/3D/ml-outputs/channel_ml_outputs.parquet"

ML_OUTPUTS_PATH <- paste0(base_dir, "/ml-outputs/", ML_OUTPUTS_FILE)

# path to the remote CONUS net parquet file
CONUS_NETWORK_FILENAME <- "conus_net.parquet"
CONUS_NETWORK_URI      <- paste0(s3_bucket, version_prefix, "/", CONUS_NETWORK_FILENAME)

### EDIT ###
