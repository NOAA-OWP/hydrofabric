### EDIT base_dir, aws_profile, and DEM_URL ###
base_dir     <- '/Users/anguswatters/Desktop/lynker-spatial'

# AWS profile to run CLI commands 
aws_profile  <- "angus-lynker"

# name of S3 bucket
s3_bucket    <- "s3://lynker-spatial/"

# location of FEMA 100 year flood plain FGB files
FEMA_S3_BUCKET <- "s3://lynker-hydrofabric/"
FEMA_S3_BUCKET_PREFIX <- "FEMA100/"
FEMA_S3_DIR <- paste0(FEMA_S3_BUCKET, FEMA_S3_BUCKET_PREFIX)

# FEMA100 year flood map FGB save location (temporary, will be deleted after processing)
FEMA_FGB_PATH <- paste0(base_dir, "/FEMA100")

# DEM URL
DEM_URL      <- "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/TIFF/USGS_Seamless_DEM_1.vrt"

# scale argument for cross_section_pts() function. 
# The percentage of the length of the transect line to try and extend a transect to see if viable Z values can be found by extending transect line
# Default setting is 50% of the original transect lines length (0.5)
EXTENSION_PCT <- 0.5

# percentage of the length each cross section that should be used as a threshold for classifying a cross section as having relief or not
# 1% of the cross sections length is the default value we are using 
# (i.e. a 100m long cross section needs a minimum of 1 meter (1%) of relief in its cross section points to be classified as "having relief")
PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF = 0.01

# Whether to collect meta data from runs to generate an output CSV (currently only being created in 02_cs_pts.R)
COLLECT_META <- TRUE

# Where should meta data CSVs be saved to? 
# Local path to save CSVs of cross section meta data during each iteration
META_PATH <- '/Users/anguswatters/Desktop/cs_meta/'
# META_PATH <-  "/local/path/to/save/cross_section_meta_data/"



# # create the directory if it does NOT exist
# if(!dir.exists(base_dir)) {
#   message(glue::glue('Base directory does not exist...\nCreating directory: {base_dir}'))
#   dir.create(base_dir)
# }

### EDIT ###
