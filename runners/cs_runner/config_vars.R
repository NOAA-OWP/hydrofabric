### EDIT base_dir, aws_profile, and DEM_URL ###
base_dir     <- '/Users/anguswatters/Desktop/lynker-spatial'

# AWS profile to run CLI commands 
aws_profile  <- "angus-lynker"

# name of S3 bucket
s3_bucket    <- "s3://lynker-spatial/"

# DEM URL
DEM_URL      <- "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/TIFF/USGS_Seamless_DEM_1.vrt"

# scale argument for cross_section_pts() function. 
# The percentage of the length of the transect line to try and extend a transect to see if viable Z values can be found by extending transect line
# Default setting is 50% of the original transect lines length (0.5)
EXTENSION_PCT <- 0.5

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
