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

# # create the directory if it does NOT exist
# if(!dir.exists(base_dir)) {
#   message(glue::glue('Base directory does not exist...\nCreating directory: {base_dir}'))
#   dir.create(base_dir)
# }

### EDIT ###
