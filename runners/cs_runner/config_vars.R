### EDIT base_dir, aws_profile, and DEM_URL ###
base_dir    <- '/Users/anguswatters/Desktop/lynker-spatial'

# AWS profile to run CLI commands 
aws_profile <- "angus-lynker"

# DEM URL
DEM_URL     <- "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/TIFF/USGS_Seamless_DEM_1.vrt"

# # create the directory if it does NOT exist
# if(!dir.exists(base_dir)) {
#   message(glue::glue('Base directory does not exist...\nCreating directory: {base_dir}'))
#   dir.create(base_dir)
# }

### EDIT ###
