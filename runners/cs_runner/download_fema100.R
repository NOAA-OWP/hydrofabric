# Running this script goes and pulls the desired FEMA100 flood fgb datasets from the lynker-hydrofabric S3 bucket then saves them into a directory within "base_dir"
# base_dir is defined within runners/workflow/root_dir.R

# NOTE: The lynker-hydrofabric S3 bucket is private at the moment

# load config variables
source("runners/cs_runner/config_vars.R")

# create FEMA100/ directory if it does NOT exist
if (!dir.exists(FEMA_FGB_PATH)) {
  message(glue::glue('FEMA100/ directory does not exist...\nCreating directory: {FEMA_FGB_PATH}'))
  dir.create(FEMA_FGB_PATH)
}


# list objects in S3 bucket, and regular expression match to nextgen_.gpkg pattern
fema_list_command <- paste0('#!/bin/bash
            # AWS S3 Bucket and Directory information
            S3_BUCKET="', FEMA_S3_DIR, '" 
            
            # Regular expression pattern to match object keys
            PATTERN=".fgb$"
            
            # AWS CLI command to list objects in the S3 bucket and use grep to filter them
            S3_OBJECTS=$(aws s3 ls "$S3_BUCKET" --profile ', aws_profile, ' | awk \'{print $4}\' | grep -E "$PATTERN")
            
            echo "$S3_OBJECTS"'
)

# ---- Get nextgen geopackages ----
# Run the script to get a list of the nextgen geopackages that matched the regular expression above
FEMA_BUCKET_KEYS <- system(fema_list_command, intern = TRUE)

# FEMA_BUCKET_OBJECTS <- paste0(FEMA_S3_BUCKET, FEMA_S3_BUCKET_PREFIX, FEMA_BUCKET_KEYS)

# Parse the selected S3 objects keys from the FEMA100 bucket directory copy them to the local destination directory if the file does NOT exist yet
for (key in FEMA_BUCKET_KEYS) {
  local_save_path <- paste0(FEMA_FGB_PATH, "/", key)
  
  if(!file.exists(local_save_path)) {
    copy_cmd <- paste0('aws s3 cp ', FEMA_S3_BUCKET, FEMA_S3_BUCKET_PREFIX, key, " ", local_save_path)
    
    message("Copying S3 object:\n", local_save_path)
    
    # system(copy_cmd)
    
    message("Download '", key, "' complete!")
    message("------------------")
  }
}