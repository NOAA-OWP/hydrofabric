# Running this script goes and pulls the desired NextGen geopackage datasets from http://www.lynker-spatial.com/, saves them into a directory within "BASE_DIR"
# BASE_DIR is defined within runners/workflow/root_dir.R

# load config variables
source("runners/cs_runner/config_vars.R")

# directory to copy nextgen bucket data too
NEXTGEN_DIR <- paste0(BASE_DIR, "/", S3_BUCKET_NEXTGEN_DIR)

# create the directory if it does NOT exist
if(!dir.exists(NEXTGEN_DIR)) {
  message("Directory does not exist at: \n\t'", NEXTGEN_DIR, "'\nCreating directory at: \n\t'", NEXTGEN_DIR, "'")
  
  dir.create(NEXTGEN_DIR)
}

# ---------------------------------------------------------------------------
# ---- List/Get the nextgen gpkgs from S3 bucket  ----
# ---------------------------------------------------------------------------

# list objects in S3 bucket, and regular expression match to nextgen_.gpkg pattern
command <- paste0('#!/bin/bash
            # AWS S3 Bucket and Directory information
            S3_BUCKET="', S3_BUCKET_NEXTGEN_DIR_URI, '"
            DESTINATION_DIR=', NEXTGEN_DIR, '
            
            # Regular expression pattern to match object keys
            PATTERN="^nextgen_[0-9][0-9][A-Za-z]*\\.gpkg$"
            
            # AWS CLI command to list objects in the S3 bucket and use grep to filter them
            S3_OBJECTS=$(aws s3 ls "$S3_BUCKET" | awk \'{print $4}\' | grep -E "$PATTERN")
            
            echo "$S3_OBJECTS"'
                  )

# ---- Get nextgen geopackages ----
# Run the script to get a list of the nextgen geopackages that matched the regular expression above
bucket_keys <- system(command, intern = TRUE)

# Parse the selected S3 objects keys and copy them to the destination directory
for (key in bucket_keys) {
  
  copy_cmd <- paste0('aws s3 cp ', S3_BUCKET_NEXTGEN_DIR_URI, key, " ", NEXTGEN_DIR, key)
  message("Copying S3 object:\n", paste0(S3_BUCKET_NEXTGEN_DIR_URI, key))
  
  system(copy_cmd)
  
  message("Download '", key, "' complete!")
  message("------------------")
}

# # ---------------------------------------------------------------------------
# # ---- List/Get reference features from S3 bucket ----
# # ---------------------------------------------------------------------------
# 
# ## Go get a list of the reference features geopackages from S3 and create a save path using the S3 file names to save reference features to local directory
# 
# # list objects in S3 bucket, and regular expression match to nextgen_.gpkg pattern
# list_ref_features <- paste0('#!/bin/bash
#             # AWS S3 Bucket and Directory information
#             S3_BUCKET="', S3_BUCKET_REF_FEATURES_URI , '"
#             
#             # Regular expression pattern to match object keys
#             PATTERN="reference_features.gpkg"
# 
#             S3_OBJECTS=$(aws s3 ls "$S3_BUCKET" | awk \'{print $4}\' | grep -E "$PATTERN")
#             
#             echo "$S3_OBJECTS"'
# )
# 
# # ---- Get a list of reference features geopackages geopackages ----
# # Run the script to get a list of the nextgen geopackages that matched the regular expression above
# ref_features <- system(list_ref_features, intern = TRUE)
# 
# ## Download reference features geopackages and save them to a local directory
# # Parse the selected S3 objects keys and copy them to the destination directory
# for (key in ref_features) {
#   # paste0(REF_FEATURES_DIR, "gpkg/")
#   copy_cmd <- paste0('aws s3 cp ', S3_BUCKET_REF_FEATURES_URI, key, ' ', paste0(REF_FEATURES_DIR, "gpkg/"), key)
#   
#   message("Copying S3 object:\n", paste0(S3_BUCKET_REF_FEATURES_URI, key))
#   system(copy_cmd)
#   
#   message("Download '", paste0(REF_FEATURES_DIR, "gpkg/", key), "' complete!")
#   message("------------------")
# }
# 
# # ---------------------------------------------------------------------------
# # ---- Get ML outputs data from S3 bucket ----
# # ---------------------------------------------------------------------------
# 
# ml_copy_cmd <- paste0('aws s3 cp ', ML_OUTPUTS_S3_URI, ' ', paste0(ML_OUTPUTS_DIR, basename(ML_OUTPUTS_S3_URI)))
# 
# message("Copying S3 object:\n", ML_OUTPUTS_S3_URI)
# system(ml_copy_cmd)
# 
# message("Download '", paste0(ML_OUTPUTS_DIR, basename(ML_OUTPUTS_S3_URI)), "' complete!")
# message("------------------")
# 
