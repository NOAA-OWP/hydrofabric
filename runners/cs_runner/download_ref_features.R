# download_ref_features.R
# This script pulls the reference features geopackage datasets from http://www.lynker-spatial.com/, saves them into a directory within "BASE_DIR"

# load config variables
source("runners/cs_runner/config_vars.R")

# ---------------------------------------------------------------------------
# ---- List/Get reference features from S3 bucket ----
# ---------------------------------------------------------------------------

# list objects in S3 bucket, and regular expression match to reference_features.gpkg pattern
list_ref_features <- paste0('#!/bin/bash
            # AWS S3 Bucket and Directory information
            S3_BUCKET="', S3_BUCKET_REF_FEATURES_URI, '"
            
            # Regular expression pattern to match object keys
            PATTERN="reference_features.gpkg"

            S3_OBJECTS=$(aws s3 ls "$S3_BUCKET" | awk \'{print $4}\' | grep -E "$PATTERN")
            
            echo "$S3_OBJECTS"'
)

# ---- Get a list of reference features geopackages ----
# Run the script to get a list of the reference features geopackages that matched the regular expression above
ref_features <- system(list_ref_features, intern = TRUE)

# Parse the selected S3 objects keys and copy them to the destination directory
for (key in ref_features) {
  copy_cmd <- paste0('aws s3 cp ', S3_BUCKET_REF_FEATURES_URI, key, ' ', paste0(REF_FEATURES_DIR, "gpkg/"), key)
  message("Copying S3 object:\n", paste0(S3_BUCKET_REF_FEATURES_URI, key))
  system(copy_cmd)
  message("Download '", paste0(REF_FEATURES_DIR, "gpkg/", key), "' complete!")
  message("------------------")
}
