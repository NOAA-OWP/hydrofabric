# Running this script goes and pulls the desired NextGen geopackage datasets from http://www.lynker-spatial.com/, saves them into a directory within "base_dir"
# base_dir is defined within runners/workflow/root_dir.R

# load config variables
source("runners/cs_runner/config_vars.R")

# name of S3 bucket
s3_bucket <- "s3://lynker-spatial/"

# nextgen bucket folder name
nextgen_bucket_folder <- "v20.1/gpkg/"

# nextgen bucket name
nextgen_prefix  <- paste0(s3_bucket, nextgen_bucket_folder)

# prerelease_prefix     <- "s3://lynker-spatial/pre-release/"

# # reference features S3 bucket prefix
# ref_features_prefix <- "s3://lynker-spatial/00_reference_features/gpkg/"

# nextgen model attributes folder in S3 bucket with parquet files
model_attr_prefix <- paste0(s3_bucket, "v20/3D/model_attributes/")

# directory to copy nextgen bucket data too
nextgen_dir <- paste0(base_dir, "/", nextgen_bucket_folder)

# create the directory if it does NOT exist
if(!dir.exists(nextgen_dir)) {
  message("Directory does not exist at: \n\t'", nextgen_dir, "'\nCreating directory at: \n\t'", nextgen_dir, "'")
  
  dir.create(nextgen_dir)
}

# model attributes directory
model_attr_dir <- paste0(base_dir, "/model_attributes/")

# create the directory if it does NOT exist
if(!dir.exists(model_attr_dir)) {
  message("Directory does not exist at: \n\t'", model_attr_dir, "'\nCreating directory at: \n\t'", model_attr_dir, "'")
  dir.create(model_attr_dir)
}

# ---------------------------------------------------------------------------
# ---- List/Get the nextgen gpkgs from S3 bucket  ----
# ---------------------------------------------------------------------------

# list objects in S3 bucket, and regular expression match to nextgen_.gpkg pattern
command <- paste0('#!/bin/bash
            # AWS S3 Bucket and Directory information
            S3_BUCKET="', nextgen_prefix, '"
            DESTINATION_DIR=', nextgen_dir, '
            
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
  
  copy_cmd <- paste0('aws s3 cp ', nextgen_prefix, key, " ", nextgen_dir, key)
  message("Copying S3 object:\n", paste0(nextgen_prefix, key))
  
  system(copy_cmd)
  
  message("Download '", key, "' complete!")
  message("------------------")
}

# ---------------------------------------------------------------------------
# ---- List/Get nextgen model attributes parquets from S3 bucket ----
# ---------------------------------------------------------------------------

# aws s3 ls s3://lynker-spatial/v20/3D/model_attributes/

# list objects in S3 bucket, and regular expression match to nextgen_.gpkg pattern
list_model_attr_cmd <- paste0('#!/bin/bash
            # AWS S3 Bucket and Directory information
            S3_BUCKET="', model_attr_prefix, '"
            
            # Regular expression pattern to match object keys
            PATTERN="^nextgen_[0-9][0-9][A-Za-z]*\\_model_attributes.parquet$"
            
            # AWS CLI command to list objects in the S3 bucket and use grep to filter them
            S3_OBJECTS=$(aws s3 ls "$S3_BUCKET" | awk \'{print $4}\' | grep -E "$PATTERN")
            
            echo "$S3_OBJECTS"'
)

# get a list of the model attributes objects in S3
model_attr_keys <- system(list_model_attr_cmd, intern = TRUE)

# Parse the selected S3 objects keys and copy them to the destination directory
for (key in model_attr_keys) {
  
  copy_cmd <- paste0('aws s3 cp ', model_attr_prefix, key, ' ', model_attr_dir, key)
  
  message("Copying S3 object:\n", paste0(model_attr_prefix, key))
  
  system(copy_cmd)
  
  message("Download '", paste0(model_attr_prefix, key), "' complete!")
  message("------------------")
}

# ---------------------------------------------------------------------------
# ---- List/Get reference features from S3 bucket ----
# ---------------------------------------------------------------------------

## Go get a list of the reference features geopackages from S3 and create a save path using the S3 file names to save reference features to local directory

# list objects in S3 bucket, and regular expression match to nextgen_.gpkg pattern
list_ref_features <- paste0('#!/bin/bash
            # AWS S3 Bucket and Directory information
            S3_BUCKET="', ref_features_prefix , '"
            
            # Regular expression pattern to match object keys
            PATTERN="reference_features.gpkg"

            S3_OBJECTS=$(aws s3 ls "$S3_BUCKET" | awk \'{print $4}\' | grep -E "$PATTERN")
            
            echo "$S3_OBJECTS"'
)

# ---- Get a list of reference features geopackages geopackages ----
# Run the script to get a list of the nextgen geopackages that matched the regular expression above
ref_features <- system(list_ref_features, intern = TRUE)

## Download reference features geopackages and save them to a local directory
# Parse the selected S3 objects keys and copy them to the destination directory
for (key in ref_features) {
  # paste0(ref_features_dir, "gpkg/")
  copy_cmd <- paste0('aws s3 cp ', ref_features_prefix, key, ' ', paste0(ref_features_dir, "gpkg/"), key)
  
  message("Copying S3 object:\n", paste0(ref_features_prefix, key))
  system(copy_cmd)
  
  message("Download '", paste0(ref_features_dir, "gpkg/", key), "' complete!")
  message("------------------")
}

# ---------------------------------------------------------------------------
# ---- Get ML outputs data from S3 bucket ----
# ---------------------------------------------------------------------------

ml_copy_cmd <- paste0('aws s3 cp ', ML_OUTPUTS_URI, ' ', paste0(ML_OUTPUTS_DIR, basename(ML_OUTPUTS_URI)))

message("Copying S3 object:\n", ML_OUTPUTS_URI)
system(ml_copy_cmd)

message("Download '", paste0(ML_OUTPUTS_DIR, basename(ML_OUTPUTS_URI)), "' complete!")
message("------------------")

