# Running this script goes and pulls the desired NextGen geopackage datasets from http://www.lynker-spatial.com/, saves them into a directory within "base_dir"
# base_dir is defined within runners/workflow/root_dir.R

library(logger)

source("runners/workflow/root_dir.R")

# name of S3 bucket
s3_bucket <- "s3://lynker-spatial/"

# nextgen bucket name
prerelease_prefix     <- "s3://lynker-spatial/pre-release/"

# nextgen model attributes folder in S3 bucket with parquet files
model_attr_prefix    <- glue::glue("{s3_bucket}v20/3D/model_attributes/")

# directory to copy nextgen bucket data too
nextgen_dir   <- glue::glue('{base_dir}/pre-release/')

# create the directory if it does NOT exist
if(!dir.exists(nextgen_dir)) {
  logger::log_info("\n\nDirectory does not exist at: \n\t'{nextgen_dir}'\nCreating directory at: \n\t'{nextgen_dir}'")
  dir.create(nextgen_dir)
}

# model attributes directory
model_attr_dir <- glue::glue('{base_dir}/model_attributes/')

# create the directory if it does NOT exist
if(!dir.exists(model_attr_dir)) {
  logger::log_info("\n\nDirectory does not exist at: \n\t'{model_attr_dir}'\nCreating directory at: \n\t'{model_attr_dir}'")
  dir.create(model_attr_dir)
}

# list objects in S3 bucket, and regular expression match to nextgen_.gpkg pattern
command <- paste0('#!/bin/bash
            # AWS S3 Bucket and Directory information
            S3_BUCKET="', prerelease_prefix, '"
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
  
  copy_cmd <- paste0('aws s3 cp ', prerelease_prefix, key, " ", nextgen_dir, key)
  logger::log_info("Copying S3 object:\n{prerelease_prefix}{key}")
  
  system(copy_cmd)
  
  logger::log_info("Download '{key}' complete!")
  logger::log_info("------------------")
}

# ---- Get nextgen model attributes parquets ----
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
  
  copy_cmd <- glue::glue("aws s3 cp {model_attr_prefix}{key} {model_attr_dir}{key}")
  # copy_cmd <- paste0('aws s3 cp ', prerelease_prefix, key, " ", nextgen_dir, key)

  logger::log_info("Copying S3 object:\n{model_attr_prefix}{key}")
  
  system(copy_cmd)
  
  logger::log_info("Download '{model_attr_prefix}{key}' complete!")

  logger::log_info("------------------")
}


