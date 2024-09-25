# load required packages
pacman::p_load(
  archive,
  hydrofabric,
  hydrofabric3D
)

# # install.packages("devtools")
# devtools::install_github("anguswg-ucsb/hydrofabric3D")

# load root directory 
source("runners/cs_runner/config_vars.R")
source("runners/cs_runner/utils.R")

sf::sf_use_s2(FALSE)

### Cross section point

# # -------------------------------------------------------------------------------------
# # ----- S3 names ------
# # -------------------------------------------------------------------------------------

# # AWS S3 bucket URI 
# S3_BUCKET_URI <- "s3://lynker-spatial/"

# # name of bucket with nextgen data
# S3_BUCKET_NAME <- "lynker-spatial"

# # the name of the folder in the S3 bucket with the nextgen data
# S3_BUCKET_NEXTGEN_DIR <- "v20.1/gpkg/"

# # full URI to the S3 bucket folder with the nextgen data 
# S3_BUCKET_NEXTGEN_DIR_URI  <- paste0(S3_BUCKET_URI, S3_BUCKET_NEXTGEN_DIR)

# # reference features S3 bucket prefix
# S3_BUCKET_REF_FEATURES_URI <- "s3://lynker-spatial/00_reference_features/gpkg/"

# # S3 prefix/folder of version run
# VERSION <- "v20.1"

# # -------------------------------------------------------------------------------------

# # -------------------------------------------------------------------------------------
# # ----- Local directories ------
# # -------------------------------------------------------------------------------------

# ### LOCAL DIRS

# # directory to copy nextgen bucket data too
# NEXTGEN_DIR      <- paste0(BASE_DIR, "/", S3_BUCKET_NEXTGEN_DIR)
# # NEXTGEN_DIR <- paste0(BASE_DIR, "/pre-release/")

# # model attributes directory
# MODEL_ATTR_DIR   <- paste0(BASE_DIR, "/model_attributes/")

# # cross-section data model data directories
# TRANSECTS_DIR    <- paste0(BASE_DIR, "/01_transects/")
# CS_PTS_DIR       <- paste0(BASE_DIR, "/02_cs_pts/")

# # final output directory with geopackages per VPU
# CS_OUTPUT_DIR    <- paste0(BASE_DIR, "/cross_sections/")

# # directory to copy nextgen bucket data too
# REF_FEATURES_DIR <- paste0(BASE_DIR, "/00_reference_features/")

# # make a directory for the ML outputs data
# ML_OUTPUTS_DIR   <- paste0(BASE_DIR, "/ml-outputs/")

# -------------------------------------------------------------------------------------
# ----- Create local directories ------
# -------------------------------------------------------------------------------------
# create directories 
dir.create(TRANSECTS_DIR, showWarnings = FALSE)
dir.create(CS_PTS_DIR,    showWarnings = FALSE)
dir.create(REF_FEATURES_DIR,     showWarnings = FALSE)
dir.create(paste0(REF_FEATURES_DIR, "gpkg/"),     showWarnings = FALSE)
dir.create(CS_OUTPUT_DIR,     showWarnings = FALSE)
dir.create(ML_OUTPUTS_DIR,     showWarnings = FALSE)

# create the directory if it does NOT exist
if(!dir.exists(NEXTGEN_DIR)) {
  message("Directory does not exist at: \n\t'", NEXTGEN_DIR, "'\nCreating directory at: \n\t'", NEXTGEN_DIR, "'")
  
  dir.create(NEXTGEN_DIR)
}

# # create the directory if it does NOT exist
# if(!dir.exists(MODEL_ATTR_DIR)) {
#   message("Directory does not exist at: \n\t'", MODEL_ATTR_DIR, "'\nCreating directory at: \n\t'", MODEL_ATTR_DIR, "'")
#   dir.create(MODEL_ATTR_DIR)
# }
# dir.create(MODEL_ATTR_DIR,  showWarnings = FALSE)

# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# ----- Get the paths / locations of reference_features data ------
# -------------------------------------------------------------------------------------

## Go get a list of the reference features geopackages from S3 and create a save path using the S3 file names to save reference features to local directory

# list objects in S3 bucket, and regular expression match to nextgen_.gpkg pattern
list_ref_features <- paste0('#!/bin/bash
            # AWS S3 Bucket and Directory information
            S3_BUCKET="', S3_BUCKET_REF_FEATURES_URI , '"
            
            # Regular expression pattern to match object keys
            PATTERN="reference_features.gpkg"

            S3_OBJECTS=$(aws s3 ls "$S3_BUCKET" | awk \'{print $4}\' | grep -E "$PATTERN")
            
            echo "$S3_OBJECTS"'
)

# --------------------------------------------------------------------------
# ---- Create empty file structure for a "new_domain" ----
# --------------------------------------------------------------------------

create_new_domain_dirs(BASE_DIR, NEW_DOMAIN_DIRNAME)

# --------------------------------------------------------------------------
# ---- Create empty file structure for a "domain_with_fema" ----
# --------------------------------------------------------------------------

create_new_domain_dirs(BASE_DIR, DOMAIN_WITH_FEMA_DIRNAME)

# --------------------------------------------------------------------------
# ---- Create empty file structure for a "new_conus_domain" ----
# --------------------------------------------------------------------------

create_new_domain_dirs(BASE_DIR, NEW_CONUS_DOMAIN_DIRNAME)

# --------------------------------------------------------------------------
# ---- Get a list of reference features geopackages ----
# --------------------------------------------------------------------------

# Run the script to get a list of the nextgen geopackages that matched the regular expression above
ref_features <- system(list_ref_features, intern = TRUE)

# ref features datasets
ref_features_keys  <- paste0(S3_BUCKET_REF_FEATURES_URI, ref_features)
ref_features_files <- paste0(REF_FEATURES_DIR, "gpkg/", ref_features)
