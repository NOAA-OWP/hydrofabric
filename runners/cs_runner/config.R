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
# LYNKER_SPATIAL_HF_S3_URI <- "s3://lynker-spatial/"

# # name of bucket with nextgen data
# LYNKER_SPATIAL_S3_BUCKET_NAME <- "lynker-spatial"

# # the name of the folder in the S3 bucket with the nextgen data
# S3_BUCKET_NEXTGEN_DIR <- "v20.1/gpkg/"

# # full URI to the S3 bucket folder with the nextgen data 
# S3_BUCKET_NEXTGEN_DIR_URI  <- paste0(LYNKER_SPATIAL_HF_S3_URI, S3_BUCKET_NEXTGEN_DIR)

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
# ---- Get a list of reference features geopackages ----
# --------------------------------------------------------------------------

# Run the script to get a list of the nextgen geopackages that matched the regular expression above
ref_features <- system(list_ref_features, intern = TRUE)

# ref features datasets
ref_features_keys  <- paste0(S3_BUCKET_REF_FEATURES_URI, ref_features)
ref_features_files <- paste0(REF_FEATURES_DIR, "gpkg/", ref_features)

# --------------------------------------------------------------------------
# ---- Create empty file structure for a "new_domain" ----
# --------------------------------------------------------------------------

create_new_domain_dirs(BASE_DIR, NEW_DOMAIN_DIRNAME, with_output = TRUE)

# --------------------------------------------------------------------------
# ---- Create empty file structure for a "domain_with_fema" ----
# --------------------------------------------------------------------------

create_new_domain_dirs(BASE_DIR, DOMAIN_WITH_FEMA_DIRNAME, with_output = TRUE)

# --------------------------------------------------------------------------
# ---- Get locations of ML parquet files in S3 ---
# --------------------------------------------------------------------------

VPU_ML_BATHYMETRY_URIS <- unlist(
  lapply(VPU_ML_BATHYMETRY_S3_DIRS, function(s3_dir) {
    s3_file <- list_s3_objects(s3_dir, ".parquet$", AWS_PROFILE)
    paste0(s3_dir, s3_file)
  })
)

# -------------------------------------------------------------------------------------
# ---- Download ML parquets for DOMAIN_WITH_FEMA ----
# -------------------------------------------------------------------------------------

# Parse the selected S3 objects keys from the FEMA100 bucket directory copy them to the local destination directory if the file does NOT exist yet
for (s3_uri in VPU_ML_BATHYMETRY_URIS) {
  # message(s3_uri)
  # s3_uri <- "s3://lynker-hydrofabric/hydrofabric/nextgen/bathymetry/multisource_river_attributes/vpuid=18/part-0.parquet"
  # s3_uri <- "s3://lynker-hydrofabric/hydrofabric/nextgen/bathymetry/multisource_river_attributes/vpuid=03W/part-0.parquet"
  # s3_uri <- "s3://lynker-hydrofabric/hydrofabric/nextgen/bathymetry/multisource_river_attributes/vpuid=21/"
  
  is_parquet <- endsWith(s3_uri, ".parquet")
  vpu_id     <- gsub(".*vpuid=([a-zA-Z0-9]+).*", "\\1", s3_uri)
  
  message("Checking S3 bucket for VPU ", vpu_id, " ML data...")
  
  if (is_parquet) {
    s3_file    <- basename(s3_uri)
    # vpu_id <- gsub(".*vpuid=([a-zA-Z0-9]+).*", "\\1", s3_uri)
    new_file_name <-     paste0(vpu_id, "_ml.parquet")
    # new_file_name <- gsub("-", "_", paste0(vpu_id, "_", s3_file))
    
    local_save_path <- paste0(DOMAIN_WITH_FEMA_ML_DIR, "/", new_file_name)
    
    if(!file.exists(local_save_path)) {
      copy_cmd <- paste0('aws s3 cp ', s3_uri,
                         " ", local_save_path, " --profile ", AWS_PROFILE)
      
      message("S3 object:\n > '", s3_uri, "'")
      message("Downloading S3 object to:\n > '", local_save_path, "'")
      # message("Copy command:\n > '", copy_cmd, "'")
      
      system(copy_cmd)
      
      message(" > '", new_file_name, "' download complete!")
      message("----------------------------------")
    } else {
      message("File already exists at:\n > '", local_save_path, "'")
    }
    
  } else {
    message("No S3 bucket ML data for VPU ", vpu_id, "...")
  }

}

# VPU_ML_BATHYMETRY_PATHS <- list.files(DOMAIN_WITH_FEMA_ML_DIR, full.names = T)
# 
# ml_outputs <- lapply(VPU_ML_BATHYMETRY_PATHS, function(prq) {
#       vpu_id     <- gsub(".*ml/([a-zA-Z0-9]+).*", "\\1", prq)
#       arrow::read_parquet(prq) %>% 
#         dplyr::mutate(vpu_id = vpu_id) 
#       }
#     ) %>% 
#   dplyr::bind_rows()

# --------------------------------------------------------------------------
# ---- Get locations of diffusive domain DEM files in S3 ----
# --------------------------------------------------------------------------

COASTAL_BATHY_DEM_S3_URIS    <- paste0(COASTAL_BATHY_DEM_S3_DIR_URI, 
                                       list_s3_objects(COASTAL_BATHY_DEM_S3_DIR_URI, ".tif$", AWS_PROFILE)
                                       )

# -------------------------------------------------------------------------------------
# ---- Download diffusive domain DEM files from S3 for DOMAIN_WITH_FEMA ----
# -------------------------------------------------------------------------------------

# Parse the selected S3 objects keys from the FEMA100 bucket directory copy them to the local destination directory if the file does NOT exist yet
for (s3_uri in COASTAL_BATHY_DEM_S3_URIS) {
  message(s3_uri)
  
  is_tif <- endsWith(s3_uri, ".tif")
  # vpu_id     <- gsub(".*vpuid=([a-zA-Z0-9]+).*", "\\1", s3_uri)
  
  message("Checking S3 bucket for DEM data...")
  
  if (is_tif) {
    
    s3_file    <- basename(s3_uri)
    # vpu_id <- gsub(".*vpuid=([a-zA-Z0-9]+).*", "\\1", s3_uri)
    # new_file_name <-     paste0(vpu_id, "_ml.parquet")
    # new_file_name <- gsub("-", "_", paste0(vpu_id, "_", s3_file))
    
    local_save_path <- paste0(DOMAIN_WITH_FEMA_DEM_DIR, "/", s3_file)
    
    if (!file.exists(local_save_path)) {
      copy_cmd <- paste0('aws s3 cp ', s3_uri,
                         " ", local_save_path, " --profile ", AWS_PROFILE)
      
      message("S3 object:\n > '", s3_uri, "'")
      message("Downloading S3 object to:\n > '", local_save_path, "'")
      # message("Copy command:\n > '", copy_cmd, "'")
      
      system(copy_cmd)
      
      message(" > '", s3_file, "' download complete!")
      message("----------------------------------")
    } else {
      message("File already exists at:\n > '", local_save_path, "'")
    }
    
  } else {
    message("No S3 bucket DEM data... ")
  }
  
}

COASTAL_BATHY_DEM_PATHS <- list.files(DOMAIN_WITH_FEMA_DEM_DIR, full.names = TRUE)




























