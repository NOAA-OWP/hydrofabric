# Running this script goes and pulls the desired FEMA100 flood fgb datasets from the lynker-hydrofabric S3 bucket then saves them into a directory within "base_dir"
# base_dir is defined within runners/workflow/root_dir.R

# NOTE: The lynker-hydrofabric S3 bucket is private at the moment

# load config variables
source("runners/cs_runner/config_vars.R")

# -------------------------------------------------------------------------------------
# ---- Create FEMA100/ directory and bounding box dir (if it does NOT exist) ----
# -------------------------------------------------------------------------------------
# create FEMA FGB directory (if not exists) 
if (!dir.exists(FEMA_FGB_PATH)) {
  message(paste0("FEMA100/ directory does not exist...\nCreating directory:\n > '", FEMA_FGB_PATH, "'"))
  dir.create(FEMA_FGB_PATH)
}

# create geojsons directory (if not exists) 
if (!dir.exists(FEMA_GEOJSON_PATH)) {
  message(paste0(FEMA_GEOJSON_PATH, " directory does not exist...\nCreating directory:\n > '", FEMA_GEOJSON_PATH, "'"))
  dir.create(FEMA_GEOJSON_PATH)
}

# create directory for cleaned FEMA geometries (if not exists) 
if (!dir.exists(FEMA_CLEAN_PATH)) {
  message(paste0(FEMA_CLEAN_PATH, " directory does not exist...\nCreating directory:\n > '", FEMA_CLEAN_PATH, "'"))
  dir.create(FEMA_CLEAN_PATH)
}

# create directory for cleaned FEMA geometries as geopackages (if not exists) 
if (!dir.exists(FEMA_GPKG_PATH)) {
  message(paste0(FEMA_GPKG_PATH, " directory does not exist...\nCreating directory:\n > '", FEMA_GPKG_PATH, "'"))
  dir.create(FEMA_GPKG_PATH)
}

# create simplified geojsons directory (if not exists)
if (!dir.exists(FEMA_SIMPLIFIED_PATH)) {
  message(paste0(FEMA_SIMPLIFIED_PATH, " directory does not exist...\nCreating directory:\n > '", FEMA_SIMPLIFIED_PATH, "'"))
  dir.create(FEMA_SIMPLIFIED_PATH)
}

# create simplified geojsons directory (if not exists)
if (!dir.exists(FEMA_DISSOLVED_PATH)) {
  message(paste0(FEMA_DISSOLVED_PATH, " directory does not exist...\nCreating directory:\n > '", FEMA_DISSOLVED_PATH, "'"))
  dir.create(FEMA_DISSOLVED_PATH)
}

# create exploded geojsons directory (if not exists)
if (!dir.exists(FEMA_EXPLODED_PATH)) {
  message(paste0(FEMA_EXPLODED_PATH, " directory does not exist...\nCreating directory:\n > '", FEMA_EXPLODED_PATH, "'"))
  dir.create(FEMA_EXPLODED_PATH)
}

# create FEMA GPKG Bounding Boxes directory (if not exists)
if (!dir.exists(FEMA_GPKG_BB_PATH)) {
  message(paste0(FEMA_GPKG_BB_PATH, " directory does not exist...\nCreating directory:\n > '", FEMA_GPKG_BB_PATH, "'"))
  dir.create(FEMA_GPKG_BB_PATH)
}


# -------------------------------------------------------------------------------------
# ---- Get list of FEMA FGB files in S3 bucket ----
# -------------------------------------------------------------------------------------

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

# -------------------------------------------------------------------------------------
# ---- Get the S3 buckets object keys for FEMA 100 FGB files ----
# -------------------------------------------------------------------------------------

# Run the script to get a list of the nextgen geopackages that matched the regular expression above
FEMA_BUCKET_KEYS <- system(fema_list_command, intern = TRUE)

# create bucket object URIs
# FEMA_BUCKET_OBJECTS <- paste0(FEMA_S3_BUCKET, FEMA_S3_BUCKET_PREFIX, FEMA_BUCKET_KEYS)

# -------------------------------------------------------------------------------------
# ---- Download FEMA 100 year FGB files from S3 ----
# -------------------------------------------------------------------------------------

# Parse the selected S3 objects keys from the FEMA100 bucket directory copy them to the local destination directory if the file does NOT exist yet
for (key in FEMA_BUCKET_KEYS[1:length(FEMA_BUCKET_KEYS)]) {
  local_save_path <- paste0(FEMA_FGB_PATH, "/", key)
  
  if(!file.exists(local_save_path)) {
    copy_cmd <- paste0('aws s3 cp ', FEMA_S3_BUCKET, FEMA_S3_BUCKET_PREFIX, key, " ", local_save_path, " --profile ", aws_profile)
    
    message("S3 object:\n > '", FEMA_S3_BUCKET, FEMA_S3_BUCKET_PREFIX, key, "'")
    message("Downloading S3 object to:\n > '", local_save_path, "'")
    # message("Copy command:\n > '", copy_cmd, "'")
    
    # system(copy_cmd)
    
    message(" > '", key, "' download complete!")
    message("----------------------------------")
  } else {
    message("File already exists at:\n > '", local_save_path, "'")
  }
}
# -------------------------------------------------------------------------------------
# ---- Run ogr2ogr to get FGB files into geojson ----
# -------------------------------------------------------------------------------------

for (key in FEMA_BUCKET_KEYS) {
  
  local_fema_path <- paste0(FEMA_FGB_PATH, "/", key)
  
  geojson_filename     <- gsub(".fgb", ".geojson", key)
  geojson_save_path    <- paste0(FEMA_GEOJSON_PATH, "/", geojson_filename)
  
  message("S3 Key: '", key, "'")
  message("Converting \n > '", key, "' to geojson '", geojson_filename, "'")
  
  ogr2ogr_command = paste0("ogr2ogr ", geojson_save_path, " ", local_fema_path)
  
  system(ogr2ogr_command)
  
  message("Saved '", geojson_filename, "' saved to: \n > '", geojson_save_path, "'")
  message()
}  

# -------------------------------------------------------------------------------------
# ---- Clean FEMA geometries (Simplify, Dissolve, Explode) ----
# -------------------------------------------------------------------------------------

# paths to FEMA 100 year flood plain files
FEMA_geojson_paths      <- list.files(FEMA_GEOJSON_PATH, full.names = TRUE)
# FEMA_BB_paths   <- list.files(FEMA_GPKG_BB_PATH, full.names = TRUE)

for (fema_file in FEMA_geojson_paths) {
  message("Fema 100 year flood plain:\n > '", basename(fema_file), "'")
  # message("Fema 100 year flood plain:\n > '", fema_file, "'")
  output_clean_filename <- gsub(".geojson", "_clean.geojson", basename(fema_file))
  output_path <- paste0(FEMA_CLEAN_PATH, "/", output_clean_filename)
  
  message("Running mapshaper 'simplify', 'dissolve', and 'explode' via CLI...")
  
  # mapshaper_command = paste0('node  --max-old-space-size=16000 /opt/homebrew/bin/mapshaper ', fema_file, ' -simplify 0.15 visvalingam -o ', output_path)
  # test_file_path <- "/Users/anguswatters/Desktop/lynker-spatial/FEMA100_simplified/Wyoming-100yr-flood_valid_clean.geojson"
  mapshaper_command = paste0('node  --max-old-space-size=16000 /opt/homebrew/bin/mapshaper ', fema_file, 
                             ' -simplify 0.15 visvalingam \\', 
                             ' -dissolve \\', 
                             ' -explode \\', 
                             ' -o ', output_path
                             )
  system(mapshaper_command)
  message("Mapshaper command: ", mapshaper_command)
  message()
}

# -------------------------------------------------------------------------------------
# ---- Convert cleaned FEMA geometries to geopackages ----
# -------------------------------------------------------------------------------------

# paths to FEMA 100 year flood plain files
FEMA_clean_paths      <- list.files(FEMA_CLEAN_PATH, full.names = TRUE)

for (fema_file in FEMA_clean_paths) {
  message("Fema 100 year flood plain:\n > '", basename(fema_file), "'")
  # message("Fema 100 year flood plain:\n > '", fema_file, "'")
  output_gpkg_filename <- gsub("_clean.geojson", "_clean.gpkg", basename(fema_file))
  output_path <- paste0(FEMA_GPKG_PATH, "/", output_gpkg_filename)
  
  message("Converting geojson files to gpkg...")
  
  message("Converting \n > '", fema_file, "' to geojson '", output_gpkg_filename, "'")
  
  # mapshaper_command = paste0('node  --max-old-space-size=16000 /opt/homebrew/bin/mapshaper ', fema_file, ' -simplify 0.15 visvalingam -o ', output_path)
  # test_file_path <- "/Users/anguswatters/Desktop/lynker-spatial/FEMA100_simplified/Wyoming-100yr-flood_valid_clean.geojson"
  ogr2ogr_command = paste0("ogr2ogr -nlt MULTIPOLYGON ", output_path, " ", fema_file)
  
  # system(ogr2ogr_command)
  
  message("Saved '", output_gpkg_filename, "' saved to: \n > '", output_path, "'")
  message()
}

# # -------------------------------------------------------------------------------------
# # ---- Apply hydrofab::clean_geometries() to cleaned FEMA geometries  ----
# # -------------------------------------------------------------------------------------
# 
# # paths to FEMA 100 year flood plain files
# FEMA_gpkg_paths      <- list.files(FEMA_GPKG_PATH, full.names = TRUE)
# 
# for (fema_file in FEMA_gpkg_paths) {
#   message("Applying final cleaning process to:\n > '", basename(fema_file), "'")
#   
#   fema <- sf::read_sf(fema_file)
#   
#   fema
#   
#   # message("Fema 100 year flood plain:\n > '", fema_file, "'")
#   output_gpkg_filename <- gsub("_clean.geojson", "_clean.gpkg", basename(fema_file))
#   output_path <- paste0(FEMA_GPKG_PATH, "/", output_gpkg_filename)
#   
#   message("Converting geojson files to gpkg...")
#   
#   message("Converting \n > '", fema_file, "' to geojson '", output_gpkg_filename, "'")
#   
#   # mapshaper_command = paste0('node  --max-old-space-size=16000 /opt/homebrew/bin/mapshaper ', fema_file, ' -simplify 0.15 visvalingam -o ', output_path)
#   # test_file_path <- "/Users/anguswatters/Desktop/lynker-spatial/FEMA100_simplified/Wyoming-100yr-flood_valid_clean.geojson"
#   ogr2ogr_command = paste0("ogr2ogr -nlt MULTIPOLYGON ", output_path, " ", fema_file)
#   
#   system(ogr2ogr_command)
#   
#   message("Saved '", output_gpkg_filename, "' saved to: \n > '", output_path, "'")
#   message()
# }

# -------------------------------------------------------------------------------------
# ---- Generate bounding box gpkg for each FEMA FGB ----
# -------------------------------------------------------------------------------------

for (key in FEMA_BUCKET_KEYS) {
 
  local_fema_path <- paste0(FEMA_FGB_PATH, "/", key)
  
  gpkg_filename   <- gsub(".fgb", "_bb.gpkg", key)
  bb_save_path    <- paste0(FEMA_GPKG_BB_PATH, "/", gpkg_filename)
  
  message("S3 Key: '", key, "'")
  message("Local FEMA file:\n > '", local_fema_path, "'")
  message("Local output FEMA bounding box file:\n > '", bb_save_path, "'")
  
  # fema <- sf::read_sf(local_fema_path)
  
  fema_bb <- 
    local_fema_path %>% 
    sf::read_sf() %>% 
    sf::st_bbox() %>% 
    sf::st_as_sfc() %>% 
    sf::st_as_sf() %>% 
    dplyr::mutate(
      fema_fgb      = key,
      fema_fgb_path = local_fema_path,
      state         = gsub("-100yr-flood_valid.fgb", "", key)
    ) %>% 
    dplyr::select(fema_fgb, fema_fgb_path, state, geometry = x) %>% 
    sf::st_transform(5070)
  
  message("Saving FEMA bounding box file:\n > '", bb_save_path, "'")
  
  sf::write_sf(fema_bb, bb_save_path)
  message()
}  
  




