# Script should be run AFTER download_fema100.R as the FEMA 100 year flood plain data needs to first be downloaded from S3
# This file will take a directory of FEMA 100 year FGB files (FEMA_FGB_PATH) the below processes to generate a cleaned, simple set of geopackages 

# Processing steps:
# - Convert FGBs to GEOJSON (via ogr2ogr)
# - Simplifies
# - Dissolves
# - Explodes
# - Convert cleaned GEOJSON to cleaned GPKGs (via ogr2ogr)
# - Apply hydrofab::clean_geometry()
# - Partition FEMA 100 geometries by VPU      # TODO still
# - Get FEMA bounding box geometries (maybe)

# load config variables
source("runners/cs_runner/config_vars.R")
source("runners/cs_runner/config.R")
source("runners/cs_runner/utils.R")

library(dplyr)
library(sf)
library(geos)

# -------------------------------------------------------------------------------------
# ---- OVERWRITE_FEMA_FILES constant logical ----
# ---- > if TRUE, processing steps will be run again 
#          and overwrite existing previously processed files
# -------------------------------------------------------------------------------------

# Default is TRUE (i.e. a fresh processing run is done from start to finish)
OVERWRITE_FEMA_FILES <- TRUE

# -------------------------------------------------------------------------------------
# ---- Create directories (if they do NOT exist) ----
# -------------------------------------------------------------------------------------

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

# create directory for FEMA geomteries partioned by VPU
if (!dir.exists(FEMA_BY_VPU_PATH)) {
  message(paste0(FEMA_BY_VPU_PATH, " directory does not exist...\nCreating directory:\n > '", FEMA_BY_VPU_PATH, "'"))
  dir.create(FEMA_BY_VPU_PATH)
}

for (VPU_SUBFOLDER in FEMA_VPU_SUBFOLDERS) {
  # create directory for FEMA geomteries by VPU
  if (!dir.exists(VPU_SUBFOLDER)) {
    message("Creating FEMA VPU subfolder...")
    message(paste0("'/", basename(VPU_SUBFOLDER), "' directory does not exist...\nCreating directory:\n > '", VPU_SUBFOLDER, "'"))
    dir.create(VPU_SUBFOLDER)
  }
}

# FEMA_VPU_SUBFOLDERS <- paste0(FEMA_BY_VPU_PATH, "/VPU_", VPU_IDS)


# create FEMA GPKG Bounding Boxes directory (if not exists)
if (!dir.exists(FEMA_GPKG_BB_PATH)) {
  message(paste0(FEMA_GPKG_BB_PATH, " directory does not exist...\nCreating directory:\n > '", FEMA_GPKG_BB_PATH, "'"))
  dir.create(FEMA_GPKG_BB_PATH)
}

# -------------------------------------------------------------------------------------
# ---- Get paths to downloaded FEMA 100 FGBs ----
# -------------------------------------------------------------------------------------

FEMA_FILENAMES  <- list.files(FEMA_FGB_PATH, full.names = FALSE)
FEMA_FILE_PATHS <- paste0(FEMA_FGB_PATH, "/", FEMA_FILENAMES)

# -------------------------------------------------------------------------------------
# ---- Run ogr2ogr to get FGB files into geojson ----
# -------------------------------------------------------------------------------------

for (file in FEMA_FILENAMES) {
  
  local_fema_path <- paste0(FEMA_FGB_PATH, "/", file)
  
  geojson_filename     <- gsub(".fgb", ".geojson", file)
  geojson_save_path    <- paste0(FEMA_GEOJSON_PATH, "/", geojson_filename)
  
  message("FEMA filename: '", file, "'")
  message("Converting \n > '", file, "' to geojson '", geojson_filename, "'")
  
  geojson_exists <- file.exists(geojson_save_path)
  
  message(" >>> '", geojson_filename, "' already exists? ", geojson_exists)
  message(" >>> Overwrite? ", OVERWRITE_FEMA_FILES)
  
  # ogr2ogr command converting FGBs to GEOJSON for mapshaper processing
  ogr2ogr_command = paste0("ogr2ogr ", geojson_save_path, " ", local_fema_path)
  
  if (OVERWRITE_FEMA_FILES) {
    system(ogr2ogr_command)
    message("Writting '", geojson_filename, "' to: \n > '", geojson_save_path, "'")
  }
  
  message()
}  

# -------------------------------------------------------------------------------------
# ---- Clean FEMA geometries (Simplify, Dissolve, Explode) ----
# -------------------------------------------------------------------------------------

# paths to FEMA 100 year flood plain files
FEMA_geojson_paths      <- list.files(FEMA_GEOJSON_PATH, full.names = TRUE)

for (file in FEMA_geojson_paths) {
  
  message("Simplify, dissolve, explode > '", basename(file), "'")
  # message("Fema 100 year flood plain:\n > '", file, "'")
  output_clean_filename <- gsub(".geojson", "_clean.geojson", basename(file))
  output_path <- paste0(FEMA_CLEAN_PATH, "/", output_clean_filename)
  
  clean_geojson_exists <- file.exists(output_path)
  message(" >>> '", output_clean_filename, "' already exists? ", clean_geojson_exists)
  message(" >>> Overwrite? ", OVERWRITE_FEMA_FILES)
  
  mapshaper_command = paste0('node  --max-old-space-size=16000 /opt/homebrew/bin/mapshaper ', file, 
                             ' -simplify 0.15 visvalingam \\', 
                             ' -dissolve \\', 
                             ' -explode \\', 
                             ' -o ', output_path
  )
  
  if (OVERWRITE_FEMA_FILES) {
    message("Running mapshaper 'simplify', 'dissolve', and 'explode' via CLI...")
    system(mapshaper_command)
    message("Writting '", output_clean_filename, "' to: \n > '", output_path, "'")
  }
  
  message()
}

# -------------------------------------------------------------------------------------
# ---- Convert cleaned FEMA geometries to geopackages ----
# -------------------------------------------------------------------------------------

# paths to FEMA 100 year flood plain files
FEMA_clean_paths      <- list.files(FEMA_CLEAN_PATH, full.names = TRUE)

for (file in FEMA_clean_paths) {
  message("Fema 100 year flood plain:\n > '", basename(file), "'")
  
  output_gpkg_filename <- gsub("_clean.geojson", "_clean.gpkg", basename(file))
  output_path          <- paste0(FEMA_GPKG_PATH, "/", output_gpkg_filename)
  
  message("Converting GEOJSON file to GPKG:\n > '", basename(file), "' > '", output_gpkg_filename, "'")
  
  # system(ogr2ogr_command)
  
  clean_gpkg_exists <- file.exists(output_path)
  
  message(" >>> '", output_gpkg_filename, "' already exists? ", clean_gpkg_exists)
  message(" >>> Overwrite? ", OVERWRITE_FEMA_FILES)
  
  ogr2ogr_command <- paste0("ogr2ogr -nlt MULTIPOLYGON ", output_path, " ", file)
  # ogr2ogr_command = paste0("ogr2ogr -nlt MULTIPOLYGON ", output_path, " ", file)
  
  if (OVERWRITE_FEMA_FILES) {
    system(ogr2ogr_command)
    message("Writting '", output_gpkg_filename, "' to: \n > '", output_path, "'")
  }
  message()
}

# # -------------------------------------------------------------------------------------
# # ---- Apply hydrofab::clean_geometries() to cleaned FEMA geometries  ----
# # -------------------------------------------------------------------------------------
# 
# paths to FEMA 100 year flood plain files
FEMA_gpkg_paths      <- list.files(FEMA_GPKG_PATH, full.names = TRUE)

for (file_path in FEMA_gpkg_paths) {
  message("Applying hydrofab::clean_geometry() to:\n > '", basename(file_path), "'")
  
  fema <-
    file_path %>%
    sf::read_sf() %>%
    sf::st_transform(5070) %>%
    sf::st_cast("POLYGON") %>%
    dplyr::mutate(
      fema_id = 1:dplyr::n()
    ) %>%
    dplyr::select(fema_id, geometry = geom)

  message(" > ", nrow(fema), " POLYGONs")
  message("Start time: ", Sys.time())

  fema_clean <- hydrofab::clean_geometry(
    catchments = fema,
    ID         = "fema_id"
    )
  
  fema_clean <-
    fema_clean %>%
    dplyr::mutate(
      source = basename(file_path),
      state  = gsub("-100yr-flood_valid_clean.gpkg", "", source)
    ) %>%
    dplyr::select(fema_id, source, state, areasqkm, geometry)

  message("End time: ", Sys.time())
  
  # geom_diff <- sf::st_difference(fema[1, ], fema_clean[1, ])
  # mapview::mapview(fema[1, ], col.regions = "red") +
  # mapview::mapview(fema_clean[1, ], col.regions = "green") +
  # mapview::mapview(geom_diff, col.regions = "white")
  
  if (OVERWRITE_FEMA_FILES) {
    message("Writting '", basename(file_path), "' to: \n > '", file_path, "'")
    sf::write_sf(
      fema_clean,
      file_path
    )
  }
  message()
  
}

# # -------------------------------------------------------------------------------------
# # ---- Partion parts of each FEMA GPKGs to the a Nextgen VPU ---- 
# # -------------------------------------------------------------------------------------

# Clean FEMA GPKG files
FEMA_CLEAN_GPKG_PATHS      <- list.files(FEMA_GPKG_PATH, full.names = TRUE)

# paths to nextgen datasets and model attribute parquet files
NEXTGEN_FILENAMES  <- list.files(nextgen_dir, full.names = FALSE)
NEXTGEN_FILE_PATHS <- paste0(nextgen_dir, NEXTGEN_FILENAMES)

for (file_path in FEMA_CLEAN_GPKG_PATHS) {
  fema_file <- basename(file_path)
  message("Partioning FEMA polygons by VPU: \n > FEMA gpkg: '", fema_file, "'")
  
  # read in fema polygons
  fema <- sf::read_sf(file_path)
  
  for (nextgen_path in NEXTGEN_FILE_PATHS) {
    nextgen_basename <- basename(nextgen_path)
    vpu              <- unlist(regmatches(nextgen_basename, gregexpr("\\d+[A-Za-z]*", nextgen_basename)))
    
    message("VPU: ", vpu)   
    message("- nextgen gpkg:\n > '", nextgen_path, "'")   
    message(" > Checking if '", fema_file, "' intersects with '", nextgen_basename, "'")
    
    # nextgen_path <- NEXTGEN_FILE_PATHS[13]
    
    # read in nextgen flowlines 
    flines <- sf::read_sf(nextgen_path, layer = "flowpaths")

    # get the FEMA polygons that intersect with the nextgen flowlines
    fema_intersect <- polygons_with_line_intersects(fema, flines)

    fema_in_nextgen <-  nrow(fema_intersect) != 0
    
    message("FEMA intersects with nextgen flowlines? ", fema_in_nextgen)

    if(fema_in_nextgen) {
      
      # create filepaths
      vpu_subfolder      <- paste0("VPU_", vpu)
      vpu_subfolder_path <- paste0(FEMA_BY_VPU_PATH, "/", vpu_subfolder)
      # vpu_subfolder_path <- FEMA_VPU_SUBFOLDERS[grepl(vpu_subfolder, FEMA_VPU_SUBFOLDERS)]
      
      fema_intersect <-
        fema_intersect %>%
        dplyr::mutate(
          vpu = vpu
        ) %>%
        dplyr::select(vpu, fema_id, source, state,
                      areasqkm, geom)
      
      # state <- gsub("-100yr-flood_valid_clean.gpkg", "", fema_file)

      fema_vpu_filename <- gsub(".gpkg", paste0("_", vpu, ".gpkg"), fema_file)
      fema_vpu_path     <- paste0(vpu_subfolder_path, "/", fema_vpu_filename)
      

      if (OVERWRITE_FEMA_FILES) {
        message("Writting '", basename(fema_vpu_filename), "' to: \n > '", fema_vpu_path, "'")
        
        sf::write_sf(
          fema_intersect,
          fema_vpu_path
        )
      }

        
    }
    message()
  }
  
  
  message(
    "--------------------------------------------------------------\n", 
    "Completed all VPU intersections for: \n > '", fema_file, "'",
    "\n--------------------------------------------------------------\n"
    )
  
}
#  
#   text = "nextgen_03W.gpkg"
#   VPU <- unlist(regmatches(text, 
#                         gregexpr("\\d+[A-Za-z]*", text)
#                         ))
# 
#   fema <-
#     file_path %>%
#     sf::read_sf() %>%
#     sf::st_transform(5070) %>%
#     sf::st_cast("POLYGON") %>%
#     dplyr::mutate(
#       fema_id = 1:dplyr::n()
#     ) %>%
#     dplyr::select(fema_id, geometry = geom)
# 
#   message(" > ", nrow(fema), " POLYGONs")
#   message("Start time: ", Sys.time())
# 
#   fema_clean <- hydrofab::clean_geometry(
#     catchments = fema,
#     ID         = "fema_id"
#     )
#   
#   fema_clean <-
#     fema_clean %>%
#     dplyr::mutate(
#       source = basename(file_path),
#       state  = gsub("-100yr-flood_valid_clean.gpkg", "", source)
#     ) %>%
#     dplyr::select(fema_id, source, state, areasqkm, geometry)
# 
#   message("End time: ", Sys.time())
#   
#   # geom_diff <- sf::st_difference(fema[1, ], fema_clean[1, ])
#   # mapview::mapview(fema[1, ], col.regions = "red") +
#   # mapview::mapview(fema_clean[1, ], col.regions = "green") +
#   # mapview::mapview(geom_diff, col.regions = "white")
#   
#   if (OVERWRITE_FEMA_FILES) {
#     message("Writting '", basename(file_path), "' to: \n > '", file_path, "'")
#     sf::write_sf(
#       fema_clean,
#       file_path
#     )
#   }
#   message()
#   
# }
# -------------------------------------------------------------------------------------
# ---- Generate bounding box gpkg for each FEMA FGB ----
# -------------------------------------------------------------------------------------

for (key in FEMA_FILENAMES) {
 
  local_fema_path <- paste0(FEMA_FGB_PATH, "/", key)
  
  gpkg_filename   <- gsub(".fgb", "_bb.gpkg", key)
  bb_save_path    <- paste0(FEMA_FGB_BB_PATH, "/", gpkg_filename)
  
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
