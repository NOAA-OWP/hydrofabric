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
# source("runners/cs_runner/config_vars.R")
# source("runners/cs_runner/config.R")
source("runners/cs_runner2/base_variables.R")
source("runners/cs_runner2/utils.R")

library(dplyr)
library(sf)
library(geos)
library(fastmap)
library(nngeo)

# ONLY runs if the main CONUS FEMA gpkg hasnt been created yet
if (!file.exists(CONUS_FEMA_GPKG_PATH)) {
  
# -------------------------------------------------------------------------------------
# ---- OVERWRITE_FEMA_FILES constant logical ----
# ---- > if TRUE, processing steps will be run again 
#          and overwrite existing previously processed files
# -------------------------------------------------------------------------------------

# Default is TRUE (i.e. a fresh processing run is done from start to finish)
OVERWRITE_FEMA_FILES  <- TRUE
DELETE_STAGING_GPKGS  <- TRUE
Sys.setenv(OGR_GEOJSON_MAX_OBJ_SIZE=0)

# -------------------------------------------------------------------------------------
# ---- Get paths to downloaded FEMA 100 FGBs ----
# -------------------------------------------------------------------------------------

FEMA_FILENAMES        <- list.files(FEMA_FGB_PATH, full.names = FALSE)
FEMA_FILE_PATHS       <- paste0(FEMA_FGB_PATH, "/", FEMA_FILENAMES)

for (file in FEMA_FILENAMES) {
  # message(file)
  
  STAGING_FILES_TO_DELETE <- c()
  
  # -------------------------------------------------------------------------------------------------------------------
  # ---- Step 1: Convert FGB to GeoJSON 
  # -------------------------------------------------------------------------------------------------------------------
 
  local_fema_path   <- paste0(FEMA_FGB_PATH, "/", file)
  geojson_filename  <- gsub(".fgb", ".geojson", file)
  geojson_save_path <- paste0(FEMA_GPKG_PATH, "/", geojson_filename)
  
  message("FEMA filename: '", file, "'")
  message("Converting \n > '", file, "' to geojson '", geojson_filename, "'")
  
  geojson_exists  <- file.exists(geojson_save_path)
  
  message(" >>> '", geojson_filename, "' already exists? ", geojson_exists)
  message(" >>> Overwrite? ", OVERWRITE_FEMA_FILES)
  
  # Step 1.1 Run FGDB to GeoJSON conversion
  ogr2ogr_command <- paste0("ogr2ogr ", geojson_save_path, " ", local_fema_path)
  
  if (OVERWRITE_FEMA_FILES || !geojson_exists) {
    system(ogr2ogr_command)
    message("Writing '", geojson_filename, "' to: \n > '", geojson_save_path, "'")
    
    STAGING_FILES_TO_DELETE <- c(STAGING_FILES_TO_DELETE, geojson_save_path)
  }
  
  # -------------------------------------------------------------------------------------------------------------------
  # ---- # Step 2: Clean GeoJSON
  # -------------------------------------------------------------------------------------------------------------------
 
  message("Simplify, dissolve, explode > '", geojson_filename, "'")
  output_clean_filename <- gsub(".geojson", "_clean.geojson", geojson_filename)
  output_clean_geojson_path     <- paste0(FEMA_GPKG_PATH, "/", output_clean_filename)
  
  clean_geojson_exists  <- file.exists(output_clean_geojson_path)
  message(" >>> '", output_clean_filename, "' already exists? ", clean_geojson_exists)
  message(" >>> Overwrite? ", OVERWRITE_FEMA_FILES)
  
  mapshaper_command = paste0('node  --max-old-space-size=16000 /opt/homebrew/bin/mapshaper ', geojson_save_path, 
                             ' -dissolve2 FLD_AR_ID \\', 
                             ' -simplify 0.1 visvalingam \\', 
                             ' -snap \\',
                             ' -o ', output_clean_geojson_path)

  # Step 2.1 Run simplify, dissolve, explode on cleaned GeoJSON
  if (OVERWRITE_FEMA_FILES || !clean_geojson_exists) {
    message("Running mapshaper 'simplify', 'dissolve', and 'explode' via CLI...")
    system(mapshaper_command)
    message("Writing '", output_clean_filename, "' to: \n > '", output_clean_geojson_path, "'")
    
    STAGING_FILES_TO_DELETE <- c(STAGING_FILES_TO_DELETE, output_clean_geojson_path)
  }
  
  # -------------------------------------------------------------------------------------------------------------------
  # ----  # Step 3: Convert cleaned GeoJSON to GeoPackage
  # -------------------------------------------------------------------------------------------------------------------
 
  message("Fema 100 year flood plain:\n > '", output_clean_filename, "'")
  output_gpkg_filename  <- gsub("_clean.geojson", "_clean.gpkg", output_clean_filename)
  output_gpkg_path      <- paste0(FEMA_GPKG_PATH, "/", output_gpkg_filename)
  
  message("Converting GEOJSON file to GPKG:\n > '", output_clean_filename, "' > '", output_gpkg_filename, "'")
  
  clean_gpkg_exists <- file.exists(output_gpkg_path)
  message(" >>> '", output_gpkg_filename, "' already exists? ", clean_gpkg_exists)
  message(" >>> Overwrite? ", OVERWRITE_FEMA_FILES)
  
  ogr2ogr_command <- paste0("ogr2ogr -nlt MULTIPOLYGON ", output_gpkg_path, " ", output_clean_geojson_path)
  
  if (OVERWRITE_FEMA_FILES || !clean_gpkg_exists) {
    system(ogr2ogr_command)
    message("Writing '", output_gpkg_filename, "' to: \n > '", output_gpkg_path, "'")
  }
  
  # -------------------------------------------------------------------------------------------------------------------
  # ---- Step 4: Apply final dissolve/snap and removal of internal boundaries in FEMA geometries  ----
  # -------------------------------------------------------------------------------------------------------------------
 
  message("Resolving internal boundaries, islands, and topology issues:\n > '", basename(output_gpkg_path), "'")
  
  fema <- sf::read_sf(output_gpkg_path)
  fema <- resolve_internal_fema_boundaries(fema, output_gpkg_path)  

  message("End time: ", Sys.time())
  
  if (OVERWRITE_FEMA_FILES) {
    message("Writting '", basename(output_gpkg_path), "' to: \n > '", output_gpkg_path, "'")
    sf::write_sf(
      fema,
      # fema_clean,
      output_gpkg_path
    )
  }
  

  
  # -------------------------------------------------------------------------------------------------------------------
  # ---- Step 5: Delete intermediary files 
  # -------------------------------------------------------------------------------------------------------------------
  
  message("Deleting intermediary files\n")
  for (delete_file in STAGING_FILES_TO_DELETE) {
    if (file.exists(delete_file)) {
      message("Deleting >>> '", delete_file, "'")
      file.remove(delete_file)
    }
    
  }

  rm(fema)

  message()
  
}

# -------------------------------------------------------------------------------------------------------------------
# ---- Apply final dissolve/snap and removal of internal boundaries in FEMA geometries  ----
# -------------------------------------------------------------------------------------------------------------------

# # paths to FEMA 100 year flood plain files
# FEMA_gpkg_paths      <- list.files(FEMA_GPKG_PATH, full.names = TRUE)

# for (file_path in FEMA_gpkg_paths) {
#   message("Resolving internal boundaries, islands, and topology issues:\n > '", basename(file_path), "'")
  
#   fema <- sf::read_sf(file_path)

#   fema <- resolve_internal_fema_boundaries(fema) 
  
#   # fema <-
#   #   fema[!sf::st_is_empty(fema), ] %>% 
#   #   sf::st_transform(5070)
  
#   # fema <-
#   #   fema %>% 
#   #   # fema[!sf::st_is_empty(fema), ] %>% 
#   #   dplyr::select(geometry = geom) %>%
#   #   add_predicate_group_id(sf::st_intersects) %>% 
#   #   sf::st_make_valid() %>% 
#   #   dplyr::group_by(group_id) %>% 
#   #   dplyr::summarise(
#   #     geometry = sf::st_combine(sf::st_union(geometry))
#   #   ) %>% 
#   #   dplyr::ungroup() %>% 
#   #   dplyr::select(-group_id) %>% 
#   #   add_predicate_group_id(sf::st_intersects) %>% 
#   #   rmapshaper::ms_dissolve(sys = TRUE, sys_mem = 16) %>% 
#   #   rmapshaper::ms_explode(sys = TRUE, sys_mem = 16) %>% 
#   #   dplyr::mutate(
#   #     fema_id = as.character(1:dplyr::n())
#   #   ) %>% 
#   #   dplyr::select(fema_id, geometry)
  
#   # # mapview::mapview(fema, color = 'cyan', col.regions = "cyan") + 
#   # # mapview::mapview(end_fema, color = 'red', col.regions = "white") 
  
#   # fema <- 
#   #   fema %>% 
#   #   dplyr::mutate(
#   #     source = basename(file_path),
#   #     state  = gsub("-100yr-flood_valid_clean.gpkg", "", source)
#   #   ) %>%
#   #   dplyr::select(fema_id, source, state, 
#   #                 # areasqkm, 
#   #                 geometry)
  
#   message("End time: ", Sys.time())
  
#   if (OVERWRITE_FEMA_FILES) {
#     message("Writting '", basename(file_path), "' to: \n > '", file_path, "'")
#     sf::write_sf(
#       # fema_clean,
#       fema,
#       file_path
#     )
#   }
#   message()
  
# }

# -------------------------------------------------------------------------------------
# ---- Partion parts of each FEMA GPKGs to a Nextgen VPU ---- 
# -------------------------------------------------------------------------------------

# Clean FEMA GPKG files
FEMA_CLEAN_GPKG_PATHS      <- list.files(FEMA_GPKG_PATH, full.names = TRUE)

# paths to nextgen datasets and model attribute parquet files
# CONUS_NEXTGEN_GPKG_PATH 

# layer_info = sf::st_layers(CONUS_NEXTGEN_GPKG_PATH)
# layer_info[layer_info$name == "flowpaths", "fields"]
# layer_inf

# flines = sf::read_sf(CONUS_NEXTGEN_GPKG_PATH, layer = "flowpaths")
# flines %>% names()
# (flines$vpuid  %>% unique()) %in% VPU_IDS
# VPU_IDS
# VPU_IDS %in% (flines$vpuid  %>% unique())

# query the unique VPU_IDS from the flowlines layer from sf::read_sf(CONUS_NEXTGEN_GPKG_PATH, layer = "flowpaths")
CONUS_VPU_IDS <- 
            CONUS_NEXTGEN_GPKG_PATH  %>% 
            sf::read_sf(query = "SELECT DISTINCT vpuid FROM flowpaths") %>%
            dplyr::pull()


# NEXTGEN_FILENAMES    <- list.files(NEXTGEN_DIR, full.names = FALSE)
# NEXTGEN_FILE_PATHS   <- paste0(NEXTGEN_DIR, NEXTGEN_FILENAMES)


for (file_path in FEMA_CLEAN_GPKG_PATHS) {
   # i = 35
  # file_path = FEMA_CLEAN_GPKG_PATHS[i]

  fema_file <- basename(file_path)
  
  message("Partioning FEMA polygons by VPU: \n > FEMA gpkg: '", fema_file, "'")
  
  # read in fema polygons
  fema <- sf::read_sf(file_path)
  
  for (vpu in CONUS_VPU_IDS) {
    # j = 8 
    # vpu = CONUS_VPU_IDS[j]
    
    # nextgen_basename <- basename(nextgen_path)
    # vpu              <- unlist(regmatches(nextgen_basename, gregexpr("\\d+[A-Za-z]*", nextgen_basename)))
    
    message("VPU: ", vpu)   
    # message("- nextgen gpkg:\n > '", nextgen_path, "'")   
    message(" > Checking if '", fema_file, "' intersects with CONUS flowpaths in VPU '", vpu, "'")
    
    # read in nextgen flowlines 
    flines <- 
            CONUS_NEXTGEN_GPKG_PATH  %>% 
            sf::read_sf(query = paste0("SELECT * FROM flowpaths WHERE vpuid = '", vpu, "'"))

    # flines <- sf::read_sf(nextgen_path, layer = "flowpaths")
    
    # get the FEMA polygons that intersect with the nextgen flowlines
    fema_intersect <- polygons_with_line_intersects(fema, flines)
    
    fema_in_nextgen <-  nrow(fema_intersect) != 0
    
    message("FEMA intersects with nextgen flowlines? ", fema_in_nextgen)
    
    if(fema_in_nextgen) {
      
      # create filepaths
      vpu_subfolder      <- paste0("vpu-", vpu)
      # vpu_subfolder_path <- paste0(FEMA_BY_VPU_PATH, "/", vpu_subfolder, "/states")
      vpu_subfolder_path <- paste0(FEMA_BY_VPU_PATH, "/", vpu_subfolder)
      
      # vpu_subfolder_path <- FEMA_VPU_SUBFOLDERS[grepl(vpu_subfolder, FEMA_VPU_SUBFOLDERS)]
      
      fema_intersect <-
        fema_intersect %>%
        dplyr::mutate(
          vpu = vpu
        ) %>%
        dplyr::select(vpu, fema_id, source, state, geom)
      
      # state <- gsub("-100yr-flood_valid_clean.gpkg", "", fema_file)
      
      # BASE_DIRS_LIST$fema_by_vpu_subsets_dirs
      
      fema_vpu_filename <- gsub(".gpkg", paste0("_", vpu, ".gpkg"), fema_file)
      fema_vpu_path     <- paste0(vpu_subfolder_path, "/subsets/", fema_vpu_filename)
      
      
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

# for (i in seq_along(BASE_DIRS_LIST$fema_by_vpu_subdirs)) {
#   vpu_subdir <- BASE_DIRS_LIST$fema_by_vpu_subdirs[i]
#   
#   message(i, " - ", vpu_subdir)
#   
#   files_to_move = list.files(vpu_subdir, pattern = ".gpkg", full.names = T)
#   
#   for (file in files_to_move) {
#     new_path <- paste0(vpu_subdir, "/subsets/",     basename(file))
#     message(" > ", file, "\n > ", new_path, "\n")
#     fs::file_move(
#       file, 
#       new_path
#     )
#   }
# }

# -------------------------------------------------------------------------------------
# ---- Loop through each VPU subfolder and merge all of the Geopackages into one---- 
# -------------------------------------------------------------------------------------
# BASE_DIRS_LIST$fema_by_vpu_subdirs
# BASE_DIRS_LIST$fema_by_vpu_subsets_dirs
# for (vpu_dir in FEMA_VPU_SUBFOLDERS) {
BASE_DIRS_LIST$fema_by_vpu_subdirs

for (i in seq_along(BASE_DIRS_LIST$fema_by_vpu_subdirs)) {
  # for (i in 1:4) {
  # i = 8
  
  vpu_dir <-  BASE_DIRS_LIST$fema_by_vpu_subdirs[i]
  # vpu_dir <-  FEMA_VPU_SUBFOLDERS[i]
  
  vpu_subset_dir <- BASE_DIRS_LIST$fema_by_vpu_subsets_dirs[i]
  
  message("Merging files in '", basename(vpu_dir), "' directory...")
  vpu_subdirs    <- list.files(vpu_subset_dir, full.names = TRUE)
  # vpu_subdirs    <- list.files(vpu_dir, full.names = TRUE)
  
  # message(paste0("\n > ", basename(vpu_subdirs), collapse = ""))
  
  # path to the merged directory where the final merged geopackge will end up
  master_name       <- paste0("fema-", gsub("VPU", "vpu", basename(vpu_dir)))
  master_gpkg_name  <- paste0(master_name, ".gpkg")
  master_filepath   <- paste0(vpu_dir, "/merged/", master_gpkg_name)
  
  # if the file already exists, remove it so we dont OVER append data to the "master file"
  if (file.exists(master_filepath)) {
    file.remove(master_filepath)
  }
  
  # fema state geopackages partioned for the specific VPU
  fema_state_gpkgs        <- list.files(vpu_subset_dir, full.names = TRUE)
  # master_output_filepath  <- paste0(vpu_dir, "/", gsub(".gpkg", "_output.gpkg", master_gpkg_name))
  # 
  # # make sure to ignore the master file if it already exists
  # fema_state_gpkgs <- fema_state_gpkgs[fema_state_gpkgs != master_filepath & fema_state_gpkgs != master_output_filepath]
  
  for(gpkg_file in fema_state_gpkgs) {
    # message(" - Appending '", basename(gpkg_file), "' to master FEMA VPU gpkg:\n  >  '",
    #         basename(gpkg_file), " > ", basename(master_filepath),
    #         "'")
    message(" > '",
            basename(gpkg_file), " > ", basename(master_filepath),
            "'")

    ogr2ogr_merge_command <- paste0("ogr2ogr -f 'gpkg' -append -nln ", master_name, " ",
                                    master_filepath,
                                    " ", gpkg_file
    )

    # if (OVERWRITE_FEMA_FILES) {
      system(ogr2ogr_merge_command)
    # }
      
  }
  
  # has_fema_state_gpkgs <- length(fema_state_gpkgs) > 0
  # 
  # if(DELETE_STAGING_GPKGS && has_fema_state_gpkgs) {
  #   message(" - Deleting individual gpkgs from '", vpu_dir, "' directory...")
  #   remove_gpkg_cmds <- paste0("rm ", fema_state_gpkgs)
  #   
  #   for (remove_cmd in remove_gpkg_cmds) {
  #     message("  >  '", remove_cmd, "'")
  #     system(remove_cmd)
  #   }
  # }
  
  # message()
  message("Merge complete!")
  message("Merged '", basename(vpu_dir), "' FEMA output geopackage:\n --> '", master_filepath, "'")
  message()
  
}
# -------------------------------------------------------------------------------------
# ----Apply simplify, dissolve, explode on the MERGED polygons  ---- 
# -------------------------------------------------------------------------------------

# # NOTE: remove past runs for testing...
# for (i in list.files(FEMA_VPU_SUBFOLDERS, full.names = T)[grepl("_output.gpkg", list.files(FEMA_VPU_SUBFOLDERS, full.names = T))]) {
#   file.remove(i)
# }
# seq_along(BASE_DIRS_LIST$fema_by_vpu_subdirs)
# rm(STAGING_FILES_TO_DELETE, i, VPU, vpu_dir)

for (i in 1:length(FEMA_VPU_SUBFOLDERS)) {

  # i = 8
  # FEMA_VPU_SUBFOLDERS
  STAGING_FILES_TO_DELETE <- c()
  
  vpu_dir    <- FEMA_VPU_SUBFOLDERS[i]
  VPU        <- basename(vpu_dir)
  
  message(i, " - Attempting to union FEMA polygons for '", VPU, "'...")
  
  # path to the merged directory where the final merged geopackage will end up
  master_name       <- paste0("fema-", gsub("VPU", "vpu", basename(vpu_dir)))
  master_gpkg_name  <- paste0(master_name, ".gpkg")
  master_filepath   <- paste0(vpu_dir, "/merged/", master_gpkg_name)
  
  master_geojson_name       <- paste0(master_name, ".geojson")
  master_geojson_filepath   <- paste0(vpu_dir, "/merged/", master_geojson_name)
  
  # updated_gpkg_name  <- gsub(".gpkg", "-output.gpkg", master_gpkg_name)
  updated_gpkg_name  <- master_gpkg_name
  updated_filepath   <- paste0(vpu_dir, "/output/", master_gpkg_name)
  
  updated_gpkg_exists  <- file.exists(updated_filepath)
  
  # remove output file if it already exists
  if (updated_gpkg_exists) {
    file.remove(updated_filepath)
  }
  
  message("VPU Merged FEMA filename: '", master_gpkg_name, "'")
  message("> Simplifying, dissolve, exploding VPU aggregated FEMA polygons... '", basename(master_filepath), "'")
  
  if(!file.exists(master_filepath)) { 
    message("No FEMA geometries in '", VPU, "'")
    message()
    next
  }
  
  message("Converting \n > '", basename(master_filepath), "' to geojson '", master_geojson_name, "'")
  
  geojson_exists  <- file.exists(master_geojson_filepath)
  
  # message(" >>> '", geojson_filename, "' already exists? ", geojson_exists)
  # message(" >>> Overwrite? ", OVERWRITE_FEMA_FILES)
  gpkg_to_geojson_cmd <- paste0("ogr2ogr -s_srs EPSG:5070 -t_srs EPSG:5070 ", master_geojson_filepath, " ", master_filepath)
  # gpkg_to_geojson_cmd <- paste0("ogr2ogr -f GEOJSON -s_srs EPSG:5070 -t_srs EPSG:5070 ", master_geojson_filepath, " ", master_filepath)
  # gpkg_to_geojson_cmd <- paste0("ogr2ogr ", master_geojson_filepath, " ", master_filepath)
  # file.remove(master_geojson_filepath)
  
  # if (OVERWRITE_FEMA_FILES || !geojson_exists) {
  system(gpkg_to_geojson_cmd)
  message("Writing '", master_geojson_name, "' to: \n > '", master_geojson_filepath, "'")
  
  STAGING_FILES_TO_DELETE <- c(STAGING_FILES_TO_DELETE, master_geojson_filepath)
  # }
  
  # master_gj <- sf::read_sf(master_geojson_filepath)
  # master_gpkg <- sf::read_sf(master_filepath)

  # Clean GeoJSON
  message("Simplify, dissolve, explode > '", master_geojson_name, "'")
  output_clean_filename      <- gsub(".geojson", "_clean.geojson", master_geojson_name)
  output_clean_geojson_path  <- paste0(vpu_dir, "/merged/", output_clean_filename)
  
  clean_geojson_exists  <- file.exists(output_clean_geojson_path)
  message(" >>> '", output_clean_filename, "' already exists? ", clean_geojson_exists)
  # message(" >>> Overwrite? ", OVERWRITE_FEMA_FILES)

  if (clean_geojson_exists) {
    file.remove(output_clean_geojson_path)
  }
  
  mapshaper_command = paste0('node  --max-old-space-size=16000 /opt/homebrew/bin/mapshaper ',
                             master_geojson_filepath, 
                             # ' -clean \\', 
                             # ' -explode \\',
                             # ' -dissolve2 \\',
                             ' -simplify 0.3 visvalingam \\', 
                             ' -snap \\',
                             ' -explode \\',
                             ' -clean \\',
                             # ' -proj EPSG:5070 \\',
                             ' -o ', output_clean_geojson_path
  )
  
  system(mapshaper_command)
  # message("Writing '", master_geojson_name, "' to: \n > '", master_geojson_filepath, "'")
  
  STAGING_FILES_TO_DELETE <- c(STAGING_FILES_TO_DELETE, output_clean_geojson_path)
  # output_clean_gpkg_filename   <- gsub(".geojson", ".gpkg", master_geojson_name)
  # output_clean_gpkg_path       <- paste0(vpu_dir, "/merged/", output_clean_gpkg_filename)
  
  # fema_vpu <- sf::read_sf(master_filepath)
  # geojson_to_gpkg_cmd <- paste0("ogr2ogr -f GPKG ", updated_filepath, " ", output_clean_geojson_path)
  geojson_to_gpkg_cmd <- paste0("ogr2ogr -nlt MULTIPOLYGON -s_srs EPSG:5070 -t_srs EPSG:5070  ", 
                                updated_filepath, " ", 
                                output_clean_geojson_path)
  # geojson_to_gpkg_cmd <- paste0("ogr2ogr ", updated_filepath, " ", output_clean_geojson_path)

  # updated_gpkg_exists  <- file.exists(updated_filepath)
  # if (updated_gpkg_exists) {
  #   file.remove(updated_filepath)
  # }
  
  
  # if (OVERWRITE_FEMA_FILES || !updated_gpkg_exists) {
  system(geojson_to_gpkg_cmd)
  message("Writing '", updated_gpkg_name, "' to: \n > '", updated_filepath, "'")
  # }

  # sf::st_layers(updated_filepath)
  
  # mapview::npts(fema)
  fema  <- 
    sf::read_sf(updated_filepath) %>%
    # sf::read_sf(output_clean_geojson_path) %>%
    # rmapshaper::ms_explode(sys=TRUE, sys_mem = 16) %>% 
    dplyr::mutate(
      vpu      = gsub("vpu-", "", VPU),
      fema_id = as.character(1:dplyr::n())
    ) %>% 
    dplyr::select(
      vpu,
      fema_id, 
      geom
    )
  
  # fema %>% 
  #   rmapshaper::ms_simplify(keep = 0.5, keep_shapes = T) %>% 
  #   dplyr::group_by(fema_id) %>% 
  #   dplyr::mutate(pts = mapview::npts(geom)) %>% 
  #   dplyr::arrange(-pts) 
  
  # remove before writting updated version
  file.remove(updated_filepath)
  
  sf::write_sf(
    fema, 
    updated_filepath, 
    append = FALSE
  )
  
  message("Deleting intermediary files\n")
  for (delete_file in STAGING_FILES_TO_DELETE) {
    if (file.exists(delete_file)) {
      message("Deleting >>> '", delete_file, "'")
      file.remove(delete_file)
    }
    
  }
}


# -------------------------------------------------------------------------------------
# ---- Store all FEMA layers in a single conus_fema.gpkg 
# -------------------------------------------------------------------------------------

FEMA_VPU_SUBFOLDERS

all_fema_vpu_layers <- list.files(BASE_DIRS_LIST$fema_by_vpu_output_dirs, full.names = TRUE)

# fema_vpu_layers <- list.files(FEMA_VPU_SUBFOLDERS, full.names = T)[grepl("_output.gpkg", list.files(FEMA_VPU_SUBFOLDERS))]

combine_gpkg_files(all_fema_vpu_layers, CONUS_FEMA_GPKG_PATH)
  

# -------------------------------------------------------------------------------------

} else {
  message("'", basename(CONUS_FEMA_GPKG_PATH), "'already exists at \n> ", CONUS_FEMA_GPKG_PATH)
}

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

# # -------------------------------------------------------------------------------------
# # ---- Union each VPU geopackage (either on state or just touching predicate) ---- 
# # -------------------------------------------------------------------------------------
# 
# for (i in 1:length(FEMA_VPU_SUBFOLDERS)) {
#   
#   vpu_dir    <- FEMA_VPU_SUBFOLDERS[i]
#   VPU        <- basename(vpu_dir)
#   
#   message(i, " - Attempting to union FEMA polygons for '", VPU, "'...")
#   
#   # path to the merged directory where the final merged geopackage will end up
#   master_name       <- paste0("fema_", gsub("VPU", "vpu", basename(vpu_dir)))
#   master_gpkg_name  <- paste0(master_name, ".gpkg")
#   master_filepath   <- paste0(vpu_dir, "/", master_gpkg_name)
#   
#   updated_gpkg_name  <- gsub(".gpkg", "_output.gpkg", master_gpkg_name)
#   updated_filepath   <- paste0(vpu_dir, "/", updated_gpkg_name)
#   
#   message("> Re-unioning and re-exploding geometries in '", basename(master_filepath), "'")
#   
#   if(!file.exists(master_filepath)) { 
#     message("No FEMA geometries in '", VPU, "'")
#     message()
#     next
#   }
#   
#   
#   fema_vpu <- sf::read_sf(master_filepath)
#   
#   # fema_vpu %>% sf::st_geometry_type() %>% unique()
#   # fema_vpu %>% mapview::npts()
#   # fema_vpu %>% sf::st_is_valid() %>% all()
#   # fema_vpu %>% 
#   #   sf::st_make_valid() %>% 
#   #   sf::st_geometry_type() %>% 
#   #   unique()
#   
#   geom_type_counts <- table(sf::st_geometry_type(fema_vpu))
#   
#   message("Geometry counts before casting all geometries to MULTIPOLYGON:")
#   for (g in seq_along(geom_type_counts)) {
#     message(" > ", names(geom_type_counts[g]), ": ", geom_type_counts[g])
#   }
#   
#   # mapview::mapview(fema_vpu, color = 'red', col.regions = 'white') +
#   #       mapview::mapview(fema_union, color = 'green', col.regions = 'white')
#   
#   # fema_vpu %>% 
#   #   sf::st_make_valid() %>% 
#   #   dplyr::filter(sf::st_geometry_type(geom) %in% c("POLYGON", "MULTIPOLYGON")) %>% 
#   #   sf::st_is_valid() %>% 
#   #   all()
#   
#   tryCatch({
#     
#     fema_vpu <- 
#       fema_vpu %>% 
#       nngeo::st_remove_holes(max_area = 200) %>% 
#       # sf::st_make_valid() %>%
#       # dplyr::filter(sf::st_geometry_type(geom) %in% c("POLYGON", "MULTIPOLYGON")) %>% 
#       add_predicate_group_id(sf::st_intersects) %>% 
#       dplyr::group_by(group_id) %>%
#       rmapshaper::ms_dissolve(sys = TRUE, sys_mem = 16) %>% 
#       rmapshaper::ms_explode(sys = TRUE, sys_mem = 16) %>% 
#       dplyr::ungroup() %>% 
#       nngeo::st_remove_holes(max_area = 200) %>% 
#       dplyr::mutate(      
#         vpu      = gsub("VPU_", "", VPU),
#         fema_id  = as.character(1:dplyr::n())
#       ) %>% 
#       dplyr::select(
#         vpu, fema_id,
#         # state, 
#         geom = geometry
#       )
#     
#   }, error = function(e) {
#     message(VPU, " threw into the following error \n ", e)
#     message(" > Cleaning ", VPU, " using a backup cleaning strategy...")
#     
#     fema_vpu <- 
#       fema_vpu %>% 
#       sf::st_make_valid() %>% 
#       dplyr::mutate(      
#         vpu      = gsub("VPU_", "", VPU),
#         fema_id  = as.character(1:dplyr::n())
#       ) %>% 
#       dplyr::select(
#         vpu, fema_id,
#         # state, 
#         geom
#       )
#   
#     })
#   
#   if (OVERWRITE_FEMA_FILES) {
#     
#     message("> Overwritting '", basename(master_filepath), "' with final clean version...")
#     
#     # union_file_path <- gsub(".gpkg", "_union.gpkg", fema_vpu_file)
#     # message("> writting '", basename(union_file_path), "' (unioned and exploded version)")
#     
#     sf::write_sf(
#       fema_vpu,
#       updated_filepath
#       # master_filepath
#       # union_file_path
#     )
#   }
#   message()
# }
