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
library(fastmap)
library(nngeo)

# devtools::install_github("anguswg-ucsb/hydrofabric3D")

# TODO: Steps that converts FGB to geojson and then geojson to gpkg can be put into a single loop
# TODO: Delete old files as needed

# -------------------------------------------------------------------------------------
# ---- OVERWRITE_FEMA_FILES constant logical ----
# ---- > if TRUE, processing steps will be run again 
#          and overwrite existing previously processed files
# -------------------------------------------------------------------------------------

# Default is TRUE (i.e. a fresh processing run is done from start to finish)
OVERWRITE_FEMA_FILES  <- TRUE
DELETE_STAGING_GPKGS  <- FALSE
# DELETE_STAGING_GPKGS  <- TRUE

# -------------------------------------------------------------------------------------
# ---- Create directories (if they do NOT exist) ----
# -------------------------------------------------------------------------------------

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
  # message(VPU_SUBFOLDER)
  
  # state_dir  = paste0(VPU_SUBFOLDER, "/states/")
  # merged_dir = paste0(VPU_SUBFOLDER, "/merged/")
  
  if (!dir.exists(VPU_SUBFOLDER)) {
    message("Creating FEMA VPU subfolder...")
    message(paste0("'/", basename(VPU_SUBFOLDER), "' directory does not exist...\n  Creating directory:\n > '", VPU_SUBFOLDER, "'"))
    dir.create(VPU_SUBFOLDER)
  }
  # if (!dir.exists(state_dir)) { 
  #   message("Creating FEMA VPU states subfolder...")
  #   message(paste0("'/", basename(state_dir), "' directory does not exist...\n  Creating directory:\n > '", state_dir, "'"))
  #   dir.create(state_dir)
  # }
  # if (!dir.exists(merged_dir)) { 
  #   message("Creating FEMA VPU merged subfolder...")
  #   message(paste0("'/", basename(merged_dir), "' directory does not exist...\n  Creating directory:\n > '", merged_dir, "'"))
  #   dir.create(merged_dir)
  # }
}

# -------------------------------------------------------------------------------------
# ---- Get paths to downloaded FEMA 100 FGBs ----
# -------------------------------------------------------------------------------------

FEMA_FILENAMES        <- list.files(FEMA_FGB_PATH, full.names = FALSE)
FEMA_FILE_PATHS       <- paste0(FEMA_FGB_PATH, "/", FEMA_FILENAMES)

for (file in FEMA_FILENAMES) {
  
  STAGING_FILES_TO_DELETE <- c()
  
  # Convert FGB to GeoJSON
  local_fema_path   <- paste0(FEMA_FGB_PATH, "/", file)
  geojson_filename  <- gsub(".fgb", ".geojson", file)
  geojson_save_path <- paste0(FEMA_GPKG_PATH, "/", geojson_filename)
  
  message("FEMA filename: '", file, "'")
  message("Converting \n > '", file, "' to geojson '", geojson_filename, "'")
  
  geojson_exists  <- file.exists(geojson_save_path)
  
  message(" >>> '", geojson_filename, "' already exists? ", geojson_exists)
  message(" >>> Overwrite? ", OVERWRITE_FEMA_FILES)
  
  ogr2ogr_command <- paste0("ogr2ogr ", geojson_save_path, " ", local_fema_path)
  
  if (OVERWRITE_FEMA_FILES || !geojson_exists) {
    system(ogr2ogr_command)
    message("Writing '", geojson_filename, "' to: \n > '", geojson_save_path, "'")
    
    STAGING_FILES_TO_DELETE <- c(STAGING_FILES_TO_DELETE, geojson_save_path)
  }
  
  # Clean GeoJSON
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
                             ' -o ', output_clean_geojson_path
  )
  
  
  if (OVERWRITE_FEMA_FILES || !clean_geojson_exists) {
    message("Running mapshaper 'simplify', 'dissolve', and 'explode' via CLI...")
    system(mapshaper_command)
    message("Writing '", output_clean_filename, "' to: \n > '", output_clean_geojson_path, "'")
    
    STAGING_FILES_TO_DELETE <- c(STAGING_FILES_TO_DELETE, output_clean_geojson_path)
  }
  
  # Convert cleaned GeoJSON to GeoPackage
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
  
  message("Deleting intermediary files\n")
  for (delete_file in STAGING_FILES_TO_DELETE) {
    if (file.exists(delete_file)) {
      message("Deleting >>> '", delete_file, "'")
      file.remove(delete_file)
    }
    
  }
  
  message()
  
}

# -------------------------------------------------------------------------------------------------------------------
# ---- Apply final dissolve/snap and removal of internal boundaries in FEMA geometries  ----
# -------------------------------------------------------------------------------------------------------------------

# paths to FEMA 100 year flood plain files
FEMA_gpkg_paths      <- list.files(FEMA_GPKG_PATH, full.names = TRUE)

for (file_path in FEMA_gpkg_paths) {
  message("Resolving internal boundaries, islands, and topology issues:\n > '", basename(file_path), "'")
  
  fema <- sf::read_sf(file_path)
  
  fema <-
    fema[!sf::st_is_empty(fema), ] %>% 
    sf::st_transform(5070)
  
  #  TODO: Snap using geos::geos_snap()
  # fema <-
  #   geos::geos_snap(
  #     geos::as_geos_geometry(fema),
  #     geos::as_geos_geometry(fema),
  #     tolerance = 1
  #     ) %>%
  #   geos::geos_make_valid()  %>%
  #   sf::st_as_sf()
  
  # TODO: we get this error when trying to use the geometry column after geos snapping
  # TODO: Error = "Error: Not compatible with STRSXP: [type=NULL]."
  # fema %>%
  # sf::st_cast("POLYGON")
  
  # TODO: Snap using sf::st_snap()
  # fema <- sf::st_snap(
  #             fema,
  #             fema,
  #             tolerance = 2
  #             )
  
  fema <-
    fema %>% 
    # fema[!sf::st_is_empty(fema), ] %>% 
    dplyr::select(geometry = geom) %>%
    add_predicate_group_id(sf::st_intersects) %>% 
    sf::st_make_valid() %>% 
    dplyr::group_by(group_id) %>% 
    dplyr::summarise(
      geometry = sf::st_combine(sf::st_union(geometry))
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-group_id) %>% 
    add_predicate_group_id(sf::st_intersects) %>% 
    rmapshaper::ms_dissolve(sys = TRUE, sys_mem = 16) %>% 
    rmapshaper::ms_explode(sys = TRUE, sys_mem = 16) %>% 
    dplyr::mutate(
      fema_id = as.character(1:dplyr::n())
    ) %>% 
    dplyr::select(fema_id, geometry)
  
  # mapview::mapview(fema, color = 'cyan', col.regions = "cyan") + 
  # mapview::mapview(end_fema, color = 'red', col.regions = "white") 
  
  fema <- 
    fema %>% 
    dplyr::mutate(
      source = basename(file_path),
      state  = gsub("-100yr-flood_valid_clean.gpkg", "", source)
    ) %>%
    dplyr::select(fema_id, source, state, 
                  # areasqkm, 
                  geometry)
  
  message("End time: ", Sys.time())
  
  if (OVERWRITE_FEMA_FILES) {
    message("Writting '", basename(file_path), "' to: \n > '", file_path, "'")
    sf::write_sf(
      # fema_clean,
      fema,
      file_path
    )
  }
  message()
  
}

# -------------------------------------------------------------------------------------
# ---- Partion parts of each FEMA GPKGs to a Nextgen VPU ---- 
# -------------------------------------------------------------------------------------

# Clean FEMA GPKG files
FEMA_CLEAN_GPKG_PATHS      <- list.files(FEMA_GPKG_PATH, full.names = TRUE)

# paths to nextgen datasets and model attribute parquet files
NEXTGEN_FILENAMES    <- list.files(nextgen_dir, full.names = FALSE)
NEXTGEN_FILE_PATHS   <- paste0(nextgen_dir, NEXTGEN_FILENAMES)

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
    
    # read in nextgen flowlines 
    flines <- sf::read_sf(nextgen_path, layer = "flowpaths")
    
    # get the FEMA polygons that intersect with the nextgen flowlines
    fema_intersect <- polygons_with_line_intersects(fema, flines)
    
    fema_in_nextgen <-  nrow(fema_intersect) != 0
    
    message("FEMA intersects with nextgen flowlines? ", fema_in_nextgen)
    
    if(fema_in_nextgen) {
      
      # create filepaths
      vpu_subfolder      <- paste0("VPU_", vpu)
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

# -------------------------------------------------------------------------------------
# ---- Loop through each VPU subfolder and merge all of the Geopackages into one---- 
# -------------------------------------------------------------------------------------

for (vpu_dir in FEMA_VPU_SUBFOLDERS) {
  # for (i in 1:4) {
  # vpu_dir = FEMA_VPU_SUBFOLDERS2[12]
  message("Merging files in '", basename(vpu_dir), "' directory...")
  # }
  
  vpu_subdirs    <- list.files(vpu_dir, full.names = TRUE)
  
  # path to the merged directory where the final merged geopackge will end up
  master_name       <- paste0("fema_", gsub("VPU", "vpu", basename(vpu_dir)))
  master_gpkg_name  <- paste0(master_name, ".gpkg")
  master_filepath   <- paste0(vpu_dir, "/", master_gpkg_name)
  
  # fema state geopackages partioned for the specific VPU
  fema_state_gpkgs <- list.files(vpu_dir, full.names = TRUE)
  
  # make sure to ignore the master file if it already exists
  fema_state_gpkgs <- fema_state_gpkgs[fema_state_gpkgs != master_filepath]
  
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
    
    if (OVERWRITE_FEMA_FILES) {
      system(ogr2ogr_merge_command)
    }
  }
  
  has_fema_state_gpkgs <- length(fema_state_gpkgs) > 0
  
  if(DELETE_STAGING_GPKGS && has_fema_state_gpkgs) {
    message(" - Deleting individual gpkgs from '", vpu_dir, "' directory...")
    # message("- Deleting individual gpkgs from 'states' directory:\n > '", states_dir, "'")
    
    remove_gpkg_cmds <- paste0("rm ", fema_state_gpkgs)
    
    for (remove_cmd in remove_gpkg_cmds) {
      message("  >  '", remove_cmd, "'")
      system(remove_cmd)
    }
  }
  
  # message()
  message("Merge complete!")
  message("Merged '", basename(vpu_dir), "' FEMA output geopackage:\n --> '", master_filepath, "'")
  message()
}

# -------------------------------------------------------------------------------------
# ---- Union each VPU geopackage (either on state or just touching predicate) ---- 
# -------------------------------------------------------------------------------------
# for (i in 5:length(FEMA_VPU_SUBFOLDERS)) {
#   
#   vpu_dir    <- FEMA_VPU_SUBFOLDERS[i]
#   VPU        <- basename(vpu_dir)
#   
#   message(i, " - Attempting to union FEMA polygons for '", VPU, "'...")
# }

for (i in 5:length(FEMA_VPU_SUBFOLDERS)) {
  
  vpu_dir    <- FEMA_VPU_SUBFOLDERS[i]
  VPU        <- basename(vpu_dir)
  
  message(i, " - Attempting to union FEMA polygons for '", VPU, "'...")
  
  # path to the merged directory where the final merged geopackage will end up
  master_name       <- paste0("fema_", gsub("VPU", "vpu", basename(vpu_dir)))
  master_gpkg_name  <- paste0(master_name, ".gpkg")
  master_filepath   <- paste0(vpu_dir, "/", master_gpkg_name)
  
  updated_gpkg_name  <- gsub(".gpkg", "_output.gpkg", master_gpkg_name)
  updated_filepath   <- paste0(vpu_dir, "/", updated_gpkg_name)
  
  message("> Re-unioning and re-exploding geometries in '", basename(master_filepath), "'")
  
  if(!file.exists(master_filepath)) { 
    message("No FEMA geometries in '", VPU, "'")
    message()
    next
  }
  
  
  fema_vpu <- sf::read_sf(master_filepath)
  
  # fema_vpu %>% sf::st_geometry_type() %>% unique()
  # fema_vpu %>% mapview::npts()
  # fema_vpu %>% sf::st_is_valid() %>% all()
  # fema_vpu %>% 
  #   sf::st_make_valid() %>% 
  #   sf::st_geometry_type() %>% 
  #   unique()
  geom_type_counts <- table(sf::st_geometry_type(fema_vpu))
  
  message("Geometry counts before casting all geometries to MULTIPOLYGON:")
  for (g in seq_along(geom_type_counts)) {
    message(" > ", names(geom_type_counts[g]), ": ", geom_type_counts[g])
  }
  
  # mapview::mapview(fema_vpu, color = 'red', col.regions = 'white') +
  #       mapview::mapview(fema_union, color = 'green', col.regions = 'white')
  
  # fema_vpu %>% 
  #   sf::st_make_valid() %>% 
  #   dplyr::filter(sf::st_geometry_type(geom) %in% c("POLYGON", "MULTIPOLYGON")) %>% 
  #   sf::st_is_valid() %>% 
  #   all()
  
  tryCatch({
    
    fema_vpu <- 
      fema_vpu %>% 
      nngeo::st_remove_holes(max_area = 200) %>% 
      sf::st_make_valid() %>%
      # dplyr::filter(sf::st_geometry_type(geom) %in% c("POLYGON", "MULTIPOLYGON")) %>% 
      add_predicate_group_id(sf::st_intersects) %>% 
      dplyr::group_by(group_id) %>% 
      rmapshaper::ms_dissolve(sys = TRUE, sys_mem = 16) %>% 
      rmapshaper::ms_explode(sys = TRUE, sys_mem = 16) %>% 
      dplyr::ungroup() %>% 
      nngeo::st_remove_holes(max_area = 200) %>% 
      dplyr::mutate(      
        vpu      = gsub("VPU_", "", VPU),
        fema_id  = as.character(1:dplyr::n())
      ) %>% 
      dplyr::select(
        vpu, fema_id,
        # state, 
        geom = geometry
      )
    
  }, error = function(e) {
    message(VPU, " threw into the following error \n ", e)
    message(" > Cleaning ", VPU, " using a backup cleaning strategy...")
    
    fema_vpu <- 
      fema_vpu %>% 
      sf::st_make_valid() %>% 
      dplyr::mutate(      
        vpu      = gsub("VPU_", "", VPU),
        fema_id  = as.character(1:dplyr::n())
      ) %>% 
      dplyr::select(
        vpu, fema_id,
        # state, 
        geom
      )
  
    })
  
    # fema_vpu2 <- 
    #   fema_vpu %>% 
    #   nngeo::st_remove_holes(max_area = 200) %>% 
    #   # sf::st_make_valid() %>% 
    #   # dplyr::filter(sf::st_geometry_type(geom) %in% c("POLYGON", "MULTIPOLYGON")) %>% 
    #   add_predicate_group_id(sf::st_intersects) %>% 
    #   dplyr::group_by(group_id) %>% 
    #   dplyr::summarise(
    #     geometry = sf::st_combine(sf::st_union(geometry))
    #   ) %>%
    #   rmapshaper::ms_explode(sys = TRUE, sys_mem = 16) %>% 
    #   dplyr::ungroup() %>% 
    #   nngeo::st_remove_holes(max_area = 200)
      
  # fema_vpu <- 
  #   fema_vpu %>% 
  #   dplyr::mutate(      
  #     vpu      = gsub("VPU_", "", VPU),
  #     fema_id  = as.character(1:dplyr::n())
  #   ) %>% 
  #   dplyr::select(
  #       vpu, fema_id,
  #       # state, 
  #       geom = geometry
  #       )
  
  
  # fema_vpu <- 
  #   fema_vpu %>% 
  #   nngeo::st_remove_holes(max_area = 200) %>%
  #   # dplyr::select(geometry = geom) %>%
  #   add_predicate_group_id(sf::st_intersects) %>% 
  #   sf::st_make_valid() %>%
  #   dplyr::group_by(group_id) %>% 
  #   dplyr::summarise(
  #     geometry = sf::st_combine(sf::st_union(geometry))
  #   ) %>% 
  #   dplyr::ungroup() %>% 
  #   dplyr::select(-group_id) %>%
  #   add_predicate_group_id(sf::st_intersects)
  # 
  # fema_vpu %>% sf::st_geometry_type() %>% unique()
  # fema_vpu %>% mapview::npts()
  # 
  # geom_type_counts <- table(sf::st_geometry_type(fema_vpu))
  # 
  # message("Geometry counts before casting all geometries to MULTIPOLYGON:")
  # for (g in seq_along(geom_type_counts)) {
  #   message(" > ", names(geom_type_counts[g]), ": ", geom_type_counts[g])
  # }
  # 
  # message("Keeping only POLYGON and MULTIPOLYGON geometries...")
  # fema_vpu <- 
  #   fema_vpu %>%
  #   dplyr::filter(sf::st_geometry_type(geometry) %in% c("POLYGON", "MULTIPOLYGON")) %>%
  #   sf::st_cast("MULTIPOLYGON") %>%
  #   sf::st_make_valid() %>%
  #   # dplyr::group_by(group_id) %>% 
  #   rmapshaper::ms_dissolve(sys = TRUE, sys_mem = 16) %>%
  #   rmapshaper::ms_explode(sys = TRUE, sys_mem = 16) %>%
  #   nngeo::st_remove_holes(max_area = 200) %>% 
  #   dplyr::mutate(
  #     fema_id = as.character(1:dplyr::n())
  #   ) %>% 
  #   dplyr::select(fema_id, geometry)
  
  # end_geom_type_counts <- table(sf::st_geometry_type(fema_vpu))
  # message("Geometry counts after all processing steps: ")
  # for (g in seq_along(end_geom_type_counts)) {
  #   message(" > ", names(end_geom_type_counts[g]), ": ", end_geom_type_counts[g])
  # }
  
  # fema_vpu2 %>% mapview::npts()
  # fema_vpu2_subset <- fema_vpu2[lengths(sf::st_intersects(fema_vpu2, fema_vpu[1:100, ]))  > 1, ]
  # mapview::mapview(fema_vpu, color = 'red', col.regions = 'white') +
  #       mapview::mapview(fema_vpu2, color = 'green', col.regions = 'white')
  #   # mapview::mapview(fema_vpu2_subset, color = 'green', col.regions = 'white')
  # # mapview::mapview(fema_vpu2[1:100, ], color = 'green', col.regions = 'white')
  
  # fema_vpu <-
  #   fema_vpu %>%
  #   # dplyr::group_by(source) %>%
  #   dplyr::mutate(
  #     # state    = tolower(gsub("-100yr-flood_valid_clean.gpkg", "", source)),
  #     vpu      = gsub("VPU_", "", VPU),
  #     fema_id  = 1:dplyr::n()
  #   ) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::select(vpu, fema_id, 
  #                 # state, 
  #                 geom = geometry)
  
  if (OVERWRITE_FEMA_FILES) {
    
    message("> Overwritting '", basename(master_filepath), "' with final clean version...")
    
    # union_file_path <- gsub(".gpkg", "_union.gpkg", fema_vpu_file)
    # message("> writting '", basename(union_file_path), "' (unioned and exploded version)")
    
    sf::write_sf(
      fema_vpu,
      updated_filepath
      # master_filepath
      # union_file_path
    )
  }
  message()
}
