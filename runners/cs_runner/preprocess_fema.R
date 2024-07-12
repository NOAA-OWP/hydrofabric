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

# TODO: Steps that converts FGB to geojson and then geojson to gpkg can be put into a single loop
# TODO: Delete old files as needed

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
  # message(VPU_SUBFOLDER)
  
  state_dir  = paste0(VPU_SUBFOLDER, "/states/")
  merged_dir = paste0(VPU_SUBFOLDER, "/merged/")
  
  if (!dir.exists(VPU_SUBFOLDER)) {
    message("Creating FEMA VPU subfolder...")
    message(paste0("'/", basename(VPU_SUBFOLDER), "' directory does not exist...\n  Creating directory:\n > '", VPU_SUBFOLDER, "'"))
    dir.create(VPU_SUBFOLDER)
  }
  
  if (!dir.exists(state_dir)) { 
      message("Creating FEMA VPU states subfolder...")
      message(paste0("'/", basename(state_dir), "' directory does not exist...\n  Creating directory:\n > '", state_dir, "'"))
      
      dir.create(state_dir)
     
     }
  
  if (!dir.exists(merged_dir)) { 
    message("Creating FEMA VPU merged subfolder...")
    message(paste0("'/", basename(merged_dir), "' directory does not exist...\n  Creating directory:\n > '", merged_dir, "'"))
    
    dir.create(merged_dir)
    
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
  # message(file)

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
  
  # start_fema <- sf::read_sf(file)
  
  # mapshaper_command = paste0('node  --max-old-space-size=16000 /opt/homebrew/bin/mapshaper ', file, 
  #                            ' -simplify 0.15 visvalingam \\', 
  #                            ' -dissolve \\', 
  #                            ' -explode \\', 
  #                            ' -o ', output_path
  # )
  
  mapshaper_command = paste0('node  --max-old-space-size=16000 /opt/homebrew/bin/mapshaper ', file, 
                             ' -dissolve2 FLD_AR_ID \\', 
                             ' -simplify 0.1 visvalingam \\', 
                             # ' -explode \\',
                             ' -snap \\',
                             ' -o ', output_path
  )
  
  if (OVERWRITE_FEMA_FILES) {
    message("Running mapshaper 'simplify', 'dissolve', and 'explode' via CLI...")
    system(mapshaper_command)
    message("Writting '", output_clean_filename, "' to: \n > '", output_path, "'")
  }
  # end_fema <- sf::read_sf(output_path)
  
  message()

  }

# -------------------------------------------------------------------------------------
# ---- Convert cleaned FEMA geojson geometries to geopackages ----
# -------------------------------------------------------------------------------------

# paths to FEMA 100 year flood plain files
FEMA_clean_paths      <- list.files(FEMA_CLEAN_PATH, full.names = TRUE)

for (file in FEMA_clean_paths) {
  message("Fema 100 year flood plain:\n > '", basename(file), "'")
  
  output_gpkg_filename <- gsub("_clean.geojson", "_clean.gpkg", basename(file))
  output_path          <- paste0(FEMA_GPKG_PATH, "/", output_gpkg_filename)
  
  message("Converting GEOJSON file to GPKG:\n > '", basename(file), "' > '", output_gpkg_filename, "'")
  
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

# # -------------------------------------------------------------------------------------------------------------------
# # ---- Apply final dissolve/snap and removal of internal boudnaries in FEMA geometries  ----
# # -------------------------------------------------------------------------------------------------------------------

# paths to FEMA 100 year flood plain files
FEMA_gpkg_paths      <- list.files(FEMA_GPKG_PATH, full.names = TRUE)

for (file_path in FEMA_gpkg_paths) {
  message("Applying hydrofab::clean_geometry() to:\n > '", basename(file_path), "'")
   
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
  # mapview::mapview(start_fema$geom, color = "red", col.regions = "red") + 
  #   mapview::mapview(end_fema$geom, color = 'limegreen', col.regions = "limegreen") + 
  #   mapview::mapview(snap_union_sf, color = 'gold', col.regions = "gold") + 
  # mapview::mapview(final_fema, color = 'white', col.regions = "white") + 
  # mapview::mapview(fin, color = 'white', col.regions = "white") 

  # message(" > ", nrow(fema), " POLYGONs")
  # message("Start time: ", Sys.time())
  # 
  # fema_clean <- hydrofab::clean_geometry(
  #   catchments = fema,
  #   ID         = "fema_id"
  #   )
  # 
  # fema_clean <-
  #   fema_clean %>%
  #   dplyr::mutate(
  #     source = basename(file_path),
  #     state  = gsub("-100yr-flood_valid_clean.gpkg", "", source)
  #   ) %>%
  #   dplyr::select(fema_id, source, state, areasqkm, geometry)
  
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
  
  # geom_diff <- sf::st_difference(fema[1, ], fema_clean[1, ])
  # mapview::mapview(fema[1, ], col.regions = "red") +
  # mapview::mapview(fema_clean[1, ], col.regions = "green") +
  # mapview::mapview(geom_diff, col.regions = "white")
  
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

# # -------------------------------------------------------------------------------------
# # ---- Partion parts of each FEMA GPKGs to the a Nextgen VPU ---- 
# # -------------------------------------------------------------------------------------

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
      vpu_subfolder_path <- paste0(FEMA_BY_VPU_PATH, "/", vpu_subfolder, "/states")
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

DELETE_STAGING_GPKGS <- FALSE 

for (vpu_dir in FEMA_VPU_SUBFOLDERS) {
 # for (i in 1:4) {
  # vpu_dir = FEMA_VPU_SUBFOLDERS[i]
  message("Merging files in '", basename(vpu_dir), "' directory...")
 # }
  
  # vpu_dir <- '/Users/anguswatters/Desktop/lynker-spatial/FEMA_BY_VPU/VPU_06'
  vpu_subdirs    <- list.files(vpu_dir, full.names = TRUE)
  
  states_dir <- vpu_subdirs[grepl(paste0(vpu_dir, "/states"), vpu_subdirs)]
  merged_dir <- vpu_subdirs[grepl(paste0(vpu_dir, "/merged"), vpu_subdirs)]
  
  # fema state geopackages partioned for the specific VPU
  fema_state_gpkgs <- list.files(states_dir, full.names = TRUE)
  
  master_name      <- paste0("fema_", gsub("VPU", "vpu", basename(vpu_dir)))
  master_gpkg_name <- paste0(master_name, ".gpkg")
  
  # path to the merged directory where the final merged geopackge will end up
  master_filepath  <- paste0(merged_dir, "/", master_gpkg_name)
  
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
  
  if(DELETE_STAGING_GPKGS) {
    message(" - Deleting individual gpkgs from '/states' directory...")
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

MERGED_DIRS <- paste0(FEMA_VPU_SUBFOLDERS, "/merged")
MERGED_DIRS
for (i in 1:length(FEMA_VPU_SUBFOLDERS)) {
  # i = 8
  vpu_dir = FEMA_VPU_SUBFOLDERS[i]
  
  VPU        <- basename(vpu_dir)
  
  message("Attempting to union FEMA polygons for '", VPU, "'...")
  
  merged_dir <- paste0(vpu_dir, "/merged")
  fema_vpu_file <- list.files(merged_dir, full.names = TRUE)
  
  has_fema_vpu_file <- ifelse(length(fema_vpu_file) > 0, TRUE, FALSE)
  # has_fema_vpu_file
  # message()
  # fema_vpu_file
# }
  if(!has_fema_vpu_file) { 
    message("No FEMA geometries in '", VPU, "'")
    message()
    next
  }

  message("> Re-unioning and re-exploding geometries in '", basename(fema_vpu_file), "'")
  
  fema_vpu_file <- fema_vpu_file[!grepl("_union.gpkg", fema_vpu_file)]

  fema_vpu <- sf::read_sf(fema_vpu_file)
  # fema_vpu <-
  #   fema_vpu %>%
  #   dplyr::group_by(source) %>%
  #   dplyr::summarise()  %>%
  #   dplyr::ungroup()
 
   mapview::npts(fema_vpu)
  fema_vpu2 <- 
    fema_vpu %>% 
    nngeo::st_remove_holes(max_area = 20) %>%
    # dplyr::select(geometry = geom) %>%
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
  
  fema_vpu2 %>% mapview::npts()
  
  fema_vpu2_subset <- fema_vpu2[lengths(sf::st_intersects(fema_vpu2, fema_vpu[1:100, ]))  > 1, ]
  
  mapview::mapview(fema_vpu[1:100, ], color = 'red', col.regions = 'white') +
    mapview::mapview(fema_vpu2_subset, color = 'green', col.regions = 'white')
      # mapview::mapview(fema_vpu2[1:100, ], color = 'green', col.regions = 'white')
    
  # message("Removing holes before dissolve...")
  fema_vpu <- nngeo::st_remove_holes(fema_vpu)
  mapview::mapview(fema_vpu[1:100, ], color = 'red', col.regions = 'white') +
  mapview::mapview(fema_vpu2[1:100, ], color = 'green', col.regions = 'white') 
  # mapview::npts(fema_vpu)
  all(sf::st_is_valid(fema_vpu))
  rmapshaper::ms_innerlines(fema_vpu)
  
  
  # message("Making valid geometries...")
  # fema_vpu <- sf::st_make_valid(fema_vpu) 
  
  # fema_vpu <- 
  #   fema_vpu %>% 
  #   sf::st_cast("MULTIPOLYGON")
  
  message("Dissolving...")
  
  # 1421 = old number of polygons
  fema_vpu2 <- rmapshaper::ms_dissolve(
                            input   = fema_vpu,
                            field   = "source",
                            sys     = TRUE,
                            sys_mem = 16
                            )
  
  message("Exploding...")
  
  # mapview::npts(fema_vpu)
  # mapview::npts(fema_vpu2)
  
  fema_vpu2 <- rmapshaper::ms_explode(
                            input   = fema_vpu2,     
                            sys     = TRUE,
                            sys_mem = 16
                            )
  # mapview::npts(fema_vpu2)
  message("Removing holes after explosion...")
  fema_vpu2 <- nngeo::st_remove_holes(fema_vpu2)
  
  fema_vpu2 <- 
    fema_vpu2 %>% 
    add_predicate_group_id(sf::st_intersects) %>% 
    dplyr::group_by(group_id) %>% 
    dplyr::summarise(
      geometry = sf::st_combine(sf::st_union(geometry))
    )
  
  mapview::mapview(fema_vpu[1:100, ], color = 'red', col.regions = 'white') +
    mapview::mapview(fema_vpu2[1:00, ], color = 'green', col.regions = 'white') 
  
  # mapview::npts(fema_vpu2)
  sf::st_is_valid(fema_vpu2) %>% all()
  
  fema_vpu2 %>% 
    sf::st_make_valid() %>% 
    sf::st_geometry_type() %>% 
    unique()
  sf::st_geometry_type(fema_vpu2) %>% unique()
  # slice_subset = 1:50
  # fema_exp_noholes[slice_subset, ]
  # mapview::mapview(  fema_vpu[1:100, ], col.regions = "dodgerblue")+ 
  # mapview::mapview(  fema_exp[slice_subset, ], col.regions = "red") +
  #   mapview::mapview(  fema_exp_noholes[slice_subset, ], col.regions = "green")
  # fema_vpu <- rmapshaper::ms_dissolve(fema_vpu,
  #                                    field = "source",
  #                                    sys = TRUE,
  #                                    sys_mem = 16
  # # )
  # fema_vpu  <- rmapshaper::ms_explode(fema_vpu,     
  #                                     sys = TRUE,
  #                                     sys_mem = 16)

  fema_vpu <-
    fema_vpu %>%
    # dplyr::group_by(source) %>%
    dplyr::mutate(
      state    = tolower(gsub("-100yr-flood_valid_clean.gpkg", "", source)),
      vpu      = gsub("VPU_", "", VPU),
      fema_id  = 1:dplyr::n()
      ) %>%
    dplyr::ungroup() %>%
    dplyr::select(vpu, fema_id, state, geom = geometry)
  
  if (OVERWRITE_FEMA_FILES) {
    union_file_path <- gsub(".gpkg", "_union.gpkg", fema_vpu_file)
    message("> Writting '", basename(union_file_path), "' (unioned and exploded version)")
    sf::write_sf(
      fema_vpu,
      union_file_path
      )
  }
  message()
}

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
