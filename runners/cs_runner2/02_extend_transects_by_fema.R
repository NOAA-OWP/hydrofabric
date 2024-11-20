
library(sf)
library(dplyr)

source("runners/cs_runner2/base_variables.R")
# source("runners/cs_runner2/utils.R")

VPU_IDS <- get_vpu_ids()

for (i in seq_along(VPU_IDS)) {
  # for (i in 21:length(VPU_IDS)) {
  
  # --------------------------------------------------------------
  # ---- Setup variables / paths ----
  # --------------------------------------------------------------
  # VERSION_DIRS_LIST
  # i = 8
  vpu                       <- VPU_IDS[i]
  CROSSWALK_ID              <- "id"
  # cs_filesnames <- get_cross_section_filenames(vpu, sep = "-")
  transect_filenames           <- get_transect_filenames(vpu, sep = "-")
 
  # TODO: Add a catch for not allowing this all to run if the base_transects do NOT exist 
  base_transects_path          <- paste0(VERSION_DIRS_LIST$transects_base_dir, "/", transect_filenames$transects_base_path)
  fema_transects_output_path   <- paste0(VERSION_DIRS_LIST$transects_fema_extended_dir, "/", transect_filenames$transects_fema_extended_path)
  
  base_transect_file_exists    <- file.exists(base_transects_path)
  fema_transects_file_exists   <- file.exists(fema_transects_output_path)
  
  do_process_transects      <- !fema_transects_file_exists || REGENERATE_TRANSECTS
  
  if (do_process_transects) {
    message("Creating VPU ", vpu, " transects:", 
            "\n - flowpaths: '",
            basename(CONUS_NEXTGEN_GPKG_PATH), "'",
            "\n - base transects: '",
            basename(base_transects_path), "'"
    )
  } else {
    message(
      "VPU ", vpu, " transects file already exists at:\n - '", fema_transects_output_path, "'", 
      "\n\n   >>> NOTE: To regenerate transects:\n\t\t- Set REGENERATE_TRANSECTS = TRUE in 'base_variables.R' \n\t\tOR\n\t\t- delete the '", 
      basename(fema_transects_output_path), "' file"
    )
    next
    
  }
  
  if (do_process_transects && fema_transects_file_exists) {
    message("Deleting old fema extended transect file:\n > '", fema_transects_output_path, "'")
    file.remove(fema_transects_output_path)
  }
  
  # read in nextgen flowlines 
  flines <- 
    CONUS_NEXTGEN_GPKG_PATH  %>% 
    sf::read_sf(query = paste0("SELECT * FROM flowpaths WHERE vpuid = '", vpu, "'")) %>% 
    hydroloom::rename_geometry("geometry")
  
  # transects
  transects <- sf::read_sf(base_transects_path) %>% 
    hydroloom::rename_geometry("geometry")
  
  # FEMA polygons
  fema <- 
    CONUS_FEMA_GPKG_PATH %>% 
    sf::read_sf(layer = get_fema_conus_layer_name(vpu)) %>% 
    rmapshaper::ms_simplify(keep_shapes = T, keep = 0.01, sys = TRUE, sys_mem = 16) %>% 
    hydroloom::rename_geometry("geometry")
  # fema <- rmapshaper::ms_simplify(fema, keep_shapes = T, keep = 0.1, sys = TRUE, sys_mem = 16)
  
  # # add mainstem to  
  # transects <- 
  #   transects  %>%
  #   dplyr::left_join(
  #       flines %>% 
  #         sf::st_drop_geometry() %>% 
  #         dplyr::select(dplyr::any_of(CROSSWALK_ID), mainstem),
  #     by = CROSSWALK_ID 
  #   )
  
  message("Extending transects out to FEMA 100yr floodplain polygon boundaries - (", Sys.time(), ")")
  
  # TODO: make sure this 3000m extension distance is appropriate across VPUs 
  # TODO: also got to make sure that this will be feasible on memory on the larger VPUs...
  extended_transects <- hydrofabric3D::extend_transects_to_polygons(
    transect_lines         = transects, 
    polygons               = fema, 
    flowlines              = flines, 
    crosswalk_id           = CROSSWALK_ID,
    grouping_id            = CROSSWALK_ID, 
    # grouping_id            = "mainstem", 
    max_extension_distance = 3000,
    reindex_cs_ids = TRUE
  )
  
  is_valid_transects <- hydrofabric3D::validate_transects(extended_transects, "id")
  is_valid_transects_against_flowlines <- hydrofabric3D::validate_transects_against_flowlines(extended_transects, flines, "id")
  
  
  
  
  
}

  
  
  
  
  
  
  
  
  
  
  
  
  