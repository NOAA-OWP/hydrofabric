
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
  vpu                       <- VPU_IDS[i]
  CROSSWALK_ID              <- "id"
  
  transect_filenames        <- get_transect_filenames(vpu, sep = "-")
  base_output_path          <- paste0(VERSION_DIRS_LIST$transects_base_dir, "/", transect_filenames$transects_base_path)
  base_transect_file_exists <- file.exists(base_output_path)
  
  do_process_transects      <- !base_transect_file_exists || REGENERATE_TRANSECTS
  # do_not_reprocess_transects <- base_transect_file_exists && !REGENERATE_TRANSECTS
  
  if (do_process_transects) {
    message("Creating VPU ", vpu, " transects:", 
            "\n - flowpaths: '",
            basename(CONUS_NEXTGEN_GPKG_PATH), "'"
    )
  } else {
    message(
      "VPU ", vpu, " transects file already exists at:\n - '", base_output_path, "'", 
      "\n\n   >>> NOTE: To regenerate transects:\n\t\t- Set REGENERATE_TRANSECTS = TRUE in 'base_variables.R' \n\t\tOR\n\t\t- delete the '", 
      basename(base_output_path), "' file"
    )
    next
    
  }
  
  # if we are going to create new transects and the transect file already exists, delete it as to not cause confusion or accidently append to an old dataset
  if (do_process_transects && base_transect_file_exists) {
    message("Deleting old transect file:\n > '", base_output_path, "'")
    file.remove(base_output_path)
  }
  
  # --------------------------------------------------------------
  # ---- Load flowpaths for VPU from conus_nextgen.gpkg ----
  # --------------------------------------------------------------
  
  # read in nextgen flowlines 
  flines <- 
    CONUS_NEXTGEN_GPKG_PATH  %>% 
    sf::read_sf(query = paste0("SELECT * FROM flowpaths WHERE vpuid = '", vpu, "'"))
  
  has_no_flowlines <- nrow(flines) == 0 
  
  if (has_no_flowlines) {
      message("Skipping VPU ", vpu, " as no flowlines were found in '", basename(CONUS_NEXTGEN_GPKG_PATH), "'..." )
      next
    }
    
  
  
  # Add an estimate bankful width based on Total downsteam drainage area (sqkm, Power law equation) 
  flines <-
    flines %>% 
    hydroloom::rename_geometry("geometry") %>% 
    hydrofabric3D::add_powerlaw_bankful_width(
      total_drainage_area_sqkm_col = "tot_drainage_areasqkm", 
      min_bf_width = 50
    ) %>% 
    dplyr::group_by(order) %>% 
      dplyr::slice(1:10) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(
        id,
        lengthkm,
        tot_drainage_areasqkm,
        bf_width,
        mainstem,
        geometry
      ) 
    # dplyr::group_by(order)
  
  # --------------------------------------------------------------
  # ---- Generate transects ---- 
  # --------------------------------------------------------------
  # flines$geometry %>% sf::st_is_empty() %>% any()
  
  # sf::write_sf(
  #   flines,
  #   "/Users/anguswatters/Desktop/wrong_cs_ids_flines_error.gpkg"
  #   # "/Users/anguswatters/Desktop/empty_geom_flines_error.gpkg"
  # )
  
  # create transect lines
  transects <- hydrofabric3D::cut_cross_sections(
    net               = flines,                        # flowlines network
    crosswalk_id      = CROSSWALK_ID,                       # Unique feature ID
    cs_widths         = flines$bf_width,     # cross section width of each "id" linestring ("hy_id")
    # cs_widths = 15,
    num               = 3,                            # number of cross sections per "id" linestring ("hy_id")
    # smooth            = FALSE,                          # smooth lines
    # densify           = NULL,  
    smooth            = TRUE,                          # smooth lines
    densify           = 3,                             # densify linestring points
    
    rm_self_intersect = TRUE,                          # remove self intersecting transects
    fix_braids        = FALSE,                         # whether to fix braided flowlines or not
    add               = TRUE                           # whether to add back the original data
  )  
  # dplyr::mutate(
    # cs_source = CS_SOURCE
  # )
  # transects$cs_lengthm 
  # transects$geometry %>% sf::st_length()
  
  # transects %>% 
  #   dplyr::mutate(
  #     new_cs_lengthm = as.numeric(sf::st_length(.))
  #   ) %>% 
  #   dplyr::filter(!dplyr::near(new_cs_lengthm, cs_lengthm, tol =   2))
  # .Machine$double.eps^1
  # crosswalk_id <- "id"
  # # reenumerate the cs_ids for each transect based on cs_measure sorting, and make sure all cross sections are correctly numbered
  # mismatches <-
  #   transects %>%
  #   sf::st_drop_geometry() %>% 
  #   dplyr::group_by(dplyr::across(dplyr::any_of(c(crosswalk_id)))) %>% 
  #   dplyr::arrange(cs_measure, .by_group = TRUE) %>% 
  #   dplyr::mutate(
  #     new_cs_id = 1:dplyr::n()
  #   ) %>% 
  #   dplyr::ungroup() %>% 
  #   dplyr::filter(cs_id != new_cs_id)
  # mismatches
  # transects %>% 
  #   dplyr::filter(id %in% c("wb-1023360", "wb-1023363"))
  # 
  # transects %>% 
  #   dplyr::filter(id %in% c("wb-1023360", "wb-1023363")) %>% 
  #   mapview::mapview()
  # 
  # # FALSE if there are any transects with different cs_ids to the newly created cs_id 
  # # Otherwise TRUE
  # has_valid_cs_ids <- !(nrow(mismatches) > 0)
  # 
  # # t2 <- transects %>% hydrofabric3D:::renumber_cs_ids("id")
  # # hydrofabric3D::validate_transects(t2, "id")
  # transects 
  is_valid_transects                    <- hydrofabric3D::validate_transects(transects, "id")
  is_valid_transects_against_flowlines  <- hydrofabric3D::validate_transects_against_flowlines(transects, flines, "id")
  
  # is_valid_transects                    <- hydrofabric3D::validate_transects(transects[c(1, 3, 55), ], "id")
  # is_valid_transects_against_flowlines  <- hydrofabric3D::validate_transects_against_flowlines(transects, flines[c(-1), ], "id")
  # is_valid_transects_against_flowlines
  # hydrofabric3D:::rm_multiflowline_intersections(transects, flines[c(-1), ])
  
  # select relevent columns
  transects <-
    transects %>% 
    dplyr::mutate(
      cs_source = CS_SOURCE
    ) %>% 
    # hydrofabric3D:::select_transects(CROSSWALK_ID)
    dplyr::select(
      dplyr::any_of(CROSSWALK_ID), 
      cs_id, 
      cs_lengthm,
      # cs_lengthm = new_cs_lengthm, 
      cs_measure,
      sinuosity,
      cs_source,
      geometry
    )
  
  if (is_valid_transects && is_valid_transects_against_flowlines) {
    message("Saving transects to:\n - filepath: '", base_output_path, "'")
      sf::write_sf(
        transects,
        base_output_path
      )
  }
  
  message(" --- VPU ", vpu, " transects generation complete! --- ")
  message()
  
}
