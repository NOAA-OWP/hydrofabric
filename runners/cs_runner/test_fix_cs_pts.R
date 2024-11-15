# Generate the flowlines layer for the final cross_sections_<VPU>.gpkg for each VPU
source("runners/cs_runner/config.R")

# # load libraries
library(hydrofabric3D)
library(dplyr)
library(sf)

# paths to nextgen datasets
NEXTGEN_FILES  <- list.files(NEXTGEN_DIR, full.names = FALSE)

# paths to nextgen datasets
transect_files <- list.files(TRANSECTS_DIR, full.names = FALSE)
transect_files <- transect_files[!grepl("updated_", transect_files)]

REF_FEATURES   <- list.files(REF_FEATURES_GPKG_DIR, full.names = FALSE)

START_CS_PTS_DIR <-  '/Users/anguswatters/Desktop/lynker-spatial/test_out'
start_cs_pts_paths     <- list.files(START_CS_PTS_DIR, full.names = T)

END_CS_PTS_DIR <-  '/Users/anguswatters/Desktop/lynker-spatial/test_fix_out'

# reference features dataframe
ref_df <- data.frame(
  vpu      = sapply(strsplit(REF_FEATURES, "_", fixed = TRUE), function(i) { i[1] }), 
  ref_file = REF_FEATURES
)

# ensure the files are in the same order and matched up by VPU
path_df <- align_files_by_vpu(
  x    = NEXTGEN_FILES,
  y    = transect_files,
  base = BASE_DIR
) %>% 
  dplyr::left_join(
    ref_df,
    by = "vpu"
  )
path_df$cs_pts_path <- start_cs_pts_paths

# path_df

# loop over the nextgen and transect datasets (by VPU) and extract point elevations across points on each transect line,
# then classify the points, and create a parquet file with hy_id, cs_id, pt_id, X, Y, Z data.
# Save parquet locally and upload to specified S3 bucket

# output_path <- paste0(BASE_DIR, "/test_out/")

for (i in 15:nrow(path_df)) {
  # i = 8
  
  start <- Sys.time()
  
  start_cs_pts_path <- path_df$cs_pts_path[i]
  
  # nextgen file and full path
  nextgen_file <- path_df$x[i]
  nextgen_path <- paste0(NEXTGEN_DIR, nextgen_file)
  
  # model attributes file and full path
  transect_file <- path_df$y[i]
  transect_path <- paste0(TRANSECTS_DIR, transect_file)
  
  # model attributes file and full path
  ref_file <- path_df$ref_file[i]
  ref_path <- paste0(REF_FEATURES_DIR, "gpkg/", ref_file)
  
  # current VPU being processed
  VPU   <- path_df$vpu[i]
  
  start <- Sys.time()
  
  message("Creating VPU ", VPU, 
          " cross section points:\n - flowpaths: '", nextgen_file,
          "'\n - start cs pts: '", basename(start_cs_pts_path), "'", 
          "'\n - transects: '", transect_file, "'", 
          "\n - waterbodies: '", ref_file, "'",
          "'\n - start time: '", start, "'"
  )
  
  ################### 
  message("Reading in starting CS points...\n > ", basename(start_cs_pts_path))
  cs_pts <- arrow::read_parquet(start_cs_pts_path)
  
  message("Reading in transects...\n > ", transect_file)
  # read in transects data
  transects <- sf::read_sf(transect_path)
  
  message("Reading in flowlines... \n > ", nextgen_file)
  # read in nextgen data
  flines <- sf::read_sf(nextgen_path, layer = "flowpaths")
  
  message("Reading in waterbodies... \n > ", ref_file)
  # read in waterbodies reference features layer
  waterbodies <- sf::read_sf(ref_path, layer = "waterbodies")
  
  # Update flowlines and transects to remove flowlines and transects that intersect with reference_features waterbodies
  feature_subsets <- wb_intersects(flines, transects, waterbodies)
  
  # replace flowlines and transects objects with updated versions in "updated_features"
  flines    <- flines[feature_subsets$valid_flowlines, ]
  transects <- transects[feature_subsets$valid_transects, ]
  
  rm(waterbodies)
  gc()
  
  start_cs_pts <- Sys.time()
  flines <- dplyr::rename(flines, hy_id = id)
  
  # test_ids <- 
  #   flines %>% 
  #   dplyr::slice(1:100) %>% 
  #   .$hy_id %>% 
  #   unique()
  # flines <-
  #   flines %>% 
  #   dplyr::filter(
  #     hy_id %in% test_ids
  #   )
  # transects <-
  #   transects %>% 
  #   dplyr::filter(
  #     hy_id %in% test_ids
  #   )
  # cs_pts <-
  #   cs_pts %>% 
  #   dplyr::filter(
  #     hy_id %in% test_ids
  #   )
  
  ids_original_cs_pts <- hydrofabric3D::add_tmp_id(cs_pts)$tmp_id 
  
  # # ------------------------------------------------------------------------
  # # ------ TESTING DATA -------
  # # ------------------------------------------------------------------------
  # system.time({
    fixed_pts <- hydrofabric3D::get_improved_cs_pts(
      cs_pts         = cs_pts,    # cross section points generated from hydrofabric3D::cross_section_pts()
      net            = flines,    # original flowline network
      # net            = flines,    # original flowline network
      transects      = transects, # original transect lines
      crosswalk_id   = "hy_id",
      points_per_cs  = NULL, 
      min_pts_per_cs = 10, # number of points per cross sections
      dem            = DEM_PATH, # DEM to extract points from
      scale          = EXTENSION_PCT, # How far to extend transects if the points need to be rechecked
      pct_of_length_for_relief = PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF, # percent of cross sections length to be needed in relief calculation to consider cross section to "have relief"
      fix_ids = FALSE,
      verbose = TRUE
    )
  # })
  # fixed_pts 
  
  # fixed_pts %>%
  #   # cs_pts %>% 
  #   sf::st_drop_geometry() %>%
  #   # dplyr::left_join(
  #   #   stream_order,
  #   #   by = "hy_id"
  #   # ) %>%
  #   dplyr::select(dplyr::any_of("hy_id"), cs_id, 
  #                 # vpu, order,
  #                 valid_banks, has_relief) %>%
  #   dplyr::group_by(
  #     dplyr::across(dplyr::any_of(c("hy_id", "cs_id")))
  #   ) %>%
  #   dplyr::slice(1) %>%
  #   dplyr::ungroup() %>%
  #   # dplyr::group_by(
  #   #   dplyr::across(dplyr::any_of(c("vpu", "order")))
  #   # ) %>%
  #   dplyr::count(valid_banks, has_relief)
  
  
  # transects %>% 
  #   sf::st_drop_geometry() %>% 
  #   dplyr::group_by(id, cs_id)
  # fixed_pts2$is_extended %>% sum()
  ids_after_fixed_pts <- hydrofabric3D::add_tmp_id(fixed_pts)$tmp_id
  # ----------------------------------------------------------------------------------------------------------------
  # ---- Cross section points parquet to S3 ----
  # ----------------------------------------------------------------------------------------------------------------
  
  # classify the cross section points
  fixed_pts <-
    fixed_pts %>% 
    # dplyr::mutate(
    #   X = sf::st_coordinates(.)[,1],
    #   Y = sf::st_coordinates(.)[,2]
    # ) %>%
    sf::st_drop_geometry() %>% 
    dplyr::select(
      hy_id, cs_id, pt_id,
      cs_lengthm,
      relative_distance,
      X, Y, Z,
      class, point_type,
      bottom, left_bank, right_bank, valid_banks, has_relief # newly added columns (03/06/2024)
    )
  
  # # Drop point geometries, leaving just X, Y, Z values
  # fixed_pts <- sf::st_drop_geometry(fixed_pts)
  
  # add Z_source column for source of elevation data
  fixed_pts <-
    fixed_pts %>%
    dplyr::mutate(
      Z_source = CS_SOURCE
    ) %>%
    dplyr::relocate(hy_id, cs_id, pt_id, cs_lengthm, relative_distance, X, Y, Z, Z_source, 
                    class, point_type, 
                    bottom, left_bank, right_bank, valid_banks, has_relief)
  
  ids_before_align <- hydrofabric3D::add_tmp_id(fixed_pts)$tmp_id
  
  message("Aligning banks and smoothing bottoms...")
  fixed_pts <- hydrofabric3D::align_banks_and_bottoms(cs_pts = fixed_pts, crosswalk_id = "hy_id")
  # fixed_pts <- hydrofabric3D::align_banks_and_bottoms(cs_pts = fixed_pts, crosswalk_id = "hy_id")
  
  ids_after_align <- hydrofabric3D::add_tmp_id(fixed_pts)$tmp_id
  
  message("Reclassifying cross section points...")
  fixed_pts <- hydrofabric3D::classify_points(
    cs_pts                    = fixed_pts, 
    crosswalk_id              = "hy_id",
    pct_of_length_for_relief  = PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF
  )
  
  ids_after_reclassify <- hydrofabric3D::add_tmp_id(fixed_pts)$tmp_id
  
  if(all(ids_original_cs_pts %in% ids_after_fixed_pts)) {
    message("All hy_id/cs_ids in ORIGINAL DEM point extraction were found in the FIXED points")
  } else {
    message(" >>> Missing hy_id/cs_ids in ORIGINAL DEM point extraction compared to the FIXED points")
  }
  
  if(all(ids_before_align %in% ids_after_align)) {
    message("All hy_id/cs_ids are kept in tact after bank alignment and bottom smoothing")
  } else {
    message(" >>> Missing hy_id/cs_ids after bank alignment and bottom smoothing")
  }
  
  if(all(ids_after_align %in% ids_after_reclassify)) {
    message("All hy_id/cs_ids are kept in tact after RECLASSIFICATION")
  } else {
    message(" >>> Missing hy_id/cs_ids after RECLASSIFICATION")
  }
  
  # all(hydrofabric3D::add_tmp_id(fixed_pts2)$tmp_id %in% hydrofabric3D::add_tmp_id(fixed_pts)$tmp_id)
  # all(hydrofabric3D::add_tmp_id(fixed_pts4)$tmp_id %in% hydrofabric3D::add_tmp_id(fixed_pts)$tmp_id)
  
  ############################################################################## 
  
  # ----------------------------------------------------------------------------------------------------------------
  # ---- Upload the cross section points parquet to S3 ----
  # ----------------------------------------------------------------------------------------------------------------
  
  # # name of file and path to save transects gpkg too
  out_file <- paste0("nextgen_fixed_", path_df$vpu[i], "_cross_sections.parquet")
  out_path <- paste0(BASE_DIR, "/test_fix_out/", out_file)
  
  message("Saving cross section points to:\n - filepath: '", out_path, "'")
  
  # save cross section points as a parquet to out_path (lynker-spatial/02_cs_pts/cs_pts_<VPU num>.parquet)
  arrow::write_parquet(fixed_pts, out_path)
  # sf::write_sf(cs_pts2, "/Users/anguswatters/Desktop/test_improve_cs_pts_classified_11.gpkg")
  # sf::write_sf(cs_pts, "/Users/anguswatters/Desktop/test_improve_cs_pts_classified_11_2.gpkg")
  
  end <- Sys.time()
  
  message("Finished cross section point generation for VPU ", VPU)
  message("- Completed at: ", end)
  message("==========================")
  
  rm(cs_pts, fixed_pts)
  gc()
  gc()
}
