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

# loop over the nextgen and transect datasets (by VPU) and extract point elevations across points on each transect line,
# then classify the points, and create a parquet file with hy_id, cs_id, pt_id, X, Y, Z data.
# Save parquet locally and upload to specified S3 bucket
for (i in 20:nrow(path_df)) {
  
  start <- Sys.time()
  
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
          "'\n - transects: '", transect_file, "'", 
          "\n - waterbodies: '", ref_file, "'",
          "'\n - start time: '", start, "'"
  )
  
  ################### 
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
  # # ------------------------------------------------------------------------
  # # ------ TESTING DATA -------
  # # ------------------------------------------------------------------------
  #   flines <-
  #     flines %>%
  #     dplyr::slice(1:3500)
  # 
  # transects <-
  #   transects %>%
  #   dplyr::filter(hy_id %in% flines$id)

  # ------------------------------------------------------------------------
  
  message("Extracting cross section points (", start_cs_pts, ")")
  # ----------------------------------------------------------------------------------------------------------------
  # ---- STEP 1: Extract cs points from DEM ----
  # ----------------------------------------------------------------------------------------------------------------
  # system.time({
  
  # get cross section point elevations
  cs_pts <- hydrofabric3D::cross_section_pts(
    
    cs             = transects,
    crosswalk_id   = "hy_id",
    points_per_cs  = NULL,
    min_pts_per_cs = 10,
    dem            = DEM_URL
  )

  # })
  
  # ----------------------------------------------------------------------------------------------------------------
  # ---- STEP 2: Remove any cross section that has ANY missing (NA) Z values, and classify the points ----
  # ----------------------------------------------------------------------------------------------------------------
  
  # sf::write_sf(cs_pts, "/Users/anguswatters/Desktop/test_improve_cs_pts_11.gpkg")
  # sf::write_sf(flines, "/Users/anguswatters/Desktop/test_improve_flines_11.gpkg")
  # sf::write_sf(transects, "/Users/anguswatters/Desktop/test_improve_transects_11.gpkg")
  
  # sf::write_sf(flines, "/Users/anguswatters/Desktop/test_improve_flines_11_2.gpkg")
  # sf::write_sf(transects, "/Users/anguswatters/Desktop/test_improve_transects_11_2.gpkg")
  # cs_pts %>% 
  #   dplyr::group_by(hy_id, cs_id) %>% 
  #   dplyr::filter(!any(is.na(Z))) %>% 
  #   dplyr::ungroup()
  # 
  # cs_pts %>% 
  #   hydrofabric3D::drop_incomplete_cs_pts("hy_id")
  
  # system.time({
  
  # STEP 2: Remove any cross section that has ANY missing (NA) Z values, and classify the points 
  cs_pts <- 
  # cs_pts2 <- 
    cs_pts %>% 
    # dplyr::group_by(hy_id, cs_id) %>% 
    # dplyr::filter(!any(is.na(Z))) %>% 
    # dplyr::ungroup() %>% 
    hydrofabric3D::drop_incomplete_cs_pts("hy_id") %>% 
    hydrofabric3D::classify_points(
      crosswalk_id             = "hy_id", 
      pct_of_length_for_relief = PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF
      )  
  
  # })
  
  ids_original_cs_pts <- hydrofabric3D::add_tmp_id(cs_pts)$tmp_id  
  # ids_original_cs_pts <- hydrofabric3D::add_tmp_id(cs_pts2)$tmp_id
  
  # sf::write_sf(cs_pts2, "/Users/anguswatters/Desktop/test_improve_cs_pts_classified_11.gpkg")
  # sf::write_sf(cs_pts, "/Users/anguswatters/Desktop/test_improve_cs_pts_classified_11_2.gpkg")
  
  
  # ----------------------------------------------------------------------------------------------------------------
  # ---- STEP 3: Try to rectify any no relief and invalid banks cross sections ----
  # ----------------------------------------------------------------------------------------------------------------
  # dplyr::rename(flines, hy_id = id)
  # profvis::profvis({
  # system.time({
    
  # # NOTE: new inplace method for improving (rectifying) any invalid cross sections where we dont have banks and relief
  # fixed_pts <- hydrofabric3D::improve_invalid_cs2(
  #   cs_pts         = cs_pts,    # cross section points generated from hydrofabric3D::cross_section_pts()
  #   net            = dplyr::rename(flines, hy_id = id),    # original flowline network
  #   # net            = flines,    # original flowline network
  #   transects      = transects, # original transect lines
  #   points_per_cs  = NULL, 
  #   min_pts_per_cs = 10, # number of points per cross sections
  #   dem            = DEM_URL, # DEM to extract points from
  #   scale          = EXTENSION_PCT, # How far to extend transects if the points need to be rechecked
  #   pct_of_length_for_relief = PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF, # percent of cross sections length to be needed in relief calculation to consider cross section to "have relief"
  #   fix_ids = FALSE,
  #   verbose = TRUE
  # )
  
  
  # system.time({
    fixed_pts <- hydrofabric3D::get_improved_cs_pts(
      cs_pts         = cs_pts,    # cross section points generated from hydrofabric3D::cross_section_pts()
      net            = dplyr::rename(flines, hy_id = id),    # original flowline network
      # net            = flines,    # original flowline network
      transects      = transects, # original transect lines
      crosswalk_id   = "hy_id",
      points_per_cs  = NULL, 
      min_pts_per_cs = 10, # number of points per cross sections
      dem            = DEM_URL, # DEM to extract points from
      scale          = EXTENSION_PCT, # How far to extend transects if the points need to be rechecked
      pct_of_length_for_relief = PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF, # percent of cross sections length to be needed in relief calculation to consider cross section to "have relief"
      fix_ids = FALSE,
      verbose = TRUE
    )
  # })
  
  # fixed_pts2$is_extended %>% sum()
  
  ids_after_fixed_pts <- hydrofabric3D::add_tmp_id(fixed_pts)$tmp_id

  # # TODO: This is taking A LOT time to process as inputs get larger, an improvement should be looked into more
  # fixed_pts <- hydrofabric3D::rectify_cs(
  # # cs_pts <- hydrofabric3D::rectify_flat_cs(
  #   cs_pts         = cs_pts,    # cross section points generated from hydrofabric3D::cross_section_pts()
  #   net            = flines,    # original flowline network
  #   transects      = transects, # original transect lines
  #   points_per_cs  = NULL, 
  #   min_pts_per_cs = 10, # number of points per cross sections
  #   dem            = DEM_URL, # DEM to extract points from
  #   scale          = EXTENSION_PCT, # How far to extend transects if the points need to be rechecked
  #   pct_of_length_for_relief = PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF, # percent of cross sections length to be needed in relief calculation to consider cross section to "have relief"
  #   fix_ids = FALSE,
  #   verbose = TRUE
  # )
  
  # get a summary dataframe and print out details message
  rectify_summary <- hydrofabric3D::rectify_summary(cs_pts, fixed_pts, verbose = TRUE)
  rectify_summary <- 
    rectify_summary %>% 
    dplyr::mutate(
      vpu = VPU
    ) %>% 
    dplyr::relocate(vpu, metric, value)
  
  readr::write_csv(rectify_summary, paste0(META_PATH, "nextgen_", VPU, "_cross_sections_metadata.csv"))
  
  # ----------------------------------------------------------------------------------------------------------------
  # ---- STEP 4: Update transects with extended transects (if exists) ----
  # ----------------------------------------------------------------------------------------------------------------
  
  # get the counts of each point type to add this data to the transects dataset
  point_type_counts <- hydrofabric3D::get_point_type_counts(fixed_pts)
  # point_type_counts <- hydrofabric3D::get_point_type_counts(fixed_pts, crosswalk_id = "hy_id")
  
  # # check the number of cross sections that were extended
  # fixed_pts$is_extended %>% table()
  message("Subsetting cross section points generated after extending transects...")
  
  # extract cross section points that have an "is_extended" value of TRUE
  extended_pts <-
    fixed_pts %>% 
    dplyr::filter(is_extended) %>% 
    hydrofabric3D::add_tmp_id()
  
  # extract transects that have a "hy_id" in the "extended_pts" dataset
  update_transects <-
    transects %>% 
    hydrofabric3D::add_tmp_id() %>% 
    dplyr::filter(tmp_id %in% unique(extended_pts$tmp_id))
  
  # if any transects were extended, update the transects dataset, and overwrite local and S3 transects geopackages
  if (nrow(update_transects) > 0) {
    message("Updating ", nrow(update_transects), " transects")
    
    update_transects <-
      update_transects %>%
      # apply extend_by_percent function to each transect line: 
      hydrofabric3D:::extend_by_percent(
        pct        = EXTENSION_PCT,
        length_col = "cs_lengthm"
      )
  
    # Filter down to ONLY points that were finalized and rectified from rectify_cs_pts()
    # remove old transects that have "tmp_id" in "extended_pts" (transects that were unchanged and are "good_to_go")
    # and then replace with old transects with the "update_transects" 
    out_transects <- 
      transects %>% 
      hydrofabric3D::add_tmp_id() %>% 
      # dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id)) %>%
      dplyr::filter(tmp_id %in% unique(hydrofabric3D::add_tmp_id(fixed_pts)$tmp_id)) %>% # Subset down to the remaining tmp_ids in the fixed points
      dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id)) %>% # remove the tmp_ids that we are going add back in with the extended versions of those tmp_ids
      dplyr::bind_rows( # bring in the new updated extended transects
        dplyr::mutate(
          update_transects,
          is_extended = TRUE
        )
      )  
    
  } else {
    
    out_transects <- 
      transects %>% 
      hydrofabric3D::add_tmp_id() %>% 
      dplyr::filter(tmp_id %in% unique(hydrofabric3D::add_tmp_id(fixed_pts)$tmp_id)) %>% 
      dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id))
  }
  
  # finalize new transects
  out_transects <- 
    out_transects %>% 
    dplyr::left_join(
      point_type_counts,
      by = c("hy_id", "cs_id")
    ) %>% 
    dplyr::left_join(
      dplyr::ungroup(
        dplyr::slice( 
          dplyr::group_by(
            dplyr::select(sf::st_drop_geometry(fixed_pts), 
                          hy_id, cs_id, bottom, left_bank, right_bank, valid_banks, has_relief), 
            hy_id, cs_id),
          1)
      ),
      by = c("hy_id", "cs_id")
    ) %>% 
    dplyr::select(hy_id, cs_source, cs_id, cs_measure, cs_lengthm, 
                  # sinuosity, 
                  is_extended, 
                  left_bank_count, right_bank_count, channel_count, bottom_count, 
                  bottom, left_bank, right_bank, valid_banks, has_relief,
                  geom
                  )
  
  # ----------------------------------------------------------------------------------------------------------------
  # ---- Re enumerate the transects & cross section points "cs_id" ----
  # ----------------------------------------------------------------------------------------------------------------
  
  # make a dataframe that has a new_cs_id column that has 
  # the cs_id renumbered to fill in any missing IDs,
  # so each hy_id has cs_ids that go from 1 - number of cross sections on hy_id
  # The dataframe below will be used to join the "new_cs_id" with 
  # the original "cs_ids" in the final cross section POINTS and UPDATED TRANSECTS output datasets
  renumbered_ids <-
    fixed_pts %>% 
    sf::st_drop_geometry() %>% 
    # dplyr::filter(hy_id %in% c("wb-2402800", "wb-2398282", "wb-2400351")) %>%
    dplyr::select(hy_id, cs_id, pt_id, cs_measure) %>% 
    dplyr::group_by(hy_id, cs_id) %>% 
    dplyr::slice(1) %>% 
    dplyr::ungroup() %>% 
    hydrofabric3D::add_tmp_id() %>% 
    dplyr::group_by(hy_id) %>% 
    dplyr::mutate(
      new_cs_id = 1:dplyr::n()
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(new_cs_id, tmp_id)
  
  # Renumber the transects to have correct CS IDs
  out_transects <- dplyr::left_join(
                      hydrofabric3D::add_tmp_id(out_transects),
                      renumbered_ids,
                      by = "tmp_id"
                    ) %>% 
                      dplyr::select(-cs_id, -tmp_id) %>% 
                      dplyr::select(hy_id, cs_source, 
                                    cs_id = new_cs_id, 
                                    cs_measure, cs_lengthm, 
                                    # sinuosity,
                                    is_extended, 
                                    left_bank_count, right_bank_count, channel_count, bottom_count, 
                                    bottom, left_bank, right_bank, valid_banks, has_relief,
                                    geometry = geom
                      )
  
  # Renumber the cross sections points to have correct CS IDs
  fixed_pts <- dplyr::left_join(
                    hydrofabric3D::add_tmp_id(fixed_pts),
                    renumbered_ids,
                    by = "tmp_id"
                  ) %>% 
                    dplyr::select(-cs_id, -tmp_id) %>% 
                    dplyr::rename(cs_id = new_cs_id)

  # ----------------------------------------------------------------------------------------------------------------
  # ---- Cross section points parquet to S3 ----
  # ----------------------------------------------------------------------------------------------------------------
  
  # classify the cross section points
  fixed_pts <-
    fixed_pts %>% 
    dplyr::mutate(
      X = sf::st_coordinates(.)[,1],
      Y = sf::st_coordinates(.)[,2]
    ) %>%
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
  fixed_pts <- hydrofabric3D::align_banks_and_bottoms(cs_pts = fixed_pts)
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
  # ---- Re upload the updated transects geopackage to S3 ----
  # ----------------------------------------------------------------------------------------------------------------
  updated_path <- gsub(transect_file, paste0("updated_", transect_file), transect_path)
  
  ## Save local and REUPLOAD TRANSECTS to S3 to update for any extended cross sections
  message("Saving updated transects to:\n - filepath: '", updated_path, "'")
  
  # save flowlines to out_path (lynker-spatial/01_transects/transects_<VPU num>.gpkg)
  sf::write_sf(
    out_transects,
    # transect_path
    updated_path
  )
  
  # command to copy transects geopackage to S3
  trans_to_s3 <- paste0("aws s3 cp ", updated_path, " ", S3_TRANSECTS_DIR, transect_file, 
                        ifelse(is.null(AWS_PROFILE), "", paste0(" --profile ", AWS_PROFILE)))
  
  message("Copy VPU ", path_df$vpu[i], " transects to S3:\n - S3 copy command:\n'", 
          trans_to_s3, 
          "'\n==========================")
  
  system(trans_to_s3, intern = TRUE)
  
  ############################################################################## 
  
  # ----------------------------------------------------------------------------------------------------------------
  # ---- Upload the cross section points parquet to S3 ----
  # ----------------------------------------------------------------------------------------------------------------
  
  # name of file and path to save transects gpkg too
  out_file <- paste0("nextgen_", path_df$vpu[i], "_cross_sections.parquet")
  out_path <- paste0(CS_PTS_DIR, out_file)
  
  message("Saving cross section points to:\n - filepath: '", out_path, "'")
  
  # save cross section points as a parquet to out_path (lynker-spatial/02_cs_pts/cs_pts_<VPU num>.parquet)
  arrow::write_parquet(fixed_pts, out_path)
  
  # command to copy cross section points parquet to S3
  copy_cs_pts_to_s3 <- paste0("aws s3 cp ", out_path, " ", S3_CS_PTS_DIR, out_file,
                              ifelse(is.null(AWS_PROFILE), "", paste0(" --profile ", AWS_PROFILE)))
  
  message("Copy VPU ", path_df$vpu[i], " cross sections to S3:\n - S3 copy command:\n'", 
          paste0("aws s3 cp ", out_path, " ", S3_CS_PTS_DIR, out_file,
                 ifelse(is.null(AWS_PROFILE), "", paste0(" --profile ", AWS_PROFILE))), 
          "'\n==========================")
  
  system(copy_cs_pts_to_s3, intern = TRUE)
  
  end <- Sys.time()
  
  message("Finished cross section point generation for VPU ", VPU)
  message("- Completed at: ", end)
  message("==========================")
  
  rm(fixed_pts)
  gc()
  gc()
}

# ###########################################################################################################################################
