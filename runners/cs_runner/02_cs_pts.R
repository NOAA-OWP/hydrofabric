
###########################################################################################################################################
################################################.     REDO EVERYTHING   #######################################################
###########################################################################################################################################

# Generate the flowlines layer for the final cross_sections_<VPU>.gpkg for each VPU
source("runners/cs_runner/config.R")

# # load libraries
library(hydrofabric3D)
library(dplyr)
library(sf)

# cross section bucket prefix
cs_pts_prefix <- paste0(s3_bucket, version_prefix, "/3D/dem-cross-sections/")
# cs_pts_prefix <- paste0(s3_bucket, "v20/3D/dem-cross-sections/")

# transect bucket prefix
transects_prefix <- paste0(s3_bucket, version_prefix, "/3D/transects/")

# paths to nextgen datasets
nextgen_files <- list.files(nextgen_dir, full.names = FALSE)

# paths to nextgen datasets
transect_files <- list.files(transects_dir, full.names = FALSE)
transect_files <- transect_files[!grepl("updated_", transect_files)]

# string to fill in "cs_source" column in output datasets
cs_source <- "hydrofabric3D"

# reference features dataframe
ref_df <- data.frame(
  vpu      = sapply(strsplit(ref_features, "_", fixed = TRUE), function(i) { i[1] }), 
  ref_file = ref_features
)

# ensure the files are in the same order and matched up by VPU
path_df <- align_files_by_vpu(
  x    = nextgen_files,
  y    = transect_files,
  base = base_dir
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
  nextgen_path <- paste0(nextgen_dir, nextgen_file)
  
  # model attributes file and full path
  transect_file <- path_df$y[i]
  transect_path <- paste0(transects_dir, transect_file)
  
  # model attributes file and full path
  ref_file <- path_df$ref_file[i]
  ref_path <- paste0(ref_features_dir, "gpkg/", ref_file)
  
  # current VPU being processed
  VPU = path_df$vpu[i]
  
  message("Creating VPU ", VPU, 
          " cross section points:\n - flowpaths: '", nextgen_file,
          "'\n - transects: '", transect_file, "'", 
          "\n - waterbodies: '", ref_file, "'",
          "'\n - start time: '", start, "'"
  )
  
  ################### 
  
  # read in transects data
  transects <- sf::read_sf(transect_path)
  
  # read in nextgen data
  flines <- sf::read_sf(nextgen_path, layer = "flowpaths")
  
  # read in waterbodies reference features layer
  waterbodies <- sf::read_sf(ref_path, layer = "waterbodies")
  
  ##### subset flowlines and transects to first 5 features for testing #####
  # flines = dplyr::slice(flines, 1:5) 
  # transects = dplyr::filter(transects, hy_id %in% unique(flines$id))
  ##### #####
  
  # system.time({
  # Update flowlines and transects to remove flowlines and transects that intersect with reference_features waterbodies
  feature_subsets <- wb_intersects(flines, transects, waterbodies)
  # })
  
  # tmp <- transects %>% 
  #   dplyr::slice(1:15000)
  # tmp_flines <- flines %>% 
  #   dplyr::filter(id %in% tmp$hy_id)
  # mapview::mapview(tmp_flines, color = "dodgerblue") + mapview::mapview(tmp, color = "red") 
  
  # replace flowlines and transects objects with updated versions in "updated_features"
  flines    <- flines[feature_subsets$valid_flowlines, ]
  transects <- transects[feature_subsets$valid_transects, ]
  
  rm(waterbodies)
  gc()
  
  #################################################################
  ##### Temporary subsetting to speed local development up ########
  #################################################################
  
  # flines <- 
  #   flines %>% 
  #   dplyr::group_by(order) %>% 
  #   dplyr::slice(1:100) %>% 
  #   dplyr::ungroup()
  # # flines <- 
  # #   flines %>% 
  # #   dplyr::slice(1:2500) 
  # 
  # transects <- 
  #   transects %>% 
  #   # dplyr::filter(hy_id %in% unique(tmp$id)) 
  #   dplyr::filter(hy_id %in% unique(flines$id)) 
  
  #################################################################
  #################################################################
  
  start_cs_pts <- Sys.time()
  message("Extracting cross section points (", start_cs_pts, ")")
  # message("Extracting cross section points (", Sys.time(),")")    
  
  # STEP 1: Extract cs points from DEM
  # system.time({
    
    # get cross section point elevations
    cs_pts <- hydrofabric3D::cross_section_pts(
      cs             = transects,
      points_per_cs  = NULL,
      min_pts_per_cs = 10,
      dem            = DEM_URL
    )
    
  # })
  
  # system.time({
  
  # STEP 2:
  # Remove any cross section that has ANY missing (NA) Z values, and classify the points 
  cs_pts <- 
    cs_pts %>% 
    # dplyr::filter(hy_id %in% c("wb-849054", "wb-845736")) %>% 
    dplyr::group_by(hy_id, cs_id) %>% 
    dplyr::filter(!any(is.na(Z))) %>% 
    dplyr::ungroup() %>% 
    hydrofabric3D::classify_points(pct_of_length_for_relief = 0.01)  
  
  # })
  
  # STEP 3: Try to rectify any no relief and invalid banks cross sections
  system.time({
    # cs_pts <- hydrofabric3D::rectify_flat_cs(
    fixed_pts <- hydrofabric3D::rectify_cs(
      cs_pts         = cs_pts,    # cross section points generated from hydrofabric3D::cross_section_pts()
      net            = flines,    # original flowline network
      transects      = transects, # original transect lines
      points_per_cs  = NULL, 
      min_pts_per_cs = 10, # number of points per cross sections
      dem            = DEM_URL, # DEM to extract points from
      scale          = EXTENSION_PCT, # How far to extend transects if the points need to be rechecked
      pct_of_length_for_relief = PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF, # percent of cross sections length to be needed in relief calculation to consider cross section to "have relief"
      fix_ids = FALSE,
      verbose = TRUE
    )
  })
  # })
  
  # get a summary dataframe and print out details message
  rectify_summary <- hydrofabric3D::rectify_summary(cs_pts, fixed_pts, verbose = TRUE)
  rectify_summary <- 
    rectify_summary %>% 
    dplyr::mutate(
      vpu = VPU
    ) %>% 
    dplyr::relocate(vpu, metric, value)
  
  readr::write_csv(rectify_summary, paste0(META_PATH, "nextgen_", VPU, "_cross_sections_metadata.csv"))
  
  # get the counts of each point type to add this data to the transects dataset
  point_type_counts <- hydrofabric3D::get_point_type_counts(fixed_pts, add = FALSE)
  
  # # check the number of cross sections that were extended
  # fixed_pts$is_extended %>% table()
  message("Subsetting cross section points generated after extending transects...")
  
  # extract cross section points that have an "is_extended" value of TRUE
  extended_pts <-
    fixed_pts %>% 
    dplyr::filter(is_extended) %>% 
    hydrofabric3D::add_tmp_id()
  # dplyr::mutate(tmp_id = paste0(hy_id, "_", cs_id)) 
  
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
      # dplyr::filter(hy_id %in% unique(extended_pts$hy_id)) %>% 
      # apply extend_by_percent function to each transect line: 
      hydrofabric3D:::extend_by_percent(
        pct        = EXTENSION_PCT,
        length_col = "cs_lengthm"
      )
    
    # # Number of transects being updated
    # if(COLLECT_META) {
    # extended_transects_count <- nrow(update_transects)
    # extended_transects_ids   <- length(unique(update_transects$tmp_id))
    # }
    
    # start_uids <- hydrofabric3D:::get_unique_tmp_ids(cs_pts)
    # end_uids <- hydrofabric3D:::get_unique_tmp_ids(fixed_pts)
    # removed_tmp_ids <- start_uids[!start_uids %in% end_uids]
    # transects %>% 
    #   hydrofabric3D::add_tmp_id() %>% 
    #   dplyr::filter(!tmp_id %in% removed_tmp_ids)
  
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
  
  # -------------------------------------------------------------------
  # ---- Re enumerate the transects & cross section points "cs_id" ----
  # -------------------------------------------------------------------
  
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
                      dplyr::select(hy_id, cs_source, cs_id = new_cs_id, 
                                    cs_measure, cs_lengthm, 
                                    # sinuosity,
                                    is_extended, 
                                    left_bank_count, right_bank_count, channel_count, bottom_count, 
                                    bottom, left_bank, right_bank, valid_banks, has_relief,
                                    geometry = geom
                      )
  
  # # fline_lengths <-  sf::st_drop_geometry(flines) %>% 
  # #   dplyr::filter(id %in% out_transects$hy_id) %>% 
  # #   dplyr::mutate(lengthm = lengthkm * 1000) %>% 
  # #   dplyr::select(hy_id = id, lengthm, lengthkm) 
  # tmp <- dplyr::left_join( out_transects, fline_lengths, by = "hy_id") %>% 
  #   dplyr::mutate(ds_distance = (cs_measure * lengthm) / 100) %>% 
  #   dplyr::select(-sinuosity) %>% 
  #   dplyr::relocate(hy_id, cs_id, cs_measure, lengthm, ds_distance, lengthkm) %>% 
  #   dplyr::rename("geometry" = geom)
  
  # Renumber the cross sections points to have correct CS IDs
  fixed_pts <- dplyr::left_join(
                    hydrofabric3D::add_tmp_id(fixed_pts),
                    renumbered_ids,
                    by = "tmp_id"
                  ) %>% 
                    dplyr::select(-cs_id, -tmp_id) %>% 
                    dplyr::rename(cs_id = new_cs_id)
  
  # mapview::mapview(transects, color = "red") +
  #  mapview::mapview(dplyr::filter(out_transects, is_extended), color = "green") +
  #  mapview::mapview(flines, color = "dodgerblue") 
  
  ###################################### 
  
  # ----------------------------------------------------------
  # ---- Cross section points parquet to S3 ----
  # ----------------------------------------------------------
  
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
      Z_source = cs_source
    ) %>%
    dplyr::relocate(hy_id, cs_id, pt_id, cs_lengthm, relative_distance, X, Y, Z, Z_source, 
                    class, point_type, 
                    bottom, left_bank, right_bank, valid_banks, has_relief)
  
  ###################################### 
  
  # ----------------------------------------------------------
  # ---- Re upload the updated transects geopackage to S3 ----
  # ----------------------------------------------------------
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
  trans_to_s3 <- paste0("aws s3 cp ", updated_path, " ", transects_prefix, transect_file, 
                        ifelse(is.null(aws_profile), "", paste0(" --profile ", aws_profile)))
  
  message("Copy VPU ", path_df$vpu[i], " transects to S3:\n - S3 copy command:\n'", 
          trans_to_s3, 
          "'\n==========================")
  
  system(trans_to_s3, intern = TRUE)
  
  ###################################### 
  ###################################### 
  
  # ----------------------------------------------------------
  # ---- Upload the cross section points parquet to S3 ----
  # ----------------------------------------------------------
  
  # name of file and path to save transects gpkg too
  out_file <- paste0("nextgen_", path_df$vpu[i], "_cross_sections.parquet")
  out_path <- paste0(cs_pts_dir, out_file)
  
  message("Saving cross section points to:\n - filepath: '", out_path, "'")
  
  # save cross section points as a parquet to out_path (lynker-spatial/02_cs_pts/cs_pts_<VPU num>.parquet)
  arrow::write_parquet(fixed_pts, out_path)
  
  # command to copy cross section points parquet to S3
  copy_cs_pts_to_s3 <- paste0("aws s3 cp ", out_path, " ", cs_pts_prefix, out_file,
                              ifelse(is.null(aws_profile), "", paste0(" --profile ", aws_profile)))
  
  message("Copy VPU ", path_df$vpu[i], " cross sections to S3:\n - S3 copy command:\n'", 
          paste0("aws s3 cp ", out_path, " ", cs_pts_prefix, out_file,
                 ifelse(is.null(aws_profile), "", paste0(" --profile ", aws_profile))), 
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
# 
# transects
# flines
# ###########################################################################################################################################
# ###########################################################################################################################################
# 
# # Generate the flowlines layer for the final cross_sections_<VPU>.gpkg for each VPU
# source("runners/cs_runner/config.R")
# 
# # # load libraries
# library(hydrofabric3D)
# library(dplyr)
# library(sf)
# 
# # cross section bucket prefix
# cs_pts_prefix <- paste0(s3_bucket, version_prefix, "/3D/dem-cross-sections/")
# # cs_pts_prefix <- paste0(s3_bucket, "v20/3D/dem-cross-sections/")
# 
# # transect bucket prefix
# transects_prefix <- paste0(s3_bucket, version_prefix, "/3D/transects/")
# 
# # paths to nextgen datasets
# nextgen_files <- list.files(nextgen_dir, full.names = FALSE)
# 
# # paths to nextgen datasets
# transect_files <- list.files(transects_dir, full.names = FALSE)
# 
# # string to fill in "cs_source" column in output datasets
# cs_source <- "hydrofabric3D"
# 
# # reference features dataframe
# ref_df <- data.frame(
#   vpu      = sapply(strsplit(ref_features, "_", fixed = TRUE), function(i) { i[1] }), 
#   ref_file = ref_features
# )
# 
# # ensure the files are in the same order and matched up by VPU
# path_df <- align_files_by_vpu(
#               x    = nextgen_files,
#               y    = transect_files,
#               base = base_dir
#             ) %>% 
#           dplyr::left_join(
#             ref_df,
#             by = "vpu"
#             )
# 
# # loop over the nextgen and transect datasets (by VPU) and extract point elevations across points on each transect line,
# # then classify the points, and create a parquet file with hy_id, cs_id, pt_id, X, Y, Z data.
# # Save parquet locally and upload to specified S3 bucket
# for (i in 1:nrow(path_df)) {
# 
#   # i = 8
#   
#   start <- Sys.time()
#   
#   # nextgen file and full path
#   nextgen_file <- path_df$x[i]
#   nextgen_path <- paste0(nextgen_dir, nextgen_file)
#   
#   # model attributes file and full path
#   transect_file <- path_df$y[i]
#   transect_path <- paste0(transects_dir, transect_file)
#   
#   # model attributes file and full path
#   ref_file <- path_df$ref_file[i]
#   ref_path <- paste0(ref_features_dir, "gpkg/", ref_file)
#   
#   # current VPU being processed
#   VPU = path_df$vpu[i]
#   
#   message("Creating VPU ", VPU, 
#           " cross section points:\n - flowpaths: '", nextgen_file,
#           "'\n - transects: '", transect_file, "'", 
#           "\n - waterbodies: '", ref_file, "'",
#           "'\n - start time: '", start, "'"
#           )
# 
#   ################### 
#   
#     # read in transects data
#     transects <- sf::read_sf(transect_path)
#   
#     # read in nextgen data
#     flines <- sf::read_sf(nextgen_path, layer = "flowpaths")
#   
#     # read in waterbodies reference features layer
#     waterbodies <- sf::read_sf(ref_path, layer = "waterbodies")
# 
#     ##### subset flowlines and transects to first 5 features for testing #####
#     # flines = dplyr::slice(flines, 1:5) 
#     # transects = dplyr::filter(transects, hy_id %in% unique(flines$id))
#     ##### #####
#     
#     # system.time({
#       # Update flowlines and transects to remove flowlines and transects that intersect with reference_features waterbodies
#       feature_subsets <- wb_intersects(flines, transects, waterbodies)
#     # })
#     
#     # Collect meta data on features and changes
#     if(COLLECT_META) {
#       
#       fline_count       <- nrow(flines)
#       transect_count    <- nrow(transects)
#       wb_count          <- nrow(waterbodies)
#       
#       fline_wb_count    <- sum(feature_subsets$valid_flowlines)
#       transect_wb_count <- sum(feature_subsets$valid_transects)
#     }
#     
#     # replace flowlines and transects objects with updated versions in "updated_features"
#     flines    <- flines[feature_subsets$valid_flowlines, ]
#     transects <- transects[feature_subsets$valid_transects, ]
#     
    # rm(waterbodies)
    # gc()
#     
#     flines <- 
#       flines %>% 
#       dplyr::group_by(order) %>% 
#       dplyr::slice(1:100) %>% 
#       dplyr::ungroup()
#     # flines <- 
#     #   flines %>% 
#     #   dplyr::slice(1:2500) 
#     
#     transects <- 
#       transects %>% 
#       # dplyr::filter(hy_id %in% unique(tmp$id)) 
#       dplyr::filter(hy_id %in% unique(flines$id)) 
#     
#     
#     start_cs_pts <- Sys.time()
#     
#     message("Extracting cross section points (", start_cs_pts,")")    
#     
#     system.time({
#   
#     # get cross section point elevations
#     cs_pts <- hydrofabric3D::cross_section_pts(
#       cs             = transects,
#       points_per_cs  = NULL,
#       min_pts_per_cs = 10,
#       dem            = DEM_URL
#       )
#     
#     })
#     
#     # Remove any cross section that has ANY missing (NA) Z values. 
#     cs_pts2 <-
#       cs_pts %>% 
#       # dplyr::filter(hy_id %in% c("wb-849054", "wb-845736")) %>% 
#       dplyr::group_by(hy_id, cs_id) %>% 
#       dplyr::filter(!any(is.na(Z))) %>% 
#       dplyr::ungroup() 
#     
#     # Stash meta data about the points 
#     pts_meta <- 
#       cs_pts2 %>% 
#       sf::st_drop_geometry() %>% 
#       dplyr::select(hy_id, cs_id, pt_id, cs_measure, is_extended)
#     
#     # classify the points
#     cs_pts2 <- hydrofabric3D::classify_points(cs_pts2, pct_of_length_for_relief = 0.01)
#     
#     system.time({
#       # cs_pts <- hydrofabric3D::rectify_flat_cs(
#       fixed_pts <- hydrofabric3D::rectify_cs(
#         cs_pts         = cs_pts2,    # cross section points generated from hydrofabric3D::cross_section_pts()
#         net            = flines,    # original flowline network
#         transects      = transects, # original transect lines
#         points_per_cs  = NULL, 
#         min_pts_per_cs = 10, # number of points per cross sections
#         dem            = DEM_URL, # DEM to extract points from
#         scale          = EXTENSION_PCT, # How far to extend transects if the points need to be rechecked
#         pct_of_length_for_relief = 0.01, # percent of cross sections length to be needed in relief calculation to consider cross section to "have relief"
#         fix_ids = FALSE,
#         verbose = TRUE
#       )
#     })
#     
#     rectify_summary <- hydrofabric3D::rectify_summary(cs_pts2, fixed_pts)
#     
#     end_cs_pts <- Sys.time()
#     message("\n ---> Completed extraction of cross section points (", end_cs_pts,")")
#     
#     if(COLLECT_META) {
#       start_cs_pts_count <- nrow(cs_pts)
#     }
#     
#     # cs_pts_time <- round(as.numeric(end_cs_pts - start_cs_pts ), 2)
#     # message("\n\n ---> Cross section point elevations processed in ", cs_pts_time)
#     
#     start_rectify <- Sys.time()
#     message("Rectifying cross section points (", start_rectify,")")
#     
#     # collect the hy_ids and number of stream orders in cs_pts
#     if(COLLECT_META) {
#     
#       cs_pts_ids       <- unique(cs_pts$hy_id)
#       start_cs_pts_ids <- length(cs_pts_ids)
#       
#       start_order_count <- 
#         flines %>% 
#         sf::st_drop_geometry() %>% 
#         dplyr::filter(id %in% cs_pts_ids) %>% 
#         dplyr::group_by(order) %>% 
#         dplyr::count() %>% 
#         tidyr::pivot_wider(names_from = order,  
#                          names_glue   = "start_order_{order}", 
#                          values_from  = n
#                          ) %>% 
#         dplyr::ungroup()
#     }
#     
#     # Remove any cross section that has ANY missing (NA) Z values. 
#     cs_pts <- 
#       cs_pts %>% 
#       # dplyr::filter(hy_id %in% c("wb-849054", "wb-845736")) %>% 
#       dplyr::group_by(hy_id, cs_id) %>% 
#       dplyr::filter(!any(is.na(Z))) %>% 
#       dplyr::ungroup()
#     # try to extend any cross sections that returned cross section points with 
#     # identical Z values within a certain threshold ("flat" cross sections)
#     
#     system.time({
#     # cs_pts <- hydrofabric3D::rectify_flat_cs(
#       fixed_pts <- hydrofabric3D::rectify_flat_cs(
#                   cs_pts         = cs_pts,    # cross section points generated from hydrofabric3D::cross_section_pts()
#                   net            = flines,    # original flowline network
#                   cs             = transects, # original transect lines
#                   points_per_cs  = NULL, 
#                   min_pts_per_cs = 10, # number of points per cross sections
#                   dem            = DEM_URL, # DEM to extract points from
#                   scale          = EXTENSION_PCT, # How far to extend transects if the points need to be rechecked
#                   threshold      = 1,    # 1 meter from bottom
#                   pct_threshold  = 0.99, # rectify if 99% points are within 1 meter from the bottom
#                   fix_ids        = FALSE
#                   )
#     })
#     
#     
#     end_rectify <- Sys.time()
#     rectify_time <- round(as.numeric(end_rectify - start_rectify ), 2)
#     
#     message("\n ---> Completed rectifying cross section points (", end_rectify,")")
#     
#     if(COLLECT_META) {
#       rectify_cs_pts_count <- nrow(fixed_pts)
#       # collect the hy_ids and number of stream orders in the RECTIFIED cs_pts
#       rectify_cs_pts_ids      <- unique(fixed_pts$hy_id)
#       rectify_cs_pts_id_count <- length(rectify_cs_pts_ids)
# 
#       rectify_order_count <- 
#         flines %>% 
#         sf::st_drop_geometry() %>% 
#         dplyr::filter(id %in% rectify_cs_pts_ids) %>% 
#         dplyr::group_by(order) %>% 
#         dplyr::count() %>% 
#         tidyr::pivot_wider(names_from = order,  
#                            names_glue   = "rectify_order_{order}", 
#                            values_from  = n
#         ) %>% 
#         dplyr::ungroup()
#     }
#     
#     rm(cs_pts)
#     gc()
#     
#     message("\n\n ---> Cross section points rectified in ", rectify_time, " (seconds?) ")
#     
#     # Remove any cross section that has ANY missing (NA) Z values. 
#     fixed_pts <- 
#       fixed_pts %>% 
#       # dplyr::filter(hy_id %in% c("wb-849054", "wb-845736")) %>% 
#       # dplyr::group_by(hy_id, cs_id) %>% 
#       # dplyr::filter(!any(is.na(Z))) %>% 
#       # dplyr::ungroup() %>% 
#       hydrofabric3D::add_tmp_id()
#     
#     # Number of cross section points after removing any cross sections that contain any NA Z values
#     if(COLLECT_META) {
#       cs_pts_na_removed_count <- nrow(fixed_pts)
#     }
#     
#     # Stash meta data about the points 
#     pts_meta <- 
#       fixed_pts %>% 
#       sf::st_drop_geometry() %>% 
#       dplyr::select(hy_id, cs_id, pt_id, cs_measure, is_extended)
#     
#     message("Classifying cross section points...")
#     
#     # # Classify points
#     # fixed_pts <- hydrofabric3D::classify_points(fixed_pts)
#     
#     # add meta data back to the points
#     fixed_pts <- 
#       fixed_pts %>% 
#       dplyr::select(-is_extended) %>% 
#       dplyr::left_join(
#         pts_meta,
#         # dplyr::select(pts_meta, -is_extended),
#         by = c("hy_id", "cs_id", "pt_id")
#         # dplyr::select(pts_meta, hy_id, cs_id, pt_id, cs_measure, is_extended)
#       )
#     
#     message("Gathering count of point types per cross section...")
#     
#     # get the counts of each point type to add this data to the transects dataset
#     point_type_counts <- hydrofabric3D::get_point_type_counts(fixed_pts, add = FALSE)
#     
#     # # check the number of cross sections that were extended
#     # fixed_pts$is_extended %>% table()
#     message("Subsetting cross section points generated after extending transects...")
#     
#     # extract cross section points that have an "is_extended" value of TRUE
#     extended_pts <-
#       fixed_pts %>% 
#       dplyr::filter(is_extended) %>% 
#       hydrofabric3D::add_tmp_id()
#       # dplyr::mutate(tmp_id = paste0(hy_id, "_", cs_id)) 
#     
#     # extract transects that have a "hy_id" in the "extended_pts" dataset
#     update_transects <-
#       transects %>% 
#       hydrofabric3D::add_tmp_id() %>% 
#       dplyr::filter(tmp_id %in% unique(extended_pts$tmp_id))
#     
#     # Number of cross section points generated from extending transects and number of tmpIDs
#     if(COLLECT_META) {
#       extended_pts_count <- nrow(extended_pts)
#       extended_pts_ids   <- length(unique(extended_pts$tmp_id))
#       extended_transects_count <- nrow(update_transects)
#       extended_transects_ids   <- length(unique(update_transects$tmp_id))
#     }
#     
#     # if any transects were extended, update the transects dataset, and overwrite local and S3 transects geopackages
#     if (nrow(update_transects) > 0) {
#       message("Updating ", nrow(update_transects), " transects")
#     
#       update_transects <-
#         update_transects %>%
#         # dplyr::filter(hy_id %in% unique(extended_pts$hy_id)) %>% 
#         # apply extend_by_percent function to each transect line: 
#         hydrofabric3D:::extend_by_percent(
#           pct        = EXTENSION_PCT,
#           length_col = "cs_lengthm"
#         )
#       
#       # # Number of transects being updated
#       # if(COLLECT_META) {
#         # extended_transects_count <- nrow(update_transects)
#         # extended_transects_ids   <- length(unique(update_transects$tmp_id))
#       # }
#       
#       # Filter down to ONLY points that were finalized and rectified from rectify_cs_pts()
#       # remove old transects that have "tmp_id" in "extended_pts" (transects that were unchanged and are "good_to_go")
#       # and then replace with old transects with the "update_transects" 
#       out_transects <- 
#         transects %>% 
#         hydrofabric3D::add_tmp_id() %>% 
#         # dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id)) %>%
#         # dplyr::filter(!tmp_id %in% )
#         dplyr::filter(tmp_id %in% unique(hydrofabric3D::add_tmp_id(fixed_pts)$tmp_id)) %>% 
#         dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id)) %>%
#         dplyr::bind_rows(
#           dplyr::mutate(
#             update_transects,
#             is_extended = TRUE
#             )
#           )  
#       
#         # dplyr::mutate(is_extended = FALSE) %>%
#         # dplyr::bind_rows(
#         #   dplyr::mutate(update_transects, is_extended = TRUE)
#         #   )  %>% 
#         # dplyr::select(-tmp_id) 
#       
#     } else {
#       
#       out_transects <- 
#         transects %>% 
#         hydrofabric3D::add_tmp_id() %>% 
#         dplyr::filter(tmp_id %in% unique(hydrofabric3D::add_tmp_id(fixed_pts)$tmp_id)) %>% 
#         dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id))
#     }
#     
#       # Number of final output transects and the number of unique tmpIDs (hy_id/cs_id , i.e. cross sections)
#       if(COLLECT_META) {
#         output_transects_count <- nrow(out_transects)
#         output_transects_ids   <- length(unique(out_transects$tmp_id))
#       }
#     
#       # finalize new transects
#       out_transects <- 
#         out_transects %>% 
#         dplyr::left_join(
#           point_type_counts,
#           by = c("hy_id", "cs_id")
#         ) %>% 
#         dplyr::select(hy_id, cs_source, cs_id, cs_measure, cs_lengthm, 
#                       # sinuosity, 
#                       is_extended, 
#                       left_bank_count, right_bank_count, channel_count, bottom_count,
#                       geom)
#       
#       # -------------------------------------------------------------------
#       # ---- Re enumerate the transects & cross section points "cs_id" ----
#       # -------------------------------------------------------------------
#       
#       # make a dataframe that has a new_cs_id column that has 
#       # the cs_id renumbered to fill in any missing IDs,
#       # so each hy_id has cs_ids that go from 1 - number of cross sections on hy_id
#       # The dataframe below will be used to join the "new_cs_id" with 
#       # the original "cs_ids" in the final cross section POINTS and UPDATED TRANSECTS output datasets
#       renumbered_ids <-
#         fixed_pts %>% 
#         sf::st_drop_geometry() %>% 
#         # dplyr::filter(hy_id %in% c("wb-2402800", "wb-2398282", "wb-2400351")) %>%
#         dplyr::select(hy_id, cs_id, pt_id, cs_measure) %>% 
#         dplyr::group_by(hy_id, cs_id) %>% 
#         dplyr::slice(1) %>% 
#         dplyr::ungroup() %>% 
#         hydrofabric3D::add_tmp_id() %>% 
#         dplyr::group_by(hy_id) %>% 
#         dplyr::mutate(
#           new_cs_id = 1:dplyr::n()
#         ) %>% 
#         dplyr::ungroup() %>% 
#         dplyr::select(new_cs_id, tmp_id)
#       
#       # Renumber the transects to have correct CS IDs
#       out_transects <- dplyr::left_join(
#                           hydrofabric3D::add_tmp_id(out_transects),
#                           renumbered_ids,
#                           by = "tmp_id"
#                         ) %>% 
#                         dplyr::select(-cs_id, -tmp_id) %>% 
#                         dplyr::select(hy_id, cs_source, cs_id = new_cs_id, 
#                                       cs_measure, cs_lengthm, 
#                                       # sinuosity,
#                                       is_extended, 
#                                       left_bank_count, right_bank_count, channel_count, bottom_count,
#                                       geometry = geom
#                                       )
#       
#       # # fline_lengths <-  sf::st_drop_geometry(flines) %>% 
#       # #   dplyr::filter(id %in% out_transects$hy_id) %>% 
#       # #   dplyr::mutate(lengthm = lengthkm * 1000) %>% 
#       # #   dplyr::select(hy_id = id, lengthm, lengthkm) 
#       # tmp <- dplyr::left_join( out_transects, fline_lengths, by = "hy_id") %>% 
#       #   dplyr::mutate(ds_distance = (cs_measure * lengthm) / 100) %>% 
#       #   dplyr::select(-sinuosity) %>% 
#       #   dplyr::relocate(hy_id, cs_id, cs_measure, lengthm, ds_distance, lengthkm) %>% 
#       #   dplyr::rename("geometry" = geom)
#       
#       # Renumber the cross sections points to have correct CS IDs
#       fixed_pts <- dplyr::left_join(
#                             hydrofabric3D::add_tmp_id(fixed_pts),
#                             renumbered_ids,
#                             by = "tmp_id"
#                           ) %>% 
#                         dplyr::select(-cs_id, -tmp_id) %>% 
#                         dplyr::rename(cs_id = new_cs_id)
#       
#     # mapview::mapview(transects, color = "red") +
#     #  mapview::mapview(dplyr::filter(out_transects, is_extended), color = "green") +
#     #  mapview::mapview(flines, color = "dodgerblue") 
# 
#     ###################################### 
# 
#     # ----------------------------------------------------------
#     # ---- Cross section points parquet to S3 ----
#     # ----------------------------------------------------------
#     
#     # classify the cross section points
#     fixed_pts <- 
#       fixed_pts %>% 
#       dplyr::mutate(
#         X = sf::st_coordinates(.)[,1],
#         Y = sf::st_coordinates(.)[,2]
#       ) %>%
#       sf::st_drop_geometry() %>% 
#       dplyr::select(
#         hy_id, cs_id, pt_id,
#         cs_lengthm,
#         relative_distance,
#         X, Y, Z,
#         class, point_type
#         )
# 
#   # Drop point geometries, leaving just X, Y, Z values
#   fixed_pts <- sf::st_drop_geometry(fixed_pts)
# 
#   # add Z_source column for source of elevation data
#   fixed_pts <-
#     fixed_pts %>%
#     dplyr::mutate(
#       Z_source = cs_source
#       ) %>%
#     dplyr::relocate(hy_id, cs_id, pt_id, cs_lengthm, relative_distance, X, Y, Z, Z_source, class)
# 
#   # Number of final output transects and the number of unique tmpIDs (hy_id/cs_id , i.e. cross sections)
#   if(COLLECT_META) {
#     output_cs_pts_count      <- nrow(fixed_pts)
#     output_cs_pts_ids        <- length(unique(hydrofabric3D::add_tmp_id(fixed_pts)$tmp_id))
#     dropped_transects_count  <- transect_count - output_transects_count
#     }
#   
#   ###################################### 
#   
#   # ----------------------------------------------------------
#   # ---- Re upload the updated transects geopackage to S3 ----
#   # ----------------------------------------------------------
#   updated_path <- gsub(transect_file, paste0("updated_", transect_file), transect_path)
#   
#   ## Save local and REUPLOAD TRANSECTS to S3 to update for any extended cross sections
#   message("Saving updated transects to:\n - filepath: '", updated_path, "'")
#   
#   # save flowlines to out_path (lynker-spatial/01_transects/transects_<VPU num>.gpkg)
#   sf::write_sf(
#     out_transects,
#     # transect_path
#     updated_path
#   )
#  
#   # command to copy transects geopackage to S3
#   trans_to_s3 <- paste0("aws s3 cp ", updated_path, " ", transects_prefix, transect_file, 
#                         ifelse(is.null(aws_profile), "", paste0(" --profile ", aws_profile)))
#   
#   message("Copy VPU ", path_df$vpu[i], " transects to S3:\n - S3 copy command:\n'", 
#           trans_to_s3, 
#           "'\n==========================")
#   
#   system(trans_to_s3, intern = TRUE)
#   
#   ###################################### 
#   ###################################### 
#   
#   # ----------------------------------------------------------
#   # ---- Upload the cross section points parquet to S3 ----
#   # ----------------------------------------------------------
#   
#   # name of file and path to save transects gpkg too
#   out_file <- paste0("nextgen_", path_df$vpu[i], "_cross_sections.parquet")
#   out_path <- paste0(cs_pts_dir, out_file)
#   
#   message("Saving cross section points to:\n - filepath: '", out_path, "'")
#   
#   # save cross section points as a parquet to out_path (lynker-spatial/02_cs_pts/cs_pts_<VPU num>.parquet)
#   arrow::write_parquet(fixed_pts, out_path)
#   
#   # command to copy cross section points parquet to S3
#   copy_cs_pts_to_s3 <- paste0("aws s3 cp ", out_path, " ", cs_pts_prefix, out_file,
#                           ifelse(is.null(aws_profile), "", paste0(" --profile ", aws_profile)))
# 
#   message("Copy VPU ", path_df$vpu[i], " cross sections to S3:\n - S3 copy command:\n'", 
#           paste0("aws s3 cp ", out_path, " ", cs_pts_prefix, out_file,
#                 ifelse(is.null(aws_profile), "", paste0(" --profile ", aws_profile))), 
#           "'\n==========================")
#   
#   system(copy_cs_pts_to_s3, intern = TRUE)
#   
#   end <- Sys.time()
#   
#   message("Finished cross section point generation for VPU ", VPU)
#   message("- Completed at: ", end)
#   message("==========================")
#   
#   if(COLLECT_META) {
#     
#     meta_df <- data.frame(
#       vpu                 = VPU,
#       start               = as.character(start),
#       end                 = as.character(end),
#       start_cs_pts        = as.character(start_cs_pts),
#       end_cs_pts          = as.character(end_cs_pts),
#       start_rectify       = as.character(start_rectify),
#       end_rectify         = as.character(end_rectify),
#       fline_count         = fline_count,
#       transect_count      = transect_count,
#       wb_count            = wb_count,
#       fline_wb_count      = fline_wb_count,
#       transect_wb_count   = transect_wb_count,
#       start_cs_pts_count  = start_cs_pts_count,
#       start_cs_pts_ids    = start_cs_pts_ids,
#       rectify_cs_pts_count       = rectify_cs_pts_count,
#       rectify_cs_pts_ids         = rectify_cs_pts_id_count,
#       extended_transects_count   = extended_transects_count,
#       extended_transects_ids     = extended_transects_ids,
#       dropped_transects          = dropped_transects_count,
#       output_transects_count     = output_transects_count,
#       output_cs_pts_count        = output_cs_pts_count,
#       output_transects_ids       = output_transects_ids,
#       output_cs_pts_ids          = output_cs_pts_ids
#     )
#     
#     order_df <- cbind(data.frame(vpu = VPU), start_order_count, rectify_order_count)
#     
#     readr::write_csv(meta_df, paste0(META_PATH, "nextgen_", VPU, "_cross_sections_metadata.csv"))
#     readr::write_csv(order_df, paste0(META_PATH, "nextgen_", VPU, "_cross_sections_streamorder.csv"))
#   }
#   
  # rm(fixed_pts)
  # gc()
  # gc()
#   
#   }
# 
