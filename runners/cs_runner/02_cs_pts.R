# Generate the flowlines layer for the final cross_sections_<VPU>.gpkg for each VPU
source("runners/cs_runner/config.R")

# # load libraries
# library(hydrofabric3D)
# # library(terrainSliceR)
# library(dplyr)
# library(sf)

# cross section bucket prefix
cs_pts_prefix <- paste0(s3_bucket, version_prefix, "/3D/dem-cross-sections/")
# cs_pts_prefix <- paste0(s3_bucket, "v20/3D/dem-cross-sections/")

# transect bucket prefix
transects_prefix <- paste0(s3_bucket, version_prefix, "/3D/transects/")

# paths to nextgen datasets
nextgen_files <- list.files(nextgen_dir, full.names = FALSE)

# paths to nextgen datasets
transect_files <- list.files(transects_dir, full.names = FALSE)

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
for (i in 1:nrow(path_df)) {
  # i = 15
  
  # nextgen file and full path
  nextgen_file <- path_df$x[i]
  nextgen_path <- paste0(nextgen_dir, nextgen_file)
  
  # model attributes file and full path
  transect_file <- path_df$y[i]
  transect_path <- paste0(transects_dir, transect_file)
  
  # model attributes file and full path
  ref_file <- path_df$ref_file[i]
  ref_path <- paste0(ref_features_dir, "gpkg/", ref_file)
  
  message("Creating VPU ", path_df$vpu[i], 
          " cross section points:\n - flowpaths: '", nextgen_file,
          "'\n - transects: '", transect_file, "'", 
          "\n - waterbodies: '", ref_file, "'"
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

    # Update flowlines and transects to remove flowlines and transects that intersect with reference_features waterbodies
    feature_subsets <- wb_intersects(flines, transects, waterbodies)

    # replace flowlines and transects objects with updated versions in "updated_features"
    flines    <- flines[feature_subsets$valid_flowlines, ]
    transects <- transects[feature_subsets$valid_transects, ]

    # get start time for log messages
    time1 <- Sys.time()

    # get cross section point elevations
    cs_pts <- hydrofabric3D::cross_section_pts(
      cs             = transects,
      points_per_cs  = NULL,
      min_pts_per_cs = 10,
      dem            = DEM_URL
      )

    # try to extend any cross sections that returned cross section points with 
    # identical Z values within a certain threshold ("flat" cross sections)
    cs_pts <- hydrofabric3D::rectify_flat_cs(
      net            = flines,
      cs             = transects,
      cs_pts         = cs_pts, 
      points_per_cs  = NULL,
      min_pts_per_cs = 10,
      dem            = DEM_URL,
      scale          = EXTENSION_PCT,
      threshold      = 0
    )

    # get end time for log messages
    time2 <- Sys.time()
    time_diff <- round(as.numeric(time2 - time1 ), 2)
    
    message("\n\n ---> Cross section point elevations processed in ", time_diff)
    
    # Remove any cross section that has ANY missing (NA) Z values. 
    cs_pts <-
      cs_pts %>%
      # dplyr::filter(hy_id %in% c("wb-849054", "wb-845736")) %>% 
      dplyr::group_by(hy_id, cs_id) %>% 
      dplyr::filter(!any(is.na(Z))) %>% 
      dplyr::ungroup()

    # # check the number of cross sections that were extended
    # cs_pts$is_extended %>% table()

    # extract cross section points that have an "is_extended" value of TRUE
    extended_pts = 
      cs_pts %>% 
      dplyr::filter(is_extended) %>% 
      dplyr::mutate(tmp_id = paste0(hy_id, "_", cs_id)) 

    # extract transects that have a "hy_id" in the "extended_pts" dataset
    update_transects = 
      transects %>% 
      dplyr::mutate(tmp_id = paste0(hy_id, "_", cs_id)) %>%
      dplyr::filter(tmp_id %in% unique(extended_pts$tmp_id))
      
    # if any transects were extended, update the transects dataset, and overwrite local and S3 transects geopackages
    if (nrow(update_transects) > 0) {
      message("Updating ", nrow(update_transects), " transects")

      update_transects =
        update_transects %>%
        # dplyr::filter(hy_id %in% unique(extended_pts$hy_id)) %>% 
        # apply extend_by_percent function to each transect line: 
        hydrofabric3D:::extend_by_percent(
          pct        = EXTENSION_PCT,
          length_col = "cs_lengthm"
        )

    # remove old transects that have "tmp_id" in "extended_pts", and replace with "update_transects"
    out_transects = 
      transects %>% 
      dplyr::mutate(tmp_id = paste0(hy_id, "_", cs_id)) %>%
      dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id)) %>%
      dplyr::bind_rows(
        update_transects
        )  %>% 
      # dplyr::mutate(is_extended = FALSE) %>%
      # dplyr::bind_rows(
      #   dplyr::mutate(update_transects, is_extended = TRUE)
      #   )  %>% 
      dplyr::select(-tmp_id) 

    # mapview::mapview(transects, color = "red") +
    #  mapview::mapview(dplyr::filter(out_transects, is_extended), color = "green") +
    #  mapview::mapview(flines, color = "dodgerblue") 

    ###################################### 
    
    ## Save local and REUPLOAD TRANSECTS to S3 to update for any extended cross sections
    message("Saving updated transects to:\n - filepath: '", transect_path, "'")
  
    # save flowlines to out_path (lynker-spatial/01_transects/transects_<VPU num>.gpkg)
    sf::write_sf(
      out_transects,
      transect_path
    )

    # command to copy transects geopackage to S3
    trans_to_s3 <- paste0("aws s3 cp ", transect_path, " ", transects_prefix, transect_file, 
                            ifelse(is.null(aws_profile), "", paste0(" --profile ", aws_profile)))

    message("Copy VPU ", path_df$vpu[i], " transects to S3:\n - S3 copy command:\n'", 
            trans_to_s3, 
            "'\n==========================")
    
    system(trans_to_s3, intern = TRUE)

    ###################################### 

    }

    # classify the cross section points
    cs_pts <-
      cs_pts %>%
      dplyr::rename(cs_widths = cs_lengthm) %>%
      hydrofabric3D::classify_points() %>%
      dplyr::mutate(
        X = sf::st_coordinates(.)[,1],
        Y = sf::st_coordinates(.)[,2]
      ) %>%
      dplyr::select(
        hy_id, cs_id, pt_id,
        cs_lengthm = cs_widths,
        relative_distance,
        X, Y, Z,
        class
        )

  # Drop point geometries, leaving just X, Y, Z values
  cs_pts <- sf::st_drop_geometry(cs_pts)

  # add Z_source column for source of elevation data
  cs_pts <-
    cs_pts %>%
    dplyr::mutate(
      Z_source = cs_source
      ) %>%
    dplyr::relocate(hy_id, cs_id, pt_id, cs_lengthm, relative_distance, X, Y, Z, Z_source, class)

  ###################################### 

  ###################################### 
  # name of file and path to save transects gpkg too
  out_file <- paste0("nextgen_", path_df$vpu[i], "_cross_sections.parquet")
  out_path <- paste0(cs_pts_dir, out_file)
  
  message("Saving cross section points to:\n - filepath: '", out_path, "'")
  
  # save cross section points as a parquet to out_path (lynker-spatial/02_cs_pts/cs_pts_<VPU num>.parquet)
  arrow::write_parquet(cs_pts, out_path)
  
  # command to copy cross section points parquet to S3
  copy_cs_pts_to_s3 <- paste0("aws s3 cp ", out_path, " ", cs_pts_prefix, out_file,
                          ifelse(is.null(aws_profile), "", paste0(" --profile ", aws_profile)))

  message("Copy VPU ", path_df$vpu[i], " cross sections to S3:\n - S3 copy command:\n'", 
          paste0("aws s3 cp ", out_path, " ", cs_pts_prefix, out_file,
                ifelse(is.null(aws_profile), "", paste0(" --profile ", aws_profile))), 
          "'\n==========================")
  
  system(copy_cs_pts_to_s3, intern = TRUE)

}
