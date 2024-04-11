# ----------------------------------------------------------------------------------------------------------------
# ---- data paths -----
# ----------------------------------------------------------------------------------------------------------------
library(dplyr)
library(hydrofabric3D)
library(sf)
library(patchwork)

# Generate the flowlines layer for the final cross_sections_<VPU>.gpkg for each VPU
source("runners/cs_runner/config.R")

# cross section bucket prefix
CS_ML_PTS_S3_PREFIX <- paste0(s3_bucket, version_prefix, "/3D/cross-sections/")
# cs_pts_prefix <- paste0(s3_bucket, "v20/3D/dem-cross-sections/")

ML_OUTPUTS_PATH <- list.files(ML_OUTPUTS_DIR, full.names = TRUE)

# paths to nextgen datasets
nextgen_files <- list.files(nextgen_dir, full.names = FALSE)

# paths to nextgen datasets
cs_files <- list.files(cs_pts_dir, full.names = FALSE)

# ensure the files are in the same order and matched up by VPU
path_df <- align_files_by_vpu(
  x    = nextgen_files,
  y    = cs_files,
  base = base_dir
) 
# dplyr::left_join(
#   ref_df,
#   by = "vpu"
# )

# ML Outputs
ml_output <- arrow::read_parquet(ML_OUTPUTS_PATH)

# cs_ml_data_path <- "/Users/anguswatters/Desktop/cs_pts_for_ml_tests/nextgen_06_cross_sections_for_ml.parquet"
# ml_output_path <- "/Users/anguswatters/Desktop/lynker-spatial/ml-outputs/channel_ml_outputs.parquet"
# conus_network_path <- 's3://lynker-spatial/v20.1/conus_net.parquet'

# for (i in 1:5000) {
#   start <- Sys.time()
#   message(i, " - time: '", start, "'")
# }

# loop over the nextgen and transect datasets (by VPU) and extract point elevations across points on each transect line,
# then classify the points, and create a parquet file with hy_id, cs_id, pt_id, X, Y, Z data.
# Save parquet locally and upload to specified S3 bucket
for (i in 1:nrow(path_df)) {
  
  # i = 8
  
  start <- round(Sys.time())
  
  
  # nextgen file and full path
  nextgen_file <- path_df$x[i]
  nextgen_path <- paste0(nextgen_dir, nextgen_file)
  
  # model attributes file and full path
  cs_file      <- path_df$y[i]
  cs_pts_path  <- paste0(cs_pts_dir, cs_file)
  
  # # model attributes file and full path
  # ref_file <- path_df$ref_file[i]
  # ref_path <- paste0(ref_features_dir, "gpkg/", ref_file)
  
  # current VPU being processed
  VPU = path_df$vpu[i]
  
  message("Augmenting DEM cross sections with ML estimated widths/depths: ", VPU, 
          " cross section points:",
          "'\n - cross section points: '", cs_file, "'", 
          "'\n - ML estimated widths/depths: '", ML_OUTPUTS_FILE, "'", 
          # "'\n - ML estimated widths/depths: '", ML_OUTPUTS_URI, "'", 
          "\n - CONUS network file: '", CONUS_NETWORK_URI, "'",
          "\n - flowpaths: '", nextgen_file,
          # "\n - waterbodies: '", ref_file, "'",
          "'\n - start time: '", start, "'"
          )
  
  # ----------------------------------------------------------------------------------------------------------------
  # ---- Read in data -----
  # ----------------------------------------------------------------------------------------------------------------
  message("Loading data...")
  
  # CONUS network parquet
  net <- 
    CONUS_NETWORK_URI %>% 
    arrow::open_dataset() %>%
    dplyr::filter(vpu == VPU) %>%
    dplyr::collect()
  
  # Cross section points parquet
  cs_pts    <- arrow::read_parquet(cs_pts_path)
  
  # TODO: Not needed 
  # fline_net <- sf::read_sf(nextgen_path, layer = "flowpaths")
  
  # ----------------------------------------------------------------------------------------------------------------
  # ---- Extract the max hydroseq "hy_id" for each flowline in the CONUS network parquet -----
  # Use this to join "hy_id" to "hf_id" in ML outputs data
  # ----------------------------------------------------------------------------------------------------------------
  
  net_subset <- 
    net %>% 
    dplyr::select(id, hf_id, hf_hydroseq) %>% 
    # dplyr::filter(id %in% unique(cs$hy_id)) %>%
    dplyr::filter(stringr::str_detect(id, "^wb-") & !is.na(id)) %>% 
    dplyr::group_by(id) %>% 
    dplyr::slice_max(hf_hydroseq, with_ties = FALSE) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(hy_id = id, 
                  hf_id,
                  hf_hydroseq
    ) 
  
  number_hyids_in_cs_pts    <- length(unique(cs_pts$hy_id ))
  number_hyids_in_conus_net <- length(unique(net_subset$hy_id))
  
  message("VPU ", VPU, ": ", 
          "\n - Number of CONUS network hy_ids: '", number_hyids_in_conus_net, "'",
          "\n - Number of cross section points hy_ids: '", number_hyids_in_cs_pts, "'"
  )
  # ----------------------------------------------------------------------------------------------------------------
  # ---- Find the stream orders for the flowlines in the network -----
  # ----------------------------------------------------------------------------------------------------------------
  
  # stream_order <- 
  #   fline_net %>% 
  #   sf::st_drop_geometry() %>% 
  #   dplyr::select(hy_id = id, 
  #                 order, 
  #                 # hydroseq, 
  #                 tot_drainage_areasqkm
  #   ) %>% 
  #   dplyr::filter(hy_id %in% net_subset$hy_id) 
  
  # ----------------------------------------------------------------------------------------------------------------
  # ---- Subset ML data to specific VPU and add "hy_id" column to ML data -----
  # ----------------------------------------------------------------------------------------------------------------
  
  
  # # ML Outputs
  # ml_output <- arrow::read_parquet(ML_OUTPUTS_PATH)
  
  # Join hy_id onto the ML outputs and then remove the rows WITHOUT matches in hy_id
  #  this should give us a (nearly) one-to-one cross walk between "hy_id" in the cross section points 
  # and "hf_id" in the ML outputs dataset
  ml_subset <- 
    ml_output %>% 
    dplyr::left_join(
      net_subset,
      by = "hf_id"
    ) %>% 
    dplyr::filter(!is.na(hy_id)) %>%
    dplyr::relocate(hy_id, hf_id) 
  
  message(round(Sys.time()), " - Joining ML width/depths estimates to cross section points...")

  # join the ML outputs data to the cross section points
  cs_pts <- 
    cs_pts %>% 
    dplyr::left_join(
      dplyr::select(ml_subset, 
                    hy_id, 
                    hf_id, 
                    owp_tw_inchan, 
                    owp_y_inchan, 
                    owp_tw_bf, 
                    owp_y_bf, 
                    owp_dingman_r),
      by = "hy_id"
    ) 
  
  message(round(Sys.time()), " - Replacing any negative width/depth estimates with cross section bottom lengths...")
  
  cs_bottom_lengths <- hydrofabric3D::get_cs_bottom_length(cs_pts)
  
  # TODO: for now we replace any negative TW values with the length of the bottom of the cross section
  # TODO: This method + the negative model output values both need to be looked into (04/05/2024)
  cs_pts <-
    cs_pts %>% 
    dplyr::left_join(
      cs_bottom_lengths, 
      by = c("hy_id", "cs_id")
    ) %>% 
    dplyr::mutate(
      owp_tw_inchan = dplyr::case_when(
        owp_tw_inchan <= 0 ~ bottom_length,
        TRUE               ~ owp_tw_inchan
      ),
      owp_tw_bf = dplyr::case_when(
        owp_tw_bf <= 0     ~ bottom_length,
        TRUE               ~ owp_tw_bf
      )
    ) %>% 
    dplyr::select(-bottom_length)
  
  # cs_pts %>% dplyr::filter(owp_tw_inchan <= 0 | owp_tw_bf <= 0)
  # dplyr::left_join( stream_order, by = "hy_id") 
  
  # extract any cross sections that didn't get matched with a "hf_id" and (or?) no ML data
  # TODO: look at this stuff with Arash (04/09/2024)
  missing_cs <- 
    cs_pts %>% 
    dplyr::filter(is.na(hf_id) | 
                  is.na(owp_tw_inchan) | is.na(owp_y_inchan) | 
                  is.na(owp_tw_bf) | is.na(owp_y_bf) |
                  is.na(owp_dingman_r)) %>% 
    hydrofabric3D::add_tmp_id()
  
  # TODO: Delete this, but time being keeping this to inspect mismatch in between "hy_id" and "hf_id" 
  readr::write_csv(
    dplyr::select(missing_cs, -tmp_id), 
    paste0(META_PATH, "nextgen_", path_df$vpu[i], "_cross_sections_missing_hf_ids.csv")
    )

  # Split the cross sections into 2 groups:
  # - "Inchannel cs" group are points with BOTH valid banks AND relief --> These get the INCHANNEL TW and Y values from the ML model
  # - "Bankful cs" group are points WITHOUT valid banks OR any relief  --> These get the BANKFUL TW and Y values from the ML model
  inchannel_cs <-
    cs_pts %>% 
    hydrofabric3D::add_tmp_id() %>% 
    dplyr::filter(!tmp_id %in% unique(missing_cs$tmp_id)) %>%  # NOTE: makes sure to remove any of the "missing" cross sections without data
    dplyr::select(-tmp_id) %>% 
    dplyr::filter(valid_banks & has_relief) %>% 
    # NOTE: temporarily rename the top widths, depths, and dingman's R columns so they 
    # work nicely with the "hydrofabric3D::add_cs_bathymetry()" function which takes a dataframe of cross section points 
    # with "TW", "DEPTH", and "DINGMAN_R" columns for each cross section
    dplyr::rename(        
      TW        = owp_tw_inchan,    
      DEPTH     = owp_y_inchan,
      DINGMAN_R = owp_dingman_r
    )
  
  bankful_cs   <- 
    cs_pts %>% 
    hydrofabric3D::add_tmp_id() %>% 
    dplyr::filter(!tmp_id %in% unique(missing_cs$tmp_id)) %>% 
    dplyr::select(-tmp_id) %>% 
    dplyr::filter(!valid_banks | !has_relief) %>% 
    dplyr::rename(
      TW        = owp_tw_bf, 
      DEPTH     = owp_y_bf,
      DINGMAN_R = owp_dingman_r
    )
  
  # sanity check that all rows are accounted for after splitting up data
  split_kept_all_rows <- nrow(cs_pts) == nrow(bankful_cs) + nrow(inchannel_cs) + nrow(missing_cs)
  # split_kept_all_rows <- nrow(cs_pts) == nrow(bankful_cs) + nrow(inchannel_cs)
  
  if (!split_kept_all_rows) {
    warning(paste0("When splitting cross section points into 'bankful' and 'inchannel' groups,", 
    "\nsome points were not put in either group.", 
    "\nLikely due to 'valid_banks' and/or 'has_relief' columns have either missing ", 
    "values or contain values other than TRUE/FALSE")
    )
  }
  message(round(Sys.time()), " - Adding cross section bathymetry using inchannel widths/depths estimates...")

  # Add bathymetry using "inchannel" estimates
  inchannel_cs <- hydrofabric3D::add_cs_bathymetry(
      cross_section_pts = inchannel_cs
    )
  
  message(round(Sys.time()), " - Adding cross section bathymetry using bankful widths/depths estimates...")

  # Add bathymetry using "bankful" estimates
  bankful_cs <- hydrofabric3D::add_cs_bathymetry(
    cross_section_pts = bankful_cs
  )
    
  # combine the inchannel and bankful cross section points back together, fill out missing values and reclassify the points
  final_cs <- dplyr::bind_rows(
                dplyr::select(
                  inchannel_cs,
                  -hf_id, -TW, -DEPTH, -DINGMAN_R, -dplyr::starts_with("owp"), -is_dem_point
                ),
                dplyr::select(
                  bankful_cs,
                  -hf_id, -TW, -DEPTH, -DINGMAN_R, -dplyr::starts_with("owp"), -is_dem_point
                ),
                dplyr::select(
                  missing_cs,
                  -hf_id, -dplyr::starts_with("owp"), -tmp_id
                )
              ) %>%
    dplyr::group_by(hy_id, cs_id) %>%
    tidyr::fill(
      c(cs_lengthm, Z_source)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      -point_type,
      -class,
      -bottom, -left_bank, -right_bank,
      -has_relief, -valid_banks
    )
  
  message(round(Sys.time()), " - Reclassifying cross section point types...")
  
  # reclassify
  final_cs <- hydrofabric3D::classify_points(final_cs)

  # final_uids <- final_cs %>% hydrofabric3D::get_unique_tmp_ids()
  # random_uids <- sample(x=final_uids, size=12)
  # cs_subset <-  dplyr::filter(hydrofabric3D::add_tmp_id(final_cs), 
                        # tmp_id %in% random_uids) 
  # hydrofabric3D::classify_points(cs_subset) %>%  hydrofabric3D::plot_cs_pts(color = "point_type")
  
  starting_uids <- hydrofabric3D::get_unique_tmp_ids(cs_pts)
  ending_uids   <- hydrofabric3D::get_unique_tmp_ids(final_cs)

  has_same_number_of_uids          <- length(starting_uids) == length(ending_uids)
  all_starting_uids_in_ending_uids <- all(starting_uids %in% ending_uids)
  all_ending_uids_in_starting_uids <- all(ending_uids %in% starting_uids)

  # throw some warnings if:
  # - the number of uids in the input is different from the output
  # - there are missing hy_id/cs_id
  if (!has_same_number_of_uids) {
    warning(paste0("The number of unique hy_id/cs_id is different in the ",
                   "starting cross section points from the final cross section points:",
                   "\n > Starting number of unique hy_id/cs_id: ", length(starting_uids),
                   "\n > Ending number of unique hy_id/cs_id: ", length(ending_uids)
                   ))
  }

  if (!all_starting_uids_in_ending_uids) {
    number_uids_not_in_ending_uids <- length(starting_uids[!starting_uids %in% ending_uids])

    # starting_uids %in% ending_uids
    warning(
      paste0("Missing hy_id/cs_id in output that are in the starting input cross section points: ",
            "\n > Number of hy_id/cs_id missing: ", number_uids_not_in_ending_uids
             )
      )

    # warning(paste0(number_uids_not_in_ending_uids, " hy_id/cs_id from the input cross section points ",
    #          "is missing from the output cross section points"))

  }
  
  # ----------------------------------------------------------------------------------------------------------------
  # ---- Upload the cross section points parquet to S3 ----
  # ----------------------------------------------------------------------------------------------------------------

  # name of file and path to save transects gpkg too
  out_file <- paste0("nextgen_", path_df$vpu[i], "_cross_sections.parquet")
  out_path <- paste0(final_dir, out_file)
  
  message(round(Sys.time()), " - Saving ML augmented cross section points to:\n - filepath: '", out_path, "'")
  
  # save cross section points as a parquet to out_path (lynker-spatial/02_cs_pts/cs_pts_<VPU num>.parquet)
  arrow::write_parquet(final_cs, out_path)
  
  s3_save_uri <- paste0(CS_ML_PTS_S3_PREFIX, out_file)
  
  # command to copy cross section points parquet to S3
  copy_cs_pts_to_s3 <- paste0("aws s3 cp ", 
                              out_path,
                              " ", 
                              s3_save_uri,
                              ifelse(is.null(aws_profile), "", paste0(" --profile ", aws_profile))
                              )
  
  message(
    "Copy VPU ", path_df$vpu[i], " ML augmented cross sections to S3:\n - S3 copy command:\n'",
    copy_cs_pts_to_s3, "'",
    "'\n=========================="
    )

  system(copy_cs_pts_to_s3, intern = TRUE)
  
  end <- round(Sys.time())
  
  message("Finished augmenting cross section points with ML for VPU ", VPU)
  message("- Completed at: ", end)
  message("==========================")
  
  rm(net, net_subset, 
     final_cs, cs_pts, 
     inchannel_cs, bankful_cs)
  gc()
  gc()
  
}
  
    
  