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

# string to fill in "cs_source" column in output datasets
cs_source <- "hydrofabric3D"

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

# loop over the nextgen and transect datasets (by VPU) and extract point elevations across points on each transect line,
# then classify the points, and create a parquet file with hy_id, cs_id, pt_id, X, Y, Z data.
# Save parquet locally and upload to specified S3 bucket
for (i in 1:nrow(path_df)) {
  
  # i = 15
  
  start <- Sys.time()
  
  
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
  
  message("Creating VPU ", VPU, 
          " cross section points:\n - flowpaths: '", nextgen_file,
          "'\n - cross section points: '", cs_file, "'", 
          # "\n - waterbodies: '", ref_file, "'",
          "'\n - start time: '", start, "'"
  )
  
  # ----------------------------------------------------------------------------------------------------------------
  # ---- Read in data -----
  # ----------------------------------------------------------------------------------------------------------------
  
  # CONUS network parquet
  net <- 
    CONUS_NETWORK_URI %>% 
    arrow::open_dataset() %>%
    dplyr::filter(vpu == VPU) %>%
    dplyr::collect()
  
  # Cross section points parquet
  cs_pts    <- arrow::read_parquet(cs_pts_path)
  
  fline_net <- sf::read_sf(nextgen_path, layer = "flowpaths")
  
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
          "\n - Number of cross section points hy_ids: '", number_hyids_in_cs_pts, "'", 
          "'\n - Number of CONUS network hy_ids: '", number_hyids_in_conus_net, "'"
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
  
  cs_bottom_lengths <- get_cs_bottom_length(cs_pts)
  
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
  # dplyr::filter(owp_tw_inchan <= 0 | owp_tw_bf <= 0) 
  # dplyr::left_join( stream_order, by = "hy_id") 
  
  missing_cs <- 
    cs_pts %>% 
    dplyr::filter(is.na(hf_id) | is.na(owp_tw_inchan) | is.na(owp_y_inchan) | is.na(owp_tw_bf) | is.na(owp_y_bf) | is.na(owp_dingman_r)) %>% 
    hydrofabric3D::add_tmp_id()
  
  # Split the cross sections into 2 groups:
  # - "Inchannel cs" group are points with BOTH valid banks AND relief --> These get the INCHANNEL TW and Y values from the ML model
  # - "Bankful cs" group are points WITHOUT valid banks OR any relief  --> These get the BANKFUL TW and Y values from the ML model
  inchannel_cs <-
    cs_pts %>% 
    hydrofabric3D::add_tmp_id() %>% 
    dplyr::filter(!tmp_id %in% unique(missing_cs$tmp_id)) %>% 
    dplyr::select(-tmp_id) %>% 
    dplyr::filter(valid_banks & has_relief) %>% 
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
  
  split_kept_all_rows <- nrow(cs_pts) == nrow(bankful_cs) + nrow(inchannel_cs) + nrow(missing_cs)
  # split_kept_all_rows <- nrow(cs_pts) == nrow(bankful_cs) + nrow(inchannel_cs)
  
  if (!split_kept_all_rows) {
    warning("When splitting cross section points into 'bankful' and 'inchannel' groups, some points were not put in either group")
  }
  
    # Add bathymetry using "inchannel" estimates
    cs_bathy_inchannel <- add_cs_bathymetry(
      cross_section_pts = inchannel_cs,
      # cross_section_pts = dplyr::slice(inchannel_cs, 1:100000),
      top_width         = "owp_tw_inchan",
      depth             = "owp_y_inchan", 
      dingman_r         = "owp_dingman_r"
    )
    
    # Add bathymetry using "bankful" estimates
    cs_bathy_bankful <- add_cs_bathymetry(
      cross_section_pts = bankful_cs,
      top_width         = "owp_tw_bf",
      depth             = "owp_y_bf", 
      dingman_r         = "owp_dingman_r"
    )
    
    
  final_cs <- dplyr::bind_rows(cs_bathy_inchannel, cs_bathy_bankful)
  
  final_cs <- dplyr::bind_rows(
                  dplyr::select(cs_bathy_inchannel, 
                                -owp_tw_bf, -owp_y_bf, -hf_id),
                  dplyr::select(cs_bathy_bankful, 
                                -owp_tw_inchan, -owp_y_inchan, -hf_id)
                ) %>% 
    dplyr::group_by(hy_id, cs_id) %>% 
    tidyr::fill(
      c(cs_lengthm, Z_source, TW, DEPTH, DINGMAN_R)
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(
      -point_type, 
      -class,
      -bottom, -left_bank, -right_bank,
      -has_relief, -valid_banks
    )
  
  final_cs <- hydrofabric3D::classify_points(final_cs)
  # final_classified %>% 
  #   dplyr::filter(!valid_banks | !has_relief) 
  # system.time({
  #   final_classified <- 
  #     final_cs %>% 
  #     hydrofabric3D::classify_points()
  # })

  final_cs <- dplyr::bind_rows(
                dplyr::relocate(
                  final_cs,
                    hy_id, cs_id, pt_id, 
                    Z, relative_distance, cs_lengthm, class, point_type, 
                    X, Y, Z_source, bottom, left_bank, right_bank, valid_banks, has_relief, 
                    TW, DEPTH, DINGMAN_R, is_dem_point
                  ),
                dplyr::relocate(  
                  dplyr::mutate(
                    dplyr::select(missing_cs, 
                                  -tmp_id, -hf_id, -owp_tw_inchan, -owp_y_inchan, -owp_tw_bf, -owp_y_bf, -owp_dingman_r),
                    TW           = NA,
                    DEPTH        = NA, 
                    DINGMAN_R    = NA,
                    is_dem_point = TRUE
                  ),
                  hy_id, cs_id, pt_id, 
                  Z, relative_distance, cs_lengthm, class, point_type, 
                  X, Y, Z_source, bottom, left_bank, right_bank, valid_banks, has_relief, 
                  TW, DEPTH, DINGMAN_R, is_dem_point
                )
              )
  # ----------------------------------------------------------------------------------------------------------------
  # ---- Upload the cross section points parquet to S3 ----
  # ----------------------------------------------------------------------------------------------------------------
  
  # name of file and path to save transects gpkg too
  out_file <- paste0("nextgen_", path_df$vpu[i], "_cross_sections.parquet")
  out_path <- paste0(final_dir, out_file)
  
  message("Saving ML augmented cross section points to:\n - filepath: '", out_path, "'")
  
  # save cross section points as a parquet to out_path (lynker-spatial/02_cs_pts/cs_pts_<VPU num>.parquet)
  arrow::write_parquet(final_cs, out_path)
  
  # command to copy cross section points parquet to S3
  copy_cs_pts_to_s3 <- paste0("aws s3 cp ", out_path, " ", cs_pts_prefix, out_file,
                              ifelse(is.null(aws_profile), "", paste0(" --profile ", aws_profile)))
  
  message("Copy VPU ", path_df$vpu[i], " cross sections to S3:\n - S3 copy command:\n'", 
          paste0("aws s3 cp ", out_path, " ", CS_ML_PTS_S3_PREFIX, out_file,
                 ifelse(is.null(aws_profile), "", paste0(" --profile ", aws_profile))), 
          "'\n==========================")
  
  system(copy_cs_pts_to_s3, intern = TRUE)
  
  end <- Sys.time()
  
  message("Finished augmenting cross section points with ML for VPU ", VPU)
  message("- Completed at: ", end)
  message("==========================")
  
  rm(fixed_pts)
  gc()
  gc()
  
}
  
    
  