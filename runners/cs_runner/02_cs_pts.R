# Generate the flowlines layer for the final cross_sections_<VPU>.gpkg for each VPU
source("runners/workflow/config.R")

# load libraries
library(terrainSliceR)
library(dplyr)
library(sf)

# name of S3 bucket
s3_bucket <- "s3://lynker-spatial/"

# transect bucket prefix
cs_pts_prefix <- glue::glue("{s3_bucket}v20/3D/dem-cross-sections/")

# paths to nextgen datasets
nextgen_files <- list.files(nextgen_dir, full.names = FALSE)

# paths to nextgen datasets
transect_files <- list.files(transects_dir, full.names = FALSE)

# string to fill in "cs_source" column in output datasets
cs_source <- "terrainSliceR"

# ensure the files are in the same order and matched up by VPU
path_df <- align_files_by_vpu(
                x    = nextgen_files,
                y    = transect_files,
                base = base_dir
              )

# loop over the nextgen and transect datasets (by VPU) and extract point elevations across points on each transect line,
# then classify the points, and create a parquet file with hy_id, cs_id, pt_id, X, Y, Z data.
# Save parquet locally and upload to specified S3 bucket
for (i in 1:nrow(path_df)) {

  # nextgen file and full path
  nextgen_file <- path_df$x[i]
  nextgen_path <- glue::glue("{nextgen_dir}{nextgen_file}")
  
  # model attributes file and full path
  transect_file <- path_df$y[i]
  transect_path <- glue::glue("{transects_dir}{transect_file}")
  
  logger::log_info("\n\nCreating VPU {path_df$vpu[i]} cross section points:\n - flowpaths: '{nextgen_file}'\n - transects: '{transect_file}'")
  
  ################### 
  
    # read in transects data
    transects <- sf::read_sf(transect_path)
  
    # read in nextgen data
    flines <- sf::read_sf(nextgen_path, layer = "flowpaths")
  
    transects <-
      transects  %>%
      dplyr::rename(lengthm = cs_lengthm)

    # get start time for log messages
    time1 <- Sys.time()
    
    # get cross section point elevations
    cs_pts <- terrainSliceR::cross_section_pts(
              cs             = transects[lengths(sf::st_intersects(transects, flines)) == 1, ],
              points_per_cs  = NULL,
              min_pts_per_cs = 10,
              dem            = '/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/TIFF/USGS_Seamless_DEM_1.vrt'
              )
    
    # get end time for log messages
    time2 <- Sys.time()
    time_diff <- round(as.numeric(time2 - time1 ), 2)
    
    logger::log_info("\n\n ---> Cross section point elevations processed in {time_diff}")
    
    # Remove any cross section that has ANY missing (NA) Z values. 
    cs_pts <-
      cs_pts %>%
      # dplyr::filter(hy_id %in% c("wb-849054", "wb-845736")) %>% 
      dplyr::group_by(hy_id, cs_id) %>% 
      dplyr::filter(!any(is.na(Z))) %>% 
      dplyr::ungroup()
    
    # classify the cross section points
    cs_pts <-
      cs_pts %>%
      dplyr::rename(cs_widths = lengthm) %>%
      terrainSliceR::classify_points() %>%
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
  
  ################### 
  
  # name of file and path to save transects gpkg too
  out_file <- glue::glue("nextgen_{path_df$vpu[i]}_cross_sections.parquet")
  out_path <- glue::glue('{cs_pts_dir}{out_file}')
  
  logger::log_info("\n\nSaving cross section points to:\n - filepath: '{out_path}'")
  
  # save cross section points as a parquet to out_path (lynker-spatial/02_cs_pts/cs_pts_<VPU num>.parquet)
  arrow::write_parquet(cs_pts, out_path)
  
  # command to copy cross section points parquet to S3
  if (!is.null(aws_profile)) {
    copy_cs_pts_to_s3 <- glue::glue("aws s3 cp {out_path} {cs_pts_prefix}{out_file} --profile {aws_profile}")
  } else {
    copy_cs_pts_to_s3 <- glue::glue("aws s3 cp {out_path} {cs_pts_prefix}{out_file}")
  }
  
  logger::log_info("\n\nCopy VPU {path_df$vpu[i]} cross sections to S3:\n - S3 copy command:\n'{copy_cs_pts_to_s3}'\n==========================")
  system(copy_cs_pts_to_s3, intern = TRUE)

}
