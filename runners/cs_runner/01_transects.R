# Generate the flowlines layer for the final cross_sections_<VPU>.gpkg for each VPU
source("runners/cs_runner/config.R")

# # # # load libraries
# library(hydrofabric3D)
# library(dplyr)
# library(sf)
# install.packages("devtools")

# transect bucket prefix
transects_prefix <- paste0(s3_bucket, version_prefix, "/3D/transects/")

# paths to nextgen datasets and model attribute parquet files
nextgen_files    <- list.files(nextgen_dir, full.names = FALSE)
model_attr_files <- list.files(model_attr_dir, full.names = FALSE)

# string to fill in "cs_source" column in output datasets
net_source <- "hydrofabric3D"

# ensure the files are in the same order and matched up by VPU
path_df <- align_files_by_vpu(
                x    = nextgen_files,
                y    = model_attr_files,
                base = base_dir
                )

# loop over each VPU and generate cross sections, then save locally and upload to S3 bucket
for(i in 1:nrow(path_df)) {
  # i = 8
  # nextgen file and full path
  nextgen_file <- path_df$x[i]
  nextgen_path <- paste0(nextgen_dir, nextgen_file)
  
  vpu <- path_df$vpu[i]

  # Get FEMA by VPU directory and files for current VPU 
  fema_vpu_dir <- paste0(FEMA_VPU_SUBFOLDERS[grepl(paste0("VPU_", vpu), basename(FEMA_VPU_SUBFOLDERS))], "/merged")
  # fema_vpu_dir <- paste0(FEMA_VPU_SUBFOLDERS[grepl(paste0("VPU_", vpu), basename(FEMA_VPU_SUBFOLDERS))], "/merged")

  vpu_fema_files <- list.files(fema_vpu_dir, full.names = TRUE)
  vpu_fema_file <- vpu_fema_files[grepl(paste0(vpu, "_union.gpkg"), vpu_fema_files)]
  
  # fema polygons and transect lines
  fema <- sf::read_sf(vpu_fema_file)
  
  # # model attributes file and full path
  # model_attr_file <- path_df$y[i]
  # model_attr_path <- paste0(model_attr_dir, model_attr_file)

  message("Creating VPU ", vpu, " transects:", 
          "\n - flowpaths: '",
          nextgen_file, "'",
           "\n - FEMA polygons: ", basename(vpu_fema_file)
          )
  # message("Creating VPU ", path_df$vpu[i], " transects:\n - flowpaths: '", nextgen_file, "'\n - model attributes: '", model_attr_file, "'")
  
  # read in nextgen data
  flines <- sf::read_sf(nextgen_path, layer = "flowpaths")
  
  # # model attributes
  # model_attrs <- arrow::read_parquet(model_attr_path)

  # # join flowlines with model atttributes
  # flines <- dplyr::left_join(
  #   flines,
  #   dplyr::select(
  #     model_attrs,
  #     id, eTW
  #   ), by = "id")

  # calculate bankfull width
  flines <-
    flines %>%
    dplyr::mutate(
      bf_width = exp(0.700    + 0.365* log(tot_drainage_areasqkm))
    ) %>%
    # dplyr::mutate( bf_width = 11 * eTW) %>%
    # dplyr::mutate( # if there are any NAs, use exp(0.700    + 0.365* log(tot_drainage_areasqkm)) equation to calculate bf_width
    #   bf_width = dplyr::case_when(
    #     is.na(bf_width) ~ exp(0.700    + 0.365* log(tot_drainage_areasqkm)),
    #     TRUE            ~ bf_width
    #   )) %>%
    dplyr::select(
      hy_id = id,
      lengthkm,
      tot_drainage_areasqkm,
      bf_width,
      geometry = geom
    )

  # flines$bf_width <- ifelse(is.na(flines$bf_width),  exp(0.700    + 0.365* log(flines$tot_drainage_areasqkm)), flines$bf_width)

  time1 <- Sys.time()

  # create transect lines
  transects <- hydrofabric3D::cut_cross_sections(
    net               = flines,                        # flowlines network
    id                = "hy_id",                       # Unique feature ID
    cs_widths         = pmax(50, flines$bf_width * 11),     # cross section width of each "id" linestring ("hy_id")
    # cs_widths         = pmax(50, flines$bf_width),     # cross section width of each "id" linestring ("hy_id")
    num               = 10,                            # number of cross sections per "id" linestring ("hy_id")
    smooth            = TRUE,                          # smooth lines
    densify           = 3,                             # densify linestring points
    rm_self_intersect = TRUE,                          # remove self intersecting transects
    fix_braids        = FALSE,                         # whether to fix braided flowlines or not
    #### Arguments used for when fix_braids = TRUE     # TODO: these methods need revision in hydrofabric3D to allow for more flexible processing for data that is NOT COMID based (i.e. hy_id)    
    # terminal_id       = NULL,
    # braid_threshold   = NULL,
    # version           = 2,
    # braid_method      = "comid",
    # precision         = 1,
    add               = TRUE                           # whether to add back the original data
  )
    
  time2 <- Sys.time()
  time_diff <- round(as.numeric(time2 - time1 ), 2)
  
  message("\n\n ---> Transects processed in ", time_diff)
  
  # name of file and path to save transects gpkg too
  out_file <- paste0("nextgen_", path_df$vpu[i], "_transects.gpkg")
  out_path <- paste0(transects_dir, out_file)
  
  message("Saving transects to:\n - filepath: '", out_path, "'")
  
  # add cs_source column and rename cs_widths to cs_lengthm
  transects <- 
    transects %>%
    dplyr::mutate(
      cs_source = net_source
    )
  
  # TODO: make sure this 3000m extension distance is appropriate across VPUs 
  # TODO: also got to make sure that this will be feasible on memory on the larger VPUs...
  transects <- hydrofabric3D::get_transect_extension_distances_to_polygons(
                                                  transect_lines         = transects, 
                                                  polygons               = fema, 
                                                  flines                 = flines, 
                                                  max_extension_distance = 3000 
                                                  )
                                                
  # save transects with only columns to be uploaded to S3 (lynker-spatial/01_transects/transects_<VPU num>.gpkg)
  sf::write_sf(
    # save dataset with only subset of columns to upload to S3
    dplyr::select(transects, 
                  hy_id,
                  cs_source,
                  cs_id,
                  cs_measure,
                  cs_lengthm,
                  # sinuosity,
                  geometry
    ),
    out_path
    )
  
  transects <- sf::read_sf(out_path)

  
  # command to copy transects geopackage to S3
  copy_to_s3 <- paste0("aws s3 cp ", out_path, " ", transects_prefix, out_file, 
                   ifelse(is.null(aws_profile), "", paste0(" --profile ", aws_profile))
                   )
  
  message("Copy VPU ", path_df$vpu[i], " transects to S3:\n - S3 copy command:\n'", 
          copy_to_s3, 
          "'\n==========================")
  
  system(copy_to_s3, intern = TRUE)
  
  message("Overwritting local copy of transects to include 'is_extended' column...\n==========================")
  # Overwrite transects with additional columns for development purposes (is_extended) to have a local copy of dataset with information about extensions
  sf::write_sf(
    dplyr::select(
      dplyr::mutate(transects, is_extended = FALSE),
      hy_id,
      cs_source,
      cs_id,
      cs_measure,
      cs_lengthm,
      # sinuosity,
      is_extended,
      geometry
    ),
    out_path
  )
}
