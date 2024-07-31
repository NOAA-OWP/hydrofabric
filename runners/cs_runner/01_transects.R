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
  
  # nextgen file and full path
  nextgen_file <- path_df$x[i]
  nextgen_path <- paste0(nextgen_dir, nextgen_file)
  
  vpu <- path_df$vpu[i]

  # Get FEMA by VPU directory and files for current VPU 
  fema_vpu_dir <- paste0(FEMA_VPU_SUBFOLDERS[grepl(paste0("VPU_", vpu), basename(FEMA_VPU_SUBFOLDERS))])
  # fema_vpu_dir <- paste0(FEMA_VPU_SUBFOLDERS[grepl(paste0("VPU_", vpu), basename(FEMA_VPU_SUBFOLDERS))], "/merged")

  vpu_fema_files <- list.files(fema_vpu_dir, full.names = TRUE)
  vpu_fema_file  <- vpu_fema_files[grepl(paste0(vpu, "_output.gpkg"), vpu_fema_files)]
  
  message("Creating VPU ", vpu, " transects:", 
          "\n - flowpaths: '",
          nextgen_file, "'",
           "\n - FEMA polygons: '", 
          basename(vpu_fema_file), "'"
          )
  
  # message("Creating VPU ", path_df$vpu[i], " transects:\n - flowpaths: '", nextgen_file, "'\n - model attributes: '", model_attr_file, "'")
  # sf::write_sf(
  #   dplyr::slice(dplyr::filter(flines, order == 2), 2), 
  #   "/Users/anguswatters/Desktop/example_flowline.gpkg"
  #   )
  
  # read in nextgen data
  flines <- sf::read_sf(nextgen_path, layer = "flowpaths")
  
  # calculate bankfull width
  flines <-
    flines %>%
    dplyr::mutate(
      bf_width = exp(0.700    + 0.365* log(tot_drainage_areasqkm))
    ) %>%
    dplyr::select(
      hy_id = id,
      lengthkm,
      tot_drainage_areasqkm,
      bf_width,
      mainstem,
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
  
  gc()
  
  time2 <- Sys.time()
  time_diff <- round(as.numeric(time2 - time1 ), 2)
  
  message("\n\n ---> Transects processed in ", time_diff)
  
  # name of file and path to save transects gpkg too
  out_file <- paste0("nextgen_", path_df$vpu[i], "_transects.gpkg")
  out_path <- paste0(transects_dir, out_file)
  
  # add cs_source column and rename cs_widths to cs_lengthm
  transects <- 
    transects %>%
    dplyr::mutate(
      cs_source = net_source
    )
  
  # ---------------------------------------------------------------------
  # --- Extend transects out to FEMA 100yr floodplains
  # ---------------------------------------------------------------------
  message("Reading in FEMA polygons...") 
  
  # fema polygons and transect lines
  fema <- sf::read_sf(vpu_fema_file)
  
  message("Simplifying FEMA polygons...")
  message(" - Number of geoms BEFORE simplifying: ", nrow(fema))
  
  # TODO: this should be a function argument OR removed, shouldn't probably forcibly and silently simplify the input polygons without user knowing..
  # keep 1% of the original points for speed
  fema <- rmapshaper::ms_simplify(fema, keep_shapes = T, keep = 0.01, sys = TRUE, sys_mem = 16)
  # fema <- rmapshaper::ms_simplify(fema, keep_shapes = T, keep = 0.1, sys = TRUE, sys_mem = 16)
  
  message(" - Number of geoms AFTER simplifying: ", nrow(fema))
  message("Extending transects out to FEMA 100yr floodplain polygon boundaries - (", Sys.time(), ")")
  
  transects <- 
    transects  %>%
    dplyr::left_join(
      dplyr::select(sf::st_drop_geometry(flines),
                    hy_id,
                    mainstem
      ),
      by = "hy_id" 
    )
  
  # TODO: make sure this 3000m extension distance is appropriate across VPUs 
  # TODO: also got to make sure that this will be feasible on memory on the larger VPUs...
  transects <- hydrofabric3D::extend_transects_to_polygons(
    transect_lines         = transects, 
    polygons               = fema, 
    flowlines              = flines, 
    crosswalk_id           = "hy_id",
    grouping_id     = "mainstem", 
    max_extension_distance = 3000 
  )
  
  message("FEMA extensions complete! - ( ", Sys.time(), " )")
  
  transects  <- dplyr::select(transects, -tmp_id)
  transects  <- hydrofabric3D::add_tmp_id(transects)
  
  transects <-
    transects %>%  
    # dplyr::select(-cs_lengthm) %>% 
    # dplyr::mutate(is_fema_extended = left_is_extended | right_is_extended) %>% 
    dplyr::select(
      hy_id, 
      cs_id, 
      cs_lengthm,
      # cs_lengthm = new_cs_lengthm, 
      cs_source,
      cs_measure,
      geometry
      # is_extended,
      # is_fema_extended,
      # geometry = geom
    )
  
  gc()
   
  # # ---------------------------------------------------------------------
  message("Saving transects to:\n - filepath: '", out_path, "'")
  
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
  
  rm(fema, transects, flines)
  gc()
}
