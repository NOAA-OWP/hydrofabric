# Generate the transects + cs_pts + cross sections layers for a single flowlines domain file and DEM file 
source("runners/cs_runner/config.R")
source("runners/cs_runner/utils.R")

# # # # load libraries
# library(hydrofabric3D)
# library(dplyr)
# library(sf)
# install.packages("devtools")

# # transect bucket prefix
# S3_TRANSECTS_DIR <- paste0(S3_BUCKET_URI, VERSION, "/3D/transects/")

# paths to NEW DOMAIN datasets 
# NEW_DOMAIN_FILES    <- list.files(NEW_DOMAIN_FLOWLINES_DIR, full.names = TRUE)

# Unique Flowline ID column name
CROSSWALK_ID <- "id"

# ---------------------------------------------------------------------
# --- Read in flowlines
# ---------------------------------------------------------------------
# sf::st_layers(DOMAIN_WITH_FEMA_SUBSET_PATH)

# Subsetting area
aoi     <- sf::read_sf(DOMAIN_WITH_FEMA_SUBSET_PATH, layer = "divides")
aoi     <- rmapshaper::ms_simplify(aoi, keep = 0.05)

id_col <- "id"
# id_col <- "divide_id"

# query the conus.gpkg for matching IDs
ids        <- unique(aoi[[id_col]])
query_ids  <- paste0(paste0("'", ids, "'"), collapse= ", ")

gpkg_layers <- sf::st_layers(DOMAIN_WITH_FEMA_FLOWLINES_PATH)
layer       <- gpkg_layers$name[gpkg_layers$name == "flowpaths"]

wkt <- 
  aoi %>% 
  rmapshaper::ms_dissolve() %>% 
  # rmapshaper::ms_explode() %>% 
  # sf::st_as_sfc() %>%
  sf::st_sf() %>%
  sf::st_geometry() %>%
  sf::st_as_text()

# wkt   <- sf::st_as_text(sf::st_geometry(aoi))
# wkt <-
#   aoi %>%
#   # sf::st_bbox() %>%
#   sf::st_as_sfc() %>%
#   sf::st_sf() %>%
#   sf::st_geometry() %>%
#   sf::st_as_text()

# length(wkt)

# read in flowlines based on IDs in AOI
flines <- sf::read_sf(DOMAIN_WITH_FEMA_FLOWLINES_PATH, layer = "flowpaths",
                      wkt_filter = wkt
                      # query = sprintf("SELECT * FROM \"%s\" WHERE %s IN (%s)", layer, id_col, query_ids),
                      # query = query
                      )

# ---------------------------------------------------------------------
# --- Split flowlines by VPU
# ---------------------------------------------------------------------

# VPUs polygons
VPU_boundaries   <- sf::st_transform(nhdplusTools::vpu_boundaries, sf::st_crs(flines))

# add a VPU ID column to each flowline
flines <- add_intersects_ids(x = flines, y = VPU_boundaries, id_col = "VPUID")

# set of unique VPUs
VPU_IDS <- unnest_ids(flines$VPUID)

# all possible FEMA dirs
AOI_FEMA_DIRS <- FEMA_VPU_SUBFOLDERS[basename(FEMA_VPU_SUBFOLDERS) %in% paste0("VPU_", VPU_IDS) ]

# flines %>% 
#   sf::st_drop_geometry() %>% 
#   dplyr::group_by(VPUID) %>% 
#   dplyr::count()
# unique(flines$VPUID)

# calculate bankfull width
flines <-
  flines %>%
  dplyr::mutate(
    bf_width = hydrofabric3D::calc_powerlaw_bankful_width(tot_drainage_areasqkm)
    # bf_width = exp(0.700    + 0.365* log(tot_drainage_areasqkm))
  ) %>%
  dplyr::select(
    dplyr::any_of(CROSSWALK_ID),
    VPUID,
    # hy_id = id,
    lengthkm,
    tot_drainage_areasqkm,
    bf_width,
    mainstem
  ) %>% 
  hydroloom::rename_geometry("geometry")

# save the flowlines subset 
DOMAIN_WITH_FEMA_FLOWLINE_SUBSET_PATH <- paste0(DOMAIN_WITH_FEMA_FLOWLINES_DIR, "/flowlines_subset.gpkg")
sf::write_sf(
  flines, 
  DOMAIN_WITH_FEMA_FLOWLINE_SUBSET_PATH
)

# split the flowlines into groups by VPU
fline_groups <- dplyr::group_split(flines, VPUID)

# loop through each VPU group of flowlines and save a local copy
for (i in seq_along(fline_groups)) {
  flowlines  <- fline_groups[[i]]
  VPU        <- unique(flowlines$VPUID)
  
  # save the flowlines subset 
  vpu_flowlines_path <- paste0(DOMAIN_WITH_FEMA_FLOWLINES_DIR, "/", gsub(", ", "_", VPU), "_flowlines.gpkg")
  
  sf::write_sf(
    flowlines, 
    vpu_flowlines_path
  )
}

for (i in seq_along(fline_groups)) {

  flowlines  <- fline_groups[[i]]
  VPU        <-  unique(flowlines$VPUID)
  
  message("(", i, "/", length(fline_groups), ")\n", 
          " > Processing flowlines in VPU group: '", VPU, "'")
  
  # unique VPUs for group
  GROUP_VPU_IDS <- unnest_ids(flowlines$VPUID)
  
  # all FEMA dirs for the current area
  GROUP_FEMA_DIRS  <- FEMA_VPU_SUBFOLDERS[basename(FEMA_VPU_SUBFOLDERS) %in% paste0("VPU_", GROUP_VPU_IDS) ]
  GROUP_FEMA_FILES <- list.files(GROUP_FEMA_DIRS, full.names = T)[grepl("_output.gpkg", list.files(GROUP_FEMA_DIRS, full.names = T))]
  # GROUP_FEMA_FILES <- list.files(GROUP_FEMA_DIRS, full.names = T)
  # GROUP_FEMA_FILES <- GROUP_FEMA_FILES[grepl("_output.gpkg", GROUP_FEMA_FILES)]
  
  # create transect lines
  transects <- hydrofabric3D::cut_cross_sections(
    net               = flowlines,                        # flowlines network
    id                = CROSSWALK_ID,                     # Unique feature ID
    cs_widths         = pmax(50, flowlines$bf_width * 11),     # cross section width of each "id" linestring ("hy_id")
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
# }
  
  # add cs_source column and rename cs_widths to cs_lengthm
  transects <- 
    transects %>%
    dplyr::mutate(
      cs_source = CS_SOURCE
    )

  # ---------------------------------------------------------------------
  # --- Extend transects out to FEMA 100yr floodplains
  # ---------------------------------------------------------------------
  message("Reading in FEMA polygons...") 
  
  # unique(flowlines$VPUID)
  
  # fema polygons and transect lines
  # fema <- sf::read_sf(vpu_fema_file)
  
  # read each FEMA geopackage into a list 
  fema <- lapply(GROUP_FEMA_FILES, function(gpkg) sf::read_sf(gpkg))
  
  fema <- 
    fema %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(
      fema_id = 1:dplyr::n()
    ) 

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
      dplyr::select(sf::st_drop_geometry(flowlines),
                    dplyr::any_of(CROSSWALK_ID),
                    mainstem
      ),
      by = CROSSWALK_ID
    )
  
  # TODO: make sure this 3000m extension distance is appropriate across VPUs 
  # TODO: also got to make sure that this will be feasible on memory on the larger VPUs...
  transects <- hydrofabric3D::extend_transects_to_polygons(
    transect_lines         = transects, 
    polygons               = fema, 
    flowlines              = flowlines, 
    crosswalk_id           = CROSSWALK_ID,
    grouping_id            = "mainstem", 
    max_extension_distance = 3000 
  )
  
  # mapview::mapview(transects, color = "green") + 
   # mapview::mapview(transects2, color = "red")
  
  transects <- 
    transects %>% 
    hydrofabric3D::add_tmp_id(x = CROSSWALK_ID) %>% 
    dplyr::mutate(is_extended = FALSE) %>%
    dplyr::select(
      dplyr::any_of(CROSSWALK_ID),
      cs_id, 
      cs_lengthm,
      cs_source,
      cs_measure,
      is_extended,
      geometry
    )
  
  out_path <- paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/",   
                     paste(GROUP_VPU_IDS, collapse = "_"), "_transects.gpkg" 
                     )
  
  message("Writting transect lines for VPU group: '", VPU, "'",
          "\n > '", out_path, "'")
  
  sf::write_sf(transects, out_path)
  
  message("Finished writting transects!")
  
}
# ---------------------------------------------------------------------
# --- Bind together transects into single dataset
# ---------------------------------------------------------------------

transect_paths <- list.files(DOMAIN_WITH_FEMA_TRANSECTS_DIR, full.names = T)
flowline_paths <- list.files(DOMAIN_WITH_FEMA_FLOWLINES_DIR, full.names = T)
flowline_paths <- flowline_paths[!basename(flowline_paths) %in% c("flowlines_subset.gpkg", "ls_conus.gpkg")]

# read each transect geopackages into a list 
transects <- lapply(transect_paths, function(gpkg) sf::read_sf(gpkg))
transects[[5]]

transects <- 
  transects %>% 
  dplyr::bind_rows()

paths_df <- data.frame(
  t = transect_paths,
  f = flowline_paths
) %>% 
  dplyr::mutate(
    vpu = gsub("_transects.gpkg", "", basename(t))
  )

for (i in 1:nrow(paths_df)) {
  
  VPU     <- paths_df$vpu[i]
  t_path  <- paths_df$t[i]
  f_path  <- paths_df$f[i]
  
  message("Creating VPU ", VPU, " cross section points:", 
          "\n - flowpaths: '", f_path, "'",
          "\n - transects: '", t_path, "'"
  )
  
  ################### 
  message("Reading in transects...\n > ", t_path)
  # read in transects data
  transects <- sf::read_sf(t_path)
  
  message("Reading in flowlines... \n > ", f_path)
  # read in nextgen data
  flines <- sf::read_sf(f_path)
  
  message("Extracting cross section points...")
  # ----------------------------------------------------------------------------------------------------------------
  # ---- STEP 1: Extract cs points from DEM ----
  # ----------------------------------------------------------------------------------------------------------------
  # system.time({
  # get cross section point elevations
  cs_pts <- hydrofabric3D::cross_section_pts(
    cs             = transects,
    crosswalk_id   = CROSSWALK_ID,
    points_per_cs  = NULL,
    min_pts_per_cs = 10,
    dem            = DEM_URL
  )
  
  # ----------------------------------------------------------------------------------------------------------------
  # ---- STEP 2: Remove any cross section that has ANY missing (NA) Z values, and classify the points  ----
  # ----------------------------------------------------------------------------------------------------------------
  # cs_pts2 %>% 
  #   dplyr::slice(1:200) %>% 
  #   dplyr::rename(hy_id = id) %>% 
  #   hydrofabric3D::plot_cs_pts(x = "pt_id", color = "point_type")
  
  cs_pts <- 
    # cs_pts2 <- 
    cs_pts %>% 
    hydrofabric3D::drop_incomplete_cs_pts(CROSSWALK_ID) %>% 
    hydrofabric3D::classify_points(
      crosswalk_id             = CROSSWALK_ID, 
      pct_of_length_for_relief = PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF
    )  
  
  # })
  
  ids_original_cs_pts <- hydrofabric3D::add_tmp_id(cs_pts, x = CROSSWALK_ID)$tmp_id  
  # ids_original_cs_pts <- hydrofabric3D::add_tmp_id(cs_pts2)$tmp_id
  
  # sf::write_sf(cs_pts2, "/Users/anguswatters/Desktop/test_improve_cs_pts_classified_11.gpkg")
  # sf::write_sf(cs_pts, "/Users/anguswatters/Desktop/test_improve_cs_pts_classified_11_2.gpkg")
  
  
  # ----------------------------------------------------------------------------------------------------------------
  # ---- STEP 3: Try to rectify any no relief and invalid banks cross sections ----
  # ----------------------------------------------------------------------------------------------------------------
  
  # system.time({
    fixed_pts <- hydrofabric3D::get_improved_cs_pts(
      cs_pts         = cs_pts,    # cross section points generated from hydrofabric3D::cross_section_pts()
      net            = flines,    # original flowline network
      # net            = flines,    # original flowline network
      transects      = transects, # original transect lines
      crosswalk_id   = CROSSWALK_ID,
      points_per_cs  = NULL, 
      min_pts_per_cs = 10, # number of points per cross sections
      dem            = DEM_URL, # DEM to extract points from
      scale          = EXTENSION_PCT, # How far to extend transects if the points need to be rechecked
      pct_of_length_for_relief = PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF, # percent of cross sections length to be needed in relief calculation to consider cross section to "have relief"
      fix_ids = FALSE,
      verbose = TRUE
    )
  # })
  
  ids_after_fixed_pts <- hydrofabric3D::add_tmp_id(cs_pts, x = CROSSWALK_ID)$tmp_id  
  
  # ----------------------------------------------------------------------------------------------------------------
  # ---- Update transects with extended transects (if exists) ----
  # ----------------------------------------------------------------------------------------------------------------
  
  out_transects <- match_transects_to_extended_cs_pts(
    transect_lines = transects, 
    fixed_cs_pts   = fixed_pts, 
    crosswalk_id   = CROSSWALK_ID
  )
  
  # ----------------------------------------------------------------------------------------------------------------
  # ---- Re enumerate the transects & cross section points "cs_id" ----
  # ----------------------------------------------------------------------------------------------------------------
  
  # fixed_pts      <- hydrofabric3D:::renumber_cs_ids(df = fixed_pts, crosswalk_id = "hy_id")
  # out_transects  <- hydrofabric3D:::renumber_cs_ids(
  #                                                   df            = dplyr::mutate(out_transects, pt_id = 1), 
  #                                                   crosswalk_id  = "hy_id"
  #                                                   ) %>% 
  #                                                   dplyr::select(-pt_id)
  
  fixed_pts      <- hydrofabric3D:::renumber_cs_ids(df = fixed_pts, crosswalk_id = CROSSWALK_ID)
  out_transects  <- hydrofabric3D:::renumber_cs_ids(df = out_transects, crosswalk_id = CROSSWALK_ID)
  
  # ----------------------------------------------------------------------------------------------------------------
  # ---- STEP 4: Update transects with extended transects (if exists) ----
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
      dplyr::any_of(CROSSWALK_ID),
      cs_id, 
      pt_id,
      cs_lengthm,
      relative_distance,
      X, Y, Z,
      class, point_type,
      bottom, left_bank, right_bank, valid_banks, has_relief # newly added columns (03/06/2024)
    )
  
  # add Z_source column for source of elevation data
  fixed_pts <-
    fixed_pts %>%
    dplyr::mutate(
      Z_source = CS_SOURCE
    ) %>%
    dplyr::relocate(
      dplyr::any_of(CROSSWALK_ID),
      cs_id, pt_id, cs_lengthm, relative_distance, X, Y, Z, Z_source, 
                    class, point_type, 
                    bottom, left_bank, right_bank, valid_banks, has_relief)
  
  ids_before_align <- hydrofabric3D::add_tmp_id(fixed_pts, x = CROSSWALK_ID)$tmp_id
  
  message("Aligning banks and smoothing bottoms...")
  fixed_pts <- hydrofabric3D::align_banks_and_bottoms(cs_pts = fixed_pts, crosswalk_id = CROSSWALK_ID)
  
  ids_after_align <- hydrofabric3D::add_tmp_id(fixed_pts, x = CROSSWALK_ID)$tmp_id
  
  message("Reclassifying cross section points...")
  
  fixed_pts <- hydrofabric3D::classify_points(
    cs_pts                    = fixed_pts, 
    crosswalk_id              = CROSSWALK_ID,
    pct_of_length_for_relief  = PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF
  )
  
  ids_after_reclassify <- hydrofabric3D::add_tmp_id(fixed_pts, x = CROSSWALK_ID)$tmp_id
  
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
  
  sf::write_sf(
    out_transects,
    paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/", VPU, "_transects.gpkg") 
  )
  
  arrow::write_parquet(
    fixed_pts,
    paste0(DOMAIN_WITH_FEMA_CS_PTS_DIR, "/", VPU, "_cs_pts.gpkg") 
  )
  
}

# TODO: Save all cross section points + transects into single geopackages (transects.gpkg, cs_pts.gpkg)


