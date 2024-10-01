# Generate the transects + cs_pts + cross sections layers for a single flowlines domain file and DEM file 
source("runners/cs_runner/config.R")
source("runners/cs_runner/utils.R")

# load libraries
# library(hydrofabric3D)
# library(dplyr)
# library(sf)

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

# bad_ids <- c("wb-14538", "wb-14686",  "wb-14687")
# flines <- 
#   flines %>% 
#   dplyr::filter(id %in% bad_ids)
# ---------------------------------------------------------------------
# --- Split flowlines by VPU
# ---------------------------------------------------------------------

# VPUs polygons
VPU_boundaries   <- sf::st_transform(nhdplusTools::vpu_boundaries, sf::st_crs(flines))

# add a VPU ID column to each flowline
flines <- add_intersects_ids(x = flines, y = VPU_boundaries, id_col = "VPUID")

# TODO: improve this, manual remap some VPUIDs so that there are 
# TODO: less small subsets of flowlines because a small bit of a different VPU is intersected 
flines <- 
  flines %>% 
  dplyr::mutate(
    VPUID = dplyr::case_when(
      VPUID == "12"     ~ "11, 12, 13",
      VPUID == "12, 13" ~ "11, 12, 13",
      VPUID == "11, 12" ~ "11, 12, 13",
      VPUID == "01"     ~ "01, 02",
      TRUE              ~ VPUID
    )
  ) 

# set of unique VPUs
VPU_IDS <- unnest_ids(flines$VPUID)
VPU_IDS

# all possible FEMA dirs
AOI_FEMA_DIRS <- FEMA_VPU_SUBFOLDERS[basename(FEMA_VPU_SUBFOLDERS) %in% paste0("VPU_", VPU_IDS) ]

flines %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(VPUID) %>%
  dplyr::count()
# unique(flines$VPUID)

# calculate bankfull width
flines <-
  flines %>%
  dplyr::mutate(
    bf_width = hydrofabric3D::calc_powerlaw_bankful_width(tot_drainage_areasqkm)
    # bf_width = exp(0.700    + 0.365* log(tot_drainage_areasqkm))
  ) %>%
  # hydrofabric3D::add_powerlaw_bankful_width("tot_drainage_areasqkm", 50) %>% 
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
  vpu_flowlines_path <- paste0(DOMAIN_WITH_FEMA_VPU_SUBSETS_DIR, "/", gsub(", ", "_", VPU), "_flowlines.gpkg")
  
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
    # cs_widths         = pmax(50, flowlines$bf_width),     # cross section width of each "id" linestring ("hy_id")
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
  
  out_path <- paste0(DOMAIN_WITH_FEMA_VPU_SUBSETS_DIR, "/",   
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
vpu_subset_paths <- list.files(DOMAIN_WITH_FEMA_VPU_SUBSETS_DIR, full.names = T)

transect_paths <- vpu_subset_paths[grepl("_transects.gpkg", vpu_subset_paths)]
flowline_paths <- vpu_subset_paths[grepl("_flowlines.gpkg", vpu_subset_paths)]

# read each transect geopackages into a list 
transects <- lapply(transect_paths, function(gpkg) sf::read_sf(gpkg))
# transects[[5]]

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
  # i = 2
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
  flowlines <- sf::read_sf(f_path)
  
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
      net            = flowlines,    # original flowline network
      # net            = flowlines,    # original flowline network
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
    paste0(DOMAIN_WITH_FEMA_VPU_SUBSETS_DIR, "/", VPU, "_transects.gpkg") 
  )
  
  arrow::write_parquet(
    fixed_pts,
    paste0(DOMAIN_WITH_FEMA_VPU_SUBSETS_DIR, "/", VPU, "_cs_pts.parquet") 
  )
  
}

# TODO: Save all cross section points + transects into single geopackages (transects.gpkg, cs_pts.gpkg)

# ---------------------------------------------------------------------------------
# ---- Merge all the transects / cs points into single files  
# ---------------------------------------------------------------------------------

vpu_subset_paths <- list.files(DOMAIN_WITH_FEMA_VPU_SUBSETS_DIR, full.names = T)

# ---------------------------------------------------------------------------------
# ---- Read all transects and save to single gpkg ----
# ---------------------------------------------------------------------------------

transect_paths <- vpu_subset_paths[grepl("_transects.gpkg", vpu_subset_paths)]

transects <- lapply(transect_paths, function(gpkg) {sf::read_sf(gpkg) }) %>% 
  dplyr::bind_rows() %>% 
  hydroloom::rename_geometry("geometry")

TRANSECTS_OUTPUT_PATH <- paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/transects.gpkg") 

sf::write_sf(
  transects,
  TRANSECTS_OUTPUT_PATH
)

# ---------------------------------------------------------------------------------
# ---- Read all CS PTs and save to single parquet ----
# ---------------------------------------------------------------------------------

cs_pts_paths   <- vpu_subset_paths[grepl("_cs_pts.parquet", vpu_subset_paths)]
cs_pts <- lapply(cs_pts_paths, function(prq) arrow::read_parquet(prq)) %>% 
  dplyr::bind_rows()

CS_PTS_OUTPUT_PATH <- paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/cs_pts.parquet")

arrow::write_parquet(
  cs_pts,
  CS_PTS_OUTPUT_PATH
)

# ---------------------------------------------------------------------------------
# ---- Inject ML predicted top widths / Dingman's R ----
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# ---- Read in ML data ----
# ---------------------------------------------------------------------------------
VPU_ML_BATHYMETRY_PATHS <- list.files(DOMAIN_WITH_FEMA_ML_DIR, full.names = T)

ML_CROSSWALK_ID <- "id"

# ml_outputs <- lapply(VPU_ML_BATHYMETRY_PATHS, function(prq) {
#   vpu_id     <- gsub(".*ml/([a-zA-Z0-9]+).*", "\\1", prq)
#   arrow::read_parquet(prq) %>% 
#     dplyr::mutate(vpu_id = vpu_id) 
# }
# ) %>% 
#   dplyr::bind_rows() %>% 
#   dplyr::select(
#     dplyr::any_of(ML_CROSSWALK_ID),
#     vpu_id,
#     owp_y_bf, owp_y_inchan, 
#     owp_tw_bf, owp_tw_inchan,
#     owp_dingman_r
#                 )
# 
# # rename ML_CROSSWALK_ID (unique ID) to match the CROSSWALK_ID in CS PTS
# # TODO: This assumes the IDs do correspond with eachother... (built from same flowlines network)
# names(ml_outputs)[names(ml_outputs) == ML_CROSSWALK_ID] = CROSSWALK_ID
# 
# # Keep only distinct ID rows
# ml_outputs <- 
  # ml_outputs %>%
  # dplyr::distinct(
  #   dplyr::across(dplyr::any_of(CROSSWALK_ID)),
  #   vpu_id,
  #   owp_y_bf, owp_y_inchan,
  #   owp_tw_bf, owp_tw_inchan,
  #   owp_dingman_r
  # )

# sf::st_layers(DOMAIN_WITH_FEMA_FLOWLINES_PATH)
# rm(ml_outputs, ml)
ml_outputs <- sf::read_sf(DOMAIN_WITH_FEMA_FLOWLINES_PATH, layer = "flowpath-attributes-ml")

ml_outputs <- 
  ml_outputs %>% 
  dplyr::select(
    dplyr::any_of(ML_CROSSWALK_ID),
    vpuid,
    owp_y_bf = YCC, 
    owp_y_inchan = Y,
    owp_tw_bf = TopWdthCC, 
    owp_tw_inchan = TopWdth,
    owp_dingman_r = dingman_r
                )
# 
# # rename ML_CROSSWALK_ID (unique ID) to match the CROSSWALK_ID in CS PTS
# # TODO: This assumes the IDs do correspond with eachother... (built from same flowlines network)
names(ml_outputs)[names(ml_outputs) == ML_CROSSWALK_ID] = CROSSWALK_ID
# 
# # Keep only distinct ID rows
ml_outputs <-
  ml_outputs %>%
  dplyr::distinct(
    dplyr::across(dplyr::any_of(CROSSWALK_ID)),
    vpu_id,
    owp_y_bf, owp_y_inchan,
    owp_tw_bf, owp_tw_inchan,
    owp_dingman_r
  )

# ---------------------------------------------------------------------------------
# ---- Read in CS PTS data ----
# ---------------------------------------------------------------------------------
CS_PTS_OUTPUT_PATH <- paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/cs_pts.parquet")

cs_pts  <- arrow::read_parquet(CS_PTS_OUTPUT_PATH)

# ---------------------------------------------------------------------------------
# ---- Join CS PTS  with ML data ---
# ---------------------------------------------------------------------------------
message(round(Sys.time()), " - Joining ML width/depths estimates to cross section points...")

# ml_outputs %>% 
#   dplyr::group_by(id) %>% 
#   dplyr::count(id) %>% 
#   dplyr::arrange(-n)

# join the ML outputs data to the cross section points
cs_pts <-
  cs_pts %>% 
  dplyr::left_join(
    dplyr::select(ml_outputs, 
                  dplyr::any_of(CROSSWALK_ID),
                  owp_tw_inchan, 
                  owp_y_inchan, 
                  owp_tw_bf, 
                  owp_y_bf, 
                  owp_dingman_r
                  ),
    by = CROSSWALK_ID
  ) 

# ---------------------------------------------------------------------------------
# ---- Fixing negative depths/widths estimates ----
# ---------------------------------------------------------------------------------
message(round(Sys.time()), " - Replacing any negative width/depth estimates with cross section bottom lengths...")

cs_bottom_lengths <- hydrofabric3D::get_cs_bottom_length(cross_section_pts = cs_pts, crosswalk_id = CROSSWALK_ID)

# TODO: for now we replace any negative TW values with the length of the bottom of the cross section
# TODO: This method + the negative model output values both need to be looked into (04/05/2024)
cs_pts <-
  cs_pts %>% 
  dplyr::left_join(
    cs_bottom_lengths, 
    by = c(CROSSWALK_ID, "cs_id")
    # by = c("hy_id", "cs_id")
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

# extract any cross sections that didn't get matched with a "hf_id" and (or?) no ML data
# TODO: look at this stuff with Arash (04/09/2024)
missing_cs <- 
  cs_pts %>% 
  dplyr::filter(
          is.na(.data[[CROSSWALK_ID]]) | 
          is.na(owp_tw_inchan) | is.na(owp_y_inchan) | 
          is.na(owp_tw_bf) | is.na(owp_y_bf) |
          is.na(owp_dingman_r)
        ) %>% 
  hydrofabric3D::add_tmp_id(x = CROSSWALK_ID)

# TODO: Delete this, but time being keeping this to inspect mismatch in between "hy_id" and "hf_id"
# readr::write_csv(
#   dplyr::select(missing_cs, -tmp_id), 
#   paste0(META_PATH, "nextgen_", path_df$vpu[i], "_cross_sections_missing_hf_ids.csv")
# )

# Split the cross sections into 2 groups:
# - "Inchannel cs" group are points with BOTH valid banks AND relief --> These get the INCHANNEL TW and Y values from the ML model
# - "Bankful cs" group are points WITHOUT valid banks OR any relief  --> These get the BANKFUL TW and Y values from the ML model
inchannel_cs <-
  cs_pts %>% 
  hydrofabric3D::add_tmp_id(x = CROSSWALK_ID) %>% 
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
  hydrofabric3D::add_tmp_id(x = CROSSWALK_ID) %>% 
  dplyr::filter(!tmp_id %in% unique(missing_cs$tmp_id)) %>% 
  dplyr::select(-tmp_id) %>% 
  dplyr::filter(!valid_banks | !has_relief) %>% 
  dplyr::rename(
    TW        = owp_tw_bf, 
    DEPTH     = owp_y_bf,
    DINGMAN_R = owp_dingman_r
  )


# sanity check that all rows are accounted for after splitting up data
split_kept_all_rows <- nrow(cs_pts) == (nrow(bankful_cs) + nrow(inchannel_cs) + nrow(missing_cs))
# split_kept_all_rows <- nrow(cs_pts) == nrow(bankful_cs) + nrow(inchannel_cs)

if (!split_kept_all_rows) {
  warning(paste0("When splitting cross section points into 'bankful' and 'inchannel' groups,", 
                 "\nsome points were not put in either group.", 
                 "\nLikely due to 'valid_banks' and/or 'has_relief' columns have either missing ", 
                 "values or contain values other than TRUE/FALSE")
  )
}

message(round(Sys.time()), " - Adding cross section bathymetry using inchannel widths/depths estimates...")
# tmp <- 
#   inchannel_cs %>% 
#   dplyr::slice(1:10000)
system.time({
  
# Add bathymetry using "inchannel" estimates
inchannel_cs <- hydrofabric3D::add_cs_bathymetry(
  cross_section_pts = inchannel_cs,
  # cross_section_pts = tmp,
  crosswalk_id      = CROSSWALK_ID
)

})



# arrow::write_parquet(inchannel_cs, "/Users/anguswatters/Desktop/test_ml_cs_pts_06.parquet")
# ml_subset %>%
#   dplyr::filter(hy_id == "wb-1005207") %>%
#   dplyr::select(owp_y_inchan, owp_tw_inchan) %>% 
#   .$owp_y_inchan
message(round(Sys.time()), " - Adding cross section bathymetry using bankful widths/depths estimates...")
system.time({
  
# Add bathymetry using "bankful" estimates
bankful_cs <- hydrofabric3D::add_cs_bathymetry(
  cross_section_pts = bankful_cs,
  # cross_section_pts = dplyr::slice(bankful_cs, 1:10000),
  crosswalk_id      = CROSSWALK_ID
)

})

# combine the inchannel and bankful cross section points back together, fill out missing values and reclassify the points
final_cs <- dplyr::bind_rows(
                dplyr::select(
                  inchannel_cs,
                  # inchannel_cs2,
                  # -hf_id, 
                  -TW, -DEPTH, -DINGMAN_R, 
                  # -is_dem_point,
                  -dplyr::starts_with("owp")
                ),
                dplyr::select(
                  bankful_cs,
                  # bankful_cs2,
                  # -hf_id, 
                  -TW, -DEPTH, -DINGMAN_R, 
                  # -is_dem_point,
                  -dplyr::starts_with("owp")
                ),
                dplyr::select(
                  dplyr::mutate(
                    missing_cs,
                    is_dem_point = FALSE
                  ),
                  # -hf_id, 
                  # -is_dem_point,
                  -dplyr::starts_with("owp"), 
                  -tmp_id
                )
              ) %>%
  dplyr::group_by(dplyr::across(dplyr::any_of(c(CROSSWALK_ID, "cs_id")))) %>% 
  # dplyr::group_by(hy_id, cs_id) %>%
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

# arrow::write_parquet(final_cs, "/Users/anguswatters/Desktop/tmp.parquet")
# final_cs <- arrow::read_parquet("/Users/anguswatters/Desktop/tmp.parquet")


message(round(Sys.time()), " - Reclassifying cross section point types...")

# reclassify
final_cs <- hydrofabric3D::classify_points(cs_pts = final_cs, 
                                            crosswalk_id = CROSSWALK_ID,
                                            pct_of_length_for_relief = PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF
                                            )

starting_uids <- hydrofabric3D::get_unique_tmp_ids(cs_pts, x = CROSSWALK_ID)
ending_uids   <- hydrofabric3D::get_unique_tmp_ids(final_cs, x = CROSSWALK_ID)

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

# tmp_path <-   paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/final_cs2.parquet")
# message("saving file", tmp_path)
# 
# final_cs <- arrow::read_parquet(tmp_path)
# 
# # save cross section points as a parquet to out_path (domain/outputs/cross_sections.parquet)
# arrow::write_parquet(
#   dplyr::select(final_cs, 
#                 -is_dem_point
#   ), 
#   paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/final_cs2.parquet")
# )

# ---------------------------------------------------------------------------------
# ---- Write final cross section points data ----
# ---------------------------------------------------------------------------------

CROSS_SECTIONS_OUTPUT_PATH <- paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/cross_sections.parquet")
message(round(Sys.time()), " - Saving ML augmented cross section points to:\n - filepath: '", CROSS_SECTIONS_OUTPUT_PATH, "'")

# save cross section points as a parquet to out_path (domain/outputs/cross_sections.parquet)
arrow::write_parquet(
  dplyr::select(final_cs, 
                -is_dem_point
  ), 
  CROSS_SECTIONS_OUTPUT_PATH
)

INTERNAL_CROSS_SECTIONS_PATH <- paste0("/Users/anguswatters/Desktop/cross_sections_is_dem_point.parquet")
arrow::write_parquet(
  final_cs,
  INTERNAL_CROSS_SECTIONS_PATH
)

# ---------------------------------------------------------------------------------
# ---- Substitue diffusive domain DEMs Z values ---- 
# ---------------------------------------------------------------------------------
final_cs <- arrow::read_parquet(CROSS_SECTIONS_OUTPUT_PATH)
# INTERNAL_CROSS_SECTIONS_PATH <- paste0("/Users/anguswatters/Desktop/cross_sections2.parquet")

bb_df <- lapply(COASTAL_BATHY_DEM_PATHS, function(path) {
                      r       <- terra::rast(path)
                      extent  <- terra::ext(r)
                      # r <- terra::project(r, "EPSG:5070")
                      # terra::set.crs(r, "EPSG:5070")
                      ext_df          <- data.frame(lapply(extent, function(i) {i}))
                      ext_df$crs      <- terra::crs(r)
                      ext_df$file     <- basename(path)
                      ext_df$path     <- path 
                      
                      return(ext_df)
              }) %>% 
  dplyr::bind_rows()
# final_cs

# INTERNAL_CROSS_SECTIONS_PATH <- paste0("/Users/anguswatters/Desktop/cross_sections_is_dem_point.parquet")
# arrow::write_parquet(
#   final_cs,
#   INTERNAL_CROSS_SECTIONS_PATH
# )

# CROSS_SECTIONS_OUTPUT_PATH <- paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/cross_sections.parquet")
# final_cs <- arrow::read_parquet(CROSS_SECTIONS_OUTPUT_PATH)
# START_PATH <- paste0("/Users/anguswatters/Desktop/cross_sections2.parquet")
# 
# arrow::write_parquet(
#   final_cs,
#   START_PATH
# )

# rm(i, cs, df, pts_subset, has_pts_in_bb, bb, updated_depths)

for (i in 1:nrow(bb_df)) {
  
  # i = 1
  
  message(i, " - Checking DEM '", bb_df$file[i], "' for CS PTS...")
  
  df <- bb_df[i, ]
  
  START_EPSG_CODE <- 5070
  
  cs <- arrow::read_parquet(START_PATH)
  cs <- 
    cs %>% 
    # dplyr::slice(1:10000) %>% 
    sf::st_as_sf(coords = c("X", "Y"), crs = START_EPSG_CODE)
  
  # convert to bounding box CRS 
  cs <- 
    cs %>% 
    sf::st_transform(df$crs[1])
  
  # create bounding box shape
  bb <- sf::st_as_sf(
          sf::st_as_sfc(
            sf::st_bbox(
                c(xmin =   df$xmin, 
                  xmax =   df$xmax, 
                  ymax =   df$ymax, 
                  ymin =   df$ymin
                  ), 
                crs = sf::st_crs(df$crs)
              )
            )
          )
  
  # get pts that are in the bounding box
  pts_subset <- sf::st_filter(
    cs, 
    bb
  )  
  
  has_pts_in_bb <- nrow(pts_subset) > 0
  
  message(nrow(pts_subset), " cs points found within '", df$file, "' DEMs bounding box")
  
  if(!has_pts_in_bb) {
    message("   > No points to update!")
    next
  }
  
  if(has_pts_in_bb) {
    
    message("   > Loading Raster")
    
    # load DEM
    dem <- terra::rast(df$path)   
    
    message("   > Extracting new cross section depth values from DEM...")
    
    updated_depths <- 
      pts_subset %>% 
      dplyr::mutate(
        Z2 = hydrofabric3D:::extract_pt_val(dem, .),
        Z_source2 = df$file
        # Z = extract_pt_val(terra::rast(dem), .)
        ) %>% 
      sf::st_drop_geometry() %>% 
      dplyr::select(
        dplyr::any_of(CROSSWALK_ID),
        cs_id, pt_id, 
        # Z,
        Z2,
        # Z_source,
        Z_source2
        )
    
    message("   > Replacing old depth values with new depth values...")
    cs <- 
      cs %>% 
      dplyr::left_join(
        updated_depths,
        by = c(CROSSWALK_ID, "cs_id", "pt_id")
      ) %>% 
      dplyr::mutate(
        Z = dplyr::case_when(
          !is.na(Z2)  ~ Z2,
          TRUE        ~ Z
        ),
        Z_source = dplyr::case_when(
          !is.na(Z_source2)  ~ Z_source2,
          TRUE               ~ Z_source
        )
      ) %>% 
      dplyr::select(-Z2, -Z_source2)
    
    message("   > Projecting CS PTs back to starting CRS (", START_EPSG_CODE, ") ...")
    cs <- 
      cs %>% 
      sf::st_transform(START_EPSG_CODE) 
    
    message("   > Dropping point geometries and preserving X / Y coordinates...")
    
    cs <- 
      cs %>% 
      dplyr::mutate(
        X = sf::st_coordinates(.)[,1],
        Y = sf::st_coordinates(.)[,2]
      ) %>%
      sf::st_drop_geometry() %>% 
      dplyr::select(
        # hy_id, 
        dplyr::any_of(CROSSWALK_ID),
        cs_id, pt_id,
        cs_lengthm,
        relative_distance,
        X, Y, 
        Z,
        class, point_type,
        Z_source,
        bottom, left_bank, right_bank, valid_banks,
        has_relief 
        
        # newly added columns (03/06/2024)
      )
    
    message("   > Overwritting original cross section points parquet with updated depth values ...")
    message("      > '", START_PATH, "'")
    
    arrow::write_parquet(cs, START_PATH)
    
  }
  
  }

cross_sections <- arrow::read_parquet(START_PATH)

cross_sections <-
  cross_sections %>% 
  dplyr::select(
          -point_type,
          -class,
          -bottom, -left_bank, -right_bank,
          -has_relief, -valid_banks
          )

# reclassify
# system.time({
cross_sections <- hydrofabric3D::classify_points(cs_pts = cross_sections, 
                                           crosswalk_id = CROSSWALK_ID,
                                           pct_of_length_for_relief = PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF
                                           )
# })

# final_cross_sections %>% 
#   dplyr::filter(id == "wb-1000") %>% 
#   dplyr::rename(hy_id = id) %>% 
#   hydrofabric3D::plot_cs_pts(color = "point_type")

# ---------------------------------------------------------------------------------
# ---- Write final cross section points data ----
# ---- Diffusive Domain DEM + FEMA + ML 
# ---------------------------------------------------------------------------------

CROSS_SECTIONS_ML_OUTPUT_PATH <- paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/cross_sections.parquet")

message(round(Sys.time()), " - Saving Diffusive DEM + FEMA + ML augmented cross section points to:\n - filepath: '", CROSS_SECTIONS_ML_OUTPUT_PATH, "'")

# sum(is.na(final_cross_sections$id))
# sum(is.na(final_cross_sections$cs_id))
# sum(is.na(final_cross_sections$pt_id))
# sum(is.na(final_cross_sections$X))

# save cross section points as a parquet to out_path (domain/outputs/cross_sections.parquet)
arrow::write_parquet(
  # dplyr::select(final_cs, 
  #               -is_dem_point
  # ), 
  cross_sections,
  CROSS_SECTIONS_ML_OUTPUT_PATH
)

# cross_sections %>% 
#   dplyr::select(id, Z_source) %>% 
#   dplyr::group_by(Z_source) %>% 
#   dplyr::count(Z_source) %>% 
#   dplyr::arrange(-n)

# CROSS_SECTIONS_ML_OUTPUT_PATH
# 
# CROSS_SECTIONS_OUTPUT_PATH
# 
# CROSS_SECTIONS_ML_OUTPUT_PATH <- paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/cross_sections.parquet")
# 
# CROSS_SECTIONS_PATH      <- paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/cross_sections.parquet")
# DEM_CROSS_SECTIONS_PATH  <- paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/dem_cross_sections.parquet")
# 
# cs <- arrow::read_parquet(CROSS_SECTIONS_PATH) %>% 
#   dplyr::rename(hy_id = id)
# 
# dem_cs <- arrow::read_parquet(DEM_CROSS_SECTIONS_PATH) %>% 
#   dplyr::rename(hy_id = id)
# 
# z_ids <- 
#   cs %>% 
#   dplyr::group_by(Z_source) %>% 
#   dplyr::slice(1:4) %>% 
#   dplyr::pull(hy_id)
# 
# # subset_ids <- z_ids[1:1]
# # z_ids <- z_ids[1:2]
# # for (i in seq_along(z_ids)) {
# for (i in seq_along(z_ids)) {
#   message(i) 
#   # i = 1
#   subset_id <- z_ids[i]
#   
#   new_cs <- 
#     cs %>% 
#     dplyr::filter(hy_id %in% subset_id)
#   
#   old_cs <- 
#     dem_cs %>% 
#     dplyr::filter(hy_id %in% subset_id)
#   
#   title_str <-  paste0(
#                 unique(new_cs$Z_source)[1], " vs. ",
#                 unique(old_cs$Z_source)[1]
#               )
#   
#   save_path = paste0(
#     "/Users/anguswatters/Desktop/z_source_plots/",
#                   gsub(".tif", "", unique(new_cs$Z_source)[1]), 
#                   "_vs_",
#                   unique(old_cs$Z_source)[1], "_", i, ".png"
#                   )
#   
#   message("making plot at \n > '", save_path, "'")
#   
#   comp_plot <- dplyr::bind_rows(
#                   new_cs,
#                   old_cs
#                 ) %>% 
#                   ggplot2::ggplot() + 
#                   ggplot2::geom_point(ggplot2::aes(x = pt_id, 
#                                                    y = Z,
#                                                    color = Z_source
#                   ),
#                   size = 3,
#                   alpha = 0.9
#                   ) + 
#                   ggplot2::labs(title = title_str) + 
#                   ggplot2::facet_wrap(hy_id~cs_id) +
#                   ggplot2::theme_bw() +
#                   ggplot2::theme(
#                     plot.title = ggplot2::element_text(size = 14),
#                     legend.text  = ggplot2::element_text(size = 14),
#                     legend.title = ggplot2::element_text(size = 14, face = "bold")
#                   )
#   
#   ggplot2::ggsave(plot = comp_plot,
#                   filename = save_path,
#                   scale = 1
#                   )
#   
#   }
# 
# 
# dplyr::bind_rows(
#   cs %>% 
#     dplyr::filter(hy_id %in% subset_ids),
#   dem_cs %>% 
#     dplyr::filter(hy_id %in% subset_ids) 
#     
# ) %>% 
#   ggplot2::ggplot() + 
#   ggplot2::geom_point(ggplot2::aes(x = pt_id, 
#                                    y = Z,
#                                    color = Z_source
#                                    
#   )) + 
#   ggplot2::labs(title = "al_nwfl_ncei_1 vs hydrofabric3D") + 
#   ggplot2::facet_wrap(hy_id~cs_id)
# 
# cs %>% 
#   dplyr::filter(hy_id %in% subset_ids) %>% 
#   ggplot2::ggplot() + 
#   ggplot2::geom_point(ggplot2::aes(x = pt_id, 
#                                    y = Z,
#                                    color = point_type
#                                    
#   )) + 
#   ggplot2::facet_wrap(hy_id~cs_id)
# 
# 
# dem_cs %>% 
#   dplyr::filter(hy_id %in% subset_ids) %>% 
#   ggplot2::ggplot() + 
#   ggplot2::geom_point(ggplot2::aes(x = pt_id, 
#                                    y = Z,
#                                    color = point_type
#                                    
#                                    )) + 
#   ggplot2::facet_wrap(hy_id~cs_id)
# 
# 
# 
# 
# 
# 













