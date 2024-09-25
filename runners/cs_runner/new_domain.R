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

# ---------------------------------------------------------------------
# --- Read in flowlines
# ---------------------------------------------------------------------
# read in nextgen data
flines <- sf::read_sf(NEW_DOMAIN_FLOWLINES_PATH, layer = "flowpaths")

VPUS           <- sf::st_transform(nhdplusTools::vpu_boundaries, sf::st_crs(flines))
vpu_flines_int <- sf::st_intersects(VPUS, flines)
vpu_aoi        <- VPUS[lengths(vpu_flines_int) > 0, ]

# VPU IDs of interest 
VPU_IDS        <- vpu_aoi$VPUID

# sf::write_sf(aoi, "all_diffusive_combined_vpus.gpkg")
# VPUS[lengths(vpu_fl_int) > 0, ] %>% mapview::mapview()
# nhdplusTools::vpu_boundaries

# flines %>% 
#   dplyr::filter(order >= 3) %>% 
# # flines %>% 
#   sf::st_transform(4326) %>% 
#   mapview::mapview()
# # ggplot2::ggplot() +
# #   ggplot2::geom_sf()
#   
# 
# flines %>% 
#   sf::st_transform(4326) %>% 
#   sf::st_bbox() %>% 
#   sf::st_as_sfc() %>% 
#   sf::st_as_sf() %>% 
#   mapview::mapview()

# ---------------------------------------------------------------------
# --- Add bankful width estimate
# ---------------------------------------------------------------------

# calculate bankfull width
flines <-
  flines %>%
  sf::st_transform(26904) %>% 
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

# ---------------------------------------------------------------------
# --- Create transect lines
# ---------------------------------------------------------------------
system.time({
  
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

})

# add cs_source column and rename cs_widths to cs_lengthm
transects <- 
  transects %>%
  dplyr::mutate(
    cs_source = CS_SOURCE
  )

# ---------------------------------------------------------------------
# --- Extract cross section points from DEM 
# ---------------------------------------------------------------------
# dem <- terra::rast(NEW_DOMAIN_DEM_PATH)

system.time({

# get cross section point elevations
cs_pts <- hydrofabric3D::cross_section_pts(
  
  cs             = transects,
  crosswalk_id   = "hy_id",
  points_per_cs  = NULL,
  min_pts_per_cs = 10,
  dem            = NEW_DOMAIN_DEM_PATH
)

})


# ---------------------------------------------------------------------
# --- Classify cross section points
# ---------------------------------------------------------------------
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

# cs_pts3 <-
#   cs_pts %>% 
#   # dplyr::group_by(hy_id, cs_id) %>% 
#   # dplyr::filter(!any(is.na(Z))) %>% 
#   # dplyr::ungroup() %>% 
#   hydrofabric3D::drop_incomplete_cs_pts("hy_id") %>% 
#   hydrofabric3D::classify_points2(
#     crosswalk_id             = "hy_id", 
#     pct_of_length_for_relief = PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF
#   )  
# 
# old_plot <- 
#   cs_pts2 %>% 
#   dplyr::slice(1:100) %>% 
#   hydrofabric3D::plot_cs_pts(color = 'point_type', size = 4)+
#   ggplot2::labs(title = "OLD CLASSES")
# 
# new_plot <- 
#   cs_pts3 %>% 
#   dplyr::slice(1:100) %>% 
#   hydrofabric3D::plot_cs_pts(color = 'point_type', size = 4) +
#   ggplot2::labs(title = "NEW CLASSES")
# 
# library(patchwork)
# old_plot / new_plot
# })

ids_original_cs_pts <- hydrofabric3D::add_tmp_id(cs_pts)$tmp_id  

# ----------------------------------------------------------------------------------------------------------------
# ---- Improve cross section points based on the bank validity & amount of relief in a cross section
# ----------------------------------------------------------------------------------------------------------------

# system.time({
fixed_pts <- hydrofabric3D::get_improved_cs_pts(
  cs_pts         = cs_pts,    # cross section points generated from hydrofabric3D::cross_section_pts()
  net            = flines,     # original flowline network
  transects      = transects, # original transect lines
  crosswalk_id   = "hy_id",
  points_per_cs  = NULL, 
  min_pts_per_cs = 10, # number of points per cross sections
  dem            = NEW_DOMAIN_DEM_PATH, # DEM to extract points from
  scale          = EXTENSION_PCT, # How far to extend transects if the points need to be rechecked
  pct_of_length_for_relief = PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF, # percent of cross sections length to be needed in relief calculation to consider cross section to "have relief"
  fix_ids = FALSE,
  verbose = TRUE
)
# })

# ----------------------------------------------------------------------------------------------------------------
# ---- Update transects with extended transects (if exists) ----
# ----------------------------------------------------------------------------------------------------------------

out_transects <- match_transects_to_extended_cs_pts(
                      transect_lines = transects, 
                      fixed_cs_pts   = fixed_pts, 
                      crosswalk_id   = "hy_id"
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

fixed_pts      <- hydrofabric3D:::renumber_cs_ids(df = fixed_pts, crosswalk_id = "hy_id")
out_transects  <- hydrofabric3D:::renumber_cs_ids(df = out_transects, crosswalk_id = "hy_id")

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
    hy_id, 
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
  dplyr::relocate(hy_id, cs_id, pt_id, cs_lengthm, relative_distance, X, Y, Z, Z_source, 
                  class, point_type, 
                  bottom, left_bank, right_bank, valid_banks, has_relief)

ids_before_align <- hydrofabric3D::add_tmp_id(fixed_pts)$tmp_id

message("Aligning banks and smoothing bottoms...")
fixed_pts <- hydrofabric3D::align_banks_and_bottoms(fixed_pts)

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

sf::write_sf(
  out_transects,
  paste0(NEW_DOMAIN_TRANSECTS_DIR, "/hi_transects.gpkg"), 
)

arrow::write_parquet(
  fixed_pts,
  paste0(NEW_DOMAIN_CS_PTS_DIR, "/hi_cs_pts.parquet"), 
)

# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------
# ---- STEP 4: Update transects with extended transects (if exists) ----
# ----------------------------------------------------------------------------------------------------------------

# get the counts of each point type to add this data to the transects dataset
point_type_counts <- hydrofabric3D::get_point_type_counts(fixed_pts)

# # check the number of cross sections that were extended
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
renumb <- hydrofabric3D:::renumber_cs_ids(fixed_pts, crosswalk_id = "hy_id")
unumb <- renumb %>% get_unique_tmp_ids()
nnumb <- fixed_pts %>% get_unique_tmp_ids()
all(unumb %in% nnumb)
all(nnumb %in% unumb)


length(unumb)
length(nnumb)

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
out_transects2 <- dplyr::left_join(
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
fixed_pts <- 
  dplyr::left_join(
      hydrofabric3D::add_tmp_id(fixed_pts),
      renumbered_ids,
      by = "tmp_id"
  ) %>% 
  dplyr::select(-cs_id, -tmp_id) %>% 
  dplyr::rename(cs_id = new_cs_id)

renumbered_ids <- 
  df %>%
  sf::st_drop_geometry() %>%
  dplyr::select(
    # hy_id, 
    dplyr::any_of(crosswalk_id),
    cs_id, pt_id, cs_measure
  ) %>%
  dplyr::group_by(dplyr::across(dplyr::any_of(c(crosswalk_id, "cs_id")))) %>% 
  # dplyr::group_by(hy_id, cs_id) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(dplyr::across(dplyr::any_of(c(crosswalk_id)))) %>% 
  # dplyr::group_by(hy_id) %>%
  dplyr::mutate(
    new_cs_id = 1:dplyr::n()
    # tmp_id = paste0(hy_id, "_", cs_id)
  ) %>%
  add_tmp_id(x = get(crosswalk_id)) %>% 
  dplyr::ungroup() %>%
  dplyr::select(new_cs_id, tmp_id)

# Join the new cs_ids back with the final output data to replace the old cs_ids
df <- dplyr::left_join(
  add_tmp_id(df, x = get(crosswalk_id)),
  # dplyr::mutate(df,tmp_id = paste0(hy_id, "_", cs_id)),
  renumbered_ids,
  by = "tmp_id"
) %>%
  dplyr::select(-cs_id, -tmp_id) %>%
  dplyr::rename("cs_id" = "new_cs_id") %>%
  dplyr::relocate(dplyr::any_of(crosswalk_id), cs_id)
# dplyr::relocate(hy_id, cs_id)

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















