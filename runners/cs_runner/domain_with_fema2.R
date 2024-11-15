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

# read in flowlines based on IDs in AOI
flines <- sf::read_sf(DOMAIN_WITH_FEMA_FLOWLINES_PATH, layer = "flowpaths",
                      wkt_filter = wkt
                      # query = sprintf("SELECT * FROM \"%s\" WHERE %s IN (%s)", layer, id_col, query_ids),
                      # query = query
                      )

# flines <- 
#   flines %>% 
#   dplyr::slice(1:1000)

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

# rm(flines2)
# unnest_ids(flines2$VPUID2)

# set of unique VPUs

VPU_IDS <- unnest_ids(flines$VPUID)
# VPU_IDS <- unnest_ids(flines2$VPUID2)
VPU_IDS

# all possible FEMA dirs
AOI_FEMA_DIRS <- FEMA_VPU_SUBFOLDERS[basename(FEMA_VPU_SUBFOLDERS) %in% paste0("VPU_", VPU_IDS) ]

flines %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(VPUID) %>%
  dplyr::count()

flines %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(order) %>%
  dplyr::count()
# unique(flines$VPUID)

# calculate bankfull width
flines <- 
  flines %>% 
  hydrofabric3D::add_powerlaw_bankful_width("tot_drainage_areasqkm", 50) %>% 
  dplyr::select(
    dplyr::any_of(CROSSWALK_ID),
    VPUID,
    # hy_id = id,
    lengthkm,
    tot_drainage_areasqkm,
    bf_width,
    order,
    mainstem
  ) %>% 
  hydroloom::rename_geometry("geometry")


sf::write_sf(
  flines,
  paste0(DOMAIN_WITH_FEMA_FLOWLINES_DIR, "/flowlines_subset.gpkg")
)

transects <- hydrofabric3D::cut_cross_sections(
  net               = flines,                        # flowlines network
  crosswalk_id      = CROSSWALK_ID,                     # Unique feature ID
  cs_widths         = flines$bf_width,
  # cs_widths         = pmax(50, flowlines$bf_width * 11),     # cross section width of each "id" linestring ("hy_id")
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

transects <- 
  transects %>% 
  dplyr::mutate(cs_source = CS_SOURCE) %>% 
  dplyr::select(id,
                cs_source,
                cs_id,
                cs_measure,
                ds_distance,
                cs_lengthm,
                sinuosity,
                geometry
  )

# paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/pre_extension_transects.gpkg")

sf::write_sf(
  transects,
  paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/pre-extension-transects.gpkg")
)

# [1] "wb-2416492" "wb-2414869" "wb-2415765" "wb-2416657" "wb-2415479" "wb-2420307" "wb-2425952" "wb-2421269" "wb-2425750"
# [10] "wb-2422028" "wb-2419212" "wb-2426092" "wb-2414971" "wb-2416120" "wb-2414960" "wb-2433259" "wb-14381"   "wb-18022"  
# [19] "wb-15173"   "wb-15160"   "wb-571"     "wb-2413739" "wb-508015"  "wb-507770"  "wb-2424752" "wb-2424545" "wb-2422849"
# [28] "wb-2417939" "wb-2425090" "wb-2422212" "wb-2422714" "wb-2419650" "wb-2419650" "wb-2416492" "wb-2414869" "wb-2415765"
# [37] "wb-2416657" "wb-2415479" "wb-2420307" "wb-2419694" "wb-2419694" "wb-2419974" "wb-2419974" "wb-2420065" "wb-2420065"
# [46] "wb-2425952" "wb-2420678" "wb-2420678" "wb-2421269" "wb-2419695" "wb-2419695" "wb-2425750" "wb-2422028" "wb-2419212"
# [55] "wb-2426091" "wb-2414971" "wb-2407999" "wb-2407999" "wb-2426131" "wb-2426131" "wb-2416120" "wb-2433525" "wb-2433525"
# [64] "wb-2433259" "wb-14381"   "wb-14576"   "wb-14576"   "wb-18022"   "wb-15173"   "wb-15160"   "wb-11097"   "wb-11097"  
# [73] "wb-14804"   "wb-14804"   "wb-571"     "wb-2408123" "wb-2408123" "wb-2413739" "wb-508015"  "wb-507770"  "wb-2424752"
# [82] "wb-2435791" "wb-2435791" "wb-2435798" "wb-2435798" "wb-2435779" "wb-2435779" "wb-2435781" "wb-2435781" "wb-2424545"
# [91] "wb-2422849" "wb-2417939" "wb-2425090" "wb-2422212" "wb-2422714"

# transects <- sf::read_sf(paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/pre-extension-transects.gpkg"))

# out_transects
# out_transects2

transects <- 
  transects %>% 
  dplyr::left_join(
    flines %>% 
      dplyr::select(id, VPUID, mainstem, tot_drainage_areasqkm) %>% 
      sf::st_drop_geometry(),
    by = "id"
  )

# sf::write_sf(
#   transects,
#   paste0("/Users/anguswatters/Desktop/tmp_trans.gpkg")
# )

GROUP_VPU_IDS <- unnest_ids(transects$VPUID)

# all FEMA dirs for the current area
GROUP_FEMA_DIRS  <- FEMA_VPU_SUBFOLDERS[basename(FEMA_VPU_SUBFOLDERS) %in% paste0("VPU_", GROUP_VPU_IDS) ]
GROUP_FEMA_FILES <- list.files(GROUP_FEMA_DIRS, full.names = T)[grepl("_output.gpkg", list.files(GROUP_FEMA_DIRS, full.names = T))]


# read each FEMA geopackage into a list 
fema <- lapply(GROUP_FEMA_FILES, function(gpkg) sf::read_sf(gpkg))
# transects$VPUID %>% unique()

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

# sf::write_sf(
#   fema,
#   paste0("/Users/anguswatters/Desktop/tmp_fema.gpkg")
# )

message(" - Number of geoms AFTER simplifying: ", nrow(fema))
message("Extending transects out to FEMA 100yr floodplain polygon boundaries - (", Sys.time(), ")")

# transects <- 
#   transects  %>%
#   dplyr::left_join(
#     dplyr::select(sf::st_drop_geometry(flines),
#                   dplyr::any_of(CROSSWALK_ID),
#                   mainstem
#     ),
#     by = CROSSWALK_ID
#   )

# test_ids <- c("wb-2416492", "wb-2414869", "wb-2415765" , "wb-2416657" , "wb-2415479",
#               "wb-2420307" , "wb-2425952" , "wb-2421269" , "wb-2425750",
#               "wb-14804",   "wb-14804",   "wb-571",     "wb-2408123", "wb-2408123"
#               )
# test_trans <- 
#   transects %>% 
#   dplyr::filter(id %in% test_ids)
# 
# test_flines <- 
#   flines %>% 
#   dplyr::filter(id %in% test_ids)

# test_fema
# test_trans$VPUID

# TODO: make sure this 3000m extension distance is appropriate across VPUs 
# TODO: also got to make sure that this will be feasible on memory on the larger VPUs...
# transects <- hydrofabric3D::extend_transects_to_polygons(
ext_transects <- hydrofabric3D::extend_transects_to_polygons(
  transect_lines         = transects,
  polygons               = fema,
  flowlines              = flines,
  # transect_lines         = test_trans, 
  # polygons               = fema, 
  # flowlines              = test_flines,
  crosswalk_id           = CROSSWALK_ID,
  grouping_id            = "mainstem", 
  max_extension_distance = 3000 
)

# sf::write_sf(
#   flines,
#   paste0("/Users/anguswatters/Desktop/tmp_flines.gpkg")
# )

# flines <- sf::read_sf( paste0("/Users/anguswatters/Desktop/tmp_flines.gpkg"))
# transects <- sf::read_sf(paste0("/Users/anguswatters/Desktop/tmp_trans.gpkg"))
# fema <- sf::read_sf( paste0("/Users/anguswatters/Desktop/tmp_fema.gpkg"))
# 
# transect_lines <- transects
# polygons <- fema
# flowlines <- flines
# bad_id   <- "wb-10813"
# mainstem <- "1977479"
# crosswalk_id           = "id"
# grouping_id            = "mainstem"
# max_extension_distance = 3000 


# hydrofabric3D:::rm_self_intersections(ext_transects)
# transects <- 
#   transects %>% 
ext_transects <-
  ext_transects %>%
  dplyr::mutate(cs_source = CS_SOURCE) %>% 
  dplyr::select(id,
                cs_source,
                cs_id,
                cs_measure,
                ds_distance,
                cs_lengthm,
                sinuosity,
                geometry
                )

# ext_transects <-
#   ext_transects %>% 
#   dplyr::select(id, cs_id, cs_lengthm, cs_measure, ds_distance, sinuosity, 
#                 # tot_drainage_areasqkm, 
#                 geometry)

sf::write_sf(
  # transects,
  ext_transects,
  paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/post-extension-transects.gpkg")
)

# get cross section point elevations
cs_pts <- hydrofabric3D::cross_section_pts(
  # cs             = transects,
  cs             = ext_transects,
  crosswalk_id   = CROSSWALK_ID,
  points_per_cs  = NULL,
  min_pts_per_cs = 10,
  dem            = DEM_PATH
)

# mapview::mapview(cs_pts) + test_trans

sf::write_sf(
  cs_pts,
  paste0(DOMAIN_WITH_FEMA_CS_PTS_DIR, "/raw-cs-pts.gpkg")
)

# DOMAIN_WITH_FEMA_CS_PTS_DIR
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

# cs_pts %>% 
#   dplyr::slice(150:250) %>% 
#   hydrofabric3D::plot_cs_pts("id", color = "point_type", size = 4)

# })

ids_original_cs_pts <- hydrofabric3D::add_tmp_id(cs_pts, x = CROSSWALK_ID)$tmp_id  
# ids_original_cs_pts <- hydrofabric3D::add_tmp_id(cs_pts2)$tmp_id

# sf::write_sf(cs_pts2, "/Users/anguswatters/Desktop/test_improve_cs_pts_classified_11.gpkg")
# sf::write_sf(cs_pts, "/Users/anguswatters/Desktop/test_improve_cs_pts_classified_11_2.gpkg")


# ----------------------------------------------------------------------------------------------------------------
# ---- STEP 3: Try to rectify any no relief and invalid banks cross sections ----
# ----------------------------------------------------------------------------------------------------------------

system.time({
fixed_pts <- hydrofabric3D::get_improved_cs_pts(
  cs_pts         = cs_pts,    # cross section points generated from hydrofabric3D::cross_section_pts()
  net            = flines,    # original flowline network
  # net = test_flines,
  # net            = flowlines,    # original flowline network
  # transects      = transects, # original transect lines
  transects      = ext_transects, # original transect lines
  crosswalk_id   = CROSSWALK_ID,
  points_per_cs  = NULL, 
  min_pts_per_cs = 10, # number of points per cross sections
  dem            = DEM_PATH, # DEM to extract points from
  scale          = EXTENSION_PCT, # How far to extend transects if the points need to be rechecked
  pct_of_length_for_relief = PCT_LENGTH_OF_CROSS_SECTION_FOR_RELIEF, # percent of cross sections length to be needed in relief calculation to consider cross section to "have relief"
  fix_ids = FALSE,
  verbose = TRUE
)
})

# sf::write_sf(
#   test_flines,
#   "/Users/anguswatters/Desktop/test_flines.gpkg"
# )
# 
# sf::write_sf(
#   test_trans,
#   "/Users/anguswatters/Desktop/test_start_trans.gpkg"
# )
# 
# sf::write_sf(
#   ext_transects,
#   "/Users/anguswatters/Desktop/test_ext_trans.gpkg"
# )
# 
# sf::write_sf(
#   cs_pts,
#   "/Users/anguswatters/Desktop/test_cs_pts.gpkg"
# )

# mapview::mapview(fixed_pts) +
#   mapview::mapview(ext_transects, color = "green") +
#   mapview::mapview(out_transects, color = "red")

fixed_pts$Z %>% is.null() %>% any()
fixed_pts$Z %>% is.na() %>% any()

ids_after_fixed_pts <- hydrofabric3D::add_tmp_id(fixed_pts, x = CROSSWALK_ID)$tmp_id  

# ----------------------------------------------------------------------------------------------------------------
# ---- Update transects with extended transects (if exists) ----
# ----------------------------------------------------------------------------------------------------------------
# ext_transects <- 
#   ext_transects %>% 
#   dplyr::mutate(cs_source = CS_SOURCE)

# end_pts <- 
#   fixed_pts %>% 
#   dplyr::filter(id == "wb-2435034", cs_id == 1) %>% 
#   # hydrofabric3D::plot_cs_pts("id")
#   dplyr::slice(
#     which.min(pt_id),
#     which.max(pt_id)
#   )
# mapview::mapview(end_pts) + tt
# tt <- 
#   ext_transects %>% 
#   dplyr::filter(id == "wb-2435034", cs_id == 1)
# hydrofabric3D:::match_transects_to_extended_cs_pts()

out_transects <- hydrofabric3D:::match_transects_to_extended_cs_pts(
  # transect_lines = transects, 
  transect_lines = ext_transects,
  fixed_cs_pts   = fixed_pts, 
  crosswalk_id   = CROSSWALK_ID
  # extension_pct  = EXTENSION_PCT
)

# sf::write_sf(
#   out_transects,
#   paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/cs-point-extension-transects.gpkg")
#   # "/Users/anguswatters/Desktop/tmp_trans_improved.gpkg"
# )

# fixed_pts2      <- hydrofabric3D:::renumber_cs_ids(df = fixed_pts, crosswalk_id = CROSSWALK_ID)
# out_transects2  <- hydrofabric3D:::renumber_cs_ids(df = out_transects, crosswalk_id = CROSSWALK_ID)

fixed_pts <- 
  fixed_pts %>% 
# fixed_pts2 <- 
  # fixed_pts %>% 
  hydrofabric3D:::renumber_cs_ids(crosswalk_id = CROSSWALK_ID) %>% 
  dplyr::group_by(id) %>% 
  dplyr::arrange(cs_id, .by_group = TRUE) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(id, cs_id) %>% 
  dplyr::arrange(pt_id, .by_group = TRUE) %>% 
  dplyr::ungroup()


out_transects <- 
  out_transects %>% 
# out_transects2 <- 
#   out_transects %>% 
  hydrofabric3D:::renumber_cs_ids(crosswalk_id = CROSSWALK_ID) %>% 
  dplyr::group_by(id) %>% 
  dplyr::arrange(cs_id, .by_group = TRUE) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(
    -left_bank_count,
    -right_bank_count, 
    -channel_count,
    -bottom_count,
    -bottom,
    -left_bank,
    -right_bank,
    -valid_banks,
    -has_relief,
    -is_extended
  ) 

sf::write_sf(
  out_transects,
  # out_transects2,
  paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/cs-point-extension-transects.gpkg")
  # "/Users/anguswatters/Desktop/tmp_trans_improved.gpkg"
)

# sf::write_sf(
#   fixed_pts2,
#   paste0(DOMAIN_WITH_FEMA_CS_PTS_DIR, "/fixed_cs_pts.gpkg")
# )

# transects %>% 
#   dplyr::group_by(id) %>% 
#   dplyr::arrange(cs_id, .by_group = TRUE) %>% 
#   dplyr::ungroup()
# 
# out_transects2 %>% 
#   # dplyr::filter(id %in% c("wb-2435034", "wb-2435045")) %>% 
#   dplyr::group_by(id) %>% 
#   dplyr::arrange(cs_id, .by_group = TRUE) %>% 
#   dplyr::ungroup()

# sf::write_sf(
#   out_transects2,
#   paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/cs-point-extension-transects.gpkg")
#   # "/Users/anguswatters/Desktop/tmp_trans_improved.gpkg"
# )
# 
# sf::write_sf(
#   fixed_pts2,
#   paste0(DOMAIN_WITH_FEMA_CS_PTS_DIR, "/fixed_cs_pts.gpkg")
# )

# trans_uids <- hydrofabric3D::get_unique_tmp_ids(out_transects, "id")
# ext_trans_uids <- hydrofabric3D::get_unique_tmp_ids(ext_transects, "id")
# fixed_pts_uids <- hydrofabric3D::get_unique_tmp_ids(fixed_pts, "id")

# all(trans_uids %in% ext_trans_uids)
# all(ext_trans_uids %in% trans_uids)
# all(trans_uids %in% fixed_pts_uids)
# all(fixed_pts_uids %in% trans_uids)
# all(ext_trans_uids %in% fixed_pts_uids)
# all(fixed_pts_uids %in% ext_trans_uids)

# ----------------------------------------------------------------------------------------------------------------
# ---- Re enumerate the transects & cross section points "cs_id" ----
# ----------------------------------------------------------------------------------------------------------------

# fixed_pts      <- hydrofabric3D:::renumber_cs_ids(df = fixed_pts, crosswalk_id = "hy_id")
# out_transects  <- hydrofabric3D:::renumber_cs_ids(
#                                                   df            = dplyr::mutate(out_transects, pt_id = 1), 
#                                                   crosswalk_id  = "hy_id"
#                                                   ) %>% 
#                                                   dplyr::select(-pt_id)
# length(out_transects_uids)
# length(fixed_pts_uids)
# 
# out_transects_uids <- hydrofabric3D::get_unique_tmp_ids(out_transects, "id")
# fixed_pts_uids <- hydrofabric3D::get_unique_tmp_ids(fixed_pts, "id")
# 
# length(out_transects_uids) == length(fixed_pts_uids) & all(out_transects_uids %in% fixed_pts_uids) & all(fixed_pts_uids %in% out_transects_uids)
# 
# fixed_pts2      <- hydrofabric3D:::renumber_cs_ids(df = fixed_pts, crosswalk_id = CROSSWALK_ID)
# out_transects2  <- hydrofabric3D:::renumber_cs_ids(df = out_transects, crosswalk_id = CROSSWALK_ID)
# 
# out_transects2_uids <- hydrofabric3D::get_unique_tmp_ids(out_transects2, "id")
# fixed_pts2_uids <- hydrofabric3D::get_unique_tmp_ids(fixed_pts2, "id")
# 
# length(out_transects2_uids) == length(fixed_pts2_uids) & all(out_transects2_uids %in% fixed_pts2_uids) & all(fixed_pts2_uids %in% out_transects2_uids)

  # paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/extended_transects.gpkg")
  
# sf::write_sf(
#   out_transects2,
#   # "/Users/anguswatters/Desktop/tmp_trans_improved2.gpkg"
#   paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/final_transects.gpkg")
# )

# object.size(out_transects2)

# classify the cross section points
fixed_pts <-
  fixed_pts %>% 
# fixed_pts2 <-
  # fixed_pts2 %>% 
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
  ) %>% 
  dplyr::mutate(
    Z_source = CS_SOURCE
  ) %>%
  dplyr::relocate(
    dplyr::any_of(CROSSWALK_ID),
    cs_id, pt_id, cs_lengthm, relative_distance, X, Y, Z, Z_source,
    class, point_type,
    bottom, left_bank, right_bank, valid_banks, has_relief
  )

# # add Z_source column for source of elevation data
# fixed_pts2 <-
#   fixed_pts2 %>%
#   dplyr::mutate(
#     Z_source = CS_SOURCE
#   ) %>%
#   dplyr::relocate(
#     dplyr::any_of(CROSSWALK_ID),
#     cs_id, pt_id, cs_lengthm, relative_distance, X, Y, Z, Z_source
#     # class, point_type, 
#     # bottom, left_bank, right_bank, valid_banks, has_relief
#     )

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

# sf::write_sf(
#   out_transects,
#   paste0(DOMAIN_WITH_FEMA_VPU_SUBSETS_DIR, "/", VPU, "_transects.gpkg") 
# )
# 

arrow::write_parquet(
  fixed_pts,
  paste0(DOMAIN_WITH_FEMA_CS_PTS_DIR, "/final-cs-pts.parquet")
)
# paste0(DOMAIN_WITH_FEMA_CS_PTS_DIR, "/final-cs-pts.parquet")


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

# CS_PTS_OUTPUT_PATH <- paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/cs_pts.parquet")
CS_PTS_OUTPUT_PATH <- paste0(DOMAIN_WITH_FEMA_CS_PTS_DIR, "/final-cs-pts.parquet")

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

starting_uids  <- hydrofabric3D::get_unique_tmp_ids(cs_pts, x = CROSSWALK_ID)
ending_uids    <- hydrofabric3D::get_unique_tmp_ids(final_cs, x = CROSSWALK_ID)

has_same_number_of_uids           <- length(starting_uids) == length(ending_uids)
all_starting_uids_in_ending_uids  <- all(starting_uids %in% ending_uids)
all_ending_uids_in_starting_uids  <- all(ending_uids %in% starting_uids)

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
# # save cross section points as a parquet to out_path (domain/outputs/cross-sections.parquet)
# arrow::write_parquet(
#   dplyr::select(final_cs,
#                 -is_dem_point
#   ),
#   paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/final_cs2.parquet")
# )

# ---------------------------------------------------------------------------------
# ---- Write final cross section points data ----
# ---------------------------------------------------------------------------------

# CROSS_SECTIONS_OUTPUT_PATH <- paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/cross-sections.parquet")
CROSS_SECTIONS_OUTPUT_PATH <- paste0(DOMAIN_WITH_FEMA_CROSS_SECTIONS_DIR,  "/cross-sections.parquet")
message(round(Sys.time()), " - Saving ML augmented cross section points to:\n - filepath: '", CROSS_SECTIONS_OUTPUT_PATH, "'")

# save cross section points as a parquet to out_path (domain/outputs/cross-sections.parquet)
arrow::write_parquet(
  dplyr::select(final_cs, 
                -is_dem_point
  ), 
  CROSS_SECTIONS_OUTPUT_PATH
)

INTERNAL_CROSS_SECTIONS_PATH <- paste0(DOMAIN_WITH_FEMA_CROSS_SECTIONS_DIR, "/cross-sections-is-dem-point.parquet")

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

# INTERNAL_CROSS_SECTIONS_PATH <- paste0("/Users/anguswatters/Desktop/cross-sections-is-dem-point.parquet")
# arrow::write_parquet(
#   final_cs,
#   INTERNAL_CROSS_SECTIONS_PATH
# )

# CROSS_SECTIONS_OUTPUT_PATH <- paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/cross-sections.parquet")
# final_cs <- arrow::read_parquet(CROSS_SECTIONS_OUTPUT_PATH)
# START_PATH <- paste0("/Users/anguswatters/Desktop/cross_sections2.parquet")
# 
# arrow::write_parquet(
#   final_cs,
#   START_PATH
# )

# rm(i, count, cs, df, tmp, aoi, new_tmp, pts_subset, dem, has_pts_in_bb, bb, updated_depths)


# CROSS_SECTIONS_OUTPUT_PATH
# START_CROSS_SECTIONS_OUTPUT_PATH <- paste0(DOMAIN_WITH_FEMA_CROSS_SECTIONS_DIR,  "/cross-sections.parquet")
START_PATH    <-  paste0(DOMAIN_WITH_FEMA_CROSS_SECTIONS_DIR,  "/cross-sections.parquet")
UPDATED_PATH  <-  paste0(DOMAIN_WITH_FEMA_CROSS_SECTIONS_DIR,  "/coastal-bathy_cross-sections.parquet")

count <- 0

for (i in 1:nrow(bb_df)) {
  
  # i = 1
  
  is_first_iter <- count == 0
  count         <- count + 1
  CURR_PATH     <- ifelse(is_first_iter, START_PATH, UPDATED_PATH)
  
  # CURR_PATH
  
  message(i, " - Checking DEM '", bb_df$file[i], "' for CS PTS...")
  
  df <- bb_df[i, ]
  
  START_EPSG_CODE <- 5070
  
  cs <- arrow::read_parquet(CURR_PATH)
  
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
    message("   > Overwritting original cross section points parquet with updated depth values ...")
    message("      > '", UPDATED_PATH, "'")
    
    cs <- 
      cs %>% 
      sf::st_transform(START_EPSG_CODE) %>% 
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
      )
    
    arrow::write_parquet(cs, UPDATED_PATH)
    
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
    
    
    # cs %>% 
    #   dplyr::slice(1:27) %>% 
    #   hydrofabric3D::plot_cs_pts(crosswalk_id = "id", color = "point_type", size = 4)
    # hydrofabric3D::classify_points(crosswalk_id = "id")
    # hydrofabric3D::plot_cs_pts(crosswalk_id = "id", color = "point_type", size = 4)
      
    # cs <- 
    #   cs %>% 
    #   # dplyr::slice(1:27) %>% 
    #   dplyr::select(-class, -point_type, 
    #                 -valid_banks, -has_relief, 
    #                 -bottom, -left_bank, -right_bank
    #                 ) %>% 
    #   hydrofabric3D::classify_points(crosswalk_id = "id") %>% 
    #   dplyr::select(
    #     # hy_id, 
    #     dplyr::any_of(CROSSWALK_ID),
    #     cs_id, pt_id,
    #     Z,
    #     relative_distance,
    #     cs_lengthm,
    #     class,
    #     point_type,
    #     X, Y, 
    #     Z_source,
    #     bottom, left_bank, right_bank, 
    #     valid_banks, has_relief 
    #     # newly added columns (03/06/2024)
    #   )
      # hydrofabric3D::plot_cs_pts(crosswalk_id = "id", color = "point_type", size = 4)
      
    # tmp <- final_cs %>% dplyr::filter(id == "wb-507785", cs_id == 4)
    # new_tmp <- cs %>% dplyr::filter(id == "wb-507785", cs_id == 4)
    
    # hydrofabric3D::plot_cs_pts(tmp, crosswalk_id = "id", color = "point_type", size = 4)
    # hydrofabric3D::plot_cs_pts(new_tmp, crosswalk_id = "id", color = "point_type", size = 4)
    
    message("   > Overwritting original cross section points parquet with updated depth values ...")
    message("      > '", UPDATED_PATH, "'")
    
    arrow::write_parquet(cs, UPDATED_PATH)
    
  }
  
}

cross_sections <- arrow::read_parquet(UPDATED_PATH)

# cross_sections %>% 
#   dplyr::select(Z_source) %>% 
#   dplyr::group_by(Z_source) %>% 
#   dplyr::count()

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

cross_sections <- 
  cross_sections %>% 
  dplyr::select(
    # hy_id, 
    dplyr::any_of(CROSSWALK_ID),
    cs_id, pt_id,
    Z,
    relative_distance,
    cs_lengthm,
    class,
    point_type,
    X, Y, 
    Z_source,
    bottom, left_bank, right_bank, 
    valid_banks, has_relief 
    # newly added columns (03/06/2024)
    )
# })

# START_PATH    <-  paste0(DOMAIN_WITH_FEMA_CROSS_SECTIONS_DIR,  "/cross-sections.parquet")
# UPDATED_PATH  <-  paste0(DOMAIN_WITH_FEMA_CROSS_SECTIONS_DIR,  "/coastal-bathy_cross-sections.parquet")

arrow::write_parquet(
  cross_sections, 
  paste0(DOMAIN_WITH_FEMA_CROSS_SECTIONS_DIR,  "/coastal-bathy_cross-sections.parquet")
  )

# -------------------------------------------------
# ---- Move final datasets to outputs/ ----
# -------------------------------------------------

select_cs_pts <- function(cs_pts, crosswalk_id = NULL) {
  
  if(is.null(crosswalk_id)) {
    # crosswalk_id  <- 'hydrofabric_id'
    stop("Please provide a valid 'crosswalk_id' which uniquely identifies each cross section in 'cs_pts'")
  }
  
  cs_pts <- 
    cs_pts %>% 
    dplyr::select(
    dplyr::any_of(c(
                    crosswalk_id,
                    "cs_id",
                    "pt_id",
                    "relative_distance",
                    "cs_lengthm",
                    "X", 
                    "Y",
                    "Z",
                    "Z_source",
                    "class", 
                    "point_type",
                    "valid_banks",
                    "has_relief"
                  )
                  )
    )
  
  return(cs_pts)
}

select_transects <- function(transects, crosswalk_id = NULL) {
  
  if(is.null(crosswalk_id)) {
    # crosswalk_id  <- 'hydrofabric_id'
    stop("Please provide a valid 'crosswalk_id' which uniquely identifies the flowline associated with each transect in 'transects'")
  }
  
  transects <- hydroloom::rename_geometry(transects, "geometry")
  
  transects <- 
    transects %>% 
    dplyr::select(
      dplyr::any_of(c(
        crosswalk_id,
        "cs_source",
        "cs_id",
        "cs_measure",
        "cs_lengthm",
        "geometry"
      )
      )
    )
  
  return(transects)
}


# arrow::read_parquet(paste0(DOMAIN_WITH_FEMA_CS_PTS_DIR, "/final-cs-pts.parquet")) %>% 
#   select_cs_pts("id")

# DEM based cross section points parquet
# - FEMA extended
# - CS based extensions
arrow::write_parquet(
  arrow::read_parquet(paste0(DOMAIN_WITH_FEMA_CS_PTS_DIR, "/final-cs-pts.parquet")) %>% 
    select_cs_pts("id"),
  # arrow::read_parquet(paste0(DOMAIN_WITH_FEMA_CS_PTS_DIR, "/final-cs-pts.parquet")),
  paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/dem_cross-sections.parquet")
)

# DEM based transect geometry gpkg 
# - FEMA extended
# - CS based extensions
sf::write_sf(
  sf::read_sf(paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/cs-point-extension-transects.gpkg")) %>% 
    select_transects("id"),
    # dplyr::rename(geometry = geom),
    paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/transects.gpkg")
)


# DEM + ML cross section points parquet 
# - FEMA extended
# - CS based extensions
# - ML based Z value updates
arrow::write_parquet(
  arrow::read_parquet(paste0(DOMAIN_WITH_FEMA_CROSS_SECTIONS_DIR, "/cross-sections.parquet")) %>% 
    select_cs_pts("id"),
  paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/ml_cross-sections.parquet")
)

# DEM + COASTAL BATHY cross section points parquet 
# - FEMA extended
# - CS based extensions
# - ML based Z value updates
# - Coastal bathy DEM Z value updates (where applicable)
arrow::write_parquet(
  arrow::read_parquet(paste0(DOMAIN_WITH_FEMA_CROSS_SECTIONS_DIR, "/coastal-bathy_cross-sections.parquet")) %>% 
    select_cs_pts("id"),
  paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/coastal-bathy_cross-sections.parquet")
)

# -------------------------------------------------
# ---- Data validation ----
# -------------------------------------------------

# paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/cs-point-extension-transects.gpkg")
# new_trans_path <- paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/transects.gpkg")

# sf::st_layers(paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/cs-point-extension-transects.gpkg"))
# sf::st_layers(paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/transects.gpkg"))

# og_transects  <- sf::read_sf(paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/cs-point-extension-transects.gpkg"))
# new_transects <- sf::read_sf(paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/transects.gpkg"))

# DOMAIN_WITH_FEMA_CROSS_SECTIONS_DIR

output_files <- list.files(DOMAIN_WITH_FEMA_OUTPUT_DIR, full.names = T)
# rm(transects, dropped_trans, cs_pts, flowlines, lengths_check, transect_lengths, rel_dist_check, wrong_lengths, mismatches, duplicate_ids, cs_pt_lengths, cs_pts_ids)

tmp_ids_list <- lapply(output_files, function(i) {
    is_gpkg    <- endsWith(i, ".gpkg")
    is_parquet <- endsWith(i, ".parquet")
    
    if(is_gpkg) {
      x <- 
        i %>% 
        sf::read_sf() %>% 
        hydrofabric3D::get_unique_tmp_ids(x = "id")
      return(x)
    }
    
    if(is_parquet) {
      x <- 
        i %>% 
        arrow::read_parquet() %>% 
        hydrofabric3D::get_unique_tmp_ids(x = "id")
      return(x)
    }
  return(NULL)
    
    
  })


tmp_ids_list

# tmp_ids_list[[1]] %in% tmp_ids_list[[2]]

res <- list()
idxs <- seq_along(tmp_ids_list)

for (i in seq_along(tmp_ids_list)) {
  # idxs <- seq_along(tmp_ids_list)
  # i = 1
  curr <- tmp_ids_list[[i]]
  
  other_idxs <- idxs[idxs != i]
  
  in_all_others <- lapply(tmp_ids_list[other_idxs], function(k) { 
      all(k %in% curr) & all(curr %in% k) 
    }) %>% 
    unlist() %>% 
    all()
  res[[i]] <- in_all_others
}

all_ids_are_matching <- 
  res %>% 
  unlist() %>% 
  all()

# validate_transects(transects, "id")
if (!all_ids_are_matching) {
  warning("Not all id/cs_id are included in all transects / cross section point datasets")
} else {
 
  message("All id / cs_id are correctly within all transects / cross section points datasets") 
}

trans_path   <- paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/transects.gpkg")
transects    <- sf::read_sf(trans_path) %>% 
                dplyr::rename(geometry = geom)

is_validated_transects           <- validate_transects(transects, "id")

is_flowline_validated_transects  <- validate_transects_against_flowlines(transects, flines, "id")

if (!(is_validated_transects && is_flowline_validated_transects)) {
  stop(trans_path, "  transects failed validity check:",
       "\n > is_validated_transects: ", is_validated_transects,
       "\n > is_flowline_validated_transects: ", is_flowline_validated_transects
  )
} else {
  message("Transects look good to go!")
}

output_files <- list.files(DOMAIN_WITH_FEMA_OUTPUT_DIR, full.names = T)
cs_pts_files <- output_files[grepl(".parquet", output_files)]

for (i in seq_along(cs_pts_files)) {
  path <- cs_pts_files[i]
  message(i, " - Checking cross section points validity for > '",   basename(path), "'")
  
  cs_pts <-arrow::read_parquet(path)
  
  is_validated_cs_pts              <- validate_cs_pts(cs_pts, crosswalk_id = "id")
  is_transect_validated_cs_pts     <- validate_cs_pts_against_transects(cs_pts, transects, crosswalk_id = "id")
  
  is_valid <- is_validated_cs_pts & is_transect_validated_cs_pts
  
  if(!is_valid) {
    stop(path, " cs_pts failed validity check:",
         "\n > is_validated_cs_pts: ", is_validated_cs_pts,
         "\n > is_transect_validated_cs_pts: ", is_transect_validated_cs_pts
         )
  } else {
    message("'",   basename(path), "' is valid and is good to go!")
  }
  
  # message("'",   basename(path), "' is validated")

  }

# og_trans <- sf::read_sf(paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/pre-extension-transects.gpkg"))
# post_trans <- sf::read_sf(paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/post-extension-transects.gpkg"))
# dplyr::filter(transects, id == "wb-10813") %>% 
#   sf::st_buffer(2500) 
# 
# bb <- 
#   dplyr::filter(transects, id == "wb-10813") %>% 
#   sf::st_buffer(5000) %>% 
#   sf::st_bbox() %>% 
#   sf::st_as_sfc() %>% 
#   sf::st_sf()
# 
#   sf::st_filter(transects, bb) 
#   sf::st_filter(transects, bb) 
#   sf::st_filter(transects, bb) 
# 
#   mapview::mapview(dplyr::filter(flines, id == "wb-10813")) +
#   mapview::mapview(    sf::st_filter(transects, bb) , color  = "green") + 
#   mapview::mapview(  sf::st_filter(og_trans, bb) , color = "red") + 
#     mapview::mapview(  sf::st_filter(post_trans, bb) , color = "gold")
# 
# transects %>% 
#   sf::st_filter(
#     bb
#   ) %>% mapview::mapview()
# 
# mapview::mapview() + mapview::mapview(dplyr::filter(flines, id == "wb-10813")) +
#   mapview::mapview(  dplyr::filter(transects, id == "wb-10813"), color  = "green") + 
#   mapview::mapview(dplyr::filter(og_trans, id == "wb-10813"), color = "red")
# 
# 
# dplyr::filter(og_trans, id == "wb-10813")
# dplyr::filter(post_trans, id == "wb-10813")
# dplyr::filter(transects, id == "wb-10813")
# 
# mapview::mapview(dplyr::filter(flines, id == "wb-10813")) +
# mapview::mapview(  dplyr::filter(transects, id == "wb-10813"), color  = "green") + 
# mapview::mapview(dplyr::filter(og_trans, id == "wb-10813"), color = "red")
# 
# transects %>% 
#   dplyr::group_by(id)  %>% 
#   dplyr::filter(cs_id == max(cs_id)) %>% 
#   dplyr::arrange(cs_id)
# 
# # cs_pts_files <- 
# cs_pts <-arrow::read_parquet(paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/dem_cross-sections.parquet"))
# 
# is_validated_transects           <- validate_transects(transects, "id")
# is_flowline_validated_transects  <- validate_transects_against_flowlines(transects, flines, "id")
# 
# is_validated_cs_pts              <- validate_cs_pts(cs_pts, crosswalk_id = "id")
# is_transect_validated_cs_pts     <- validate_cs_pts_against_transects(cs_pts, transects, crosswalk_id = "id")

# sf::st_layers(paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/transects.gpkg"))_

validate_transects_self_intersections <- function(transects) {
  return(
    nrow(transects) == nrow(hydrofabric3D:::rm_self_intersections(transects))
  )
}

validate_transects_cs_id_enumeration <- function(transects, crosswalk_id = NULL) {
  
  # reenumerate the cs_ids for each transect based on cs_measure sorting, and make sure all cross sections are correctly numbered
  mismatches <-
    transects %>%
    sf::st_drop_geometry() %>% 
    dplyr::group_by(dplyr::across(dplyr::any_of(c(crosswalk_id)))) %>% 
    dplyr::arrange(cs_measure, .by_group = TRUE) %>% 
    dplyr::mutate(
      new_cs_id = 1:dplyr::n()
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(cs_id != new_cs_id)
  
  # FALSE if there are any transects with different cs_ids to the newly created cs_id 
  # Otherwise TRUE
  has_valid_cs_ids <- !(nrow(mismatches) > 0)
  
  return(
    has_valid_cs_ids
  )
  
}

validate_transects_cs_length <- function(transects, crosswalk_id = NULL) {
  
  # re calculate transect geometry length and compare to cs_lengthm column
  wrong_lengths <-
    transects %>%
    dplyr::mutate(
      new_cs_length = as.numeric(sf::st_length(.)) 
    ) %>% 
    dplyr::filter(cs_lengthm != new_cs_length)
  
  # FALSE if there are any transects with different cs_lengthm than the freshly calculated new_cs_length
  has_correct_lengths <- !(nrow(wrong_lengths) > 0)
  
  return(
    has_correct_lengths
  )
  
}

validate_transects_unique_ids <- function(transects, crosswalk_id = NULL) {
    
  duplicate_ids <- 
    transects %>%
    sf::st_drop_geometry() %>% 
    hydrofabric3D::add_tmp_id(x = crosswalk_id) %>% 
    dplyr::select(tmp_id) %>% 
    dplyr::group_by(tmp_id) %>% 
    dplyr::count() %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(n > 1)
    
  # FALSE if there are ANY rows in the duplicate_ids dataframe above 
  # (i.e. a count of greater than 1 for any tmp_id (<crosswalk_id>_cs_id))
  has_unique_ids <- !(nrow(duplicate_ids) > 0)
  
  return(
    has_unique_ids
  )
  
}

validate_transects_cs_measure <- function(transects) {
  
  min_cs_measure <-
    transects %>%
    sf::st_drop_geometry() %>% 
    dplyr::pull(cs_measure) %>% 
    min()
  
  max_cs_measure <-
    transects %>%
    sf::st_drop_geometry() %>% 
    dplyr::pull(cs_measure) %>% 
    max()
  
  # cs_measure should always be:
  # greater than or equal to 0 AND 
  # less than or equal to 100 as its a percentage along a flowline
  has_valid_cs_measure <- min_cs_measure >= 0 & max_cs_measure <= 100
  
  return(
    has_valid_cs_measure
  )
  
}

validate_transects_has_complete_geometries <- function(transects, crosswalk_id = NULL) {
  
  has_empty_geoms <- 
    transects %>% 
    sf::st_is_empty() %>% 
    any()
  
  return(
    !has_empty_geoms
  )
  
}

validate_transects_has_crs <- function(transects) {
  
  missing_crs <- 
    transects %>% 
    sf::st_crs() %>% 
    is.na()
    
  return(
    !missing_crs
  )
  
}

validate_transects <- function(transects, 
                               crosswalk_id = NULL
                               ) {
  
  # # standardize geometry name
  # transects <- hydroloom::rename_geometry(transects, "geometry")
  
  REQUIRED_COLS <- c(crosswalk_id, "cs_id", "cs_source", "cs_measure", "cs_lengthm", "geometry")
  
  # validate dataframe has all correct columns  
  has_all_valid_cols         <- hydrofabric3D:::validate_df(
                                          x = transects, 
                                          cols = REQUIRED_COLS,
                                          obj_name = "transects"
                                          )  
  
  # Validate every flowline (id) has a cs_id of 1:number of transects
  has_valid_cs_ids           <- validate_transects_cs_id_enumeration(transects, crosswalk_id = crosswalk_id)
  
  # Validate there are no self intersections
  has_no_self_intersections  <- validate_transects_self_intersections(transects)
  
  # validate the cs_lengthm column equals the actual transect geometry length
  has_correct_lengths        <- validate_transects_cs_length(transects, crosswalk_id = crosswalk_id)
  
  # validate no duplicate id / cs_id combos
  has_unique_ids             <- validate_transects_unique_ids(transects, crosswalk_id = crosswalk_id)
  
  # validate cs measure is never greater than 100 (i think)
  has_valid_cs_measure       <- validate_transects_cs_measure(transects)
  
  # make sure transects have no empty geometries 
  has_complete_geoemetries   <- validate_transects_has_complete_geometries(transects)
  
  # make sure transects have a CRS
  has_crs                    <- validate_transects_has_crs(transects)
  
  # if everything is TRUE, return true, otherwise return FALSE (or throw an error...?)
  is_validated_transects <- all(
                              c(
                                has_all_valid_cols,
                                has_valid_cs_ids,
                                has_no_self_intersections,
                                has_correct_lengths,
                                has_unique_ids,
                                has_valid_cs_measure,
                                has_complete_geoemetries,
                                has_crs
                              )
                            )
  
  return(is_validated_transects)
  
}

# all ids in transects are in flowlines
validate_transects_ids_in_flowlines <- function(transects, flowlines, crosswalk_id = NULL) {
  
  # flowlines <- flines
  
  transect_ids <- 
    transects %>% 
    sf::st_drop_geometry() %>% 
    dplyr::select(dplyr::any_of(crosswalk_id)) %>% 
    dplyr::pull(dplyr::any_of(crosswalk_id)) %>% 
    unique()
  
  flowline_ids <-
    flowlines %>% 
    sf::st_drop_geometry() %>% 
    dplyr::select(dplyr::any_of(crosswalk_id)) %>% 
    dplyr::pull(dplyr::any_of(crosswalk_id)) %>% 
    unique()
  
  all_transect_ids_in_flowline_ids <- all(transect_ids %in% flowline_ids)
  
  return(
    all_transect_ids_in_flowline_ids
  )
  
}

# make sure no transect crosses more than a single flowline, a single time
validate_transects_flowline_intersections <- function(transects, flowlines) {
  
  # flowlines <- flines
  return (
    nrow(transects) == nrow(hydrofabric3D:::rm_multiflowline_intersections(transects, flowlines)) 
    )
}

# validate 2 SF objects have the same CRS
validate_same_crs <- function(x, y) {
  
  return (
    sf::st_crs(x) == sf::st_crs(y)
  )
  
}

validate_transects_against_flowlines <- function(transects, 
                                                 flowlines,  
                                                 crosswalk_id = NULL
                                                ) {
  
  # # standardize geometry name
  # transects <- hydroloom::rename_geometry(transects, "geometry")
  
  REQUIRED_COLS <- c(crosswalk_id, "cs_id", "cs_source", "cs_measure", "cs_lengthm", "geometry")
  
  # validate dataframe has all correct columns  
  has_all_valid_cols         <- hydrofabric3D:::validate_df(
                                                x = transects, 
                                                cols = REQUIRED_COLS,
                                                obj_name = "transects"
                                              )  
  
  # all ids in transects are in flowlines
  all_transect_ids_in_flowline_ids  <- validate_transects_ids_in_flowlines(transects, flowlines, crosswalk_id = crosswalk_id)
  
  # transects only intersects a single flowline, a single time
  has_valid_flowline_intersects     <- validate_transects_flowline_intersections(transects, flowlines)
  
  # transects and flowlines have the same CRS
  has_same_crs                      <- validate_same_crs(transects, flowlines)
  
  
  # if everything is TRUE, return true, otherwise return FALSE (or throw an error...?)
  is_flowline_validated_transects <- all(
                                        c(
                                          has_all_valid_cols,
                                          all_transect_ids_in_flowline_ids,
                                          has_valid_flowline_intersects,
                                          has_same_crs
                                        )
                                      )
  
  return(is_flowline_validated_transects)
  
}

# validate_transects_against_cs_pts <- function(
#                                           transects, 
#                                           cs_pts,  
#                                           crosswalk_id = NULL
#                                           ) {
#   
#   # # standardize geometry name
#   # transects <- hydroloom::rename_geometry(transects, "geometry")
#   
#   REQUIRED_COLS <- c(crosswalk_id, "cs_id", "cs_source", "cs_measure", "cs_lengthm", "geometry")
#   
#   # validate dataframe has all correct columns  
#   has_all_valid_cols         <- hydrofabric3D:::validate_df(
#     x = transects, 
#     cols = REQUIRED_COLS,
#     obj_name = "transects"
#   )  
# }

validate_cs_pts_cs_id_enumeration <- function(cs_pts, crosswalk_id = NULL) {
  
  # reenumerate the cs_ids for each transect based on cs_measure sorting, and make sure all cross sections are correctly numbered
  mismatches <-
    cs_pts %>%
    # dplyr::slice(1:150) %>% 
    sf::st_drop_geometry() %>% 
    dplyr::select(dplyr::any_of(crosswalk_id), cs_id) %>% 
    dplyr::group_by(dplyr::across(dplyr::any_of(c(crosswalk_id, "cs_id")))) %>% 
    dplyr::slice(1) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(dplyr::across(dplyr::any_of(c(crosswalk_id)))) %>% 
    dplyr::arrange(cs_id, .by_group = TRUE) %>% 
    dplyr::mutate(
      new_cs_id = 1:dplyr::n()
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(cs_id != new_cs_id)
  
  # FALSE if there are any transects with different cs_ids to the newly created cs_id 
  # Otherwise TRUE
  has_valid_cs_ids <- !(nrow(mismatches) > 0)
  
  return(
    has_valid_cs_ids
  )
  
}

validate_cs_pts_pt_id_enumeration <- function(cs_pts, crosswalk_id = NULL) {
  
  # reenumerate the pt_ids to make sure the pt_ids are valid values of 1:number of points in cross section
  mismatches <-
    cs_pts %>%
    sf::st_drop_geometry() %>% 
    dplyr::select(dplyr::any_of(crosswalk_id), cs_id, pt_id) %>% 
    dplyr::group_by(dplyr::across(dplyr::any_of(c(crosswalk_id, "cs_id")))) %>% 
    dplyr::arrange(pt_id, .by_group = TRUE) %>% 
    dplyr::mutate(
      new_pt_id = 1:dplyr::n()
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(pt_id != new_pt_id)
  
  # FALSE if there are any cs_pts with  pt_ids different from the newly created new_pt_id 
  # Otherwise TRUE
  has_valid_pt_ids <- !(nrow(mismatches) > 0)
  
  return(
    has_valid_pt_ids
  )
  
}

validate_cs_pts_relative_distance <- function(cs_pts, crosswalk_id = NULL) {
  
  # make sure relative distance is greater than or equal to 0
  min_relative_distance <-
    cs_pts %>%
    # dplyr::slice(1:50) %>% 
    sf::st_drop_geometry() %>% 
    dplyr::pull(relative_distance) %>% 
    min()
  
  
  # reenumerate the pt_ids to make sure the pt_ids are valid values of 1:number of points in cross section
  rel_dist_check <-
    cs_pts %>%
    # dplyr::slice(1:50) %>% 
    sf::st_drop_geometry() %>% 
    dplyr::select(dplyr::any_of(crosswalk_id), cs_id, pt_id, relative_distance, cs_lengthm) %>% 
    dplyr::group_by(dplyr::across(dplyr::any_of(c(crosswalk_id, "cs_id")))) %>% 
    dplyr::summarise(
      cs_lengthm = max(cs_lengthm),
      # min_rel_dist = min(relative_distance),
      max_rel_dist = max(relative_distance)
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
        # TODO: as long as the lengths are within 1 meter, thats equal
        is_valid_relative_dist = abs(cs_lengthm - max_rel_dist) <= 1
        # approx_equal_lengths = all.equal(cs_lengthm, cs_pts_lengthm, tolerance = 0.01)
      ) 
    
    # dplyr::filter(max_rel_dist > cs_lengthm)
  
  # relative distance is always greater than or equal to 0 and less than the cross sections length
  has_valid_relative_dist_min        <- min_relative_distance >= 0
  has_valid_relative_dist_maximums   <- all(rel_dist_check$is_valid_relative_dist)
  
  return(
    has_valid_relative_dist_min && has_valid_relative_dist_maximums
  )
  
}

validate_cs_pts_point_types <- function(cs_pts) {
  
  # make sure only "left_bank", "right_bank", "channel", and "bottom" values exist in cs_pts point_type column
  valid_point_types  <- c("left_bank", "right_bank", "channel", "bottom")
  
  # unique point types in cs_pts
  unique_point_types <- unique(cs_pts$point_type)
  
  has_only_valid_point_types <- all(unique_point_types %in% valid_point_types)
  
  return(
    has_only_valid_point_types
  )
  
}

validate_cs_pts <- function(
    cs_pts,  
    crosswalk_id = NULL
) {
  
  # # standardize geometry name
  # transects <- hydroloom::rename_geometry(transects, "geometry")
  
  REQUIRED_COLS <- c(crosswalk_id, "cs_id", "pt_id", 
                     "relative_distance", "cs_lengthm", "X", "Y", "Z", "Z_source",
                     "class", "point_type", "valid_banks", "has_relief"
  )
  
  # validate dataframe has all correct columns  
  has_all_valid_cols         <- hydrofabric3D:::validate_df(
    x = cs_pts, 
    cols = REQUIRED_COLS,
    obj_name = "cs_pts"
  )  
  
  # make sure valid cs_ids
  has_valid_cs_pts_cs_ids       <- validate_cs_pts_cs_id_enumeration(cs_pts, crosswalk_id = crosswalk_id)
  
  # make sure valid pt_ids
  has_valid_cs_pts_pt_ids       <- validate_cs_pts_pt_id_enumeration(cs_pts, crosswalk_id = crosswalk_id)
  
  # check cs_pts have only valid relative_distance values
  has_valid_relative_distances  <- validate_cs_pts_relative_distance(cs_pts, crosswalk_id = crosswalk_id)
  
  has_valid_point_types         <- validate_cs_pts_point_types(cs_pts)
  
  # if everything is TRUE, return true, otherwise return FALSE (or throw an error...?)
  is_validated_cs_pts <- all(
    c(
      has_all_valid_cols,
      has_valid_cs_pts_cs_ids,
      has_valid_cs_pts_pt_ids,
      has_valid_relative_distances,
      has_valid_point_types
    )
  )
  
  return(is_validated_cs_pts)
  
}

# validate all cs_pts id/cs_ids are in the transects
validate_cs_pt_ids_in_transects <- function(cs_pts, transects, crosswalk_id = NULL) {
  
  # flowlines <- flines
  cs_pts_ids <-
    cs_pts %>% 
    sf::st_drop_geometry() %>% 
    hydrofabric3D::get_unique_tmp_ids(x = crosswalk_id)
  
  transect_ids <-
    transects %>% 
    sf::st_drop_geometry() %>% 
    hydrofabric3D::get_unique_tmp_ids(x = crosswalk_id)
  
  
  all_cs_pts_ids_in_transects <- all(cs_pts_ids %in% transect_ids)
  all_transect_ids_in_cs_pts  <- all(transect_ids %in% cs_pts_ids)
  same_number_of_ids          <- length(cs_pts_ids) == length(transect_ids)  
  
  is_valid_cs_pts_ids <- all(
    c(
      all_cs_pts_ids_in_transects,
      all_transect_ids_in_cs_pts,
      same_number_of_ids
    )
  )
  
  return(
    is_valid_cs_pts_ids
  )
  
}

validate_cs_pts_length_against_transects <- function(cs_pts, transects, crosswalk_id = NULL) {
  
  cs_pt_lengths <-
    cs_pts %>% 
    sf::st_drop_geometry() %>% 
    dplyr::select(dplyr::any_of(crosswalk_id), cs_id, cs_lengthm) %>% 
    dplyr::group_by(dplyr::across(dplyr::any_of(c(crosswalk_id, "cs_id")))) %>% 
    dplyr::slice(1) %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(
      cs_pts_lengthm = cs_lengthm
    )
  
  transect_lengths <-
    transects %>% 
    sf::st_drop_geometry() %>% 
    dplyr::select(dplyr::any_of(crosswalk_id), cs_id, cs_lengthm) %>% 
    dplyr::group_by(dplyr::across(dplyr::any_of(c(crosswalk_id, "cs_id")))) %>% 
    dplyr::slice(1) %>% 
    dplyr::ungroup()
  
  lengths_check <- 
    dplyr::left_join(
      transect_lengths,
      cs_pt_lengths,
      by = c(crosswalk_id, "cs_id")
    ) %>% 
    dplyr::mutate(
      # TODO: as long as the lengths are within 1 meter, thats equal
      approx_equal_lengths = abs(cs_lengthm - cs_pts_lengthm) <= 1
      # approx_equal_lengths = all.equal(cs_lengthm, cs_pts_lengthm, tolerance = 0.01)
    ) 
  
  all_lengths_are_equal <- all(lengths_check$approx_equal_lengths)
  
  return(
    all_lengths_are_equal
  )
  
}

validate_cs_pts_against_transects <- function(
    cs_pts,  
    transects, 
    crosswalk_id = NULL
) {
  
  # # standardize geometry name
  # transects <- hydroloom::rename_geometry(transects, "geometry")
  
  # make sure all id/cs_id combos are in both transects and cs_pts
  has_valid_cs_pts_ids  <- validate_cs_pt_ids_in_transects(cs_pts, transects, crosswalk_id = crosswalk_id)
  
  # make sure cs_lengthm matches from transects to cs_pts
  has_matching_lengths  <- validate_cs_pts_length_against_transects(cs_pts, transects, crosswalk_id = crosswalk_id)
  
  # if everything is TRUE, return true, otherwise return FALSE (or throw an error...?)
  is_transect_validated_cs_pts <- all(
    c(
      has_valid_cs_pts_ids,
      has_matching_lengths
    )
  )
  
  return(is_transect_validated_cs_pts)
  
}

# og_transects   <- sf::read_sf(paste0(DOMAIN_WITH_FEMA_TRANSECTS_DIR, "/cs-point-extension-transects.gpkg"))
# new_transects  <- sf::read_sf(paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/transects.gpkg"))

# new_transects %>% 
#   dplyr::arrange(cs_lengthm) %>% 
#   dplyr::slice(1:1000) %>% 
#   mapview::mapview()
# 
# hydrofabric3D:::rm_self_intersections(new_transects)

# --------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------
# 
og_trans <- sf::read_sf( "/Users/anguswatters/Desktop/lynker-spatial/domain_with_fema/transects/pre_extension_transects.gpkg")

transects <-  sf::read_sf(output_files[[4]])

transects2 <- hydrofabric3D:::rm_self_intersections(transects)

old_ids <- transects %>% hydrofabric3D::get_unique_tmp_ids("id")
new_ids <- transects2 %>% hydrofabric3D::get_unique_tmp_ids("id")
diff_tmp_ids <- old_ids[!old_ids %in% new_ids]
strsplit(diff_ids, "_")

diff_ids <- lapply(diff_tmp_ids, function(i) {
  strsplit(i, "_")[[1]][1]
}) %>% unlist()

diff_trans <- 
  transects %>% 
  # hydrofabric3D::add_tmp_id(x = "id") %>% 
  dplyr::filter(id %in% diff_ids)
  # sf::st_buffer(50) %>% 
  # sf::st_bbox() %>% 
  # sf::st_as_sfc() %>% 
  # sf::st_sf()

og_diff_trans <- 
  og_trans %>% 
  # hydrofabric3D::add_tmp_id(x = "id") %>% 
  dplyr::filter(id %in% diff_ids)

mapview::mapview(diff_trans, color = "green") +
  mapview::mapview(og_diff_trans, color = "red")

old_trans <- transects[lengths(sf::st_intersects(transects, diff_bb)) > 0, ]

# final_cross_sections %>% 
#   dplyr::filter(id == "wb-1000") %>% 
#   dplyr::rename(hy_id = id) %>% 
#   hydrofabric3D::plot_cs_pts(color = "point_type")

# ---------------------------------------------------------------------------------
# ---- Write final cross section points data ----
# ---- Diffusive Domain DEM + FEMA + ML 
# ---------------------------------------------------------------------------------

CROSS_SECTIONS_ML_OUTPUT_PATH <- paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/cross-sections.parquet")

message(round(Sys.time()), " - Saving Diffusive DEM + FEMA + ML augmented cross section points to:\n - filepath: '", CROSS_SECTIONS_ML_OUTPUT_PATH, "'")
# transects <- 
# sum(is.na(final_cross_sections$id))
# sum(is.na(final_cross_sections$cs_id))
# sum(is.na(final_cross_sections$pt_id))
# sum(is.na(final_cross_sections$X))

# save cross section points as a parquet to out_path (domain/outputs/cross-sections.parquet)
arrow::write_parquet(
  # dplyr::select(final_cs, 
  #               -is_dem_point
  # ), 
  cross_sections,
  CROSS_SECTIONS_ML_OUTPUT_PATH
)


























