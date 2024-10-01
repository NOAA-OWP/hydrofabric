# Generate the transects + cs_pts + cross sections layers for a single flowlines domain file and DEM file 
source("runners/cs_runner/config.R")
source("runners/cs_runner/utils.R")

# load libraries
# library(hydrofabric3D)
# library(dplyr)
# library(sf)

# Unique Flowline ID column name
CROSSWALK_ID <- "id"

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
  # dplyr::slice(1:1000)

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
  # dplyr::slice(1:1000)

# # Replace any topwidth values that are GREATER than the actual cross section length (meters)
# bankful_cs2 <- hydrofabric3D:::fix_oversized_topwidths(
#   cross_section_pts = bankful_cs2,
#   crosswalk_id = CROSSWALK_ID
# )
# bankful_cs %>% 
#   dplyr::filter(TW >= cs_lengthm) %>% 
#   dplyr::select(TW)

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
# system.time({
  
  # Add bathymetry using "inchannel" estimates
  inchannel_cs <- hydrofabric3D::add_cs_bathymetry(
    cross_section_pts = inchannel_cs,
    # cross_section_pts = tmp,
    crosswalk_id      = CROSSWALK_ID
  )
  
# })

  gc()


# arrow::write_parquet(inchannel_cs, "/Users/anguswatters/Desktop/test_ml_cs_pts_06.parquet")
# ml_subset %>%
#   dplyr::filter(hy_id == "wb-1005207") %>%
#   dplyr::select(owp_y_inchan, owp_tw_inchan) %>% 
#   .$owp_y_inchan
message(round(Sys.time()), " - Adding cross section bathymetry using bankful widths/depths estimates...")
# system.time({
  
  # Add bathymetry using "bankful" estimates
  bankful_cs <- hydrofabric3D::add_cs_bathymetry(
    cross_section_pts = bankful_cs,
    # cross_section_pts = dplyr::slice(bankful_cs, 1:10000),
    crosswalk_id      = CROSSWALK_ID
  )
  
# })

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
tmp_path <-   paste0(DOMAIN_WITH_FEMA_OUTPUT_DIR, "/final_cs.parquet")
message("saving file", tmp_path)

arrow::write_parquet(
  final_cs,
  tmp_path
)

