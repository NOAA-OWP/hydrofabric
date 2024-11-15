# Create an empty file structure 
# base_dir: character, top level directory path
# Directory tree:
# base_dir/
#   └── lynker-spatial/
#     ├── hydrofabric/
#     ├── dem/
#         ├── vrt/
#         ├── tif/
#     ├── cs-extension-polygons/
create_local_hydrofabric_base_dirs <- function(base_dir) {
  
  
  # build paths
  hydrofabric_dir <- paste0(base_dir, "/hydrofabric")
  
  # DEM dirs
  dem_dir                    <- file.path(base_dir, "dem")
  dem_vrt_dir                <- file.path(dem_dir, "vrt")
  dem_tif_dir                <- file.path(dem_dir, "tif")
  
  # polygons for transect extensions
  cs_extension_polygons_dir     <- file.path(base_dir, "cs-extension-polygons")
  
  # FEMA data 
  fema_dir                      <- file.path(base_dir, "fema")
  
  fema_fgb_dir     <- file.path(fema_dir, "fema-fgb")
  fema_geojson_dir <- file.path(fema_dir, "fema-geojson")
  fema_clean_dir   <- file.path(fema_dir, "fema-clean")
  fema_gpkg_dir    <- file.path(fema_dir, "fema-gpkg")
  
  # BY VPU folders 
  VPU_IDS               <- get_vpu_ids()
  fema_by_vpu_dir       <- file.path(fema_dir, "fema-by-vpu")
  fema_by_vpu_subdirs   <- paste0(fema_by_vpu_dir, "/vpu-", VPU_IDS)
  
  # create base directories
  create_if_not_exists(base_dir)
  create_if_not_exists(hydrofabric_dir)
  
  # DEM dirs
  create_if_not_exists(dem_dir)
  create_if_not_exists(dem_vrt_dir)
  create_if_not_exists(dem_tif_dir)
  
  # extension polygons
  create_if_not_exists(cs_extension_polygons_dir)
  
  
  # Create FEMA folders
  create_if_not_exists(fema_dir)
  create_if_not_exists(fema_fgb_dir)
  create_if_not_exists(fema_geojson_dir)
  create_if_not_exists(fema_clean_dir)
  create_if_not_exists(fema_gpkg_dir)
  create_if_not_exists(fema_by_vpu_dir)
  
  for (path in fema_by_vpu_subdirs) {
    create_if_not_exists(path)
  }
  
}

get_vpu_ids <- function() {
  VPU_IDS              <- c('01', '02', '03N', '03S', '03W', '04', '05', '06', '07', '08', '09', 
                            '10L', '10U', '11', '12', '13', '14', '15', '16', '17', '18', '20', '21')
  # VPU_IDS              <- sf::st_drop_geometry(nhdplusTools::get_boundaries())$VPUID
  
  return(VPU_IDS)
  
}

# # Base directory for local file storage
# BASE_DIR    <- '/Volumes/T7SSD/lynker-spatial'
# base_dir <- BASE_DIR
# # FEMA100 year flood map FGB save location (temporary, will be deleted after processing)
# FEMA_FGB_PATH        <- file.path(BASE_DIRS_LIST$fema_dir, "fema_fgb")
# FEMA_GEOJSON_PATH    <- file.path(BASE_DIRS_LIST$fema_dir, "fema_geojson")
# FEMA_CLEAN_PATH      <- file.path(BASE_DIRS_LIST$fema_dir, "fema_clean")
# FEMA_GPKG_PATH       <- file.path(BASE_DIRS_LIST$fema_dir, "fema_gpkg")
# FEMA_BY_VPU_PATH     <- file.path(BASE_DIRS_LIST$fema_dir, "FEMA_BY_VPU")
# 
# 
# VPU_IDS              <- c('01', '02', '03N', '03S', '03W', '04', '05', '06', '07', '08', '09', 
#                           '10L', '10U', '11', '12', '13', '14', '15', '16', '17', '18', '20', '21')
# # VPU_IDS              <- sf::st_drop_geometry(nhdplusTools::get_boundaries())$VPUID
# 
# paste0("'", sf::st_drop_geometry(nhdplusTools::get_boundaries())$VPUID, "'", collapse = ", ")
# FEMA_VPU_SUBFOLDERS  <- paste0(FEMA_BY_VPU_PATH, "/VPU_", VPU_IDS)

# Create an empty file structure for a new version within a specified base_dir
# base_dir: character, top level directory path
# Directory tree:
# base_dir/
#   └── lynker-spatial/
#     ├── hydrofabric/
    #     ├── version_number/
    #         ├── network/
    #         ├── transects/
    #         ├── cross-sections/
          #         ├── dem/
          #         ├── dem-ml/
          #         ├── dem-coastal-bathy/
          #         ├── dem-points/
create_new_version_dirs <- function(base_dir, version, with_output = FALSE) {
  # version = "v3.0"
  # base_dir <- BASE_DIR
  
  # build paths
  hydrofabric_dir   <- paste0(base_dir, "/hydrofabric")
  version_base_dir  <- paste0(hydrofabric_dir, "/", version)
  
  # polygons for transect extensions
  ml_dir                     <- paste0(version_base_dir, "/ml")
  
  # reference features 
  ref_features_dir           <- paste0(version_base_dir, "/reference-features")
  
  # conus network gpkg
  network_dir                <- paste0(version_base_dir, "/network")
  
  # transects
  transects_dir              <- paste0(version_base_dir, "/transects")
  
  # cross sections dirs
  cross_sections_dir                 <- paste0(version_base_dir, "/cross-sections")
  cross_sections_dem_dir             <- paste0(cross_sections_dir, "/dem")
  cross_sections_ml_dir              <- paste0(cross_sections_dir, "/dem-ml")
  cross_sections_coastal_bathy_dir   <- paste0(cross_sections_dir, "/dem-coastal-bathy")
  cross_sections_dem_pts_dir         <- paste0(cross_sections_dir, "/dem-points")
  
  if(with_output) {
    output_dir       <- paste0(version_base_dir, "/outputs")
  }
  
  # create version BASE dir
  create_if_not_exists(version_base_dir)
  
  # CONUS dir
  create_if_not_exists(network_dir)
  
  # ML data
  create_if_not_exists(ml_dir)
  
  # reference features data
  create_if_not_exists(ref_features_dir)
  
  # transects
  create_if_not_exists(transects_dir)
  
  # CS pts
  create_if_not_exists(cross_sections_dir)
  create_if_not_exists(cross_sections_dem_dir)
  create_if_not_exists(cross_sections_ml_dir)
  create_if_not_exists(cross_sections_coastal_bathy_dir)
  create_if_not_exists(cross_sections_dem_pts_dir)
  
  if(with_output) {
    create_if_not_exists(output_dir)
  }
  
}

create_if_not_exists <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Created directory: '", dir_path, "'\n")
  }
}

# get a list of top level directories for main directory
get_base_dir_paths <- function(base_dir) {
  # base_dir = BASE_DIR
  # version = "v3.0"
  
  hydrofabric_dir  <- file.path(base_dir, "hydrofabric")
  
  dem_dir          <- file.path(base_dir, "dem")
  dem_vrt_dir      <- file.path(base_dir, "dem", "vrt")
  dem_tif_dir      <- file.path(base_dir, "dem", "tif")
  
  cs_extension_polygons_dir <- file.path(base_dir, "cs-extension-polygons")
  
  # FEMA data 
  
  fema_dir         <- file.path(base_dir, "fema")
  
  fema_fgb_dir     <- file.path(fema_dir, "fema-fgb")
  fema_geojson_dir <- file.path(fema_dir, "fema-geojson")
  fema_clean_dir   <- file.path(fema_dir, "fema-clean")
  fema_gpkg_dir    <- file.path(fema_dir, "fema-gpkg")
  
  # BY VPU folders 
  VPU_IDS               <- get_vpu_ids()
  fema_by_vpu_dir       <- file.path(fema_dir, "fema-by-vpu")
  fema_by_vpu_subdirs   <- paste0(fema_by_vpu_dir, "/vpu-", VPU_IDS)
  
  return(
      list(
      hydrofabric_dir = hydrofabric_dir,
      dem_dir = dem_dir,
      dem_vrt_dir = dem_vrt_dir, 
      dem_tif_dir = dem_tif_dir,
      cs_extension_polygons_dir = cs_extension_polygons_dir,
      fema_dir = fema_dir,
      fema_fgb_dir = fema_fgb_dir,
      fema_geojson_dir = fema_geojson_dir,
      fema_clean_dir = fema_clean_dir,
      fema_gpkg_dir = fema_gpkg_dir,
      fema_by_vpu_dir = fema_by_vpu_dir,
      fema_by_vpu_subdirs = fema_by_vpu_subdirs
    )
  )
}

# get list of a specific directories in a version directory
get_version_base_dir_paths <- function(base_dir, version) {
  # base_dir = BASE_DIR
  # version = "v3.0"
  
  hydrofabric_dir   <- file.path(base_dir, "hydrofabric")
  
  version_base_dir  <- file.path(hydrofabric_dir, version)
  
  # polygons for transect extensions
  ml_dir                     <- file.path(version_base_dir, "ml")
  
  # reference features 
  ref_features_dir           <- file.path(version_base_dir, "reference-features")
  
  # conus network gpkg
  network_dir                <- file.path(version_base_dir, "network")
  
  # transects
  transects_dir              <- file.path(version_base_dir, "transects")
  
  # cross sections dirs
  cross_sections_dir                 <- file.path(version_base_dir, "cross-sections")
  cross_sections_dem_dir             <- file.path(cross_sections_dir, "dem")
  cross_sections_ml_dir              <- file.path(cross_sections_dir, "dem-ml")
  cross_sections_coastal_bathy_dir   <- file.path(cross_sections_dir, "dem-coastal-bathy")
  cross_sections_dem_pts_dir         <- file.path(cross_sections_dir, "dem-points")
  
  return(
    list(
      hydrofabric_dir    = hydrofabric_dir,
      version_base_dir   = version_base_dir,
      ref_features_dir   = ref_features_dir, 
      network_dir        = network_dir,
      ml_dir = ml_dir,
      transects_dir      = transects_dir,
      cross_sections_dir = cross_sections_dir,
      cross_sections_dem_dir     = cross_sections_dem_dir,
      cross_sections_dem_pts_dir = cross_sections_dem_pts_dir,
      cross_sections_ml_dir      = cross_sections_ml_dir,
      cross_sections_coastal_bathy_dir = cross_sections_coastal_bathy_dir
    )
  )
}

list_s3_objects <- function(s3_bucket, pattern = NULL, aws_profile = NULL) {
  
  profile_option <- if (!is.null(aws_profile)) paste0("--profile ", aws_profile) else ""
  
  if (is.null(pattern) || pattern == "") {
    grep_command <- ""  # no filtering if empty or NULL
  } else {
    grep_command <- paste0(" | grep -E \"", pattern, "\"")  # grep if a pattern is given
  }
  
  cmd <- paste0(
    '#!/bin/bash\n',
    'S3_BUCKET="', s3_bucket, '"\n',
    'PATTERN="', pattern, '"\n',
    'S3_OBJECTS=$(aws s3 ls "$S3_BUCKET" ', profile_option, ' | awk \'{print $4}\' | grep -E "$PATTERN")\n',
    'echo "$S3_OBJECTS"'
  )
  # cmd <- paste0(
  #   '#!/bin/bash\n',
  #   'S3_BUCKET="', s3_bucket, '"\n',
  #   'S3_OBJECTS=$(aws s3 ls "$S3_BUCKET" ', profile_option, ' | awk \'{print $4}\'', grep_command, ')\n',
  #   'echo "$S3_OBJECTS"'
  # )
  ls_output <- system(cmd, intern = TRUE)
  return(ls_output)
}

download_tiles <- function(tile_paths, output_dir) {
  # output_dir <- DEM_TIF_DIR
  # tile_paths
  # error_tiles <- 
  tif_save_paths       <- paste0(output_dir, "/", basename(tile_paths))
  
  error_tiles <- data.frame(
    tile   = basename(tile_paths),
    status = TRUE
  )
  
  for (i in seq_along(tile_paths)) {
    # i = 1
    tile_path         <- tile_paths[i]
    tif_save_path     <- tif_save_paths[i]
    
    message("[", i, "]", 
            "\n > Tile: ", basename(tile_path),
            "\n > Output path: ", tif_save_path
    )
    
    download_tif_cmd    <- paste0("curl -o ", tif_save_path, " ", tile_path)
    
    tryCatch({
      
      tif_download_output <- system(download_tif_cmd, intern = TRUE)
      message("   > Succesfully downloaded tile: ", basename(tile_path))
      
    }, error = function(e) {
      
      message("Error downloading tile: ", basename(tile_path))
      message("ERROW below: \n  ", e)
      
      error_tiles[error_tiles$tile == basename(tile_path), ] <- FALSE
      
    })
    
  }
  
  return(error_tiles)
  
}

# Given 2 character vectors of filenames both including VPU strings after a "nextgen_" string, match them together to
# make sure they are aligned and in the same order
# x is a character vector of file paths with a VPU ID preceeded by a "nextgen_" string 
# y is a character vector of file paths with a VPU ID preceeded by a "nextgen_" string 
# base is a character vector of the base directory of the files. Defaults to NULL
# Returns a dataframe with VPU, x, and y columns
align_files_by_vpu <- function(
    x, 
    y, 
    base = NULL
) {

  # Regular expression pattern to match numeric pattern after "nextgen_" and remove everything after the ending period
  regex_pattern <- "nextgen_(\\d+[A-Za-z]?).*"
  
  # path dataframe for X filepaths
  x_paths <- data.frame(x = x)
  
  # path dataframe for Y filepaths
  y_paths <- data.frame(y = y)
  
  # generate VPU IDs based on file path regular expression matching with "regex_pattern" above
  x_paths$vpu <- gsub(regex_pattern, "\\1", x_paths$x)
  y_paths$vpu <- gsub(regex_pattern, "\\1", y_paths$y)
  
  # match paths based on VPU column
  matched_paths <- dplyr::left_join(
    x_paths,
    y_paths,
    by = "vpu"
  ) 
  
  # reorder columns
  matched_paths <- dplyr::relocate(matched_paths, vpu, x, y)
  
  if(!is.null(base)) {
    matched_paths$base_dir <- base
  }
  
  return(matched_paths)
  
}

# -------------------------------------------------------------------------------------
# FEMA processing functions: 
# -------------------------------------------------------------------------------------

resolve_internal_fema_boundaries <- function(fema, source_file = "") {
#   message("Resolving internal boundaries, islands, and topology issues:\n > '", basename(file_path), "'")
  
#   fema <- sf::read_sf(file_path)
  
  fema <-
    fema[!sf::st_is_empty(fema), ] %>% 
    sf::st_transform(5070)
  
  fema <-
    fema %>% 
    dplyr::select(geometry = geom) %>%
    add_predicate_group_id(sf::st_intersects) %>% 
    sf::st_make_valid() %>% 
    dplyr::group_by(group_id) %>% 
    dplyr::summarise(
      geometry = sf::st_combine(sf::st_union(geometry))
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-group_id) %>% 
    add_predicate_group_id(sf::st_intersects) %>% 
    rmapshaper::ms_dissolve(sys = TRUE, sys_mem = 16) %>% 
    rmapshaper::ms_explode(sys = TRUE, sys_mem = 16) %>% 
    dplyr::mutate(
      fema_id = as.character(1:dplyr::n())
    ) %>% 
    dplyr::select(fema_id, geometry)
  
  fema <- 
    fema %>% 
    dplyr::mutate(
      source = basename(source_file),
      state  = gsub("-100yr-flood_valid_clean.gpkg", "", source)
    ) %>%
    dplyr::select(fema_id, source, state, geometry)
  
  return(fema)
}

add_predicate_group_id <- function(polys, predicate) {
  # GROUP BY SPATIAL PREDICATES
  # ----------------------------------------- 
  # predicate = sf::st_touches
  # polys <- sf_df
  # ----------------------------------------- 
  
  
  relations <- predicate(polys)
  
  relations <- lapply(seq_along(relations), function(i) { as.character(sort(unique(c(relations[i][[1]], i)))) })
  
  group_ids_map <- fastmap::fastmap()
  ids_to_groups <- fastmap::fastmap()
  
  group_id <- 0
  
  for (i in seq_along(relations)) {
    
    predicate_ids <- relations[i][[1]]
    
    # message("(", i, ") - ", predicate_ids)
    # message("Start Group ID: ", group_id)
    
    id_group_check <- ids_to_groups$has(predicate_ids)
    
    if(any(id_group_check)) {
      
      known_groups  <- ids_to_groups$mget(predicate_ids)
      known_group   <- known_groups[unname(sapply(known_groups , function(kg) {
        !is.null(kg)
      }))][[1]]
      
      # message("IDs part of past group ID > '", known_group, "'")
      
      past_group_ids     <- group_ids_map$get(known_group)[[1]]
      updated_group_ids  <- as.character(
        sort(as.numeric(unique(c(past_group_ids, predicate_ids))))
      )
      
      group_ids_map$set(known_group, list(updated_group_ids))
      
      new_ids <- predicate_ids[!predicate_ids %in% past_group_ids]
      
      # message("Adding ", new_ids, " to seen set...")
      
      # add any newly added IDs to the seen map
      for (seen_id in new_ids) {
        # message(seen_id)
        ids_to_groups$set(as.character(seen_id), as.character(group_id))
      }
      
    } else {
      # get a new group ID number
      group_id <- group_id + 1    
      # message("IDs form NEW group > '", group_id, "'")
      
      # create a new key in the map with the predicate IDs list as the value
      group_ids_map$set(as.character(group_id), list(predicate_ids))
      
      # message("Adding ", predicate_ids, " to seen set...")
      
      # add each predicate ID to the map storing the seen indexes and their respecitve group IDs 
      for (seen_id in predicate_ids) {
        # message(seen_id)
        ids_to_groups$set(as.character(seen_id), as.character(group_id))
      }
    }
    # message("End group ID: ", group_id, "\n") 
  }
  
  group_ids   <- group_ids_map$as_list() 
  
  grouping_df <- lapply(seq_along(group_ids), function(i) {
    # i = 2
    grouping  <- group_ids[i] 
    group_id  <- names(grouping)
    indices   <- grouping[[1]][[1]]
    
    data.frame(
      index      = as.numeric(indices),
      group_id   = rep(group_id, length(indices))   
    )
    
  }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::arrange(i) 
  
  # count up the number of IDs for each group, well use this to determine which group 
  # to put any indices that had MULTIPLE groups they were apart of (use the group with the most other members)
  group_id_counts <- 
    grouping_df %>% 
    dplyr::group_by(group_id) %>% 
    dplyr::count() %>% 
    # dplyr::arrange(-n) %>% 
    dplyr::ungroup()
  
  # select the IDs with the most other members
  grouping_df <- 
    grouping_df %>% 
    dplyr::left_join(
      group_id_counts, 
      by = 'group_id'
    ) %>% 
    dplyr::group_by(index) %>% 
    dplyr::slice_max(n, with_ties = FALSE) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-n) %>% 
    dplyr::arrange(-index) 
  
  polys$group_id <- grouping_df$group_id
  
  return(polys)
  
}

# Update flowlines and transects to remove flowlines and transects that intersect with reference_features waterbodies
# flowlines: flowlines linestring sf object
# trans: transects linestring sf object
# waterbodies: waterbodies polygon sf object
# Returns a list of length 2 with logical vectors that subsets the "flowlines" and "transects" sf objects to remove flowlines and transects that intersect waterbodies
### Returns a list of length 2 with updated "flowlines" and "transects" sf objects
wb_intersects <- function(flowlines, trans, waterbodies) {
  
  ########  ########  ########  ########  ########  ########
  
  flowlines_geos <- geos::as_geos_geometry(flowlines)
  wbs_geos <- geos::as_geos_geometry(waterbodies)
  
  # temporary ID for transects that is the "hy_id", underscore, "cs_id", used for subsetting in future steps
  trans$tmp_id <- paste0(trans$hy_id, "_", trans$cs_id)
  
  message("Checking flowlines against waterbodies...")
  
  # create an index between flowlines and waterbodies 
  wb_index <- geos::geos_intersects_matrix(flowlines_geos, wbs_geos)
  
  # remove any flowlines that cross more than 1 waterbody
  to_keep  <- flowlines[lengths(wb_index) == 0, ]
  to_check <- flowlines[lengths(wb_index) != 0, ]
  
  # subset transects to the hy_ids in "to_check" set of flowlines
  trans_check <- trans[trans$hy_id %in% unique(to_check$id), ]
  # trans_check <- trans_geos[trans$hy_id %in% unique(to_check$id)]
  
  # check where the transects linestrings intersect with the waterbodies
  trans_geos_check <- geos::as_geos_geometry(trans_check)
  
  message("Checking transects against waterbodies (v2) ...")
  wb_trans_index <- geos::geos_intersects_matrix(trans_geos_check, wbs_geos)                    # (NEW METHOD)
  # wb_trans_index <- geos::geos_intersects_any(trans_geos_check, wbs_geos[unlist(wb_index)])   # (OLD METHOD)
  
  # sum(lengths(wb_trans_index) == 0)
  # length(wb_trans_index)
  
  # within the transects lines that are on a flowline that crosses a waterbody, 
  # check if any of these transects line DO NOT CROSS A WATERBODY AT ALL
  trans_keep <- trans_check[lengths(wb_trans_index) == 0, ]                       # (NEW METHOD)
  # trans_keep <- trans_check[!wb_trans_index, ]                                  # (OLD METHOD)
  
  # preserve any flowlines that CROSS A WATERBODY BUT ALSO HAVE A TRANSECT LINE that does NOT cross any waterbodies
  to_check <- to_check[to_check$id %in% unique(trans_keep$hy_id), ]
  
  # update flowlines to keep with flowlines that intersect a waterbody BUT STILL,
  # have transects that are NOT in the waterbody
  to_keep <- dplyr::bind_rows(to_keep, to_check)
  
  # 'tmp_ids' of transects that are being checked and also the transects within trans_check 
  # that were determined to be valid (are being kept)
  check_ids <- unique(trans_check$tmp_id)
  keep_ids <- unique(trans_keep$tmp_id)
  
  # logical vectors of which flowlines/transects to keep (KEEP == TRUE)
  # - Remove any transects that are on flowlines that cross a waterbody AND the transect crosses the waterbody too.
  # - Keep original transects that are not on flowlines that intersect waterbodies AND 
  # also the transects that do NOT intersect waterbodies but are on a flowline that DOES intersect a waterbody
  valid_flowlines <- flowlines$id %in% to_keep$id
  valid_transects <- trans$tmp_id %in% dplyr::filter(trans,
                                                     !tmp_id %in% check_ids[!check_ids %in% keep_ids])$tmp_id
  
  # return alist of updated flowlines and transects 
  return(
    list(
      "valid_flowlines" =  valid_flowlines,
      "valid_transects" =  valid_transects
    )
  )
  
  # # within the transects lines that are on a flowline that crosses a waterbody, 
  # # check if any of these transects line DO NOT CROSS A WATERBODY AT ALL
  # trans_keep <- trans_check[!trans_wb_index, ]
  # # trans_keep <- trans_check[lengths(trans_wb_index2) == 0, ]
  # 
  # # preserve any flowlines that CROSS A WATERBODY BUT ALSO HAVE A TRANSECT LINE that does NOT cross any waterbodies
  # to_check <- to_check[to_check$id %in% unique(trans_keep$hy_id), ]
  # 
  # # update flowlines to keep with flowlines that intersect a waterbody BUT STILL,
  # # have transects that are NOT in the waterbody
  # to_keep <- dplyr::bind_rows(to_keep, to_check)
  # 
  # # 'tmp_ids' of transects that are being checked and also the transects within trans_check 
  # # that were determined to be valid (are being kept)
  # check_ids <- unique(trans_check$tmp_id)
  # keep_ids <- unique(trans_keep$tmp_id)
  # 
  # # logical vectors of which flowlines/transects to keep (KEEP == TRUE)
  # # - Remove any transects that are on flowlines that cross a waterbody AND the transect crosses the waterbody too.
  # # - Keep original transects that are not on flowlines that intersect waterbodies AND 
  # # also the transects that do NOT intersect waterbodies but are on a flowline that DOES intersect a waterbody
  # valid_flowlines <- flowlines$id %in% to_keep$id
  # valid_transects <- trans$tmp_id %in% dplyr::filter(trans,
  #                                                    !tmp_id %in% check_ids[!check_ids %in% keep_ids])$tmp_id
  # 
  # # return alist of updated flowlines and transects 
  # return(
  #   list(
  #     "valid_flowlines" =  valid_flowlines,
  #     "valid_transects" =  valid_transects
  #   )
  # )
}

add_intersects_ids <- function(x, y, id_col) {
  # make sure the crs are tjhe same
  y <- sf::st_transform(y, sf::st_crs(x))
  
  # Perform the intersection
  intersections <- sf::st_intersects(x, y)
  
  # add the intersected values to the first dataframe
  x[[id_col]] <- unlist(lapply(intersections, function(idx) {
    if (length(idx) > 0) {
      paste0(unlist(y[[id_col]][idx]), collapse = ", ")
    } else {
      NA
    }
  }))
  
  return(x)
}

unnest_ids <- function(ids) {
  return(
    unique(unlist(strsplit(unique(ids), ", ")))
  )
}

#' Get the polygons that interesect with any of the linestring geometries
#' This is just a wrapper around geos::geos_intersects_matrix. Takes in sf dataframes, uses geos, then outputs sf dataframes
#' @param polygons polygon sf object. Default is NULL
#' @param lines linestring sf object. Default is NULL.
#'
#' @return sf dataframe of polygons that intersect with the linestrings
polygons_with_line_intersects <- function(polygons = NULL, lines = NULL) {
  
  if (is.null(polygons)) {
    stop("NULL 'polygons' argument, provide an sf dataframe of POLYGON or MULTIPOLYGON geometries")
  }
  
  if (is.null(lines)) {
    stop("NULL 'lines' argument, provide an sf dataframe of LINESTRING or MULTILINESTRING geometries")
  }
  
  # Convert the SF geometries to geos geometries
  polygons_geos   <- geos::as_geos_geometry(polygons)
  lines_geos      <- geos::as_geos_geometry(lines)
  
  # create an index between the polygons and linestrings
  lines_index <-  geos::geos_intersects_matrix(polygons_geos, lines_geos)
  
  # get the polygons that have atleast 1 intersection with the 'lines'
  polygons_with_lines <- polygons[lengths(lines_index) != 0, ]
  
  return(polygons_with_lines)
}

# TODO: DElete these NEW DOMAIN functions...
# Create an empty file structure 
# base_dir: character, top level directory path
# domain_dirname: character, name of the intended new domain directory, if folder exists, then the required subdirectories are created (if they DO NOT exist)

# Directory tree:
# base_dir/
#   └── domain_dirname/
#     ├── flowlines/
#     ├── dem/
#     ├── transects/
#     ├── cross_sections/
#     └── cs_pts/
create_new_domain_dirs <- function(base_dir, domain_dirname, with_output = FALSE) {
  
  # build paths
  domain_dir         <- paste0(base_dir, "/", domain_dirname)
  flowlines_dir      <- paste0(domain_dir, "/flowlines")
  domain_subset_dir  <- paste0(domain_dir, "/domain_subset")
  dem_dir            <- paste0(domain_dir, "/dem")
  transects_dir      <- paste0(domain_dir, "/transects")
  cross_sections_dir <- paste0(domain_dir, "/cross_sections")
  cs_pts_dir         <- paste0(domain_dir, "/cs_pts")
  vpu_subsets_dir    <- paste0(domain_dir, "/vpu-subsets")
  
  if(with_output) {
    output_dir       <- paste0(domain_dir, "/outputs")
  }
  
  # create directories
  create_if_not_exists(domain_dir)
  create_if_not_exists(flowlines_dir)
  create_if_not_exists(domain_subset_dir)
  create_if_not_exists(dem_dir)
  create_if_not_exists(transects_dir)
  create_if_not_exists(cross_sections_dir)
  create_if_not_exists(cs_pts_dir)
  create_if_not_exists(vpu_subsets_dir)
  
  if(with_output) {
    create_if_not_exists(output_dir)
  }
  
}

# get path strings for a domain dir (based of a base dir and domain dirname)
# NOTE: this does NOT guarentee that these folders exist, 
# NOTE: it just gets the paths if they were created by create_new_domain_dirs() 
get_new_domain_paths <- function(base_dir, domain_dirname, with_output = FALSE) {
  
  # build paths
  domain_dir         <- paste0(base_dir, "/", domain_dirname)
  flowlines_dir      <- paste0(domain_dir, "/flowlines")
  domain_subset_dir  <- paste0(domain_dir, "/domain_subset")
  dem_dir            <- paste0(domain_dir, "/dem")
  transects_dir      <- paste0(domain_dir, "/transects")
  cross_sections_dir <- paste0(domain_dir, "/cross_sections")
  cs_pts_dir         <- paste0(domain_dir, "/cs_pts")
  vpu_subsets_dir    <- paste0(domain_dir, "/vpu-subsets")
  output_dir         <- ifelse(with_output, paste0(domain_dir, "/outputs"), NA)
  
  
  # named list of file paths
  return(
    list(
      base_dir           = base_dir, 
      domain_dir         = domain_dir,
      flowlines_dir      = flowlines_dir,
      domain_subset_dir  = domain_subset_dir,
      dem_dir            = dem_dir,
      transects_dir      = transects_dir,
      cross_sections_dir = cross_sections_dir,
      cs_pts_dir         = cs_pts_dir,
      vpu_subsets_dir    = vpu_subsets_dir,
      output_dir         = output_dir
    )
  )
  
}

download_3dep_vrt <- function(base_dir) {
  
  ## Cmd+A/Cmd+C from: http://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Elevation/13/TIFF/current/
  ### paste w/ `datapasta::vector_paste_vertical()`
  ### Some manual cleaning of header and footer mess...
  ### Reason? Un-scrapable page, and no index.gpkg...
  
  t <- c(
    "0              n06e162/",
    "0              n06e163/",
    "0              n07e134/",
    "0              n07e151/",
    "0              n07e152/",
    "0              n07e158/",
    "0              n08e134/",
    "0              n08e151/",
    "0              n08e152/",
    "0              n08e158/",
    "0              n09e134/",
    "0              n10e138/",
    "0              n14e144/",
    "0              n15e145/",
    "0              n16e145/",
    "0              n18w065/",
    "0              n18w066/",
    "0              n18w067/",
    "0              n18w068/",
    "0              n19w065/",
    "0              n19w066/",
    "0              n19w067/",
    "0              n19w068/",
    "0              n19w156/",
    "0              n20w155/",
    "0              n20w156/",
    "0              n20w157/",
    "0              n21w156/",
    "0              n21w157/",
    "0              n21w158/",
    "0              n22w157/",
    "0              n22w158/",
    "0              n22w159/",
    "0              n22w160/",
    "0              n22w161/",
    "0              n23w160/",
    "0              n23w161/",
    "0              n25w081/",
    "0              n25w082/",
    "0              n25w083/",
    "0              n26w081/",
    "0              n26w082/",
    "0              n26w098/",
    "0              n26w099/",
    "0              n27w081/",
    "0              n27w082/",
    "0              n27w083/",
    "0              n27w098/",
    "0              n27w099/",
    "0              n27w100/",
    "0              n28w081/",
    "0              n28w082/",
    "0              n28w083/",
    "0              n28w097/",
    "0              n28w098/",
    "0              n28w099/",
    "0              n28w100/",
    "0              n28w101/",
    "0              n29w081/",
    "0              n29w082/",
    "0              n29w083/",
    "0              n29w090/",
    "0              n29w096/",
    "0              n29w097/",
    "0              n29w098/",
    "0              n29w099/",
    "0              n29w100/",
    "0              n29w101/",
    "0              n29w104/",
    "0              n30w081/",
    "0              n30w082/",
    "0              n30w083/",
    "0              n30w084/",
    "0              n30w085/",
    "0              n30w086/",
    "0              n30w089/",
    "0              n30w090/",
    "0              n30w091/",
    "0              n30w092/",
    "0              n30w093/",
    "0              n30w094/",
    "0              n30w095/",
    "0              n30w096/",
    "0              n30w097/",
    "0              n30w098/",
    "0              n30w099/",
    "0              n30w100/",
    "0              n30w101/",
    "0              n30w102/",
    "0              n30w103/",
    "0              n30w104/",
    "0              n30w105/",
    "0              n31w082/",
    "0              n31w083/",
    "0              n31w084/",
    "0              n31w085/",
    "0              n31w086/",
    "0              n31w087/",
    "0              n31w088/",
    "0              n31w089/",
    "0              n31w090/",
    "0              n31w091/",
    "0              n31w092/",
    "0              n31w093/",
    "0              n31w094/",
    "0              n31w095/",
    "0              n31w096/",
    "0              n31w097/",
    "0              n31w098/",
    "0              n31w099/",
    "0              n31w100/",
    "0              n31w101/",
    "0              n31w102/",
    "0              n31w103/",
    "0              n31w104/",
    "0              n31w105/",
    "0              n31w106/",
    "0              n31w107/",
    "0              n32w081/",
    "0              n32w082/",
    "0              n32w083/",
    "0              n32w084/",
    "0              n32w085/",
    "0              n32w086/",
    "0              n32w087/",
    "0              n32w088/",
    "0              n32w089/",
    "0              n32w090/",
    "0              n32w091/",
    "0              n32w092/",
    "0              n32w093/",
    "0              n32w094/",
    "0              n32w095/",
    "0              n32w096/",
    "0              n32w097/",
    "0              n32w098/",
    "0              n32w099/",
    "0              n32w100/",
    "0              n32w101/",
    "0              n32w102/",
    "0              n32w103/",
    "0              n32w104/",
    "0              n32w105/",
    "0              n32w106/",
    "0              n32w107/",
    "0              n32w108/",
    "0              n32w109/",
    "0              n32w110/",
    "0              n32w111/",
    "0              n32w112/",
    "0              n32w113/",
    "0              n32w114/",
    "0              n33w080/",
    "0              n33w081/",
    "0              n33w082/",
    "0              n33w083/",
    "0              n33w084/",
    "0              n33w085/",
    "0              n33w086/",
    "0              n33w087/",
    "0              n33w088/",
    "0              n33w089/",
    "0              n33w090/",
    "0              n33w091/",
    "0              n33w092/",
    "0              n33w093/",
    "0              n33w094/",
    "0              n33w095/",
    "0              n33w096/",
    "0              n33w097/",
    "0              n33w098/",
    "0              n33w099/",
    "0              n33w100/",
    "0              n33w101/",
    "0              n33w102/",
    "0              n33w103/",
    "0              n33w104/",
    "0              n33w105/",
    "0              n33w106/",
    "0              n33w107/",
    "0              n33w108/",
    "0              n33w109/",
    "0              n33w110/",
    "0              n33w111/",
    "0              n33w112/",
    "0              n33w113/",
    "0              n33w114/",
    "0              n33w115/",
    "0              n33w116/",
    "0              n33w117/",
    "0              n33w118/",
    "0              n33w119/",
    "0              n34w078/",
    "0              n34w079/",
    "0              n34w080/",
    "0              n34w081/",
    "0              n34w082/",
    "0              n34w083/",
    "0              n34w084/",
    "0              n34w085/",
    "0              n34w086/",
    "0              n34w087/",
    "0              n34w088/",
    "0              n34w089/",
    "0              n34w090/",
    "0              n34w091/",
    "0              n34w092/",
    "0              n34w093/",
    "0              n34w094/",
    "0              n34w095/",
    "0              n34w096/",
    "0              n34w097/",
    "0              n34w098/",
    "0              n34w099/",
    "0              n34w100/",
    "0              n34w101/",
    "0              n34w102/",
    "0              n34w103/",
    "0              n34w104/",
    "0              n34w105/",
    "0              n34w106/",
    "0              n34w107/",
    "0              n34w108/",
    "0              n34w109/",
    "0              n34w110/",
    "0              n34w111/",
    "0              n34w112/",
    "0              n34w113/",
    "0              n34w114/",
    "0              n34w115/",
    "0              n34w116/",
    "0              n34w117/",
    "0              n34w118/",
    "0              n34w119/",
    "0              n34w120/",
    "0              n34w121/",
    "0              n35w076/",
    "0              n35w077/",
    "0              n35w078/",
    "0              n35w079/",
    "0              n35w080/",
    "0              n35w081/",
    "0              n35w082/",
    "0              n35w083/",
    "0              n35w084/",
    "0              n35w085/",
    "0              n35w086/",
    "0              n35w087/",
    "0              n35w088/",
    "0              n35w089/",
    "0              n35w090/",
    "0              n35w091/",
    "0              n35w092/",
    "0              n35w093/",
    "0              n35w094/",
    "0              n35w095/",
    "0              n35w096/",
    "0              n35w097/",
    "0              n35w098/",
    "0              n35w099/",
    "0              n35w100/",
    "0              n35w101/",
    "0              n35w102/",
    "0              n35w103/",
    "0              n35w104/",
    "0              n35w105/",
    "0              n35w106/",
    "0              n35w107/",
    "0              n35w108/",
    "0              n35w109/",
    "0              n35w110/",
    "0              n35w111/",
    "0              n35w112/",
    "0              n35w113/",
    "0              n35w114/",
    "0              n35w115/",
    "0              n35w116/",
    "0              n35w117/",
    "0              n35w118/",
    "0              n35w119/",
    "0              n35w120/",
    "0              n35w121/",
    "0              n36w076/",
    "0              n36w077/",
    "0              n36w078/",
    "0              n36w079/",
    "0              n36w080/",
    "0              n36w081/",
    "0              n36w082/",
    "0              n36w083/",
    "0              n36w084/",
    "0              n36w085/",
    "0              n36w086/",
    "0              n36w087/",
    "0              n36w088/",
    "0              n36w089/",
    "0              n36w090/",
    "0              n36w091/",
    "0              n36w092/",
    "0              n36w093/",
    "0              n36w094/",
    "0              n36w095/",
    "0              n36w096/",
    "0              n36w097/",
    "0              n36w098/",
    "0              n36w099/",
    "0              n36w100/",
    "0              n36w101/",
    "0              n36w102/",
    "0              n36w103/",
    "0              n36w104/",
    "0              n36w105/",
    "0              n36w106/",
    "0              n36w107/",
    "0              n36w108/",
    "0              n36w109/",
    "0              n36w110/",
    "0              n36w111/",
    "0              n36w112/",
    "0              n36w113/",
    "0              n36w114/",
    "0              n36w115/",
    "0              n36w116/",
    "0              n36w117/",
    "0              n36w118/",
    "0              n36w119/",
    "0              n36w120/",
    "0              n36w121/",
    "0              n36w122/",
    "0              n37w076/",
    "0              n37w077/",
    "0              n37w078/",
    "0              n37w079/",
    "0              n37w080/",
    "0              n37w081/",
    "0              n37w082/",
    "0              n37w083/",
    "0              n37w084/",
    "0              n37w085/",
    "0              n37w086/",
    "0              n37w087/",
    "0              n37w088/",
    "0              n37w089/",
    "0              n37w090/",
    "0              n37w091/",
    "0              n37w092/",
    "0              n37w093/",
    "0              n37w094/",
    "0              n37w095/",
    "0              n37w096/",
    "0              n37w097/",
    "0              n37w098/",
    "0              n37w099/",
    "0              n37w100/",
    "0              n37w101/",
    "0              n37w102/",
    "0              n37w103/",
    "0              n37w104/",
    "0              n37w105/",
    "0              n37w106/",
    "0              n37w107/",
    "0              n37w108/",
    "0              n37w109/",
    "0              n37w110/",
    "0              n37w111/",
    "0              n37w112/",
    "0              n37w113/",
    "0              n37w114/",
    "0              n37w115/",
    "0              n37w116/",
    "0              n37w117/",
    "0              n37w118/",
    "0              n37w119/",
    "0              n37w120/",
    "0              n37w121/",
    "0              n37w122/",
    "0              n37w123/",
    "0              n38w076/",
    "0              n38w077/",
    "0              n38w078/",
    "0              n38w079/",
    "0              n38w080/",
    "0              n38w081/",
    "0              n38w082/",
    "0              n38w083/",
    "0              n38w084/",
    "0              n38w085/",
    "0              n38w086/",
    "0              n38w087/",
    "0              n38w088/",
    "0              n38w089/",
    "0              n38w090/",
    "0              n38w091/",
    "0              n38w092/",
    "0              n38w093/",
    "0              n38w094/",
    "0              n38w095/",
    "0              n38w096/",
    "0              n38w097/",
    "0              n38w098/",
    "0              n38w099/",
    "0              n38w100/",
    "0              n38w101/",
    "0              n38w102/",
    "0              n38w103/",
    "0              n38w104/",
    "0              n38w105/",
    "0              n38w106/",
    "0              n38w107/",
    "0              n38w108/",
    "0              n38w109/",
    "0              n38w110/",
    "0              n38w111/",
    "0              n38w112/",
    "0              n38w113/",
    "0              n38w114/",
    "0              n38w115/",
    "0              n38w116/",
    "0              n38w117/",
    "0              n38w118/",
    "0              n38w119/",
    "0              n38w120/",
    "0              n38w121/",
    "0              n38w122/",
    "0              n38w123/",
    "0              n38w124/",
    "0              n39w075/",
    "0              n39w076/",
    "0              n39w077/",
    "0              n39w078/",
    "0              n39w079/",
    "0              n39w080/",
    "0              n39w081/",
    "0              n39w082/",
    "0              n39w083/",
    "0              n39w084/",
    "0              n39w085/",
    "0              n39w086/",
    "0              n39w087/",
    "0              n39w088/",
    "0              n39w089/",
    "0              n39w090/",
    "0              n39w091/",
    "0              n39w092/",
    "0              n39w093/",
    "0              n39w094/",
    "0              n39w095/",
    "0              n39w096/",
    "0              n39w097/",
    "0              n39w098/",
    "0              n39w099/",
    "0              n39w100/",
    "0              n39w101/",
    "0              n39w102/",
    "0              n39w103/",
    "0              n39w104/",
    "0              n39w105/",
    "0              n39w106/",
    "0              n39w107/",
    "0              n39w108/",
    "0              n39w109/",
    "0              n39w110/",
    "0              n39w111/",
    "0              n39w112/",
    "0              n39w113/",
    "0              n39w114/",
    "0              n39w115/",
    "0              n39w116/",
    "0              n39w117/",
    "0              n39w118/",
    "0              n39w119/",
    "0              n39w120/",
    "0              n39w121/",
    "0              n39w122/",
    "0              n39w123/",
    "0              n39w124/",
    "0              n40w075/",
    "0              n40w076/",
    "0              n40w077/",
    "0              n40w078/",
    "0              n40w079/",
    "0              n40w080/",
    "0              n40w081/",
    "0              n40w082/",
    "0              n40w083/",
    "0              n40w084/",
    "0              n40w085/",
    "0              n40w086/",
    "0              n40w087/",
    "0              n40w088/",
    "0              n40w089/",
    "0              n40w090/",
    "0              n40w091/",
    "0              n40w092/",
    "0              n40w093/",
    "0              n40w094/",
    "0              n40w095/",
    "0              n40w096/",
    "0              n40w097/",
    "0              n40w098/",
    "0              n40w099/",
    "0              n40w100/",
    "0              n40w101/",
    "0              n40w102/",
    "0              n40w103/",
    "0              n40w104/",
    "0              n40w105/",
    "0              n40w106/",
    "0              n40w107/",
    "0              n40w108/",
    "0              n40w109/",
    "0              n40w110/",
    "0              n40w111/",
    "0              n40w112/",
    "0              n40w113/",
    "0              n40w114/",
    "0              n40w115/",
    "0              n40w116/",
    "0              n40w117/",
    "0              n40w118/",
    "0              n40w119/",
    "0              n40w120/",
    "0              n40w121/",
    "0              n40w122/",
    "0              n40w123/",
    "0              n40w124/",
    "0              n40w125/",
    "0              n41w073/",
    "0              n41w074/",
    "0              n41w075/",
    "0              n41w076/",
    "0              n41w077/",
    "0              n41w078/",
    "0              n41w079/",
    "0              n41w080/",
    "0              n41w081/",
    "0              n41w082/",
    "0              n41w083/",
    "0              n41w084/",
    "0              n41w085/",
    "0              n41w086/",
    "0              n41w087/",
    "0              n41w088/",
    "0              n41w089/",
    "0              n41w090/",
    "0              n41w091/",
    "0              n41w092/",
    "0              n41w093/",
    "0              n41w094/",
    "0              n41w095/",
    "0              n41w096/",
    "0              n41w097/",
    "0              n41w098/",
    "0              n41w099/",
    "0              n41w100/",
    "0              n41w101/",
    "0              n41w102/",
    "0              n41w103/",
    "0              n41w104/",
    "0              n41w105/",
    "0              n41w106/",
    "0              n41w107/",
    "0              n41w108/",
    "0              n41w109/",
    "0              n41w110/",
    "0              n41w111/",
    "0              n41w112/",
    "0              n41w113/",
    "0              n41w114/",
    "0              n41w115/",
    "0              n41w116/",
    "0              n41w117/",
    "0              n41w118/",
    "0              n41w119/",
    "0              n41w120/",
    "0              n41w121/",
    "0              n41w122/",
    "0              n41w123/",
    "0              n41w124/",
    "0              n41w125/",
    "0              n42w070/",
    "0              n42w071/",
    "0              n42w072/",
    "0              n42w073/",
    "0              n42w074/",
    "0              n42w075/",
    "0              n42w076/",
    "0              n42w077/",
    "0              n42w078/",
    "0              n42w079/",
    "0              n42w080/",
    "0              n42w081/",
    "0              n42w082/",
    "0              n42w083/",
    "0              n42w084/",
    "0              n42w085/",
    "0              n42w086/",
    "0              n42w087/",
    "0              n42w088/",
    "0              n42w089/",
    "0              n42w090/",
    "0              n42w091/",
    "0              n42w092/",
    "0              n42w093/",
    "0              n42w094/",
    "0              n42w095/",
    "0              n42w096/",
    "0              n42w097/",
    "0              n42w098/",
    "0              n42w099/",
    "0              n42w100/",
    "0              n42w101/",
    "0              n42w102/",
    "0              n42w103/",
    "0              n42w104/",
    "0              n42w105/",
    "0              n42w106/",
    "0              n42w107/",
    "0              n42w108/",
    "0              n42w109/",
    "0              n42w110/",
    "0              n42w111/",
    "0              n42w112/",
    "0              n42w113/",
    "0              n42w114/",
    "0              n42w115/",
    "0              n42w116/",
    "0              n42w117/",
    "0              n42w118/",
    "0              n42w119/",
    "0              n42w120/",
    "0              n42w121/",
    "0              n42w122/",
    "0              n42w123/",
    "0              n42w124/",
    "0              n42w125/",
    "0              n43w071/",
    "0              n43w072/",
    "0              n43w073/",
    "0              n43w074/",
    "0              n43w075/",
    "0              n43w076/",
    "0              n43w077/",
    "0              n43w078/",
    "0              n43w079/",
    "0              n43w080/",
    "0              n43w081/",
    "0              n43w082/",
    "0              n43w083/",
    "0              n43w084/",
    "0              n43w085/",
    "0              n43w086/",
    "0              n43w087/",
    "0              n43w088/",
    "0              n43w089/",
    "0              n43w090/",
    "0              n43w091/",
    "0              n43w092/",
    "0              n43w093/",
    "0              n43w094/",
    "0              n43w095/",
    "0              n43w096/",
    "0              n43w097/",
    "0              n43w098/",
    "0              n43w099/",
    "0              n43w100/",
    "0              n43w101/",
    "0              n43w102/",
    "0              n43w103/",
    "0              n43w104/",
    "0              n43w105/",
    "0              n43w106/",
    "0              n43w107/",
    "0              n43w108/",
    "0              n43w109/",
    "0              n43w110/",
    "0              n43w111/",
    "0              n43w112/",
    "0              n43w113/",
    "0              n43w114/",
    "0              n43w115/",
    "0              n43w116/",
    "0              n43w117/",
    "0              n43w118/",
    "0              n43w119/",
    "0              n43w120/",
    "0              n43w121/",
    "0              n43w122/",
    "0              n43w123/",
    "0              n43w124/",
    "0              n43w125/",
    "0              n44w069/",
    "0              n44w070/",
    "0              n44w071/",
    "0              n44w072/",
    "0              n44w073/",
    "0              n44w074/",
    "0              n44w075/",
    "0              n44w076/",
    "0              n44w077/",
    "0              n44w078/",
    "0              n44w079/",
    "0              n44w080/",
    "0              n44w081/",
    "0              n44w083/",
    "0              n44w084/",
    "0              n44w085/",
    "0              n44w086/",
    "0              n44w087/",
    "0              n44w088/",
    "0              n44w089/",
    "0              n44w090/",
    "0              n44w091/",
    "0              n44w092/",
    "0              n44w093/",
    "0              n44w094/",
    "0              n44w095/",
    "0              n44w096/",
    "0              n44w097/",
    "0              n44w098/",
    "0              n44w099/",
    "0              n44w100/",
    "0              n44w101/",
    "0              n44w102/",
    "0              n44w103/",
    "0              n44w104/",
    "0              n44w105/",
    "0              n44w106/",
    "0              n44w107/",
    "0              n44w108/",
    "0              n44w109/",
    "0              n44w110/",
    "0              n44w111/",
    "0              n44w112/",
    "0              n44w113/",
    "0              n44w114/",
    "0              n44w115/",
    "0              n44w116/",
    "0              n44w117/",
    "0              n44w118/",
    "0              n44w119/",
    "0              n44w120/",
    "0              n44w121/",
    "0              n44w122/",
    "0              n44w123/",
    "0              n44w124/",
    "0              n44w125/",
    "0              n45w067/",
    "0              n45w068/",
    "0              n45w069/",
    "0              n45w070/",
    "0              n45w071/",
    "0              n45w072/",
    "0              n45w073/",
    "0              n45w074/",
    "0              n45w075/",
    "0              n45w076/",
    "0              n45w077/",
    "0              n45w083/",
    "0              n45w084/",
    "0              n45w085/",
    "0              n45w086/",
    "0              n45w087/",
    "0              n45w088/",
    "0              n45w089/",
    "0              n45w090/",
    "0              n45w091/",
    "0              n45w092/",
    "0              n45w093/",
    "0              n45w094/",
    "0              n45w095/",
    "0              n45w096/",
    "0              n45w097/",
    "0              n45w098/",
    "0              n45w099/",
    "0              n45w100/",
    "0              n45w101/",
    "0              n45w102/",
    "0              n45w103/",
    "0              n45w104/",
    "0              n45w105/",
    "0              n45w106/",
    "0              n45w107/",
    "0              n45w108/",
    "0              n45w109/",
    "0              n45w110/",
    "0              n45w111/",
    "0              n45w112/",
    "0              n45w113/",
    "0              n45w114/",
    "0              n45w115/",
    "0              n45w116/",
    "0              n45w117/",
    "0              n45w118/",
    "0              n45w119/",
    "0              n45w120/",
    "0              n45w121/",
    "0              n45w122/",
    "0              n45w123/",
    "0              n45w124/",
    "0              n45w125/",
    "0              n46w068/",
    "0              n46w069/",
    "0              n46w070/",
    "0              n46w071/",
    "0              n46w072/",
    "0              n46w073/",
    "0              n46w074/",
    "0              n46w075/",
    "0              n46w084/",
    "0              n46w085/",
    "0              n46w086/",
    "0              n46w087/",
    "0              n46w088/",
    "0              n46w089/",
    "0              n46w090/",
    "0              n46w091/",
    "0              n46w092/",
    "0              n46w093/",
    "0              n46w094/",
    "0              n46w095/",
    "0              n46w096/",
    "0              n46w097/",
    "0              n46w098/",
    "0              n46w099/",
    "0              n46w100/",
    "0              n46w101/",
    "0              n46w102/",
    "0              n46w103/",
    "0              n46w104/",
    "0              n46w105/",
    "0              n46w106/",
    "0              n46w107/",
    "0              n46w108/",
    "0              n46w109/",
    "0              n46w110/",
    "0              n46w111/",
    "0              n46w112/",
    "0              n46w113/",
    "0              n46w114/",
    "0              n46w115/",
    "0              n46w116/",
    "0              n46w117/",
    "0              n46w118/",
    "0              n46w119/",
    "0              n46w120/",
    "0              n46w121/",
    "0              n46w122/",
    "0              n46w123/",
    "0              n46w124/",
    "0              n46w125/",
    "0              n47w068/",
    "0              n47w069/",
    "0              n47w070/",
    "0              n47w071/",
    "0              n47w084/",
    "0              n47w085/",
    "0              n47w086/",
    "0              n47w087/",
    "0              n47w088/",
    "0              n47w089/",
    "0              n47w090/",
    "0              n47w091/",
    "0              n47w092/",
    "0              n47w093/",
    "0              n47w094/",
    "0              n47w095/",
    "0              n47w096/",
    "0              n47w097/",
    "0              n47w098/",
    "0              n47w099/",
    "0              n47w100/",
    "0              n47w101/",
    "0              n47w102/",
    "0              n47w103/",
    "0              n47w104/",
    "0              n47w105/",
    "0              n47w106/",
    "0              n47w107/",
    "0              n47w108/",
    "0              n47w109/",
    "0              n47w110/",
    "0              n47w111/",
    "0              n47w112/",
    "0              n47w113/",
    "0              n47w114/",
    "0              n47w115/",
    "0              n47w116/",
    "0              n47w117/",
    "0              n47w118/",
    "0              n47w119/",
    "0              n47w120/",
    "0              n47w121/",
    "0              n47w122/",
    "0              n47w123/",
    "0              n47w124/",
    "0              n47w125/",
    "0              n48w068/",
    "0              n48w069/",
    "0              n48w070/",
    "0              n48w087/",
    "0              n48w088/",
    "0              n48w089/",
    "0              n48w090/",
    "0              n48w091/",
    "0              n48w092/",
    "0              n48w093/",
    "0              n48w094/",
    "0              n48w095/",
    "0              n48w096/",
    "0              n48w097/",
    "0              n48w098/",
    "0              n48w099/",
    "0              n48w100/",
    "0              n48w101/",
    "0              n48w102/",
    "0              n48w103/",
    "0              n48w104/",
    "0              n48w105/",
    "0              n48w106/",
    "0              n48w107/",
    "0              n48w108/",
    "0              n48w109/",
    "0              n48w110/",
    "0              n48w111/",
    "0              n48w112/",
    "0              n48w113/",
    "0              n48w114/",
    "0              n48w115/",
    "0              n48w116/",
    "0              n48w117/",
    "0              n48w118/",
    "0              n48w119/",
    "0              n48w120/",
    "0              n48w121/",
    "0              n48w122/",
    "0              n48w123/",
    "0              n48w124/",
    "0              n48w125/",
    "0              n49w089/",
    "0              n49w090/",
    "0              n49w091/",
    "0              n49w092/",
    "0              n49w093/",
    "0              n49w094/",
    "0              n49w095/",
    "0              n49w096/",
    "0              n49w097/",
    "0              n49w098/",
    "0              n49w099/",
    "0              n49w100/",
    "0              n49w101/",
    "0              n49w102/",
    "0              n49w103/",
    "0              n49w104/",
    "0              n49w105/",
    "0              n49w106/",
    "0              n49w107/",
    "0              n49w108/",
    "0              n49w109/",
    "0              n49w110/",
    "0              n49w111/",
    "0              n49w112/",
    "0              n49w113/",
    "0              n49w114/",
    "0              n49w115/",
    "0              n49w116/",
    "0              n49w117/",
    "0              n49w118/",
    "0              n49w119/",
    "0              n49w120/",
    "0              n49w121/",
    "0              n49w122/",
    "0              n49w123/",
    "0              n49w124/",
    "0              n49w125/",
    "0              n50w095/",
    "0              n50w096/",
    "0              n50w097/",
    "0              n50w098/",
    "0              n50w099/",
    "0              n50w100/",
    "0              n50w101/",
    "0              n50w107/",
    "0              n50w108/",
    "0              n50w122/",
    "0              n50w123/",
    "0              n50w124/",
    "0              n52e177/",
    "0              n52e178/",
    "0              n52e179/",
    "0              n52w174/",
    "0              n52w176/",
    "0              n52w177/",
    "0              n52w178/",
    "0              n52w179/",
    "0              n52w180/",
    "0              n53e172/",
    "0              n53e173/",
    "0              n53e174/",
    "0              n53e175/",
    "0              n53e177/",
    "0              n53e178/",
    "0              n53e179/",
    "0              n53w169/",
    "0              n53w170/",
    "0              n53w171/",
    "0              n53w172/",
    "0              n53w173/",
    "0              n53w174/",
    "0              n53w175/",
    "0              n53w176/",
    "0              n53w177/",
    "0              n54e172/",
    "0              n54e173/",
    "0              n54w167/",
    "0              n54w168/",
    "0              n54w169/",
    "0              n54w170/",
    "0              n55w131/",
    "0              n55w132/",
    "0              n55w133/",
    "0              n55w134/",
    "0              n55w160/",
    "0              n55w161/",
    "0              n55w162/",
    "0              n55w163/",
    "0              n55w164/",
    "0              n55w165/",
    "0              n55w166/",
    "0              n55w167/",
    "0              n56w130/",
    "0              n56w131/",
    "0              n56w132/",
    "0              n56w133/",
    "0              n56w134/",
    "0              n56w135/",
    "0              n56w156/",
    "0              n56w157/",
    "0              n56w159/",
    "0              n56w160/",
    "0              n56w161/",
    "0              n56w162/",
    "0              n56w163/",
    "0              n56w164/",
    "0              n57w131/",
    "0              n57w132/",
    "0              n57w133/",
    "0              n57w134/",
    "0              n57w135/",
    "0              n57w136/",
    "0              n57w153/",
    "0              n57w154/",
    "0              n57w155/",
    "0              n57w156/",
    "0              n57w157/",
    "0              n57w158/",
    "0              n57w159/",
    "0              n57w160/",
    "0              n57w161/",
    "0              n57w162/",
    "0              n57w170/",
    "0              n57w171/",
    "0              n58w133/",
    "0              n58w134/",
    "0              n58w135/",
    "0              n58w136/",
    "0              n58w137/",
    "0              n58w153/",
    "0              n58w154/",
    "0              n58w155/",
    "0              n58w156/",
    "0              n58w157/",
    "0              n58w158/",
    "0              n58w159/",
    "0              n58w170/",
    "0              n58w171/",
    "0              n59w134/",
    "0              n59w135/",
    "0              n59w136/",
    "0              n59w137/",
    "0              n59w138/",
    "0              n59w139/",
    "0              n59w152/",
    "0              n59w153/",
    "0              n59w154/",
    "0              n59w155/",
    "0              n59w156/",
    "0              n59w157/",
    "0              n59w158/",
    "0              n59w159/",
    "0              n59w160/",
    "0              n59w161/",
    "0              n59w162/",
    "0              n59w163/",
    "0              n60w135/",
    "0              n60w136/",
    "0              n60w137/",
    "0              n60w138/",
    "0              n60w139/",
    "0              n60w140/",
    "0              n60w141/",
    "0              n60w142/",
    "0              n60w143/",
    "0              n60w144/",
    "0              n60w145/",
    "0              n60w146/",
    "0              n60w147/",
    "0              n60w148/",
    "0              n60w149/",
    "0              n60w150/",
    "0              n60w151/",
    "0              n60w152/",
    "0              n60w153/",
    "0              n60w154/",
    "0              n60w155/",
    "0              n60w156/",
    "0              n60w157/",
    "0              n60w158/",
    "0              n60w159/",
    "0              n60w160/",
    "0              n60w161/",
    "0              n60w162/",
    "0              n60w163/",
    "0              n60w164/",
    "0              n60w165/",
    "0              n60w166/",
    "0              n60w167/",
    "0              n60w168/",
    "0              n61w140/",
    "0              n61w141/",
    "0              n61w142/",
    "0              n61w143/",
    "0              n61w144/",
    "0              n61w145/",
    "0              n61w146/",
    "0              n61w147/",
    "0              n61w148/",
    "0              n61w149/",
    "0              n61w150/",
    "0              n61w151/",
    "0              n61w152/",
    "0              n61w153/",
    "0              n61w154/",
    "0              n61w155/",
    "0              n61w156/",
    "0              n61w157/",
    "0              n61w158/",
    "0              n61w159/",
    "0              n61w160/",
    "0              n61w161/",
    "0              n61w162/",
    "0              n61w163/",
    "0              n61w164/",
    "0              n61w165/",
    "0              n61w166/",
    "0              n61w167/",
    "0              n61w168/",
    "0              n61w173/",
    "0              n61w174/",
    "0              n62w142/",
    "0              n62w143/",
    "0              n62w144/",
    "0              n62w145/",
    "0              n62w146/",
    "0              n62w147/",
    "0              n62w148/",
    "0              n62w149/",
    "0              n62w150/",
    "0              n62w151/",
    "0              n62w152/",
    "0              n62w153/",
    "0              n62w154/",
    "0              n62w155/",
    "0              n62w156/",
    "0              n62w157/",
    "0              n62w158/",
    "0              n62w159/",
    "0              n62w160/",
    "0              n62w161/",
    "0              n62w162/",
    "0              n62w163/",
    "0              n62w164/",
    "0              n62w165/",
    "0              n62w166/",
    "0              n62w167/",
    "0              n63w142/",
    "0              n63w143/",
    "0              n63w144/",
    "0              n63w145/",
    "0              n63w146/",
    "0              n63w147/",
    "0              n63w148/",
    "0              n63w149/",
    "0              n63w150/",
    "0              n63w151/",
    "0              n63w152/",
    "0              n63w153/",
    "0              n63w154/",
    "0              n63w155/",
    "0              n63w156/",
    "0              n63w157/",
    "0              n63w158/",
    "0              n63w159/",
    "0              n63w160/",
    "0              n63w161/",
    "0              n63w162/",
    "0              n63w163/",
    "0              n63w164/",
    "0              n63w165/",
    "0              n63w166/",
    "0              n63w167/",
    "0              n63w170/",
    "0              n64w141/",
    "0              n64w142/",
    "0              n64w143/",
    "0              n64w144/",
    "0              n64w145/",
    "0              n64w146/",
    "0              n64w147/",
    "0              n64w148/",
    "0              n64w149/",
    "0              n64w150/",
    "0              n64w151/",
    "0              n64w152/",
    "0              n64w153/",
    "0              n64w154/",
    "0              n64w155/",
    "0              n64w156/",
    "0              n64w157/",
    "0              n64w158/",
    "0              n64w159/",
    "0              n64w160/",
    "0              n64w161/",
    "0              n64w162/",
    "0              n64w163/",
    "0              n64w164/",
    "0              n64w165/",
    "0              n64w169/",
    "0              n64w170/",
    "0              n64w171/",
    "0              n64w172/",
    "0              n65w141/",
    "0              n65w142/",
    "0              n65w143/",
    "0              n65w144/",
    "0              n65w145/",
    "0              n65w146/",
    "0              n65w147/",
    "0              n65w148/",
    "0              n65w149/",
    "0              n65w150/",
    "0              n65w151/",
    "0              n65w152/",
    "0              n65w153/",
    "0              n65w154/",
    "0              n65w155/",
    "0              n65w156/",
    "0              n65w157/",
    "0              n65w158/",
    "0              n65w159/",
    "0              n65w160/",
    "0              n65w161/",
    "0              n65w162/",
    "0              n65w163/",
    "0              n65w164/",
    "0              n65w165/",
    "0              n65w166/",
    "0              n65w167/",
    "0              n66w141/",
    "0              n66w142/",
    "0              n66w143/",
    "0              n66w144/",
    "0              n66w145/",
    "0              n66w146/",
    "0              n66w147/",
    "0              n66w148/",
    "0              n66w149/",
    "0              n66w150/",
    "0              n66w151/",
    "0              n66w152/",
    "0              n66w153/",
    "0              n66w154/",
    "0              n66w155/",
    "0              n66w156/",
    "0              n66w157/",
    "0              n66w158/",
    "0              n66w159/",
    "0              n66w160/",
    "0              n66w161/",
    "0              n66w162/",
    "0              n66w163/",
    "0              n66w164/",
    "0              n66w165/",
    "0              n66w166/",
    "0              n66w167/",
    "0              n66w168/",
    "0              n66w169/",
    "0              n67w141/",
    "0              n67w142/",
    "0              n67w143/",
    "0              n67w144/",
    "0              n67w145/",
    "0              n67w146/",
    "0              n67w147/",
    "0              n67w148/",
    "0              n67w149/",
    "0              n67w150/",
    "0              n67w151/",
    "0              n67w152/",
    "0              n67w153/",
    "0              n67w154/",
    "0              n67w155/",
    "0              n67w156/",
    "0              n67w157/",
    "0              n67w158/",
    "0              n67w159/",
    "0              n67w160/",
    "0              n67w161/",
    "0              n67w162/",
    "0              n67w163/",
    "0              n67w164/",
    "0              n67w165/",
    "0              n67w166/",
    "0              n67w167/",
    "0              n67w168/",
    "0              n68w141/",
    "0              n68w142/",
    "0              n68w143/",
    "0              n68w144/",
    "0              n68w145/",
    "0              n68w146/",
    "0              n68w147/",
    "0              n68w148/",
    "0              n68w149/",
    "0              n68w150/",
    "0              n68w151/",
    "0              n68w152/",
    "0              n68w153/",
    "0              n68w154/",
    "0              n68w155/",
    "0              n68w156/",
    "0              n68w157/",
    "0              n68w158/",
    "0              n68w159/",
    "0              n68w160/",
    "0              n68w161/",
    "0              n68w162/",
    "0              n68w163/",
    "0              n68w164/",
    "0              n68w165/",
    "0              n68w166/",
    "0              n69w141/",
    "0              n69w142/",
    "0              n69w143/",
    "0              n69w144/",
    "0              n69w145/",
    "0              n69w146/",
    "0              n69w147/",
    "0              n69w148/",
    "0              n69w149/",
    "0              n69w150/",
    "0              n69w151/",
    "0              n69w152/",
    "0              n69w153/",
    "0              n69w154/",
    "0              n69w155/",
    "0              n69w156/",
    "0              n69w157/",
    "0              n69w158/",
    "0              n69w159/",
    "0              n69w160/",
    "0              n69w161/",
    "0              n69w162/",
    "0              n69w163/",
    "0              n69w164/",
    "0              n69w165/",
    "0              n69w166/",
    "0              n69w167/",
    "0              n70w141/",
    "0              n70w142/",
    "0              n70w143/",
    "0              n70w144/",
    "0              n70w145/",
    "0              n70w146/",
    "0              n70w147/",
    "0              n70w148/",
    "0              n70w149/",
    "0              n70w150/",
    "0              n70w151/",
    "0              n70w152/",
    "0              n70w153/",
    "0              n70w154/",
    "0              n70w155/",
    "0              n70w156/",
    "0              n70w157/",
    "0              n70w158/",
    "0              n70w159/",
    "0              n70w160/",
    "0              n70w161/",
    "0              n70w162/",
    "0              n70w163/",
    "0              n70w164/",
    "0              n70w165/",
    "0              n71w143/",
    "0              n71w144/",
    "0              n71w145/",
    "0              n71w146/",
    "0              n71w147/",
    "0              n71w148/",
    "0              n71w149/",
    "0              n71w150/",
    "0              n71w151/",
    "0              n71w152/",
    "0              n71w153/",
    "0              n71w154/",
    "0              n71w155/",
    "0              n71w156/",
    "0              n71w157/",
    "0              n71w158/",
    "0              n71w159/",
    "0              n71w160/",
    "0              n71w161/",
    "0              n71w162/",
    "0              n71w163/",
    "0              n71w164/",
    "0              n72w155/",
    "0              n72w156/",
    "0              n72w157/",
    "0              n72w158/",
    "0              s14w170/",
    "0              s14w171/"
  )
  
  
  # sub out HTML copy pattern with vsi URL
  t2 <- gsub("0              ",
             "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13/TIFF/current/",
             t)
  
  # add file paths following NED scheme
  t3 <- paste0(t2,
               "USGS_13_",
               basename(t2),
               ".tif")
  
  # Write table to data-raw
  write.table(t3,
              paste0(base_dir, "/ned_list_USGS_13.txt"),
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
  
  
}
# dem_base_dir <- "/Users/anguswatters/Desktop/transects_paper/data/dem"
# download_3dep_vrt(dem_base_dir)
# 
# # Create meta data object of three NED resoruces
# ned <- data.frame(rbind(
#   # c(id        = "USGS_1",
#   #   URL       = "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1",
#   #   varname   = "30m elevation",
#   #   long_name = "30m (1 arcsec) National Elevation Dataset",
#   #   units     = "m"),
#   # 
#   # c(id        = "USGS_2",
#   #   URL       = "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/2",
#   #   varname   = "60m elevation",
#   #   long_name = "60m (2 arcsec) National Elevation Dataset Alaska",
#   #   units     = "m"),
#   
#   c(id        = "USGS_13",
#     URL       = "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13",
#     varname   = "10m elevation",
#     long_name = "10m (1/3th arcsec) National Elevation Dataset",
#     units     = "m")
# ))
# 
# 
# 
# 
# # Loop over the three resolutions
# for (i in 1:length(ned)) {
#   i = 1
#   
#   # Define output text file path
#   txt_file  <- paste0(dem_base_dir, "/ned_list_", ned$id[i], "_2.txt")
#   
#   # Define output VRT path
#   vrt_file <- paste0(paste0(dem_base_dir, "/ned_", ned$id[i], ".vrt"))
#   
#   # If VRT does NOT exist, build VRT
#   if (!file.exists(vrt_file)) {
#     
#     # read the corresponding index.gpkg
#     files <- sf::read_sf(ned$domain_url[i])
#     DEM_UR
#     DEM_PATH        <- "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/TIFF/USGS_Seamless_DEM_1.vrt"
#     # Build full HTTPS paths to "./current/"
#     files <- c(file.path(ned$URL[i], "TIFF/current", gsub("[.]/", "", files$location)))
#     
#     files <- "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/TIFF/USGS_Seamless_DEM_1.vrt"
#     
#     # write list of files to text file
#     write.table(files, txt_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
#     
#     # build VRT from text file input using GDAL system call ...
#     system(paste("gdalbuildvrt -input_file_list", txt_file, vrt_file))
#   }
#   
#   logger::log_info("Finished ", ned$id[i], "...")
# }

# 
# match_transects_to_extended_cs_pts <- function(transect_lines, fixed_cs_pts, crosswalk_id, extension_pct = 0.5 ) {
#   
#   # transect_lines = transects
#   # fixed_cs_pts   = fixed_pts
#   # crosswalk_id   = CROSSWALK_ID
#   
#   fixed_cs_pts <- nhdplusTools::rename_geometry(fixed_cs_pts, "geometry")
#   transect_lines    <- nhdplusTools::rename_geometry(transect_lines, "geometry")
#   
#   # get the counts of each point type to add this data to the transect_lines dataset
#   point_type_counts <- hydrofabric3D::get_point_type_counts(classified_pts = fixed_cs_pts, 
#                                                             crosswalk_id = crosswalk_id)
#   # Check the number of cross sections that were extended
#   message("Subsetting cross section points generated after extending transect_lines...")
#   
#   # extract cross section points that have an "is_extended" value of TRUE
#   extended_pts <- 
#     fixed_cs_pts %>%
#     dplyr::filter(is_extended) %>%
#     hydrofabric3D::add_tmp_id(x = crosswalk_id)
#   
#   # extended_pts %>% 
#   #   get_unique_tmp_ids() %>% 
#   #   length()
#   
#   # extract transect_lines that have a "crosswalk_id" in the "extended_pts" dataset
#   update_transect_lines <- 
#     transect_lines %>%
#     hydrofabric3D::add_tmp_id(x = crosswalk_id) %>%
#     dplyr::filter(tmp_id %in% unique(extended_pts$tmp_id))
#   
#   cs_pt_uids    <- unique(hydrofabric3D::add_tmp_id(fixed_cs_pts, x = crosswalk_id)$tmp_id)
#   
#   # If any transect_lines were extended, update the transect_lines dataset, and overwrite local and S3 transect_lines geopackages
#   if (nrow(update_transect_lines) > 0) {
#     message("Updating ", nrow(update_transect_lines), " transect_lines")
#     
#     
#     # update_transect_lines <- 
#     #   update_transect_lines %>% 
#     #   dplyr::rename(hy_id := !!sym(crosswalk_id)) 
#     # 
#     update_transect_lines <- 
#       update_transect_lines %>%
#       # apply extend_by_percent function to each transect line:
#       hydrofabric3D:::extend_by_percent(
#         crosswalk_id = crosswalk_id,
#         pct        = extension_pct,
#         length_col = "cs_lengthm"
#       )
#     
#     update_transect_lines <- hydroloom::rename_geometry(update_transect_lines, "geometry")
#     
#     # update_transect_lines <- 
#     #   update_transect_lines %>%  
#     #   dplyr::rename(!!sym(crosswalk_id) := hy_id)
#     
#     # cs_pt_uids    <- unique(hydrofabric3D::add_tmp_id(fixed_cs_pts, x = get(crosswalk_id))$tmp_id)
#     # transect_uids <- unique(hydrofabric3D::add_tmp_id(transect_lines, x = get(crosswalk_id))$tmp_id)
#     
#     # Filter down to ONLY points that were finalized and rectified from rectify_cs_pts()
#     # Remove old transect_lines that have "tmp_id" in "extended_pts" (transect_lines that were unchanged and are "good_to_go")
#     # and then replace with old transect_lines with the "update_transect_lines"
#     out_transect_lines <-
#       transect_lines %>%
#       hydrofabric3D::add_tmp_id(x = crosswalk_id) %>%
#       dplyr::filter(tmp_id %in% cs_pt_uids) %>% 
#       dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id)) %>%
#       dplyr::bind_rows(
#         dplyr::mutate(update_transect_lines, is_extended = TRUE)
#       )
#     
#     # transect_lines %>% 
#     #   hydrofabric3D::add_tmp_id(x = "hy_id") %>%
#     #   # dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id)) %>%
#     #   dplyr::filter(tmp_id %in% unique(hydrofabric3D::add_tmp_id(fixed_pts, x = "hy_id")$tmp_id)) %>% # Subset down to the remaining tmp_ids in the fixed points
#     #   dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id)) %>% # remove the tmp_ids that we are going add back in with the extended versions of those tmp_ids
#     #   dplyr::bind_rows( # bring in the new updated extended transect_lines
#     #     dplyr::mutate(
#     #       update_transect_lines,
#     #       is_extended = TRUE
#     #     )
#     #   )  
#   } else {
#     # If no transect_lines were extended
#     out_transect_lines <- 
#       transect_lines %>%
#       hydrofabric3D::add_tmp_id(x = crosswalk_id) %>%
#       dplyr::filter(tmp_id %in% cs_pt_uids) %>% 
#       # dplyr::filter(tmp_id %in% unique(hydrofabric3D::add_tmp_id(fixed_cs_pts, x = get(crosswalk_id))$tmp_id)) %>%
#       dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id))
#   }
#   
#   # Finalize new transect_lines
#   out_transect_lines <- 
#     out_transect_lines %>%
#     dplyr::left_join(
#       point_type_counts, 
#       by = c(crosswalk_id, "cs_id")
#     ) %>%
#     dplyr::left_join(
#       dplyr::ungroup(
#         dplyr::slice(
#           dplyr::group_by(
#             dplyr::select(sf::st_drop_geometry(fixed_cs_pts),
#                           dplyr::any_of(crosswalk_id), 
#                           cs_id, bottom, left_bank, right_bank, valid_banks, has_relief
#             ),
#             dplyr::across(dplyr::any_of(c(crosswalk_id, "cs_id")))
#           ),
#           1
#         )
#       ),
#       by = c(crosswalk_id, "cs_id")
#     ) %>%
#     dplyr::select(
#       dplyr::any_of(crosswalk_id),
#       cs_source, cs_id, cs_measure, cs_lengthm,
#       # sinuosity,
#       is_extended,
#       left_bank_count, right_bank_count, channel_count, bottom_count,
#       bottom, left_bank, right_bank, valid_banks, has_relief,
#       geometry
#     ) %>% 
#     dplyr::mutate(
#       is_extended = ifelse(is.na(is_extended), FALSE, is_extended)
#     )  
#   
#   return(out_transect_lines)
# }

# utility function for getting transects extended and 
# matching cross section points that went through "get_improved_cs_pts()" and that were extended for improvement
# returns the extended version of the transects 
# match_transects_to_extended_cs_pts <- function(transect_lines, fixed_cs_pts, crosswalk_id) {
#   
#   # transect_lines = transects
#   # fixed_cs_pts   = fixed_pts
#   # crosswalk_id   = CROSSWALK_ID
#   
#   fixed_cs_pts <- nhdplusTools::rename_geometry(fixed_cs_pts, "geometry")
#   transect_lines    <- nhdplusTools::rename_geometry(transect_lines, "geometry")
#   
#   # get the counts of each point type to add this data to the transect_lines dataset
#   point_type_counts <- hydrofabric3D::get_point_type_counts(classified_pts = fixed_cs_pts, 
#                                                             crosswalk_id = crosswalk_id)
#   # Check the number of cross sections that were extended
#   message("Subsetting cross section points generated after extending transect_lines...")
#   
#   # extract cross section points that have an "is_extended" value of TRUE
#   extended_pts <- 
#     fixed_cs_pts %>%
#     dplyr::filter(is_extended) %>%
#     hydrofabric3D::add_tmp_id(x = crosswalk_id)
#   
#   # extended_pts %>% 
#   #   get_unique_tmp_ids() %>% 
#   #   length()
#   
#   # extract transect_lines that have a "crosswalk_id" in the "extended_pts" dataset
#   update_transect_lines <- 
#     transect_lines %>%
#     hydrofabric3D::add_tmp_id(x = crosswalk_id) %>%
#     dplyr::filter(tmp_id %in% unique(extended_pts$tmp_id))
#   
#   cs_pt_uids    <- unique(hydrofabric3D::add_tmp_id(fixed_cs_pts, x = crosswalk_id)$tmp_id)
#   
#   # If any transect_lines were extended, update the transect_lines dataset, and overwrite local and S3 transect_lines geopackages
#   if (nrow(update_transect_lines) > 0) {
#     message("Updating ", nrow(update_transect_lines), " transect_lines")
#     
#     
#     update_transect_lines <- 
#       update_transect_lines %>% 
#       dplyr::rename(hy_id := !!sym(crosswalk_id)) 
#     
#     update_transect_lines <- 
#       update_transect_lines %>%
#       # apply extend_by_percent function to each transect line:
#       hydrofabric3D:::extend_by_percent(
#         pct        = EXTENSION_PCT,
#         length_col = "cs_lengthm"
#       )
#     
#     update_transect_lines <- hydroloom::rename_geometry(update_transect_lines, "geometry")
#     
#     update_transect_lines <- 
#       update_transect_lines %>%  
#       dplyr::rename(!!sym(crosswalk_id) := hy_id)
# 
#     # cs_pt_uids    <- unique(hydrofabric3D::add_tmp_id(fixed_cs_pts, x = get(crosswalk_id))$tmp_id)
#     # transect_uids <- unique(hydrofabric3D::add_tmp_id(transect_lines, x = get(crosswalk_id))$tmp_id)
#     
#     # Filter down to ONLY points that were finalized and rectified from rectify_cs_pts()
#     # Remove old transect_lines that have "tmp_id" in "extended_pts" (transect_lines that were unchanged and are "good_to_go")
#     # and then replace with old transect_lines with the "update_transect_lines"
#     out_transect_lines <-
#       transect_lines %>%
#       hydrofabric3D::add_tmp_id(x = crosswalk_id) %>%
#       dplyr::filter(tmp_id %in% cs_pt_uids) %>% 
#       dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id)) %>%
#       dplyr::bind_rows(
#         dplyr::mutate(update_transect_lines, is_extended = TRUE)
#       )
#     
#     # transect_lines %>% 
#     #   hydrofabric3D::add_tmp_id(x = "hy_id") %>%
#     #   # dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id)) %>%
#     #   dplyr::filter(tmp_id %in% unique(hydrofabric3D::add_tmp_id(fixed_pts, x = "hy_id")$tmp_id)) %>% # Subset down to the remaining tmp_ids in the fixed points
#     #   dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id)) %>% # remove the tmp_ids that we are going add back in with the extended versions of those tmp_ids
#     #   dplyr::bind_rows( # bring in the new updated extended transect_lines
#     #     dplyr::mutate(
#     #       update_transect_lines,
#     #       is_extended = TRUE
#     #     )
#     #   )  
#   } else {
#     # If no transect_lines were extended
#     out_transect_lines <- 
#       transect_lines %>%
#       hydrofabric3D::add_tmp_id(x = crosswalk_id) %>%
#       dplyr::filter(tmp_id %in% cs_pt_uids) %>% 
#       # dplyr::filter(tmp_id %in% unique(hydrofabric3D::add_tmp_id(fixed_cs_pts, x = get(crosswalk_id))$tmp_id)) %>%
#       dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id))
#   }
#   
#   # Finalize new transect_lines
#   out_transect_lines <- 
#     out_transect_lines %>%
#     dplyr::left_join(
#       point_type_counts, 
#       by = c(crosswalk_id, "cs_id")
#     ) %>%
#     dplyr::left_join(
#       dplyr::ungroup(
#         dplyr::slice(
#           dplyr::group_by(
#             dplyr::select(sf::st_drop_geometry(fixed_cs_pts),
#                           dplyr::any_of(crosswalk_id), 
#                           cs_id, bottom, left_bank, right_bank, valid_banks, has_relief
#             ),
#             dplyr::across(dplyr::any_of(c(crosswalk_id, "cs_id")))
#           ),
#           1
#         )
#       ),
#       by = c(crosswalk_id, "cs_id")
#     ) %>%
#     dplyr::select(
#       dplyr::any_of(crosswalk_id),
#       cs_source, cs_id, cs_measure, cs_lengthm,
#       # sinuosity,
#       is_extended,
#       left_bank_count, right_bank_count, channel_count, bottom_count,
#       bottom, left_bank, right_bank, valid_banks, has_relief,
#       geometry
#     ) %>% 
#     dplyr::mutate(
#       is_extended = ifelse(is.na(is_extended), FALSE, is_extended)
#     )  
#   
#   return(out_transect_lines)
# }



