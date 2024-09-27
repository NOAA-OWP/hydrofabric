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
  
  create_if_not_exists <- function(dir_path) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
      message("Created directory: '", dir_path, "'\n")
    }
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


# utility function for getting transects extended and 
# matching cross section points that went through "get_improved_cs_pts()" and that were extended for improvement
# returns the extended version of the transects 
match_transects_to_extended_cs_pts <- function(transect_lines, fixed_cs_pts, crosswalk_id) {
  
  # transect_lines = transects
  # fixed_cs_pts   = fixed_pts
  # crosswalk_id   = CROSSWALK_ID
  
  fixed_cs_pts <- nhdplusTools::rename_geometry(fixed_cs_pts, "geometry")
  transect_lines    <- nhdplusTools::rename_geometry(transect_lines, "geometry")
  
  # get the counts of each point type to add this data to the transect_lines dataset
  point_type_counts <- hydrofabric3D::get_point_type_counts(classified_pts = fixed_cs_pts, 
                                                            crosswalk_id = crosswalk_id)
  # Check the number of cross sections that were extended
  message("Subsetting cross section points generated after extending transect_lines...")
  
  # extract cross section points that have an "is_extended" value of TRUE
  extended_pts <- 
    fixed_cs_pts %>%
    dplyr::filter(is_extended) %>%
    hydrofabric3D::add_tmp_id(x = crosswalk_id)
  
  # extended_pts %>% 
  #   get_unique_tmp_ids() %>% 
  #   length()
  
  # extract transect_lines that have a "crosswalk_id" in the "extended_pts" dataset
  update_transect_lines <- 
    transect_lines %>%
    hydrofabric3D::add_tmp_id(x = crosswalk_id) %>%
    dplyr::filter(tmp_id %in% unique(extended_pts$tmp_id))
  
  cs_pt_uids    <- unique(hydrofabric3D::add_tmp_id(fixed_cs_pts, x = crosswalk_id)$tmp_id)
  
  # If any transect_lines were extended, update the transect_lines dataset, and overwrite local and S3 transect_lines geopackages
  if (nrow(update_transect_lines) > 0) {
    message("Updating ", nrow(update_transect_lines), " transect_lines")
    
    
    update_transect_lines <- 
      update_transect_lines %>% 
      dplyr::rename(hy_id := !!sym(crosswalk_id)) 
    
    update_transect_lines <- 
      update_transect_lines %>%
      # apply extend_by_percent function to each transect line:
      hydrofabric3D:::extend_by_percent(
        pct        = EXTENSION_PCT,
        length_col = "cs_lengthm"
      )
    
    update_transect_lines <- hydroloom::rename_geometry(update_transect_lines, "geometry")
    
    update_transect_lines <- 
      update_transect_lines %>%  
      dplyr::rename(!!sym(crosswalk_id) := hy_id)

    # cs_pt_uids    <- unique(hydrofabric3D::add_tmp_id(fixed_cs_pts, x = get(crosswalk_id))$tmp_id)
    # transect_uids <- unique(hydrofabric3D::add_tmp_id(transect_lines, x = get(crosswalk_id))$tmp_id)
    
    # Filter down to ONLY points that were finalized and rectified from rectify_cs_pts()
    # Remove old transect_lines that have "tmp_id" in "extended_pts" (transect_lines that were unchanged and are "good_to_go")
    # and then replace with old transect_lines with the "update_transect_lines"
    out_transect_lines <-
      transect_lines %>%
      hydrofabric3D::add_tmp_id(x = crosswalk_id) %>%
      dplyr::filter(tmp_id %in% cs_pt_uids) %>% 
      dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id)) %>%
      dplyr::bind_rows(
        dplyr::mutate(update_transect_lines, is_extended = TRUE)
      )
    
    # transect_lines %>% 
    #   hydrofabric3D::add_tmp_id(x = "hy_id") %>%
    #   # dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id)) %>%
    #   dplyr::filter(tmp_id %in% unique(hydrofabric3D::add_tmp_id(fixed_pts, x = "hy_id")$tmp_id)) %>% # Subset down to the remaining tmp_ids in the fixed points
    #   dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id)) %>% # remove the tmp_ids that we are going add back in with the extended versions of those tmp_ids
    #   dplyr::bind_rows( # bring in the new updated extended transect_lines
    #     dplyr::mutate(
    #       update_transect_lines,
    #       is_extended = TRUE
    #     )
    #   )  
  } else {
    # If no transect_lines were extended
    out_transect_lines <- 
      transect_lines %>%
      hydrofabric3D::add_tmp_id(x = crosswalk_id) %>%
      dplyr::filter(tmp_id %in% cs_pt_uids) %>% 
      # dplyr::filter(tmp_id %in% unique(hydrofabric3D::add_tmp_id(fixed_cs_pts, x = get(crosswalk_id))$tmp_id)) %>%
      dplyr::filter(!tmp_id %in% unique(extended_pts$tmp_id))
  }
  
  # Finalize new transect_lines
  out_transect_lines <- 
    out_transect_lines %>%
    dplyr::left_join(
      point_type_counts, 
      by = c(crosswalk_id, "cs_id")
    ) %>%
    dplyr::left_join(
      dplyr::ungroup(
        dplyr::slice(
          dplyr::group_by(
            dplyr::select(sf::st_drop_geometry(fixed_cs_pts),
                          dplyr::any_of(crosswalk_id), 
                          cs_id, bottom, left_bank, right_bank, valid_banks, has_relief
            ),
            dplyr::across(dplyr::any_of(c(crosswalk_id, "cs_id")))
          ),
          1
        )
      ),
      by = c(crosswalk_id, "cs_id")
    ) %>%
    dplyr::select(
      dplyr::any_of(crosswalk_id),
      cs_source, cs_id, cs_measure, cs_lengthm,
      # sinuosity,
      is_extended,
      left_bank_count, right_bank_count, channel_count, bottom_count,
      bottom, left_bank, right_bank, valid_banks, has_relief,
      geometry
    ) %>% 
    dplyr::mutate(
      is_extended = ifelse(is.na(is_extended), FALSE, is_extended)
    )  
  
  return(out_transect_lines)
}



