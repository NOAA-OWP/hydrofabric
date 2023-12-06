# load required packages
pacman::p_load(
  archive,
  hydrofabric,
  hydrofabric3D
  # terrainSliceR
)

# load root directory 
source("runners/cs_runner/config_vars.R")

sf::sf_use_s2(FALSE)

### Cross section point 

### S3 names

# name of S3 bucket
s3_bucket <- "s3://lynker-spatial/"

# name of bucket with nextgen data
nextgen_bucket <- "lynker-spatial"

# reference features S3 bucket prefix
ref_features_prefix <- "s3://lynker-spatial/00_reference_features/gpkg/"

# S3 prefix/folder of version run
version_prefix <- "v20.1"
# version_prefix <- "v20"

### LOCAL DIRS
# directory to copy nextgen bucket data too
nextgen_dir <- paste0(base_dir, "/pre-release/")

# model attributes directory
model_attr_dir <- paste0(base_dir, "/model_attributes/")

# cross-section data model data directories
transects_dir <- paste0(base_dir, "/01_transects/")
cs_pts_dir <- paste0(base_dir, "/02_cs_pts/")

# final output directory with geopackages per VPU
final_dir <- paste0(base_dir, "/cross_sections/")

# directory to copy nextgen bucket data too
ref_features_dir <- paste0(base_dir, "/00_reference_features/")

# create directories 
dir.create(transects_dir, showWarnings = FALSE)
dir.create(cs_pts_dir,    showWarnings = FALSE)
dir.create(ref_features_dir,     showWarnings = FALSE)
dir.create(paste0(ref_features_dir, "gpkg/"),     showWarnings = FALSE)
dir.create(final_dir,     showWarnings = FALSE)
# dir.create(model_attr_dir,  showWarnings = FALSE)

## Go get a list of the reference features geopackages from S3 and create a save path using the S3 file names to save reference features to local directory

# list objects in S3 bucket, and regular expression match to nextgen_.gpkg pattern
list_ref_features <- paste0('#!/bin/bash
            # AWS S3 Bucket and Directory information
            S3_BUCKET="', ref_features_prefix , '"
            
            # Regular expression pattern to match object keys
            PATTERN="reference_features.gpkg"

            S3_OBJECTS=$(aws s3 ls "$S3_BUCKET" | awk \'{print $4}\' | grep -E "$PATTERN")
            
            echo "$S3_OBJECTS"'
)

# ---- Get a list of reference features geopackages ----

# Run the script to get a list of the nextgen geopackages that matched the regular expression above
ref_features <- system(list_ref_features, intern = TRUE)

# ref features datasets
ref_features_keys <- paste0(ref_features_prefix, ref_features)
ref_features_files <- paste0(ref_features_dir, "gpkg/", ref_features)

###
### UTILITY FUNCTION FOR MATCHING FILES BASED ON VPU STRING ###
###

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
  ########  ########  ########  ########  ########  ########
  # if(type == 1) {
    
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
    
    # trans[trans$hy_id %in% unique(to_check$id), ] %>% nrow()
    # 
    # trans_geos[trans$hy_id %in% unique(to_check$id)]
    
    # check where the transects linestrings intersect with the waterbodies
    trans_geos_check <- geos::as_geos_geometry(trans_check)
    
    message("Checking transects against waterbodies...")
    trans_wb_index <- geos::geos_intersects_any(
      trans_geos_check,  
      wbs_geos[unlist(wb_index)]
    )
  
    # within the transects lines that are on a flowline that crosses a waterbody, 
    # check if any of these transects line DO NOT CROSS A WATERBODY AT ALL
    trans_keep <- trans_check[!trans_wb_index, ]
    # trans_keep <- trans_check[lengths(trans_wb_index2) == 0, ]
  
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
  # }
  ########  ########  ########  ########  ########  ########
  ########  ########  ########  ########  ########  ########
  # if(type == 2) {
  # 
  #   # temporary ID for transects that is the "hy_id", underscore, "cs_id", used for subsetting in future steps
  #   trans$tmp_id <- paste0(trans$hy_id, "_", trans$cs_id)
  #   
  #   # trans$order <- 1:nrow(trans)
  #   
  #   message("Checking flowlines against waterbodies...")
  #   
  #   # create an index between flowlines and waterbodies 
  #   wb_index <- sf::st_intersects(flowlines, waterbodies)
  #   
  #   # remove any flowlines that cross more than 1 waterbody
  #   to_keep <- flowlines[lengths(wb_index) == 0, ]
  #   to_check <- flowlines[lengths(wb_index) != 0, ]
  #   
  #   # subset transects to the hy_ids in "to_check" set of flowlines
  #   trans_check <- trans[trans$hy_id %in% unique(to_check$id), ]
  #   
  #   message("Checking transects against waterbodies...")
  #   
  #   # check where the transects linestrings intersect with the waterbodies
  #   trans_wb_index <- sf::st_intersects(trans_check, 
  #                                       waterbodies[unlist(wb_index), ]
  #   )
  #   # trans_wb_index <- sf::st_intersects(trans_check, wbs)
  #   
  #   # within the transects lines that are on a flowline that crosses a waterbody, 
  #   # check if any of these transects line DO NOT CROSS A WATERBODY AT ALL
  #   trans_keep <- trans_check[lengths(trans_wb_index) == 0, ]
  #   
  #   # preserve any flowlines that CROSS A WATERBODY BUT ALSO HAVE A TRANSECT LINE that does NOT cross any waterbodies
  #   to_check <- to_check[to_check$id %in% unique(trans_keep$hy_id), ]
  #   
  #   # update flowlines to keep with flowlines that intersect a waterbody BUT STILL,
  #   # have transects that are NOT in the waterbody
  #   to_keep <- dplyr::bind_rows(to_keep, to_check)
  #   
  #   # 'tmp_ids' of transects that are being checked and also the transects within trans_check 
  #   # that were determined to be valid (are being kept)
  #   check_ids <- unique(trans_check$tmp_id)
  #   keep_ids <- unique(trans_keep$tmp_id)
  #   
  #   # logical vectors of which flowlines/transects to keep (KEEP == TRUE)
  #   # - Remove any transects that are on flowlines that cross a waterbody AND the transect crosses the waterbody too.
  #   # - Keep original transects that are not on flowlines that intersect waterbodies AND 
  #   # also the transects that do NOT intersect waterbodies but are on a flowline that DOES intersect a waterbody
  #   valid_flowlines <- flowlines$id %in% to_keep$id
  #   valid_transects <- trans$tmp_id %in% dplyr::filter(trans,
  #                                                      !tmp_id %in% check_ids[!check_ids %in% keep_ids])$tmp_id
  #   
  #   # return alist of updated flowlines and transects 
  #   return(
  #     list(
  #       "valid_flowlines" =  valid_flowlines,
  #       "valid_transects" =  valid_transects
  #     )
  #   )
  #   
  #   # # - Remove any transects that are on flowlines that cross a waterbody AND the transect crosses the waterbody too.
  #   # # - Keep original transects that are not on flowlines that intersect waterbodies AND 
  #   #     # also the transects that do NOT intersect waterbodies but are on a flowline that DOES intersect a waterbody
  #   # trans <-
  #   #   trans %>%
  #   #   dplyr::filter(
  #   #     !tmp_id %in% check_ids[!check_ids %in% keep_ids]
  #   #     ) %>%
  #   #   dplyr::select(-tmp_id)
  #   # 
  #   # # return alist of updated flowlines and transects 
  #   # return(
  #   #   list(
  #   #     "flowlines" = to_keep,
  #   #     "transects" = trans
  #   #     )
  #   # )
  # }
}
