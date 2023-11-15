# load required packages
pacman::p_load(
  archive,
  terrainSliceR
)

# load root directory 
source("runners/cs_runner/config_vars.R")

sf::sf_use_s2(FALSE)

# name of bucket with nextgen data
nextgen_bucket <- "lynker-spatial"

# directory to copy nextgen bucket data too
nextgen_dir <- paste0(base_dir, "/pre-release/")

# model attributes directory
model_attr_dir <- paste0(base_dir, "/model_attributes/")

# cross-section data model data directories
transects_dir <- paste0(base_dir, "/01_transects/")
cs_pts_dir <- paste0(base_dir, "/02_cs_pts/")

# final output directory with geopackages per VPU
final_dir <- paste0(base_dir, "/cross_sections/")

# create directories 
dir.create(transects_dir, showWarnings = FALSE)
dir.create(cs_pts_dir,    showWarnings = FALSE)
dir.create(final_dir,     showWarnings = FALSE)
# dir.create(model_attr_dir,  showWarnings = FALSE)


##### UTILITY FUNCTION FOR MATCHING FILES BASED ON VPU STRING ######
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
