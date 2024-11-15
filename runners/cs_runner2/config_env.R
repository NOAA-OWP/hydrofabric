# load required packages
pacman::p_load(
  archive,
  hydrofabric,
  hydrofabric3D,
  dplyr,
  sf
)

# # install.packages("devtools")
# devtools::install_github("anguswg-ucsb/hydrofabric3D")

# load root directory 
source("runners/cs_runner2/base_variables.R")
source("runners/cs_runner2/utils.R")

sf::sf_use_s2(FALSE)

# create empty base directories 
create_local_hydrofabric_base_dirs(base_dir = BASE_DIR)

# create a new version directory
create_new_version_dirs(base_dir = BASE_DIR, 
                        version = VERSION, 
                        with_output = TRUE)

