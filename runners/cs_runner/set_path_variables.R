# Generate the list of file paths for locally stored nextgen datasets:
# - NEXTGEN_FILES --> (list NEXTGEN_DIR files that were downloaded via download_nextgen.R)


# load config variables
source("runners/cs_runner/config_vars.R")

# NOTE: SET VARIABLE FOR REST OF PROCESSING 
# paths to nextgen datasets
NEXTGEN_FILES   <- list.files(NEXTGEN_DIR, full.names = FALSE)
