# Running this script goes and pulls the desired NextGen geopackage datasets from http://www.lynker-spatial.com/, saves them into a directory within "BASE_DIR"
# BASE_DIR is defined within runners/workflow/root_dir.R

# load config variables
source("runners/cs_runner2/base_variables.R")
# source("runners/cs_runner/config_vars.R")

# ---------------------------------------------------------------------------
# ---- Download conus_nextgen.gpkg 
# ---------------------------------------------------------------------------

copy_cmd <- paste0('aws s3 cp ', CONUS_NEXTGEN_S3_URI, " ", CONUS_NEXTGEN_GPKG_PATH)
message("Copying S3 object:\n", CONUS_NEXTGEN_S3_URI)

if (!file.exists(CONUS_NEXTGEN_GPKG_PATH)) {
  tryCatch({
      system(copy_cmd)
      message("Download '", basename(CONUS_NEXTGEN_GPKG_PATH), "' complete!")
      message("------------------")
  }, error = function(e) {
      message("Error downloading conus_nextgen.gpkg")
      message(e)
      stop()
  })

} else {
  message("conus_nextgen.gpkg file already exists at\n > '", CONUS_NEXTGEN_GPKG_PATH, "'")
}
