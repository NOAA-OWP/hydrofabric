# Running this script goes and pulls the desired NextGen geopackage datasets from http://www.lynker-spatial.com/, saves them into a directory within "BASE_DIR"
# BASE_DIR is defined within runners/workflow/root_dir.R

# load config variables
source("runners/cs_runner2/base_variables.R")
# source("runners/cs_runner/config_vars.R")

# ---------------------------------------------------------------------------
# ---- Download conus_nextgen.gpkg 
# ---------------------------------------------------------------------------

copy_cmd <- paste0('aws s3 cp ', CONUS_ML_S3_URI, " ", CONUS_ML_PARQUET_PATH)
message("Copying S3 object:\n", CONUS_ML_S3_URI)

if (!file.exists(CONUS_ML_PARQUET_PATH)) {
  tryCatch({
    system(copy_cmd)
    message("Download '", basename(CONUS_ML_PARQUET_PATH), "' complete!")
    message("------------------")
  }, error = function(e) {
    message("Error downloading conus_nextgen.gpkg")
    message(e)
    stop()
  })
  
} else {
  message("'", basename(CONUS_ML_PARQUET_PATH), "' file already exists at\n > '", CONUS_ML_PARQUET_PATH, "'")
}
