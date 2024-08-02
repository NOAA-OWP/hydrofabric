# download_ml_outputs.R
# This script pulls the ML outputs data from http://www.lynker-spatial.com/, saves them into a directory within "BASE_DIR"

# load config variables
source("runners/cs_runner/config_vars.R")

# ---------------------------------------------------------------------------
# ---- Get ML outputs data from S3 bucket ----
# ---------------------------------------------------------------------------

ml_copy_cmd <- paste0('aws s3 cp ', ML_OUTPUTS_S3_URI, ' ', paste0(ML_OUTPUTS_DIR, basename(ML_OUTPUTS_S3_URI)))

message("Copying S3 object:\n", ML_OUTPUTS_S3_URI)
system(ml_copy_cmd)
message("Download '", paste0(ML_OUTPUTS_DIR, basename(ML_OUTPUTS_S3_URI)), "' complete!")
message("------------------")
