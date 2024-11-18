### Run this file to have all runner scripts run in order

# downloads nextgen datasets 
source("runners/cs_runner2/config_env.R")

# downloads datasets
# - Nextgen data
# - Reference features (for waterbody filtering)
# - ML outputs
# - FEMA 100 year floodplain polygons (FGBs)
# - 3DEP DEM VRT
source("runners/cs_runner2/download_conus_nextgen.R")
source("runners/cs_runner2/download_conus_ref_features.R")
source("runners/cs_runner2/download_fema.r")
source("runners/cs_runner2/process_fema.R")
source("runners/cs_runner2/download_dem_from_vrt.R")
source("runners/cs_runner2/download_ml_outputs.R")

# generate and upload transects datasets 
source("runners/cs_runner2/01_transects.R")

# generate and upload cross sections points datasets 
source("runners/cs_runner2/02_cs_pts.R")

# Apply machine learning topwidths and depths estimates to DEM cross section points
source("runners/cs_runner2/03_inject_ml.R")
