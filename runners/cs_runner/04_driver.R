### Run this file to have all runner scripts run in order

# downloads nextgen datasets 
source("runners/cs_runner/config.R")

# downloads datasets
# - Nextgen data
# - Reference features (for waterbody filtering)
# - ML outputs
# - FEMA 100 year floodplain polygons (FGBs)
source("runners/cs_runner/download_nextgen.R")
source("runners/cs_runner/download_ref_features.R")
source("runners/cs_runner/download_ml_outputs.R")
source("runners/cs_runner/download_fema100.R")

# simplify, dissolve, FEMA polygons and partition FEMA polygons by VPU 
source("runners/cs_runner/partition_fema_by_vpu.R")

# generate and upload transects datasets 
source("runners/cs_runner/01_transects.R")

# generate and upload cross sections points datasets 
source("runners/cs_runner/02_cs_pts.R")

# Apply machine learning topwidths and depths estimates to DEM cross section points
source("runners/cs_runner/03_inject_ml.R")
