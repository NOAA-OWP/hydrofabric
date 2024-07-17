### Run this file to have all runner scripts run in order

# downloads nextgen datasets 
source("runners/cs_runner/config.R")

# downloads nextgen datasets 
source("runners/cs_runner/download_nextgen.R")

# download FEMA100 year FGBs
source("runners/cs_runner/download_fema100.R")

# simplify, dissolve, FEMA polygons and partition FEMA polygons by VPU 
source("runners/cs_runner/partition_fema_by_vpu.R")

# generate and upload transects datasets 
source("runners/cs_runner/01_transects.R")

# generate and upload cross sections points datasets 
source("runners/cs_runner/02_cs_pts.R")

# Apply machine learning topwidths and depths estimates to DEM cross section points
source("runners/cs_runner/02_cs_pts.R")
