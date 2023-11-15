### Run this file to have all runner scripts run in order

# downloads nextgen datasets 
source("runners/cs_runner/config.R")

# downloads nextgen datasets 
source("runners/cs_runner/download_nextgen.R")

# generate and upload transects datasets 
source("runners/cs_runner/01_transects.R")

# generate and upload cross sections points datasets 
source("runners/cs_runner/02_cs_pts.R")
