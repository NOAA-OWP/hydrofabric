
# Build National Topology -------------------------------------------------
library(sf); library(arrow)
"/Volumes/Transcend/ngen/CONUS-hydrofabric/05_nextgen/conus.gpkg" %>% 
  read_sf('network') %>% 
  write_parquet("data/conus_net.parquet")
