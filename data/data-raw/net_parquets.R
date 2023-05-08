
# Build National Topology -------------------------------------------------
library(sf); library(arrow)
"/Volumes/Transcend/ngen/CONUS-hydrofabric/05_nextgen/conus.gpkg" %>% 
  read_sf('network_lookup') %>% 
  write_parquet("data/conus_net.parquet")

"/Volumes/Transcend/ngen/CONUS-hydrofabric/05_nextgen/conus.gpkg" %>% 
  read_sf('hydrolocation_lookup') %>% 
  write_parquet("data/conus_hl.parquet")