
# Build National Topology -------------------------------------------------
library(sf); library(arrow)
vaa = get_vaa("hydroseq") %>% 
  rename(hf_id = comid, hf_hydroseq = hydroseq)

xx = "/Volumes/Transcend/ngen/CONUS-hydrofabric/05_nextgen/conus.gpkg" %>% 
  read_sf('network') %>% 
  left_join(vaa, by = "hf_id")
  
write_parquet(xx, "data/conus_net.parquet")
