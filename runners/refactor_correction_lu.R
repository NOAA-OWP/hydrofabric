source("runners/config.R")

old = "11172023"

old_rds = glue("runners/data/refactored_{old}_topology.rds")
new_rds = glue("runners/data/refactored_{ref_date}_topology.rds")


df = data.frame(old = list.files(glue('{base}/refactored_{old}'), 
                 pattern = ".gpkg$",
                 full.names = TRUE),
                 new = list.files(glue('{base}/refactored_{ref_date}'), 
                                 pattern = ".gpkg$",
                                 full.names = TRUE)) %>% 
  mutate(vpu = gsub("refactor_", "", gsub(".gpkg", "", basename(old))))


new_topo = old_topo = list()

for(i in 1:nrow(df)){
  message(i)
   old_topo[[i]] = 
       read_sf(df$old[i], 'lookup_table') %>% 
         mutate(vpu = df$vpu[i])
   
   new_topo[[i]] = 
     read_sf(df$new[i], 'lookup_table') %>% 
     mutate(vpu = df$vpu[i])
     
}


saveRDS(bind_rows(old_topo), old_rds)
saveRDS(bind_rows(new_topo), new_rds)

o = readRDS(old_rds) %>% 
  rename(id = reconciled_ID,
         old_mainstem = LevelPathID,
         VPU = vpu) %>% 
  mutate(mainstem = NULL)

n = readRDS(new_rds) %>% 
  rename(new_id = reconciled_ID,
         new_mainstem = LevelPathID,
         new_VPU = vpu) %>% 
  mutate(mainstem = NULL)


lu = full_join(o, n,  by = c("NHDPlusV2_COMID"),  relationship = "many-to-many")

saveRDS(lu, "runners/data/refactor_lu.rds")

