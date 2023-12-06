source("runners/config.R")

if(FIX){ pipeline = filter(pipeline, !is.na(corrected_refactor)) }
global = TRUE
# Produce Nextgen Fabrics -------------------------------------------------

for(i in 1:nrow(pipeline)){
  
  #if(!file.exists(pipeline$nextgen[i])){
    
    o    = apply_nexus_topology(gpkg        = ifelse(global, pipeline$uniform_global[i], pipeline$uniform[i]),
                                nexus_prefix = nexus_prefix,
                                terminal_nexus_prefix = terminal_nexus_prefix,
                                coastal_nexus_prefix = coastal_nexus_prefix,
                                internal_nexus_prefix = internal_nexus_prefix,
                                catchment_prefix = catchment_prefix,
                                waterbody_prefix = waterbody_prefix,
                                vpu         = pipeline$vpus[i],
                                export_gpkg = pipeline$nextgen[i]) %>%
      add_flowpath_attributes() %>%
      add_lake_attributes(lake_path = lake_path) %>%
      append_style(layer_names = c("nexus", "hydrolocations", "flowpaths", "divides", "lakes"))

    bas = read_sf(o, "divides")
    bas = clean_geometry(bas, ID = "divide_id")
    write_sf(bas, o, "divides")
  #}
}

national_merge(gpkg = pipeline$nextgen, outfile = conus_gpkg )

# Build National Topology -------------------------------------------------
vaa = get_vaa(c('hydroseq', 'areasqkm'), updated_network = TRUE) %>%
  rename(hf_id = comid, hf_hydroseq = hydroseq, hf_areasqkm = areasqkm)

r = data.table::fread(ms_lookup) %>% 
  select(mainstem = lp_mainstem, ref_mainstem_uri = ref_mainstem_id) %>% 
  mutate(ref_mainstem_uri = glue('https://geoconnex.us/ref/mainstems/{ref_mainstem_uri}'))

net = left_join(read_sf(conus_gpkg, 'network'), vaa, by = "hf_id") %>% 
  left_join(r, by = "mainstem")

write_parquet(net, conus_net)

hl = read_sf(conus_gpkg, "hydrolocations")
write_sf(hl, glue('{base_fgb}/hydrolocations.fgb'), overwrite = TRUE)

nexus = read_sf(conus_gpkg, "nexus")
write_sf(nexus, glue('{base_fgb}/nexus.fgb'), overwrite = TRUE)

fps = read_sf(conus_gpkg, "flowpaths")
write_sf(fps, glue('{base_fgb}/flowpaths.fgb'), overwrite = TRUE)

divides = read_sf(conus_gpkg, "divides")
write_sf(divides, glue('{base_fgb}/divides.fgb'), overwrite = TRUE)
