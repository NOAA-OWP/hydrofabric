source("runners/config.R")

camels_data  = read.csv(camels) |>
  select(gauge_id, COMID) |>
  mutate(outfile = glue("{base_camels}/Gage_{gauge_id}.gpkg")) %>% 
  mutate(id = 1:n())

vaa = get_vaa(c('hydroseq', 'areasqkm'), updated_network = TRUE) %>%
  rename(hf_id = comid, hf_areasqkm = areasqkm)

#l = st_layers(list.files(base_gpkg, full.names = TRUE)[1])

for(i in 1:nrow(camels_data)){
  
  message(crayon::yellow(i, "of", nrow(camels_data)))
  
  if(!file.exists(camels_data$outfile[i])){
    
    hf = subset_network(comid = camels_data$COMID[i],
                        base_dir = glue("{base}/{version}"),
                        outfile = camels_data$outfile[i],
                        lyrs  = c(
                          "divides",
                          "nexus",
                          "flowpaths",
                          "network",
                          "hydrolocations",
                          "flowpath_attributes"
                        ))
    
    net = read_sf(hf, "network") %>% 
      left_join(vaa, by = "hf_id") %>% 
      select(divide_id, hf_id, hf_areasqkm, areasqkm) %>% 
      filter(complete.cases(.)) %>% 
      mutate(r = hf_areasqkm / areasqkm) %>% 
      select(divide_id, hf_id, r)
    
    ## EXRTACT HydroATLAS vars
    ha = open_dataset(glue('{base}/hydroatlas_vars.parquet')) %>% 
      filter(hf_id %in% net$hf_id) %>% 
      collect() %>% 
      right_join(net, by = "hf_id") %>% 
      select(-hydroatlas_id, -hf_id)
    
    cols = names(ha)[!names(ha) %in% names(net)]
    
    tmap = mutate(ha, across(dplyr::any_of(cols), ~ .x * r)) %>% 
      group_by(divide_id) %>% 
      summarise(across(dplyr::any_of(cols), ~ sum(.x))) %>% 
      ungroup()
    
    write_sf(tmap, hf, "hydroatlas_atts")
  }
}


st_layers(hf)
