source("config.R")

# Extract and Build POIs --------------------------------------------------

for(i in 1:length(vpus)){
  
  VPU = vpus[i]
  
  refactored_gpkg = get_hydrofabric(VPU = VPU, type = "refactor", dir = base_refactored)
  reference_gpkg  = get_hydrofabric(VPU = VPU, type = "reference", dir = base_reference)
  
  fl = read_hydrofabric(reference_gpkg, "flowlines")[[1]] %>% 
    st_transform(4326) 

  poi_layer = grep("POIs_*", st_layers(reference_gpkg)$name, value = TRUE)

  # Community POIs  ----
  
  hl  = read_sf(reference_gpkg, poi_layer) %>% 
    st_transform(4326) %>% 
    mutate(X = st_coordinates(.)[,1],
           Y = st_coordinates(.)[,2]) %>% 
    st_drop_geometry() %>%
    select(hf_id = COMID, X, Y, poi_id = id,  any_of(paste0("Type_", community_hl_types))) %>%
    mutate_at(vars(matches("Type_")), as.character) %>% 
    tidyr::pivot_longer(-c(poi_id, hf_id, X, Y)) %>%
    filter(!is.na(value)) %>%
    mutate(hl_reference = gsub("Type_", "", name), hl_link = as.character(value)) %>% 
    select(hf_id, poi_id, hl_reference, hl_link, X, Y) %>% 
    tidyr::separate_longer_delim(hl_link, ",")
  

  # Coastal Model ----
  
  coastal = read.csv(coastal_gages) %>% 
    select(hl_link = SITE_NO, lat = LAT_NHD, lon = LON_NHD) %>% 
    st_as_sf(coords = c("lon", "lat"), crs = 4326)  %>% 
    mutate(hl_reference = "CoastalGage") %>% 
    mutate(hl_link = as.character(hl_link)) %>% 
    rename_geometry("geometry") %>% 
    mutate(X = st_coordinates(.)[,1], Y = st_coordinates(.)[,2])
  
  coastal$hf_id = fl$COMID[st_nearest_feature(coastal, fl)]
  
  # RouteLink ----
  
  hs = get_vaa("hydroseq")
  
  lu = read_sf(refactored_gpkg, "lookup_table") %>% 
    select(comid = member_COMID, hy_id = reconciled_ID) %>% 
    mutate(comid = as.numeric(comid))
  
  rl = open_dataset(rl_file) %>% 
    select(comid, NHDWaterbodyComID) %>% 
    collect() %>% 
    filter(!is.na(NHDWaterbodyComID)) %>% 
    left_join(hs, by = "comid") %>% 
    filter(!is.na(hydroseq)) %>% 
    group_by(NHDWaterbodyComID) %>% 
    slice_min(hydroseq) %>% 
    ungroup() %>%  
    select(hf_id = comid, hl_link = NHDWaterbodyComID) %>% 
    mutate(hl_reference = "RL_WBOut", 
           hl_link = as.character(hl_link)) %>% 
    distinct() %>% 
    filter(hf_id %in% fl$COMID)
  
  
  # AHPS ----
  a = read_sf(fim_ahps) %>% 
    filter(substring(HUC8, 1,2) == substring(VPU, 1,2)) %>% 
    filter(!usgs_site_code %in% filter(hl, hl_reference == "Gages")$hl_link) %>% 
    st_transform(st_crs(fl))
  
  a$hf_id = fl$COMID[st_nearest_feature(a, fl)]
  
  a$hf_id == a$nwm_feature_id
  
  a2 = a %>% 
    st_transform(4326) %>% 
    mutate(X = st_coordinates(.)[,1],
           Y = st_coordinates(.)[,2]) %>% 
    st_drop_geometry() %>% 
    select(nws_lid, usgs_site_code, hf_id, X, Y) %>% 
    tidyr::pivot_longer(-c(hf_id, X, Y)) %>% 
    filter(complete.cases(.)) %>% 
    mutate(name = ifelse(name == "usgs_site_code", "Gages", name)) %>% 
    mutate(name = ifelse(name == "nws_lid", "AHPS", name)) %>% 
    rename(hl_link = value, hl_reference = name)
  
  # Merge and Build ----
  
  tmp = bind_rows(hl, st_drop_geometry(coastal), rl, a2) %>% 
    group_by(hf_id) %>% 
    mutate(hl_id = paste0(VPU, cur_group_id())) %>% 
    ungroup() 
  
  sub = filter(fl, COMID %in% tmp$hf_id)
  
  outlets = get_node(sub)
  outlets$hf_id = sub$COMID
  
  lu = read_sf(refactored_gpkg, "lookup_table") %>% 
    select(hf_id = NHDPlusV2_COMID, ID = reconciled_ID, mainstem = LevelPathID, member_COMID) %>% 
    group_by(hf_id) %>% 
    slice_min(as.numeric(member_COMID)) %>% 
    ungroup() %>% 
    select(-member_COMID)

  hydrolocations = full_join(tmp, outlets, by = "hf_id") %>%
    left_join(lu, by = "hf_id") %>% 
    mutate(hl_position = "outflow") %>% 
    st_as_sf()
  
  hydrolocations = mutate(hydrolocations,
              X = ifelse(is.na(X), st_coordinates(hydrolocations)[,1], X),
              Y = ifelse(is.na(X), st_coordinates(hydrolocations)[,1], Y)) %>% 
    select(hl_link, hl_reference, hf_id, ID,  X, Y) %>% 
    mutate(VPUID = VPU)

  write_sf(hydrolocations, glue("{base_hl}/hl_{VPU}.gpkg"))
  
}

l = bind_rows(lapply(list.files(base_hl, full.names = TRUE), read_sf))

l$hl_id = 1:nrow(l)

write_sf(l, full_hl)
