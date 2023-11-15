# fs::file_copy('/Users/mjohnson/Library/Application Support/org.R-project.R/R/ngen.hydrofab/RouteLink_nwm_v2_2_3.fst',
#               glue("{base}/RouteLink_nwm_v2_2_3.fst"))
# 
# arrow::write_parquet(tt,  glue("{base}/RouteLink_nwm_v2_2_3.parquet"))

pacman::p_load(hydrofabric, glue, arrow) 
devtools::load_all()

vpus = vpu_boundaries$VPUID[1:21]

base = '/Volumes/MyBook/nextgen'
hldir = glue("{base}/hydrolocations")
dir.create(hldir)

# Extract and Build POIs --------------------------------------------------

for(i in 1:length(vpus)){
  
  VPU = vpus[i]
  
  refactored_gpkg = get_hydrofabric(VPU = VPU, 
                                    type = "refactor",
                                    dir = glue("{base}/refactored"),
                                    overwrite = FALSE)
  
  fl = read_hydrofabric(refactored_gpkg, "flowlines")[[1]] %>% 
    st_transform(4326) 
  
  type =  c('HUC12', 'Gages', 'TE', 'NID', "WBOut")
  poi_layer = grep("POIs_*", st_layers(refactored_gpkg)$name, value = TRUE)
  
  # Community POIs  ----
  
  hl  = read_sf(refactored_gpkg, poi_layer) %>% 
    st_drop_geometry() %>%
    select(hy_id = ID, paste0("Type_", type)) %>%
    mutate_at(vars(matches("Type_")), as.character) %>% 
    pivot_longer(-hy_id) %>%
    filter(!is.na(value)) %>%
    mutate(hl_reference = gsub("Type_", "", name), hl_link = as.character(value)) %>% 
    select(hy_id, hl_reference, hl_link) %>% 
    tidyr::separate_longer_delim(hl_link, ",")
  

  # Coastal Model ----
  
  pts = read.csv(glue('{base}/GAGE_SUMMARY.csv')) %>% 
    select(hl_link = SITE_NO, lat = LAT_NHD, lon = LON_NHD) %>% 
    st_as_sf(coords = c("lon", "lat"), crs = 4326)  %>% 
    mutate(hl_reference = "CoastalGage") %>% 
    mutate(hl_link = as.character(hl_link)) %>% 
    rename_geometry("geometry")
  
  pts$hy_id = fl$ID[st_nearest_feature(pts, fl)]
  
  # RouteLink ----
  
  hs = get_vaa("hydroseq")
  
  lu = read_sf(refactored_gpkg, "lookup_table") %>% 
    select(comid = member_COMID, hy_id = reconciled_ID) %>% 
    mutate(comid = as.numeric(comid))
  
  topo = st_drop_geometry(fl) %>% 
    select(hy_id = ID, toID)
  
  topo2 = left_join(topo, select(topo, hy_id = toID, fromID = hy_id), by = "hy_id")
  
  rl = read_parquet(glue("{base}/RouteLink_nwm_v2_2_3.parquet")) %>% 
    filter(!is.na(NHDWaterbodyComID)) %>% 
    select(comid, to, NHDWaterbodyComID) %>% 
    left_join(hs) %>% 
    group_by(NHDWaterbodyComID, by = "comid") %>% 
    left_join(lu, by = "comid") %>% 
    left_join(topo2, by = "hy_id", relationship = "many-to-many") %>% 
    mutate(WBIn = !hy_id %in% toID,
           WBOut = hydroseq == min(hydroseq)) %>% 
    ungroup() %>% 
    filter(complete.cases(.)) 
  
  #TODO: need to actually add these pre-refactor
  #
  # WBIn = filter(rl, WBIn) %>% 
  #   select(hy_id = fromID, hl_link = NHDWaterbodyComID) %>% 
  #   mutate(hl_reference = "RL_WBIn", hl_link = as.character(hl_link)) %>% 
  #   distinct()
  #   
  WBOut = filter(rl, WBOut) %>% 
    select(hy_id, hl_link = NHDWaterbodyComID) %>% 
    mutate(hl_reference = "RL_WBOut", 
           hl_link = as.character(hl_link)) %>% 
    distinct()
  
  # Merge and Build ----
  
  tmp = bind_rows(hl, st_drop_geometry(pts), WBOut) %>% 
    group_by(hy_id) %>% 
    mutate(hl_id = paste0(VPU, cur_group_id())) %>% 
    ungroup() 
  
  sub = filter(fl, ID %in% tmp$hy_id)
  
  outlets = get_node(sub)
  outlets$hy_id = sub$ID
  
  hydrolocations = full_join(tmp, outlets, by = "hy_id") %>%
    mutate(hl_position = "outflow")
  
  write_sf(hydrolocations, glue("{hldir}/hl_{VPU}.gpkg"))
  
}
