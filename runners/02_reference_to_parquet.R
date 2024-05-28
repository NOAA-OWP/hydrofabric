source("runners/config.R")

ref_gpkg   <- '/Users/mjohnson/Downloads/reference_CONUS.gpkg'
ref_poi <- '/Users/mjohnson/Downloads/reference_CONUS_poigeom.gpkg'
ref_div   <- '/Users/mjohnson/Downloads/reference_catchments.gpkg'


ref_net_out    <- glue('/Users/mjohnson/hydrofabric/v2.2/reference/conus_network')


ref_div_out    <- glue('/Users/mjohnson/hydrofabric/v2.2/reference/conus_divides')
ref_fl_out     <- glue('/Users/mjohnson/hydrofabric/v2.2/reference/conus_flowlines')
ref_hl_out     <- glue('/Users/mjohnson/hydrofabric/v2.2/reference/conus_hydrolocations')

# Reference Fabric --------------------------------------------------------

fl = as_sqlite(ref_gpkg, 'reference_flowline') %>% 
  select(id = comid, 
         toid = tocomid, 
         poi_id = POI_ID, 
         terminalpa, 
         mainstemlp = levelpathi, 
         vpuid, 
         reachcode, 
         frommeas,
         tomeas,
         lengthkm, 
         areasqkm,
         streamorde, 
         totdasqkm, 
         hydroseq, 
         dnhydroseq,
         geom) %>% 
  read_sf_dataset_sqlite()

g = nhdplusTools::get_node(fl, position = "end") %>%
  mutate(outlet_X = st_coordinates(.)[,1], 
         outlet_Y = st_coordinates(.)[,2])

fl = bind_cols(fl, st_drop_geometry(g))

div = as_sqlite(ref_div, 'reference_catchments') %>%
  select(divide_id = featureid, vpuid, areasqkm, geom) %>%
  read_sf_dataset_sqlite()

poi_geom = read_sf(ref_poi, "poi_geometry") %>% 
  mutate(X = st_coordinates(.)[,1], 
         Y = st_coordinates(.)[,2]) %>% 
  st_drop_geometry() %>% 
  distinct() %>% 
  select(-hy_id)

hl = as_sqlite(ref_gpkg, 'poi_data') %>% 
  select(hy_id, hl_link, hl_reference, nat_poi_id, vpuid) %>% 
  collect() %>% 
  mutate(hl_reference = gsub("Type_", "", hl_reference),
         hl_uri = paste0(hl_reference, "-",hl_link)) %>% 
  filter(hl_reference %in% community_hl_types) %>% 
  left_join(poi_geom,
            by = c('nat_poi_id', 'vpuid')) %>%
  select(id = hy_id,
         poi_id = nat_poi_id,
         starts_with('hl_'),
         vpuid, X, Y)

filter(hl, hl_reference == "Gages", hl_link == "06752260") %>% 
  st_as_sf(coords = c("X", "Y"), crs = 5070) %>% 
  mapview::mapview()
  
net = st_drop_geometry(div) %>% 
  mutate(id = divide_id, vpuid = NULL) %>% 
  full_join(st_drop_geometry(select(fl, -areasqkm)), by = "id") %>% 
  mutate(hf_id = id,  topo = "fl-fl") %>% 
  left_join(select(hl, -vpuid, -X, -Y),
            by = "poi_id",
            relationship = "many-to-many")

mutate(div,
       has_flowline = ifelse(divide_id %in% fl$id, TRUE, FALSE),
       id = ifelse(has_flowline, divide_id, NA)) %>% 
  group_by(vpuid) %>% 
  write_sf_dataset(path  = ref_div_out,  
                   version = 2.6,
                   hf_version = "2.2")

mutate(fl,
       has_divide = ifelse(id %in% div$divide_id, TRUE, FALSE),
       divide_id = ifelse(has_divide, id, NA)) %>% 
  group_by(vpuid) %>% 
  write_sf_dataset(path  = ref_fl_out,  
                   version = 2.6,
                   hf_version = "2.2")

group_by(net, vpuid) %>% 
  arrow::write_dataset(ref_net_out, 
                       version = 2.6)

group_by(hl, vpuid) %>% 
  arrow::write_dataset(path  = ref_hl_out,  
                   version = 2.6)


system("aws s3 sync /Users/mjohnson/hydrofabric/v2.2/reference s3://lynker-spatial/hydrofabric/v2.2/reference/")












# Refactor Fabric --------------------------------------------------------

rfc_gpkg   <- glue("{base}/refactor/refactor_CONUS.gpkg")
rfc_net    <- glue("{base}/refactor/conus_network")
rfc_div    <- glue("{base}/refactor/conus_divides")
rfc_fl     <- glue("{base}/refactor/conus_flowlines")


ref = as_sqlite(ref_gpkg, 'reference_network') %>% 
  select(hf_id = comid, 
         mainstemlp = levelpathi, 
         hf_hydroseq = hydroseq) %>% 
  collect() 

lu = as_sqlite(rfc_gpkg, "lookup_table") %>% 
  select(hf_id = nhdplusv2_comid, 
         id = nat_id,
         member_comid) %>% 
  collect()
  
lu2 = distinct(left_join(lu, ref, by = "hf_id")) 

lu3 = group_by(lu2, id) %>% 
  arrange(hf_hydroseq) %>% 
  #rank order member comid from outlet to upstream
  mutate(member_comid = paste(member_comid, collapse = ','))  %>% 
  ungroup()

lu4 = arrange(lu3, id, hf_hydroseq) %>% 
  filter(!duplicated(id)) %>% 
  mutate(hydroseq = rank(hf_hydroseq), hf_id = NULL, hf_hydroseq = NULL)

rm(lu); rm(lu2); rm(lu3); rm(ref); gc()

divide = as_sqlite(rfc_gpkg, 'refactored_divides') %>% 
  collect() %>% 
  rename(id = nat_id, areasqkm = inc_areasqkm, vpuid = vpu) %>% 
  st_as_sf() %>% 
  st_set_crs(5070)

group_by(divide, vpuid) %>% 
  write_sf_dataset(path  = rfc_div,  format = "parquet", hf_version = version,  hive_style = TRUE)

flowpath = as_sqlite(rfc_gpkg, 'reconciled_flowpaths') %>% 
  collect() %>% 
  rename(id = nat_id, toid = to_nat_id, vpuid = vpu) %>% 
  st_as_sf() %>% 
  st_set_crs(5070)

fps = left_join(flowpath, lu4, by = "id") %>% 
  select(id, toid, vpuid, lengthkm, totdasqkm, mainstemlp, hydroseq)

group_by(fps, vpuid) %>% 
  write_sf_dataset(path  = rfc_fl,  format = "parquet", hf_version = version, hive_style = TRUE)

net = st_drop_geometry(divide) %>% 
  select(id, areasqkm) %>% 
  mutate(divide_id = id, vpuid = NULL) %>% 
  full_join(st_drop_geometry(fps), by = "id")

group_by(net, vpuid) %>% 
  arrow::write_dataset(ref_net_out)

gc()
