source("runners/config.R")

ref_gpkg   <- glue('{base}/reference/reference_CONUS.gpkg')
ref_net    <- glue('{base}/reference/conus_network')
ref_div    <- glue('{base}/reference/conus_divides')
ref_fl     <- glue('{base}/reference/conus_flowlines')

rfc_gpkg   <- glue("{base}/refactor/refactor_CONUS.gpkg")
rfc_net    <- glue("{base}/refactor/conus_network")
rfc_div    <- glue("{base}/refactor/conus_divides")
rfc_fl     <- glue("{base}/refactor/conus_flowlines")

# Reference Fabric --------------------------------------------------------

fl = as_sqlite(ref_gpkg, 'reference_flowline') %>%
  select(id = comid, 
         toid = tocomid, 
         poi_id, 
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
  collect() %>%
  st_as_sf()

group_by(fl, vpuid) %>% 
  write_sf_dataset(path  = ref_fl,  format = "parquet", hf_version = version, hive_style = TRUE)

g = nhdplusTools::get_node(fl, position = "end") %>%
  mutate(outlet_X = st_coordinates(.)[,1], outlet_Y = st_coordinates(.)[,2])

net = bind_cols(st_drop_geometry(fl), st_drop_geometry(g))

div = as_sqlite(ref_gpkg, 'reference_catchment') %>%
  select(divide_id = featureid, vpuid, areasqkm, geom) %>%
  collect() %>%
  st_as_sf()

group_by(div, vpuid) %>% 
  write_sf_dataset(path  = ref_div,  format = "parquet", hf_version = version,  hive_style = TRUE)

net = st_drop_geometry(div) %>% 
  mutate(id = divide_id, 
         vpuid = NULL) %>% 
  full_join(net, by = "id") %>% 
  mutate(hf_id = id, topo = "fl-fl")

group_by(net, vpuid) %>% 
  arrow::write_dataset(ref_net)

# Refactor Fabric --------------------------------------------------------

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
  full_join(select(st_drop_geometry(flowpath)), by = "id")




outlet_hf
outlet_X
outlet_Y

group_by(net, vpuid) %>% 
  arrow::write_dataset(rfc_net)

gc()
