source("runners/config.R")

s <- '/Users/mikejohnson/hydrofabric/v2.2'

ref_gpkg   <- glue('{s}/reference_CONUS.gpkg')
# ref_poi    <- glue('/Users/mjohnson/Downloads/reference_CONUS_poigeom.gpkg'
ref_div    <- glue('{s}/reference_catchments.gpkg')

ref_net_out    <- glue('{s}/reference/conus_network')
ref_div_out    <- glue('{s}/reference/conus_divides')
ref_fl_out     <- glue('{s}/reference/conus_flowlines')
ref_hl_out     <- glue('{s}/reference/gfv2_hydrolocations')

# Reference Fabric --------------------------------------------------------

fl = as_sqlite(ref_gpkg, 'reference_flowline') %>% 
  select(id = comid, 
         toid = tocomid, 
         #poi_id = POI_ID, 
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
         outlet_Y = st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

g2 = nhdplusTools::get_node(fl, position = "start") %>%
  mutate(inlet_X = st_coordinates(.)[,1], 
         inlet_Y = st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

fl = bind_cols(fl, bind_cols(g, g2))

div = as_sqlite(ref_div, 'reference_catchments') %>%
  select(divide_id = featureid, vpuid, areasqkm, geom) %>%
  read_sf_dataset_sqlite()

div = st_set_crs(div, st_crs(fl))

poi_geom = read_sf(ref_gpkg, "poi_geometry") %>% 
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
            by = "id",
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
  arrow::write_dataset(ref_net_out, version = 2.6)

group_by(hl, vpuid) %>% 
  arrow::write_dataset(path  = ref_hl_out,  
                   version = 2.6)

system(glue("aws s3 sync {s}/reference s3://lynker-spatial/hydrofabric/v2.2/reference/"))

# Hawaii ------------------------------------------------------------------

ref_gpkg_hi   <- '/Users/mjohnson/hydrofabric/reference_20.gpkg'

ref_net_out_hi    <- glue('/Users/mjohnson/hydrofabric/v2.2/reference/hawaii_network')
ref_div_out_hi   <- glue('/Users/mjohnson/hydrofabric/v2.2/reference/hawaii_divides')
ref_fl_out_hi     <- glue('/Users/mjohnson/hydrofabric/v2.2/reference/hawaii_flowlines')
ref_hl_out_hi     <- glue('/Users/mjohnson/hydrofabric/v2.2/reference/hawaii_hydrolocations')

# Reference Fabric --------------------------------------------------------
as_sqlite(ref_gpkg_hi)

fl = as_sqlite(ref_gpkg_hi, 'nhd_flowline') %>% 
  select(id = COMID, 
         toid = toCOMID, 
         #poi_id = POI_ID, 
         terminalpa = TerminalPa, 
         mainstemlp = LevelPathI, 
         vpuid = VPUID, 
         reachcode = REACHCODE, 
         frommeas = FromMeas,
         tomeas = ToMeas,
         lengthkm = LENGTHKM, 
         areasqkm = AreaSqKM,
         streamorde = StreamOrde, 
         totdasqkm = TotDASqKM, 
         hydroseq = Hydroseq, 
         dnhydroseq = DnHydroseq,
         geom) %>% 
  read_sf_dataset_sqlite()

g = nhdplusTools::get_node(fl, position = "end") %>%
  mutate(outlet_X = st_coordinates(.)[,1], 
         outlet_Y = st_coordinates(.)[,2])

fl = bind_cols(fl, st_drop_geometry(g))

div = as_sqlite(ref_gpkg_hi, 'nhd_catchment') %>%
  select(divide_id = FEATUREID,  areasqkm = AreaSqKM, geom) %>%
  mutate(vpuid = 20) %>% 
  read_sf_dataset_sqlite()

hl = read_sf(ref_gpkg_hi, "final_POIS") %>% 
  mutate(X = st_coordinates(.)[,1], 
         Y = st_coordinates(.)[,2]) %>% 
  st_drop_geometry() %>% 
  mutate(across(starts_with('Type_'), as.character)) %>% 
  tidyr::pivot_longer(-c(COMID, X,Y), values_to = "hl_link", names_to = "hl_reference") %>% 
  tidyr::drop_na(hl_link) %>% 
  mutate(hl_reference = gsub("Type_", "", hl_reference),
         hl_uri = paste0(hl_reference, "-",hl_link),
         vpuid = 20) %>%
  filter(hl_reference %in% c(community_hl_types, "Elev")) %>% 
  group_by(X,Y) %>% 
  mutate(poi_id = cur_group_id()) %>% 
  ungroup() %>% 
  select(id = COMID,
         poi_id,
         starts_with('hl_'),
         vpuid, X, Y)
  
net = st_drop_geometry(div) %>% 
  mutate(id = divide_id, vpuid = NULL) %>% 
  full_join(st_drop_geometry(select(fl, -areasqkm)), by = "id") %>% 
  mutate(hf_id = id,  topo = "fl-fl") %>% 
  left_join(select(hl, -vpuid, -X, -Y),
            by = "id",
            relationship = "many-to-many")

mutate(div,
       has_flowline = ifelse(divide_id %in% fl$id, TRUE, FALSE),
       id = ifelse(has_flowline, divide_id, NA)) %>% 
  write_sf_dataset(path  = ref_div_out_hi,  
                   version = 2.6,
                   hf_version = "2.2")

mutate(fl,
       has_divide = ifelse(id %in% div$divide_id, TRUE, FALSE),
       divide_id = ifelse(has_divide, id, NA)) %>% 
  write_sf_dataset(path  = ref_fl_out_hi,  
                   version = 2.6,
                   hf_version = "2.2")

net %>% 
  arrow::write_dataset(ref_net_out_hi, 
                       version = 2.6)

hl %>% 
  arrow::write_dataset(path  = ref_hl_out_hi,  
                       version = 2.6)


system("aws s3 sync {s}/reference s3://lynker-spatial/hydrofabric/v2.2/reference/")
system("aws s3 sync s3://lynker-spatial/hydrofabric/v2.2/reference/ {s}/reference")


open_dataset('s3://lynker-spatial/hydrofabric/v2.2/reference/hawaii_divides') %>% 
  read_sf_dataset() %>% 
  mapview::mapview()






