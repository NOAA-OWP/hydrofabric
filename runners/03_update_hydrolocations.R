devtools::load_all('.')
library(RNetCDF)
library(dataRetrieval)
library(data.table)

source("runners/config.R")

# Hydrolocations ----------------------------------------------------------

# Purpose --> This code constructs a conus_hl (conus hydrolocation) file from a set of inputs,
# https://github.com/NOAA-OWP/hydrofabric/issues/40

# Hydrolocations have the following schema
# > poi_id: Unique Identifier to aggregate same physical POIs
# > hl_source: Source of hl
# > hl_reference: Type of Hydrolocation
# > hl_link: Identifier from source dataset
# 
# > X: Longitude
# > Y: Latitude
#
# > hf_id: ID from reference hydrofabric
# > hf_source: Source of reference hydrofabric

# Reference Fabric --------------------------------------------------------
# fl = as_sqlite(ref_gpkg, 'reference_flowline') %>%
#   select(comid, tocomid, poi_id, levelpathi, vpuid, hydroseq,  geom) %>%
#   collect() %>%
#   st_as_sf()
# 
# g = nhdplusTools::get_node(fl, position = "end") %>%
#   mutate(outlet_X = st_coordinates(.)[,1],
#          outlet_Y = st_coordinates(.)[,2])
# 
# ref = bind_cols(st_drop_geometry(fl), st_drop_geometry(g))
# 
# write_parquet(ref, '/Volumes/MyBook/TNC/v2.2/reference/conus_network.parquet')

source = '/Users/mjohnson/hydrofabric'

ref_net  <- collect(open_dataset(glue('{source}/v2.2/reference/conus_network')))

# ref_gpkg = '/Volumes/MyBook/TNC/v2.2/reference/reference_CONUS.gpkg'

# GFv20 Reference -------------------------------------------------------------------

# meta = sf::read_sf(ref_gpkg, 'poi_data') %>% 
#   filter(hl_reference %in% paste0('Type_', community_hl_types)) %>% 
#   select(hf_id = hy_id, 
#          hl_id = poi_id, 
#          hl_link, 
#          hl_reference)  %>% 
#   mutate(hf_source = "reference_features", 
#          hl_source = "GFv20", 
#          hl_reference = gsub("Type_", "", hl_reference))
# 
# 
# community_pois = select(sf::read_sf(ref_gpkg, 'poi_geometry'),  
#                         hf_id = hy_id, 
#                         hl_id = poi_id, ) %>% 
#   mutate(X = st_coordinates(.)[ ,1],
#          Y = st_coordinates(.)[ ,2]) %>% 
#   st_drop_geometry() %>% 
#   right_join(meta, by = c('hf_id', 'hl_id'), relationship = "many-to-many") %>% 
#   # TMP for now
#   mutate(hl_id = NA) %>% 
#   relocate(all_of(schema)) 

gfv20 = open_dataset(glue("{source}/v2.2/reference/conus_hydrolocations")) %>% 
  collect() %>% 
  mutate(hf_id = id,
         hf_source = "reference_features", 
         hl_source = "GFv20",
         usgs_poi_id = NULL) %>% 
  select(all_of(schema)) %>% 
  relocate(all_of(schema)) %>% 
  distinct() 

# RouteLink ---------------------------------------------------------------

rl_file = glue('/Users/mjohnson/hydrofabric/RouteLink_CONUS_309.nc')

nc = open.nc(rl_file)

var = c('link', 'NHDWaterbodyComID')

rl_pois <- lapply(1:length(var), FUN = function(x){
  var.get.nc(nc, var[x])}) %>% 
  bind_cols(.name_repair = "unique_quiet") |>
  setNames(c("comid", "NHDWaterbodyComID")) |>
  filter(NHDWaterbodyComID > 0) |>
  left_join(select(ref_net, comid = hf_id, hydroseq), by = "comid") |>
  filter(!is.na(hydroseq)) |> 
  tidyr::drop_na() |> 
  group_by(NHDWaterbodyComID) |>
  arrange(hydroseq)|>
  filter(row_number()==1 | row_number()==n()) |> 
  mutate(WBOut_rl = ifelse(hydroseq == min(hydroseq), comid, NA),
         WBIn_rl  = ifelse(hydroseq == max(hydroseq), comid, NA)) |>
  ungroup() |>
  select(hl_link = NHDWaterbodyComID, WBOut_rl, WBIn_rl) |> 
  pivot_longer(-hl_link, values_to = "hf_id", names_to = "hl_reference") |> 
  tidyr::drop_na() |> 
  mutate(poi_id = NA,  hf_source = "reference_features", hl_source = "nwm_v3.0.9_routelink") |> 
  left_join(distinct(select(ref_net, hf_id, X = outlet_X, Y = outlet_Y)), by = 'hf_id') %>%
  select(all_of(schema)) %>% 
  relocate(all_of(schema)) %>% 
  distinct() 

# NWM Lakes  --------------------------------------------------------------

nc  <- open.nc(glue('{source}/LAKEPARM_CONUS_216.nc'))

var <- c("lake_id")

lk <- lapply(1:length(var), FUN = function(x){
  var.get.nc(nc, var[x])}) %>% 
  bind_cols(.name_repair = "unique_quiet") |>
  setNames(c("hl_link")) |> 
  mutate(hl_reference = "nwmlake",
         poi_id = NA, 
         hf_source = "reference_features",
         hl_source = "nwm_v2.1.6_LAKEPARAM")

nwmlake_pois <- distinct(select(rl_pois, hl_link, hf_id, hl_reference, X, Y))  %>% 
  filter(hl_reference == 'WBOut_rl') |>
  select(-hl_reference) |>
  right_join(lk, by = "hl_link") %>% 
  tidyr::drop_na(X) %>% 
  select(all_of(schema)) %>% 
  relocate(all_of(schema)) %>% 
  distinct() 

# NWM Reservoirs  --------------------------------------------------------------

f <- list.files(glue('/Users/mjohnson/hydrofabric'), 
                pattern = "reservoir", 
                full.names = TRUE)

res <- list()

for(i in 1:length(f)){
  nc <- open.nc(f[i])
  
  type <- paste0("nwm_v3.0.9_", gsub(".nc", "", basename(f[i])))
  
  usgs_var <- c("usgs_lake_id", 'usgs_gage_id')
  
  usgs_lake <- lapply(1:length(usgs_var), FUN = function(x){
    var.get.nc(nc, usgs_var[x])}) %>% 
    bind_cols(.name_repair = "unique_quiet") %>% 
    setNames(c("hl_link", "hl_id")) %>% 
    mutate(hl_source = type,
           hl_reference = "usgs_gage") %>% 
    select(starts_with("hl_")) %>% 
    left_join(distinct(select(filter(rl_pois, hl_reference == "WBOut_rl"), 
                              hl_link, X, Y, starts_with("hf_"))), 
              by = "hl_link") %>% 
    mutate(hl_link = hl_id,
           poi_id = NA) 
  
  rfc_var <- c("rfc_lake_id", 'rfc_gage_id')
  
  rfc_lake <- lapply(1:length(rfc_var), FUN = function(x){
    var.get.nc(nc, rfc_var[x])}) %>% 
    bind_cols(.name_repair = "unique_quiet") %>% 
    setNames(c("hl_link", "hl_id")) %>% 
    mutate(hl_source = type,
           hl_reference = "rfc_gage") %>% 
    select(starts_with("hl_")) %>% 
    left_join(distinct(select(filter(rl_pois, hl_reference == "WBOut_rl"), 
                              hl_link, X, Y, starts_with("hf_"))), 
              by = "hl_link") %>% 
    mutate(hl_link = hl_id,
           hl_id = NA) 
  
  usace_var <- c("usace_lake_id", 'usace_gage_id')
  
  usace_lake <- lapply(1:length(usace_var), FUN = function(x){
    var.get.nc(nc, usace_var[x])}) %>% 
    bind_cols(.name_repair = "unique_quiet") %>% 
    setNames(c("hl_link", "hl_id")) %>% 
    mutate(hl_source = type,
           hl_reference = "usace_gage") %>% 
    select(starts_with("hl_")) %>% 
    left_join(distinct(select(filter(rl_pois, hl_reference == "WBOut_rl"), 
                              hl_link, X, Y, starts_with("hf_"))), 
              by = "hl_link") %>% 
    mutate(hl_link = hl_id,
           poi_id = NA) %>% 
    select(all_of(schema)) %>% 
    distinct() 
  
  
  res[[i]] <- bind_rows(usgs_lake, rfc_lake, usace_lake)
}

nwm_res <- bind_rows(res) %>% 
  select(all_of(schema)) %>% 
  relocate(all_of(schema)) %>% 
  distinct() 


# Calibration Gages -------------------------------------------------------

cal <- read.csv(glue('{source}/gage_list.txt')) %>% 
  mutate(site_no = paste0("0", site_no))

tmp <- filter(gfv20, hl_reference == "Gages") %>% 
  filter(hl_link %in% cal$site_no)

stopifnot(nrow(tmp) == nrow(cal))

calib_poi <- tmp |>
  mutate(hl_id = NA,
         hl_source = "noaaowp",
         hl_reference = "nwm_v3.0.7_calibration_gage") |>
  select(all_of(schema)) %>% 
  relocate(all_of(schema)) %>% 
  distinct() 


# Coastal Gages -----------------------------------------------------------

coastal_network <- read.csv(glue('{source}/GAGE_SUMMARY.csv')) |>
  mutate(hl_link = sprintf('%08s', SITE_NO), 
         hf_id = NA,
         hl_id = NA,
         hl_source = "coastal",
         hl_reference = "coastal_gage",
         hf_source = "reference_features") |>
  select(starts_with(c("hl_", "hf_")), Y = LAT_NHD, X = LON_NHD) 

for(i in 1:nrow(coastal_network)){
  coastal_network$hf_id[i] <- as.numeric(
    findNLDI(location = c(coastal_network$X[i],
                          coastal_network$Y[i]))$origin$comid
  )
}

coastal_gage_poi <- coastal_network %>% 
  st_as_sf(coords = c("X", "Y"), crs = 4326) %>% 
  st_transform(5070) %>% 
  mutate(X = st_coordinates(.)[,1], 
         Y =  st_coordinates(.)[,2],
         poi_id = NA) %>% 
  st_drop_geometry() %>% 
  select(all_of(schema)) %>% 
  relocate(all_of(schema)) %>% 
  distinct() 

# Coastal Domain ----------------------------------------------------------

coastal_domain <- read_sf(glue('{source}/coastal_domain.gpkg')) |>
  st_transform(5070) |>
  st_cast("MULTILINESTRING")

m <-  st_intersects(st_transform(vpu_boundaries[1:21,], 5070), coastal_domain)

v = vpu_boundaries[1:21,]$VPUID[lengths(m) > 0]

domain_pois <- list()

for(i in 1:length(v)){
  
  fl_tmp <-  open_dataset(ref_fl) %>% 
    dplyr::filter(vpuid == !!v[i]) %>% 
    read_sf_dataset()

  d <- st_intersection(coastal_domain, AOI::bbox_get(fl_tmp))
  
  touches <- st_intersects(fl_tmp, d)
  
  fl_map <- fl_tmp[lengths(touches) > 0, ]
  
  domain_pois[[i]] <- st_intersection(coastal_domain, fl_map) |>
    select(hf_id = id, domain_id) |>
    st_collection_extract("POINT") |>
    st_cast("POINT") |>
    distinct() |>
    group_by(hf_id) |>
    mutate(hl_link = paste(domain_id[1], hf_id, sep = "-")) |>
    slice(1) %>% 
    ungroup() |>
    mutate(
      poi_id = NA,
      hl_source = "coastal",
      hl_reference = domain_id,
      hf_source = "reference_features",
      domain_id = NULL,
    ) %>% 
    mutate(X = st_coordinates(.)[,1],
           Y = st_coordinates(.)[,2]) |> 
    st_drop_geometry()
  
  message("VPU", v[i], " (", i, "/", length(v), ")")
}

coastal_domain_pois <- bind_rows(domain_pois) %>% 
  select(all_of(schema)) %>% 
  relocate(all_of(schema)) %>% 
  distinct() 

# lid --------------------------------------------------------------------

lid <- read_sf(glue('{source}/nws_lid.gpkg')) |>
  mutate(hf_id = as.numeric(nwm_feature_id)) |>
  select(nws_lid, usgs_site_code, hf_id) |>
  pivot_longer(-c(hf_id, geom), names_to = "hl_reference", values_to = "hl_link") |>
  filter(!is.na(hl_link)) |>
  mutate(
    poi_id = NA,
    hl_source = "fim",
    hf_source = "reference_features",
  )  %>%  
  mutate(X = st_coordinates(.)[,1],
         Y = st_coordinates(.)[,2]) |> 
  st_drop_geometry() |>
  select(all_of(schema)) %>% 
  relocate(all_of(schema)) %>% 
  distinct() 

# Final  ------------------------------------------------------------------

hl = rbindlist(list(rl_pois,
                    nwmlake_pois,
                    nwm_res,
                    calib_poi,
                    coastal_gage_poi,
                    coastal_domain_pois,
                    lid
)) %>% 
  select(-X, -Y) %>% 
  left_join(distinct(select(ref_net, X = outlet_X, Y = outlet_Y, hf_id = id)), 
            by = "hf_id",
            relationship = "many-to-many") %>% 
  bind_rows(gfv20) %>% 
  group_by(X, Y) %>% 
  mutate(count = n(),
         poi_id = cur_group_id()) %>% 
  ungroup() %>% 
  left_join(distinct(select(ref_net, hf_id, mainstemlp, vpuid)), 
            by = "hf_id",
            relationship = "many-to-many") %>% 
  distinct()

# Write Hydrolocation Table
hl %>% 
  group_by(vpuid) %>% 
  arrow::write_dataset(glue("{source}/conus_hl"), version = 2.6)

aws s3 sync /Users/mjohnson/hydrofabric/conus_hl s3://lynker-spatial/hydrofabric/conus_hl

open_dataset(glue("{source}/conus_hl")) %>% collect()


