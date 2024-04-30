devtools::load_all('.')
library(httr)
library(RNetCDF)
library(dataRetrieval)
library(data.table)

# Purpose --> This code constructs a conus_hl (conus hydrolocation) file from a set of inputs,
# https://github.com/NOAA-OWP/hydrofabric/issues/40

map_poi = function(df){
  df %>% 
    sf::st_as_sf(coords = c("X", "Y"), crs = 4326) %>% 
    mapview::mapview()
}

# Hydrolocations ----------------------------------------------------------
# Hydrolocations have the following schema
# > hl_id: Unique Identifier to aggregate same physical POIs
# > hl_source: Source of hl
# > hl_reference: Type of Hydrolocation
# > hl_link: Identifier from soure dataset
# > X: Longitude
# > Y: Latiture
# > hf_id: ID from reference hydrofabric
# > hf_source: Source of reference hydrofabric

# Target Schema
schema <- c("hl_id", "hl_source", "hl_reference", "hl_link", "X", "Y",
            "hf_id", "hf_source")

# Reference Fabric --------------------------------------------------------
ref_base <- "/Volumes/MyBook/Reference Features - targets/"
ref_net  <- read_parquet(glue('{ref_base}reference_features.parquet'))

# GFv20 -------------------------------------------------------------------

sb_base <- "https://www.sciencebase.gov/catalog/file/get/60be0e53d34e86b93891012b?f=__disk__"

GFv20_POIs_type <- glue('{sb_base}88/be/47/88be47f9e53f44eda6fdb8f1e4e0c79f85e33d06')
GFv20_POIs      <- glue('{sb_base}9b/05/b0/9b05b0e75e91d1692adcb1c2398e75e751865023')
gfv20_type      <- glue("{base}/hydrolocations/GFv20_POIs_type.csv")
gfv20_geo       <- glue("{base}/hydrolocations/GFv20_POIs.geojson")

try( GET(url = GFv20_POIs_type, 
         write_disk(gfv20_type, overwrite = TRUE)), 
     silent = TRUE )

try( GET(url = GFv20_POIs, 
         write_disk(gfv20_geo, overwrite = TRUE)), 
     silent = TRUE )

meta = fread(gfv20_type) %>% 
  mutate(Type_Gages = stringi::stri_pad_left(Type_Gages, 8, "0"),
         Type_HUC12 = stringi::stri_pad_left(Type_HUC12, 12, "0")) %>%  
  select(hf_id = nhdpv2_COMID, 
         hl_id = provider_id,  
         any_of(paste0("Type_", community_hl_types)), Y, X) |>
  mutate_at(vars(matches("Type_")), as.character) |>
  pivot_longer(-c(hl_id, hf_id, X, Y), names_to = "hl_reference", values_to = "hl_link") |>
  filter(!is.na(hl_link)) |>
  distinct() %>% 
  mutate(hf_source = "reference_features", hl_source = "GFv20", hl_reference = gsub("Type_", "", hl_reference))

community_pois =  select(read_sf(gfv20_geo), hl_id = provider_id) %>% 
  mutate(X = st_coordinates(.)[ ,1],
         Y = st_coordinates(.)[ ,2]) %>% 
  st_drop_geometry() %>% 
  right_join(meta, by = "hl_id") %>% 
  # TMP for now
  mutate(hl_id = NA) %>% 
  #left_join(select(ref_net, hf_id = comid, hf_hydroseq = hydroseq), by = 'hf_id') %>% 
  relocate(all_of(schema)) 

# RouteLink ---------------------------------------------------------------

rl_file = glue('{base}/nwm_v3.0.7/conus/RouteLink_CONUS.nc')

nc = open.nc(rl_file)

var = c('link', 'NHDWaterbodyComID')

rl_pois <- lapply(1:length(var), FUN = function(x){
  var.get.nc(nc, var[x])}) %>% 
  bind_cols(.name_repair = "unique_quiet") |>
  setNames(c("comid", "NHDWaterbodyComID")) |>
  filter(NHDWaterbodyComID > 0) |>
  left_join(get_vaa("hydroseq", updated_network = TRUE), by = "comid") |>
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
  mutate(hl_id = NA,  hf_source = "reference_features", hl_source = "nwm_v3.0.7_routelink") |> 
  left_join(select(ref_net, hf_id = comid, X, Y), by = 'hf_id') %>%
  distinct() |>
  relocate(all_of(schema))

# NWM Lakes  --------------------------------------------------------------

nc  <- open.nc(lake_path)

var <- c("lake_id", "lat", "lon")

lk <- lapply(1:length(var), FUN = function(x){
  var.get.nc(nc, var[x])}) %>% 
  bind_cols(.name_repair = "unique_quiet") |>
  setNames(c("hl_link", "Y", "X")) |> 
  mutate(hl_reference = "nwmlake",
         hl_id = NA, 
         hf_source = "reference_features",
         hl_source = "nwm_v2.1.6_LAKEPARAM")

nwmlake_pois <- distinct(select(rl, hl_link, hf_id, hl_reference))  %>% 
  filter(hl_reference == 'WBOut_rl') |>
  select(-hl_reference) |>
  right_join(lk, by = "hl_link") |>
  relocate(all_of(schema))

# NWM Reservoirs  --------------------------------------------------------------

f <- list.files(glue('{base}/nwm_v3.0.7/conus'), 
                pattern = "reservoir", 
                full.names = TRUE)

res <- list()

for(i in 1:length(f)){
  nc <- open.nc(f[i])
  
  type <- paste0("nwm_v3.0.7_", gsub(".nc", "", basename(f[i])))
  
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
           hl_id = NA) 
  
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
           hl_id = NA) %>% 
    relocate()
  
  res[[i]] <- bind_rows(usgs_lake, rfc_lake, usace_lake)
}

nwm_res <- relocate(bind_rows(res), all_of(schema))

# Calibration Gages -------------------------------------------------------

cal <- read.csv(glue('{base}/hydrolocations/gage_list.txt')) %>% 
  mutate(site_no = stringi::stri_pad_left(site_no, 8, "0") )

tmp <- filter(community_pois, hl_reference == "Gages") %>% 
  filter(hl_link %in% cal$site_no)

stopifnot(nrow(tmp) == nrow(cal))

calib_poi <- tmp |>
  mutate(hl_id = NA,
         hl_source = "noaaowp",
         hl_reference = "nwm_v3.0.7_calibration_gage") |>
  relocate(all_of(schema))


# Coastal Gages -----------------------------------------------------------

coastal_network <- read.csv(glue('{base}/hydrolocations/GAGE_SUMMARY.csv')) |>
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


coastal_gage_poi <- relocate(coastal_network, all_of(schema))

# Coastal Domain ----------------------------------------------------------

coastal_domain <- read_sf(glue('{base}/hydrolocations/coastal_domain.gpkg')) |>
  st_transform(5070) |>
  st_cast("MULTILINESTRING")

m <-  st_intersects(st_transform(vpu_boundaries[1:21,], 5070), coastal_domain)

df <-  data.frame(gpkgs = list.files(ref_base, pattern = ".gpkg", full.names = TRUE)) |>
  mutate(vpu = sapply(strsplit(basename(gpkgs), "_"), `[[`, 1)) |> 
  filter(vpu %in% vpu_boundaries[1:21,]$VPUID[lengths(m) > 0])

domain_pois <- list()

for(i in 1:nrow(df)){
  
  fl_tmp <- read_sf(df$gpkgs[i], "flowlines")
  
  d <- st_intersection(coastal_domain, AOI::bbox_get(fl_tmp))
  
  touches <- st_intersects(fl_tmp, d)
  
  fl_map <- fl_tmp[lengths(touches) > 0, ]
  
  domain_pois[[i]] <- st_intersection(coastal_domain, fl_map) |>
    select(hf_id = COMID, domain_id) |>
    st_collection_extract("POINT") |>
    st_cast("POINT") |>
    distinct() |>
    group_by(domain_id) |>
    mutate(hl_link = paste(domain_id, df$vpu[i], 1:n(), sep = "-")) |>
    ungroup() |>
    mutate(
      hl_id = NA,
      hl_source = "coastal",
      hl_reference = domain_id,
      hf_source = "reference_features",
      domain_id = NULL,
    ) |> 
    st_transform(4326) |>
    mutate(X = st_coordinates(.)[,1],
           Y = st_coordinates(.)[,2]) |> 
    st_drop_geometry()
  
}

coastal_domain_pois <-  relocate(bind_rows(domain_pois), all_of(schema))

# AHPS --------------------------------------------------------------------

ahps <- read_sf(glue('{base}/hydrolocations/nws_lid.gpkg')) |>
  mutate(hf_id = as.numeric(nwm_feature_id)) |>
  select(nws_lid, usgs_site_code, hf_id) |>
  pivot_longer(-c(hf_id, geom), names_to = "hl_reference", values_to = "hl_link") |>
  filter(!is.na(hl_link)) |>
  mutate(
    hl_id = NA,
    hl_source = "fim",
    hf_source = "reference_features",
  ) |> 
  st_transform(4326) |>
  mutate(X = st_coordinates(.)[,1],
         Y = st_coordinates(.)[,2]) |> 
  st_drop_geometry() |>
  relocate(all_of(schema))

# Final  ------------------------------------------------------------------

hl = rbindlist(list(community_pois,
                    rl_pois,
                    nwmlake_pois,
                    nwm_res,
                    calib_poi,
                    coastal_gage_poi,
                    coastal_domain_pois,
                    ahps
)) |>
  group_by(X, Y) |> 
  mutate(poi_id = cur_group_id(),
         count = n(),
         hl_id = NULL) %>% 
  ungroup() 

unlink(full_hl)

# Write Table
write_sf(select(hl, -count), full_hl, "properties")

# Write spatial location: POIs
select(hl, poi_id, X, Y, count) |> 
  distinct() |> 
  st_as_sf(coords = c("X", "Y"), 
           remove = FALSE,  
           crs = 4326) |>
  write_sf(full_hl, "pois")

# Write spatial location: HL's
st_as_sf(
  select(hl, -count),
  coords = c("X", "Y"),
  remove = FALSE,
  crs = 4326) |>
  write_sf(full_hl, "hydrolocations")

table(hl$hl_reference)
