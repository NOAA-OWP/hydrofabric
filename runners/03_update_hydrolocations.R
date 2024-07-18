source("runners/config.R")
devtools::load_all('.')

schema


# Hydrolocations ----------------------------------------------------------

# Purpose --> This code constructs a conus_hl (conus hydrolocation) file from a set of inputs,
# https://github.com/NOAA-OWP/hydrofabric/issues/40

# Hydrolocations have the following schema
# > poi_id: Unique Identifier to aggregate same physical POIs

# > hl_source: Source of hl
# > hl_reference: Type of Hydrolocation
# > hl_link: Identifier from source dataset
# > hl_position: is the hydrolocation on the inlet or outlet of the hf_id
#
# > X: Longitude (EPSG: 5070)
# > Y: Latitude  (EPSG: 5070)

# > hf_id: ID from reference hydrofabric
# > hf_source: Source of reference hydrofabric

# Reference Fabric --------------------------------------------------------

s <- '/Users/mikejohnson/hydrofabric/v2.2/reference'

hl_inventory <- '/Users/mikejohnson/hydrofabric/hydrolocation-inventory'

system(
  glue("aws s3 sync s3://lynker-spatial/hydrofabric/hydrolocation-inventory/ {hl_inventory}")
)

f <- list.files(hl_inventory, full.names = TRUE)

ref_net <- collect(open_dataset(glue('{s}/conus_network')))

# 1. GFv20 Reference -------------------------------------------------------------------

gfv20 <- open_dataset(glue("{s}/gfv2_hydrolocations")) |>
  collect() %>%
  mutate(hf_id = id,
         hf_source = "reference_features",
         hl_source = "GFv20") %>%
  distinct() %>%
  mutate(hl_position = case_when(hl_reference == 'XXX' ~ "inlet", TRUE ~ "outlet")) |>
  select(all_of(schema)) |>
  relocate(all_of(schema))

# 2. Calibration Gages -------------------------------------------------------

cal <- grep("gage_list", f, value = TRUE) |> 
  read.csv()  |> 
  mutate(site_no = paste0("0", site_no))

tmp <- filter(gfv20, hl_reference == "Gages", hl_link %in% cal$site_no)

## --
stopifnot(nrow(tmp) == nrow(cal))
## --

calib_poi <- mutate(tmp,
                    hl_id = NA,
                    hl_source = "noaaowp",
                    hl_reference = "nwm_v3.0.7_calibration_gage") |>
  select(all_of(schema)) %>%
  relocate(all_of(schema)) %>%
  distinct()

# 3. lid --------------------------------------------------------------------

lid <- read_sf(grep("nws_lid", f, value = TRUE)) |>
  mutate(hf_id = as.numeric(nwm_feature_id)) |>
  select(nws_lid, usgs_site_code, hf_id) |>
  pivot_longer(-c(hf_id, geom), names_to = "hl_reference", values_to = "hl_link") |>
  filter(!is.na(hl_link)) |>
  mutate(hl_source = "fim",
         hf_source = "reference_features",
         hl_position = "outlet",)  %>%
  mutate(X = st_coordinates(.)[, 1], Y = st_coordinates(.)[, 2]) |>
  st_drop_geometry() |>
  distinct() |>
  left_join(select(gfv20, poi_id, hl_link), by = "hl_link") |>
  select(all_of(schema)) %>%
  relocate(all_of(schema))

# 4. RouteLink ---------------------------------------------------------------
nc <- open.nc(grep("RouteLink", f, value = TRUE))

var <- c('link', 'NHDWaterbodyComID')

rl_pois <- lapply( 1:length(var), FUN = function(x) { var.get.nc(nc, var[x])} ) |>
  bind_cols(.name_repair = "unique_quiet") |>
  setNames(c("hf_id", "hl_link")) |>
  filter(hl_link > 0) |>
  distinct() |>
  left_join(distinct(select(ref_net, hf_id, toid, hydroseq)), by = "hf_id") |>
  filter(!is.na(hydroseq)) |>
  tidyr::drop_na() |>
  group_by(hl_link) |>
  arrange(hydroseq) |>
  mutate(hl_link = as.character(hl_link)) |>
  distinct() %>%
  left_join(distinct(select(
    ref_net, hf_id, X = outlet_X, Y = outlet_Y
  )))

rl_outlets = filter(rl_pois, hydroseq == min(hydroseq)) |>
  ungroup()  |>
  select(hl_link, hf_id) |>
  mutate(
    hf_source = "reference_features",
    hl_source = "nwm_v3.0.9_routelink",
    hl_reference = "WBOut",
    hl_position = "outlet",
    poi_id = NA
  )

rl_inlets = filter(rl_pois, !hf_id %in% toid) |>
  ungroup() %>%
  left_join(select(ref_net, fromid = id, hf_id = toid), by = "hf_id") |>
  select(hl_link, hf_id = fromid)  |>
  filter(!is.na(hf_id)) %>%
  mutate(
    hf_source = "reference_features",
    hl_source = "nwm_v3.0.9_routelink",
    hl_reference = "WBIn",
    hl_position = "outlet",
    poi_id = NA
  )

rl_wbs = bind_rows(rl_outlets, rl_inlets) %>%
  left_join(distinct(select(
    ref_net, hf_id, X = outlet_X, Y = outlet_Y
  ))) %>%
  select(all_of(schema)) %>%
  relocate(all_of(schema)) %>%
  distinct()

# 5. NWM Lakes  --------------------------------------------------------------

nc  <- open.nc(grep("LAKEPARM", f, value = TRUE))

var <- c("lake_id")

lk <- lapply(
  1:length(var),
  FUN = function(x) {
    var.get.nc(nc, var[x])
  }
) %>%
  bind_cols(.name_repair = "unique_quiet") |>
  setNames(c("hl_link")) |>
  mutate(
    hl_reference = "nwmlake",
    hl_link = as.character(hl_link),
    poi_id = NA,
    hf_source = "reference_features",
    hl_source = "nwm_v2.1.6_LAKEPARAM"
  )

nwmlake_pois <- distinct(select(rl_wbs, hl_link, hf_id, hl_reference, hl_position, X, Y))  %>%
  # Pin these lakes to the same outlet in the RL
  filter(hl_reference == 'WBOut') |>
  select(-hl_reference) |>
  inner_join(lk, by = "hl_link") %>%
  #lakes that are present here but not in RL are Canadian
  tidyr::drop_na(X) %>%
  select(all_of(schema)) %>%
  relocate(all_of(schema)) %>%
  distinct()

# 6. NWM Reservoirs  --------------------------------------------------------------

res_file <- grep("reservoir", f, value = TRUE)

res <- list()

master = distinct(select(
  filter(rl_wbs, hl_reference == "WBOut"),
  hl_link,
  hl_position,
  X,
  Y,
  starts_with("hf_")
))

for (i in 1:length(res_file)) {
  nc <- open.nc(res_file[i])
  
  type <- paste0("nwm_v3.0.9_", gsub(".nc", "", basename(res_file[i])))
  
  usgs_var <- c("usgs_lake_id", 'usgs_gage_id')
  
  usgs_lake <- lapply(
    1:length(usgs_var),
    FUN = function(x) {
      var.get.nc(nc, usgs_var[x])
    }
  ) %>%
    bind_cols(.name_repair = "unique_quiet") %>%
    setNames(c("hl_link", "hl_id")) %>%
    mutate(hl_source = type, hl_reference = "usgs_gage") %>%
    mutate(hl_link = as.character(hl_link)) |> 
    left_join(master, by = "hl_link") %>%
    mutate(hl_link = hl_id, poi_id = NA)
  
  rfc_var <- c("rfc_lake_id", 'rfc_gage_id')
  
  rfc_lake <- lapply(
    1:length(rfc_var),
    FUN = function(x) {
      var.get.nc(nc, rfc_var[x])
    }
  ) %>%
    bind_cols(.name_repair = "unique_quiet") %>%
    setNames(c("hl_link", "hl_id")) %>%
    mutate(hl_source = type, hl_reference = "rfc_gage") %>%
    mutate(hl_link = as.character(hl_link)) |> 
    left_join(master, by = "hl_link") %>%
    mutate(hl_link = hl_id, hl_id = NA)
  
  usace_var <- c("usace_lake_id", 'usace_gage_id')
  
  usace_lake <- lapply(
    1:length(usace_var),
    FUN = function(x) {
      var.get.nc(nc, usace_var[x])
    }
  ) %>%
    bind_cols(.name_repair = "unique_quiet") %>%
    setNames(c("hl_link", "hl_id")) %>%
    mutate(hl_source = type, hl_reference = "usace_gage") %>%
    mutate(hl_link = as.character(hl_link)) |> 
    left_join(master, by = "hl_link") %>%
    mutate(hl_link = hl_id, poi_id = NA) %>%
    select(all_of(schema)) %>%
    distinct()
  
  
  res[[i]] <- bind_rows(usgs_lake, rfc_lake, usace_lake)
}

nwm_res <- bind_rows(res) %>%
  select(all_of(schema)) %>%
  relocate(all_of(schema)) %>%
  distinct()

# 7. Coastal Gages -----------------------------------------------------------

coastal_gages <- read.csv(grep("GAGE_SUMMARY", f, value = TRUE)) |>
  mutate(
    hl_link = sprintf('%08s', SITE_NO),
    hf_id = NA,
    hl_id = NA,
    hl_source = "coastal",
    hl_reference = "coastal_gage",
    hl_position = "outlet",
    hf_source = "reference_features"
  ) |>
  select(starts_with(c("hl_", "hf_")), Y = LAT_NHD, X = LON_NHD)

coastal_domain <- read.csv(grep("coastal", f, value = TRUE)) |>
  mutate(
    hl_link = hl_link,
    hf_id = NA,
    hl_id = NA,
    hl_source = "coastal",
    hl_reference = "coastal_domain",
    hl_position = "outlet",
    hf_source = "reference_features"
  ) |>
  select(starts_with(c("hl_", "hf_")), Y = lat, X = long)

coastal = bind_rows(coastal_gages, coastal_domain)

for (i in 1:nrow(coastal)) {
  coastal$hf_id[i] <- 
    as.numeric(findNLDI(location = c( coastal$X[i], coastal$Y[i]))$origin$comid)
}

coastal_poi <- coastal_network %>%
  st_as_sf(coords = c("X", "Y"), crs = 4326) %>%
  st_transform(5070) %>%
  mutate(X = st_coordinates(.)[, 1],
         Y =  st_coordinates(.)[, 2],
         poi_id = NA) %>%
  st_drop_geometry() %>%
  select(all_of(schema)) %>%
  relocate(all_of(schema)) %>%
  distinct()

# Coastal Domain ----------------------------------------------------------

# coastal_domain <- read_sf(grep("coastal_domain", f, value = TRUE)) |>
#   st_transform(5070) |>
#   st_cast("MULTILINESTRING")
#
# ref_fl     <- glue('{s}/conus_flowlines')
#
# m <-  st_intersects(st_transform(vpu_boundaries[1:21,], 5070), coastal_domain)
#
# v = vpu_boundaries[1:21,]$VPUID[lengths(m) > 0]
#
# domain_pois <- list()
#
# for(i in 1:length(v)){
#
#   fl_tmp <-  open_dataset(ref_fl) %>%
#     dplyr::filter(vpuid == !!v[i]) %>%
#     read_sf_dataset()
#
#   d <- st_intersection(coastal_domain, AOI::bbox_get(fl_tmp))
#
#   touches <- st_intersects(fl_tmp, d)
#
#   fl_map <- fl_tmp[lengths(touches) > 0, ]
#
#   domain_pois[[i]] <- st_intersection(coastal_domain, fl_map) |>
#     select(hf_id = id, domain_id) |>
#     st_collection_extract("POINT") |>
#     st_cast("POINT") |>
#     distinct() |>
#     group_by(hf_id) |>
#     mutate(hl_link = paste(domain_id[1], hf_id, sep = "-")) |>
#     slice(1) %>%
#     ungroup() |>
#     mutate(
#       poi_id = NA,
#       hl_source = "coastal",
#       hl_reference = domain_id,
#       hl_position = "outlet",
#       hf_source = "reference_features",
#       domain_id = NULL,
#     ) %>%
#     mutate(X = st_coordinates(.)[,1],
#            Y = st_coordinates(.)[,2]) |>
#     st_drop_geometry()
#
#   message("VPU", v[i], " (", i, "/", length(v), ")")
# }
#
# coastal_domain_pois <- bind_rows(domain_pois) %>%
#   select(all_of(schema)) %>%
#   relocate(all_of(schema)) %>%
#   distinct()

# Final  ------------------------------------------------------------------
hl_fin = rbindlist(
  list(
    rl_wbs,
    nwmlake_pois,
    nwm_res,
    calib_poi,
    coastal_poi,
    lid,
    gfv20
  )
) |> 
  group_by(X, Y, hf_id) |>
  arrange(-poi_id) |>
  mutate(poi_id = zoo::na.locf(poi_id, na.rm=FALSE),
         gid = cur_group_id()) |> 
  ungroup()

gid  = unique(hl_fin$gid)
nids = unique(gid[!gid %in% unique(hl_fin$poi_id)])

hl_fin_complete = filter(hl_fin, is.na(poi_id)) |>
  group_by(gid) %>%
  mutate(
    gid = cur_group_id(),
    poi_id = nids[gid]) %>%
  ungroup() |> 
  bind_rows(filter(hl_fin, !is.na(poi_id))) %>%
  group_by(poi_id) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  left_join(distinct(select(ref_net, hf_id, mainstemlp, vpuid)), by = "hf_id", relationship = "many-to-many") %>%
  distinct() %>%
  mutate(hl_uri = paste0(hl_reference, "-", hl_link)) 

# Write Hydrolocation Table
hl_fin_complete %>%
  group_by(vpuid) %>%
  arrow::write_dataset(glue("{s}/conus_hydrolocations"), version = 2.6)

system(glue("aws s3 sync {s}/conus_hydrolocations s3://lynker-spatial/hydrofabric/v2.2/reference/conus_hydrolocations"))

