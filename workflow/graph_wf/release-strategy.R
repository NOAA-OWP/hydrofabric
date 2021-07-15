# -------------------------------------------------------------------------
data = readRDS("data/sample_network.rds")
dir = "data"; 
version = "beta"
# -------------------------------------------------------------------------

cat = data$cat
fl = data$flowpath

net = build_node_net(fl)

nexus_prefix     = "nex-"
catchment_prefix = "cat-"
waterbody_prefix = "wb-"


dir.create(file.path(dir, version), showWarnings = FALSE)

cfile <- file.path(dir, version,  "catchment_data.geojson")
wfile <- file.path(dir, version,  "flowpath_data.geojson")
nfile <- file.path(dir, version, "nexus_data.geojson")

cfe <- file.path(dir, version, "catchment_edge_list.json")
wfe <- file.path(dir, version, "waterbody_edge_list.json")
cwe <- file.path(dir, version, "crosswalk-mapping.json")
 
gpkg <- file.path(dir, version, paste0('ngen-', version,  ".gpkg"))

fline = net$network %>%
  select(ID, from, to, lengthkm, member_comid, levelpath, toType, fromType) 

nexus <-   net$node

catchment_edge_list = st_drop_geometry(fline) %>% 
  mutate(id = paste0(catchment_prefix, ID), 
         #to = paste0(ifelse(toType == "nex", "", paste0(toType, "-")), nexus_prefix, to),
         #from = paste0(ifelse(fromType == "nex", "", paste0(fromType, "-")), nexus_prefix, from)) %>% 
         to = paste0(nexus_prefix, to),
         from = paste0(nexus_prefix, from)) %>%
  select(id, to, from) 

waterbody_edge_list = st_drop_geometry(fline) %>% 
  mutate(id = paste0(waterbody_prefix, ID), 
         to = paste0(nexus_prefix, to),
         from = paste0(nexus_prefix, from)) %>%
         # to = paste0(ifelse(toType == "nex", "", paste0(toType, "-")), nexus_prefix, to),
         # from = paste0(ifelse(fromType == "nex", "", paste0(fromType, "-")), nexus_prefix, from)) %>% 
  select(id, to, from) 

names(catchment_edge_list) <- tolower(names(catchment_edge_list))
names(waterbody_edge_list) <- tolower(names(waterbody_edge_list))

###############################################################

catchment_data = cat %>%
  mutate(areasqkm = as.numeric(set_units(st_area(.), "km2"))) %>% 
  select(ID, areasqkm) %>%
  mutate(ID = paste0(catchment_prefix, ID)) %>%
  left_join(catchment_edge_list, by = c("ID" = "id"))

flowpath_data = select(fline, ID, lengthkm, main_id = levelpath) %>%
  mutate(ID = paste0(catchment_prefix, ID)) %>%
  left_join(catchment_edge_list, by = c("ID" = "id"))

nexus_data =  mutate(nexus, 
                     ID = paste0(nexus_prefix, nexID)) %>% 
  select(ID)

###############################################################

usgs = readRDS("data/usgs_rc.rds") %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  select(COMID, siteID, name) %>%
  st_transform(5070)

oo = st_filter(usgs, catchment_data)

what_nwis_data <- dataRetrieval::whatNWISdata(siteNumber = gsub("USGS-", "", oo$siteID))

nwis_sites <- filter(what_nwis_data, parm_cd == "00060" & data_type_cd == "uv") %>%
  st_as_sf(coords = c("dec_long_va", "dec_lat_va"),
           crs = 4269) %>%
  st_transform(5070) %>% 
  st_join(cat) %>% 
  select(site_no, ID)

oo = left_join(select(st_drop_geometry(fl), ID, member_comid), st_drop_geometry(nwis_sites))

nhd_crosswalk <- lapply(1:nrow(oo), function(x) {
  list(COMID = strsplit(oo$member_comid[x], ",")[[1]],
       site_no = jsonlite::unbox(oo$site_no[x]))
})

names(nhd_crosswalk) = paste0(catchment_prefix, oo$ID)

for(i in 1:length(nhd_crosswalk)){
  if(is.na(nhd_crosswalk[[i]]$site_no)){
    nhd_crosswalk[[i]]$site_no <- NULL
  }
}

#######################################
write_geojson <- function(x, y) {
  names(x) <- tolower(names(x))
  unlink(y)
  write_sf(st_make_valid(st_transform(x, 4326)), y, 
           layer_options = c("ID_FIELD=id", "ID_TYPE=String"))
}

write_geojson(catchment_data, cfile)
write_geojson(flowpath_data, wfile)
write_geojson(x = nexus_data, nfile)

jsonlite::write_json(catchment_edge_list, cfe, pretty = TRUE)
jsonlite::write_json(waterbody_edge_list, wfe, pretty = TRUE)
jsonlite::write_json(nhd_crosswalk, cwe, pretty = TRUE, auto_unbox = FALSE)

write_sf(catchment_data, gpkg, 'cat')
write_sf(flowpath_data,  gpkg, 'fl')
write_sf(nexus_data, gpkg, 'nex')

zip(zipfile = file.path(dir, version, paste0('ngen-',version, ".zip")), 
    files   = dir(file.path(dir, version), full.names = TRUE))




