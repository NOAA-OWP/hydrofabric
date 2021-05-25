source('workflow/utils.R')
#############################################################################
path       <- "workflow/cache/ngen_01a-2.gpkg"
out_path   <- "workflow/cache/ngen_01a-3.gpkg"
#############################################################################

cat = read_sf(path, "catchments")

fl  = read_sf(path, "flowpaths") %>%
  left_join(st_drop_geometry(cat))

#############################################################################

shorts = filter(fl, stream_order == 1) %>%
  filter(lengthkm < .6)

filter(shorts, toID == 3526)

ll = st_intersects(shorts, fl)

mapping  <- data.frame(shortID = rep(shorts$ID, sapply(ll, length)),
                       toID    = rep(shorts$toID,sapply(ll, length)),
                       intID   = fl$ID[unlist(ll)],
                       so      = fl$stream_order[unlist(ll)])  %>%
  mutate(toID_match = ifelse(toID == intID, TRUE,FALSE)) %>%
  filter(shortID != intID)

m = mapping %>%
  group_by(shortID) %>%
  arrange(-toID_match, -so) %>%
  slice(1) %>%
  ungroup()

#################################################
cat_dt = mutate(data.table(cat), ind = 1:n())

for (i in 1:nrow(m)) {
  short_fl = filter(cat_dt,  ID == m$shortID[i])
  int_flow = filter(cat_dt,  ID == m$intID[i])
  new = st_union(short_fl$geom, int_flow$geom)
  cat_dt[int_flow$ind, geom := st_as_sf(new)]
}

new_cat = filter(st_as_sf(cat_dt), !ID %in% m$shortID) %>%
  mutate(areasqkm = as.numeric(st_area(.)/1e6))

######
#Drop FL but merge member COMID into the merged one!

fl_dt = mutate(data.table(fl), ind = 1:n())

for (i in 1:nrow(m)) {
  short_fl = filter(fl_dt,  ID == m$shortID[i])
  int_flow = filter(fl_dt,  ID == m$intID[i])
  new = st_line_merge(st_union(short_fl$geom, int_flow$geom))
  fl_dt[int_flow$ind, member_COMID := paste0(short_fl$member_COMID, int_flow$member_COMID)]
  #fl_dt[int_flow$ind, geom := st_as_sf(new)]
}

new_fl = st_as_sf(fl_dt) %>%
  filter(!ID %in% m$shortID) %>%
  flowpaths_to_linestrings() %>%
  mutate(lengthkm = as.numeric(st_length(.)/1e3))
#################################################
#################################################
write_sf(new_cat, "workflow/cache/ngen_tmp.gpkg", "catchments")
write_sf(new_fl, "workflow/cache/ngen_tmp.gpkg", "flowpaths")
#################################################
#################################################

stout_fl = new_fl %>%
  left_join(select(st_drop_geometry(new_cat), ID, areasqkm)) %>%
  filter(areasqkm < 3, stream_order == 1)

ll = st_intersects(stout_fl, new_fl)

mapping  <- data.frame(stoutID = rep(stout_fl$ID, sapply(ll, length)),
                 toID          = rep(stout_fl$toID,sapply(ll, length)),
                 stoutArea     = rep(stout_fl$areasqkm,sapply(ll, length)),
                 intID         = new_fl$ID[unlist(ll)],
                 so            = new_fl$stream_order[unlist(ll)],
                 intArea = new_fl$areasqkm[unlist(ll)])  %>%
  mutate(toID_match = ifelse(toID == intID, TRUE,FALSE),
         new_area = stoutArea + intArea) %>%
  filter(stoutID != intID, new_area  < 15) %>%
  arrange(stoutID)

#  Merge short ID with intID
# if multiple options prioritize toIDmatch = TRUE
# if all are FALSE, prioritize the larger stream order
m = mapping %>%
  group_by(stoutID) %>%
  arrange(-so,-toID_match, -new_area) %>%
  slice(1) %>%
  ungroup()


#################################################

cat_dt = data.table(mutate(new_cat, ind = 1:n()))

m1 = filter(m, intID %in% stoutID)
ints = unique(m1$intID)

for (i in 1:length(ints)) {
  short_fls = filter(cat_dt,  ID == m1$stoutID[i])
  int_flow = filter(cat_dt,   ID == m1$intID[i])
  new = st_union(short_fl$geom, int_flow$geom)
  cat_dt[int_flow$ind, geom := st_as_sf(new)]
}


m2 = m %>% filter(!intID %in% stoutID)

for (i in 1:nrow(m2)) {
  short_fl = filter(cat_dt,  ID == m2$stoutID[i])
  int_flow = filter(cat_dt,  ID == m2$intID[i])
  new = st_union(short_fl$geom, int_flow$geom)
  cat_dt[int_flow$ind, geom := st_as_sf(new)]
}


new_cat_2 = st_as_sf(cat_dt) %>%
  filter(!ID %in% m$stoutID) %>%
  mutate(areasqkm = as.numeric(st_area(.)/1e6)) %>%
  select(-ind)

################################################

fl_dt = data.table(new_fl) %>%
  mutate(ind = 1:n())

for (i in 1:nrow(m)) {
  short_fl = filter(fl_dt,  ID == m$stoutID[i])
  int_flow = filter(fl_dt,  ID == m$intID[i])
  #new = st_line_merge(st_union(short_fl$geom, int_flow$geom))
  fl_dt[int_flow$ind, member_COMID := paste0(short_fl$member_COMID, int_flow$member_COMID)]
}

new_fl_2 = st_as_sf(fl_dt) %>%
  filter(!ID %in% m$stoutID) %>%
  flowpaths_to_linestrings() %>%
  mutate(lengthkm = as.numeric(st_length(.)/1e3))


#################################################

make_plot(new_fl_2, new_cat_2, "First Order Dissolve HF") %>%
  ggsave(filename = "workflow/cache/img/04-ngen-firstorder-dissolve.png",
         units = "in", height = 4, width = 8)

message("Dropped: ", nrow(cat) - nrow(new_cat_2), " features")

if(all(validate_network(new_fl_2, new_cat_2))){
  unlink(out_path)
  write_sf(new_cat_2, out_path, "catchments")
  write_sf(new_fl_2, out_path, "flowpaths")
} else {
  message("Error in: ", basename(rstudioapi::getSourceEditorContext()$path))
}
