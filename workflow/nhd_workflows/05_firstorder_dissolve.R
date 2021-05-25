source('workflow/utils.R')
#############################################################################
path       <- "workflow/nhd_workflows/cache/ngen_01a-3.gpkg"
out_path   <- "workflow/nhd_workflows/cache/ngen_01a-4.gpkg"
min_area   <- 3
max_area   <- 15
min_length <- 0.6

cat = read_sf(path, "catchments")
fl  = read_sf(path, "flowpaths")

fl$inlet = !fl$comid %in%fl$tocomid
s = nhdplusTools::get_node(fl, "start")
inlets = s[!fl$comid %in%fl$tocomid,]
write_sf(inlets, out_path, "inlets")

cat = left_join(cat, st_drop_geometry(select(fl, comid, tocomid)), by = "comid")

#############################################################################
# merge short
shorts = filter(fl, order == 1, lengthkm < min_length)

cat_dt = mutate(data.table(cat), ind = 1:n())
fl_dt  = mutate(data.table(fl),  ind = 1:n())

removed_ids = list()

for (i in 1:nrow(shorts)) {
  # Merge Catchments into their toID ...
  short_cat = filter(cat_dt,  comid == shorts$comid[i])
  int_cat   = filter(cat_dt,  comid %in% shorts$tocomid[i])
  # Aggregate catchment geometries
  if(nrow(short_cat) != 0 & nrow(int_cat) != 0 ){
    new_geom = st_union(short_cat$geom, int_cat$geom)
    cat_dt[int_cat$ind, geom := st_as_sf(new_geom)]
  }

  # Aggregate flowlines, member COMIDs, and member IDs
  short_fl = filter(fl_dt,  comid == shorts$comid[i])
  int_fl   = filter(fl_dt,  comid == shorts$tocomid[i])

  if(nrow(short_fl) != 0 & nrow(int_fl) != 0 ){
    new_geom = st_union(short_fl$geom, int_fl$geom)
    fl_dt[int_fl$ind, geom := st_as_sf(new_geom)]
    fl_dt[int_fl$ind, member_COMID := paste0(short_fl$member_COMID, int_fl$member_COMID)]
  }

  # Collect the IDs to remove...
  removed_ids[[i]] = shorts$comid[i]
}

removed_ids = unlist(removed_ids)

new_cat = filter(st_as_sf(cat_dt), !comid %in% removed_ids) %>%
  mutate(areasqkm = as.numeric(st_area(.)/1e6)) %>%
  select(comid, areasqkm)

new_fl =  filter(st_as_sf(fl_dt), !comid %in% removed_ids) %>%
  mutate(lengthkm = as.numeric(st_length(.)/1e3)) %>%
  select(-ind)

new_cat = left_join(new_cat, st_drop_geometry(select(fl, comid, order, tocomid)), by = "comid")

#################################################
write_sf(new_cat, "workflow/cache/tester.gpkg", "cats")
write_sf(new_fl, "workflow/cache/tester.gpkg", "fl")
#################################################

stouts = filter(new_cat, areasqkm <= min_area) %>%
  filter(order == 1)  %>%
  mutate(int_area = new_cat$areasqkm[match(tocomid, new_cat$comid)],
         merge_area = areasqkm + int_area) %>%
  filter(merge_area <= max_area | areasqkm <= 1)

cat_dt2 = mutate(data.table(new_cat), ind = 1:n())

removed_ids2 = list()

for (i in 1:nrow(stouts)) {
  # Merge Catchments into their toID ...
  short_cat = filter(cat_dt2,  comid == stouts$comid[i])
  int_cat   = filter(cat_dt2,  comid == stouts$tocomid[i])

  # Aggregate catchment geometries
  if(nrow(short_cat) != 0 & nrow(int_cat) != 0 ){
    new_geom = st_union(short_cat$geom, int_cat$geom)
    cat_dt2[int_cat$ind, geom := st_as_sf(new_geom)]
  }

  # Collect the IDs to remove...
  removed_ids2[[i]] = stouts$comid[i]
}

removed_ids2 = unlist(removed_ids2)

new_cat2 = filter(st_as_sf(cat_dt2), !comid %in% removed_ids2) %>%
  mutate(areasqkm = as.numeric(st_area(.)/1e6)) %>%
  select(comid, areasqkm, order, tocomid)

new_fl2 =  filter(new_fl, !comid %in% removed_ids2) %>%
  mutate(lengthkm = as.numeric(st_length(.)/1e3))

stouts = filter(new_cat2, areasqkm <= min_area) %>%
  filter(order == 1)  %>%
  mutate(int_area = new_cat2$areasqkm[match(tocomid, new_cat2$comid)],
         merge_area = areasqkm + int_area) %>%
  filter(merge_area <= max_area | areasqkm <= 1)

cat_dt3 = mutate(data.table(new_cat2), ind = 1:n())

removed_ids3 = list()

for (i in 1:nrow(stouts)) {
  # Merge Catchments into their toID ...
  short_cat = filter(cat_dt3,  comid == stouts$comid[i])
  int_cat   = filter(cat_dt3,  comid == stouts$tocomid[i])

  # Aggregate catchment geometries
  if(nrow(short_cat) != 0 & nrow(int_cat) != 0 ){
    new_geom = st_union(short_cat$geom, int_cat$geom)
    cat_dt3[int_cat$ind, geom := st_as_sf(new_geom)]
  }

  # Collect the IDs to remove...
  removed_ids3[[i]] = stouts$comid[i]
}

removed_ids3 = unlist(removed_ids3)

new_cat3 = filter(st_as_sf(cat_dt3), !comid %in% removed_ids3) %>%
  mutate(areasqkm = as.numeric(st_area(.)/1e6)) %>%
  select(comid, areasqkm, order, tocomid)

new_fl3 =  filter(new_fl2, !comid %in% removed_ids3) %>%
  mutate(lengthkm = as.numeric(st_length(.)/1e3))

#################################################
make_plot(new_fl3, new_cat3, "First Order Dissolve HF") %>%
  ggsave(filename = "workflow/nhd_workflows/cache/img/04-ngen-firstorder-dissolve.png",
         units = "in", height = 4, width = 8)

message("Dropped: ", nrow(cat) - nrow(new_cat2), " features")

unlink(out_path)
write_sf(new_cat3, out_path, "catchments")
write_sf(new_fl3, out_path, "flowpaths")

