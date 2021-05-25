source('workflow/utils.R')
#############################################################################
path       <- "workflow/nhd_workflows/cache/ngen_01a-3.gpkg"
out_path   <- "workflow/nhd_workflows/cache/ngen_01a-4.gpkg"
min_area   <- 3
max_area   <- 15
min_length <- 0.6

cat = read_sf(path, "catchments")

fl  = read_sf(path, "flowpaths") %>%
  left_join(st_drop_geometry(cat), by = "comid") %>%
  left_join(select(st_drop_geometry(.), to_lp = levelpathi, tocomid = comid)) %>%
  mutate(inlet = !comid %in% tocomid)

shorts = filter(fl, (levelpathi != to_lp) & inlet) %>%
  filter(areasqkm <= min_area | lengthkm <= min_length) %>%
  arrange(areasqkm)

cat_dt = mutate(data.table(cat), ind = 1:n())
fl_dt  = mutate(data.table(fl),  ind = 1:n())

removed_ids = list()

for (i in 1:nrow(shorts)) {
  # Merge Catchments into their tocomid ...
  short_cat = filter(cat_dt,  comid == shorts$comid[i])
  int_cat   = filter(cat_dt,  comid %in% shorts$tocomid[i])

  if((short_cat$areasqkm + int_cat$areasqkm) <= max_area){
    # Aggregate catchment geometries
    if(nrow(short_cat) != 0 & nrow(int_cat) != 0 ){
      new_geom = st_union(short_cat$geom, int_cat$geom)
      cat_dt[int_cat$ind, geom := st_as_sf(new_geom)]
    }

    # Aggregate flowlines, member COMIDs, and member IDs
    short_fl = filter(fl_dt,  comid == shorts$comid[i])
    int_fl   = filter(fl_dt,  comid == shorts$tocomid[i])

    # Remove flowlines from network
    if(nrow(short_fl) != 0 & nrow(int_fl) != 0 ){
      #new_geom = st_union(short_fl$geom, int_fl$geom)
      #fl_dt[int_fl$ind, geom := st_as_sf(new_geom)]
      fl_dt[int_fl$ind, member_COMID := paste0(short_fl$member_COMID, int_fl$member_COMID)]
    }

    # Collect the IDs to remove...
    removed_ids[[i]] = shorts$comid[i]
  }
}

removed_ids = unlist(removed_ids)

new_cat = filter(st_as_sf(cat_dt), !comid %in% removed_ids) %>%
  mutate(areasqkm = as.numeric(st_area(.)/1e6)) %>%
  select(comid, areasqkm)

new_fl =  filter(st_as_sf(fl_dt), !comid %in% removed_ids) %>%
  mutate(lengthkm = as.numeric(st_length(.)/1e3)) %>%
  select(-ind, -areasqkm) %>%
  left_join(st_drop_geometry(new_cat), by = "comid")

#################################################
#################################################

shorts = new_fl %>%
  filter(inlet, areasqkm <= 1 | lengthkm <= min_length)

cat_dt2 = mutate(data.table(new_cat), ind = 1:n())
fl_dt2  = mutate(data.table(new_fl),  ind = 1:n())

removed_ids = list()

for (i in 1:nrow(shorts)) {
  # Merge Catchments into their tocomid ...
  short_cat = filter(cat_dt2,  comid == shorts$comid[i])
  int_cat   = filter(cat_dt2,  comid %in% shorts$tocomid[i])

  # Aggregate catchment geometries
    if(nrow(short_cat) != 0 & nrow(int_cat) != 0 ){
      new_geom = st_union(short_cat$geom, int_cat$geom)
      cat_dt2[int_cat$ind, geom := st_as_sf(new_geom)]
    }

    # Aggregate flowlines, member COMIDs, and member IDs
    short_fl = filter(fl_dt2,  comid == shorts$comid[i])
    int_fl   = filter(fl_dt2,  comid == shorts$tocomid[i])

    # Remove flowlines from network
    if(nrow(short_fl) != 0 & nrow(int_fl) != 0 ){
      #new_geom = st_union(short_fl$geom, int_fl$geom)
      #fl_dt[int_fl$ind, geom := st_as_sf(new_geom)]
      fl_dt2[int_fl$ind, member_COMID := paste0(short_fl$member_COMID, int_fl$member_COMID)]
    }

    # Collect the IDs to remove...
    removed_ids[[i]] = shorts$comid[i]
  }


removed_ids = unlist(removed_ids)

new_cat2 = filter(st_as_sf(cat_dt2), !comid %in% removed_ids) %>%
  mutate(areasqkm = as.numeric(st_area(.)/1e6)) %>%
  select(comid, areasqkm)

new_fl2 = filter(st_as_sf(fl_dt2), !comid %in% removed_ids) %>%
  mutate(lengthkm = as.numeric(st_length(.)/1e3)) %>%
  select(-ind, -areasqkm) %>%
  left_join(st_drop_geometry(new_cat2), by = "comid")


mutate(lengthkm = as.numeric(st_length(.)/1e3))

#################################################

make_plot(new_fl3, new_cat3, "Headwater Aggregation HF") %>%
  ggsave(filename = "workflow/nhd_workflows/cache/img/04-ngen-firstorder-dissolve.png",
         units = "in", height = 4, width = 8)

message("Dropped: ", nrow(cat) - nrow(new_cat2), " features")

unlink(out_path)
write_sf(new_cat2, out_path, "catchments")
write_sf(new_fl2, out_path, "flowpaths")

