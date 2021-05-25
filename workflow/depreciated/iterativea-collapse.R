source('workflow/utils.R')
#############################################################################
path       <- "workflow/cache/ngen_01a-2.gpkg"
out_path   <- "workflow/cache/ngen_01a-3.gpkg"
min_area   <- 3
max_area   <- 15
min_length <- 0.6

## ISSUE ID: 5385

intersection_mapping = function(sub, full){

  ll = st_intersects(sub, full)

  data.frame(ID = rep(sub$ID, sapply(ll, length)),
             toID          = rep(sub$toID,sapply(ll, length)),
             IDarea        = rep(sub$areasqkm,sapply(ll, length)),
             IDso          = rep(sub$stream_order,sapply(ll, length)),
             intID         = full$ID[unlist(ll)],
             intSO          = full$stream_order[unlist(ll)],
             intArea = full$areasqkm[unlist(ll)])  %>%
    mutate(toID_match = ifelse(toID == intID, TRUE,FALSE),
           new_area = IDarea + intArea) %>%
    filter(ID != intID) %>%
    group_by(ID) %>%
    mutate(n = 1:n()) %>%
    ungroup()
}

dissolve_network = function(m1, cat_dt, fl_dt, min_area = 3, max_area = 15, filter_fun = NULL){

  if(nrow(m1) != 0){
    for(i in 1:nrow(m1)) {
      og   = filter(cat_dt, ID == m1$toID[i])
      kill = filter(cat_dt, ID == m1$ID[i])
      new_cat_geom = st_union(st_as_sf(rbindlist(list(og,kill)))) %>% st_make_valid()
      new_member_ids = paste(c(og$member_COMID, kill$member_COMID), collapse = ",")
      fl_dt[og$ind, member_COMID := new_member_ids]
      cat_dt[og$ind, geom := st_as_sf(new_cat_geom)]
    }

    new_cat = filter(st_as_sf(cat_dt), !ID %in% m1$ID) %>%
      mutate(areasqkm = as.numeric(st_area(.)/1e6)) %>%
      select(-ind)

    new_fl =  filter(st_as_sf(fl_dt), !ID %in% m1$ID) %>%
      select(-ind, -areasqkm) %>%
      left_join(st_drop_geometry(new_cat), by = "ID")

    cat_dt = data.table(mutate(new_cat, ind = 1:n()))
    fl_dt  = data.table(mutate(new_fl,  ind = 1:n()))

    m1 = filter_fun(new_fl)
  }

  return(list(m1 = m1, fl_dt = fl_dt, cat_dt = cat_dt ))

}

cat = read_sf(path, "catchments")

fl  = read_sf(path, "flowpaths") %>%
  left_join(st_drop_geometry(cat)) %>%
  flowpaths_to_linestrings()

#############################################################################
# merge short flowlines into the downstream catchment
shorts = filter(fl, stream_order == 1) %>%
  filter(lengthkm < min_length)

m = intersection_mapping(shorts, fl) %>%
  group_by(ID) %>%
  arrange(-toID_match, -intSO) %>%
  slice(1) %>%
  ungroup()

cat_dt = mutate(data.table(cat), ind = 1:n())
fl_dt  = mutate(data.table(fl),  ind = 1:n())

for (i in 1:nrow(m)) {
  short_fl = filter(cat_dt,  ID == m$ID[i])
  int_flow = filter(cat_dt,  ID == m$toID[i])
  new = st_union(short_fl$geom, int_flow$geom)
  cat_dt[int_flow$ind, geom := st_as_sf(new)]

  short_fl = filter(fl_dt,  ID == m$ID[i])
  int_flow = filter(fl_dt,  ID == m$toID[i])
  new = st_line_merge(st_union(short_fl$geom, int_flow$geom))
  fl_dt[int_flow$ind, geom := st_as_sf(new)]
  fl_dt[int_flow$ind, member_COMID := paste0(short_fl$member_COMID, int_flow$member_COMID)]
}

new_cat = filter(st_as_sf(cat_dt), !ID %in% m$ID) %>%
  mutate(areasqkm = as.numeric(st_area(.)/1e6)) %>%
  select(-ind)

new_fl =  filter(st_as_sf(fl_dt), !ID %in% m$ID) %>%
  flowpaths_to_linestrings() %>%
  mutate(lengthkm = as.numeric(st_length(.)/1e3)) %>%
  select(-ind, -areasqkm) %>%
  left_join(st_drop_geometry(new_cat))

#################################################
#################################################
validate_network(new_fl, new_cat)
write_sf(new_fl, "workflow/cache/test.gpkg", "fl")
write_sf(new_cat, "workflow/cache/test.gpkg", "cat")
# merge 1st orders that flow into other first orders
first_orders = function(fl){
  fl %>%
    filter(areasqkm < min_area, stream_order == 1) %>%
    intersection_mapping(fl) %>%
    filter(n == 1, IDso == 1, intSO == 1, new_area < max_area) %>%
    filter(new_area < max_area) %>%
    group_by(intID) %>%
    slice_min(new_area) %>%
    ungroup() %>%
    filter(!ID %in% intID)
}

obj = list(m1 = first_orders(new_fl),
          cat_dt = data.table(mutate(new_cat, ind = 1:n())),
          fl_dt  = data.table(mutate(new_fl,  ind = 1:n())))

len = 1000000000000
len_last = -1

while(len != len_last){
  len_last = len
  obj = dissolve_network(obj$m1, obj$cat_dt, obj$fl_dt, min_area, max_area, first_orders)
  dim(obj$cat_dt)
  dim(obj$fl_dt)
  len = nrow(obj$m1)
  print(len)
}

out_cat =  st_as_sf(obj$cat_dt) %>%
  mutate(areasqkm = as.numeric(st_area(.)/1e6)) %>%
  select(-ind)

out_fl = st_as_sf(obj$fl_dt) %>%
  select(-ind, -areasqkm)  %>%
  left_join(st_drop_geometry(out_cat))

#################################################
write
#################################################
validate_network(out_fl, out_cat)

t = filter(tmp, intID == 5382)

filter(out_cat, ID %in% t$ID) %>%
  mapview() +
  filter(out_fl, ID %in% t$ID) +
  mapview(filter(out_fl, ID  == 5382), color = "red") +
  mapview(filter(out_cat, ID == 5382), color = "red")
# m1 = out_fl %>%
#   filter(areasqkm < min_area, stream_order == 1) %>%
#   intersection_mapping(out_fl) %>%
#   filter(new_area < max_area) %>%
#   group_by(intID) %>%
#   slice_min(new_area) %>%
#   ungroup() %>%
#   filter(!ID %in% intID)

other_orders = function(fl){
  tmp = fl %>%
    filter(areasqkm < min_area, stream_order == 1) %>%
    intersection_mapping(fl) %>%
    filter(new_area < max_area) %>%
    group_by(intID) %>%
    arrange(-intSO, new_area) %>%
    slice(1) %>%
    ungroup() %>%
    filter(!ID %in% intID)
}

obj = list(m1     = other_orders(out_fl),
           cat_dt = data.table(mutate(out_cat, ind = 1:n())),
           fl_dt  = data.table(mutate(out_fl,  ind = 1:n()))
)

len = 1000000000000
len_last = -1

while(len != len_last){
  len_last = len
  obj = dissolve_network(m1 = obj$m1,
                         cat_dt = obj$cat_dt,
                         fl_dt = obj$fl_dt,
                         min_area, max_area, filter_fun = other_orders)
  len = nrow(obj$m1)
  print(dim(obj$cat_dt))
  print(dim(obj$fl_dt))
  print(len)
}

filter(obj$m1, ID == 5385)
filter(obj$m1, intID == 5385)

out_cat2 =  st_as_sf(obj$cat_dt) %>%
  mutate(areasqkm = as.numeric(st_area(.)/1e6)) %>%
  select(-ind)

out_fl2 = st_as_sf(obj$fl_dt) %>%
  select(-ind, -areasqkm)  %>%
  left_join(st_drop_geometry(out_cat2))

filter(out_cat2, ID == 5382) %>%
  mapview()

#################################################
#################################################
make_plot(out_fl, out_cat, "First Order Dissolve HF") %>%
  ggsave(filename = "workflow/cache/img/04-ngen-firstorder-dissolve.png",
         units = "in", height = 4, width = 8)

message("Dropped: ", nrow(cat) - nrow(new_cat_2), " features")

if(all(validate_network(out_fl, out_cat))){
  unlink(out_path)
  write_sf(out_cat, out_path, "catchments")
  write_sf(out_fl, out_path, "flowpaths")
} else {
  message("Error in: ", basename(rstudioapi::getSourceEditorContext()$path))
}
#################################################
#################################################
