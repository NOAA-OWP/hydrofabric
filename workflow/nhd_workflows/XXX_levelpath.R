source('workflow/utils.R')

#############################################################################
path         <- "workflow/nhd_workflows/cache/ngen_01a-0.gpkg"
out_path     <- "workflow/nhd_workflows/cache/ngen_01a-1.gpkg"
ideal_size   <-  10   # in km2
#############################################################################

# TO DISSOLVE LEVELPATHS WE NEED:
#   areasqkm
#   lengthkm
#   levelpath
#   hydroseq
# 
# RECORD:
#   ngen_id
#   areasqkm
#   lengthkm
#   member_comid



test = 2239178
# Preprocessing .... ------------------------------------------------------
cat = read_sf(path, "catchment")

fl  = read_sf(path, "flowline") %>%
  mutate(member_COMID = comid, areasqkm = NULL ) %>% 
  right_join(st_drop_geometry(cat), by = "comid") %>%
  mutate(areasqkm = ifelse(is.na(areasqkm), 0, areasqkm)) %>% 
  st_cast("LINESTRING")

# fl_net = sfnetworks::as_sfnetwork(select(fl,comid)) %>% 
#   sfnetworks::activate(edges) %>% 
#   filter(comid == 721640)
# 
# new_id_mapping = data.frame(ngen = 1:nrow(fl),
#                             comid = fl$comid,
#                             to_comid = fl$to_comid,
#                             from_comid = fl$from_comid) %>% 
#   mutate(to_comid = match(to_comid, comid),
#          from_comid = match(from_comid, comid)) %>% 
#   mutate(to_comid = ifelse(is.na(to_comid), -1, to_comid),
#          from_comid = ifelse(is.na(to_comid), -1, from_comid))
#   
#                             to_ngen = fl$to_comid,
#                             from_ngen = fl$from_comid,)


  

singles = fl %>%
  group_by(levelpath) %>%
  mutate(n = n()) %>%
  filter(n == 1)  %>%
  mutate(member_COMID = as.character(comid)) %>% 
  ungroup() %>% 
  dplyr::select(comid, member_COMID, to_comid, from_comid, order, hydroseq)


multis = fl %>%
  st_drop_geometry()  %>%
  group_by(levelpath) %>%
  mutate(n = n())     %>%
  filter(n > 1)       %>%
  arrange(-hydroseq)  %>%
  mutate(ind = cs_group(.data$areasqkm, ideal_size)) %>%
  arrange(hydroseq)   %>%
  ungroup() 

mapping = multis %>%
  group_by(levelpath, ind) %>%
  arrange(hydroseq) %>%
  summarize(
    member_COMID = paste(comid, collapse = ","),
    comid        = first(comid),
    hydroseq     = first(hydroseq),
    to_comid     = first(to_comid),
    from_comid   = last(from_comid),
    order = max(order)
  ) %>%
  ungroup()

tab1 =  multis %>%
  select(comid, levelpath, ind) %>% 
  left_join(select(fl, comid), by = "comid") %>%
  group_by(levelpath, ind) %>% 
  mutate(n1 = n()) %>% 
  ungroup()

to_merge = filter(tab1, n1 > 1) %>% 
  setDT()

no_need_to_merge = filter(tab1, n1 == 1) %>% 
  select(levelpath, ind, geom) %>% 
  setDT()

fl_xx = to_merge[, .(geom   = st_union(geom)), by = list(levelpath, ind)]

out_fl = data.table::rbindlist(list(fl_xx, no_need_to_merge)) %>% 
  st_as_sf() %>%
  left_join(mapping, by = c("levelpath", 'ind')) %>%
  bind_rows(singles) %>% 
  mutate(lengthkm = as.numeric(set_units(st_length(.), "km")))

######################

single_cat   = filter(cat, comid %in% singles[['comid']]) %>% 
  setDT()
no_merge_cat = filter(cat, comid %in% filter(tab1, n1 == 1)[['comid']]) %>% 
  setDT()
merge_cat    = inner_join(cat, mutate(filter(tab1, n1 > 1), geom = NULL), by = "comid") %>% 
  setDT()

xx = merge_cat[, .(geom = st_union(geom)), by = list(levelpath, ind)]

merged = st_as_sf(xx) %>%
  left_join(dplyr::select(mapping, comid, ind, levelpath), by = c("levelpath", 'ind')) %>% 
  dplyr::select(comid) %>% 
  mutate(areasqkm = as.numeric(set_units(st_area(.), "km2"))) %>% 
  dplyr::select(comid,areasqkm,geom) %>% 
  st_as_sf() %>% 
  st_collection_extract("POLYGON") %>% 
  setDT()


out_cat = data.table::rbindlist(list(merged, no_merge_cat, single_cat)) %>% 
  st_as_sf() %>%
  mutate(areasqkm = as.numeric(set_units(st_area(.), "km2")))

#--------------------------------------------------
#############################################################################

make_plot(out_fl, out_cat, "LevelPath Dissolve HF") %>%
  ggsave(
    filename = "workflow/nhd_workflows/cache/img/02-ngen-levelpath-dissolve.png",
    units = "in",
    height = 4,
    width =  8
  )

message("Dropped: ", nrow(fl) - nrow(out_fl), " features")

unlink(out_path)
write_sf(out_cat, out_path, "catchment")
write_sf(out_fl, out_path, "flowline")

