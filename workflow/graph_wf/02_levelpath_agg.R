source('workflow/utils.R')

#############################################################################
path         <- "workflow/graph_wf/cache/ngen_01a-0.gpkg"
out_path     <- "workflow/graph_wf/cache/ngen_lp.gpkg"
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

# Preprocessing .... ------------------------------------------------------
cat = read_sf(path, "catchment") %>% 
  mutate(areasqkm = as.numeric(set_units(st_area(.), "km2")))

fl  = read_sf(path, "flowline") %>% 
  left_join(select(st_drop_geometry(cat), areasqkm, ID)) %>% 
  filter(levelpath == 2302807)

filter(fl, member_comid == "7731675")
table(st_geometry_type(fl))

singles = fl %>%
  group_by(levelpath) %>%
  mutate(n = n()) %>%
  filter(n == 1)  %>%
  ungroup() %>% 
  select(ID, member_comid)

table(st_geometry_type(singles))

multis = fl %>%
  filter(!ID %in% singles$ID) %>% 
  st_drop_geometry()  %>%
  group_by(levelpath) %>%
  mutate(n = n())     %>%
  filter(n > 1)       %>%
  arrange(-hydroseq)  %>%
  mutate(ind = cs_group(.data$areasqkm, ideal_size)) %>%
  ungroup() 

filter(multis, member_comid == "7731675")$ind

mapping = multis %>%
  group_by(levelpath, ind) %>%
  summarize(member_comid = paste(member_comid, collapse = ","), n = n()) %>%
  ungroup() 

tab1 =  multis %>%
  select(ID, levelpath, ind) %>% 
  left_join(select(fl, ID), by = "ID") %>%
  group_by(levelpath, ind) %>% 
  mutate(n1 = n(), ID = ID[1]) %>% 
  ungroup()

no_need_to_merge = filter(tab1, n1 == 1) %>% 
  bind_rows(singles) %>% 
  st_as_sf()

to_merge = st_as_sf(filter(tab1, n1 > 1)) %>% 
  setDT()

mapview(filter(st_as_sf(to_merge), ind == 2)) + fl

fl_xx = to_merge[, .(geom   = st_line_merge(st_union(geom))), by = list(levelpath, ind)]


table(st_geometry_type(st_as_sf(fl_xx)))

out_fl = st_as_sf(fl_xx) %>% 
  filter(st_geometry_type(.) == "MULTILINESTRING")
  
mapview(out_fl[4,])
  
  
  inner_join(select(mapping, ID, member_comid), by = "ID") %>% 
  flowpaths_to_linestrings()

table(st_geometry_type(out_fl))


%>% 
  bind_rows(no_need_to_merge) %>% 
  mutate(lengthkm = as.numeric(set_units(st_length(.), "km"))) %>% 
  flowpaths_to_linestrings()
  
######################

single_cat   = filter(cat, ID %in% singles$ID) %>% 
  select(ID)

no_merge_cat = inner_join(cat, filter(mutate(tab1, geom = NULL), n1 == 1), by = "ID") %>% 
  select(levelpath, ind, geom)

merge_cat    = setDT(inner_join(cat, filter(mutate(tab1, geom = NULL), n1 > 1), by = "ID"))


cat_xx = merge_cat[, .(geom = st_union(geom)), by = list(levelpath, ind)]

merged = do.call(rbind, list(st_as_sf(cat_xx), no_merge_cat)) %>% 
  left_join(mapping, by = c("levelpath", 'ind')) %>% 
  dplyr::select(ID) 

out_cat = bind_rows(merged, st_as_sf(single_cat)) %>% 
  mutate(areasqkm = as.numeric(set_units(st_area(.), "km2"))) 

dim(out_cat)
dim(out_fl)
#--------------------------------------------------
#############################################################################

message("Dropped: ", nrow(fl) - nrow(out_fl), " features")

xx = prep_ngen(fl = out_fl, cat = out_cat, ID_col = "ID")

unlink(out_path)
write_sf(xx$nex, out_path,  "nex")
write_sf(xx$fl, out_path, "flowline")
write_sf(xx$cat, out_path, "catchment")

# make_plot(out_fl, out_cat, "LevelPath Dissolve HF") %>%
#   ggsave(
#     filename = "workflow/nhd_workflows/cache/img/02-ngen-levelpath-dissolve.png",
#     units = "in",
#     height = 4,
#     width =  8
#   )

