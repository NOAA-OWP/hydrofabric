gpkg = "data/rf01.gpkg"

cat = read_sf(gpkg, layer = "catch")
fl  = read_sf(gpkg, layer = "flowpath")

# issues: there are silver of the pigtails reamining! can be ideitifed by
#          duplicate IDs, need to clean up and fill holes...

explored_cat = rmapshaper::ms_explode(cat) %>%
  mutate(area = as.numeric(st_area(.)))

frags = explored_cat %>%
  filter(area %in% c(450, 900)) %>% # These come from the 30m DEM so are either full or half cells
  rmapshaper::ms_dissolve() %>%
  rmapshaper::ms_explode()

non_frags = explored_cat %>%
  filter(!area %in% c(450, 900)) %>%
  select(ID)


ints = st_intersection(frags, select(cat, ID)) %>%
  st_collection_extract("POLYGON") %>%
  mutate(area = as.numeric(st_area(.))) %>%
  st_drop_geometry() %>%
  group_by(rmapshaperid, ID) %>%
  summarise(area = sum(area)) %>%
  ungroup() %>%
  group_by(rmapshaperid) %>%
  slice_max(area, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(frags) %>%
  select(-area, -rmapshaperid)

tmp_cat = bind_rows(ints, non_frags) %>%
  st_as_sf() %>%
  group_by(ID) %>%
  summarise()

# Interior Ring Removal and Fill ------------------------------------------
# example = 15192 (outside), 15193 (interior)

inters = st_intersects(tmp_cat)

interior_rings2 = data.frame(ID = tmp_cat$ID,
                             n  = lengths(inters)) %>%
  filter(n == 2)

int_list = inters[which(tmp_cat$ID %in% interior_rings2$ID)]
lengths(int)

all = tmp_cat$ID[unlist(int_list)]

exteriors = all[!all %in% interior_rings2$ID]

outs = list()

for(i in 1:nrow(interior_rings2)){
  outs[[i]] = tmp_cat[int_list[[i]],] %>%
    mutate(area = st_area(.), id = 1) %>%
    arrange(-area) %>%
    group_by(id) %>%
    summarize(ID = ID[1])
}

new_exteriors = bind_rows(outs) %>%
  select(-id) %>%
  group_by(ID) %>%
  summarise()

new_cat = tmp_cat %>%
  filter(!ID %in% all) %>%
  bind_rows(new_exteriors)

write_sf(new_cat, "data/hopefull_2.gpkg")




