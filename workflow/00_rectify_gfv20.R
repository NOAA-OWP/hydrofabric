source('workflow/utils.R')

#############################################################################
path     <- "workflow/cache/gfv20/01a.gpkg"
out_path <- "workflow/cache/ngen_01a-0.gpkg"
#############################################################################

# Catchments only carry around geometry and ID
in_cat <- read_sf(path, layer = "divides") %>%
  st_set_crs(5070) %>%
  select(ID) %>%
  mutate(areasqkm = as.numeric(st_area(.))/1e6)

in_fl = read_sf(path, layer = "reconciled") %>%
  mutate(lengthkm = as.numeric(st_length(.))/1e3, LENGTHKM = NULL) %>%
  st_transform(5070)

make_plot(in_fl, in_cat, "Raw GFv2.0 HF") %>%
  ggsave(filename = "workflow/cache/img/00-raw-gf-v20.png",
         units = "in", height = 4, width = 8)

#############################################################################

cat <- in_cat %>%
  ms_explode() %>%
  filter(!duplicated(.)) %>%
  mutate(area = as.numeric(st_area(.)))

#write_sf(cat, "workflow/cache/all-cat.gpkg")

ids <- filter(cat, duplicated(cat$ID))

cat_no_problem <- filter(cat, !ID %in% ids$ID)

#write_sf(cat_no_problem, "workflow/cache/no-prob.gpkg")

challenges = filter(cat, ID %in% ids$ID) %>%
  mutate(tmpID = 1:n())

#write_sf(challenges, "workflow/cache/chal.gpkg")

base_cats = challenges %>%
  group_by(ID) %>%
  slice_max(area) %>%
  bind_rows(cat_no_problem)

#write_sf(base_cats, "workflow/cache/base.gpkg")

fragments = filter(challenges, !tmpID %in% base_cats$tmpID)

message(nrow(fragments), " fragments to clean...")

frags = fragments %>%
  ms_dissolve() %>%
  ms_explode() %>%
  mutate(area = as.numeric(st_area(.))) %>%
  st_make_valid()

#write_sf(frags, "workflow/cache/frags.gpkg")

message(nrow(frags), " Consoldated fragmnets...")

ints = st_intersection(frags, st_make_valid(base_cats)) %>%
  st_collection_extract("LINESTRING") %>%
  mutate(l = st_length(.)) %>%
  group_by(rmapshaperid) %>%
  slice_max(l, with_ties = FALSE)

#write_sf(t, "workflow/cache/ints.gpkg")

tj = right_join(frags, select(st_drop_geometry(ints), ID, rmapshaperid), by = "rmapshaperid") %>%
  bind_rows(base_cats) %>%
  group_by(ID) %>%
  mutate(n = n())

tt = filter(tj, n > 1) %>%
  summarise() %>%
  bind_rows(filter(tj, n == 1))

cat = tt %>%
  select(ID) %>%
  rename(geometry = geom) %>%
  #ms_simplify(.9) %>%
  mutate(areasqkm = as.numeric(st_area(.)/1e6))

in_fl$stream_order = nhdplusTools::get_streamorder(in_fl)
#
# validate_network(fl, cat)
#
# ### FIXING THE EXTRA CATCHMENT PROBLEM!!
#
# #Flowlines without a catchment...
# homeless =  fl %>% filter(!ID %in% cat$ID)
# flow_to_homeless = filter(fl, toID %in% homeless$ID)
# homeless_flows_to = filter(fl, ID %in% homeless$toID)
# keystone_homeless = filter(homeless, ID %in% flow_to_homeless$toID)
# inlet_homeless = filter(homeless, !ID %in% keystone_homeless$ID)


fl2  = inner_join(st_as_sf(in_fl), st_drop_geometry(cat), by = 'ID')


#--------------------------------------------------
#############################################################################

make_plot(fl2, cat, "Fragment Dissolve HF") %>%
  ggsave(filename = "workflow/cache/img/01-ngen-frag-dissolve.png",
         units = "in", height = 4, width = 8)


#if(all(validate_network(fl2, cat))){
unlink(out_path)
write_sf(cat, out_path, "catchments")
write_sf(fl2, out_path, "flowpaths")
#} else {
#  message("Error in: ", basename(rstudioapi::getSourceEditorContext()$path))
#}




