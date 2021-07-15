source('workflow/utils.R')
#############################################################################
path     <- "workflow/nhd_workflows/cache/ngen_01a-0.gpkg"
out_path <- "workflow/nhd_workflows/cache/ngen_01a-1.gpkg"
#############################################################################

# Catchments only carry around geometry and ID
in_cat <- read_sf(path, layer = "catchment") 
in_fl = read_sf(path, layer = "flowline") 

#############################################################################

cat <- in_cat %>%
  ms_explode() %>%
  filter(!duplicated(.)) %>%
  mutate(area = as.numeric(st_area(.)))

#write_sf(cat, "workflow/cache/all-cat.gpkg")

ids <- filter(cat, duplicated(cat$comid))

cat_no_problem <- filter(cat, !comid %in% ids$comid)

#write_sf(cat_no_problem, "workflow/cache/no-prob.gpkg")

challenges = filter(cat, comid %in% ids$comid) %>%
  mutate(tmpID = 1:n())

#write_sf(challenges, "workflow/cache/chal.gpkg")

base_cats = challenges %>%
  group_by(comid) %>%
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

message(nrow(frags), " Consoldated fragments...")

ints = st_intersection(frags, st_make_valid(base_cats)) %>%
  st_collection_extract("LINESTRING") %>%
  mutate(l = st_length(.)) %>%
  group_by(rmapshaperid) %>%
  slice_max(l, with_ties = FALSE)

#write_sf(t, "workflow/cache/ints.gpkg")

tj = right_join(frags,
                dplyr::select(st_drop_geometry(ints), comid, rmapshaperid),
                by = "rmapshaperid") %>%
  bind_rows(base_cats) %>%
  group_by(comid) %>%
  mutate(n = n())

tt = filter(tj, n > 1) %>%
  summarise() %>%
  bind_rows(filter(tj, n == 1))

cat = tt %>%
  dplyr::select(comid) %>%
  rename(geometry = geom) %>%
  #ms_simplify(.9) %>%
  mutate(areasqkm = as.numeric(set_units(st_area(.), "km2")))

#--------------------------------------------------
#############################################################################

make_plot(in_fl, cat, "Fragment Dissolve HF") %>%
  ggsave(filename = "workflow/nhd_workflows/cache/img/01-ngen-frag-dissolve.png",
         units = "in", height = 4, width = 8)

unlink(out_path)
write_sf(cat, out_path, "catchments")
write_sf(in_fl, out_path, "flowpaths")





