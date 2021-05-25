source('workflow/utils.R')

#############################################################################
path       <- "workflow/cache/ngen_01a-1.gpkg"
out_path   <- "workflow/cache/ngen_01a-2.gpkg"
frag_def   <-  1000 # m2
#############################################################################

# Preprocessing .... ------------------------------------------------------

cat = read_sf(path, layer = "catchments")
fl  = read_sf(path, layer = "flowpaths")

# issues: there are silver of the pigtails remaining! can be ideitifed by
#          duplicate IDs, need to clean up and fill holes...
exploded_cat = rmapshaper::ms_explode(cat) %>%
  mutate(area = as.numeric(st_area(.)))

frags = exploded_cat %>%
  # These come from the 30m DEM so are either full or half cells
  # The bulk then are either 450, or 900. But there are a few other
  # fragments that we can clean with this step
  filter(area  < frag_def) %>%
  # Dissolve internal edges
  ms_dissolve() %>%
  # explode external edges
  ms_explode()

non_frags = exploded_cat %>%
  # non-frags are "ideal" catchments
  filter(area >= frag_def) %>%
  select(ID)

cat_valid  = st_make_valid(cat)

ints = st_intersection(frags, select(cat_valid, ID)) %>%
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
  select(-area, -rmapshaperid) %>%
  st_as_sf() %>%
  group_by(ID) %>%
  summarise()

# flowlines  --------------------------------------------------------------

int_lines = filter(fl, ID %in% ints$ID) %>%
  group_by(ID) %>%
  summarize(outlet = which.min(Hydroseq),
            LevelPathID = LevelPathID[1],
            ID = min(ID),
            toID = toID[outlet],
            Hydroseq = Hydroseq[outlet],
            member_COMID = paste(member_COMID, collapse = ",")) %>%
  ungroup() %>%
  select(-outlet) %>%
  st_line_merge()

tmp_lines = bind_rows(int_lines, filter(fl, ID %in% non_frags$ID) ) %>%
  group_by(ID) %>%
  mutate(n = n())

dup_lines = filter(tmp_lines, n > 1) %>%
  st_as_sf() %>%
  summarize(outlet = which.min(Hydroseq),
            LevelPathID = LevelPathID[1],
            ID = min(ID),
            toID = toID[outlet],
            Hydroseq = Hydroseq[outlet],
            member_COMID = paste(member_COMID, collapse = ",")) %>%
  ungroup()

all_lines = filter(tmp_lines, n == 1) %>%
  st_as_sf() %>%
  bind_rows(dup_lines) %>%
  select(-n)

# cats  --------------------------------------------------------------

tmp_cat = bind_rows(ints, non_frags) %>%
  group_by(ID) %>%
  mutate(n = n())

dups = filter(tmp_cat, n > 1) %>%
  st_as_sf() %>%
  summarise()

all_cats = filter(tmp_cat, n == 1) %>%
  st_as_sf() %>%
  bind_rows(dups) %>%
  select(-n)

all(validate_network(all_lines, all_cats))

# Interior Ring Removal and Fill ------------------------------------------
# example = 15192 (outside), 15193 (interior)

inters = st_intersects(all_cats)

interior_rings = data.frame(ID = all_cats$ID,
                            n = lengths(inters)) %>%
  filter(n == 2)

int_list = inters[which(all_cats$ID %in% interior_rings$ID)]

all = all_cats$ID[unlist(int_list)]

exteriors = all[!all %in% interior_rings$ID]

outs = list()

for(i in 1:nrow(interior_rings)){
  outs[[i]] = all_cats[int_list[[i]],] %>%
    mutate(area = st_area(.), id = 1) %>%
    arrange(-area) %>%
    group_by(id) %>%
    summarize(ID = ID[1])
}

new_exteriors = bind_rows(outs) %>%
  select(-id) %>%
  group_by(ID) %>%
  summarise()

new_cat = all_cats %>%
  filter(!ID %in% all) %>%
  bind_rows(new_exteriors)

new_fl = all_lines %>%
  filter(ID %in% new_cat$ID)

#--------------------------------------------------
message("Dropped: ", nrow(cat) - nrow(new_cat), " features")

if(all(validate_network(new_fl, new_cat))){
  unlink(out_path)
  write_sf(new_cat, out_path, "catchments")
  write_sf(new_fl, out_path, "flowpaths")
} else {
  message("Error in: ", basename(rstudioapi::getSourceEditorContext()$path))
}




