library(dplyr)
library(sf)

#############################################################################
path = "/Users/mjohnson/github/gfv2.0/workspace/cache/01a.gpkg"
out_path = "data/ngen_01a-0.gpkg"
#############################################################################

# Notes: odd projections....
#        some duplicated features
#        overly complex spatial topologies
#        more flowpaths then catchments

cat = read_sf(gpkg, layer = "divides") %>%
  st_set_crs(5070) %>%
  filter(!duplicated(.)) %>%
  rmapshaper::ms_simplify(.9) %>%
  mutate(area = as.numeric(st_area(.)/1e6))

which(duplicated(cat$ID))

fl  = read_sf(gpkg, layer = "reconciled") %>%
  st_set_crs(4326) %>%
  st_transform(5070) %>%
  filter(!duplicated(.))

fl  = filter(fl, ID %in% cat$ID)
cat = filter(cat, ID %in% fl$ID)

#--------------------------------------------------
unlink(out_path)
write_sf(cat, out_path, "catchments")
write_sf(fl, out_path, "flowpaths")
