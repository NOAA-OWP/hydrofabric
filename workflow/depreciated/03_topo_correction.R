source('workflow/utils.R')

#############################################################################
path       <- "workflow/cache/ngen_01a-4.gpkg"
out_path   <- "workflow/cache/ngen_01a-5.gpkg"
#############################################################################
cat = read_sf(path, "catchments")
fl  = read_sf(path, "flowpaths")

nodes <- fl %>%
  st_coordinates() %>%
  as_tibble() %>%
  rename(edgeID = L2) %>%
  group_by(edgeID) %>%
  slice(c(1, n())) %>%
  ungroup() %>%
  mutate(start_end = rep(c('start', 'end'), times = n()/2)) %>%
  mutate(xy = paste(.$X, .$Y)) %>%
  group_by(xy)

nodes$nodeID =  group_indices(nodes)

nodes = nodes %>%
    ungroup() %>%
    select(-xy)%>%
    distinct(nodeID, .keep_all = TRUE) %>%
    select(-c(edgeID, start_end)) %>%
    st_as_sf(coords = c('X', 'Y'), crs = st_crs(fl)) %>%
    mutate(n = lengths(st_intersects(., fl)))

nn = filter(nodes, n > 1)
write_sf(nn, "workflow/cache/nodes2.gpkg")

#############################################################################
message("Adjusting to ", nrow(nn), " nodes...")

ll = st_intersects(nn, fl)

mapping  <- data.frame(nodeID  = rep(nn$nodeID, sapply(ll, length)),
                       intID   = fl$ID[unlist(ll)])  %>%
  filter(nodeID != intID)

cat2 = mutate(cat, ind = as.integer(1:n())) %>%
  data.table()

for (x in 1:nrow(nn)) {
  tmp_cat = filter(cat2, ID %in% filter(mapping, nodeID %in% nn$nodeID[x])$intID) %>%
    st_as_sf()
  cat2[tmp_cat$ind, geom := st_geometry(st_snap(st_geometry(tmp_cat),
                                             st_geometry(nn[x,]), tolerance = 60))]

}

out =  st_as_sf(cat2)

# THIS HAS TO BE DONE THIS WAY. st_make_valid on 'out' fails.
inv = which(!sf::st_is_valid(out))
invalid = out[inv,] %>% st_make_valid()
out$geom[inv] = invalid$geom

write_sf(out, "workflow/cache/cats2.gpkg")

out = out %>%
  mutate(areasqkm = as.numeric(st_area(.)/1e6)) %>%
  select(-ind) %>%
  st_cast("POLYGON")

unlink(out_path)
write_sf(out, out_path, "catchments")
write_sf(fl, out_path, "flowpaths")

