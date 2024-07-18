library(RNetCDF)
library(hydrofabric)

base    = "/Volumes/MyBook/conus-hydrofabric/v20.1"
gpkg    = glue("{base}/conus.gpkg")
network = glue("{base}/conus_net.parquet")
nc      = open.nc('/Users/mjohnson/Downloads/diffusive_natural_xs.nc')

vars = c("x", "y", "link", "feature_id", "n", "z",  "xid_d")

raw = lapply(
  1:length(vars),
  FUN = function(x) {
    var.get.nc(nc, vars[x], unpack = TRUE)
  }
) %>%
  bind_cols() %>%
  setNames(vars) %>%
  st_as_sf(coords = c("x", "y"), crs  = 6349)

transects = raw %>%
  group_by(link) %>%
  arrange(xid_d) %>%
  slice(c(1, n())) %>%
  summarize(link = link[1]) %>%
  st_cast("LINESTRING")

fps = read_sf(gpkg, "flowpaths")

net = read_parquet(network) %>%
  select(id, link = hf_id, vpu, hf_hydroseq, hydroseq, mainstem)

transects = left_join(transects, net, by = "link") %>%
  distinct() %>%
  st_transform(fps)

subset = st_filter(fps, transects)

map = st_join(transects, subset) %>%
  st_drop_geometry() %>%
  select(link, touches = id)

transects = left_join(transects, map, by = "link", relationship = "many-to-many") %>%
  filter(id == touches)  %>%
  rename(hy_id = id) %>%
  mutate(cs_source = "ras",
         cs_length = st_length(.)) %>%
  group_by(hy_id) %>%
  mutate(cs_id = rank(-rank(hf_hydroseq))) %>%
  ungroup() %>%
  select(hy_id, mainstem, cs_source, cs_id, link)

subset = filter(subset, id %in% transects$hy_id)

map2 = select(st_drop_geometry(transects), hy_id, cs_id, link)

cs_pts = left_join(raw, map2, by = "link") %>%
  mutate(
    pt_id = feature_id + 1,
    X = st_coordinates(.)[, 1],
    Y = st_coordinates(.)[, 2],
    Z_source  = "ras",
    roughness = n
  ) %>%
  group_by(hy_id, cs_id) %>%
  mutate(
    mind = min(xid_d),
    relative_dist = xid_d - mind,
    pt_measure = 100 * round(relative_dist / max(relative_dist), 4)
  ) %>%
  ungroup() %>%
  select(hy_id,
         cs_id,
         pt_id,
         X,
         Y,
         Z = z,
         Z_source,
         roughness,
         relative_dist,
         pt_measure)

div = filter(read_sf(gpkg, "divides"), id %in% subset$id)

unlink("data/troute_test.gpkg")
write_sf(transects, "data/troute_test.gpkg", 'transects')
write_sf(cs_pts, "data/troute_test.gpkg", 'cs_pts')
write_sf(select(subset, hy_id = id, divide_id),
         "data/troute_test.gpkg",
         'flowpaths')
write_sf(div, "data/troute_test.gpkg", 'divides')

write_parquet(st_drop_geometry(cs_pts))

write