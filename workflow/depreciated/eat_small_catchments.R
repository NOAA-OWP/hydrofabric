source('workflow/utils.R')

# Rules:
  # 1. Merge small into the downstream catchment/flowpath

#############################################################################
path       <- "workflow/cache/ngen_01a-2.gpkg"
out_path   <- "workflow/cache/ngen_01a-3.gpkg"
min_size   <- 3 #km2
#############################################################################

# Preprocessing .... ------------------------------------------------------

cat = read_sf(path, layer = "catchments") %>%
  mutate(area = as.numeric(st_area(.) / 1e6))

fl  = read_sf(path, layer = "flowpaths") %>%
  mutate(length = as.numeric(st_length(.) / 1e3))

to_small  = filter(cat, area < min_size)
message("Need to dissolve ", nrow(to_small), " catchments")

# Will always dissolve into the downstream "toID"

snodes = nhdplusTools::get_node(fl, position = "start") %>%
  mutate(fpID = fl$ID, type = "start")

enodes = nhdplusTools::get_node(fl, position = "end") %>%
  mutate(fpID = fl$ID, type = "end")

# t = bind_rows(snodes, enodes) %>%
#   mutate(X = st_coordinates(.)[,1], Y = st_coordinates(.)[,2]) %>%
#   st_drop_geometry() %>%
#   arrange(fpID, type) %>%
#   mutate(g = group_indices(.,X,Y))

mapview(filter(fl, ID == 674)) +
  t[1:5,]


nexi = bind_rows(snodes, enodes)

lwgeom::st_split(fl, nexi)


t = bind_rows(snodes, enodes) %>%
  select(fpID, type) %>%
  st_intersection(select(fl, ID)) %>%
  filter(fpID != ID)


tmids = c(674,669,675,71)#filter(nexi, g == 43)$fpID

mapview(filter(cat, ID %in% tmids)) +
  mapview(filter(fl, ID %in% tmids))



nexi = bind_rows(snodes, enodes) %>%
  mutate(X = st_coordinates(.)[,1], Y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  arrange(fpID) %>%
  mutate(g = group_indices(.,X,Y)) %>%
  group_by(g) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  #filter(count != 1)  %>%
  arrange(g)

st_as_sf(nexi, coords = c("X", "Y"), crs =  5070) %>%
  write_sf(enodes,"workflow/cache/nodes.gpkg")

filter()

t2

tmids = c(675, 674, 669, 678, 677, 473)#filter(nexi, g == 43)$fpID

mapview(filter(cat, ID %in% tmids)) +
  mapview(filter(fl, ID %in% tmids))



filter(nexi, g == 25477)


spatial_nodes = nexi %>%
  tidyr::pivot_wider(id_cols = c(X, Y), names_from = type, values_from = fpID) %>%
  mutate(nex_id == 1:n())



# Merging... --------------------------------------------------------------

for(i in 1:nrow(to_small)){
  # Find all catchments that get flow from the small offender
  t1 = spatial_nodes[which(sapply(spatial_nodes$fromFP, function(y) to_small$ID[i] %in% y)),]
  id = filter(cat, ID %in% unlist(t1$toFP)) %>%
      slice_max(area) %>%
      pull(ID)

  if(length(id) == 0){
    id = filter(cat, ID %in% unlist(t1$fromFP)) %>%
      slice_max(area) %>%
      pull(ID)
  }

  cat$ID[which(cat$ID == to_small$ID[i])] = id
  fl$ID[which(fl$ID == to_small$ID[i])] = id

  message(i)
}

dup_cat = filter(cat, )

list1 = spatial_nodes$toFP
t = data.frame(toFP = unlist(list1), index = rep(seq(length(list1)), lengths(list1)))
list2 = spatial_nodes$fromFP
t2 = data.frame(fromFP = unlist(list2), index = rep(seq(length(list2)), lengths(list2)))


indexes = list()
for(i in 1:nrow(to_small)){

  tmp = t2$index[which(t2$fromFP == to_small$ID[i])]
  inds = t$toFP[which(t$index == tmp)]

if(length(inds) == 0){
  tmp = t$index[which(t$toFP == to_small$ID[i])]
  inds = t2$fromFP[which(t2$index == tmp)]
}

indexes[[i]] = inds
}

tm = unlist(indexes)
length(tm)


tt = rep(1:nrow(spatial_nodes), each = lengths(spatial_nodes$toFP))
unlist(spatial_nodes$toFP)


unlist(t1$toFP)


t1 = spatial_nodes[which(sapply(spatial_nodes$fromFP, function(y) to_small$ID %in% y)),]


tt = which(sapply(to_small$ID, function(y)  spatial_nodes$fromFP %in% y))
tt[1:5]


oo = spatial_nodes$toFP[tt]
oo[1]

