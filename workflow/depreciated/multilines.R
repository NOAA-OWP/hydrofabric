path       <- "workflow/cache/ngen_01a-5.gpkg"
out_path   <- "workflow/cache/ngen_01a-6.gpkg"
## ISSUES 11065, disconected from network
#############################################################################
cat = read_sf(path, "catchments")
fl  = read_sf(path, "flowpaths")

fl2 = flowpaths_to_linestrings(fl)

mapping  <- data.frame(flID    = rep(fl$ID, sapply(ll, length)),
                       intID   = fl$ID[unlist(ll)],
                       so      = fl$stream_order[unlist(ll)]) %>%
  filter(flID != intID)

#########################

ii = which(st_geometry_type(fl2) == "MULTILINESTRING")
ml_fl = fl[ii,]
unlink("workflow/cache/fls.gpkg")
write_sf(ml_fl, "workflow/cache/fls.gpkg", "raw")

for(i in 1:nrow(ml_fl)){

  test = ml_fl[i,]
  ints = filter(mapping, flID == test$ID)
  poss = filter(fl, ID %in% ints$intID)

  exp = ms_explode(test) %>%
    mutate(l = st_length(.)) %>%
    select(l, ID) %>%
    mutate(core =  lengths(st_intersects(., st_buffer(poss,1))))

  edges = node_filter(fls = exp)

  #mapview(edges) + test

  new = st_union(edges$geom)

# mapview(new) + test

if(st_geometry_type(new) == "MULTILINESTRING"){
  new = st_line_merge(new)
}

  ml_fl$geom[i] = new

}

iii = which(st_geometry_type(ml_fl) == "MULTILINESTRING")
dangles = ml_fl[iii,]



# if(st_geometry_type(new) == "MULTILINESTRING"){
#
#   d = st_geometry(st_difference(exp, new))
#
#   n = d[lengths(st_touches(d, ms_explode(new))) == 2, ]
#   new2 = st_union(new, n)
#
#   fls = st_as_sf(ms_explode(new2) ) %>%
#     filter(!duplicated(.))
#
#   nodes <- fls %>%
#     st_coordinates() %>%
#     as_tibble() %>%
#     rename(edgeID = L1) %>%
#     group_by(edgeID) %>%
#     slice(c(1, n())) %>%
#     ungroup() %>%
#     mutate(start_end = rep(c('start', 'end'), times = n()/2)) %>%
#     mutate(xy = paste(.$X, .$Y)) %>%
#     group_by(xy)
#
#   nodes$nodeID =  group_indices(nodes)
#
#   tmp = fls %>%
#     mutate(from =  filter(nodes, start_end == 'start')$nodeID,
#            to   =  filter(nodes, start_end == 'end')$nodeID,
#            l = st_length(.))  %>%
#     group_by(to) %>%
#     mutate(nTO = n()) %>%
#     ungroup() %>%
#     group_by(from) %>%
#     mutate(nFROM = n()) %>%
#     ungroup()
#
#   new = filter(tmp, !(nTO == 2 & nFROM == 2)) %>%
#     st_union()
#
#
#   if(st_geometry_type(new) == "MULTILINESTRING"){
#     new = st_line_merge(new)
#   }
#
#   if(st_geometry_type(new) == "MULTILINESTRING"){
#     solo     = as.numeric(which(table(c(tmp$from, tmp$to)) == 1))
#
#     if(length(solo) == max(tmp$from, tmp$to)){
#       pt = filter(nodes, nodeID == solo[2]) %>% st_as_sf(coords = c("X", "Y"), crs = 5070)
#       new = st_snap(new, st_geometry(pt),  tolerance = 60)
#     } else {
#       solo = solo[!solo %in% c(1, max(tmp$from, tmp$to))]
#       new = filter(tmp, !to %in% solo, !from %in% solo ) %>%
#         st_union()
#     }
#   }
#
#   if(st_geometry_type(new) == "MULTILINESTRING"){
#     new = st_line_merge(new)
#   }
#
# }

}

write_sf(ml_fl, "workflow/cache/fls.gpkg", "fixed", append = T)

iii = which(st_geometry_type(ml_fl) == "MULTILINESTRING")
ml_fl2 = fl[ii,]
dangles = ml_fl2[iii,]



fls = ms_explode(dangles[1,])

mapview(ml_fl[339,])
ml_fl$ID[iii]

