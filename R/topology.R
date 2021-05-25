topology_doctor = function(graph){

  nexi <- graph$nodes
  fl   =  graph$edges
  cat = graph$cat
  n = filter(nexi, count > 1)


  for (x in 1:nrow(n)) {
    tmp_fl = filter(fl, to %in% n$nodeID[x] | from %in% n$nodeID[x])
    tmp_cat = filter(cat, ID %in% tmp_fl$catID )
    st_geometry(cat[cat$ID %in% tmp_cat$ID,]) = st_geometry(st_snap(tmp_cat, n[x,], tolerance = 60))
  }

  cat = st_difference(cat)

  o <- suppressWarnings({
    st_intersection(cat, fl) %>%
      filter(ID != catID, as.numeric(st_length(.)) > 1) %>%
      st_collection_extract("LINESTRING")
  })

  for (x in 1:nrow(o)) {

    to_cut   <- cat[cat$ID == o$ID[x], ]
    to_merge <- cat[cat$ID == o$catID[x], ]

    spliter  <- st_split(to_cut, st_union(fl)) %>%
      st_collection_extract("POLYGON") %>%
      mutate(area = as.numeric(st_area(.))) %>%
      arrange(area)

    oddball <- bind_rows(to_merge, spliter[1, ]) %>%
      st_snap_to_grid(size = 1e-3) %>%
      st_union()

    the_rest <- spliter[2:nrow(spliter), ] %>%
      st_snap_to_grid(size = 1e-3) %>%
      st_union

    pts_fl        <- st_cast(st_geometry(fl[fl$ID == o$ID.1[x],]), "POINT")
    pts_fl1       <- pts_fl[1]
    pts_fl2       <- pts_fl[length(pts_fl)]

    pts_cut       <- st_cast(st_geometry(o[x,]), "POINT")

    magic_pt      <- if(pts_fl1 %in% pts_cut) { pts_fl1} else { pts_fl2 }

    pt_to_move    <- pts_cut[!pts_cut %in% magic_pt, ]

    oddball_pts   <- st_cast(oddball, "POINT")
    id            <- which.min(abs(st_distance(oddball_pts, pt_to_move )))
    oddball_pts[id] <- magic_pt
    oddball_poly  <- st_combine(oddball_pts) %>%
      st_cast("POLYGON")

    the_rest_pts  <- st_cast(the_rest, "POINT")
    id            <- which.min(abs(st_distance(the_rest_pts, pt_to_move)))
    the_rest_pts[id] <- magic_pt
    the_rest_poly <- st_combine(the_rest_pts) %>% st_cast("POLYGON")

    st_geometry(cat[cat$ID == to_cut$ID,])   <- the_rest_poly
    st_geometry(cat[cat$ID == to_merge$ID,])    <- oddball_poly
  }

  merge_fl = left_join(fl, st_drop_geometry(select(cat, -areasqkm)), by = "ID") %>%
    mutate(length = as.numeric(set_units(st_length(.), 'km'))) %>%
    mutate(ID = catID)
  merge_cat = cat %>% mutate(areasqkm = as.numeric(set_units(st_area(.), 'km2')))

  nexus_flow_graph(merge_fl, merge_cat)

}
