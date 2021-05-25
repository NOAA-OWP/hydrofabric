collapse_junctions = function(graph, og_catchments, threshold = 1, max_size = 15){

  candidates = filter(graph$nodes, count == 2)

  paths_to_collapse = filter(graph$edges, from %in% candidates$nodeID | to %in% candidates$nodeID) %>%
    filter(lengthkm < threshold)

  message("Trying to collapse: ", nrow(paths_to_collapse), " junction")

  for(i in 1:nrow(paths_to_collapse)){

    line      = paths_to_collapse[i,]
    junction2 = filter(candidates, nodeID %in% c(line$to, line$from))
    junction3 = filter(graph$nodes, nodeID %in% c(line$to, line$from), count == 3)

    og_nhd_catchment   = og_catchments[unlist(st_intersects(junction3, og_catchments)),] %>%
      st_buffer(.0001)

    paths        = filter(graph$edges, from %in% c(line$to, line$from) | to %in% c(line$to, line$from))
    path_to_join = filter(paths, ID != line$ID,  from == junction2$nodeID | to == junction2$nodeID)
    all_touching = filter(graph$edges,  from == junction2$nodeID | to == junction2$nodeID |
                            from == junction3$nodeID | to == junction3$nodeID)

    all_cat = filter(graph$cat, ID %in% all_touching$ID)

    st_geometry(graph$edges[which(graph$edges$ID == path_to_join$ID),]) =
      st_union(st_geometry(line), st_geometry(path_to_join)) %>%
      st_line_merge()

    graph$edges = filter(graph$edges, ID != line$ID)

    st_geometry(graph$cat[which(graph$cat$ID == path_to_join$catID),])  =
      st_union(c(st_geometry(og_nhd_catchment),
                 st_geometry(filter(graph$cat, ID == path_to_join$catID))))

    t = st_difference(st_geometry(filter(graph$cat, ID == line$catID)),
                      st_geometry(og_nhd_catchment)) %>%
      st_cast("POLYGON") %>%
      st_as_sf() %>%
      slice_max(st_area(.))

    st_geometry(graph$cat[which(graph$cat$ID == line$catID),])  = st_geometry(t)
  }



  sub_graph = nexus_flow_graph(mutate(graph$edges, ID = catID), graph$cat)

  ######

  #CHECK: Good up to here
  mapview(sub_graph)

  candidates = filter(sub_graph$nodes, count == 2) %>%
    mutate(area = NA, catID = NA)

  for(i in 1:nrow(candidates)){
    tmp  = filter(sub_graph$edges, from %in% candidates$nodeID[i] | to %in% candidates$nodeID[i]) %>%
      left_join(st_drop_geometry(sub_graph$cat), by = "ID" ) %>%
      summarize(a = sum(areasqkm), catID = toJSON(catID))

    candidates$area[i] = tmp$a
    candidates$catID[i] = tmp$catID
  }

  h = filter(candidates, area <= max_size)

  for(i in 1:nrow(h)){

    cats_to_merge = filter(sub_graph$cat, ID %in% fromJSON(h$catID[i])) %>%
      mutate(catID = ID, ID = NULL) #%>%
      #left_join(st_drop_geometry(sub_graph$edges), by = "catID")
    fl_to_merge   = filter(sub_graph$edges, catID %in% cats_to_merge$catID)
#
#     mapview(graph$cat)  +
#       st_geometry(graph$edges[which(graph$edges$ID == lower_fl),]) +
#       st_union(st_geometry(filter(graph$edges, ID %in% fromJSON(h$catID[i]))))

    lower_cat       = slice_min(fl_to_merge, hydroseq_min, n= 1)
    lower_fl        = slice_min(fl_to_merge, hydroseq_min, n= 1)
    upper_cat       = slice_max(fl_to_merge, hydroseq_min, n= 1)
    upper_fl        = slice_max(fl_to_merge, hydroseq_min, n= 1)

    st_geometry(sub_graph$cat[which(sub_graph$cat$ID == lower_cat$catID),]) =
      st_union(st_geometry(cats_to_merge))
    sub_graph$cat  = filter(sub_graph$cat, ID != upper_cat$catID)

    replace_id = which(graph$edges$ID == lower_fl$ID)

    st_geometry(sub_graph$edges[replace_id,]) =
      st_line_merge(st_union(fl_to_merge))

    sub_graph$edges$comids[replace_id] = toJSON(unique(c(fromJSON(lower_fl$comids),  fromJSON(upper_fl$comids))))

    sub_graph$edges  = filter(sub_graph$edge, ID != upper_fl$ID)

    sub_graph$cat   = sub_graph$cat   %>%
      mutate(areasqkm = as.numeric(set_units(st_area(.), "km2")))
    sub_graph$edges = sub_graph$edges %>%
      mutate(lengthkm = as.numeric(set_units(st_length(.), "km")))
  }



  fl = multi_to_line(sub_graph$edges) %>%
       mutate(ID = catID)

  out = nexus_flow_graph(fl = fl, cat = sub_graph$cat)

}

