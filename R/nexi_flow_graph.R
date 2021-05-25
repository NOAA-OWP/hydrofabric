find_nexi = function(fl){
  nexus  <- st_set_geometry(fl, find_node(fl, position = NULL)) %>%
    mutate(type = "to")
  inlets  <- st_set_geometry(fl, find_node(fl, position = 1)) %>%
    mutate(type = "from")

  do.call(rbind, list(inlets, nexus)) %>%
    mutate(X = st_coordinates(.)[,1], Y = st_coordinates(.)[,2]) %>%
    mutate(xy = paste(.$X, .$Y)) %>%
    mutate(nodeID = paste0("nex-", group_indices(., factor(xy, levels = unique(xy))))) %>%
    dplyr::select(-xy)
}

nexus_flow_graph = function(fl, cat){

  nexi = find_nexi(fl)

  start = fl %>%
    lwgeom::st_snap_to_grid(size = 1-5) %>%
    lwgeom::st_split(nexi) %>%
    st_collection_extract("LINESTRING") %>%
    mutate(catID = ID, ID = 1:n()) %>%
    mutate(lengthkm = as.numeric(set_units(st_length(.), "km"))) %>%
    select(hydroseq_min, comids, lengthkm, catID, ID)

  nexi = find_nexi(start)

  bounds = nexi %>%
    group_by(nodeID) %>%
    arrange(type) %>%
    filter(n() == 1)  %>%
    mutate(type = ifelse(type == "from", "inlet", "outlet"))

  nexi$nodeID[which(nexi$nodeID == filter(bounds, type == "outlet")$nodeID)] = "nex-0"

  map    <- st_intersects(nexi, start)
  grow   <-  list()

  for(i in 1:nrow(nexi)){

    contributing = start[map[[i]], ]

    grow[[i]] = filter(nexi, ID %in% contributing$ID) %>%
      st_drop_geometry() %>%
      tidyr::pivot_wider(id_cols = ID,
                         names_from = type,
                         values_from = nodeID)
  }


  nodes = do.call(rbind, grow) %>%
    left_join(st_drop_geometry(select(start, ID, catID)), by = 'ID') %>%
    filter(!duplicated(.)) %>%
    arrange(ID) %>%
    mutate(from =
             ifelse(from %in% filter(bounds, type == "inlet")$nodeID,
                    paste0("cat-", bounds$ID[which(bounds$catID %in% catID)]), from),
           to =
             ifelse(to %in% filter(bounds, type == "outlet")$nodeID, "cat-0", to))

  nodes_sf = filter(nexi, nodeID %in% c(nodes$to, nodes$from)) %>%
    dplyr::select(nodeID) %>%
    group_by(nodeID) %>%
    mutate(count = n()) %>%
    distinct(nodeID, .keep_all = TRUE) %>%
    arrange(nodeID)

  edges = left_join(start, select(nodes, -catID), by = "ID") %>%
    multi_to_line()

  return(list(nodes = nodes_sf, edges = edges, cat = cat))
}


