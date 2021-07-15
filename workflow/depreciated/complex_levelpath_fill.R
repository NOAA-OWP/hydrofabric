function(lpID, network){
  # Base level path
  lp = filter(network, levelpath == lpID)  
  
  # Find current top and tail
  c = build_node_net(lp)
  og_head  = filter(c$network, fromType == "hw")
  og_tail  = filter(c$network, toType == "term")
  
  # find all network elements touching these nodes in 2 degrees
  candidate = st_filter(network, st_filter(network, c$node)) 
  
  # build a new flow network from the candidate flowlines
  c = build_node_net(candidate)
  
  # REMOVE INTERIOR LOOPS
  new_lp = c$network %>% 
    mutate(l = as.numeric(st_length(.))) %>% 
    group_by(node1 = pmin(to, from),
             node2 = pmax(to, from)) %>% 
    slice_min(l, n = 1)
  
  c = build_node_net(new_lp)
  
  # Break candidate FPs at the graph nodes and add new temp ID
  breaks = st_collection_extract(lwgeom::st_split(c$network, c$node),"LINESTRING") %>% 
    mutate(new_id = 1:n())
  
  # Set while loop for dropping 'HW' nodes
  count  = 1e9
  new_lp = breaks
  
  while(count != nrow(new_lp)){
    count = nrow(new_lp)
    c = build_node_net(new_lp)
    new_lp = filter(c$network, fromType != "hw" | ID %in% og_head$ID) %>% 
      filter(toType != "term" | ID %in% og_tail$ID)
  }
  
  # Check if resulting flow path is continous linestring
  cond = !nrow(st_as_sf(st_cast(st_line_merge(st_union(new_lp)), "LINESTRING"))) == 1
  tmp = new_lp
  
  # If not, start droping inflow paths
  while(cond){
    tmp = st_cast(st_line_merge(st_union(tmp)), "LINESTRING") %>%
      st_as_sf() %>%
      mutate(ID = 1:n(), l = as.numeric(st_length(.)))
    
    if(nrow(tmp) > 1){
      tmp = slice_max(tmp, l, n = (nrow(tmp)-1))
      cond = ifelse(nrow(tmp) > 1, TRUE, FALSE)
    } else {
      cond = FALSE
    }
  }
  
  #will erase `y` (tmp) from `x` (new_lp)
  erase_frags = st_difference(new_lp, st_union(st_combine(tmp)))
  
  new_lp = filter(new_lp, !new_id %in% erase_frags$new_id)
  c = build_node_net(new_lp)
  
  # Use resulting new LP to 
  final_lp = c$network %>% 
    mutate(levelpath = as.numeric(lpID))
  
  #mapview(final_lp, color = 'red') + frags
  frags = filter(breaks, !new_id %in% final_lp$new_id & ID %in% final_lp$ID)
  
  if(nrow(frags) > 0){
    extras = st_filter(network, breaks)
    st_geometry_type(extras) %>% table()
    ids = list()
    
    for(j in 1:nrow(frags)){
      # ID the flowpath to merge into
      merge_here = filter(st_filter(extras, frags[j,]), !ID %in% final_lp$ID)
      # Find that ID in extras and merge the geom
      extras$geom[extras$ID == merge_here$ID] <- build_flow_line(merge_here$geom, 
                                                                 frags$geom[j])
      # record the ID
      ids[[j]] = merge_here$ID
    } 
    new_merges = filter(extras, ID %in% c(unlist(ids)))
  } else {
    new_merges = NULL
  }
  
  improvements = bind_rows(final_lp, new_merges)
  bind_rows(improvements, filter(network, !ID %in% improvements$ID))
}