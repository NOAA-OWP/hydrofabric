source('workflow/utils.R')
path         <- "workflow/graph_wf/cache/ngen_01a-0.gpkg"
## INGEST
fl  = read_sf(path, "flowline")

# BUILD OUT FLOW GRAPH
net = build_node_net(fl)
# ISOLATE NODES AND NETWORK
nodes = net$node
network = net$network %>% 
  mutate(lengthkm = as.numeric(set_units(st_length(.), "km")))


##### FIND LEVELPATHS WITH HOLES!! ####
# CAST TO sp FOR rgeos SPEED 
SPDF =  as_Spatial(network)
# set IDS, build sp object
rownames(SPDF@data) <- sapply(slot(SPDF, "lines"), function(x) slot(x, "ID"))  
tmp <- rgeos::gLineMerge(SPDF,
                         byid = TRUE,
                         id = network$levelpath)
ids <- data.frame(ID=sapply(slot(tmp, "lines"), function(x) slot(x, "ID")))
rownames(ids)  <- ids$ID

# Cast to SP, line merge geometries, rename default ID column
out = st_as_sf(sp::SpatialLinesDataFrame(tmp,ids)) %>% 
  st_line_merge() %>% 
  rename('levelpath' = ID)

# WHICH ARE MULTILINSTRINGS?
runner = out$levelpath[st_geometry_type(out) == "MULTILINESTRING"]

# FIX NETWORK -------------------------------------------------------------
# i = 6
# i = which(runner == 2238467)
lpID = runner[100]

fill_level_path = function(lpID, network){
    # Base level path
    lp = filter(network, levelpath == lpID) 
    # Find current top and tail
    og_head  = slice_min(lp,hydroseq)
    og_tail  = slice_max(lp,hydroseq)

    h1 = nhdplusTools::get_node(og_head, "start")
    h2 = nhdplusTools::get_node(og_head, "end")
    head = ifelse(lengths(st_intersects(h1, lp)) > 
                    lengths(st_intersects(h2, lp)),
           h2, h1)
    
    t1 = nhdplusTools::get_node(og_tail, "start")
    t2 = nhdplusTools::get_node(og_tail, "end")
    tail = ifelse(lengths(st_intersects(t1, lp)) > 
                    lengths(st_intersects(t2, lp)),
               t2, t1)
    
    candidate = st_filter(network, st_union(st_buffer(lp, 250)))
   
    t = st_line_merge(st_union(candidate)) %>% 
      st_cast("LINESTRING") %>% 
      st_as_sf()

    net = as_sfnetwork(t, directed = FALSE) %>%
      activate("edges") %>%
      filter(!edge_is_multiple()) %>%
      filter(!edge_is_loop()) %>% 
      convert(
        to_spatial_shortest_paths,
        from = head[[1]], to = tail[[1]],
        weights = edge_length()
      )
    
    new_lp = activate(net, "edges") %>% 
      st_as_sf() 
      
    st_anti_filter = function(.x, .y, .predicate = st_intersects) {
      filter(.x, lengths(.predicate(.x, .y)) == 0)
    }
    
    c = build_node_net(candidate)
    
    breaks = st_collection_extract(lwgeom::st_split(candidate, c$node),"LINESTRING") %>% 
      mutate(tmpLength = round(as.numeric(set_units(st_length(.), "km")) / lengthkm, 2),
             tmpLength = ifelse(tmpLength == 1.00, NA, tmpLength),
             newID = paste0(, ifelse(!is.na(tmpLength),tmpLength,""))) 
  
    new_lp_full = st_filter(breaks, new_lp, .predicate = st_covered_by) %>% 
      mutate(levelpath = lpID) %>% 
      mutate(member_comid = paste(member_comid, newID, sep = ","))
    
    new_lp_full$member_comid
  
    frags = st_anti_filter(breaks,  new_lp, .predicate = st_within) %>% 
      st_filter(new_lp, .predicate = st_touches) %>% 
      filter(ID %in% new_lp_full$ID) 

    new_lp_full$member_comid
    frags$member_comid
    
    corrections = list()
    
    if(nrow(frags) > 0){
      for(j in 1:nrow(frags)){
      cand = filter(c$network, ID == frags$ID[j])
      
      non_lp_connection = filter(c$network, to %in% c(cand$from) ) %>% 
        filter(!ID %in% new_lp_full$ID)
      # mapview(non_lp_connection, color = "Red") + cand + new_lp_full
      non_lp_connection$geom = build_flow_line(frags$geom[j],
                                               non_lp_connection$geom)
      
      non_lp_connection$member_comid = paste(non_lp_connection$member_comid,
                                             frags$newID[j], sep = ",")

      
      corrections[[j]] = non_lp_connection
      } 
    } else {
      corrections = NULL
    }

    corr = bind_rows(corrections) 

    new_lp_full = new_lp_full %>% 
      mutate(levelpath = as.numeric(levelpath))

    improvements = bind_rows(new_lp_full, corr) 
    net = bind_rows(improvements, filter(network, !ID %in% improvements$ID))
}
    

full_net = filter(net, levelpath %in% improvements$levelpath) 
cat2     = filter(cat, ID %in% full_net$ID)

out = list(cat = cat2, flowpath = full_net)
saveRDS(out, "data/sample_network.rds")

    ## BREAKS
    
    # ## FRAGS
    # 
    # ## RE-ID
    par(mfrow = c(1,3))
    plot(candidate$geom, col = "red")
    plot(lp$geom, col  = "green", add = TRUE)

    plot(candidate$geom, col = "red")
    plot(lp$geom, col  = "green", add = TRUE)
    plot(head[[1]], col  = "blue", add = TRUE)
    plot(tail[[1]], col  = "purple", add = TRUE)

    plot(candidate$geom, col = "red")
    plot(head[[1]], col  = "blue", add = TRUE)
    plot(tail[[1]], col  = "purple", add = TRUE)
    plot(net, add = TRUE)
}

raw_nhd = 

refactor_hy = function(){}


