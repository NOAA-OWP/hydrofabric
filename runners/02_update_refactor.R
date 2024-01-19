source("runners/config.R")

tmp = filter(pipeline, !is.na(corrected_refactor))

if (nrow(tmp) > 0) {

  for (i in 1:nrow(tmp)) {
    
    gpkg = tmp$refactored_gpkg[i]
    message("Fixing: ", gpkg)
    
# 1. Remove Flowlines --------------------------------------------------------
    
    rm = filter(remove_fp_ids, VPU %in% tmp$vpus[i]) %>% 
      tidyr::separate_longer_delim(id, delim = ",")
    
    fps = read_sf(gpkg, "refactored_flowpaths") %>%
      filter(!ID %in% rm$id)

# 3. Mainstem Update ---------------------------------------------------------

    ram = filter(reassign_mainstem, VPU %in% tmp$vpus[i])
    mod = filter(fps, ID %in%  ram$id)
    
    mod$LevelPathID = ram$mainstem[match(mod$ID, ram$id)]
    
    fps = bind_rows(filter(fps, !ID %in% ram$id), mod)
    
# 4. Redigitize --------------------------------------------------------------
    rd =  filter(redigitize_flowpaths, VPU %in% tmp$vpus[i])
    mod = filter(fps, ID %in%  rd$id)
    mod = st_reverse(mod)
    fps = bind_rows(filter(fps, !ID %in% rd$id), mod)
    
    fps$LENGTHKM = add_lengthkm(fps)
    
# Merge Flowlines ---------------------------------------------------------

    mfps =  filter(flowlines_to_merge, VPU %in% tmp$vpus[i]) %>% 
      tidyr::separate_longer_delim(to_merge, delim = ",")
    
    u = unique(mfps$id)
    ll = list()
    
    if(length(u) > 0){
      for(k in 1:length(u)){
        tmp_mfps = filter(mfps, id == u[k])
        mod = filter(fps, ID %in%  tmp_mfps$to_merge)
        
        tmp_fl = mod %>% 
          dplyr::summarise(ID = as.numeric(tmp_mfps$id[1]), 
                           toID = as.numeric(tmp_mfps$toid[1]),
                           TotDASqKM = max(TotDASqKM),
                           member_COMID = paste(member_COMID, collapse = ","),
                           LevelPathID = LevelPathID[1], 
                           refactor_ID = refactor_ID[1])
        
        st_geometry(tmp_fl) = st_line_merge(st_combine(mod))
        
        tmp_fl$LENGTHKM = add_lengthkm(tmp_fl)
        ll[[k]] = tmp_fl
      }
      
      new = bind_rows(ll)
      fps = bind_rows(filter(fps, !ID %in% mfps$to_merge), new)
      
      fps$LENGTHKM = add_lengthkm(fps)
    }
    
    
# Snap Nodes --------------------------------------------------------------
      
  sn = filter(snap_nodes, VPU %in% tmp$vpus[i])
    
  if(nrow(sn) > 0){
    for(i in 1:nrow(sn)){
      new_end_node = get_node(filter(fps, ID == sn$start_node_of[i]), "start")
      old_line = filter(fps, ID == sn$end_node_of[i])
      
      ll = rbind(st_coordinates(old_line)[, c('X', 'Y')], 
                 st_coordinates(new_end_node)) %>% 
        st_linestring() %>% 
        st_sfc(crs = st_crs(old_line))
      
      st_geometry(old_line) <- st_geometry(ll)
      
      fps = bind_rows(filter(fps, !ID %in% old_line$ID), old_line)
    }
  }
     
# 2. Topo Update -------------------------------------------------------------
    
    ut = filter(update_topo, VPU %in% tmp$vpus[i])
    mod = filter(fps, ID %in% ut$id)
    
    mod$toID = ut$toid[match(mod$ID, ut$id)]
    
    fps = bind_rows(filter(fps, !ID %in% ut$id), mod)
     
    write_sf(fps, tmp$corrected_refactor[i], "refactored_flowpaths")


  

# 4. Merge Divides --------------------------------------------------------------
        
    div = read_sf(gpkg, "refactored_divides")
    
    if(tmp$vpus[i] %in% divides_to_merge$VPU){
     
      map = divides_to_merge %>% 
        filter(VPU %in% tmp$vpus[i]) %>% 
        tidyr::separate_longer_delim(to_merge, delim = ",")
      
      u = unique(map$id)
      
      ll = list()
      
      for(j in 1:length(u)){
        map2 = filter(map, id %in% u[j])
        
        mod = filter(div, ID %in% map2$to_merge) %>% 
          mutate(ID = as.numeric(u[j]))
    
        ll[[j]] = union_polygons(mod, "ID") %>% 
          mutate(member_COMID =  paste(mod$member_COMID, collapse = ","),
                 areasqkm = add_areasqkm(.),
                 rpu = mod$rpu[1])
        
      }
      
      div = bind_rows(rename_geometry(filter(div, !ID %in% map$to_merge), "geometry"), 
                rename_geometry(bind_rows(ll), "geometry"))
      
    }
    
# 4. Remove Divides --------------------------------------------------------------
    
    
    if(tmp$vpus[i] %in% remove_divide_ids$VPU){
      
      map = remove_divide_ids %>% 
        filter(VPU %in% tmp$vpus[i]) %>% 
        tidyr::separate_longer_delim(to_merge, delim = ",")
      
      u = unique(map$id)
      
      ll = list()
      
      for(j in 1:length(u)){
        map2 = filter(map, id %in% u[j])
        
        mod = filter(div, ID %in% map2$to_merge) %>% 
          mutate(ID = as.numeric(u[j]))
        
        ll[[j]] = union_polygons(mod, "ID") %>% 
          mutate(member_COMID =  paste(mod$member_COMID, collapse = ","),
                 areasqkm = add_areasqkm(.),
                 rpu = mod$rpu[1])
        
      }
      
      div = bind_rows(rename_geometry(filter(div, !ID %in% map$to_merge), "geometry"), 
                      rename_geometry(bind_rows(ll), "geometry"))
      
    }
    
    
    write_sf(div, tmp$corrected_refactor[i], "refactored_divides")
    
    
  
  }
}
