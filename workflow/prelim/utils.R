library(sf)
library(dplyr)
library(rmapshaper)
library(mapview)
library(data.table)
library(ggplot2)
library(units)
library(sfnetworks)

flowpaths_to_linestrings = function(fp){
  bool = (st_geometry_type(st_geometry(fp)) == "MULTILINESTRING")
  multis = fp[bool, ]
  if(nrow(multis) > 0){
    st_geometry(multis) = st_line_merge(st_geometry(multis))
  }
  singles = fp[!bool, ]
  bind_rows(multis, singles)
}

validate_network = function(fp, div){
  c("same-number" = nrow(fp) == nrow(div),
    "same-IDs"   = sum(fp$ID %in% div$ID) == nrow(div),
    #"in-flows-intact"    = sum(fp$toID %in% fp$ID | is.na(fp$toID)) == nrow(fp),

  "no-na-levelpaths" = sum(is.na(fp$LevelPathID)) == 0,
  "no-na-fl-ids" = sum(is.na(fp$ID)) == 0,
  "no-na-hydroseq" =sum(is.na(fp$Hydroseq) ) == 0,
  "no-na-member-comid" =sum(is.na(fp$member_COMID)) == 0,
  "no-na-cat-ids" =    sum(is.na(cat$ID)) == 0)

}

build_node_net = function(fl){

  use_sf()
  initGRASS(
    gisBase = grass_path,
    home = tempdir(),
    gisDbase = tempdir(),
    mapset = "PERMANENT",
    location = basename(tempfile()),
    override = TRUE
  )

  writeVECT(
    SDF = fl,
    vname = 'flowpaths',
    v.in.ogr_flags = c("o", 'overwrite')
  )

  # Execute v.clean
  execGRASS(
    cmd = 'v.clean',
    input = 'flowpaths',
    output = 'fp_cleaned',
    tool = 'break',
    flags = c('overwrite', 'c')
  )

  # Read back into R
  edges <- readVECT('fp_cleaned') %>%
    select(ID) %>%
    rename(geometry = geom) %>%
    mutate(edgeID = 1:n())  %>%
    st_set_crs(5070)

  nodes <- edges %>%
    st_coordinates() %>%
    as_tibble() %>%
    rename(edgeID = L1) %>%
    group_by(edgeID) %>%
    slice(c(1, n())) %>%
    ungroup() %>%
    mutate(start_end = rep(c('start', 'end'), times = n()/2)) %>%
    mutate(xy = paste(.$X, .$Y)) %>%
    group_by(xy) %>%
    mutate(n = n())

  nodes$nodeID =  group_indices(nodes)

  nodes = nodes %>%
    ungroup() %>%
    select(-xy)

  edges = edges %>%
    mutate(from =  filter(nodes, start_end == 'start')$nodeID,
           to   =  filter(nodes, start_end == 'end')$nodeID) %>%
    select(ID, edgeID, to, from) %>%
    st_drop_geometry()

  nodes <- nodes %>%
    distinct(nodeID, .keep_all = TRUE) %>%
    select(-c(edgeID, start_end)) %>%
    st_as_sf(coords = c('X', 'Y')) %>%
    st_set_crs(st_crs(fl))

  nodes$n = lengths(st_intersects(nodes, fl))

  return(list(nodes = nodes, edges = edges))

}

cs_group <- function(x, threshold) {
  cumsum <- 0
  group <- 1
  result <- numeric()
  for (i in 1:length(x)) {
    cumsum <- cumsum + x[i]
    if (cumsum > threshold) {
      group <- group + 1
      cumsum <- x[i]
    }
    result = c(result, group)
  }
  return (result)
}

make_plot = function(fl, cat, title, min = 3, ideal = 10, max = 15){
  library(patchwork)

  ggplot(cat, aes(x = areasqkm)) +
    stat_ecdf(color = "red", alpha = .5, lwd = 1) +
    geom_vline(xintercept = min, col = "gray20") +
    geom_vline(xintercept = ideal, col = "darkred") +
    geom_vline(xintercept = max, col = "gray20") +
    xlim(0, max+5) +
    labs(title = title,
         subtitle = paste(
           format(nrow(cat), big.mark=",", scientific=FALSE), "basins/flowlines")) +
    theme_bw() +
    ggplot(fl, aes(x = lengthkm)) +
    stat_ecdf(color = "blue", alpha = .5, lwd = 1) +
    geom_vline(xintercept = .6) +
    geom_vline(xintercept = 5) +
    xlim(0, 15) +
    theme_bw()
}


node_filter = function(fls) {

  oo = as_sfnetwork(fls) %>%
    activate(edges) %>%
    st_as_sf() %>%
    select(to,from,core, l) %>%
    mutate(ID =1:n() )

  inlets  = filter(oo, !from %in% to, core == 0)
  outlets = filter(oo, !to %in% from, core == 0)

  paths   = filter(oo, !ID %in% c(inlets$ID, outlets$ID))

  out = summarize(paths) %>%
    st_line_merge()

  if (st_geometry_type(out) != "MULTILINESTRING") {
    return(out)
  } else {

    non_purge = paths %>% filter(core > 0)

    purgable = paths %>%
      filter(core == 0) %>%
      group_by(to, from) %>%
      slice_max(l) %>%
      ungroup()

    if(nrow(purgable) == 1){
      out = summarize(bind_rows(non_purge, purgable)) %>%
        st_line_merge()

      if (st_geometry_type(out) != "MULTILINESTRING") {
        return(out)
      }
    } else {
      return (NA)
    }

  }

}
#   non_purge = paths %>% filter(core > 0)
#   purgable = paths %>% filter(core == 0)
#
#   if(nrow(purgeable) > 0){
#
#
#   }
#   oo2 = as_sfnetwork(purgable) %>%
#     activate("edges") %>%
#     mutate(weight = edge_length()) %>%
#     convert(to_spatial_shortest_paths, from = 1, to = 8) %>%
#     activate(edges) %>%
#     st_as_sf()
#
#   bind_rows(oo2, non_purge) %>%
#     summarise() %>%
#     st_line_merge() %>%
#     mapview()
#
#   mapview(oo2) + non_purge
#
#   plot(oo2)# %>% plot()
#
#   p = purgable %>%
#     group_by(to) %>%
#     slice_max(l, n = 1) %>%
#     ungroup()
#
#   inlets  = filter(p, !from %in% to)
#   outlets = filter(p, !to %in% from)
#
#   paths   = filter(oo, !ID %in% c(inlets$ID, outlets$ID))
#
#   mapview(p)
#   mapview(non_purge, zcol = "core") + p
#     # %>%
#     # group_by(from) %>%
#     # slice_min(l, n = 1) %>%
#     # ungroup() %>%
#     # group_by(from, to) %>%
#     # slice_max(l) %>%
#     # ungroup() %>%
#     # summarize()
#
#   if (st_geometry_type(oo) == "MULTILINESTRING") {
#     ls = oo %>%
#       st_line_merge() %>%
#       st_cast("LINESTRING")
#     if (nrow(ls) == 2) {
#       nodes <- ls %>%
#         st_coordinates() %>%
#         as_tibble() %>%
#         rename(edgeID = L1) %>%
#         group_by(edgeID) %>%
#         slice(c(1, n())) %>%
#         ungroup() %>%
#         mutate(start_end = rep(c('start', 'end'), times = n() / 2)) %>%
#         mutate(xy = paste(.$X, .$Y)) %>%
#         group_by(xy) %>%
#         st_as_sf(coords = c("X", "Y"), crs = st_crs(fls))
#
#       ind = st_touches(nodes, fls) %>%
#         unlist() %>%
#         table()
#
#       ind = ind[ind > 1]
#
#       if(length(ind) > 0){
#         ind = which.max(ind ) %>%
#           names() %>%
#           as.numeric()
#
#         bind_rows(ls, fls[ind, ]) %>%
#           summarize()
#       } else {
#         ls[which.max(st_length(ls)),]
#       }
#     } else {
#       return(ls)
#     }
#   } else {
#     return(oo)
#   }
# }
#
#
# fls = dangles[6,] %>%
#   ms_explode()
#
#
# mapview(fls)
#
# poss = filter(fl, ID %in% ml$intID)
#
# oo = as_sfnetwork(fls) %>%
#   mutate(topo = tidygraph::node_topo_order())
#
# mapview(activate(oo, 'edges') %>% st_as_sf()) +
# mapview(activate(oo, 'nodes') %>% st_as_sf())
#   st_as_sf()  %>%
#   mutate(l = st_length(.)) %>%
#   select(to, from, l)
#
#
#
#
# st_touches(fls)
#
# mapview(fls)
#
# longest = slice_max(oo, l)
#
# longest = bind_rows(longest,
#                     filter(oo, from == longest$to),
#                     filter(oo, to   == longest$from))
#
# longest = bind_rows(longest,
#                     filter(oo, from == 4),
#                     filter(oo, to   == 3))
#
# mapview(longest) + fls
#
# cont = bind_rows(slice_max(oo, l), filter(oo, from == longest$to),
#                  US = filter(oo, to   == longest$from))
#
# mapview(longest, color = "red") + fls
# which(table(c(oo$from, oo$to)) >2)
#
# which(oo$to == oo$from)
#
# while(oo$from[2] != oo$to[1]){
#   oo = oo[-c(2:(which(oo$from == oo$to[1])-1)),]
# }
#
# st_line_merge(fls) %>%
#   mapview()
#
#
#
# mapview(oo)

intersection_mapping = function(sub, full){

  ll = st_intersects(sub, full)

  data.frame(ID = rep(sub$ID, sapply(ll, length)),
             toID          = rep(sub$toID,sapply(ll, length)),
             IDarea        = rep(sub$areasqkm,sapply(ll, length)),
             IDso          = rep(sub$stream_order,sapply(ll, length)),
             intID         = full$ID[unlist(ll)],
             intSO          = full$stream_order[unlist(ll)],
             intArea = full$areasqkm[unlist(ll)])  %>%
    mutate(toID_match = ifelse(toID == intID, TRUE,FALSE),
           new_area = IDarea + intArea) %>%
    filter(ID != intID) %>%
    group_by(ID) %>%
    mutate(n = 1:n()) %>%
    ungroup()
}

dissolve_network = function(m1, cat_dt, fl_dt, min_area = 3, max_area = 15, filter_fun = NULL){

  if(nrow(m1) != 0){
    for(i in 1:nrow(m1)) {
      og   = filter(cat_dt, ID == m1$toID[i])
      kill = filter(cat_dt, ID == m1$ID[i])
      new_cat_geom = st_union(st_as_sf(rbindlist(list(og,kill)))) %>% st_make_valid()
      new_member_ids = paste(c(og$member_COMID, kill$member_COMID), collapse = ",")
      fl_dt[og$ind, member_COMID := new_member_ids]
      cat_dt[og$ind, geom := st_as_sf(new_cat_geom)]
    }

    new_cat = filter(st_as_sf(cat_dt), !ID %in% m1$ID) %>%
      mutate(areasqkm = as.numeric(st_area(.)/1e6)) %>%
      select(-ind)

    new_fl =  filter(st_as_sf(fl_dt), !ID %in% m1$ID) %>%
      select(-ind, -areasqkm) %>%
      left_join(st_drop_geometry(new_cat), by = "ID")

    cat_dt = data.table(mutate(new_cat, ind = 1:n()))
    fl_dt  = data.table(mutate(new_fl,  ind = 1:n()))

    m1 = filter_fun(new_fl)
  }

  return(list(m1 = m1, fl_dt = fl_dt, cat_dt = cat_dt ))
}

get_nhd_crosswalk <- function(x, catchment_prefix = "catchment_",
                              network_order = NULL, sites = data.frame(local_id = character(0),
                                                                       site_no = character(0))) {

  nhd_crosswalk <- st_drop_geometry(x) %>%
    select(.data$ID, .data$member_COMID) %>%
    mutate(member_COMID = strsplit(.data$member_COMID, ",")) %>%
    unnest(cols = c("member_COMID")) %>%
    mutate(local_id = paste0(catchment_prefix, .data$ID)) %>%
    select(.data$local_id, COMID = .data$member_COMID)

  if(!is.null(network_order)) {
    outlet_comid <- dplyr::mutate(nhd_crosswalk, outlet_COMID = as.integer(.data$COMID)) %>%
      left_join(network_order, by = c("outlet_COMID" = "COMID")) %>%
      group_by(.data$local_id) %>%
      filter(.data$Hydroseq == min(.data$Hydroseq)) %>%
      select(-.data$Hydroseq, -.data$COMID) %>%
      ungroup()

    nhd_crosswalk <- dplyr::left_join(nhd_crosswalk,
                                      outlet_comid,
                                      by = "local_id")

    nhd_crosswalk <- dplyr::left_join(nhd_crosswalk,
                                      sites,
                                      by = "local_id")

    new_names <- unique(nhd_crosswalk$local_id)

    nhd_crosswalk <- lapply(unique(nhd_crosswalk$local_id), function(x, df) {

      df_sub <- df[df$local_id == x, ]

      out <- list(COMID = df_sub$COMID)

      if(any(!is.na(df_sub$site_no))) {
        out$site_no <- unique(df_sub$site_no[!is.na(df_sub$site_no)])
      }

      out$outlet_COMID <- unique(df_sub$outlet_COMID)

      out

    }, df = nhd_crosswalk)

    names(nhd_crosswalk) <- new_names

  }
  nhd_crosswalk
}
  
system.file("shape/nc.shp", package="sf")
  
prep_ngen = function(fl, cat, ID_col = "ID"){
  
  length(fl[[ID_col]])
  sort(fl[[ID_col]]) %>% length()
  sum(duplicated(fl$ID))
  
  id_map = data.frame(
    old_id = sort(fl[[ID_col]]),
    new_id = 1:nrow(fl)
  )
  
    fl_net = as_sfnetwork(fl)
    nodes = st_as_sf(mutate(activate(fl_net,"nodes"), nexID = 1:n()))
    fl    = st_as_sf(mutate(activate(fl_net,"edges"))) %>% 
      rename(old_id = !!ID_col) %>% 
      mutate(ID = id_map$new_id[match(old_id, id_map$old_id)],
             old_id = NULL)
    
    cat    = cat  %>% 
      rename(old_id = !!ID_col) %>% 
      mutate(ID = id_map$new_id[match(old_id, id_map$old_id)],
             old_id = NULL)

    hw   = fl$from[!fl$from %in% fl$to]
    term = fl$to[!fl$to %in% fl$from]
    
    nodes$type = ifelse(nodes$nexID %in% hw, "hw", NA)
    nodes$type = ifelse(nodes$nexID %in% term, "term", nodes$type)
    nodes$type = ifelse(is.na(nodes$type), "nex", nodes$type)
    
    return(list(nex = nodes,
                fl = fl,
                cat  = cat))
}


build_flow_line = function(geom1, geom2){
  g = st_union(c(geom1, geom2))
  
  if(st_geometry_type(g) == "MULTILINESTRING"){
    g = st_line_merge(g)
  }
  
  g
}


build_node_net = function(fl, add_type = TRUE){
  net = suppressWarnings({ sfnetworks::as_sfnetwork(fl)  })
  # network    = st_as_sf(mutate(activate(net,"edges"))) 
  network = net %>% 
  activate('edges') %>% 
  filter(!edge_is_multiple()) %>%
  filter(!edge_is_loop()) %>% 
  st_as_sf()
  
  nodes = st_as_sf(mutate(activate(net,"nodes"), nexID = 1:n()))
    
  
  if(add_type){
    id_map = network %>% 
      # Identify fromID from fromNODE
      left_join(select(st_drop_geometry(network), fromID  = ID, tmpto = to),   
                by = c("from" = "tmpto")) %>% 
      # Identify toID from toNODE
      left_join(select(st_drop_geometry(network), toID    = ID, tmpfrom = from),
                by = c("to" = "tmpfrom"))
    
    hw   = id_map$from[!id_map$from %in% id_map$to]
    term = id_map$to[!id_map$to %in% id_map$from]
    nodes$type = ifelse(nodes$nexID %in% hw, "hw", NA)
    nodes$type = ifelse(nodes$nexID %in% term, "term", nodes$type)
    nodes$type = ifelse(is.na(nodes$type), "nex", nodes$type)
    
    network$fromType = left_join(select(network, from), 
                                 st_drop_geometry(nodes), by = c('from' = "nexID"))$type
    
    network$toType = left_join(select(network, to), 
                               st_drop_geometry(nodes), by = c('to' = "nexID"))$type
  }
  
  list(node = nodes, network = network)
}