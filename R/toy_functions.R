make_hy = function(flowpaths, catchments, ID = "comid"){
  out = full_join(og_cat, as.data.frame(og_fl), by = ID) %>%
    rename("cat_geom" = geom.x, "fl_geom" = geom.y)
}

realize = function(x, type = "cat"){
  st_set_geometry(x, paste0(type, "_geom"))
}

update_measures = function(hy){
  hy %>%
    activate("cat") %>%
    mutate(areasqkm = set_units(st_area(.), "km2")) %>%
    activate("fl") %>%
    mutate(lengthkm = set_units(st_length(.), "km"))
}


plot_hy = function(hy){
  plot(st_geometry(activate(hy, "cat")), lwd = .5)
  plot(st_geometry(activate(hy, "fl")), lwd = .5, add = TRUE, col = "blue")
}


multi_to_line = function(obj){
  for (i in 1:nrow(obj)) {
    if (st_geometry_type(st_geometry(obj)[i]) == "MULTILINESTRING") {
      st_geometry(obj)[i] <- st_line_merge(st_geometry(obj)[i])
    }
  }
  obj
}


multi_to_poly = function(obj){
  for (i in 1:nrow(obj)) {
    if (st_geometry_type(st_geometry(obj)[i]) == "MULTIPOLYGON") {
      st_geometry(obj)[i] <- st_combine(st_cast(st_geometry(obj)[i], "POLYGON"))
    }
    message(i)
  }
  obj
}

make_graph = function(hy){
  sfnetworks::as_sfnetwork(activate(hy, "fl")) %>%
    sfnetworks::activate("nodes") %>%
    mutate(topo_order = tidygraph::node_topo_order())
}

