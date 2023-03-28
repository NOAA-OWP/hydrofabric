subset_hf = function(gpkg,
                     hy_id,
                     comid,
                     poi_id,
                     origin  = NULL,
                     export_gpkg){
  
  
  tmp = select(read_sf(gpkg, "network"), id, toid)
  n   = st_layers(gpkg)$name
  
  crs = st_layers(gpkg)$crs
  
  ids = unique(c(unlist(get_sorted(tmp, outlets = origin))))
  
  db <- DBI::dbConnect(RSQLite::SQLite(), gpkg)
  on.exit(DBI::dbDisconnect(db))
  
  for(j in 1:length(n)){
    hydrofab::hyaggregate_log("INFO", glue("Subsetting: {n[j]} ({j}/{length(n)})"))
    t = tbl(db, n[j]) %>%
      filter(if_any(contains('id'), ~ . %in% ids)) %>%
      collect()
    
    if(!any(is.na(as.character(crs[[j]])))){
      t = st_as_sf(t, crs = crs[[j]])
    }
    
    write_sf(t, export_gpkg, n[j])
  }
  export_gpkg
}



