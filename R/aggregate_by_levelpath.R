find_node = function(fl, position){
  tmp = data.frame(st_coordinates(fl))

  if(ncol(tmp) == 4){
    tmp %>%
      group_by(L2) %>%
      slice(ifelse(is.null(position), n(), position)) %>%
      ungroup() %>%
      st_as_sf(coords = c("X", "Y"), crs = st_crs(fl)) %>%
      st_geometry()
  } else {
    tmp %>%
      group_by(L1) %>%
      slice(ifelse(is.null(position), n(), position)) %>%
      ungroup() %>%
      st_as_sf(coords = c("X", "Y"), crs = st_crs(fl)) %>%
      st_geometry()
  }
}

aggregate_by_levelpath = function(fl, cat, ideal_size, max_size, max_length){

fl = left_join(fl, st_drop_geometry(cat), by = "comid")

inlets     <- st_set_geometry(fl, find_node(fl, position = 1)) %>%
  filter(!fromnode %in% tonode)

out        <- list()

for (i in 1:nrow(inlets)) {
  this.lp <- filter(fl, levelpathi == inlets$levelpathi[i]) %>% arrange(-hydroseq)
  values  <- this.lp$areasqkm
  values2 <- this.lp$lengthkm
  indexes <- c()
  count   <- 1

  while (length(values) > 0) {
    v = cumsum(values)
    v = ifelse(v > max_size, v + 1e9, v)

    v2 = cumsum(values2)

    index = min(which.min(abs(v - ideal_size)), sum(v2 < max_length))

    inds    <- 1:index
    values  <- values[-inds]
    values2 <- values2[-inds]
    indexes <- c(indexes, rep(count, length(inds)))
    count   <- count + 1
  }

  all_from  <- lapply(1:nrow(this.lp), function(x){ fromJSON(this.lp$fromCOMID[x]) })
  # all_to  <- lapply(1:nrow(this.lp), function(x){ fromJSON(this.lp$toCOMID[x]) })
  # all_comid  <- lapply(1:nrow(this.lp), function(x){ this.lp$comid[x] })

  out[[i]] <-  suppressMessages({
    group_by(this.lp, ind = indexes) %>%
    summarise(levelpathi  = levelpathi[1],
              hydroseq_min = min(hydroseq),
              streamorde   = max(streamorde),
              toCOMID      = toJSON(unlist(Map(fromJSON, toCOMID))),
              fromCOMID    = toJSON(unique(unlist(all_from))),
              comids       = toJSON(unique(comid)))
  })
}


new_fl <- do.call(rbind, out) %>%
  mutate(lengthkm = as.numeric(set_units(st_length(.), "km")),
         ID = 1:n(), levelpathi,
         ind = NULL) %>%
  multi_to_line()

catchments = list()

suppressWarnings({
  suppressMessages({
    for(i in 1:nrow(new_fl)){
    catchments[[i]] = filter(cat, comid %in% fromJSON(new_fl$comids[i])) %>%
      summarize(ID = new_fl$ID[i])
   }
  })
})

new_cat <- suppressMessages({
do.call(rbind, catchments) %>%
 # multi_to_poly() %>%
  mutate(areasqkm = as.numeric(set_units(st_area(.), "km2"))) %>%

  rmapshaper::ms_simplify(.9)
})

return(list(fl = new_fl, cat = new_cat))
}
