purge_flowlines = function(fl, cat){

  new_fl <- fl %>%
    group_by(ID) %>%
    slice_max(streamorde) %>%
    # summarise() %>%
    # slice_min(hydroseq_min, n = 1) %>%
    ungroup() %>%
    multi_to_line() %>%
    mutate(lengthkm = as.numeric(set_units(st_length(.), "km")))

  return(list(fl = new_fl, cat = cat))
}

reduce_size = function(fl, cat, min_size = 3, max_size = 15){

  merge_cat <- left_join(cat, st_drop_geometry(fl), by = "ID")

  to_try   <-  merge_cat %>%
    filter(areasqkm < min_size) %>%
    arrange(areasqkm)

  message("Will try to merge: ", nrow(to_try))

  if(nrow(to_try) > 0){

    for(i in 1:nrow(to_try)){
      this_cat   <- filter(merge_cat, ID == to_try$ID[i])
      toCOMIDs   <- fromJSON(this_cat$toCOMID)
      fromCOMIDs <- fromJSON(this_cat$fromCOMID)

      all <- lapply(1:nrow(merge_cat), function(x) { fromJSON(merge_cat$comids[x])})

      ds_cat <- merge_cat %>%
        filter(grepl(paste(toCOMIDs, collapse = "|"), all)) %>%
        filter(ID != this_cat$ID) %>%
        mutate(type = "DS") %>%
        filter(hydroseq_min == max(hydroseq_min, na.rm =TRUE))

      us_cat <- suppressMessages({ merge_cat %>%
        filter(grepl(paste(fromCOMIDs, collapse = "|"), all)) %>%
        filter(ID != this_cat$ID) %>%
        mutate(type = "US") %>%
        st_filter(this_cat) %>%
        filter(!ID %in% ds_cat$ID)
      })

      # The goal here is to identify the upstream flowpath
      # that contributes to a catchment with an area less then a
      # user defined threshold. To do this the "target" fl is
      # intersected with all "contender" flowlines.
      # The intersection identifies how many flowlines contribute to each junction
      # We remove any contender flowpaths with the same level path as
      # the target flow line as these were handled in our aggregation by 'levelpath'.
      # We then identify the maximum junctions at each contending flowpath and the
      # contenders are ranked my "size" of junction. We prefer smaller junctions

      if(nrow(us_cat) > 0){

        pts = do.call(rbind,
                      list(filter(fl, ID %in% us_cat$ID),
                           filter(fl, ID %in% this_cat$ID))) %>%
              st_intersection() %>%
              st_collection_extract("POINT") %>%
              filter(!levelpathi %in% this_cat$levelpathi) %>%
              group_by(ID) %>%
              slice_max(n.overlaps, n = 1) %>%
              ungroup() %>%
              arrange(n.overlaps) %>%
              select(ID, n.overlaps) %>%
              left_join(st_drop_geometry(select(merge_cat,
                                                areasqkm, ID)), "ID") %>%
              mutate(joined_area = areasqkm + this_cat$areasqkm) %>%
              filter(joined_area < max_size) %>%
              slice(1)

        if(nrow(pts) > 0){
          us_cat = filter(us_cat, ID == pts$ID)#, streamorde <= this_cat$streamorde)
        }   else {
          us_cat = NULL
        }

      } else {

        us_cat = NULL

      }



      topo <- do.call(rbind, list(ds_cat, us_cat))

      if (this_cat$streamorde != 1) {
        # If stream order is not 1
        # then only allow merging with same order stream downstream,
        # Upstreamjoinin
        topo <- filter(topo, streamorde <= this_cat$streamorde) %>%
          mutate(joined_area = areasqkm + this_cat$areasqkm) %>%
          filter(joined_area < max_size)

        if (nrow(topo) > 1) {
          to_merge <- filter(topo, type == "DS")
        } else {
          to_merge <- topo
        }
      } else {
        to_merge <- filter(topo, type == "DS")
      }

      new_area <- (to_merge$areasqkm + this_cat$areasqkm)
      if(length(new_area) != 1){ new_area = 1e9}

      if (new_area < max_size) {
        keep_id <- which(merge_cat$ID == to_merge$ID)
        lp = c(this_cat$levelpathi,to_merge$levelpathi)[which.max(c(this_cat$streamorde,
                                                                    to_merge$streamorde))]

        merge_cat$comids[keep_id]       <- toJSON(c(fromJSON(this_cat$comids),
                                                    fromJSON(to_merge$comids)))
        merge_cat$toCOMID[keep_id]      <- toJSON(c(fromJSON(this_cat$toCOMID),
                                                    fromJSON(to_merge$toCOMID)))
        merge_cat$fromCOMID[keep_id]    <- toJSON(c(fromJSON(this_cat$fromCOMID),
                                                    fromJSON(to_merge$fromCOMID)))
        st_geometry(merge_cat[keep_id,]) <- st_union(st_geometry(this_cat),
                                                    st_geometry(to_merge))
        merge_cat$levelpathi[keep_id]   <- lp

        # Needed for recursive checking!
        merge_cat$areasqkm[keep_id] = new_area

        merge_cat <- filter(merge_cat, ID != this_cat$ID)

        fl        <- mutate(fl, ID = ifelse(ID == this_cat$ID, to_merge$ID, ID))
      }
    }
  }

  merge_fl = left_join(select(fl, ID, hydroseq_min, streamorde),
                       st_drop_geometry(select(merge_cat, -lengthkm, -areasqkm, -hydroseq_min, -streamorde)), by = "ID") %>%
    mutate(lengthkm = as.numeric(set_units(st_length(.), 'km')))

  merge_cat = select(merge_cat, ID) %>%
    mutate(areasqkm = as.numeric(set_units(st_area(.), 'km2'))) %>%
    select(ID, areasqkm)

  return(list(fl = merge_fl, cat = merge_cat ))
}

aggregate_by_size = function(fl, cat, min_size, max_size){
  try1      = reduce_size(fl, cat, min_size, max_size)
  purge_fl  = purge_flowlines(fl = try1$fl, cat = try1$cat)
  #try2      = reduce_size(fl = purge_fl$fl, cat = purge_fl$cat, min_size, max_size)
  #purge_fl2 = purge_flowlines(try2$fl, try2$cat)
  return(purge_fl)
}
