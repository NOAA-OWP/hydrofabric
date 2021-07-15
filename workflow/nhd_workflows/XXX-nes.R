source('workflow/utils.R')

# Intermidiates
path       <- "workflow/nhd_workflows/cache/ngen_01a-1.gpkg"
out_path   <- "workflow/nhd_workflows/cache/ngen_01a-2.gpkg"
min_size   <-  3   # in km2
min_length <- .6
#############################################################################

# TO DISSOLVE LEVELPATHS WE NEED:
  #   Headwater make (e.g. from not in to)
  #   areasqkm
  #   lengthkm
  #   to_id
  #   from_id

# RECORD:
#   ngen_id
#   areasqkm
#   lengthkm
#   member_comid


cat = read_sf(path, "catchment") %>% 
  mutate(ind = 1:n()) %>% 
  st_make_valid()

fl  = read_sf(path, "flowline") %>% 
  inner_join(st_drop_geometry(cat), by = "comid") %>% 
  dplyr::select(comid, to_comid, from_comid, hydroseq, 
                levelpath, member_COMID, order, areasqkm,
                lengthkm) %>% 
  mutate(ind = 1:n())


# ints = small basins that are not headwaters...

ints = filter(fl, from_comid != -1) %>% 
  filter(areasqkm <= min_size | lengthkm <= min_length)

# try to merge into to_id, but, if terminal, merge into from_id
t = data.frame(merge = ints$comid, 
               into_these = ifelse(ints$to_comid == -1, 
                                   ints$from_comid, 
                                   ints$to_comid)) 

cat_dt = setDT(as.data.frame(cat))
fl_dt  = setDT(as.data.frame(fl))
remove_ids <- list()

system.time({
  for (i in 1:nrow(t)) {
    small_cat = cat_dt[comid == t$merge[i],]
    into_cat  = cat_dt[comid == t$into_these[i],]

    if(nrow(small_cat) != 0 & nrow(into_cat) != 0 ){
      new_geom = c(small_cat$geom, into_cat$geom) %>%  st_union()
      cat_dt[into_cat$ind, geom := st_as_sf(new_geom)]
      remove_cats[[i]]   = small_cat$comid

      small_fl = fl_dt[comid == t$merge[i],]
      into_fl  = fl_dt[comid == t$into_these[i],]
    #setkey(to_merge, hydroseq)
    
    # Aggregate catchment geometries
      new_geom = c(small_fl$geom, into_fl$geom) %>%  st_union()
    
    fl_dt[into_fl$ind, `:=`( geom       =  new_geom,
        #to_comid   = into_fl$to_comid,
        from_comid  = small_fl$from_comid,
        #hydroseq   = into_fl$hydroseq,
        #order      = into_fl$order,
        #levelpath  = into_fl$levelpath,
        member_COMID = paste(c(into_fl$member_COMID, small_fl$member_COMID), collapse = ",")
      )]
    
    remove_ids[[i]]   = t$merge[i]
    }
  }
})
  
new_cat = filter(st_as_sf(cat_dt), !comid %in% unlist(remove_ids)) %>%
  mutate(areasqkm = as.numeric(st_area(.)/1e6)) %>%
  dplyr::select(comid, areasqkm)

new_fl =  filter(st_as_sf(fl_dt), !comid %in% unlist(remove_ids)) %>%
  mutate(lengthkm = as.numeric(st_length(.)/1e3)) %>%
  dplyr::select(-ind)

#--------------------------------------------------
#############################################################################

make_plot(new_fl, new_cat, "LevelPath Dissolve HF") %>%
  ggsave(
    filename = "workflow/nhd_workflows/cache/img/03-ngen-interior-dissolve.png",
    units = "in",
    height = 4,
    width =  8
  )

message("Dropped: ", nrow(fl) - nrow(new_fl), " features")

unlink(out_path)
write_sf(new_cat, out_path, "catchment")
write_sf(new_fl, out_path, "flowline")



