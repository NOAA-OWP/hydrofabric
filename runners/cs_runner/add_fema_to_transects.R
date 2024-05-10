library(dplyr) 
library(lwgeom)

# generate the flowlines layer for the final cross_sections_<vpu>.gpkg for each vpu
source("runners/cs_runner/config.r")

# transect bucket prefix
transects_prefix <- paste0(s3_bucket, version_prefix, "/3d/transects/")

# paths to nextgen datasets and model attribute parquet files
nextgen_files   <- list.files(nextgen_dir, full.names = FALSE)
# fema_files      <- list.files(fema_fgb_path, full.names = FALSE)
# fema_bb_files   <- list.files(fema_fgb_bb_path, full.names = FALSE)
transects_files <- list.files(transects_dir, full.names = FALSE)
transects_files <- transects_files[!grepl("updated", transects_files)]


net_source <- "hydrofabric3d"

# ensure the files are in the same order and matched up by vpu
path_df <- align_files_by_vpu(
  x    = nextgen_files,
  y    = transects_files,
  base = base_dir
)

path_df

# loop over each vpu and generate cross sections, then save locally and upload to s3 bucket
# for(i in 1:nrow(path_df)) {
  
  i = 8
  
  # nextgen file and full path
  nextgen_file <- path_df$x[i]
  nextgen_path <- paste0(nextgen_dir, nextgen_file)
  
  transect_file <- path_df$y[i]
  transect_path <- paste0(transects_dir, transect_file)
  
  vpu <- path_df$vpu[i]
  transect_path
  
  # # model attributes file and full path
  # model_attr_file <- path_df$y[i]
  # model_attr_path <- paste0(model_attr_dir, model_attr_file)
  
  message("creating vpu ", path_df$vpu[i], "\n - transects: ", transect_file, "\n - flowpaths: '", nextgen_file, "'")
  # message("creating vpu ", path_df$vpu[i], " transects:\n - flowpaths: '", nextgen_file, "'\n - model attributes: '", model_attr_file, "'")
  
  # read in nextgen data
  flines <- sf::read_sf(nextgen_path, layer = "flowpaths")

  fema_vpu_dir <- paste0(FEMA_VPU_SUBFOLDERS[grepl(paste0("VPU_", vpu), basename(FEMA_VPU_SUBFOLDERS))], "/merged")
  fema_vpu_dir
  vpu_fema_files <- list.files(fema_vpu_dir, full.names = TRUE)
  # vpu_fema_file1 <- vpu_fema_files[grepl("_union.gpkg", vpu_fema_files)]
  vpu_fema_file <- vpu_fema_files[grepl(paste0(vpu, ".gpkg"), vpu_fema_files)]
  vpu_fema_file
  # fema1 <- sf::read_sf(vpu_fema_file1)
  # fema2 <- sf::read_sf(vpu_fema_filev2)  
  fema <- sf::read_sf(vpu_fema_file)
  # fema %>%  mapview::npts() 
  # fema2 %>%  mapview::npts() 
  
  fema <- 
    fema %>%
    rmapshaper::ms_simplify(keep = 0.10)
   
  fema %>%  mapview::npts()
  
  transects <- sf::read_sf(transect_path)
  
  transects_geos <- geos::as_geos_geometry(transects)
  fema_geos      <- geos::as_geos_geometry(fema)     

  length(fema_geos)
  
  # fema_geos  
  transects_fema_matrix <- geos::geos_intersects_matrix(transects_geos, fema_geos) 
  fema_transects_matrix <- geos::geos_intersects_matrix(fema_geos, transects_geos) 
  
  # transects_fema_matrix
  # get the polygons that have atleast 1 intersection with the 'lines'
  # transects_fema_matrix[[557]]
  # fema_tmp <- fema[transects_fema_matrix[[557]], ]
  # trans_tmp <- transects[557, ]
  
 # fema_dissolve <-  rmapshaper::ms_dissolve(fema_tmp)
 # fema_simple <- rmapshaper::ms_simplify(fema_tmp, keep = 1) %>% 
   # rmapshaper::ms_dissolve(field = "state")

  # mapview::mapview(trans_tmp, color = "green") + 
  #   mapview::mapview(fema_tmp[1, ], col.region = "red") + 
  #   mapview::mapview(fema_tmp[2, ],  col.region = "dodgerblue") + 
  #   mapview::mapview(fema_simple,  col.region = "yellow")

  fema_transects_matrix
  
  trans_fema    <- transects[lengths(transects_fema_matrix) != 0, ]
  # fema_polygons <- fema[lengths(fema_transects_matrix) != 0, ]
  fema_polygons <- fema_geos[lengths(fema_transects_matrix) != 0]
  
  
  fema_lines <- 
    fema_polygons %>% 
    sf::st_as_sf() %>% 
    sf::st_cast("MULTILINESTRING") %>% 
    geos::as_geos_geometry() %>% 
    geos::geos_simplify_preserve_topology(25)
  
  # length(fema_geos)
  # nrow(fema)
  # fema_geos     <- geos::as_geos_geometry(fema_polygons)
  
  tmp_trans <- trans_fema[1:5000, ]
  
  # tmp_centroid <- sf::st_centroid(tmp_trans)
  # sf::st_segmentize()
  # split_trans <- lwgeom::st_split(tmp_trans, tmp_centroid)
  # # tmp_trans
  max_extension_distance <- 2500
  
  # unlist(tmp_trans$geom)
  # unlist(sf::st_segmentize(tmp_trans, 4)$geom)
  min_segmentation <- min(tmp_trans$cs_lengthm %/% 2)
  
  # which.min(tmp_trans$cs_lengthm %/% 2)
  
  segmented_trans <- sf::st_segmentize(tmp_trans, min_segmentation)
  
  # unlist(segmented_trans$geom)
  unique(lengths(segmented_trans$geom))
  length(lengths(segmented_trans$geom))
  lengths(segmented_trans$geom)  
  
  # left_trans   <- lwgeom::st_linesubstring(segmented_trans, 0, 0.50)
  # right_trans  <- lwgeom::st_linesubstring(segmented_trans, 0.50, 1)
  
  # mapview::mapview(left_trans, col.regions = "dodgerblue") +
    # mapview::mapview(tmp_trans, color = "red") +
    #   mapview::mapview(tmp_trans[42, ], color = "yellow") +
    # mapview::mapview(right_trans, color = "dodgerblue") +
    # mapview::mapview(left_trans, color = "green")
  
  
  left_trans <- 
    segmented_trans %>% 
    lwgeom::st_linesubstring(0, 0.50) %>% 
    dplyr::mutate(
      partition         = "left",
      partition_lengthm = as.numeric(sf::st_length(geom))
    ) %>% 
    hydrofabric3D::add_tmp_id() %>% 
    dplyr::select(tmp_id, hy_id, cs_source, cs_id, cs_measure,  
                  cs_lengthm, is_extended, partition, partition_lengthm, geom)
  
  # Find the distances from the right side of transect lines 
  right_trans <- 
    segmented_trans %>% 
    lwgeom::st_linesubstring(0.50, 1) %>% 
    dplyr::mutate(
      partition         = "right",
      partition_lengthm = as.numeric(sf::st_length(geom))
    ) %>% 
    hydrofabric3D::add_tmp_id() %>%
    dplyr::select(tmp_id, hy_id, cs_source, cs_id, cs_measure,  
                  cs_lengthm, is_extended, partition, partition_lengthm, geom)
 
  # convert the transect geometries to geos types
  # get the fema polygon indices for the transect halves that are completely within a fema polygon
  # add the fema polygons index as a column to the transect dataframes
  left_trans_geos   <- geos::as_geos_geometry(left_trans)
  right_trans_geos  <- geos::as_geos_geometry(right_trans)
  
  left_within_matrix  <- geos::geos_within_matrix(left_trans_geos, fema_polygons)
  right_within_matrix <- geos::geos_within_matrix(right_trans_geos, fema_polygons)
  
  left_within_vect    <- lapply(left_within_matrix, function(i) { if(length(i) > 0) { c(i) } else { c(NA) } })
  right_within_vect   <- lapply(right_within_matrix, function(i) { if(length(i) > 0) { c(i) } else { c(NA) } })
  
  # add the fema polygon indexes as columns
  left_trans$left_fema_index    <- left_within_vect
  right_trans$right_fema_index  <- right_within_vect
  
  # add boolean columns whether the transect is fully within the FEMA polygons
  left_trans <- 
    left_trans %>% 
    dplyr::mutate(
      left_is_within_fema = dplyr::case_when(
        !is.na(left_fema_index) ~ TRUE,
        TRUE                    ~ FALSE
      )
    ) %>% 
    dplyr::select(tmp_id, hy_id, cs_source, cs_id, cs_measure,  
                  cs_lengthm, is_extended, partition, partition_lengthm,
                  left_fema_index, left_is_within_fema, 
                  geom
    )
   
  right_trans <- 
    right_trans %>% 
    dplyr::mutate(
      right_is_within_fema = dplyr::case_when(
        !is.na(right_fema_index) ~ TRUE,
        TRUE               ~ FALSE
      )
    ) %>% 
    dplyr::select(tmp_id, hy_id, cs_source, cs_id, cs_measure,  
                  cs_lengthm, is_extended, partition, partition_lengthm,
                  right_fema_index, right_is_within_fema, 
                  geom
    ) 
 # left_trans
 #  
 # sf::st_within(left_trans, fema_polygons)
 # geos::as_geos_geometry(fema_polygons)
 # geos::as_geos_geometry(left_trans)
 # geos::geos_within(geos::as_geos_geometry(left_trans), geos::as_geos_geometry(fema_polygons))
  
  left_trans_geos     <- geos::as_geos_geometry(left_trans)
  left_within_matrix  <- geos::geos_within_matrix(left_trans_geos, fema_polygons)

  left_within_vect    <- lapply(left_within_matrix, function(i) {
                          if(length(i) > 0) { c(i) } else { c(NA) } }
                          )
  left_trans
  
  
  left_trans <- 
    left_trans %>% 
    dplyr::mutate(
      left_is_within_fema = dplyr::case_when(
        !is.na(left_fema_index) ~ TRUE,
        TRUE               ~ FALSE
      )
    ) %>% 
    dplyr::select(tmp_id, hy_id, cs_source, cs_id, cs_measure,  
                  cs_lengthm, is_extended, partition, partition_lengthm,
                  left_fema_index, left_is_within_fema, 
                  geom
                  )
  left_trans
  # left_trans$geom2 <- left_trans_geos 
  
  # left_trans <- 
  #   left_trans %>% 
  #   sf::st_drop_geometry() %>%
  #   dplyr::mutate(
  #     geom = geos::as_geos_geometry(left_trans_geos)
  #     # geom2 = left_trans_geos
  #   )
  
  left_trans_geos
  
  # ----------------------------------------------------------------------------------------------------------------
  # Loop over every left and right halfs of transects and 
  # if they are fully within FEMA polygons, get the minimum extension distance required for the transect to meet the FEMA polygon boundary
  # ----------------------------------------------------------------------------------------------------------------
  
  left_ids          <- left_trans$tmp_id
  
  left_fema_indexes <- left_trans$left_fema_index
  left_fema_bool    <- left_trans$left_is_within_fema
  
  # preallocate vector that stores the extension. distances
  left_extension_dists <- vctrs::vec_c(rep(0, length(left_ids)))
  
  # all_equal_length_vects <- all(length(left_ids) == length(left_fema_indexes) && length(left_ids) == length(left_fema_bool))
  # 1:length(left_ids)
  extension_count = 0
  
  for(i in 1:length(left_ids)) {
    # i = 1
    tmp_id                 <- left_ids[i]
    is_within_fema_polygon <- left_fema_bool[i]
    fema_index             <- left_fema_indexes[i]

    message("Transect: '", tmp_id, "' - (", i, ")")
    # if(is_within_fema_polygon) {
    #   break
    # }
    # fema_index <- left_trans$left_fema_index[i]
    # is_within_fema_polygon = ifelse(!is.na(left_fema_index), TRUE, FALSE)
    
    if(is_within_fema_polygon) {
      
      message("- Left side of transect intersects with FEMA")
      message("\t > FEMA index: ", fema_index)
      extension_count = extension_count + 1
      message("\t > extension_count: ", extension_count)
      
      trans_geom <- left_trans_geos[i]
      index_vect <- sort(unlist(fema_index))
      
      # mapview::mapview(sf::st_as_sf(left_trans_geos[i]), color = "red") + sf::st_as_sf(fema_lines[index_vect])
      
      # fema_lines[index_vect]
      
      dist_to_fema <- hydrofabric3D:::geos_bs_distance(
                          distances    = 1:5000,
                          line         = trans_geom,
                          geoms_to_cut = fema_lines[index_vect],
                          direction    = "both"
                        )
      
      left_extension_dists[i] <- dist_to_fema
      
    }
    message()
  }
  
  left_trans_geos
  vctrs::vec_c(
    vctrs::vec_c(
      left_trans$hy_id
  ),
  vctrs::vec_cast(left_trans$cs_id)
  )
  
 #  
 #  # mapview::mapview(left_trans[477, ]) + fema_polygons[1888, ]
 # left_trans <- 
 #    left_trans %>%
 #    dplyr::mutate(
 #      fema_index = unlist(sf::st_within(., fema_polygons))
 #    ) %>% 
 #    dplyr::relocate(fema_index)
 #  
  # left_trans_geos <- geos::as_geos_geometry(left_trans)
  
  sort(na.omit(unlist(unique(left_trans$fema_index))))
  fema_polygons[na.omit(unlist(unique(left_trans$fema_index)))]
  
  fema_polygons %>% na.omit()
  
  # NOTE: sorting the fema polygon indices (not sure if necessary)
  left_fema <- fema_polygons[sort(na.omit(unlist(unique(left_trans$fema_index))))]
  
  left_fema
  
  left_fema %>% plot()
  
  geos::geos_make_linestring(geom = left_fema)
  left_trans$fema_index
  sort(na.omit(unlist(unique(left_trans$fema_index))))
  
  # left_fema <- fema_polygons[na.omit(unlist(unique(left_trans$fema_index)))]
  left_fema_lines <- 
    left_fema %>% 
    sf::st_as_sf() %>% 
    sf::st_cast("MULTILINESTRING") %>% 
    geos::as_geos_geometry() %>% 
    geos::geos_simplify_preserve_topology(25)
  
  left_fema %>% 
    sf::st_as_sf() %>% 
    mapview::npts()
  
  geos::geos_simplify_preserve_topology(left_fema_lines, 50) %>% 
    geos::geos_num_coordinates() %>%
    sum()
  
  geos::geos_num_coordinates(left_fema_lines) %>% sum()
  
  geos::geos_simplify_preserve_topology(left_fema_lines, 1) %>% 
    .[3] %>% 
    plot()
  left_fema_lines[3] %>% plot()
  
  
  left_extension_dists <- lapply(seq_along(left_trans_geos), function(i) {
    
    hydrofabric3D:::geos_bs_distance(
      distances    = 1:2000,
      line         = left_trans_geos[i],
      geoms_to_cut = left_fema_lines,
      direction    = "head"
    )
  }) %>% 
    unlist()
  
  left_trans$left_head_extension_dist <- left_extension_dists
  
  # Find the distances from the right side of transect lines 
  right_trans <- 
    segmented_trans %>% 
    lwgeom::st_linesubstring(0.50, 1) %>% 
    dplyr::mutate(
      partition         = "right",
      partition_lengthm = as.numeric(sf::st_length(geom))
    ) %>% 
    dplyr::select(hy_id, cs_source, cs_id, cs_measure,  
                  cs_lengthm, is_extended, partition, partition_lengthm, geom)
  right_trans_geos <- geos::as_geos_geometry(right_trans)
  right_within_matrix <- geos::geos_within_matrix(right_trans_geos, geos::as_geos_geometry(fema_polygons))
 
  right_within_vect <- lapply(right_within_matrix, function(i) {
                      if(length(i) > 0) { c(i) } else { c(NA_real_) } }
                      )
 right_trans$fema_index <- right_within_vect
# 
#  right_trans <- 
#     right_trans %>%
#     dplyr::mutate(
#       fema_index = unlist(sf::st_within(., fema_polygons))
#     ) %>% 
#     dplyr::relocate(fema_index)
#   
  
  right_fema <- fema_polygons[unique(right_trans$fema_index),] 
  right_fema_lines <- 
    right_fema %>% 
    sf::st_cast("LINESTRING")
   
  right_extension_dists <- lapply(seq_along(right_trans_geos), function(i) {
    
    hydrofabric3D:::geos_bs_distance(
      distances    = 1:2000,
      line         = right_trans_geos[i],
      geoms_to_cut = right_fema_lines,
      direction    = "tail"
    )
  }) %>% 
    unlist()
  
  right_trans$right_tail_extension_dists <- right_extension_dists
  right_trans$right_tail_extension_dists
  unlist(sf::st_within(left_trans, fema_polygons))
  
  left_trans <- 
    left_trans %>%
    dplyr::mutate(
      fema_index = unlist(sf::st_within(., fema_polygons))
    ) %>% 
  dplyr::relocate(fema_index)
  
  left_fema <- fema_polygons[unique(left_trans$fema_index),]
  left_fema_lines <- 
    left_fema %>% 
    sf::st_cast("LINESTRING")
  
  geos::as_geos_geometry(left_fema)
  mapview::mapview(left_fema, col.regions = "dodgerblue") +
    mapview::mapview(left_ls, color = "yellow") +
    mapview::mapview(left_trans, color = "red") +
    mapview::mapview(right_trans, color = "green")
  
  left_trans_geos <- geos::as_geos_geometry(left_trans)
  left_trans_geos %>% plot()
  left_trans_geos[2]
  
  length(left_trans_geos)
  
  left_extension_dists <- lapply(seq_along(left_trans_geos), function(i) {
    
   hydrofabric3D:::geos_bs_distance(
      distances    = 1:2000,
      line         = left_trans_geos[i],
      geoms_to_cut = left_fema_lines,
      direction    = "head"
    )
  }) %>% 
    unlist()
  
  left_extensions  <- geos::geos_empty()
  
  for (i in 1:length(left_extension_dists)) {
    dist = left_extension_dists[i]
    geos_line <- left_trans_geos[i]
    message(glue::glue("i: {i}\ndist: {dist}"))
    
    extended <- hydrofabric3D::geos_extend_line(
      geos_line,
      dist,
      "head"
    )
    left_extensions <- vctrs::vec_c(left_extensions, extended)
    
  }
  # index for only valid transects
  # is_valid <- !geos::geos_is_empty(left_extensions)
  left_extensions <- left_extensions[!geos::geos_is_empty(left_extensions)]
  # !geos::geos_is_empty(left_extensions)
  
  new_left_trans <- 
    left_trans %>% 
    sf::st_drop_geometry() %>% 
    dplyr::mutate(
      geom = sf::st_as_sfc(left_extensions)
    ) %>% 
    sf::st_as_sf()
  # geos::sf
  mapview::mapview(left_fema, col.regions = "dodgerblue") +
    mapview::mapview(left_ls, color = "yellow") +
    mapview::mapview(left_trans, color = "red") +
    mapview::mapview(new_left_trans, color = "green")
  mapply(function(geom, dist) {
      hydrofabric3D::geos_extend_line(geom, dist, "head") 
    },
  left_trans_geos,
  left_extension_dists
  )
  
  left_extensions <- lapply(seq_along(left_trans_geos), function(i) {
      
    extend_dist <- hydrofabric3D:::geos_bs_distance(
        distances    = 1:2000,
        line         = left_trans_geos[i],
        geoms_to_cut = left_fema_lines,
        direction    = "head"
      )
      
    hydrofabric3D::geos_extend_line(left_trans_geos[i], extend_dist, "head") 
    
    })
  # unlist(left_extensions)
  
  left_extension_dists <- lapply(seq_along(left_trans_geos), function(i) {
    hydrofabric3D:::geos_bs_distance(
      distances    = 1:2000,
      line         = left_trans_geos[i],
      geoms_to_cut = left_fema_lines,
      direction    = "head"
    )
  }) 
  
  distance_to_extend <- 
    hydrofabric3D:::geos_bs_distance(
      distances    = 1:1500,
      line         = left_trans_geos[1],
      geoms_to_cut = left_fema_lines,
      direction    = "head"
    )
  
  extended <- hydrofabric3D::geos_extend_line(left_trans_geos[1], distance_to_extend, "head") %>% 
    sf::st_as_sf()
  mapview::mapview(left_fema, col.regions = "dodgerblue") +
    mapview::mapview(left_ls, color = "yellow") +
    mapview::mapview(left_trans, color = "red") +
    mapview::mapview(extended, color = "green")
    # mapview::mapview(right_trans, color = "green")
  left_trans$geom %>% sf::st_length()
  plot(segmented_trans$geom, col = "red", lwd =5)
  plot(left_trans$geom, col = "green", lwd=3, add = TRUE)
  plot(right_trans$geom, col = "blue", lwd=3, add = TRUE)
  unlist(left_trans$geom)
   unlist(right_trans$geom)
  
  unlist(tmp_trans$geom )
  unlist(split_trans$geom )
  split_trans %>% 
    sf::st_collection_extract("LINESTRING")
  split_trans$geom
  
  tmp_trans$geom
  nngeo::st_segments(tmp_trans) %>% 
    .$result %>% 
    plot()
  
  mapview::mapview(tmp_trans) + tmp_centroid
  geos::geos_clip_by_rect()
  
  transects_with_
  
  lengths(transects_fema_matrix)
  mapview::mapview(transects_with_fema, color = "green") + fema_with_transects
  unique(hydrofabric3D::add_tmp_id(transects)$tmp_id)[1:30]
  transects  %>% 
    hydrofabric3D::add_tmp_id(transects) %>% 
    .$tmp_id %>% 
    unique() %>% .[1:30]
  trans_subset <- 
    transects %>% 
    hydrofabric3D::add_tmp_id() %>% 
    dplyr::filter(tmp_id %in%  unique(hydrofabric3D::add_tmp_id(transects)$tmp_id)[1:30])
  
  fema_subset <- 
    fema %>% 
    dplyr::filter(fema_id == "1268")

  extended <- hydrofabric3D:::extend_by_length(trans_subset,    rep(500, nrow(trans_subset)))
  extended
  clipped_trans <- rmapshaper::ms_clip(extended, fema)
  
  rmapshaper::ms_clip(extended, fema_subset)
  mapview::mapview(trans_subset, color = "red") + 
    mapview::mapview(extended, color = "yellow") +
    # mapview::mapview(  sf::st_difference(extended, fema_subset), color = "green") +
    mapview::mapview(clipped_trans,  color = "green") +
    fema 
  
   rep(50, nrow(trans_subset))
  extended <- hydrofabric3D:::extend_by_length(trans_subset,    rep(50, nrow(trans_subset)))
  
  
  hydrofabric3D::geos_extend_line(trans_subset, 50) %>% 
    sf::st_as_sf() %>% mapview::mapview()
  
  
  
  
  
  
  
  
  
  
  
  
  
 
  
  
       
  
     
  
     
  
     
  
    
  
    