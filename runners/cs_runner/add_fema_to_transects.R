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
  transects <- sf::read_sf(transect_path)
  
  # fema %>%  mapview::npts() 
  # fema2 %>%  mapview::npts() 
  
  # Give a set of transecct linestrings and poylgons and get the minimum distance to extend each transect line (from both directions, to try and reach the edge of a "polygons")
  # internal function for extending transect lines out to FEMA 100 year flood plain polygons
  # transect_lines, set of Sf linestrigns to extend (only if the transect lines are ENTIRELLY within a polygons)
  # polygons, set of sf polygons that transect lines should be exteneded 
  # max_extension_distance numeric, maximum distance (meters) to extend a transect line in either direction to try and intersect one of the "polygons"
  get_transect_extension_distances_to_polygons <- function(transect_lines, polygons, max_extension_distance) {
    
    # transect_lines <- transects
    # polygons <- fema
    # max_extension_distance <- 2500
    
    # keep 10% of the original points for speed
    polygons <- rmapshaper::ms_simplify(polygons, keep = 0.10)
     
    # fema %>%  mapview::npts()
    
    # transects <- sf::read_sf(transect_path)
    
    transects_geos  <- geos::as_geos_geometry(transects)
    polygons_geos   <- geos::as_geos_geometry(polygons)     
  
    # length(polygons_geos)
    
    # polygons_geos  
    transects_polygons_matrix <- geos::geos_intersects_matrix(transects_geos, polygons_geos) 
    polygons_transects_matrix <- geos::geos_intersects_matrix(polygons_geos, transects_geos) 
    
    # transects_polygons_matrix
    # get the polygons that have atleast 1 intersection with the 'lines'
    # transects_polygons_matrix[[557]]
    # fema_tmp <- fema[transects_polygons_matrix[[557]], ]
    # trans_tmp <- transects[557, ]
    
   # fema_dissolve <-  rmapshaper::ms_dissolve(fema_tmp)
   # fema_simple <- rmapshaper::ms_simplify(fema_tmp, keep = 1) %>% 
     # rmapshaper::ms_dissolve(field = "state")
  
    # mapview::mapview(trans_tmp, color = "green") + 
    #   mapview::mapview(fema_tmp[1, ], col.region = "red") + 
    #   mapview::mapview(fema_tmp[2, ],  col.region = "dodgerblue") + 
    #   mapview::mapview(fema_simple,  col.region = "yellow")
  
    # polygons_transects_matrix
    
    # subset the transects and polygons to only those with intersections
    intersect_transects  <- transects[lengths(transects_polygons_matrix) != 0, ]
    intersect_polygons   <- polygons_geos[lengths(polygons_transects_matrix) != 0]
    
    # Convert our intersecting polygons to LINESTRINGS b/c we DON'T NEED polygons to calculate extension distances from our transect lines
    # This can be done with just linestrings (not sure if this is actually more performent but I'm pretty sure it is....)
    intersect_lines <- 
      intersect_polygons %>% 
      sf::st_as_sf() %>% 
      sf::st_cast("MULTILINESTRING") %>% 
      geos::as_geos_geometry() %>% 
      geos::geos_simplify_preserve_topology(25)
    
    # length(polygons_geos)
    # nrow(fema)
    # polygons_geos     <- geos::as_geos_geometry(intersect_polygons)
    
    # tmp_trans <- intersect_transects[1:5000, ]
    
    # tmp_centroid <- sf::st_centroid(tmp_trans)
    # sf::st_segmentize()
    # split_trans <- lwgeom::st_split(tmp_trans, tmp_centroid)
    # # tmp_trans
    # max_extension_distance <- 2500
    
    # unlist(tmp_trans$geom)
    # unlist(sf::st_segmentize(tmp_trans, 4)$geom)
    
    # use half of the shortest transect line as the segmentation length for all transects (ensures all transects will have a midpoint...?)
    # TODO: Double check this logic.
    min_segmentation <- min(intersect_transects$cs_lengthm %/% 2)
    
    # which.min(intersect_transects$cs_lengthm %/% 2)
    
    # make each transect line have way more segments so we can take a left and right half of each transect line
    segmented_trans <- sf::st_segmentize(intersect_transects, min_segmentation)
    
    # unlist(segmented_trans$geom)
    unique(lengths(segmented_trans$geom))
    length(lengths(segmented_trans$geom))
    lengths(segmented_trans$geom)  
    
    # left_trans   <- lwgeom::st_linesubstring(segmented_trans, 0, 0.50)
    # right_trans  <- lwgeom::st_linesubstring(segmented_trans, 0.50, 1)
    
    # mapview::mapview(left_trans, col.regions = "dodgerblue") +
      # mapview::mapview(intersect_transects, color = "red") +
      #   mapview::mapview(intersect_transects[42, ], color = "yellow") +
      # mapview::mapview(right_trans, color = "dodgerblue") +
      # mapview::mapview(left_trans, color = "green")
    
    # Seperate the transect lines into LEFT and RIGHT halves
    # We do this so we can check if a side of a transect is ENTIRELY WITHIN a polygon. 
    # If the half is entirely within a polygon, 
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
    left_trans_geos     <- geos::as_geos_geometry(left_trans)
    right_trans_geos    <- geos::as_geos_geometry(right_trans)
    
    left_within_matrix  <- geos::geos_within_matrix(left_trans_geos, intersect_polygons)
    right_within_matrix <- geos::geos_within_matrix(right_trans_geos, intersect_polygons)
    
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
    
    left_distances <- calc_extension_distances(
      geos_geoms             = left_trans_geos,
      ids                    = left_trans$tmp_id,
      lines_to_cut           = intersect_lines,
      lines_to_cut_indices   = left_trans$left_fema_index,
      direction              = "head",
      max_extension_distance = max_extension_distance
    )
    
    right_distances <- calc_extension_distances(
      geos_geoms             = right_trans_geos,
      ids                    = right_trans$tmp_id,
      lines_to_cut           = intersect_lines,
      lines_to_cut_indices   = right_trans$right_fema_index,
      direction              = "tail",
      max_extension_distance = max_extension_distance
    )  
    
    left_trans$left_distance   <- left_distances
    right_trans$right_distance <- right_distances
    
    extensions_by_id <- dplyr::left_join(
                          sf::st_drop_geometry(
                            dplyr::select(left_trans, 
                                          hy_id, cs_id, left_distance)
                            ),
                          sf::st_drop_geometry(
                            dplyr::select(right_trans, 
                                        hy_id, cs_id, 
                                        right_distance)
                            ),
                          by = c("hy_id", "cs_id")
                          )
    transects_with_distances <- 
      transects %>% 
      dplyr::left_join(
                extensions_by_id,
                by = c("hy_id", "cs_id")
            ) %>% 
      hydrofabric3D::add_tmp_id()
    
    fline_id_array <- flines$id
    
    # Convert the net object into a geos_geometry
    flines_geos <- geos::as_geos_geometry(flines)
    
    transect_hy_id_array <- transects_with_distances$hy_id
    transect_cs_id_array <- transects_with_distances$cs_id
    
    transect_geoms  <- geos::as_geos_geometry(transects_with_distances$geom)
   
    left_distances  <- transects_with_distances$left_distance
    right_distances <- transects_with_distances$right_distance
    
    # preallocate vector that stores the extension. distances
    new_transects <- vctrs::vec_c(rep(geos::geos_empty(), length(transect_ids)))
    # new_transects <- geos::geos_empty()
    # measures  <- vctrs::vec_c()
    
    for (i in seq_along(transect_ids)) {
      message("i: ", i)
      i = 1
      
      current_trans <- transect_geoms[i]
      
      current_hy_id <- transect_hy_id_array[i]
      current_cs_id <- transect_cs_id_array[i]
      
      current_fline <- flines_geos[fline_id_array == current_hy_id]
    
      # current_fline <- 
      
       
      left_distance_to_extend <- left_distances[i]
      left_distance_to_extend
      right_distance_to_extend <- right_distances[i]
      right_distance_to_extend
      
      no_extension_required <- (left_distance_to_extend == 0 && right_distance_to_extend == 0)
      message("Transect tmp_id: ", curr_tmp_id, " - (", i, ")")
      
      if(no_extension_required) {
        message("Skipping -left/right extension are both 0")
        next
      }
      
      trans
      left_trans
      curr_tmp_id
      
      curr_fema_index <- 
        left_trans %>% 
        dplyr::filter(hy_id == current_hy_id, cs_id == current_cs_id) %>% 
        .$left_fema_index %>% .[[1]]
      
      left_extended_trans  <- hydrofabric3D::geos_extend_line(current_trans, left_distance_to_extend, "head")
      right_extended_trans <- hydrofabric3D::geos_extend_line(current_trans, right_distance_to_extend, "tail")
      
      # neighbor_transects <- 
      left_intersects_fline <- geos::geos_intersection(
        left_extended_trans,
        current_fline
      )
      
      left_intersects_fline_once = geos::geos_type(fline_intersects) == "point"
      left_intersects_fline_once
      
      
      fline_intersects
      geos::geos_type(fline_intersects) == "point"
      
      geos::geos_intersects(current_trans, )
      
      # if(
      #   geos::geos_type(fline_intersects)
      #   
      #   
      # )
      mapview::mapview(sf::st_as_sf(trans), color = "green") +
       mapview::mapview(sf::st_as_sf(left_extended_trans), color = "red")     +
        mapview::mapview(sf::st_as_sf(right_extended_trans), color = "red")     +
        mapview::mapview(      sf::st_as_sf(intersect_polygons[[curr_fema_index]]), col.regions = "dodgerblue")   
      
      
      
        
      
    }
    
    
    left_trans$left_extension_distance
    
    length(right_distances)
    
    
    # Calculate the minimum distance a line would need to extend to reach the boundary of the polygon/line that the input geometries are entirely within 
    calc_extension_distances <- function(geos_geoms, ids, lines_to_cut, lines_to_cut_indices, direction = "head", max_extension_distance = 2500) {
      #####   #####   #####   #####   #####
      # geos_geoms   <- left_trans_geos
      # ids          <- left_trans$tmp_id
      # lines_to_cut <- intersect_lines
      # lines_to_cut_indices <- left_trans$left_fema_index
      # direction = "head"
      # max_extension_distance = 2500
      # geos_geoms             = left_trans_geos
      # ids                    = left_trans$tmp_id
      # lines_to_cut           = intersect_lines
      # lines_to_cut_indices   = left_trans$left_fema_index
      # direction              = "head"
      # max_extension_distance = max_extension_distance
      
      
      
      #####   #####   #####   #####   #####
      
      if (!direction %in% c("head", "tail")) {
          stop("Invalid 'direction' value, must be one of 'head' or 'tail'")
      }
      
      # preallocate vector that stores the extension. distances
      extension_dists <- vctrs::vec_c(rep(0, length(ids)))

      # extension_dists <- vector(mode = "numeric", length = nrow(trans_data))
      for (i in seq_along(ids)) {
        # i = 118
        curr_id           <- ids[i]
        is_within_polygon <- any(!is.na(lines_to_cut_indices[[i]]))
        polygon_index     <- lines_to_cut_indices[[i]]
        # any(is_within_polygon)
        message("Transect: '", curr_id, "' - (", i, ")")
        
        if (is_within_polygon) {
          message("- Side of transect intersects with FEMA")
          message("\t > FEMA index: ", polygon_index)
          
          curr_geom  <- geos_geoms[[i]]
          index_vect <- sort(unlist(polygon_index))
          
           distance_to_extend <- hydrofabric3D:::geos_bs_distance(
            distances    = 1:max_extension_distance,
            line         = curr_geom,
            geoms_to_cut = lines_to_cut[index_vect],
            direction    = direction
          )
          
          extension_dists[i] <- distance_to_extend
        }
        
      }
           
      return(extension_dists)
    }
    
    
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
        
        # mapview::mapview(sf::st_as_sf(left_trans_geos[i]), color = "red") + sf::st_as_sf(intersect_lines[index_vect])
        
        # intersect_lines[index_vect]
        
        dist_to_fema <- hydrofabric3D:::geos_bs_distance(
                            distances    = 1:max_extension_distance,
                            line         = trans_geom,
                            geoms_to_cut = intersect_lines[index_vect],
                            direction    = "head"
                          )
        
        left_extension_dists[i] <- dist_to_fema
        
      }
      message()
    }
    
    left_trans$left_extension_dist <- left_extension_dists
  }
  tmp <- left_trans[1:5, ]
  tmp
  tmp_extended <- hydrofabric3D:::extend_by_length(tmp, tmp$left_extension_dist)
  
  tmp_extended
  
  fema_idx <- unique(unlist(dplyr::select(tmp, left_fema_index)$left_fema_index))
  
  mapview::mapview(dplyr::select(tmp, -left_fema_index), color = "red") +  
    mapview::mapview(dplyr::select(tmp_extended, -left_fema_index), color = "green")  +
    sf::st_as_sf(intersect_polygons[fema_idx])
  
  length(left_extension_dists)
  
  left_trans_geos
  vctrs::vec_c(
    vctrs::vec_c(
      left_trans$hy_id
  ),
  vctrs::vec_cast(left_trans$cs_id)
  )
  
 #  
 #  # mapview::mapview(left_trans[477, ]) + intersect_polygons[1888, ]
 # left_trans <- 
 #    left_trans %>%
 #    dplyr::mutate(
 #      fema_index = unlist(sf::st_within(., intersect_polygons))
 #    ) %>% 
 #    dplyr::relocate(fema_index)
 #  
  # left_trans_geos <- geos::as_geos_geometry(left_trans)
  
  sort(na.omit(unlist(unique(left_trans$fema_index))))
  fema_polygons[na.omit(unlist(unique(left_trans$fema_index)))]
  
  fema_polygons %>% na.omit()
  
  # note: sorting the fema polygon indices (not sure if necessary)
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
    sf::st_cast("multilinestring") %>% 
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
    
    hydrofabric3d:::geos_bs_distance(
      distances    = 1:2000,
      line         = left_trans_geos[i],
      geoms_to_cut = left_fema_lines,
      direction    = "head"
    )
  }) %>% 
    unlist()
  
  left_trans$left_head_extension_dist <- left_extension_dists
  
  # find the distances from the right side of transect lines 
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
                      if(length(i) > 0) { c(i) } else { c(na_real_) } }
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
    sf::st_cast("linestring")
   
  right_extension_dists <- lapply(seq_along(right_trans_geos), function(i) {
    
    hydrofabric3d:::geos_bs_distance(
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
    sf::st_cast("linestring")
  
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
    
   hydrofabric3d:::geos_bs_distance(
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
    
    extended <- hydrofabric3d::geos_extend_line(
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
      hydrofabric3d::geos_extend_line(geom, dist, "head") 
    },
  left_trans_geos,
  left_extension_dists
  )
  
  left_extensions <- lapply(seq_along(left_trans_geos), function(i) {
      
    extend_dist <- hydrofabric3d:::geos_bs_distance(
        distances    = 1:2000,
        line         = left_trans_geos[i],
        geoms_to_cut = left_fema_lines,
        direction    = "head"
      )
      
    hydrofabric3d::geos_extend_line(left_trans_geos[i], extend_dist, "head") 
    
    })
  # unlist(left_extensions)
  
  left_extension_dists <- lapply(seq_along(left_trans_geos), function(i) {
    hydrofabric3d:::geos_bs_distance(
      distances    = 1:2000,
      line         = left_trans_geos[i],
      geoms_to_cut = left_fema_lines,
      direction    = "head"
    )
  }) 
  
  distance_to_extend <- 
    hydrofabric3d:::geos_bs_distance(
      distances    = 1:1500,
      line         = left_trans_geos[1],
      geoms_to_cut = left_fema_lines,
      direction    = "head"
    )
  
  extended <- hydrofabric3d::geos_extend_line(left_trans_geos[1], distance_to_extend, "head") %>% 
    sf::st_as_sf()
  mapview::mapview(left_fema, col.regions = "dodgerblue") +
    mapview::mapview(left_ls, color = "yellow") +
    mapview::mapview(left_trans, color = "red") +
    mapview::mapview(extended, color = "green")
    # mapview::mapview(right_trans, color = "green")
  left_trans$geom %>% sf::st_length()
  plot(segmented_trans$geom, col = "red", lwd =5)
  plot(left_trans$geom, col = "green", lwd=3, add = true)
  plot(right_trans$geom, col = "blue", lwd=3, add = true)
  unlist(left_trans$geom)
   unlist(right_trans$geom)
  
  unlist(trans_fema$geom )
  unlist(split_trans$geom )
  split_trans %>% 
    sf::st_collection_extract("linestring")
  split_trans$geom
  
  tmp_trans$geom
  nngeo::st_segments(tmp_trans) %>% 
    .$result %>% 
    plot()
  
  mapview::mapview(tmp_trans) + tmp_centroid
  geos::geos_clip_by_rect()
  
  transects_with_
  
  lengths(transects_polygons_matrix)
  mapview::mapview(transects_with_fema, color = "green") + fema_with_transects
  unique(hydrofabric3d::add_tmp_id(transects)$tmp_id)[1:30]
  transects  %>% 
    hydrofabric3d::add_tmp_id(transects) %>% 
    .$tmp_id %>% 
    unique() %>% .[1:30]
  trans_subset <- 
    transects %>% 
    hydrofabric3d::add_tmp_id() %>% 
    dplyr::filter(tmp_id %in%  unique(hydrofabric3d::add_tmp_id(transects)$tmp_id)[1:30])
  
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
  
  
  
  
  
  
  
  
  
  
  
  
  
 
  
  
       
  
     
  
     
  
     
  
    
  
    