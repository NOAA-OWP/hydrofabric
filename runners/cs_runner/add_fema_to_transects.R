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
  
 
  fema_vpu_dir <- paste0(FEMA_VPU_SUBFOLDERS[grepl(paste0("VPU_", vpu), basename(FEMA_VPU_SUBFOLDERS))], "/merged")
  vpu_fema_files <- list.files(fema_vpu_dir, full.names = TRUE)
  # vpu_fema_file1 <- vpu_fema_files[grepl("_union.gpkg", vpu_fema_files)]
  
  vpu_fema_file <- vpu_fema_files[grepl(paste0(vpu, ".gpkg"), vpu_fema_files)]
  # vpu_fema_file <- vpu_fema_files[grepl(paste0(vpu, "_union.gpkg"), vpu_fema_files)]
  vpu_fema_file
  
  # fema polygons and transect lines
  fema <- sf::read_sf(vpu_fema_file)
  # mapview::npts(fema)
  
  transects <- sf::read_sf(transect_path)
  
  # read in nextgen flowlines data
  flines <- sf::read_sf(nextgen_path, layer = "flowpaths")

  # fema$fema_id %>% unique() %>% length()
  
  # # union then explode FEMA polygons 
  # fema <- 
  #   fema %>%
  #   sf::st_union()
  # 
  # fema <- rmapshaper::ms_explode(fema)
  # # fema %>% mapview::npts()
  # 
  # # reassign IDs and change geometry column name
  # fema <- 
  #   fema %>% 
  #   sf::st_as_sf() %>% 
  #   dplyr::mutate(fema_id = 1:dplyr::n()) %>% 
  #   dplyr::select(fema_id, geom = x)
  
  # sf::st_union() %>% 
  #   rmapshaper::ms_explode() %>% 
  #   sf::st_as_sf() %>% 
  #   dplyr::mutate(fema_id = 1:dplyr::n()) %>% 
  #   dplyr::select(fema_id, geom = x)
  # Given 2 geos_geometry point geometries, create a line between the 2 points
  # start: geos_geoemtry, point
  # end: geos_geoemtry, point
  # Returns geos_geometry linestring
  make_line_from_start_and_end_pts <- function(start, end, line_crs) {
        Y_start <- geos::geos_y(start)
        X_start <- geos::geos_x(start)
        Y_end   <- geos::geos_y(end)
        X_end   <- geos::geos_x(end)
        
        # make the new transect line from the start and points 
        geos_ls <- geos::geos_make_linestring(x = c(X_start, X_end),
                                                 y = c(Y_start, Y_end), 
                                                 crs = line_crs)
                                              
        return(geos_ls)                                              
        
        
      }
  
  # Give a set of transecct linestrings and poylgons and get the minimum distance to extend each transect line (from both directions, to try and reach the edge of a "polygons")
  # internal function for extending transect lines out to FEMA 100 year flood plain polygons
  # transect_lines, set of Sf linestrigns to extend (only if the transect lines are ENTIRELLY within a polygons)
  # polygons, set of sf polygons that transect lines should be exteneded 
  # max_extension_distance numeric, maximum distance (meters) to extend a transect line in either direction to try and intersect one of the "polygons"
  get_transect_extension_distances_to_polygons <- function(transect_lines, 
                                                           polygons, 
                                                           flines,
                                                           max_extension_distance) {
    
    ###    ###    ###    ###    ###    ###    ###
    ###    ###    ###    ###    ###    ###    ###
    transect_lines <- transects
    polygons <- fema
    # # flines <- flines
    # max_extension_distance <- 3000
    max_extension_distance = 3500
    # ###    ###    ###    ###    ###    ###    ###
    ###    ###    ###    ###    ###    ###    ###
    
    # keep 10% of the original points for speed
    polygons <- rmapshaper::ms_simplify(polygons, keep_shapes = F, keep = 0.10)
    
    mapview::npts(fema)
    mapview::npts(polygons)
    # transects <- sf::read_sf(transect_path)
    
     # polygons
    transects_geos  <- geos::as_geos_geometry(transects)
    polygons_geos   <- geos::as_geos_geometry(polygons)     
    
    # geos::geos_union( polygons_geos, polygons_geos) %>% length()
    # polygons_geos %>% length()
    # polygons
    # polygons_geos  
    transects_polygons_matrix <- geos::geos_intersects_matrix(transects_geos, polygons_geos) 
    polygons_transects_matrix <- geos::geos_intersects_matrix(polygons_geos, transects_geos) 

    # subset the transects and polygons to only those with intersections
    intersect_transects  <- transects[lengths(transects_polygons_matrix) != 0, ]
    intersect_polygons   <- polygons_geos[lengths(polygons_transects_matrix) != 0]
    
    # Convert our intersecting polygons to LINESTRINGS b/c we DON'T NEED polygons to calculate extension distances from our transect lines
    # This can be done with just linestrings (not sure if this is actually more performent but I'm pretty sure it is....)
    intersect_lines <- 
      intersect_polygons %>% 
      # geos::geos_make_valid() %>% 
      sf::st_as_sf() %>% 
      # sf::st_union() %>% 
      # rmapshaper::ms_explode() %>% 
      # sf::st_as_sf() %>% 
      # dplyr::mutate(fema_id = 1:dplyr::n()) %>% 
      # dplyr::select(fema_id, geom = x) %>% 
      sf::st_cast("MULTILINESTRING") %>% 
      geos::as_geos_geometry() %>% 
      geos::geos_simplify_preserve_topology(20)
    
    # mapview::npts(sf::st_as_sf(intersect_lines))

    
 #    intersect_polygons %>% 
 #      geos::geos_make_valid() %>% 
 #      geos::geos_is_valid() %>% all()
 # is.null(intersect_lines$geometry )
    # use half of the shortest transect line as the segmentation length for all transects (ensures all transects will have a midpoint...?)
    # TODO: Double check this logic.
    min_segmentation <- min(intersect_transects$cs_lengthm %/% 2)
    
    # which.min(intersect_transects$cs_lengthm %/% 2)
    
    # make each transect line have way more segments so we can take a left and right half of each transect line
    segmented_trans <- sf::st_segmentize(intersect_transects, min_segmentation)
    
    # unlist(segmented_trans$geom)
    # unique(lengths(segmented_trans$geom))
    # length(lengths(segmented_trans$geom))
    # lengths(segmented_trans$geom)  
    
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
    
    # max_extension_distance = 3000
    # which(transects_with_distances$hy_id == "wb-1003839")
    
    left_trans[which(left_trans$hy_id == "wb-1003839"), ]$cs_lengthm
    right_trans[which(left_trans$hy_id == "wb-1003839"), ]$cs_lengthm
    
    # which(right_trans$hy_id == "wb-1003839")
    
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
    
    left_trans$left_distance    <- left_distances
    right_trans$right_distance  <- right_distances
    
    extensions_by_id %>% 
      dplyr::filter(hy_id == "wb-1003839")
    
    # distance to extend LEFT and/or RIGHT for each hy_id/cs_id
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
    
     ######### ######## ######## ######## ####### 
    # TODO: This is temporary !!!!!!
    fema_index_df <- dplyr::left_join(
      sf::st_drop_geometry(
        dplyr::select(left_trans, 
                      hy_id, cs_id, left_distance, left_fema_index)
      ),
      sf::st_drop_geometry(
        dplyr::select(right_trans, 
                      hy_id, cs_id, 
                      right_distance, right_fema_index)
      ),
      by = c("hy_id", "cs_id")
    )
    ######## ######## ######## ######## ########
    
    # foi <- sf::st_as_sf(intersect_polygons[fema_uids]) %>% dplyr::mutate(
    #   fema_id = fema_uids
    # )
    # 
    # polygons_to_merge <- sf::st_as_sf(intersect_polygons[fema_uids]) %>% dplyr::mutate(
    #   fema_id = fema_uids
    # ) %>% 
    #   dplyr::filter(fema_id %in% c(1563, 1566, 1567, 590))
    # merged_polygon <- 
    #   polygons_to_merge %>% 
    #   sf::st_union()
    # merged_polygon %>% 
    #   rmapshaper::ms_explode()
    # mapview::mapview(foi, col.regions = "dodgerblue") +
    #   mapview::mapview(polygons_to_merge, col.regions = "yellow") +
    #   mapview::mapview(merged_polygon, col.regions = "green") +
    #   mapview::mapview(toi, color = "red") +
    #   mapview::mapview(og_trans, color = "green")
    # polygons %>% 
      # dplyr::filter(fema_id )
    ######## ######## ######## ######## ########
    
    # TODO: Add left/right extension distancces to transect data
    # TODO: this can ultimately just be the "transects" variable, dont need to make new "transects_with_distances" variable
    transects_with_distances <- 
      transects %>% 
      dplyr::left_join(
                extensions_by_id,
                by = c("hy_id", "cs_id")
            ) %>% 
      dplyr::mutate(
        left_distance = dplyr::case_when(
          is.na(left_distance) ~ 0,
          TRUE                 ~ left_distance
        ),
        right_distance = dplyr::case_when(
          is.na(right_distance) ~ 0,
          TRUE                  ~ right_distance
        )
        # left_distance = dplyr::case_when(
        #   left_distance == 0 ~ NA,
        #   TRUE               ~ left_distance
        # ),
        # right_distance = dplyr::case_when(
        #   right_distance == 0 ~ NA,
        #   TRUE                ~ right_distance
        # )
      ) %>% 
      hydrofabric3D::add_tmp_id()
    
    transects_with_distances %>% 
      dplyr::filter(hy_id == "wb-1003839")
    # transects
    # flines %>% 
    #   dplyr::arrange(-tot_drainage_areasqkm) 
    # flines %>% 
    #   dplyr::arrange(-hydroseq) %>% 
    #   dplyr::filter(hydroseq == min(hydroseq) | hydroseq == max(hydroseq)) %>%
    #   # dplyr::filter(hydroseq == max(hydroseq)) %>% 
    #   mapview::mapview()
    
    # extensions_by_id
    transects_with_distances %>% 
      dplyr::filter(!tmp_id %in% hydrofabric3D::add_tmp_id(extensions_by_id)$tmp_id)
    
    ######## ######## ######## ######## ######## ######## ######## ######## 
    left_extended_flag[1:20]
    right_extended_flag[1:20]
    both_extended_flag[1:20]
    
    range_vect <- 1:500
    
    fema_uids <- unique(c(unlist(fema_index_df[range_vect, ]$left_fema_index),     unlist(fema_index_df[range_vect, ]$right_fema_index)))
    fema_uids <- fema_uids[!is.na(fema_uids)]
fema_uids
    foi <- sf::st_as_sf(intersect_polygons[fema_uids]) %>% dplyr::mutate(
      fema_id = fema_uids
    )
    # ooi <- sf::st_as_sf()
    # toi <- sf::st_as_sf(new_transects[1:20])
     toi <- sf::st_as_sf(transect_geoms[range_vect])
    toi
    og_trans <- transects_with_distances[range_vect, ]
    mapview::mapview(foi, col.regions = "dodgerblue") +
      mapview::mapview(toi, color = "red") +
      mapview::mapview(og_trans, color = "green")
    ######## ######## ### ##### ######## ######## ######## ######## ######## 
    
    fline_id_array       <- flines$id
    
    # Convert the net object into a geos_geometry
    flines_geos          <- geos::as_geos_geometry(flines)
    
    transect_hy_id_array <- transects_with_distances$hy_id
    transect_cs_id_array <- transects_with_distances$cs_id
    
    transect_geoms       <- geos::as_geos_geometry(transects_with_distances$geom)
   
    left_distances       <- transects_with_distances$left_distance
    right_distances      <- transects_with_distances$right_distance
    
    # preallocate vector that stores the extension. distances
    new_transects <- vctrs::vec_c(rep(geos::geos_empty(), length(transect_hy_id_array)))
    
    left_extended_flag   <- rep(FALSE, length(transect_hy_id_array))   
    right_extended_flag  <- rep(FALSE, length(transect_hy_id_array))
    both_extended_flag   <- rep(FALSE, length(transect_hy_id_array))
    
    # new_transects <- geos::geos_empty()
    # # measures  <- vctrs::vec_c()
    # transects_with_distances[1:20, ]
    # transects[1:20, ]
    
    # number of geometries that will be iterated over, keeping this variable to reference in message block  
    total <- length(transect_hy_id_array)

    # output a message every ~10% intervals
    message_interval <- total %/% 5
    number_of_skips = 0
    
    for (i in seq_along(transect_hy_id_array)) {
      # message("i: ", i)
      # if(i > 2000) {
      #   break
      # }
      # i = 1
      # Check if the iteration is a multiple of 100
      if (i %% message_interval == 0) {
        
        # get the percent complete
        percent_done <- round(i/total, 2) * 100
        
        # Print the message every "message_interval"
        # if(verbose) { 
          message(i, " > ", percent_done, "% ") 
          message("Number of skips: ", number_of_skips)
          # }
        # message("Iteration ", i, " / ", length(extended_trans), 
        #         " - (", percent_done, "%) ")
        
      }
      # which(transects_with_distances$hy_id == "wb-1003839")
      # i = 9587
      # get the current transect, hy_id, cs_id, flowline, and extension distances
      current_trans <- transect_geoms[i]
      
      current_hy_id <- transect_hy_id_array[i]
      current_cs_id <- transect_cs_id_array[i]
      
      current_fline <- flines_geos[fline_id_array == current_hy_id]
       
      left_distance_to_extend <- left_distances[i]
      right_distance_to_extend <- right_distances[i]
      
      no_extension_required <- (left_distance_to_extend == 0 && right_distance_to_extend == 0)
      # no_extension_required <- is.na(left_distance_to_extend) && is.na(right_distance_to_extend)
      # message("Transect tmp_id: ", curr_tmp_id, " - (", i, ")")
      
      if(no_extension_required) {
        # message("Skipping -left/right extension are both 0")
        number_of_skips = number_of_skips + 1
        
        next
      }
      
      # transects_with_distances %>% 
      #   dplyr::filter(hy_id == "wb-1003839")
      # transects
#      curr_fema_index <- 
#        left_trans %>% 
#        dplyr::filter(hy_id == current_hy_id, cs_id == current_cs_id) %>% 
#        .$left_fema_index %>% .[[1]]
      
      # message("Extending transect line left and right")
      # extend the lines
      left_extended_trans  <- hydrofabric3D::geos_extend_line(current_trans, left_distance_to_extend, "head")
      right_extended_trans <- hydrofabric3D::geos_extend_line(current_trans, right_distance_to_extend, "tail")
      
      # neighbor_transects <- 
      # message("Checking left and right intersections with flowline...")
      
      # Check that the extended transect lines only intersect the current flowline once
      left_intersects_fline <- geos::geos_intersection(
        left_extended_trans,
        # current_fline
        flines_geos
      )
      
      right_intersects_fline <- geos::geos_intersection(
        right_extended_trans,
        # current_fline
        flines_geos
      )
      
      # sum(geos::geos_type(left_intersects_fline) == "point")
      # sum(geos::geos_type(right_intersects_fline) == "point")
      
      
      # which(geos::geos_type(left_intersects_fline) == "point")
      
      # transects %>% 
      #   dplyr::filter(hy_id == current_hy_id) %>% 
      #   sf::st_length()
      # current_hy_id
      # 
      # sf::st_as_sf(current_trans) %>% 
      #   sf::st_length()
      # tmp <- sf::read_sf(transect_path)
      # 
      # tmp %>% 
      #   dplyr::filter(hy_id == current_hy_id) %>% 
      #   sf::st_length()
      
      # mapview::mapview(sf::st_as_sf(flines_geos[which(geos::geos_type(left_intersects_fline) == "point")])) +
      #   mapview::mapview(sf::st_as_sf(current_trans), color = "red") +
      # mapview::mapview(sf::st_as_sf(left_extended_trans), color = "green") + 
      # mapview::mapview(sf::st_as_sf(right_extended_trans), color = "green")
      # Define conditions to decide which version of the transect to use
      
      # 1. Use transect with extension in BOTH directions
      # 2. Use transect with LEFT extension only
      # 3. Use transect with RIGHT extension only
      
      # left_intersects_fline_once  <- geos::geos_type(left_intersects_fline) == "point"
      # right_intersects_fline_once <- geos::geos_type(right_intersects_fline) == "point"
      left_intersects_fline_once  <- sum(geos::geos_type(left_intersects_fline) == "point") == 1 && sum(geos::geos_type(left_intersects_fline) == "multipoint") == 0 
      right_intersects_fline_once <- sum(geos::geos_type(right_intersects_fline) == "point") == 1 && sum(geos::geos_type(right_intersects_fline) == "multipoint") == 0 
      
      # sum(geos::geos_type(left_intersects_fline) == "point") == 1
      # sum(geos::geos_type(right_intersects_fline) == "point") == 1
      # sum(geos::geos_type(left_intersects_fline) == "multipoint") == 0 
      
      
      # # TODO: Consider doing the opppsite of these conditions (i.e. "left_intersects_other_transects" = TRUE) 
      # left_does_not_intersect_other_transects  <- !any(geos::geos_intersects(left_extended_trans, transect_geoms[-i]))
      # right_does_not_intersect_other_transects <- !any(geos::geos_intersects(right_extended_trans, transect_geoms[-i]))  
      # 
      # use_left_extension  <- left_intersects_fline_once && left_does_not_intersect_other_transects
      # use_right_extension <- right_intersects_fline_once && right_does_not_intersect_other_transects
      # use_both_extensions <- use_left_extension && use_right_extension
      
      
      # TODO: This is the opposite phrasing of these conditions, i think this is clearer to read
      left_intersects_other_transects  <- any(geos::geos_intersects(left_extended_trans, transect_geoms[-i]))
      right_intersects_other_transects <- any(geos::geos_intersects(right_extended_trans, transect_geoms[-i]))  
      
      # # make sure the extended transects don't hit any of the newly extended transects
      # # NOTE: I think this could be just done with a single transect list that starts with the original transects and if an update happens then we replace that transect
      # left_intersects_new_transects  <- any(geos::geos_intersects(left_extended_trans, new_transects))
      # right_intersects_new_transects <- any(geos::geos_intersects(right_extended_trans, new_transects))
      
      # make TRUE/FALSE flags stating which transect should we use
      # - BOTH extensions
      # - LEFT ONLY extensions
      # - RIGHT only extensions
      use_left_extension  <- left_intersects_fline_once && !left_intersects_other_transects
      use_right_extension <- right_intersects_fline_once && !right_intersects_other_transects
      use_both_extensions <- use_left_extension && use_right_extension
      
      # merged_trans <- geos::geos_union(left_extended_trans, right_extended_trans)
      # sf::st_union(sf::st_cast(sf::st_as_sf(merged_trans), "LINESTRING"))
      # mapview::mapview(sf::st_as_sf(merged_trans), color = "green") +
      #   mapview::mapview(sf::st_as_sf(left_start), col.region = "red") +
      #   mapview::mapview(sf::st_as_sf(left_end), col.region = "red") + 
      #   mapview::mapview(sf::st_as_sf(right_start), col.region = "dodgerblue") +
      #   mapview::mapview(sf::st_as_sf(right_end), col.region = "dodgerblue")
      
      # if(use_both_extensions) {
      
      # Get the start and end of both extended tranects
      left_start  <- geos::geos_point_start(left_extended_trans)
      left_end    <- geos::geos_point_end(left_extended_trans)
      right_start <- geos::geos_point_start(right_extended_trans)
      right_end   <- geos::geos_point_end(right_extended_trans)
      
      # }
      # Extend in BOTH directions
      if(use_both_extensions) {
        # message("Extend direction: BOTH")
        start  <- left_start
        end    <- right_end
        
      # extend ONLY the left side
      } else if(use_left_extension && !use_right_extension) {
         # message("Extend direction: LEFT")       
        start  <- left_start
        end    <- left_end
        
      # Extend ONLY the right side
      } else if(!use_left_extension && use_right_extension) {
         # message("Extend direction: RIGHT")       
        start  <- right_start
        end    <- right_end
        
      # DO NOT extend either direction
      } else {
        # message("No extension")   
        # TODO: Really dont need to do anything 
        # TODO: in this scenario because we just use the original transect line
        start  <- left_end
        end    <- right_start
      }
      
      # trans_needs_extension <- use_left_extension || use_right_extension
      
      # if(trans_needs_extension) {
        
        line_crs      <- wk::wk_crs(current_trans)
        updated_trans <- make_line_from_start_and_end_pts(start, end, line_crs)
        
      # }
        
        if(use_left_extension) {
          left_extended_flag[i]  <- TRUE
        }
        
        if(use_right_extension) {
          right_extended_flag[i] <- TRUE
        }
        
        if(use_both_extensions) {
          both_extended_flag[i] <- TRUE
        }
        
        # new_transects[i] <- updated_trans
        transect_geoms[i] <- updated_trans
        
        # start %>% class()
      
     
    }
    
    transects2 <- transects 
      # dplyr::mutate(
      #   new_cs_lengthm = as.numeric(sf::st_length(geom))
      # ) %>% 
      # dplyr::relocate(hy_id, cs_id, cs_lengthm, new_cs_lengthm)
    
    
    # Update the "transects_to_extend" with new geos geometries ("geos_list")
    sf::st_geometry(transects2) <- sf::st_geometry(sf::st_as_sf(transect_geoms))
    
    transects2 <- 
      transects2 %>% 
      dplyr::mutate(
        new_cs_lengthm = as.numeric(sf::st_length(geom))
      ) %>% 
      dplyr::relocate(hy_id, cs_id, cs_lengthm, new_cs_lengthm)
    
    # transects2 %>% 
    #   dplyr::filter(
    #     new_cs_lengthm > cs_lengthm
    #   )
    # 
    
    transects2$left_is_extended  <- left_extended_flag
    transects2$right_is_extended <- right_extended_flag
    
    transects2 %>% 
      dplyr::filter(left_is_extended, right_is_extended)
    
    any_extended <- 
      transects2 %>% 
      dplyr::filter(left_is_extended | right_is_extended)
    
    any_flines <- 
      flines %>% 
      dplyr::filter(id %in% any_extended$hy_id)
    
    # left_only_extended <- 
    #   transects2 %>% 
    #   dplyr::filter(left_is_extended, !right_is_extended)
    # 
    # left_only_flines <- 
    #   flines %>% 
    #   dplyr::filter(id %in% left_only_extended$hy_id)
    # 
    # right_only_extended <- 
    #   transects2 %>% 
    #   dplyr::filter(!left_is_extended, right_is_extended)
    # 
    # right_only_flines <- 
    #   flines %>% 
    #   dplyr::filter(id %in% right_only_extended$hy_id)
  
    # left_fema_polygons <- 
      # left_trans %>% 
    # ----------------------------------------------------------------
    # ------- Subset data for mapping ----------- 
    # ----------------------------------------------------------------
    extended_for_map <- 
      any_extended %>% 
      dplyr::slice(1:5000)
    
    og_transects_for_map <- 
      transects %>% 
      hydrofabric3D::add_tmp_id() %>% 
      dplyr::filter(tmp_id %in% hydrofabric3D::add_tmp_id(extended_for_map)$tmp_id)
    
    fema_indexes_in_aoi <-
      dplyr::bind_rows(
        sf::st_drop_geometry(
          dplyr::rename(left_trans, fema_index = left_fema_index)
        ),
        sf::st_drop_geometry(     
          dplyr::rename(right_trans, 
                        fema_index = right_fema_index)
        )
      ) %>% 
      hydrofabric3D::add_tmp_id() %>% 
      dplyr::filter(tmp_id %in% hydrofabric3D::add_tmp_id(extended_for_map)$tmp_id) %>% 
      # dplyr::filter(
      #   # tmp_id %in% hydrofabric3D::add_tmp_id(left_only_extended)$tmp_id |
      #   # tmp_id %in% hydrofabric3D::add_tmp_id(right_only_extended)$tmp_id
      #  
      #   tmp_id %in%  unique(hydrofabric3D::add_tmp_id(dplyr::filter(transects2, left_is_extended, right_is_extended))$tmp_id)
      #   
      #   ) %>% 
      # dplyr::filter(left_is_within_fema | right_is_within_fema) %>% 
      # dplyr::slice(1:200) %>%
      .$fema_index %>% 
      unlist() %>% 
      na.omit() %>% 
      unique() %>% 
      sort()
      # length()
    
    sf::st_as_sf(intersect_polygons[fema_indexes_in_aoi])

    # hydrofabric3D::add_tmp_id(left_only_extended)$tmp_id
    # transects_with_distances
    # %>% 
    mapview::mapview( sf::st_as_sf(intersect_polygons[fema_indexes_in_aoi]),  col.regions = "lightblue") + 
      mapview::mapview(any_flines, color = "dodgerblue") + 
      mapview::mapview(og_transects_for_map, color = "green") + 
      mapview::mapview(extended_for_map, color = "red") 

      # mapview::mapview(left_only_flines, color = "dodgerblue") + 
      #     mapview::mapview(right_only_flines, color = "dodgerblue") + 
    # mapview::mapview(left_only_extended, color = "red") + 
    #      mapview::mapview(right_only_extended, color = "green")
    # transects$cs_lengthm  <- length_list
    # ----------------------------------------------------------------
    # ----------------------------------------------------------------
    # ----------------------------------------------------------------  
  
      # make the new transect line from the start and points 
      final_line <- geos::geos_make_linestring(x = c(X_start, X_end),
                                               y = c(Y_start, Y_end), 
                                               crs = wk::wk_crs(current_trans)
                                               )
      geos::geos_make_collection(start, type_id = "LINESTRING")
      
      
      left_start <- geos::geos_point_start(left_extended_trans)
      right_end  <- geos::geos_point_end(right_extended_trans)
      
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
        mapview::mapview(      sf::st_as_sf(intersect_polygons[[curr_fema_index]]), col.regions = "dodgerblue")  +  
      mapview::mapview(      sf::st_as_sf(final_line), color = "yellow")   
      
      
      
      
        
      
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
      # geos_geoms             = left_trans_geos
      # ids                    = left_trans$tmp_id
      # lines_to_cut           = intersect_lines
      # lines_to_cut_indices   = left_trans$left_fema_index
      # direction              = "head"
      # max_extension_distance = max_extension_distance
      # 
      
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
  
  
  
  
  
  
  
  
  
  
  
  
  
 
  
  
       
  
     
  
     
  
     
  
    
  
    