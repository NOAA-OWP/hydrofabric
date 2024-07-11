#' Get the polygons that interesect with any of the linestring geometries
#' This is just a wrapper around geos::geos_intersects_matrix. Takes in sf dataframes, uses geos, then outputs sf dataframes
#' @param polygons polygon sf object. Default is NULL
#' @param lines linestring sf object. Default is NULL.
#'
#' @return sf dataframe of polygons that intersect with the linestrings
polygons_with_line_intersects <- function(polygons = NULL, lines = NULL) {
  
  if (is.null(polygons)) {
    stop("NULL 'polygons' argument, provide an sf dataframe of POLYGON or MULTIPOLYGON geometries")
  }
  
  if (is.null(lines)) {
    stop("NULL 'lines' argument, provide an sf dataframe of LINESTRING or MULTILINESTRING geometries")
  }
  
  # Convert the SF geometries to geos geometries
  polygons_geos   <- geos::as_geos_geometry(polygons)
  lines_geos      <- geos::as_geos_geometry(lines)
  
  # create an index between the polygons and linestrings
  lines_index <-  geos::geos_intersects_matrix(polygons_geos, lines_geos)
  
  # get the polygons that have atleast 1 intersection with the 'lines'
  polygons_with_lines <- polygons[lengths(lines_index) != 0, ]
  
  return(polygons_with_lines)
}

add_predicate_group_id <- function(polys, predicate) {
  # GROUP BY SPATIAL PREDICATES
  # ----------------------------------------- 
  # predicate = sf::st_touches
  # polys <- sf_df
  # ----------------------------------------- 
  
  
  relations <- predicate(polys)
  
  relations <- lapply(seq_along(relations), function(i) { as.character(sort(unique(c(relations[i][[1]], i)))) })
  
  group_ids_map <- fastmap::fastmap()
  ids_to_groups <- fastmap::fastmap()
  
  group_id <- 0
  
  for (i in seq_along(relations)) {
    
    predicate_ids <- relations[i][[1]]
    
    # message("(", i, ") - ", predicate_ids)
    # message("Start Group ID: ", group_id)
    
    id_group_check <- ids_to_groups$has(predicate_ids)
    
    if(any(id_group_check)) {
      
      known_groups  <- ids_to_groups$mget(predicate_ids)
      known_group   <- known_groups[unname(sapply(known_groups , function(kg) {
        !is.null(kg)
      }))][[1]]
      
      # message("IDs part of past group ID > '", known_group, "'")
      
      past_group_ids     <- group_ids_map$get(known_group)[[1]]
      updated_group_ids  <- as.character(
        sort(as.numeric(unique(c(past_group_ids, predicate_ids))))
      )
      
      group_ids_map$set(known_group, list(updated_group_ids))
      
      new_ids <- predicate_ids[!predicate_ids %in% past_group_ids]
      
      # message("Adding ", new_ids, " to seen set...")
      
      # add any newly added IDs to the seen map
      for (seen_id in new_ids) {
        # message(seen_id)
        ids_to_groups$set(as.character(seen_id), as.character(group_id))
      }
      
    } else {
      # get a new group ID number
      group_id <- group_id + 1    
      # message("IDs form NEW group > '", group_id, "'")
      
      # create a new key in the map with the predicate IDs list as the value
      group_ids_map$set(as.character(group_id), list(predicate_ids))
      
      # message("Adding ", predicate_ids, " to seen set...")
      
      # add each predicate ID to the map storing the seen indexes and their respecitve group IDs 
      for (seen_id in predicate_ids) {
        # message(seen_id)
        ids_to_groups$set(as.character(seen_id), as.character(group_id))
      }
    }
    # message("End group ID: ", group_id, "\n") 
  }
  
  group_ids   <- group_ids_map$as_list() 
  
  grouping_df <- lapply(seq_along(group_ids), function(i) {
    # i = 2
    grouping  <- group_ids[i] 
    group_id  <- names(grouping)
    indices   <- grouping[[1]][[1]]
    
    data.frame(
      index      = as.numeric(indices),
      group_id   = rep(group_id, length(indices))   
    )
    
  }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::arrange(i) 
  
  # count up the number of IDs for each group, well use this to determine which group 
  # to put any indices that had MULTIPLE groups they were apart of (use the group with the most other members)
  group_id_counts <- 
    grouping_df %>% 
    dplyr::group_by(group_id) %>% 
    dplyr::count() %>% 
    # dplyr::arrange(-n) %>% 
    dplyr::ungroup()
  
  # select the IDs with the most other members
  grouping_df <- 
    grouping_df %>% 
    dplyr::left_join(
      group_id_counts, 
      by = 'group_id'
    ) %>% 
    dplyr::group_by(index) %>% 
    dplyr::slice_max(n, with_ties = FALSE) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-n) %>% 
    dplyr::arrange(-index) 
  
  polys$group_id <- grouping_df$group_id
  
  return(polys)
  
}

