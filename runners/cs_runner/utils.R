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
