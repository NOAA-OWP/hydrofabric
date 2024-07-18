#' @title Read gpkg as SQLite
#' @param gpkg path
#' @param lyr layer of GPKG. If NULL (default the database tables are printed)
#' @param ignore pattern for layers to be ignored description
#' @returns an S4 object that inherits from DBIConnection. This object is used to communicate with the database engine.
#' @export 

as_sqlite = function(gpkg, lyr= NULL, ignore = "gpkg_|rtree_|sqlite_"){
  
  db = dbConnect(RSQLite::SQLite(), gpkg)
  
  if(is.null(lyr)){
    tbls = dbListTables(db)
    tbls = tbls[!grepl(ignore, tbls)]
    print(tbls)
    dbDisconnect(db)
  } else {
    if(lyr %in% dbListTables(db)){
      db = tbl(db, lyr)
    } else {
      stop(lyr, ' not in gpkg.', call. = FALSE)
      dbDisconnect(db)
    }
  }
  db
}
  
#' @title Extract spatial data from an SQLite connection
#' @param tbl  remote temporary table
#' @returns an sf object (or data.frame is non spatial)
#' @export 

read_sf_dataset_sqlite = function(tbl){
  
  srs = dbConnect(RSQLite::SQLite(), tbl$src$con@dbname) %>% 
    tbl("gpkg_spatial_ref_sys") %>% 
    collect() %>% 
    slice_tail(n = 1)
  
  d = collect(tbl)
  
  if(any(c("geom", "geometry") %in% names(d))){
    st_as_sf(d, crs = srs$definition)
  } else {
    warning("no simple features geometry column present")
    d
  }
}