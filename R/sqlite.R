#' @title Read gpkg as SQLite
#' @param gpkg path
#' @param lyr layer of GPKG. If NULL (default the database tables are printed)
#' @returns an S4 object that inherits from DBIConnection. This object is used to communicate with the database engine.
#' @export 

as_sqlite = function(gpkg, lyr= NULL){
  
  db = dbConnect(RSQLite::SQLite(), gpkg)
  
  if(is.null(lyr)){
    print(dbListTables(db))
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
