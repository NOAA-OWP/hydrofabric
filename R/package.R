#' @keywords internal
"_PACKAGE"

#' Hydrofabric Data Model
#' @family catalog

"hf_dm"

NULL

#' @importFrom DBI dbConnect dbDisconnect dbListTables
#' @importFrom RSQLite SQLite
#' @importFrom arrow read_parquet write_parquet open_dataset
#' @importFrom glue  glue

#' @rawNamespace import(climateR, except = c(plot))
#' @rawNamespace import(dplyr, except = c(intersect, union))

#' @import hydrofab
#' @import ngen.hydrofab
#' @import nhdplusTools
#' @import zonal
#' @import terra
#' @import sf

#' @export
arrow::read_parquet

#' @export
arrow::write_parquet

#' @export
arrow::open_dataset

#' @export
glue::glue