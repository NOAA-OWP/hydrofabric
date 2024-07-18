#' @title Read QML
#' @description Reads in a QML file 
#' @param qml_file character path
#' @return QML contents
#' @export

read_qml = function (qml_file)  {
  paste(readLines(qml_file), collapse = "\n")
}

#' @title Create style row
#' @param gpkg_path GPKG path
#' @param layer_name layer name
#' @param style_name style name
#' @param style_qml style QML 
#' @return data.frame
#' @export

create_style_row = function (gpkg_path,
                             layer_name,
                             style_name,
                             style_qml) {
  geom_col <-
    sf::st_read( gpkg_path,
      query = paste0(
        "SELECT column_name from gpkg_geometry_columns ", "WHERE table_name = '",
        layer_name, "'"),
      quiet = TRUE)[1, 1]
  
  data.frame(
    f_table_catalog = "",
    f_table_schema = "",
    f_table_name = layer_name,
    f_geometry_column = geom_col,
    styleName = style_name,
    styleQML = style_qml,
    styleSLD = "",
    useAsDefault = TRUE,
    description = "Generated for hydrofabric",
    owner = "",
    ui = NA,
    update_time = Sys.time()
  )
}

#' Append style to GPKG
#'
#' @param gpkg_path GPKG path
#' @param qml_dir Directory path
#' @param layer_names layer names to populate
#' @return gpkg path
#' @export

append_style = function (gpkg_path,
                         qml_dir = system.file("qml", package = "hydrofabric"),
                         layer_names) {
  
  f = list.files(qml_dir, full.names = TRUE)
  
  good_layers = gsub(".qml", "", basename(f))
  
  layer_names = layer_names[layer_names %in% good_layers]
  
  files = grep(paste(layer_names, collapse = "|"), f, value = TRUE)
  
  styles <- sapply(files, read_qml)
  style_names <-
    sapply(layer_names, paste0, "__hydrofabric_style")
  style_rows <-
    do.call(
      rbind,
      mapply(
        create_style_row,
        layer_names,
        style_names,
        styles,
        MoreArgs = list(gpkg_path = gpkg_path),
        SIMPLIFY = FALSE
      )
    )
  if ("layer_styles" %in% sf::st_layers(gpkg_path)$name) {
    try(sf::st_delete(gpkg_path, "layer_styles"), silent = TRUE)
  }
  
  if(!is.null(style_rows)){
    sf::st_write(
      style_rows,
      gpkg_path,
      layer = "layer_styles",
      append = FALSE,
      quiet = TRUE
    )
  }
  
  return(gpkg_path)
}