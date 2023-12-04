is.url <- function(x) {
  grepl("www.|http:|https:", x)
}

read_qml = function (qml_file)  {
  paste(readLines(qml_file), collapse = "\n")
}

create_style_row = function (gpkg_path,
                             layer_name,
                             style_name,
                             style_qml) {
  geom_col <-
    sf::st_read(
      gpkg_path,
      query = paste0(
        "SELECT column_name from gpkg_geometry_columns ",
        "WHERE table_name = '",
        layer_name,
        "'"
      ),
      quiet = TRUE
    )[1,
      1]
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
  sf::st_write(
    style_rows,
    gpkg_path,
    layer = "layer_styles",
    append = FALSE,
    quiet = TRUE
  )
  return(gpkg_path)
}

#' Access Hydrofabric Network
#' @param VPU Vector Processing Unit
#' @param base_s3 the base hydrofabric directory to access in Lynker's s3
#' @param cache_dir should data be cached to a local directory? Will speed up multiple subsets in the same region
#' @param cache_overwrite description. Should a cached file be overwritten
#' @return file path
#' @export

get_fabric = function(VPU,
                      base_s3 = 's3://lynker-spatial/v20/',
                      cache_dir = NULL,
                      cache_overwrite = FALSE) {
  Key <- NULL
  
  xx = aws.s3::get_bucket_df(
    bucket = dirname(base_s3),
    prefix = glue::glue("{basename(base_s3)}/gpkg"),
    region = 'us-west-2'
  ) |>
    dplyr::filter(grepl(basename(base_s3), Key) &
                    grepl(paste0(VPU, ".gpkg$"), Key)) |>
    dplyr::filter(!grepl("[.]_", Key))
  
  if (!is.null(cache_dir)) {
    dir.create(cache_dir,
               recursive = TRUE,
               showWarnings = FALSE)
    gpkg = glue::glue("{cache_dir}/{basename(xx$Key)}")
    if (cache_overwrite) {
      unlink(gpkg)
    }
    temp = FALSE
  } else {
    gpkg = tempfile(fileext  = ".gpkg")
    temp = TRUE
  }
  
  if (!file.exists(gpkg)) {
    aws.s3::save_object(
      bucket = xx$Bucket,
      object = xx$Key,
      file = gpkg,
      region = 'us-west-2'
    )
  }
  
  return(gpkg)
  
}

#' Subset Hydrofabric Network
#' @param id hydrofabric id (relevant only to nextgen fabrics)
#' @param comid NHDPlusV2 COMID
#' @param hl_uri hydrolocation URI (relevant only to nextgen fabrics)
#' @param nldi_feature list with names 'featureSource' and 'featureID' where 'featureSource' is derived from the "source" column of the response of dataRetrieval::get_nldi_sources() and the 'featureID' is a known identifier from the specified 'featureSource'.
#' @param xy Location given as vector of XY in CRS 4326 (long, lat)
#' @param bbox a numeric vector of length four, with xmin, ymin, xmax and ymax values
#' @param base_s3 the base hydrofabric directory to access in Lynker's s3
#' @param base_dir the base hydrofabric directory
#' @param lyrs layers to extract. Default is all possible in the hydrofabric GPKG data model
#' @param areasqkm the maximum upstream area to subset
#' @param pathlengthkm the maximum upstream path length to subset
#' @param ms_pathlengthkm the maximum upstream mainstem path length to subset
#' @param outfile file path to write to. Must have ".gpkg" extension
#' @param cache_dir should data be cached to a local directory? Will speed up multiple subsets in the same region
#' @param cache_overwrite description. Should a cached file be overwritten
#' @return file path (outfile) or list of features
#' @export

subset_network = function(id = NULL,
                          comid = NULL,
                          hl_uri = NULL,
                          nldi_feature = NULL,
                          xy = NULL,
                          bbox = NULL,
                          base_s3 = 's3://lynker-spatial/v20/',
                          base_dir = NULL,
                          lyrs  = c(
                            "divides",
                            "nexus",
                            "flowpaths",
                            "network",
                            "hydrolocations",
                            "reference_flowline",
                            "reference_catchment",
                            "refactored_flowpaths",
                            "refactored_divides"
                          ),
                          areasqkm = NULL,
                          pathlengthkm = NULL,
                          ms_pathlengthkm = NULL,
                          outfile = NULL,
                          cache_dir = NULL,
                          qml_dir = system.file("qml", package = "hydrofabric"),
                          cache_overwrite = FALSE) {
  Key <-
    hf_hydroseq <-
    hf_id  <- hydroseq <- member_COMID <- toid  <- vpu <- NULL
  
  lookup <-c(id = "ID", id = "COMID", toid = "toID", toid = "toCOMID")
  
  if (!is.null(bbox)) {
    stopifnot(length(bbox) == 4)
    x = sf::st_as_sfc(sf::st_bbox(bbox, crs = 4326))
    vpuid = sf::st_filter(sf::st_transform(vpu_boundaries, 4326), x)$VPUID
  } else {
    if (!is.null(base_dir)) {
      base = data.frame(file = list.files(base_dir,
                                          full.names = TRUE,
                                          recursive = TRUE)) |>
        dplyr::mutate(base = basename(file))
      
      if ("conus_net.parquet" %in% base$base) {
        net = filter(base, base == "conus_net.parquet") %>%
          dplyr::pull(file) %>%
          arrow::open_dataset() |>
          dplyr::select(id, toid, hf_id, hl_uri, hf_hydroseq, hydroseq, vpu) |>
          dplyr::collect() |>
          dplyr::distinct()
      } else {
        stop("conus_net.parquet not found in ", base_dir)
      }
    } else {
      net = tryCatch({
        arrow::open_dataset(glue::glue(base_s3, "conus_net.parquet")) |>
          dplyr::select(id, toid, hf_id, hl_uri, hf_hydroseq, hydroseq, vpu) |>
          dplyr::collect() |>
          dplyr::distinct()
      } , error = function(e) {
        NULL
      })
    }
    
    
    if (!is.null(id) & !is.null(net)) {
      comid = dplyr::filter(net, id == !!id | toid == !!id) |>
        dplyr::slice_max(hf_hydroseq) |>
        dplyr::pull(hf_id)
    }
    
    if (!is.null(nldi_feature)) {
      if (length(nldi_feature) == 1) {
        if (is.url(nldi_feature)) {
          comid = sf::read_sf(nldi_feature)$nhdpv2_comid
        }
      } else {
        comid = nhdplusTools::get_nldi_feature(nldi_feature)$comid
      }
    }
    
    if (!is.null(xy)) {
      comid = nhdplusTools::discover_nhdplus_id(point = sf::st_sfc(sf::st_point(c(xy[1], xy[2])), crs = 4326))
    }
    
    if (!is.null(hl_uri) & !is.null(net)) {
      if (is.url(hl_uri)) {
        comid = sf::read_sf(hl_uri)$nhdpv2_comid
      } else {
        origin = dplyr::filter(net, hl_uri == !!hl_uri) |>
          dplyr::slice_max(hf_hydroseq) |>
          dplyr::pull(toid) |>
          unique()
      }
    }
    
    if (!is.null(comid) & !is.null(net)) {
      origin = dplyr::filter(net, hf_id == comid) |>
        dplyr::arrange(hf_hydroseq, hydroseq) |>
        dplyr::slice(1) |>
        dplyr::pull(id) |>
        unique()
    } else if (is.null(net)) {
      origin = comid
    }
    
    if (is.null(origin) | length(origin) > 1) {
      stop("Single origin not found")
      print(origin)
    }
    
    if (is.null(net)) {
      xx = suppressMessages({
        nhdplusTools::get_nhdplus(comid = comid)
      })
      
      v = nhdplusTools::vpu_boundaries
      
      vpuid = v$VPUID[which(lengths(sf::st_intersects(
        sf::st_transform(v, sf::st_crs(xx)), xx
      )) > 0)]
    } else {
      vpuid = unique(dplyr::pull(dplyr::filter(net, id == origin |
                                                 toid == origin), vpu))
    }
  }
  
  
  if (!is.null(base_dir)) {
    gpkg = dplyr::filter(base, grepl(vpuid, base) & grepl("gpkg", base))$file
  } else {
    gpkg = get_fabric(
      VPU = vpuid,
      base_s3 = base_s3,
      cache_dir = cache_dir,
      cache_overwrite = cache_overwrite
    )
  }
  
  if (length(gpkg) > 1 & !is.null(bbox)) {
    time = file.info(gpkg)$ctime
    order = rank(-rank(time))
    
    warning(
      "Multi-file subsets are not supported. You asked for:",
      paste0("\n\t", gpkg[order], " (", time[order], ")"),
      "\n\n",
      "Choosing the newest file only."
    )
    gpkg = gpkg[which.max(time)]
  }
  
  lyrs = lyrs[lyrs %in% sf::st_layers(gpkg)$name]
  
  db <- DBI::dbConnect(RSQLite::SQLite(), gpkg)
  on.exit(DBI::dbDisconnect(db))
  
  if (!is.null(bbox)) {
    return(subset_bbox(
      gpkg,
      bbox = bbox,
      lyrs = lyrs,
      outfile = outfile,
      qml_dir = qml_dir
    ))
    
  } else {
    
    if (!is.null(net)) {
      sub_net = dplyr::distinct(dplyr::select(dplyr::filter(net, vpu == vpuid), id, toid))
    } else {
      sub_net =     dplyr::tbl(db, lyrs[grepl("flowline|flowpath", lyrs)]) |>
        dplyr::select(dplyr::any_of(
          c(
            "id",
            "toid",
            'COMID',
            'toCOMID',
            "ID",
            "toID",
            "member_COMID"
          )
        )) |>
        dplyr::collect() |>
        dplyr::rename(dplyr::any_of(lookup))
      
      if ("member_COMID" %in% names(sub_net)) {
        origin =     dplyr::filter(sub_net, grepl(origin, member_COMID)) |>
          dplyr::pull(id)
        
        sub_net =     dplyr::select(sub_net, -member_COMID)
      }
    }
    
    message("Starting from: `",  origin, "`")
    
    tmap = suppressWarnings({
      nhdplusTools::get_sorted(dplyr::distinct(sub_net), outlets = origin)
    })
    
    if (grepl("nex", utils::tail(tmap$id, 1))) {
      tmap = utils::head(tmap, -1)
    }
  
    if(any(!is.null(areasqkm), 
           !is.null(pathlengthkm), 
           !is.null(ms_pathlengthkm))){
      
      print(nrow(tmap))

      tmap = area_length_filter(tmap = tmap,
                                gpkg = gpkg,
                                areasqkm = areasqkm,
                                pathlengthkm = pathlengthkm,
                                ms_pathlengthkm = ms_pathlengthkm)
      print(nrow(tmap))
    }
    
    ids = unique(c(unlist(tmap)))
    hydrofabric = list()
    
    for (j in 1:length(lyrs)) {
      message(glue::glue("Subsetting: {lyrs[j]} ({j}/{length(lyrs)})"))
      
      l  = sf::st_layers(gpkg)
      ind = which(l$name %in% lyrs[j])
      crs = l$crs[[ind]]
      
      t =     dplyr::tbl(db, lyrs[j]) |>
        dplyr::filter(dplyr::if_any(dplyr::any_of(
          c('COMID',  'FEATUREID', 'divide_id', 'id', 'ds_id', "ID")
        ), ~ . %in% !!ids)) |>
        dplyr::collect()
      
      if (all(!any(is.na(as.character(crs))), nrow(t) > 0)) {
        if (any(c("geometry", "geom") %in% names(t))) {
          t = sf::st_as_sf(t, crs = crs)
        } else {
          t = t
        }
      }
      
      if (!is.null(outfile)) {
        sf::write_sf(t, outfile, lyrs[j])
      } else {
        hydrofabric[[lyrs[j]]] = t
      }
    }
    
    if (is.null(cache_dir)) {
      unlink(gpkg)
    }
    
    if (!is.null(outfile)) {
      outfile = append_style(outfile, qml_dir = qml_dir, layer_names = lyrs)
      return(outfile)
    } else {
      return(hydrofabric)
    }
  }
}

#' Subset a Hydrofabric by Bounding Box
#' @param gpkg a path to a gpkg to subset
#' @inheritParams subset_network
#' @return file path (outfile) or list of features
#' @export

subset_bbox = function(gpkg,
                       bbox,
                       lyrs  = c(
                         "divides",
                         "nexus",
                         "flowpaths",
                         "network",
                         "hydrolocations",
                         "reference_flowline",
                         "reference_catchment",
                         "refactored_flowpaths",
                         "refactored_divides"),
                       outfile = NULL,
                       qml_dir = system.file("qml", package = "hydrofabric")) {
  lookup <-
    c(
      id = "ID",
      id = "COMID",
      toid = "toID",
      toid = "toCOMID"
    )
  
  x = sf::st_as_sfc(sf::st_bbox(bbox, crs = 4326))
  l = sf::st_layers(gpkg)
  lyrs = lyrs[lyrs %in% l$name]
  
  poly_ind = which(grepl("poly", tolower(l$geomtype)))
  
  db <- DBI::dbConnect(RSQLite::SQLite(), gpkg)
  on.exit(DBI::dbDisconnect(db))
  
  tmap = sf::read_sf(gpkg,
                     l$name[poly_ind],
                     wkt_filter = sf::st_as_text(sf::st_transform(x, l$crs[[poly_ind]]))) |>
    sf::st_drop_geometry() |>
    dplyr::select(dplyr::any_of(
      c("id",
        "toid",
        'COMID',
        'toCOMID',
        "ID",
        "toID",
        "member_COMID")
    )) |>
    dplyr::rename(dplyr::any_of(lookup))
  
  if (nrow(tmap) == 0) {
    stop('No data is provided bbox')
  }
  
  ids = unique(c(unlist(tmap)))
  hydrofabric = list()
  
  for (j in 1:length(lyrs)) {
    message(glue::glue("Subsetting: {lyrs[j]} ({j}/{length(lyrs)})"))
    
    ind = which(l$name %in% lyrs[j])
    crs = l$crs[[ind]]
    
    t = dplyr::tbl(db, lyrs[j]) |>
      dplyr::filter(dplyr::if_any(dplyr::any_of(
        c('COMID',  'FEATUREID', 'divide_id', 'id', 'ds_id', "ID")
      ), ~ . %in% !!ids)) |>
      dplyr::collect()
    
    if (all(!any(is.na(as.character(crs))), nrow(t) > 0)) {
      if (any(c("geometry", "geom") %in% names(t))) {
        t = sf::st_as_sf(t, crs = crs)
      } else {
        t = t
      }
    }
    
    if (!is.null(outfile)) {
      sf::write_sf(t, outfile, lyrs[j])
    } else {
      hydrofabric[[lyrs[j]]] = t
    }
  }
  
  if (!is.null(outfile)) {
    outfile = append_style(outfile, qml_dir = qml_dir, layer_names = lyrs)
    return(outfile)
  } else {
    return(hydrofabric)
  }
}


#' Area Length Thresholds
#' @param tmap Topologic Map
#' @param gpkg a geopackage
#' @inheritParams subset_network
#' @return data.frame

area_length_filter = function(tmap,
                              gpkg,
                              areasqkm = NULL,
                              pathlengthkm = NULL,
                              ms_pathlengthkm = NULL) {
  l = sf::st_layers(gpkg)
  lyr = l$name[which(grepl("line", tolower(l$geomtype)))]
  
  fps = sf::read_sf(gpkg, lyr) |>
    dplyr::select(id, lengthkm, areasqkm, mainstem)
  
  x = dplyr::left_join(tmap, sf::st_drop_geometry(fps), by = "id") |>
    dplyr::mutate(
      lengthkm = dplyr::coalesce(lengthkm, 0),
      areasqkm = dplyr::coalesce(areasqkm, 0)
    )
  
  x = x[nrow(x):1, ]
  
  print(nrow(x))
  
  if (!is.null(ms_pathlengthkm)) {
    opts = dplyr::filter(x, toid == utils::tail(x, 1)$toid)
    
    ms = dplyr::filter(x, mainstem %in% opts$mainstem) |>
         dplyr::group_by(mainstem) |>
         dplyr::summarise(sum = sum(lengthkm)) |>
         dplyr::slice_max(sum) |>
         dplyr::pull(mainstem)
    
    print(ms)
    o = dplyr::filter(x, mainstem == ms) |>
        dplyr::mutate(cl = cumsum(lengthkm)) |>
        dplyr::filter(cl <= ms_pathlengthkm)
    
  } else if (!is.null(pathlengthkm)) {
    o = dplyr::mutate(x, cl  = cumsum(lengthkm)) |>
        dplyr::filter(cl <= pathlengthkm)
  } else if (!is.null(areasqkm)) {
    o = dplyr::mutate(x, ca = cumsum(areasqkm)) |>
        dplyr::filter(ca <= !!areasqkm)
  }
  
  print(nrow(o))
  sub  = x[1:which(x$id == utils::tail(o, 1)$toid), ]
  print(nrow(sub))
  dplyr::filter(x, toid %in% unlist(dplyr::select(sub, id, toid))) |>
    dplyr::select(id, toid)

}

