#' @title Is URL
#' @description Identifies if a character string is a URL by searching 
#' for www, http, or https
#' @param x character string
#' @return boolean
#' @export

is.url <- function(x) {
  grepl("www.|http:|https:", x)
}

.subset = function(source = 's3://lynker-spatial/hydrofabric',
                   type = "reference",
                   hf_version = "2.2",
                   vpuid  = NULL,
                   outlet = NULL,
                   topo = NULL ){
 
    hook <- source("{source}/v{hf-version}/{type}/conus")
    
    d = open_dataset(glue("{hook}_network")) 
    
    net  = collect(select(filter(d, vpuid == vpuid), id, toid))
    
    subset = nhdplusTools::sort_network(net, outlets = outlet)
    
    ll = unique(as.vector(unlist(subset[-nrow(subset),])))
    
    fl = open_dataset(glue("{hook}_flowlines")) %>% 
      filter(vpuid == vpuid) %>% 
      filter(comid %in% ll) %>% 
      read_sf_dataset() %>% 
      st_set_crs(5070)
    
    div = open_dataset(glue("{hook}_divides")) %>% 
      filter(vpuid == vpuid) %>% 
      filter(featureid %in% ll) %>% 
      read_sf_dataset() %>% 
      st_set_crs(5070)
    
    net = open_dataset(glue("{hook}_network")) %>% 
      filter(vpuid == vpuid) %>% 
      filter(featureid %in% ll) %>% 
      read_sf_dataset() %>% 
      st_set_crs(5070)
    
    if(is.null(outfile)){
      list(divides = div,
           flowpaths = fl,
           network = net)
    } else {
      write_sf(fl, outfile, "flowpaths")
      write_sf(div, outfile, "divides")
      write_sf(net, outfile, "network")
      return(outfile)
    }

}


# subset_network = function(id = NULL,
#                           comid = NULL,
#                           hl_uri = NULL,
#                           nldi_feature = NULL,
#                           xy = NULL,
#                           bbox = NULL,
#                           type = "reference",
#                           version = "2.2", 
#                           source = "s3://lynker-spatial/hydrofabric",
#                           base_dir = NULL,
#                           lyrs  = c(
#                             "divides",
#                             "nexus",
#                             "flowpaths",
#                             "network",
#                             "hydrolocations",
#                             "flowpath_attributes",
#                             "reference_flowline",
#                             "reference_catchment",
#                             "refactored_flowpaths",
#                             "refactored_divides"
#                           ),
#                           areasqkm = NULL,
#                           pathlengthkm = NULL,
#                           ms_pathlengthkm = NULL,
#                           outfile = NULL,
#                           cache_dir = NULL,
#                           qml_dir = system.file("qml", package = "hydrofabric"),
#                           cache_overwrite = FALSE) {
#   Key <-
#     hf_hydroseq <-
#     hf_id  <- hydroseq <- member_COMID <- toid  <- vpu <- divide_id <- NULL
#   
#   lookup <-c(id = "ID", id = "COMID", toid = "toID", toid = "toCOMID")
#   
#   if (!is.null(bbox)) {
#     stopifnot(length(bbox) == 4)
#     x = sf::st_as_sfc(sf::st_bbox(bbox, crs = 4326))
#     vpuid = sf::st_filter(sf::st_transform(nhdplusTools::vpu_boundaries, 4326), x)$VPUID
#   } else {
#     if (!is.null(base_dir)) {
#       base = data.frame(file = list.files(base_dir,
#                                           full.names = TRUE,
#                                           recursive = TRUE)) |>
#         dplyr::mutate(base = basename(file))
#       
#       if ("conus_net.parquet" %in% base$base) {
#         net = filter(base, base == "conus_net.parquet") %>%
#           dplyr::pull(file) %>%
#           arrow::open_dataset() |>
#           dplyr::select(id, toid, divide_id, hf_id, hl_uri, hf_hydroseq, hydroseq, vpu) |>
#           dplyr::collect() |>
#           dplyr::distinct()
#       } else {
#         stop("conus_net.parquet not found in ", base_dir)
#       }
#     } else {
#       net = tryCatch({
#         arrow::open_dataset(glue::glue(base_s3, "conus_net.parquet")) |>
#           dplyr::select(id, divide_id, toid, hf_id, hl_uri, hf_hydroseq, hydroseq, vpu) |>
#           dplyr::distinct() |>
#           dplyr::collect() 
#           
#       } , error = function(e) {
#         NULL
#       })
#     }
#     
#     if (!is.null(id) & !is.null(net)) {
#       comid = dplyr::filter(net, id == !!id | toid == !!id | divide_id == !!id) |>
#         dplyr::slice_max(hf_hydroseq) |>
#         dplyr::pull(hf_id) |>
#         unique()
#     }
#     
#     if (!is.null(nldi_feature)) {
#       if (length(nldi_feature) == 1) {
#         if (is.url(nldi_feature)) {
#           comid = sf::read_sf(nldi_feature)$nhdpv2_comid
#         }
#       } else {
#         comid = nhdplusTools::get_nldi_feature(nldi_feature)$comid
#       }
#     }
#     
#     if (!is.null(xy)) {
#       comid = nhdplusTools::discover_nhdplus_id(point = sf::st_sfc(sf::st_point(c(xy[1], xy[2])), crs = 4326))
#     }
#     
#     if (!is.null(hl_uri) & !is.null(net)) {
#       if (is.url(hl_uri)) {
#         comid = sf::read_sf(hl_uri)$nhdpv2_comid
#       } else {
#         origin = dplyr::filter(net, hl_uri == !!hl_uri) |>
#           dplyr::slice_max(hf_hydroseq) |>
#           dplyr::pull(toid) |>
#           unique()
#       }
#     }
#     
#     if (!is.null(comid) & !is.null(net)) {
#       origin = dplyr::filter(net, hf_id == comid) |>
#         dplyr::arrange(hf_hydroseq, hydroseq) |>
#         dplyr::slice(1) |>
#         dplyr::pull(id) |>
#         unique()
#     } else if (is.null(net)) {
#       origin = comid
#     }
#     
#     if (is.null(origin) | length(origin) > 1) {
#       stop("Single origin not found")
#       print(origin)
#     }
#     
#     if (is.null(net)) {
#       xx = suppressMessages({
#         nhdplusTools::get_nhdplus(comid = comid)
#       })
#       
#       v = nhdplusTools::vpu_boundaries
#       
#       vpuid = v$VPUID[which(lengths(sf::st_intersects(
#         sf::st_transform(v, sf::st_crs(xx)), xx
#       )) > 0)]
#     } else {
#       vpuid = unique(dplyr::pull(dplyr::filter(net, id == origin |
#                                                  toid == origin), vpu))
#     }
#   }
#   
#   
#   if (!is.null(base_dir)) {
#     gpkg = dplyr::filter(base, grepl(paste0("_", vpuid, ".gpkg"), base))$file
#   } else {
#     gpkg = get_fabric(
#       VPU = vpuid,
#       base_s3 = base_s3,
#       cache_dir = cache_dir,
#       cache_overwrite = cache_overwrite
#     )
#   }
#   
#   if (length(gpkg) > 1 & !is.null(bbox)) {
#     time = file.info(gpkg)$ctime
#     order = rank(-rank(time))
#     
#     warning(
#       "Multi-file subsets are not supported. You asked for:",
#       paste0("\n\t", gpkg[order], " (", time[order], ")"),
#       "\n\n",
#       "Choosing the newest file only."
#     )
#     gpkg = gpkg[which.max(time)]
#   }
#   
#   lyrs = lyrs[lyrs %in% sf::st_layers(gpkg)$name]
#   
#   db <- DBI::dbConnect(RSQLite::SQLite(), gpkg)
#   on.exit(DBI::dbDisconnect(db))
#   
#   if (!is.null(bbox)) {
#     return(subset_bbox(
#       gpkg,
#       bbox = bbox,
#       lyrs = lyrs,
#       outfile = outfile,
#       qml_dir = qml_dir
#     ))
#     
#   } else {
#     
#     if (!is.null(net)) {
#       sub_net = dplyr::distinct(dplyr::select(dplyr::filter(net, vpu == vpuid), id, toid))
#     } else {
#       sub_net =     dplyr::tbl(db, lyrs[grepl("flowline|flowpath", lyrs)]) |>
#         dplyr::select(dplyr::any_of(
#           c(
#             "id",
#             "toid",
#             'COMID',
#             'toCOMID',
#             "ID",
#             "toID",
#             "member_COMID"
#           )
#         )) |>
#         dplyr::collect() |>
#         dplyr::rename(dplyr::any_of(lookup))
#       
#       if ("member_COMID" %in% names(sub_net)) {
#         origin =     dplyr::filter(sub_net, grepl(origin, member_COMID)) |>
#           dplyr::pull(id)
#         
#         sub_net =     dplyr::select(sub_net, -member_COMID)
#       }
#     }
#     
#     message("Starting from: `",  origin, "`")
#     
#     tmap = suppressWarnings({
#       nhdplusTools::get_sorted(dplyr::distinct(sub_net), outlets = origin)
#     })
#     
#     if (grepl("nex", utils::tail(tmap$id, 1))) {
#       tmap = utils::head(tmap, -1)
#     }
#   
#     if(any(!is.null(areasqkm), 
#            !is.null(pathlengthkm), 
#            !is.null(ms_pathlengthkm))){
#       
#       print(nrow(tmap))
# 
#       tmap = area_length_filter(tmap = tmap,
#                                 gpkg = gpkg,
#                                 areasqkm = areasqkm,
#                                 pathlengthkm = pathlengthkm,
#                                 ms_pathlengthkm = ms_pathlengthkm)
#       print(nrow(tmap))
#     }
#     
#     ids = unique(c(unlist(tmap)))
#     hydrofabric = list()
#     
#     for (j in 1:length(lyrs)) {
#       message(glue::glue("Subsetting: {lyrs[j]} ({j}/{length(lyrs)})"))
#       
#       l  = sf::st_layers(gpkg)
#       ind = which(l$name %in% lyrs[j])
#       crs = l$crs[[ind]]
#       
#       t =     dplyr::tbl(db, lyrs[j]) |>
#         dplyr::filter(dplyr::if_any(dplyr::any_of(
#           c('COMID',  'FEATUREID', 'divide_id', 'id', 'ds_id', "ID")
#         ), ~ . %in% !!ids)) |>
#         dplyr::collect()
#       
#       if (all(!any(is.na(as.character(crs))), nrow(t) > 0)) {
#         if (any(c("geometry", "geom") %in% names(t))) {
#           t = sf::st_as_sf(t, crs = crs)
#         } else {
#           t = t
#         }
#       }
#       
#       if (!is.null(outfile)) {
#         sf::write_sf(t, outfile, lyrs[j])
#       } else {
#         hydrofabric[[lyrs[j]]] = t
#       }
#     }
#     
#     if (is.null(cache_dir)) {
#       if(is.null(base_dir)){
#         unlink(gpkg)
#       }
#     }
#     
#     if (!is.null(outfile)) {
#       outfile = append_style(outfile, qml_dir = qml_dir, layer_names = lyrs)
#       return(outfile)
#     } else {
#       return(hydrofabric)
#     }
#   }
# }

#' Subset a Hydrofabric by Bounding Box
#' @param gpkg a path to a gpkg to subset
#' @inheritParams get_subset
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
#' @inheritParams get_subset
#' @return data.frame

area_length_filter = function(tmap,
                              gpkg,
                              areasqkm = NULL,
                              pathlengthkm = NULL,
                              ms_pathlengthkm = NULL) {
  
  lengthkm <- mainstem <- toid <-  cl <-  ca <-  NULL
  
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


#' Enhance Hydrofabric Network
#' @param gpkg hydrofabric id (relevant only to nextgen fabrics)
#' @param base_s3 the base hydrofabric directory to access in Lynker's s3
#' @param lyrs layers to extract. Default is "forcing_weights" and "model_attributes"
#' @return file path of gpkg
#' @export

add_parameters = function(gpkg = NULL,
                          lyrs = c("forcing_weights", "model_attributes"),
                          base_s3 = "s3://lynker-spatial/v20.1"){
  
  divide_id <-  NULL
  
  div = read_sf(gpkg, "divides")
  
  if("forcing_weights" %in% lyrs){
    open_dataset(glue('{base_s3}/forcing_weights.parquet')) |>
      filter(divide_id %in% !!div$divide_id) |>
      collect() |>
      write_sf(gpkg, "forcing_weights")
  }
  
  if("model_attributes" %in% lyrs){
    open_dataset(glue('{base_s3}/model_attributes.parquet')) |>
      filter(divide_id %in% !!div$divide_id) |>
      collect() |>
      write_sf(gpkg, "model_attributes")
  }
  
  return(gpkg)
}

