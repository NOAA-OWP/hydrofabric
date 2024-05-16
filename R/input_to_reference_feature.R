# A function to find the most terminal feature given ids or coordinates
#' @title Find a Reference Feature
#' @param network table from network file default is NULL. datatype: dataframe e.g., conus_network
#' @param id hydrofabric id. datatype: string / vector of strings e.g., 'wb-10026' or c('wb-10026', 'wb-10355') 
#' @param comid NHDPlusV2 COMID. datatype: int / vector of int e.g., 61297116 or c(61297116 , 6129261) 
#' @param hl_id hydrolocation id. datatype: int / vector of int e.g., 01236 or c(01236 , 01244) 
#' @param hl_uri hydrolocation URI. datatype: string / vector of string / a url e.g., HUC12-010100100101 or c(HUC12-010100100101 , HUC12-010100110104) 
#' @param poi_id POI identifier. datatype: int / vector of int e.g., 266387 or c(266387, 266745)
#' @param nldi_feature list with names 'featureSource' and 'featureID' where 'featureSource' is derived from the "source" column of the response of dataRetrieval::get_nldi_sources() and the 'featureID' is a known identifier from the specified 'featureSource'. datatype: a url e.g., 'https://labs.waterdata.usgs.gov/api/nldi/linked-data/census2020-nhdpv2'
#' @param xy Location given as vector of XY and CRS (e.g., 4326) (longitude, latitude, crs)
#' @inheritParams get_vpu_fabric


input_to_reference_feature = function(id = NULL, 
                                      comid = NULL,  
                                      hl_id = NULL, 
                                      hl_uri = NULL, 
                                      poi_id = NULL, 
                                      nldi_feature = NULL, 
                                      xy = NULL, 
                                      type = "reference",
                                      version = "2.2", 
                                      source = "s3://lynker-spatial/hydrofabric"
                                      ) {
  
  
  # ____ NOTES ________ # 
  # 
  # id = 101
  # network_dir = 's3://lynker-spatial/hydrofabric/v2.2/reference/conus_network/'
  # if(!is.null(id)){
  #   net = open_dataset(network_dir) %>% 
  #     filter(comid == !!id) %>% 
  #     select(vpuid, id) %>% 
  #     distinct() %>% 
  #     collect()
  # }
  # 
  # if(!is.null(hl_id)){
  #   open_dataset('/Volumes/MyBook/conus-hydrofabric/v2.2/conus_hl') %>% 
  #     filter(hl_id == !!hl_id) %>% 
  #     select(vpuid, poi_id) %>% 
  #     distinct() %>% 
  #     collect()
  # }
  # 
  # if(!is.null(poi_id)){
  #   open_dataset(network_dir) %>% 
  #     filter(poi_id == !!poi_id) %>% 
  #     select(vpuid, comid) %>% 
  #     distinct() %>% 
  #     collect()
  # }
  # 
  # if(net$topo == "fl-fl"){
  #   outlet = net$id
  # } else if(net$topo == "fl-nex"){
  #   outlet = net$toid
  # }
  # 
  # return(vpu = net$vpuid, outlet = outlet)
  # Initialize as Null variables an default varaibles
  # 
  toid <- divide_id <- hf_hydroseq <- hf_id <- hydroseq <- vpu <- NULL

  # Network present -------------------------------------------------
  # If network file is present return origin such that it includes id, toid, hf_id, and vpu 
  # as dataframe it works with single variable input and vector 

  # Given an id and network file 
  if (!is.null(id) & !is.null(net)) {
    # Check if input is a string or a vector
    if (!is.character(id) && !is.vector(id)) {
        stop("id input must be a single string or a vector of text elements.")
    }
    origin <- dplyr::filter(net, id %in% !!id) |>
              group_by(id) %>%
              mutate(max_hf_hydroseq = max(hf_hydroseq, na.rm = TRUE)) %>%
              filter(hf_hydroseq == max_hf_hydroseq) %>%
              distinct(id, .keep_all = TRUE) %>%
              select(id, toid, hf_id, vpu) %>%
              as.data.frame()
    if (is.null(origin)) {
      stop("Single origin not found")
    }
    return(origin)
  }

  # Given an comid and network file 
  if (!is.null(comid) & !is.null(net)) {
    # Cast to int
    comid <- as.integer(comid)
    origin <- dplyr::filter(net, comid %in% !!hf_id) |>
              group_by(id) %>%
              mutate(max_hf_hydroseq = max(hf_hydroseq, na.rm = TRUE)) %>%
              filter(hf_hydroseq == max_hf_hydroseq) %>%
              distinct(id, .keep_all = TRUE) %>%
              select(id, toid, hf_id, vpu) %>%
              as.data.frame()
    if (is.null(origin)) {
      stop("Single origin not found")
    }
    return(origin)
  }

  # Given an hl_id and network file 
  if (!is.null(hl_id) & !is.null(net)) {
    origin <- dplyr::filter(net, hl_id %in% !!hl_id) |>
              group_by(id) %>%
              mutate(max_hf_hydroseq = max(hf_hydroseq, na.rm = TRUE)) %>%
              filter(hf_hydroseq == max_hf_hydroseq) %>%
              distinct(id, .keep_all = TRUE) %>%
              select(id, toid, hf_id, vpu) %>%
              as.data.frame()
    if (is.null(origin)) {
      stop("Single origin not found")
    }
    return(origin)
  }

  # Given an hl_id and network file 
  if (!is.null(hl_id) & !is.null(net)) {
    origin <- dplyr::filter(net, hl_id %in% !!hl_id) |>
              group_by(id) %>%
              mutate(max_hf_hydroseq = max(hf_hydroseq, na.rm = TRUE)) %>%
              filter(hf_hydroseq == max_hf_hydroseq) %>%
              distinct(id, .keep_all = TRUE) %>%
              select(id, toid, hf_id, vpu) %>%
              as.data.frame()
    if (is.null(origin)) {
      stop("Single origin not found")
    }
    return(origin)
  }

  # Given a hl_url and network file 
  if (!is.null(hl_uri) & !is.null(net)) {
    # If url  and comid not present then retrive the comid first and keep the most terminal
    if (is.url(hl_uri)) {
      if (is.null(comid)) {
          # Retrive the comid for features
          comid = unique(suppressWarnings(sf::read_sf(nldi_feature) %>%
                  select(one_of("COMID", "comid", "nhdpv2_comid")) %>%
                  unlist(use.names = FALSE)))        
          comid <- comid[comid != ""]
      }
      # Filter net file to those comids and keep the most terminal
      origin <- net[net$hf_id %in% comid, ] |>
                group_by(id) %>%
                mutate(max_hf_hydroseq = max(hf_hydroseq, na.rm = TRUE)) %>%
                filter(hf_hydroseq == max_hf_hydroseq) %>%
                distinct(id, .keep_all = TRUE) %>%
                select(id, toid, hf_id, vpu) %>%
                as.data.frame()
      if (is.null(origin)) {
        stop("Single origin not found")
      }
      return(origin)
    } else {
        origin <- dplyr::filter(net, hl_uri %in% !!hl_uri) |>
                  group_by(id) %>%
                  mutate(max_hf_hydroseq = max(hf_hydroseq, na.rm = TRUE)) %>%
                  filter(hf_hydroseq == max_hf_hydroseq) %>%
                  distinct(id, .keep_all = TRUE) %>%
                  select(id, toid, hf_id, vpu) %>%
                  as.data.frame()
    }
    if (is.null(origin)) {
      stop("Single origin not found")
    }
    return(origin)
  }
  
  # Given an poi id
  if (!is.null(poi_id) & !is.null(net)) {
    # Cast to int
    poi_id <- as.integer(poi_id)
    
    # Read network and grab comids
    hook = glue("{s3}/v{version}/{type}/conus")
    net_data = arrow::open_dataset(glue("{hook}_network")) 
    comids <- net_data %>%
              filter(poi_id %in% poi_id) %>%
              select(poi_id, comid)%>%
              as.data.frame()

    filtered_net <- net %>%
                    filter(hf_id %in% comids$comid)

    # Add poi_id column
    merged_df <- merge(filtered_net, comids, by.x = "hf_id", by.y = "comid")

    origin <- merged_df |>
              group_by(poi_id) %>%
              mutate(max_hf_hydroseq = max(hf_hydroseq, na.rm = TRUE)) %>%
              filter(hf_hydroseq == max_hf_hydroseq) %>%
              distinct(id, .keep_all = TRUE) %>%
              as.data.frame()%>%
              select(id, toid, hf_id, vpu) 
    return(origin)
  }

  # Given an nldi_feature
  if (!is.null(nldi_feature) & !is.null(net)) {
    # Retrive comid using nldi_feature
    if (length(nldi_feature) == 1) {
      # Pull url metadata if given a url
      if (is.url(nldi_feature)) {
          # comid column name can change depending on the dataset being pulled check for all
          comid = unique(suppressWarnings(sf::read_sf(nldi_feature) %>%
                  select(one_of("COMID", "comid", "nhdpv2_comid")) %>%
                  unlist(use.names = FALSE)))
          comid <- comid[comid != ""]
      }
    # Otherwise use nldi feature id to pull metadata
    } else {
    comid = unique(suppressWarnings(nhdplusTools::get_nldi_feature(nldi_feature) %>%
                select(one_of("COMID", "comid", "nhdpv2_comid")) %>%
                unlist(use.names = FALSE)))
    comid <- comid[comid != ""]
    }
    origin <- dplyr::filter(net, comid %in% !!hf_id) |>
              group_by(id) %>%
              mutate(max_hf_hydroseq = max(hf_hydroseq, na.rm = TRUE)) %>%
              filter(hf_hydroseq == max_hf_hydroseq) %>%
              distinct(id, .keep_all = TRUE) %>%
              select(id, toid, hf_id, vpu) %>%
              as.data.frame()
    if (is.null(origin)) {
      stop("Single origin not found")
    }
    return(origin)
  }

  # Given a coordiantes 
  if (!is.null(xy) & !is.null(net)) {
    # cast coordiantes to numeric and crs to int
    xy[1:2] <- as.numeric(xy[1:2])
    xy[3] <- as.integer(xy[3])
    comid = nhdplusTools::discover_nhdplus_id(point = sf::st_sfc(sf::st_point(c(xy[1], xy[2])), crs = xy[3]))
    origin <- dplyr::filter(net, comid %in% !!hf_id) |>
              group_by(id) %>%
              mutate(max_hf_hydroseq = max(hf_hydroseq, na.rm = TRUE)) %>%
              filter(hf_hydroseq == max_hf_hydroseq) %>%
              distinct(id, .keep_all = TRUE) %>%
              select(id, toid, hf_id, vpu) %>%
              as.data.frame()
    if (is.null(origin)) {
      stop("Single origin not found")
    }
    return(origin)
  }
  
  # Network not present -------------------------------------------------
  # If network file is not present return origin such that it includes comid and vpu 
  # as dataframe it works with single variable input and vector 

  if (is.null(net)) {
    # Check if no comid is provided then compute one
    if (is.null(comid)) {
      # If nldi_feature is provided
      if (!is.null(nldi_feature)) {
        # Retrive the comid for features
        comid = unique(suppressWarnings(sf::read_sf(nldi_feature) %>%
                select(one_of("COMID", "comid", "nhdpv2_comid")) %>%
                unlist(use.names = FALSE)))        
        comid <- comid[comid != ""]
      # Else if poi is provided
      } else if (!is.null(poi_id)) {
        # Cast to int
        poi_id <- as.integer(poi_id)
        
        # Read network and grab comids
        hook = glue("{s3}/v{version}/{type}/conus")
        net_data = arrow::open_dataset(glue("{hook}_network")) 
        poi_df <- net_data %>%
                  filter(poi_id %in% poi_id) %>%
                  as.data.frame()

        comid <- poi_df |>
              group_by(poi_id) %>%
              mutate(max_hf_hydroseq = max(hydroseq, na.rm = TRUE)) %>%
              filter(hydroseq == max_hf_hydroseq) %>%
              dplyr::pull(comid) |>
              unique()
      # Then xy is provided
      } else {
        # cast coordiantes to numeric and crs to int
        xy[1:2] <- as.numeric(xy[1:2])
        xy[3] <- as.integer(xy[3])
        comid = nhdplusTools::discover_nhdplus_id(point = sf::st_sfc(sf::st_point(c(xy[1], xy[2])), crs = xy[3]))
      }
    }

    nhd_plus = suppressMessages({
      nhdplusTools::get_nhdplus(comid = comid)
    })
    
    vpu_bounds = nhdplusTools::vpu_boundaries
    
    vpuid = vpu_bounds$VPUID[which(lengths(sf::st_intersects(
      sf::st_transform(vpu_bounds, sf::st_crs(nhd_plus)), nhd_plus
    )) > 0)]

    origin <- data.frame(hf_id = comid, vpu = vpuid)
    origin$id <- NaN
    origin$toid <- NaN
    return(origin)
  } 
}

