source("runners/config.R")
library(tidyr)
library(AOI)

# Purpose --> This code constructs a conus_hl (conus hydrolocation)file for a set of inputs,
# These include:
#   Reference Fabric POIs
#   Coastal doamin (1) gages and (2) boundaries
#   AHPS gages from FIM
#   Routelink information from the NWM

# Extract and Build POIs --------------------------------------------------
#   Work across VPU set

for (i in seq_along(vpus)) {
  VPU <- vpus[i]
  
  refactored_gpkg <- get_hydrofabric(VPU = VPU, type = "refactor", dir = base_refactored)
  reference_gpkg  <- get_hydrofabric(VPU = VPU, type = "reference", dir = base_reference)
  
  # These are our source flowlines from the reference fabric
  fl <- read_hydrofabric(reference_gpkg, "flowlines")[[1]] |>
    st_transform(4326)
  
  # Identify the name of the POI layer in the reference fabric files
  poi_layer <- grep("POIs_*", st_layers(reference_gpkg)$name, value = TRUE)
  
  # Community POIs  ----
  # Extract POIs from reference fabric and pivot wide table to long representation
  hl  <- read_sf(reference_gpkg, poi_layer[!grepl("tmp", poi_layer)]) |>
    #filter(Type_NID == 'CA00775') %>% 
    st_transform(4326) %>%
    select(hf_id = COMID, poi_id = id,  any_of(paste0("Type_", community_hl_types))) |>
    mutate_at(vars(matches("Type_")), as.character) |>
    pivot_longer(-c(poi_id, hf_id, geom), names_to = "hl_reference", values_to = "hl_link") |>
    filter(!is.na(hl_link)) |>
    distinct() |>
    mutate(hl_reference = gsub("Type_", "", hl_reference),
           source = "reference_pois") |>
    rename_geometry("geometry")
  
  # Coastal Model ----
  # Coastal Gages and domain edges are added!
  # Starting with the Gages
  coastal_network <- read.csv(coastal_gages) |>
    select(hl_link = SITE_NO, lat = LAT_NHD, lon = LON_NHD) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326)  |>
    st_join(st_transform(select(vpu_boundaries, VPUID), 4326)) |>
    filter(VPUID == VPU) |>
    mutate(hl_reference = "CoastalGage",
           hl_link = sprintf("%08s",  hl_link),
           source = "coastal_csv",
           VPUID = NULL) |>
    rename_geometry("geometry") %>% 
    st_join(select(fl, hf_id = COMID), join = st_nearest_feature)

  # Moving to domain intersection
  coastal_domain <- read_sf("runners/data/coastal_domain.gpkg") |>
    st_transform(st_crs(fl)) |>
    st_cast("MULTILINESTRING") |>
    st_intersection(bbox_get(fl))
  
  touches <- st_intersects(fl, coastal_domain)
  
  fl2 <- fl[lengths(touches) > 0, ]
  
  if (nrow(fl2) == 0) {
    coastal_domain_pois <- NULL
  } else {
    coastal_domain_pois <- st_intersection(coastal_domain, fl2) |>
      select(hf_id = COMID, domain_id) |>
      st_collection_extract("POINT") |>
      st_cast("POINT") |>
      distinct() |>
      group_by(domain_id) |>
      mutate(hl_link = paste0(domain_id, "-", 1:n())) |>
      ungroup() |>
      mutate(
        hl_reference = paste0("CoastalDomain-", domain_id),
        domain_id = NULL,
        source = "coastal_domain"
      ) |>
      rename_geometry("geometry")
  }
  
  #Consolidating Coastal POIs
  coastal = bind_rows(coastal_domain_pois, coastal_network)
  
  # RouteLink ----
  
  lu <- read_sf(refactored_gpkg, "lookup_table") |>
    select(comid = member_COMID, hy_id = reconciled_ID) |>
    mutate(comid = as.numeric(comid))
  
  rl <- open_dataset(rl_file) |>
          select(comid, NHDWaterbodyComID) |>
          filter(comid %in% fl$COMID,
                 !is.na(NHDWaterbodyComID)) |>
        collect() |>
        left_join(get_vaa("hydroseq"), by = "comid") |>
        filter(!is.na(hydroseq)) |>
        group_by(NHDWaterbodyComID) |>
        # Max HF is the outlet
        slice_max(hydroseq) |>
        ungroup() |>
        select(hf_id = comid, hl_link = NHDWaterbodyComID) |>
        mutate(
          hl_reference = "RL_WBOut",
          hl_link = as.character(hl_link),
          source = "routelink"
        ) 
       
  rl = filter(fl, COMID %in% rl$hf_id) |> 
    select(hf_id = COMID) %>% 
    st_set_geometry(st_geometry(get_node(.,"end"))) |> 
    right_join(rl) %>% 
    rename_geometry("geometry") #%>% 
    #filter(!hf_id %in% filter(hl, hl_reference == "WBOut")$hf_id)
  
  # AHPS ----
  ahps <- read_sf(fim_ahps) |>
    filter(substring(HUC8, 1, 2) == substring(VPU, 1, 2)) |>
    filter(!usgs_site_code %in% filter(hl, hl_reference == "Gages")$hl_link) |>
    st_join(st_transform(select(vpu_boundaries, "VPUID"), 5070)) |>
    filter(VPUID == VPU) |>
    st_transform(4326) %>%
    mutate(hf_id = nwm_feature_id) |>
    select(nws_lid, usgs_site_code, hf_id) |>
    pivot_longer(-c(hf_id, geom)) %>%
    filter(!is.na(value)) |>
    mutate(
      name = ifelse(name == "usgs_site_code", "Gages", name),
      name = ifelse(name == "nws_lid", "AHPS", name),
      source = "ahps",
      hf_id = as.numeric(hf_id)
    ) |>
    rename(hl_link = value, hl_reference = name) %>% 
    rename_geometry("geometry")
  
  # Merge and Build ----
  
  tmp1 <-
    bind_rows(hl,
              coastal,
              rl,
              ahps) %>%
    mutate(X = st_coordinates(.)[,1],
           Y = st_coordinates(.)[,2]) |>
    group_by(hf_id, X, Y) |>
    mutate(n = n())
  
   u <- unique(filter(tmp1, n > 2)$hf_id)
   
   tmp <- st_join(tmp1, select(filter(fl, COMID %in% u), COMID)) %>% 
     filter(hf_id == COMID) |>
     bind_rows(filter(tmp1, !hf_id %in% u)) %>% 
     mutate(n = NULL, COMID = NULL) %>% 
     group_by(hf_id) |>
     mutate(hl_id = paste0(VPU, cur_group_id())) |>
     ungroup() |>
     distinct()
    
  lu <- read_sf(refactored_gpkg, "lookup_table") |>
    select(hf_id = NHDPlusV2_COMID,
           ID    = reconciled_ID,
           mainstem = LevelPathID,
           member_COMID) |>
    group_by(hf_id) |>
    slice_max(as.numeric(member_COMID)) |>
    ungroup() |>
    select(-member_COMID)

    left_join(tmp, lu, by = "hf_id") |>
    mutate(hl_position = "outflow", 
           VPUID = VPU, 
           hf_source = "nhdplusv2") |>
    select(hl_id, poi_id, hl_reference, hl_link, hl_source = source, hl_position, 
           hf_id, hf_source, ID, mainstem, X, Y, VPUID) %>% 
    write_sf(glue("{base_hl}/hl_{VPU}.gpkg"))
  
}


l <- bind_rows(lapply(list.files(base_hl, full.names = TRUE), read_sf))

table(l$hl_reference) %>% 
  hist()

library(ggplot2)
ggplot() + 
  geom_bar(data = l, aes(x = reorder(hl_reference, -table(hl_reference)[hl_reference]))) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "", title = paste0('Hydrolocations: ', length(unique(l$hl_id)), ", Indexed featrues: ", nrow(l)))

write_sf(l, full_hl)
