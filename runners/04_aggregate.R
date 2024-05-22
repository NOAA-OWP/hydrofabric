source("runners/config.R")
## TASK 1: build out uniform catchment distribution ----

base       <- "/Volumes/MyBook/TNC"
version    <- "2.2"
base_rfc   <- glue("{base}/v{version}/refactor")
full_hl = '/Volumes/MyBook/conus-hydrofabric/v2.2/conus_hl'
rfc_div    <- glue("{base_rfc}/conus_divides")
rfc_fl     <- glue("{base_rfc}/conus_flowlines")

base_ref   <- glue("{base}/v{version}/reference")
ref_gpkg   <- glue('{base_ref}/reference_CONUS.gpkg')
ref_net    <- glue('{base_ref}/conus_network')
ref_div    <- glue('{base_ref}/conus_divides')
ref_fl     <- glue('{base_ref}/conus_flowlines')

cw = read_parquet(huc12_cw)

for (i in 1:nrow(pipeline)) {
  
  vpu = pipeline$vpus[i]
  
  divides = open_dataset(rfc_div) %>% 
    filter(vpuid == vpu) %>% 
    sfarrow::read_sf_dataset()
  
  flowpaths = open_dataset(rfc_fl) %>% 
    filter(vpuid == vpu) %>% 
    sfarrow::read_sf_dataset()

  hl_sub = open_dataset('/Volumes/MyBook/conus-hydrofabric/v2.2/conus_hl') %>% 
    filter(vpuid == vpu) %>% 
    select(poi_id, hf_id, geometry) %>% 
    sfarrow::read_sf_dataset() %>% 
    st_set_crs(5070) %>% 
    left_join(select(st_drop_geometry(flowpaths), id, hf_id),
              by = "hf_id",
              relationship = "many-to-many")
  
  gpkg = aggregate_to_distribution(
    divide = divides,
    flowpath = flowpaths, 
    vpu                    = vpu,
    outfile                = glue("{base}/v{version}/refactor/{vpu}_test.gpkg"),
    hydrolocations         = hl_sub,
    overwrite = TRUE)
  
  reference_divides = open_dataset(ref_div) %>% 
    filter(vpuid == vpu) %>% 
    sfarrow::read_sf_dataset() %>% 
    st_set_crs(5070)
  
  gpkg = add_nonnetwork_divides(gpkg, 
                                huc12 = cw, 
                                reference_divides =  reference_divides)
  
  message("Finished ", i, " of ", nrow(pipeline))
}


