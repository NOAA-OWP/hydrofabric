source("runners/config.R")
## TASK 1: build out uniform catchment distribution ----

base       <- "/Users/mikejohnson/hydrofabric"
version    <- "2.2"

base_rfc   <- glue("{base}/v{version}/refactor")
  full_hl    <- glue("{base}/v{version}/conus_hydrolocations")
  rfc_div    <- glue("{base_rfc}/conus_divides")
  rfc_fl     <- glue("{base_rfc}/conus_flowlines")
  rfc_net    <- glue('{base_rfc}/conus_network')

base_ref   <- glue("{base}/v{version}/reference")
  ref_div    <- glue('{base_ref}/conus_divides')

cw = read_parquet('/Volumes/MyBook/conus-hydrofabric/huc12_nhdplusv2_cw.parquet')

for (i in 1:nrow(pipeline)) {
  
  vpu = pipeline$vpus[i]
  
  divides = open_dataset(rfc_div) %>% 
    filter(vpuid == vpu) %>% 
    read_sf_dataset()
  
  network = open_dataset(rfc_net) %>% 
    filter(vpuid == vpu) %>% 
    select(id, member_comid) |> 
    distinct() |> 
    collect() |> 
    group_by(id) |> 
    mutate(member_comid = paste(member_comid, collapse = ",")) |> 
    distinct() |> 
    ungroup()
  
  flowpaths = open_dataset(rfc_fl) %>% 
    filter(vpuid == vpu) %>% 
    read_sf_dataset() |> 
    left_join(network, by = "id")

  hl_sub = open_dataset(rfc_net) %>% 
    filter(vpuid == vpu) %>% 
    distinct() |> 
    collect()

  gpkg = aggregate_to_distribution(
    divide = divides,
    flowpath = flowpaths, 
    vpu                    = vpu,
    outfile                = glue("/Users/mikejohnson/hydrofabric/v2.2/tmp-aggregate/{vpu}.gpkg"),
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




