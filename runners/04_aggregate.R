source("runners/config.R")
## TASK 1: build out uniform catchment distribution ----
cw = read_parquet(huc12_cw)

if(FIX){ pipeline = filter(pipeline, !is.na(corrected_refactor)) }

for (i in 1:nrow(pipeline)) {
  
  hl = filter(read_sf(full_hl), VPUID == pipeline$vpus[i]) %>% 
    mutate(hl_position = "outflow")
  
  gpkg = ifelse(is.na(pipeline$corrected_refactor[i]), 
                pipeline$refactored_gpkg[i], 
                pipeline$corrected_refactor[i])
    
  gpkg = aggregate_to_distribution(
    gpkg                   = gpkg,
    vpu                    = pipeline$vpus[i],
    divide                 = 'refactored_divides',
    outfile                = pipeline$uniform[i],
    hydrolocations         = hl,
    overwrite = TRUE)
  
  gpkg = add_nonnetwork_divides(gpkg, 
                                huc12 = cw, 
                                reference_gpkg = pipeline$reference_gpkg[i])
  
  message("Finished ", i, " of ", nrow(pipeline))
}


