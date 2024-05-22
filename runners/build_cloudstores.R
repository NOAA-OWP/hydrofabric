ref_dir <- '/Volumes/MyBook/conus-hydrofabric/reference_05012024'


# Runner ------------------------------------------------------------------

vars    <- c("catchments", "flowlines", "waterbodies")
vpus     <- nhdplusTools::vpu_boundaries$VPUID[1:21]

dir.create(glue('{ref_dir}/geoparquet'), showWarnings = FALSE)

for(i in 1:length(vars)){
  dir.create(glue('{ref_dir}/geoparquet/{vars[i]}'), showWarnings = FALSE)
}

for(i in 1:length(vpus)) {
  for (j in 1:length(vars)) {
    system(
      glue(
        'ogr2ogr -f parquet {ref_dir}/geoparquet/{vars[j]}/{vpus[i]}_{vars[j]}.parquet {ref_dir}/gpkg/{vpus[i]}_reference_features.gpkg {vars[j]}'
      )
    )
  }
}

n = arrow::read_parquet('/Volumes/Transcend/reference_geometries/enhd_nhdplusatts.parquet') |>
  select(comid, tocomid, hydroseq, levelpathi, vpuid)

# comid
# featureid
# tocomid
# hf_id = comid
# hl_uri = poi_id
# hf_hydroseq = hydroseq
# hydroseq = hydroseq
# vpu = hydroseq
# 
library(remotes)
install_github("r-spatial/sf", configure.args = "--with-proj-lib=$(brew --prefix)/lib/")

nets = list()

for(i in 1:length(vpus)){
  
  n2 = filter(n, vpuid == vpus[i])
  nets[[i]] <- arrow::open_dataset(glue('{ref_dir}/geoparquet/catchments/{vpus[i]}_catchments.parquet')) %>% 
    select(featureid) %>% 
    collect() %>%
    mutate(comid = ifelse(featureid %in% n2$comid, featureid, NA)) %>% 
    full_join(n2, by = "comid") %>% 
    mutate(vpuid == vpus[i])
}

network = bind_rows(nets)

arrow::write_parquet(network, glue::glue('{ref_dir}/conus_net.parquet'))


