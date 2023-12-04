# 4608*4
# 3841 *4
# 
# gdalwarp -t_srs '+proj=lcc +lat_0=40.000008 +lon_0=-97 +lat_1=30 +lat_2=60 +x_0=0 +y_0=0 +R=6370000 +units=m +no_defs' -te -2303999.62876143 -1920000.70008381 2304000.37123857 1919999.29991619 -ts 18432 15364 -r average /Volumes/MyBook/3DEP/usgs_1.vrt /Volumes/MyBook/3DEP/usgs_1_250mm.tif
# 
# gdaldem slope /Volumes/MyBook/3DEP/usgs_1_250mm.tif /Volumes/MyBook/3DEP/usgs_250m_slope.tif
# gdaldem aspect /Volumes/MyBook/3DEP/usgs_1_250mm.tif /Volumes/MyBook/3DEP/usgs_250m_aspect.tif
# 
# gdalwarp -t_srs '+proj=lcc +lat_0=40.000008 +lon_0=-97 +lat_1=30 +lat_2=60 +x_0=0 +y_0=0 +R=6370000 +units=m +no_defs' -te -2303999.62876143 -1920000.70008381 2304000.37123857 1919999.29991619 -ts 18432 15364 -r average /Volumes/MyBook/nlcd/nlcd_2021_impervious_l48_20230630/nlcd_2021_impervious_l48_20230630.img /Volumes/MyBook/imperv_250m.tif
# 
# #https://zenodo.org/records/4460354
# gdalwarp -t_srs '+proj=lcc +lat_0=40.000008 +lon_0=-97 +lat_1=30 +lat_2=60 +x_0=0 +y_0=0 +R=6370000 +units=m +no_defs' -te -2303999.62876143 -1920000.70008381 2304000.37123857 1919999.29991619 -ts 18432 15364 -r average /Volumes/MyBook/CONUS_TWI_epsg5072_30m_unmasked.tif /Volumes/MyBook/twi_250m.tif

source('runners/config.R')
library(powerjoin)

for(i in 1:nrow(pipeline)){
  
  if(!file.exists(pipeline$cfe[i])){
    o    = add_cfe_noahowp_attributes(gpkg        = pipeline$nextgen[i],
                                      nwm_dir     = nwm_dir,
                                      outfile     = pipeline$cfe[i])
  }

  message(i)
}

for(i in 1:nrow(pipeline)){
  
  if(!file.exists(pipeline$atts[i])){
    message("---- VPU: ", pipeline$vpu[i],  "----")
    message("---- Data Preparation ----")
    div = read_sf(pipeline$nextgen[i], 'divides')
    
    message("---- Weight Grid ----")
    w = weight_grid(rast(dem_file), div, ID = "divide_id")
    
    message("---- mean variables ----")
    x1 = execute_zonal(data = c(rast(dem_file), rast(slope_file), rast(imp_file)), 
                              w = w,
                              ID = "divide_id", 
                              join = FALSE)
    names(x1) = c("divide_id", "elevation_mean", "slope_mean", "impervious_mean")
    
    message("---- c. mean variables ----")
    x2 = execute_zonal(data = rast(aspect_file), 
                              w = w,
                              ID = "divide_id", 
                              fun = circular_mean,
                              join = FALSE)
    names(x2) = c("divide_id", "aspect_c_mean")
    
    message("---- TWI ----")
    x3 = execute_zonal(data = rast(twi_file), 
                              w = w,
                              ID = "divide_id", 
                              fun = equal_population_distribution,
                              groups = 4,
                              join = FALSE)
    names(x3) = c("divide_id", "twi_dist_4")
    
    x4 = st_centroid(div$geom) %>% 
      st_transform(4326) %>% 
      st_coordinates() %>% 
      as.data.frame() %>% 
      mutate(divide_id = div$divide_id)
    
    message("---- Old variables ----")
    vars = read_parquet(pipeline$cfe[i])
    
    message("---- Write out data ----")
    out = power_full_join(list(x1,x2,x3,x4, vars), by = "divide_id")
    write_parquet(out, pipeline$atts[i])
    
  }
  
}
  
unlink(glue("{base}/cfe"), recursive = TRUE)
unlink(glue('{base}/{version}/fgb/hydrolocations.fgb'))
unlink(glue('{base}/{version}/fgb/nexus.fgb'))
unlink(glue('{base}/{version}/fgb/flowpaths.fgb'))
unlink(glue('{base}/{version}/fgb/divides.fgb'))

open_dataset(glue('{base}/{version}/model_attributes')) %>% 
  collect() %>% 
  write_parquet(glue('{base}/{version}/model_attributes.parquet'))

  dir.create(glue('{base}/{version}/fgb/'), recursive = TRUE, showWarnings = FALSE)

hl = read_sf(glue('{base}/{version}/conus.gpkg'), "hydrolocations")
  write_sf(hl, glue('{base}/{version}/fgb/hydrolocations.fgb'))
  
nexus = read_sf(glue('{base}/{version}/conus.gpkg'), "nexus")
  write_sf(nexus, glue('{base}/{version}/fgb/nexus.fgb'))
  
fps = read_sf(glue('{base}/{version}/conus.gpkg'), "flowpaths")
  write_sf(fps, glue('{base}/{version}/fgb/flowpaths.fgb'))
  
divides = read_sf(glue('{base}/{version}/conus.gpkg'), "divides")
  write_sf(divides, glue('{base}/{version}/fgb/divides.fgb'))

  
  