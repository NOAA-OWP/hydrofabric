# Dominant Soil Type [ISLTYP]

vars = c('bexp', 'cwpvt', 'dksat', 'ISLTYP', 'IVGTYP', 
         'mfsno', 'mp', 'psisat', 'refkdt', 'slope', 
         'smcmax', 'smcwlt', 'vcmx25')


soil = ngen.hydrofab::correct_nwm_spatial('/Volumes/MyBook/conus-hydrofabric/nwm_v3.0.7/conus/soilproperties_CONUS_FullRouting.nc')
wrf  = ngen.hydrofab::correct_nwm_spatial('/Volumes/MyBook/conus-hydrofabric/nwm_v3.0.7/conus/wrfinput_CONUS.nc')

all = c(soil, wrf)
all = all[[ grepl(paste(vars, collapse = "|"), names(all)) ]] 

names(all) = gsub("_Time=1", "", names(all))

for(i in 1:nlyr(all)){
  message(i)
  writeRaster(all[[i]], glue('/Volumes/MyBook/conus-hydrofabric/nwm_v3.0.7/conus/tif/{names(all[[i]])}.tif'), gdal = c("TILED=YES", 
                                                                                                                       "COPY_SRC_OVERVIEWS=YES", 
                                                                                                                       "COMPRESS=DEFLATE"))
}


library(hydrofabric)

wrf  = correct_nwm_spatial('/Volumes/MyBook/conus-hydrofabric/nwm_v3.0.7/conus/wrfinput_CONUS.nc')

AOI = AOI::aoi_get(state = "FL")

crop(wrf[['ISLTYP_Time=1']], project(vect(AOI), crs(wrf))) %>% 
  classify(data.frame(from = 16, to = 100)) %>% 
  plot()

wrf[['ISLTYP_Time=1']] %>% 
  classify(data.frame(from = 16, to = 100)) %>% 
  plot()


s = rast('/Volumes/MyBook/conus-hydrofabric/nwm_v3.0.7/conus/tif/ISLTYP.tif') 

tt = terra::classify(s, select(soil_class, sub, generic))
plot(tt)



crop(s, project(vect(AOI), crs(lu))) %>% 
  classify(data.frame(from = 16, to = 100)) %>% 
  plot()
# Soil Generic

soil_class = tribble(
~type, ~sub,  ~generic,
"land", 1,  1,
"land", 2,  1,
"land", 3,  1,
"land", 4,  1,
"land", 5,  1,
"land", 6,  1,
"land", 7,  1,
"land", 8,  1,
"land", 9,  1,
"land", 10,  1,
"land", 11,  1,
"land", 12,  1,
"land", 13,  1,
"land", 15,  1,
"land", 17,  1,
"land", 18,  1,
"land", 19,  1,

"water", 14,  14,

"ice", 16,  16
)

# Land cover Generic

tribble(
  ~type, ~sub,  ~generic,
  "urban", 1,  1,
  
  "agricultural", 2,  2,
  "agricultural", 3,  2,
  "agricultural", 4,  2,
  "agricultural", 5,  2,
  "agricultural", 6,  2,
  
  "grassland", 7,  7,
  "grassland", 8,  7,
  "grassland", 9,  7,
  "grassland", 10, 7,
  
  "forest", 11,  11,
  "forest", 12,  11,
  "forest", 13,  11,
  "forest", 14,  11,
  "forest", 15,  11,
  
  "water", 16,  16,
  
  "wetland", 17,  17,
  "wetland", 18,  17,
  
  "barren", 19,  19,
  "barren", 25,  19,
  "barren", 26,  19,
  "barren", 27,  19,
  
  "tundra", 20,  20,
  "tundra", 21,  21,
  "tundra", 22,  22,
  "tundra", 23,  23,
  
  "ice", 24,  24
)
