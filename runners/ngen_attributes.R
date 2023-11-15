library(terra)

r = rast('/Users/mjohnson/Downloads/ngen_gridded_data.nc')

template = list(ext = ext(-2303999.62876143, 2304000.37123857, 
                          -1920000.70008381, 1919999.29991619), crs = "PROJCS[\"Sphere_Lambert_Conformal_Conic\",GEOGCS[\"GCS_Sphere\",DATUM[\"D_Sphere\",SPHEROID[\"Sphere\",6370000.0,0.0]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Lambert_Conformal_Conic\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",-97.0],PARAMETER[\"standard_parallel_1\",30.0],PARAMETER[\"standard_parallel_2\",60.0],PARAMETER[\"latitude_of_origin\",40.000008],UNIT[\"Meter\",1.0]];-35691800 -29075200 126180232.640845;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision")

n =names(r)
layers = sapply(strsplit(n, "_"),"[[",1)
n2 = gsub('_Time=1', "", n)

for(i in 1:length(layers)){
  
  d = r[[grepl(layers[i], n)]]
  names(d) = gsub('_Time=1', "", names(d))
  o = glue::glue("nwm/{layers[i]}.tif")
  terra::ext(d) <- template$ext
  terra::crs(d) <- template$crs
  
  terra::writeRaster(d, o, filetype = "COG", overwrite = T)
  message(i, " of ", length(layers))
}