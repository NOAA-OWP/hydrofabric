outdir = '/Volumes/MyBook/conus-hydrofabric/v20.1/forcing_weights'
version = '/Volumes/MyBook/conus-hydrofabric/v20.1/gpkg'
df = data.frame(f = list.files(version, pattern = ".gpkg$", full.names = TRUE))

dir.create(outdir)
     
df$vpus = gsub(".gpkg", "", gsub("nextgen_", "", basename(df$f)))

df$outfile = glue('/Volumes/MyBook/conus-hydrofabric/v20.1/forcing_weights/forcing_weights_{df$vpus}.parquet')

fg = terra::rast('/Users/mjohnson/Downloads/nwm.t00z.medium_range.forcing.f001.conus.nc')[[4]]

for(i in 1:nrow(df)){
  message(i)
  w   = weight_grid(fg, 
                    read_sf(df$f[i], "divides"), 
                    ID = "divide_id")
  
  w$grid_id = "medium_range.forcing.conus"
  write_parquet(w, df$outfile[i])
}


all = list()

for(j in 1:nrow(df)){
  all[[j]] = arrow::read_parquet(df$outfile[j])
}

all2 = dplyr::bind_rows(all)

length(unique(all2$divide_id))
all2[1,]

arrow::write_parquet(all2, '/Volumes/MyBook/conus-hydrofabric/v20.1/forcing_weights.parquet')

data.frame(
  grid_id = "medium_range.forcing.conus",
  schema = "climateR",
  X1        = terra::xmin(fg) ,
  Xn           = terra::xmax(fg) ,
  Y1           = terra::ymin(fg) ,
  Yn          = terra::ymax(fg) ,
  resX         = terra::res(fg)[2],
  resY       = terra::res(fg)[2],
  ncols      = terra::ncol(fg)  ,
  nrows = terra::nrow(fg),
  crs = terra::crs(fg, proj = TRUE)
) %>%
  jsonlite::write_json("/Volumes/MyBook/conus-hydrofabric/v20.1/forcing_grids.json")


