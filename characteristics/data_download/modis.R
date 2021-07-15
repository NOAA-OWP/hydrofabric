# Downloaded and corrected from here:
# http://book.ecosens.org/modis-sinusoidal-grid-download/
# modis_grid = file.path('/Volumes/Transcend/land-cover-tests/MODIS/modis_grid/modis_sinusoidal_grid_world.shp') %>%
#   rgdal::readOGR(verbose = FALSE) %>%
#   st_as_sf() %>%
#   st_set_crs(modis_proj)
# 
# write_sf(modis_grid, "characteristics/data/modis_grid.gpkg")


# MODIS -------------------------------------------------------------------
library(sf)

# LAI:
# NDVI:

base.url = 'https://e4ftl01.cr.usgs.gov/MOLT/'
type = "MOD15A2H.006"
AOI = AOI::aoi_get(state = "conus")
startDate = "2000-02-01"
endDate = "2020-12-31"
dir = '/Volumes/Transcend/ngen/climate/'

dates = seq.Date(as.Date(startDate), as.Date(endDate), by = "m")

modis_proj = "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext"



reg_bbox_modis_proj   <-  st_transform(AOI, st_crs(modis_grid))
ints  <-  st_filter(modis_grid, reg_bbox_modis_proj)
# What MODIS tiles intersect with the REGION DOMAIN
tiles <-  paste0("h", sprintf("%02d", as.numeric(ints$h)), "v", sprintf("%02d", as.numeric(ints$v)))

dates = paste0(base.url, type) %>% 
  read_html() %>%
  html_nodes("a") %>%
  html_attr('href')

dates = lubridate::ymd(dates) 
dates = dates[dplyr::between(dates, as.Date(startDate), as.Date(endDate))] %>% na.omit()

fs::dir_create(file.path('/Volumes/Transcend/ngen/climate/','MODIS/', type,  dates ))


for(i in 1:length(dates)){
  path = paste0('https://e4ftl01.cr.usgs.gov/MOLT/', type, '/', gsub("-", ".", dates[i]))
  
  files = path %>% 
    read_html() %>%
    html_nodes("a") %>%
    html_attr('href')
  
  f1 = grep(".hdf$", files, value = TRUE)
  f2 = grep(paste(tiles, collapse = "|"), f1, value = TRUE)
  
  urls = file.path(path, f2)
  
  tmp = file.path('/Volumes/Transcend/ngen/climate','MODIS', type,  dates[i], f2 )
  
  message("Downloading: ", path)
  lapply(1:length(urls), function(x){
    if(!file.exists(tmp[x])){
      httr::GET(urls[x],
                write_disk(tmp[x], overwrite = TRUE),
                progress(),
                config(netrc = TRUE, netrc_file = getNetrcPath()),
                set_cookies("LC" = "cookies"))
    }
  })
}



gdalinfo_raw <- sf::gdal_utils(util = "info", tmp[1], quiet  = F)

rr = strsplit(gdalinfo_raw, "\\n")[[1]]

subdataset_rawnames <-  rr[grepl(glob2rx("*SUBDATASET*NAME*"), rr)]
subdataset_names <- sapply(X = seq(length(subdataset_rawnames)),
                           FUN = function(X) {
                             split1 <- strsplit(subdataset_rawnames[X], "=")
                             return(gsub("\"", "", split1[[1]][2]))
                           })

dd = gsub(tmp[1],"[file]",subdataset_names)

stringr::str_replace(dd, "\\[file\\]", tmp)

gsub("\\[.*?\\]", tmp, rep(dd[1], length(tmp)))

files = paste0("HDF4_EOS:EOS_GRID:",tmp,":MOD_Grid_monthly_1km_VI:1 km monthly NDVI")

unlink('/Users/mikejohnson/Desktop/2001-01-01/test2.tif')
sf::gdal_utils(util = "warp",
               source = files,
               destination  = '/Users/mikejohnson/Desktop/2001-01-01/test2.tif',
               options = c("-of", 'GTiff')
)

}
