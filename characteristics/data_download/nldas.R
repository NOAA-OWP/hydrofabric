ts   = "M"
type = "FORA"
startDate = "2000-01-01"
endDate = "2020-12-31"
var  = "APCP"
ext  = "tif"
dir  = '/Volumes/Transcend/ngen/climate/NLDAS/'

get_nldas = function(
  type = "FORA",
  ts   = "M",
  startDate = "2000-01-01",
  endDate = "2020-12-31",
  var  = "APCP",
  ext  = "tif",
  dir  = '/Volumes/Transcend/ngen/climate/NLDAS/'){
  
  date = seq.Date(as.Date(startDate), as.Date(endDate), by = "m")
  
  bbox = "25%2C-125%2C53%2C-67"
  type2 = paste0('NLDAS_', type, "0125_", ts)
  date2 = gsub("-",  "", date)
  year = lubridate::year(date)
  month = sprintf("%02s", lubridate::month(date))
  ext = ifelse(ext == "nc", ".nc4", ".tif")
  
  url = paste0('https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=',
               '%2Fdata',
               '%2FNLDAS',
               '%2F',type2,'.002',
               '%2F', year,
               '%2F',type2,'.A',year,month,'.002','.grb',
               '&FORMAT=dGlmLw',
               '&BBOX=', bbox,
               '&LABEL=',type2,'.A',year,month,'.002','.grb.SUB', ext,
               '&SHORTNAME=',type2,
               '&SERVICE=L34RS_LDAS',
               '&VERSION=1.02',
               '&DATASET_VERSION=002',
               '&VARIABLES=', var)
  
  tmp = file.path(dir, var, paste0(type2, ".", date, ext))
  fs::dir_create(dirname(tmp))
  
  lapply(1:length(url), function(x){
    if(!file.exists(tmp[x])){
      httr::GET(url[x],
                write_disk(tmp[x], overwrite = TRUE),
                progress(),
                config(netrc = TRUE, netrc_file = getNetrcPath()),
                set_cookies("LC" = "cookies"))
    }
  })
  
}


get_nldas(
  type = "FORA",
  ts   = "M",
  startDate = "2000-01-01",
  endDate = "2020-12-31",
  var  = "APCP",
  ext  = "tif",
  dir  = '/Volumes/Transcend/ngen/climate/NLDAS/')

get_nldas(
  type = "FORA",
  ts   = "M",
  startDate = "2000-01-01",
  endDate = "2020-12-31",
  var  = "PEVAP",
  ext  = "tif",
  dir  = '/Volumes/Transcend/ngen/climate/NLDAS/')

get_nldas(
  type = "VIC",
  ts   = "M",
  startDate = "2000-01-01",
  endDate = "2020-12-31",
  var  = "SNOWC",
  ext  = "tif",
  dir  = '/Volumes/Transcend/ngen/climate/NLDAS/')
