pacman::p_load(hydrofabric, glue) 
devtools::load_all()

vpus = vpu_boundaries$VPUID[1:21]

base = '/Volumes/MyBook/nextgen'

ref = glue("{base}/reference")
refac = glue("{base}/refactored")

ow = TRUE

dir.create(ref, showWarnings = FALSE); dir.create(refac, showWarnings = FALSE)

# Refreshed last on: 
  # 06-27-2023 
  # ...
  # ...
  # ...

for(i in 15:length(vpus)){
  
  source('secret/sb_tools.R')
  
  get_hydrofabric(VPU = vpus[i], 
                  type = "reference", 
                  dir =  ref, 
                  overwrite = ow)
  
  get_hydrofabric(VPU = vpus[i], 
                  type = "refactor",  
                  dir = refac, 
                  overwrite = ow)
  
  message(vpus[i])
}

library(dplyr)
sf::sf_use_s2(FALSE)

tmp = tempfile(fileext = ".geojson")
httr::GET("https://earth-info.nga.mil/php/download.php?file=hydrobasins_level2", httr::write_disk(tmp))

xx = sf::read_sf('/Users/mjohnson/Downloads/hydrobasins_level2.geojson')

xx2 = st_make_valid(ms_simplify(xx))


