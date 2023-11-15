## TASK 2: Assign Globally Unique Identifiers ----

pacman::p_load(hydrofabric)
devtools::load_all()

vpus = vpu_boundaries$VPUID[1:21]

base = '/Volumes/MyBook/nextgen'

## TASK 1: build out uniform catchment distribution ----
process = data.frame(
  vpus = vpus,
  outfiles = glue("{base}/uniform/uniform_{vpus}.gpkg"),
  global = glue("{base}/global_uniform/uniform_{vpus}.gpkg")
)

unlink(process$global)

gs_file = 'https://code.usgs.gov/wma/nhgf/reference-hydrofabric/-/raw/04cd22f6b5f3f53d10c0b83b85a21d2387dfb6aa/workspace/cache/rpu_vpu_out.csv'

modifications = read.csv(gs_file) %>% 
  filter(VPUID != toVPUID) %>% 
  rename(from = COMID, to = toCOMID)

meta = assign_global_identifiers(gpkgs = process$outfiles, 
                                 outfiles = process$global,
                                 modifications = modifications)

## TASK 3: Assign Globally Unique Identifiers
## 
for(i in 1:nrow(process)){
  try(append_style(process$global[i], layer_names = c("flowpaths", "divides", "hydrolocations")), silent = TRUE)
}