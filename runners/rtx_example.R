# Install -----------------------------------------------------------------
# install.packages("remotes") 
# remotes::install_github("NOAA-OWP/hydrofabric")


# Setup -------------------------------------------------------------------

library(hydrofabric)

make_map = function(file, pois) {
  hf = read_hydrofabric(file)
  mapview::mapview(hf$catchments) + hf$flowpaths + pois
}

## NOTE: What is glue?
x <- "This is hydrofabric"
y <- "version 2.2"

glue("{x} {y}")

## NOTE: Finding help/docs
?read_hydrofabric

### ---- Data Stores  ---- ###
using_local_example <- "/Volumes/MyBook/conus-hydrofabric"
using_s3_example    <- "s3://lynker-spatial/hydrofabric"

# Geoparquet store
fs::dir_tree(glue("{using_local_example}/v2.2/reference/conus_divides"))
open_dataset(glue("{using_local_example}/v2.2/reference/conus_divides"))

# Parquet store
fs::dir_tree(glue("{using_local_example}/v2.2/reference/conus_network"))
open_dataset(glue("{using_local_example}/v2.2/reference/conus_network"))

# Replicated in s3
system(glue('aws s3 ls {using_s3_example}/v2.2/reference/conus_divides/'))
open_dataset(glue('{using_s3_example}/v2.2/reference/conus_divides/'))

# Getting data
### Defered Evaluation
open_dataset(glue('{using_s3_example}/v2.2/reference/conus_network/')) %>% 
  filter(id == 101) %>% 
  select(id, toid)

### Extract from s3
open_dataset(glue('{using_s3_example}/v2.2/reference/conus_network/')) %>% 
  filter(id == 101) %>% 
  select(id, toid) %>% 
  collect()

### Extract from local
open_dataset(glue('{using_local_example}/v2.2/reference/conus_network/')) %>% 
  filter(id == 101) %>% 
  select(id, toid) %>% 
  collect()


### ---- Sample outfiles for today ---- ###
reference_file  <- "vignettes/tutorial/poudre.gpkg"
refactored_file <- "vignettes/tutorial/refactored.gpkg"
aggregated_file <- "vignettes/tutorial/aggregated.gpkg"
nextgen_file    <- "vignettes/tutorial/nextgen.gpkg"




# Get Reference Fabric ----------------------------------------------------
## ---  Define starting feature by source and ID
## https://waterdata.usgs.gov/monitoring-location/06752260
## https://reference.geoconnex.us/collections/gages/items?provider_id=06752260

(gage <- list(featureSource = "nwis", featureID = "06752260"))

# Use get_subset to build a reference subset
get_subset(
  nldi_feature = gage,
  source  = using_local_example,
  type = "reference",
  hf_version = "2.2",
  outfile = reference_file,
  overwrite = TRUE
)


## ---  Get some Points of Interest
hf = read_hydrofabric(reference_file)

pois = open_dataset(glue("{using_local_example}/v2.2/conus_hl")) %>%
  filter(hl_source == 'GFv20') %>%
  collect() %>%
  st_as_sf(coords = c("X", "Y"), crs = 5070) %>%
  st_filter(hf$catchments)

make_map(reference_file, pois)

# Build a Refactored Fabric -----------------------------------------------

refactored = refactor(
  reference_file,
  split_flines_meters = 10000,
  collapse_flines_meters = 1000,
  collapse_flines_main_meters = 1000,
  pois = pois,
  fac = '/vsis3/lynker-spatial/gridded-resources/fac.vrt',
  fdr = '/vsis3/lynker-spatial/gridded-resources/fdr.vrt',
  outfile = refactored_file
)

make_map(refactored_file, pois)



# Build an Aggregated Network ---------------------------------------------

hydrolocations = read_sf(refactored_file, 'lookup_table') %>%
  inner_join(pois, by = c("NHDPlusV2_COMID" = "hf_id")) %>%
  select(poi_id, NHDPlusV2_COMID, id = reconciled_ID) %>%
  distinct()

aggregate_to_distribution(
  gpkg = refactored_file,
  hydrolocations = hydrolocations,
  ideal_size_sqkm = 10,
  min_length_km = 1,
  min_area_sqkm = 3,
  outfile = aggregated_file,
  overwrite = TRUE
)

make_map(aggregated_file, pois)


# Generate a NextGen Network ----------------------------------------------
unlink(nextgen_file)

apply_nexus_topology(aggregated_file, export_gpkg = nextgen_file)

hf = read_hydrofabric(nextgen_file)

mapview::mapview(hf$catchments) + hf$flowpaths + read_sf(nextgen_file, "nexus")

# Populate some data for NOAHOWP / CFE------------------------------------------

div = read_sf(nextgen_file, "divides")

noah_owp_vars = c(
  "bexp",
  "cwpvt",
  "dksat",
  "ISLTYP",
  "IVGTYP",
  "mfsno",
  "mp",
  "psisat",
  "refkdt",
  "slope",
  "smcmax",
  "smcwlt",
  "vcmx25"
)

urls = glue("/vsis3/lynker-spatial/gridded-resources/nwm/conus/{noah_owp_vars}.tif")

# Full Domain
rast(urls)

# Extracted Domain
grids = dap(URL = urls, AOI = div)

# Summarize data
noah_owp_data = execute_zonal(grids, div, ID = "divide_id", 
                              fun = "mean",
                              join = FALSE)
write_parquet(noah_owp_data,
              'vignettes/tutorial/nextgen.parquet')
# View
plot(noah_owp_data['mean.smcwlt_soil_layers_stag.3'])

# Maybe find some other data
## - Can use any catalog or shortcut function
gm = getGridMET(
  AOI = div,
  varname = "tmmx",
  startDate = "2024-03-01",
  endDate = "2024-03-15"
)

summary_data = execute_zonal(gm, div, ID = "divide_id")

## ERROR bad request!
getNLDAS(
  AOI = div,
  varname = "tmmx",
  startDate = "2024-03-01"
)

# Fix and get NLDAS data 
(nldas = getNLDAS(
  AOI = div,
  varname = "tmp2m",
  startDate = "2024-03-01"
))

plot(nldas$tmp2m)


# Weight Grids ------------------------------------------------------------

fg = terra::rast('/Users/mjohnson/Downloads/nwm.t00z.medium_range.forcing.f001.conus.nc')[[4]]

w = zonal::weight_grid(fg, div, ID = "divide_id")

w$grid_id = "medium_range.forcing.conus"

write_parquet(w, 'vignettes/tutorial/nextgen_grid_weights.parquet')

# Populate Flowpath Attributes --------------------------------------------
# https://www.nco.ncep.noaa.gov/pmb/codes/nwprod/
latest_nwm_version()

get_routelink_path()

add_flowpath_attributes(nextgen_file)

# Schema, showing as_sqlite utility
as_sqlite(nextgen_file)

# Data
as_sqlite(nextgen_file, 'flowpath_attributes')

# Future!!
open_dataset(glue('{using_s3_example}/v2.2/reference/routelink_ls/'))






# Make it pretty :)  ------------------------------------------------------

append_style(nextgen_file, layer_names = c("divides", "flowpaths", "nexus"))
