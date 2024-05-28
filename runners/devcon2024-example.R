
# Load Needed Libaries ----------------------------------------------------

library(hydrofabric)
library(powerjoin)

# Helper Function ---------------------------------------------------------

make_map = function(file, pois) {
  hf = read_hydrofabric(file)
  mapview::mapview(hf$catchments) + hf$flowpaths + pois
}

### ---- Sample out files and source for today ---- ###
fs::dir_create("tutorial")

# Setup Needed File Structure ---------------------------------------------

source    <- '/Users/mjohnson/hydrofabric/'

reference_file  <- "tutorial/poudre.gpkg"
refactored_file <- "tutorial/refactored.gpkg"
aggregated_file <- "tutorial/aggregated.gpkg"

nextgen_file       <- "tutorial/poudre_ng.gpkg"
model_atts_file    <- "tutorial/poudre_ng_attributes.parquet"
model_weights_file <- "tutorial/poudre_ng_weights.parquet"

# Extract Subset of Reference Fabric --------------------------------------

get_subset(
  hl_uri = "Gages-06752260",
  source  = using_local_example,
  type = "reference",
  hf_version = "2.2",
  lyrs = c("divides", "flowlines", "network"),
  outfile = reference_file,
  overwrite = TRUE
)

hf = read_hydrofabric(reference_file)

# Pull out community POIs for Subset --------------------------------------

pois = open_dataset(glue("{source}/v2.2/conus_hl")) %>%
  filter(hl_source == 'GFv20', 
         vpuid %in% unique(hf$flowpaths$vpuid),
         hf_id %in% hf$flowpaths$id) %>%
  collect() %>%
  st_as_sf(coords = c("X", "Y"), crs = 5070)

make_map(reference_file, pois)


# Build Refactored Fabric -------------------------------------------------

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


# Remap Hydrolocations ----------------------------------------------------

hydrolocations = read_sf(refactored_file, 'lookup_table') %>%
  inner_join(pois, by = c("NHDPlusV2_COMID" = "hf_id")) %>%
  select(poi_id, NHDPlusV2_COMID, id = reconciled_ID) %>%
  distinct()


# Aggregate Network -------------------------------------------------------

aggregate_to_distribution(
  gpkg = refactored_file,
  hydrolocations = hydrolocations,
  ideal_size_sqkm = 10,
  min_length_km = 1,
  min_area_sqkm = 3,
  outfile = aggregated_file,
  overwrite = TRUE )

make_map(aggregated_file, pois)


# Apply NextGen Topology and Structure ------------------------------------

unlink(nextgen_file)
apply_nexus_topology(aggregated_file, export_gpkg = nextgen_file)

make_map(nextgen_file, read_sf(nextgen_file, "nexus"))

# Enriching the Network

# Enrich your new network --------------------------------------------------

# LS Hooks for gridded data
vsi <- "/vsis3/lynker-spatial/gridded-resources"

# Divides for NextGen Fabric
div <- read_sf(nextgen_file, "divides")

# NOAH OWP Variables (Grid -> Polygon)
nom_vars <- c("bexp", "dksat", "psisat", "smcmax", "smcwlt")

# Read the first layer of each from the remote COG Store
# No need oto crop here as only the connection is being established
r = rast(glue("{vsi}/nwm/conus/{nom_vars}.tif"), 
         lyrs = seq(1,length(nom_vars)*4, by = 4))

# Get Mode Beta parameter
modes = execute_zonal(r[[1]], 
                    fun = mode,
                    div, ID = "divide_id", 
                    join = FALSE)  %>% 
    setNames(gsub("fun.", "", names(.)))

# Get Geometric Mean of Saturated soil hydraulic conductivity, 
# and matric potential
gm = execute_zonal(r[[2:3]], 
                    fun = geometric_mean,
                    div, ID = "divide_id", 
                    join = FALSE)  %>% 
    setNames(gsub("fun.", "", names(.)))

# Get Mean Saturated value of soil moisture and Wilting point soil moisture
m = execute_zonal(r[[4:5]], 
                    fun = "mean",
                    div, ID = "divide_id", 
                    join = FALSE)  %>% 
    setNames(gsub("mean.", "", names(.)))


# Merge all tables into one
d1 <- power_full_join(list(modes, gm, m), by = "divide_id")

# GWBUCKET Rescale Data  (Polygon --> Polygon) --------------------------------

crosswalk <- as_sqlite(nextgen_file, "network") |>
    select(hf_id, divide_id) |>
    collect()

d2 <- open_dataset(glue("{source}/v2.2/reference/conus_routelink")) |>
    select(hf_id , starts_with("gw_")) |>
    inner_join(mutate(crosswalk, hf_id = as.integer(hf_id)), by = "hf_id") |>
    group_by(divide_id) |>
    collect() |>
    summarize(
      gw_Coeff = round(weighted.mean(gw_Coeff, 
                                     w = gw_Area_sqkm, 
                                     na.rm = TRUE), 9),
      gw_Zmax_mm  = round(weighted.mean(gw_Zmax_mm,  
                                        w = gw_Area_sqkm, 
                                        na.rm = TRUE), 9),
      gw_Expon = mode(floor(gw_Expon))
    )

# Downscaleing Variables
# 
# 
d3 <- st_centroid(div) |>
  st_transform(4326) |>
  st_coordinates() |>
  data.frame() |>
  mutate(divide_id = div$divide_id)

dem_vars <- c("elev", "slope", "aspect")

r  <- rast(glue('{vsi}/250m_grids/usgs_250m_{dem_vars}.tif'))

d4 <- execute_zonal(r[[1:2]], 
                    div, ID = "divide_id", 
                    join = FALSE) |>
    setNames(c("divide_id", "elevation_mean", " slope"))

d5 <- execute_zonal(r[[3]], 
                     div, ID = "divide_id", fun = circular_mean, 
                     join = FALSE) |>
    setNames(c("divide_id", "aspect_c_mean"))


model_attributes <- power_full_join(list(d1, d2, d3, d4, d5), by = "divide_id")
  
# Build a Forcing Weights File --------------------------------------------

type = "medium_range.forcing"

w = weight_grid(rast(glue('{vsi}/{type}.tif')), div, ID = "divide_id") |> 
  mutate(grid_id = type)

# Write it all out --------------------------------------------------------
write_parquet(w, model_weights_file)
write_parquet(model_attributes, model_atts_file)





# Enhance the flow network ------------------------------------------------


crosswalk <- as_sqlite(nextgen_file, "network") |>
    select(hf_id, id, divide_id, hydroseq, poi_id) |>
    filter(!is.na(poi_id)) %>% 
    collect() %>% 
    slice_min(hydroseq)

open_dataset(glue("{source}/v2.2/reference/conus_routelink/")) |>
    select(hf_id, starts_with("ml_")) 

(cs <- open_dataset(glue("{source}/v2.2/reference/conus_routelink/")) |>
    select(hf_id, ml_y_bf_m, ml_tw_bf_m, ml_r) %>% 
    inner_join(mutate(crosswalk, hf_id = as.integer(hf_id)), by = "hf_id") |>
    collect() %>% 
    summarise(TW = mean(ml_tw_bf_m),
              r = mean(ml_r),
              Y = mean(ml_y_bf_m),
              poi_id = poi_id[1]))

bathy = AHGestimation::cross_section(r = cs$r, TW = cs$TW, Ymax = cs$Y) 

plot(bathy$x, bathy$Y, type = "l", 
     ylab = "Releative distance (m)", 
     xlab = "Depth (m)", 
     main = glue("Average XS at POI: {cs$poi_id}"))

library(plotly)

crosswalk <- as_sqlite(nextgen_file, "network") |>
    select(hf_id, id, toid, divide_id, hydroseq, poi_id) |>
    collect() %>% 
    slice_max(hydroseq)

cw = open_dataset(glue('{source}/v2.1.1/nextgen/conus_network')) %>% 
  semi_join(crosswalk, by = "hf_id") %>% 
  collect() 

message(sum(cw$lengthkm), " kilometers of river")

open_dataset(glue('{source}/v2.1.1/nextgen/conus_xs')) %>% 
  filter(vpuid %in% unique(cw$vpuid), hf_id %in% unique(cw$id)) %>% 
  group_by(hf_id, cs_id) %>% 
  collect() %>% 
  mutate(uid = cur_group_id()) %>% 
  plot_ly(x = ~X, y = ~Y, z = ~Z,  split = ~as.factor(uid),
          type = 'scatter3d', mode = 'markers+lines',
          line = list(width = 3), marker = list(size = 2)) %>% 
  layout(list(aspectmode='manual',
              aspectratio = list(x=100, y=100, z=1)),
              showlegend = FALSE)

add_flowpath_attributes(nextgen_file, source = source)

as_sqlite(nextgen_file, 'flowpath_attributes') %>% 
  collect() %>% 
  head()



# Add some symbology ------------------------------------------------------



append_style(nextgen_file, layer_names = c("divides", "flowpaths", "nexus"))

