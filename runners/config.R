## NextGen Hydrofabric Config File ##
## 
## This file is used to populate scripts that are part of the Makefile pipeline
library(tidyr)
library(RNetCDF)
library(dataRetrieval)
library(data.table)
#library(AOI)

# ---- Set up software ---- #
options(scipen = 999)
dev_mode = TRUE
FIX = FALSE
ref_date = "04052024"

if(dev_mode){
  message("DEVMODE: ON")
  message("FIX ONLY: ", ifelse(FIX, "ON", "OFF"))
  devtools::load_all()
  # Either install and update your own dev versions or set dev_mode to FALSE
  devtools::load_all(glue::glue('{dirname(getwd())}/ngen.hydrofab'))
  devtools::load_all(glue::glue('{dirname(getwd())}/hfsubsetR'))
  # devtools::load_all(glue::glue('{dirname(getwd())}/zonal'))
  # devtools::load_all(glue('{dirname(getwd())}/zonal'))
} else {
  message("DEVMODE: OFF")
  message("FIX ONLY: ", ifelse(FIX, "ON", "OFF"))
  library(hydrofabric)
}

sf::sf_use_s2(FALSE)

# Runner Parameters -------------------------------------------------------

# 1. Runners will build and access data from a defined directory. Put that here:
dir      <- "/Volumes/MyBook"
# 2. Define the HF version you are creating:
version  <- "2.2"

# POI  config: 
schema <- c("poi_id", 
            "hl_source", "hl_reference", "hl_link", "hl_position",
            "X", "Y",
            "hf_id", "hf_source")

base_ref <- "/Volumes/MyBook/TNC/v2.2/reference"

# Desired community POI types
community_hl_types    <-  c('HUC12', 'Gages', 'TE', 'NID', "resops", "hilarri", "WBOut", "WBIn", "AR", "Term")

# 3. Define your processing VPUS
vpus <- vpu_boundaries$VPUID[1:21]

# Refactor and Aggregate parameters
ideal_size_sqkm <- 10
min_length_km   <- 1 
min_area_sqkm   <- 3

# NextGen Prefixes
nexus_prefix          <- "nex-"
terminal_nexus_prefix <- "tnx-"
coastal_nexus_prefix  <- "cnx-"
internal_nexus_prefix <- "inx-"
catchment_prefix      <- "cat-"
waterbody_prefix      <- "wb-"

# Directory for NWM data
#  Download needed files for add_cfe_noahowp_attributes(...)
#  From: https://www.nco.ncep.noaa.gov/pmb/codes/nwprod/{latest_version}/parm/domain/
#   soilproperties_CONUS_FullRouting.nc
#   wrfinput_CONUS.nc
#   GWBUCKPARM_CONUS_FullRouting.nc
nwm_dir       <- "/Volumes/Transcend/nwmCONUS-v216"
# Someday this might be available here (https://www.nco.ncep.noaa.gov/pmb/codes/nwprod/nwm.v3.0.7/parm/domain/), 
# but for now I am not allowed to share it per OWP. So hard path it is!
lake_path     <- glue('{nwm_dir}/LAKEPARM_CONUS.nc')

# Directory Structure -----------------------------------------------------

# # <<< BASE >>>> 
# ├── conus_hl.gpkg
# ├── reference
# │   ├── reference_01.gpkg
# │   ├── ...
# ├── refactored
# │   ├── refactor_01.gpkg
# │   ├── ....
# ├── hydrolocations
# │   ├── hl_01.gpkg
# │   ├── ...
# ├── corrected_refactor
# │   ├── refactor_01.gpkgh
# │   ├── ...
# ├── uniform
# │   ├── uniform_01.gpkg
# │   ├── ...
# ├── global_uniform
# │   ├── uniform_01.gpkg
# │   ├── ...
# ├── VERSION
# │   ├── conus.gpkg
# │   ├── conus_net.parquet
# │   ├── model_attributes.parquet
# │   ├── camels
# │   │   ├── Gages-01169000.gpkg
# │   │   ├── ...
# │   ├── fgb
# │   │   ├── divides.fgb
# │   │   ├── ...
# │   ├── gpkg
# │   │   ├── nextgen_01.gpkg
# │   │   ├── ...
# │   ├── model_attributes
# │   │   ├── nextgen_01.parquet
# │   │   ├── ...

base     <- glue("{dir}/conus-hydrofabric/v{version}")
dir.create(base, recursive = TRUE, showWarnings = FALSE)

base_reference <- glue("{base}/reference")
dir.create(base_reference, showWarnings = FALSE)

base_refactored = glue("{base}/refactored_{ref_date}")
dir.create(base_refactored, showWarnings = FALSE)

#base_corrected = glue("{base}/corrected_refactor_{ref_date}")
base_corrected = glue("{base}/corrected_refactor")
dir.create(base_corrected, showWarnings = FALSE)

base_hl = glue("{base}/hydrolocations")
dir.create(base_hl, showWarnings = FALSE)

base_fgb   = glue('{base}/{version}/fgb')
dir.create(base_fgb, recursive = TRUE, showWarnings = FALSE)

base_gpkg = glue('{base}/{version}/gpkg')
dir.create(base_gpkg, recursive = TRUE, showWarnings = FALSE)

base_cfe = glue('{base}/cfe')
dir.create(base_cfe, showWarnings = FALSE)

base_atts = glue('{base}/{version}/model_attributes')
dir.create(base_atts, showWarnings = FALSE)

base_uniform = glue('{base}/uniform')
dir.create(base_uniform, showWarnings = FALSE)

base_global_uniform = glue('{base}/global_uniform')
dir.create(base_global_uniform, showWarnings = FALSE)

base_camels = glue('{base}/{version}/camels')
dir.create(base_camels, showWarnings = FALSE)

# Input Datasets ----------------------------------------------------------------

# Gage Summary data is shared by the coastal team
coastal_gages = glue('{base}/GAGE_SUMMARY.csv')

# Downloaded from FIM ESIP bucket
# Source: s3://noaa-nws-owp-fim/hand_fim/inputs/ahps_sites/nws_lid.gpkg
# I will ask if I can disseminate this in Lynker Spatial
fim_ahps      = glue('{base}/nws_lid.gpkg')

community_hl_types    <-  c('HUC12', 'Gages', 'TE', 'NID', "resops", "hilarri", "WBOut", "WBIn", "AR")

full_hl    <- glue("{base}/conus_hl.gpkg")

poi_fgb     <- gsub("hl.gpkg", "poi.fgb", full_hl) 

hl_parquet <- gsub(".gpkg", "", full_hl)

# Converted from latest RL file 
# Source: https://www.nco.ncep.noaa.gov/pmb/codes/nwprod/{latest_version}/parm/domain/
# Converted version: s3://lynker-spatial/tabular-resources/RouteLink_nwm_v2_2_3.parquet
rl_file       = glue("{base}/RouteLink_nwm_v2_2_3.parquet")

# Downloaded from USGS
# https://www.sciencebase.gov/catalog/item/57eaa10fe4b09082500db04e
# Stored here: s3://lynker-spatial/tabular-resources/huc12_nhdplusv2_cw.parquet
huc12_cw      = glue("{base}/huc12_nhdplusv2_cw.parquet")

gs_file       = 'https://code.usgs.gov/wma/nhgf/reference-hydrofabric/-/raw/04cd22f6b5f3f53d10c0b83b85a21d2387dfb6aa/workspace/cache/rpu_vpu_out.csv'

ms_lu         = 'https://github.com/internetofwater/ref_rivers/releases/download/v2.1/mainstem_lookup.csv'

camels        <- 'https://lynker-spatial.s3.amazonaws.com/tabular-resources/camels_compiled.csv'

hydroatlas    <- 's3://lynker-spatial/hydroATLAS/hydroatlas_vars.parquet'

pipeline = data.frame(
  vpus            = vpus,
  uniform         = glue("{base_uniform}/uniform_{vpus}.gpkg"),
  uniform_global  = glue("{base_global_uniform}/uniform_{vpus}.gpkg"),
  nextgen         = glue('{base_gpkg}/nextgen_{vpus}.gpkg'),
  cfe             = glue("{base_cfe}/cfe_noahowp_{vpus}.parquet"),
  atts            = glue("{base_atts}/nextgen_{vpus}.parquet")#,
  # refactored_gpkg = sapply(1:length(vpus), FUN = \(x){ get_hydrofabric(vpus[x], type = "refactor",  dir = base_refactored) }),
  # reference_gpkg  = sapply(1:length(vpus), FUN = \(x){ get_hydrofabric(vpus[x], type = "reference", dir = base_reference) })
)


# Data Definition ---------------------------------------------------------

# Define:
## Target SRS (t_srs)
## Target extent (te)
## Target grid resolution (ts)
## Resampling method (r)

target = "-t_srs '+proj=lcc +lat_0=40.000008 +lon_0=-97 +lat_1=30 +lat_2=60 +x_0=0 +y_0=0 +R=6370000 +units=m +no_defs' -te -2303999.62876143 -1920000.70008381 2304000.37123857 1919999.29991619 -ts 18432 15364 -r average"

# "Raw files" are used for processing here. They are not needed in 99% of cases!!! 
# We provide all the data here: s3://lynker-spatial/gridded-resources/250m_grids/

# Sourced from: 
# '/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/TIFF/USGS_Seamless_DEM_1.vrt'
# If you feel you need to run this yourself be a good citizen and download the vrt and requisite data to avoid
# exploiting USGS resources. That said, you shouldn't need to s
raw_dem = '/Volumes/MyBook/3DEP/usgs_1.vrt'
# Sourced from: https://www.mrlc.gov/
raw_imp_file <-'/Volumes/MyBook/nlcd/nlcd_2021_impervious_l48_20230630/nlcd_2021_impervious_l48_20230630.img'
# Source Data: https://zenodo.org/records/4460354
raw_twi_file  <- '/Volumes/MyBook/CONUS_TWI_epsg5072_30m_unmasked.tif'



# ------ Processing ----- #
dem_file      = glue("{dir}/250/usgs_1_250mm.tif")

if(!file.exists(dem_file)){
  glue("gdalwarp {target} {raw_dem} {dem_file}")
}

slope_file    <- glue("{dir}/250m_grids/usgs_250m_slope.tif")

if(!file.exists(slope_file)){
  glue("gdaldem slope {dem_file} {slope_file}")
}

aspect_file   <- glue("{dir}/250m_grids/usgs_250m_aspect.tif")

if(!file.exists(slope_file)){
  glue("gdaldem aspect {dem_file} {aspect_file}")
}

imp_file      <- glue("{dir}/250m_grids/imperv_250m.tif")

if(!file.exists(imp_file)){
  glue("gdalwarp {target} {raw_imp_file} {imp_file}")
}

twi_file      <- glue("{dir}/250m_grids/twi_250m.tif")

if(!file.exists(twi_file)){
  glue("gdalwarp {target} {raw_twi_file} {twi_file}")
}


# Built Datasets ----------------------------------------------------------------
#These are files that will be built in the base directory as the workflow progresses:

conus_net  <- glue("{base}/{version}/conus_net.parquet")

conus_gpkg <- glue("{base}/{version}/conus.gpkg")

atts       <- glue("{base}/{version}/model_attributes.parquet")

# Refactor Corrections ----------------------------------------------------

update_topo = tribble(
  ~VPU, ~id,  ~toid,
  # https://code.usgs.gov/wma/nhgf/reference-hydrofabric/-/issues/141
  "10U", 10010387,  10010388,
  # https://code.usgs.gov/wma/nhgf/reference-hydrofabric/-/issues/143
  "10L", 10023804,  10023793,
  # https://github.com/NOAA-OWP/hydrofabric/issues/30
  "01", 10040692,  10018009,
  "01", 10018009,  10018008,
  # https://github.com/NOAA-OWP/hydrofabric/issues/31
  "01", 10012415,  10012226,
  "01", 10012226,  10012125,
  "01", 10040737,  10002480,
  "01", 10003900,  10004000,
  "01", 10003943,  10003900,
  "10L", 10000635, 10001228,
  "10U", 10035086, 10035063,
  "10U", 10125872, 10125873,
  "10U", 10179149, 10179145,
  "10U", 10179831, 10179145
)

remove_fp_ids = tribble(
  ~VPU, ~id,
  # https://code.usgs.gov/wma/nhgf/reference-hydrofabric/-/issues/141
  "10U", '10010386',
  #https://github.com/NOAA-OWP/hydrofabric/issues/29
  "07", '10049093',
  # https://github.com/NOAA-OWP/hydrofabric/issues/31
  "01", '10016831',
  "01", '10003899',
  "10L", '10013259,10013257,10013260,10013261,10013262,10013263,10013270,10013269,10013254,10013256,10013253,10013255,10013329,10013129,10013258',
  "10L", "10000441",
  "10U", "10035132",
  "10U", "10152753,10152754,10126589,10125886,10125883,10125884,10125885,10126574,10125887,10126657,10126458",
  "10U", "10179150,10179151,10179147,10179146,10179153,10179148,10179152",
  "05", "10111481"
)

remove_divide_ids = tribble(
  ~VPU, ~divide_id,
  # https://github.com/NOAA-OWP/hydrofabric/issues/31
  "05", "10111481"
)

redigitize_flowpaths = tribble(
  ~VPU, ~id,  
  # https://github.com/NOAA-OWP/hydrofabric/issues/30
  "01", 10018009
)


snap_nodes =  tribble(
  ~VPU, ~end_node_of, ~start_node_of,
  "01", 10003900, 10004000
)
  

divides_to_merge = tribble(
    ~VPU, ~to_merge, ~id, ~toid,
    # https://github.com/NOAA-OWP/hydrofabric/issues/29
    "07", "10049093,10047655", '10047655', '10049122',
    #https://github.com/NOAA-OWP/hydrofabric/issues/31
    "01", "10040781,10040780,10012125", "10012125", '10012500',
    "01", "10016831,10016708", "10016708", '10016707',
    "01", "10040738,10040737", "10040737", '10002480',
    "10L", "10013329,10013129", "10013251", "10013134",
    "10L", "10013254,10013255,10013256", '10013252', '10013251',
    "10L", "10013259,10013260,10013258,10013261,10013269,10013262,10013263,10013270",'10013333', '100132552',
    "10U", '10179150,10179151', '10179145', '10179945',
    "10U", '10179945,10179153,10179152', '10179945', '10179156'
)

flowlines_to_merge = 
  tribble(
    ~VPU, ~to_merge, ~id, ~toid,
    #https://github.com/NOAA-OWP/hydrofabric/issues/31
    "01",  "10040781,10040780,10012125", "10012125", '10012500',
    "01",  "10040738,10040737", "10040737", '10002480',
    "16",  "0046993,10047551", "10047551", "NA",
    "10U", "10179144,10179145", '10179145', '10179945',
    "10U", "10179945,10179154,10179155", '10179945', '10179156',
    "10U", '10093159,10093160,10093161,10093162,10095778','10095778','NA'
)

reassign_mainstem =
  tribble(
    ~VPU,  ~id, ~mainstem,
    # https://github.com/NOAA-OWP/hydrofabric/issues/30
    "01", 10018008, 2111164,
    #https://github.com/NOAA-OWP/hydrofabric/issues/31
    "10L", 10001228, 1135798
)

VPUS = c(
    update_topo$VPU,
    remove_fp_ids$VPU,
    redigitize_flowpaths$VPU,
    divides_to_merge$VPU,
    flowlines_to_merge$VPU,
    reassign_mainstem$VPU,
    snap_nodes$VPU
  )

pipeline$corrected_refactor = ifelse(vpus %in% unique(VPUS), glue("{base_corrected}/corrected_refactor_{vpus}.gpkg"), NA)
