## NextGen Hydrofabric Config File ##
options(scipen = 999)
dev_mode = TRUE
FIX = TRUE

if(dev_mode){
  message("DEVMODE: ON")
  message("FIX ONLY: ", ifelse(FIX, "ON", "OFF"))
  devtools::load_all()
  devtools::load_all(glue('/Users/mjohnson/github/ngen-hydrofab'))
  # devtools::load_all(glue('{dirname(getwd())}/hydrofab'))
  # devtools::load_all(glue('{dirname(getwd())}/zonal'))
} else {
  message("DEVMODE: OFF")
  message("FIX ONLY: ", ifelse(FIX, "ON", "OFF"))
  library(hydrofabric)
}

sf::sf_use_s2(FALSE)

# Runner Parameters -------------------------------------------------------

base = "/Volumes/MyBook/conus-hydrofabric"
version = "v20.1"

vpus = vpu_boundaries$VPUID[1:21]

ideal_size_sqkm = 10
min_length_km = 1 
min_area_sqkm = 3

nexus_prefix = "nex-"
terminal_nexus_prefix = "tnx-"
coastal_nexus_prefix = "cnx-"
internal_nexus_prefix = "inx-"
catchment_prefix = "cat-"
waterbody_prefix = "wb-"

community_hl_types =  c('HUC12', 'Gages', 'TE', 'NID', "WBOut")

# Input Datasets ----------------------------------------------------------------

coastal_gages = glue('{base}/GAGE_SUMMARY.csv')

fim_ahps      = '/Users/mjohnson/Downloads/nws_lid.gpkg'

rl_file       = glue("{base}/RouteLink_nwm_v2_2_3.parquet")

huc12_cw      = glue("{base}/huc12_nhdplusv2_cw.parquet")

gs_file       = 'https://code.usgs.gov/wma/nhgf/reference-hydrofabric/-/raw/04cd22f6b5f3f53d10c0b83b85a21d2387dfb6aa/workspace/cache/rpu_vpu_out.csv'

ms_lu         = 'https://github.com/internetofwater/ref_rivers/releases/download/v2.1/mainstem_lookup.csv'

dem_file      = "/Volumes/MyBook/3DEP/usgs_1_250mm.tif"

slope_file    = "/Volumes/MyBook/3DEP/usgs_250m_slope.tif"

aspect_file   = "/Volumes/MyBook/3DEP/usgs_250m_aspect.tif"

imp_file      = '/Volumes/MyBook/imperv_250m.tif'

twi_file      = '/Volumes/MyBook/twi_250m.tif'

lake_path     = '/Volumes/Transcend/nwmCONUS-v216/LAKEPARM_CONUS.nc'

camels        = '/Users/mjohnson/github/hydrofabric_attributes/data/camels_compiled.csv'

hydroatlas    = glue('{base}/hydroatlas_vars.parquet')

nwm_dir       = "/Volumes/Transcend/nwmCONUS-v216"

# Built Datasets ----------------------------------------------------------------

full_hl   = glue("{base}/conus_hl.gpkg")

conus_net = glue("{base}/{version}/conus_net.parquet")

conus_gpkg = glue("{base}/{version}/conus.gpkg")

atts = glue("{base}/{version}/model_attributes.parquet")

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

base_reference_features = glue("{base}/reference")
  dir.create(base_reference_features, showWarnings = FALSE)

base_reference = glue("{base}/reference")
  dir.create(base_reference, showWarnings = FALSE)

base_refactored = glue("{base}/refactored")
  dir.create(base_refactored, showWarnings = FALSE)
  
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

pipeline = data.frame(
  vpus            = vpus,
  uniform         = glue("{base_uniform}/uniform_{vpus}.gpkg"),
  uniform_global  = glue("{base_global_uniform}/uniform_{vpus}.gpkg"),
  nextgen         = glue('{base_gpkg}/nextgen_{vpus}.gpkg'),
  cfe             = glue("{base_cfe}/cfe_noahowp_{vpus}.parquet"),
  atts            = glue("{base_atts}/nextgen_{vpus}.parquet"),
  refactored_gpkg = sapply(1:length(vpus), FUN = \(x){ get_hydrofabric(vpus[x], type = "refactor",  dir = base_refactored) }),
  reference_gpkg  = sapply(1:length(vpus), FUN = \(x){ get_hydrofabric(vpus[x], type = "reference", dir = base_reference) })
)

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
  "10U", "10179150,10179151,10179147,10179146,10179153,10179148,10179152"
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
    "07", 10018008, 2111164,
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


