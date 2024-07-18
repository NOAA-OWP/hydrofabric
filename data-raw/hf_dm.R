library(dm)

hf_dm = list(
  flowlines = data.frame(
    id = integer(1L),
    toid = integer(1L),
    mainstem = integer(1L),
    order = integer(1L),
    hydroseq = integer(1L),
    lengthkm = numeric(1L),
    areasqkm = numeric(1L),
    tot_drainage_areasqkm = numeric(1L),
    divide_id = integer(1L),
    has_divide = logical(1),
    geometry = double(1L)
  ),
  
  flowpath_attributes = data.frame(
    id = integer(1L),
    rl_gages = character(1L),
    NHDWaterbodyComid = character(1L),
    Qi = numeric(1L),
    MusK = numeric(1L),
    MusX = numeric(1L),
    n = numeric(1L),
    So = numeric(1L),
    ChSlp = numeric(1L),
    BtmWdth = numeric(1L),
    Kchan = numeric(1L),
    nCC = numeric(1L),
    TopWdthCC = numeric(1L),
    TopWdth = numeric(1L),
    length_m = numeric(1L)
  ),
  
  
  lakes = data.frame(
    id           = integer(1L),
    toid  = integer(1L),
    hl_id = integer(1L),
    hl_reference = character(1L),
    hl_link      = character(1L),
    hl_uri       = character(1L),
    Dam_Length  = numeric(1L) ,
    ifd  = numeric(1L),
    LkArea  = numeric(1L),
    LkMxE  = numeric(1L),
    OrificeA  = numeric(1L),
    OrificeC  = numeric(1L),
    OrificeE   = numeric(1L),
    time = numeric(1L),
    WeirC = numeric(1L),
    WeirE = numeric(1L),
    WeirL = numeric(1L),
    geometry = double(1L)
  ),
  
  
  divides = data.frame(
    divide_id = integer(1L),
    id   = integer(1L),
    toid = integer(1L),
    ds_id = integer(1L),
    areasqkm = numeric(1L),
    type = character(1L),
    has_flowline = logical(1),
    geometry = double(1L)
  ),
  
  nexus = data.frame(
    id      = integer(1L),
    toid    = integer(1L),
    hl_id   = character(1L),
    hl_uri      = character(1L),
    type    = character(1L),
    geometry = double(1L)
  ),
  
  # Can have 1:many, NOT! part of topology
  hydrolocations = data.frame(
    hl_id = integer(1L),
    id           = integer(1L),
    hl_reference = character(1L),
    hl_link      = character(1L),
    hl_uri       = character(1L),
    hl_position  = character(1L),
    geometry     = numeric(1L)
  ),
  
  network = data.frame(
    # Topology
    id         = integer(1L),
    toid       = integer(1L),
    # Associations
    divide_id  = integer(1L),
    ds_id = integer(1L),
    mainstem   = integer(1L),
    hydroseq   = integer(1L),
    wb_id      = integer(1L),
    hl_id      = integer(1L),
    hl_uri     = character(1L),
    # Reference
    hf_source  = integer(1),
    hf_id      = numeric(1L),
    #hf_id_part = integer(1L),
    # Description
    lengthkm = numeric(1L),
    areasqkm = numeric(1L),
    tot_drainage_areasqkm = numeric(1L),
    type = character(1L),
    vpu = character(1L)
  ),
  
  WB = data.frame(
    id    = integer(1L),
    wb_id = integer(1L),
    wb_area = numeric(1L),
    wb_source = numeric(1L),
    geometry = numeric(1L)
  )
)


meta = tibble::tribble(
  ~ Attribute,
  ~ Description,
  "id",
  "A hydrofabric specfic, globaly unique flowline identifier",
  "hf_source",
  "The source dataset for the hydrofabric development (e.g. NHDPlusV2)",
  "hf_id",
  "The unique identifier in the source hydrofabric (hf_source)",
  #"hf_id_part",     "If the original hydrofabric identifier was split, the sub part. Parts increase from outlet to inlet",
  "divide_id",
  "A hydrofabric specfic, globaly unique divide identifier",
  "mainstem",
  "A nationally unique identifier assigned to the set of flowlines that compose a stream from its headwater to its mouth",
  "hl_id",
  "A hydrofabric specifc, globaly unique, hydrologic location identifier",
  "hl_reference",
  "The hydrologic location type (e.g. Gage, NID)",
  "hl_link",
  "The unique identifier in the source hydrologic location dataset.",
  "hl_uri",
  "A comma seperated conncatination of hl_reference and hl_link that mirrors thse used in the NLDI",
  "hl_position",
  "Position of hydrolocation on flowpath (inflow, outflow, along)",
  "wb_id",
  "Water body Common Identifier from wb_source",
  "ds_id",
  "Most downstream adjacent divide. Only applicable to internal catchments",
  
  "toid",
  "The identifier of the directly downstream flowpath/flowline",
  "lengthkm",
  "The length in kilometers of the flowpath element",
  "areasqkm",
  "The area of the incremental divides",
  "tot_drainage_areasqkm",
  "The total upstream area contributing to the feature",
  
  "order",
  "Strahler stream order",
  "hydroseq",
  "VPU based hydrologic sort. Increases from downstream to upstream",
  
  "geometry",
  "Simple Features Geometry (POINT, LINESTRING, POLYGON)",
  "type",
  "Type of network feature (network, internal, coastal, connector)",
  
  "wb_area",
  "Waterbody area",
  "wb_source",
  "Waterbody source",
  
  "has_divide",
  "Does an abstract catchment have a divide realization",
  "has_flowline",
  "Does an abstract catchment have a flowline/flowpath realization",
  "vpu",
  "A processing unit used to segment large scale hydrofabrics",
  
  
  'rl_gages' ,
  'NHD Gage Event ID from SOURCE_FEA field in Gages feature class',
  'NHDWaterbodyComid' ,
  'ComID of NHDWaterbody feature associated using spatial join (intersection) between NHDFlowline_network and Waterbodies"',
  'Qi' ,
  'Length weighted Initial flow in link (CMS)',
  'MusK' ,
  'Length weighted  Muskingum routing time (s)',
  'MusX' ,
  'Length weighted  Muskingum weighting coefficient',
  'n' ,
  "Length weighted Manning's roughness" ,
  'So' ,
  "Slope computed from the aggregated flow network",
  'ChSlp' ,
  "Length weighted Channel side slope" ,
  'BtmWdth' ,
  "Length weighted Bottom width of channel (m)" ,
  'Kchan' ,
  "Length weighted channel conductivity",
  'nCC' ,
  "Length weighted Compound Channel Manning's n",
  'TopWdthCC' ,
  "Compound Channel Top Width (m)",
  'TopWdth' ,
  "Length weighted Top Width (m)",
  'length_m' ,
  "Length computed from the aggregated flow network",
  
  'Dam_Length'  ,
  "Dam Length (m)" ,
  'ifd'  ,
  "Initial fraction water depth",
  'LkArea'  ,
  "Lake area (sq. km)",
  'LkMxE'  ,
  "Maximum lake elevation (m ASL)",
  'OrificeA'  ,
  "Orifice cross-sectional area (sq. m)",
  'OrificeC'  ,
  "Orifice coefficient",
  'OrificeE'   ,
  "Orifice elevation (m ASL)",
  'time' ,
  "time of measurement",
  'WeirC' ,
  "Weir coefficient" ,
  'WeirE' ,
  "Weir elevation (m ASL)" ,
  'WeirL' ,
  "Weir length (m)"
)

x = unique(unlist(sapply(hf_dm, names)))


if (!all(length(x[!x %in% meta$Attribute]) == 0, length(meta$Attribute[!meta$Attribute %in% x]) == 0)) {
  stop()
}

dm = dm::dm(
  network = hf_dm$network,
  flowlines = hf_dm$flowlines,
  hydrolocations = hf_dm$hydrolocations,
  divides = hf_dm$divides,
  WB = hf_dm$WB,
  lakes = hf_dm$lakes,
  nexus = hf_dm$nexus,
  flowpath_attributes = hf_dm$flowpath_attributes
) %>%
  dm_add_pk(flowlines, id)  %>%
  dm_add_pk(hydrolocations, hl_id)  %>%
  dm_add_pk(WB, wb_id)  %>%
  dm_add_pk(divides, divide_id)  %>%
  dm_add_pk(nexus, id)  %>%
  dm_add_pk(flowpath_attributes, id)  %>%
  dm_add_pk(lakes, id)  %>%
  dm_set_colors(
    red = flowlines,
    red = divides,
    red = hydrolocations,
    gray = WB,
    red  = network,
    blue = flowpath_attributes,
    blue = nexus,
    blue = lakes
  )

hf_dm = list(dm = dm, meta = meta)
usethis::use_data(hf_dm, overwrite = TRUE)
