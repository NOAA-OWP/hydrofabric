gpkg = '/Users/mikejohnson/hydrofabric/conus_nextgen.gpkg'

subset_fabric <- get_subset(gpkg = '/Users/mikejohnson/hydrofabric/conus_nextgen.gpkg', 
                            comid = 101,
                            outfile = "inst/extdata/comid-101-subset.gpkg")

subset_fabric <- get_subset(gpkg = '/Users/mikejohnson/hydrofabric/CONUS/v2.2/reference/ls_reference.gpkg', 
                            comid = 101,
                            lyrs = c("divides", "flowpaths"),
                            outfile = "inst/extdata/reference-comid-101-subset.gpkg")


read_sf('/Users/mikejohnson/hydrofabric/CONUS/v2.2/reference/community_hydrolocations.gpkg') |> 
  st_as_sf(coords = c("X", "Y"), crs = 5070) |> 
  st_filter(read_sf("inst/extdata/reference-comid-101-subset.gpkg", "divides"))
