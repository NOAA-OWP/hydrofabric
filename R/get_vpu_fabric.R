#' Extract VPU GPKG from Lynker-Spatial
#' Data from a Parquet Dataset will be extracted for a VPU. 
#' Optionally this data can be written to a GPKG.
#' @param vpu VPU ID
#' @param type hydrofabric type
#' @param version hydrofabric version
#' @param source hydrofabric source (local root directory or s3 link)
#' @param outfile Data will be written to a file if supplied. 
#' Must be a gpkg file path. 
#' If NULL a list object of divides, flowpaths, and a network table will be returned.
#' @return
#' @export

get_vpu_fabric = function(vpu = "01", 
                   type = "reference",
                   version = "2.2", 
                   source = "s3://lynker-spatial/hydrofabric",
                   outfile = NULL){

  vpuid <- NULL
   
  fl = open_dataset(glue('{source}/v{version}/{type}/conus_flowlines/')) %>% 
    filter(vpuid == vpu) %>% 
    sfarrow::read_sf_dataset() %>% 
    st_set_crs(5070)
  
  div = open_dataset(glue('{source}/v{version}/{type}/conus_divides/')) %>% 
    filter(vpuid == vpu) %>% 
    sfarrow::read_sf_dataset() %>% 
    st_set_crs(5070)
  
  net = open_dataset(glue('{source}/v{version}/{type}/conus_network/')) %>% 
    filter(vpuid == vpu) %>% 
    collect() 
  
  if(is.null(outfile)){
    list(divides = div,
         flowpaths = fl,
         network = net)
  } else {
    write_sf(fl, outfile, "flowpaths")
    write_sf(div, outfile, "divides")
    write_sf(net, outfile, "network")
    return(outfile)
  }
}

x = get_vpu_fabric("01", "reference")
x = get_vpu_fabric("01", "reference", outfile = "test01.gpkg")
