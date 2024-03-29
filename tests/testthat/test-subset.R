nldi_feature = list('featureSource' = "nwissite",  featureID = "USGS-05428500")
pt = st_coordinates(get_nldi_feature(nldi_feature))

test_that("Reference Subset Works", {
  
  s3 = 's3://lynker-spatial/01_reference/'
  
  expect_error(subset_network(id = "wb-1282898", base_s3 = s3, cache_dir = "data"))
  
  expect_error(subset_network(hl_id = "USGS-05428500", base_s3 = s3, cache_dir = "data"))
  
  xx_comid = subset_network(comid = 13293750, base_s3 = s3, 
                            cache_dir = "data", 
                            cache_overwrite = TRUE)
    expect_equal(length(xx_comid), 2)
    expect_equal(nrow(xx_comid$reference_flowline), 169)
    expect_equal(nrow(xx_comid$reference_catchment), 168)
  
  xx_nldi = subset_network(nldi_feature = nldi_feature,  base_s3 = s3, cache_dir = "data")
    expect_equal(length(xx_nldi), 2)
    expect_equal(nrow(xx_nldi$reference_flowline), 169)
    expect_equal(nrow(xx_nldi$reference_catchment), 168)
  
  xx_loc = subset_network(xy = pt, base_s3 = s3)
    expect_equal(length(xx_loc), 2) 
    expect_equal(nrow(xx_loc$reference_flowline), 169)
    expect_equal(nrow(xx_loc$reference_catchment), 168)
  
})


test_that("Refactor Subset Works", {
  
  s3 = 's3://lynker-spatial/02_refactored/'
  
  expect_error(subset_network(id = "wb-1282898", base_s3 = s3, cache_dir = "data"))
  
  expect_error(subset_network(hl_id = "USGS-05428500", base_s3 = s3, cache_dir = "data"))
  
  xx_comid = subset_network(comid = 13293750, base_s3 = s3, 
                            cache_dir = "data", 
                            cache_overwrite = F)
  expect_equal(length(xx_comid), 2)
  expect_equal(nrow(xx_comid$refactored_flowpaths), 121)
  expect_equal(nrow(xx_comid$refactored_flowpaths), 121)
  
  xx_nldi = subset_network(nldi_feature = nldi_feature,  base_s3 = s3, cache_dir = "data")
  expect_equal(length(xx_nldi), 2)
  expect_equal(nrow(xx_nldi$refactored_flowpaths), 121)
  expect_equal(nrow(xx_nldi$refactored_flowpaths), 121)
  
  xx_loc = subset_network(xy = pt, base_s3 = s3)
  expect_equal(length(xx_loc), 2) 
  expect_equal(nrow(xx_loc$refactored_flowpaths), 121)
  expect_equal(nrow(xx_loc$refactored_flowpaths), 121)
  
})



test_that("NextGen Subset Works", {
  
  s3 = 's3://lynker-spatial/v20.1/'
  
  xx_id = subset_network(id = "nex-1089012", base_s3 = s3, cache_dir = "data", cache_overwrite = TRUE)
    expect_equal(length(xx_id), 5)
    expect_equal(nrow(xx_id$flowpaths), 58)
    expect_equal(nrow(xx_id$divides), 58)
  
  xx_hl = subset_network(hl_uri = "Gages-05428500", base_s3 = s3, cache_dir = "data")
    expect_equal(length(xx_hl), 5)
    expect_equal(nrow(xx_hl$flowpaths), 58)
    expect_equal(nrow(xx_hl$divides), 58)
  
  xx_comid = subset_network(comid = 13293750, base_s3 = s3, cache_dir = "data")
    expect_equal(length(xx_comid), 5)
    expect_equal(nrow(xx_comid$flowpaths), 58)
    expect_equal(nrow(xx_comid$divides), 58)
  
  xx_nldi = subset_network(nldi_feature = nldi_feature, base_s3 = s3, cache_dir = "data")
    expect_equal(length(xx_nldi), 5)
    expect_equal(nrow(xx_nldi$flowpaths), 58)
    expect_equal(nrow(xx_nldi$divides), 58)
  
  xx_loc = subset_network(xy = pt, base_s3 = s3, cache_dir = "data")
    expect_equal(length(xx_loc), 5) 
    expect_equal(nrow(xx_loc$flowpaths), 58)
    expect_equal(nrow(xx_loc$divides), 58)
  
})


subset_network(nldi_feature = 'https://reference.geoconnex.us/collections/gages/items?provider_id=06752260')
subset_network(hl_uri = 'https://reference.geoconnex.us/collections/gages/items?provider_id=06752260')
