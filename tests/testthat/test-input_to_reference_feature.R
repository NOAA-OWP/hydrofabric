# Initialize constant variables
source = 's3://lynker-spatial/hydrofabric'
type = "reference"
hf_version = "2.2"
outfile = 'subset_test.gpkg'

# Test with an comid and network file 
test_that("Terminal point check for comid passed", {

    comid <- 718648
  
    expect_error(get_subset(outfile = outfile, comid = comid, source = source, type = type, hf_version = hf_version))

    subset_flowpaths <- st_read("subset_test.gpkg", layer = "flowpaths")
    subset_divides <- st_read("subset_test.gpkg", layer = "divides")
    subset_network <- st_read("subset_test.gpkg", layer = "network")
    
    # Flowlines
    expect_equal(length(subset_flowpaths), 10)
    expect_equal(nrow(subset_flowpaths), 5)
    expect_equal(unique(subset_flowpaths$vpuid), "01")
    expect_equal(min(subset_flowpaths$mainstemlp), 1923133)
    expect_equal(max(subset_flowpaths$mainstemlp), 1923137) 

    # Divides
    expect_equal(length(subset_divides), 4)
    expect_equal(nrow(subset_divides), 5)
    expect_equal(unique(subset_divides$vpuid), "01")
    expect_equal(min(subset_divides$divide_id), 718618)
    expect_equal(max(subset_divides$divide_id), 718654) 

    # Network
    expect_equal(length(subset_network), 15)
    expect_equal(nrow(subset_network), 5)
    expect_equal(unique(subset_network$vpuid), "01")
    expect_equal(min(subset_network$hf_id), 718618)
    expect_equal(max(subset_network$hf_id), 718654) 
})

# Test with an id and network file 
test_that("Terminal point check for id passed", {

    id <- 923
  
    expect_error(get_subset(outfile = outfile, id = id, source = source, type = type, hf_version = hf_version))

    subset_flowpaths <- st_read("subset_test.gpkg", layer = "flowpaths")
    subset_divides <- st_read("subset_test.gpkg", layer = "divides")
    subset_network <- st_read("subset_test.gpkg", layer = "network")
    
    # Flowlines
    expect_equal(length(subset_flowpaths), 10)
    expect_equal(nrow(subset_flowpaths), 2)
    expect_equal(unique(subset_flowpaths$vpuid), "01")
    expect_equal(min(subset_flowpaths$mainstemlp), 1819851)
    expect_equal(max(subset_flowpaths$mainstemlp), 1819851) 

    # Divides
    expect_equal(length(subset_divides), 4)
    expect_equal(nrow(subset_divides), 2)
    expect_equal(unique(subset_divides$vpuid), "01")
    expect_equal(min(subset_divides$divide_id), 915)
    expect_equal(max(subset_divides$divide_id), 1359) 

    # Network
    expect_equal(length(subset_network), 15)
    expect_equal(nrow(subset_network), 2)
    expect_equal(unique(subset_network$vpuid), "01")
    expect_equal(min(subset_network$hf_id), 915)
    expect_equal(max(subset_network$hf_id), 1359) 
})

# Test with an nldi_feature
test_that("Terminal point check for nldi_feature passed", {

    nldi_feature <- list(featureSource = "nwissite", featureID = "USGS-08279500")
  
    expect_error(get_subset(outfile = outfile, nldi_feature = nldi_feature, source = source, type = type, hf_version = hf_version))

    subset_flowpaths <- st_read("subset_test.gpkg", layer = "flowpaths")
    subset_divides <- st_read("subset_test.gpkg", layer = "divides")
    subset_network <- st_read("subset_test.gpkg", layer = "network")
    
    # Flowlines
    expect_equal(length(subset_flowpaths), 10)
    expect_equal(nrow(subset_flowpaths), 3364)
    expect_equal(unique(subset_flowpaths$vpuid), "13")
    expect_equal(min(subset_flowpaths$mainstemlp), 1)
    expect_equal(max(subset_flowpaths$mainstemlp), 38435) 

    # Divides
    expect_equal(length(subset_divides), 4)
    expect_equal(nrow(subset_divides), 3354)
    expect_equal(unique(subset_divides$vpuid), "13")
    expect_equal(min(subset_divides$divide_id), 17029298)
    expect_equal(max(subset_divides$divide_id), 943020360) 

    # Network
    expect_equal(length(subset_network), 15)
    expect_equal(nrow(subset_network), 3364)
    expect_equal(unique(subset_network$vpuid), "13")
    expect_equal(min(subset_network$hf_id), 17029298)
    expect_equal(max(subset_network$hf_id), 943020360) 
})

# Test with an poi_id and network file 
test_that("Terminal point check for poi_id passed", {

    poi_id <- 74719
  
    expect_error(get_subset(outfile = outfile, poi_id = poi_id, source = source, type = type, hf_version = hf_version))

    subset_flowpaths <- st_read("subset_test.gpkg", layer = "flowpaths")
    subset_divides <- st_read("subset_test.gpkg", layer = "divides")
    subset_network <- st_read("subset_test.gpkg", layer = "network")
    
    # Flowlines
    expect_equal(length(subset_flowpaths), 10)
    expect_equal(nrow(subset_flowpaths), 1083)
    expect_equal(unique(subset_flowpaths$vpuid), "08")
    expect_equal(min(subset_flowpaths$mainstemlp), 473130)
    expect_equal(max(subset_flowpaths$mainstemlp), 1583826) 

    # Divides
    expect_equal(length(subset_divides), 4)
    expect_equal(nrow(subset_divides), 1078)
    expect_equal(unique(subset_divides$vpuid), "08")
    expect_equal(min(subset_divides$divide_id), 653948)
    expect_equal(max(subset_divides$divide_id), 938020449) 

    # Network
    expect_equal(length(subset_network), 15)
    expect_equal(nrow(subset_network), 1083)
    expect_equal(unique(subset_network$vpuid), "08")
    expect_equal(min(subset_network$hf_id), 74719)
    expect_equal(max(subset_network$hf_id), 938020449) 
})

# Test with an hl_uri and network file 
test_that("Terminal point check for hl_uri passed", {

    hl_uri <- 'HUC12-010100070603'
  
    expect_error(get_subset(outfile = outfile, hl_uri = hl_uri, source = source, type = type, hf_version = hf_version))

    subset_flowpaths <- st_read("subset_test.gpkg", layer = "flowpaths")
    subset_divides <- st_read("subset_test.gpkg", layer = "divides")
    subset_network <- st_read("subset_test.gpkg", layer = "network")
    
    # Flowlines
    expect_equal(length(subset_flowpaths), 10)
    expect_equal(nrow(subset_flowpaths), 130)
    expect_equal(unique(subset_flowpaths$vpuid), "01")
    expect_equal(min(subset_flowpaths$mainstemlp), 1921949)
    expect_equal(max(subset_flowpaths$mainstemlp), 1922079) 

    # Divides
    expect_equal(length(subset_divides), 4)
    expect_equal(nrow(subset_divides), 127)
    expect_equal(unique(subset_divides$vpuid), "01")
    expect_equal(min(subset_divides$divide_id), 4288361)
    expect_equal(max(subset_divides$divide_id), 4292703) 

    # Network
    expect_equal(length(subset_network), 15)
    expect_equal(nrow(subset_network), 130)
    expect_equal(unique(subset_network$vpuid), "01")
    expect_equal(min(subset_network$hf_id), 4288361)
    expect_equal(max(subset_network$hf_id), 4292703) 
})