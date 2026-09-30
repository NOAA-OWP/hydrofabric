test_that("every domain but Hawaii and PRVI credits TDX and HydroSHEDS", {
  withr::local_envvar(HF_DATA_ATTRIBUTION = "")
  owed <- function(v) {
    s <- hydrofabric:::.hf_attribution(v)
    c(nga = grepl("NGA", s, fixed = TRUE),
      hs  = grepl("HydroSHEDS", s, fixed = TRUE))
  }
  for (v in c("01", "04", "17", "19", "22")) expect_true(all(owed(v)), info = v)
  # Sub-VPU tiles resolve to their parent domain rather than falling through.
  expect_true(all(owed("19a")))
  # Hawaii and PRVI are NHDPlus HR throughout: a US Government work, no credit owed.
  for (v in c("20", "21")) expect_false(any(owed(v)), info = v)
  # An unknown or absent vpuid credits both: over-crediting is cosmetic,
  # under-crediting is a licence breach.
  expect_true(all(owed(NA_character_)))
  expect_true(all(owed("")))
})

test_that("an operator override replaces the derived attribution", {
  withr::local_envvar(HF_DATA_ATTRIBUTION = "Custom notice.")
  expect_equal(hydrofabric:::.hf_attribution("01"), "Custom notice.")
})

test_that("a grid credit is not repeated when policy already carries it", {
  withr::local_envvar(HF_DATA_ATTRIBUTION = "")
  f <- tempfile(fileext = ".gpkg")
  con <- DBI::dbConnect(RSQLite::SQLite(), f)
  DBI::dbExecute(con, "CREATE TABLE gpkg_contents (table_name TEXT)")
  DBI::dbDisconnect(con)
  hydrofabric:::.hf_stamp_outputs(f, vpuid = "01",
    grid_paths = "rasters/hydrosheds_v2/hydrosheds_conus_dir.vrt")
  con <- DBI::dbConnect(RSQLite::SQLite(), f)
  a <- jsonlite::fromJSON(DBI::dbGetQuery(con,
    "SELECT metadata FROM gpkg_metadata WHERE mime_type='application/json'")$metadata[1])$attribution
  lic <- DBI::dbGetQuery(con,
    "SELECT metadata FROM gpkg_metadata WHERE md_standard_uri LIKE '%spdx%'")$metadata[1]
  DBI::dbDisconnect(con)
  expect_equal(lengths(regmatches(a, gregexpr("HydroSHEDS", a, fixed = TRUE))), 1L)
  expect_equal(lic, "CC-BY-SA-4.0")
})

test_that("the data license is uniform across the product", {
  withr::local_envvar(HF_DATA_LICENSE = "")
  # TDX-Hydro's ShareAlike term compels it where TDX is in the lineage; the rest
  # take it for consistency. Credit varies by source, the license does not.
  for (v in c("01", "04", "19", "19a", "20", "21", "22", NA_character_))
    expect_equal(hydrofabric:::.hf_data_license(v), "CC-BY-SA-4.0",
                 info = ifelse(is.na(v), "NA", v))
})

test_that("an operator override replaces the default license", {
  withr::local_envvar(HF_DATA_LICENSE = "CC0-1.0")
  expect_equal(hydrofabric:::.hf_data_license("20"), "CC0-1.0")
})

test_that("fac_source names the grid a fallback value came from", {
  expect_equal(hydrofabric:::.fac_source_label("rasters/hydrosheds_v2/hydrosheds_conus_acc.vrt"),
               "hydrosheds")
  expect_equal(hydrofabric:::.fac_source_label("rasters/hydrosheds_v2/na_acc_1s_v2.vrt"),
               "hydrosheds")
  expect_equal(hydrofabric:::.fac_source_label("rasters/gedtm30/gedtm_fac.vrt"), "gedtm30")
  # An unrecognised grid is still reported as not-primary rather than mislabelled.
  expect_equal(hydrofabric:::.fac_source_label("rasters/other/mystery.tif"), "fallback")
  expect_equal(hydrofabric:::.fac_source_label(""), "fallback")
  expect_equal(hydrofabric:::.fac_source_label(NULL), "fallback")
})
