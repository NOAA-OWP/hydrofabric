# HydroLAKES is CC BY 4.0, so a build that reads HydroLAKES polygons must credit
# them. The credit is keyed on the `source` column of the lake set the build read,
# and carried from the VPU products into the merged domain product.

.lake_fixture <- function(source) {
  f <- withr::local_tempfile(fileext = ".gpkg", .local_envir = parent.frame())
  sf::st_write(sf::st_sf(lake_id = 1L, source = source,
    geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 3338)),
    f, "lakes", quiet = TRUE)
  f
}

test_that("a HydroLAKES lake set is credited", {
  f <- .lake_fixture("HydroLAKES v1.0")
  expect_identical(hydrofabric:::.hf_lakes_attribution(f),
                   hydrofabric:::.HF_ATTR_HYDROLAKES)
})

test_that("a lake set from another source carries no HydroLAKES credit", {
  f <- .lake_fixture("NHDPlus HR")
  expect_length(hydrofabric:::.hf_lakes_attribution(f), 0L)
  expect_length(hydrofabric:::.hf_lakes_attribution(NULL), 0L)
  expect_length(hydrofabric:::.hf_lakes_attribution("/no/such/lakes.gpkg"), 0L)
})

test_that("the HydroLAKES credit is carried into a merged product", {
  g <- withr::local_tempfile(fileext = ".gpkg")
  sf::st_write(sf::st_sf(id = 1L, geometry = sf::st_sfc(sf::st_point(c(0, 0)),
    crs = 3338)), g, "flowpaths", quiet = TRUE)
  hydrofabric:::.hf_stamp_outputs(g, vpuid = "19",
    lakes_path = .lake_fixture("HydroLAKES v1.0"))
  expect_true(hydrofabric:::.HF_ATTR_HYDROLAKES %in%
                hydrofabric:::.hf_grid_attribution_inherited(g))
})
