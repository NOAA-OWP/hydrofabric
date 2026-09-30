test_that("a grid that does not cover the domain is rejected", {
  skip_if_not_installed("terra")
  ref <- tempfile(fileext = ".gpkg")
  # A reference whose extent sits in one place ...
  con <- DBI::dbConnect(RSQLite::SQLite(), ref)
  DBI::dbExecute(con, "CREATE TABLE gpkg_contents (table_name TEXT, min_x REAL,
                       min_y REAL, max_x REAL, max_y REAL, srs_id INTEGER)")
  DBI::dbExecute(con, "INSERT INTO gpkg_contents VALUES ('flowpaths', 0, 0, 10, 10, 4326)")
  DBI::dbExecute(con, "CREATE TABLE gpkg_spatial_ref_sys (srs_id INTEGER,
                       organization TEXT, organization_coordsys_id INTEGER)")
  DBI::dbExecute(con, "INSERT INTO gpkg_spatial_ref_sys VALUES (4326, 'EPSG', 4326)")
  DBI::dbDisconnect(con)

  covering <- tempfile(fileext = ".tif")
  elsewhere <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::rast(terra::ext(-5, 20, -5, 20), resolution = 1,
    crs = "EPSG:4326", vals = 1), covering)
  terra::writeRaster(terra::rast(terra::ext(100, 120, 40, 60), resolution = 1,
    crs = "EPSG:4326", vals = 1), elsewhere)

  expect_silent(suppressMessages(
    check_grid_coverage(ref, covering, covering, strict = TRUE)))
  expect_error(check_grid_coverage(ref, elsewhere, elsewhere, strict = TRUE),
               "must cover the domain")

  # A fallback that cannot reach the domain is reported, not fatal: it never fires.
  expect_silent(suppressMessages(
    check_grid_coverage(ref, covering, covering, elsewhere, elsewhere, strict = TRUE)))

  cov <- suppressMessages(check_grid_coverage(ref, covering, covering, strict = FALSE))
  expect_equal(unname(cov[["fdr"]]), 1)
})
