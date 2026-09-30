# Gage drainage-area reporting.
#
# This was a tolerance test that raised a flag. The flag depended on how many
# reference reaches the aggregation bundled into a flowpath rather than on where
# the gage sat: smaller aggregates made gages fail a test they had passed while
# their delineated area was unchanged. The numbers are published; the judgement
# is the reader's.

test_that("both drainage areas and their ratio are written, with no flag", {
  f <- tempfile(fileext = ".gpkg"); on.exit(unlink(f))
  pt <- sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(1, 1)), crs = 5070)
  net <- data.frame(
    hl_reference = c("nwis-1", "nwis-1", "nwis-2"),
    flowpath_id  = c("fp-1", "fp-1", "fp-2"),
    total_dasqkm = c(2, 2, 200),
    flowline_total_dasqkm = c(1000, 3, 199), stringsAsFactors = FALSE)
  sf::st_write(sf::st_sf(hl_reference = c("nwis-1", "nwis-2"), geometry = pt),
    f, "hydrolocations", quiet = TRUE)
  con <- DBI::dbConnect(RSQLite::SQLite(), f)
  DBI::dbWriteTable(con, "network", net, overwrite = TRUE)
  DBI::dbDisconnect(con)

  qc <- hydrofabric:::.report_poi_drainage_area(f)
  expect_false("poi_da_flag" %in% names(qc))
  expect_false("poi_da_resnap_flowpath" %in% names(qc))
  a <- qc[qc$hl_reference == "nwis-1", ]
  expect_equal(a$poi_da_ref_sqkm, 1000)   # largest member reference DA
  expect_equal(a$poi_da_ngen_sqkm, 2)     # delineated DA of the flowpath
  expect_equal(a$poi_da_ratio, 0.998)     # |1000 - 2| / 1000, reported not judged
})

test_that("a reference area below one DEM cell yields NA, not a huge ratio", {
  # A reference DA smaller than the grid the fabric is delineated on is absent,
  # not small; dividing by it would turn the product's own area into a verdict.
  f <- tempfile(fileext = ".gpkg"); on.exit(unlink(f))
  sf::st_write(sf::st_sf(hl_reference = "nwis-3",
    geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 5070)),
    f, "hydrolocations", quiet = TRUE)
  con <- DBI::dbConnect(RSQLite::SQLite(), f)
  DBI::dbWriteTable(con, "network", data.frame(hl_reference = "nwis-3",
    flowpath_id = "fp-3", total_dasqkm = 300,
    flowline_total_dasqkm = 1e-6, stringsAsFactors = FALSE), overwrite = TRUE)
  DBI::dbDisconnect(con)
  qc <- hydrofabric:::.report_poi_drainage_area(f)
  expect_true(is.na(qc$poi_da_ratio))
  expect_equal(qc$poi_da_ngen_sqkm, 300)
})
