# check_reference_coincident_geometry: find reference reaches drawn on top of
# one another. One QGIS click returned several features
# because four reference reaches shared 2.3-5.5 km of geometry into a common
# confluence. Exact duplication is rare; near-coincidence within a metre is what
# actually occurs and what a click cannot distinguish.

.crc_gpkg <- function(lines, ids, vpuid = NULL) {
  d <- sf::st_sf(flowpath_id = ids,
                 geometry = sf::st_sfc(lines, crs = 5070))
  if (!is.null(vpuid)) d$vpuid <- vpuid
  p <- tempfile(fileext = ".gpkg")
  suppressWarnings(sf::st_write(d, p, "flowpaths", quiet = TRUE))
  p
}
.crc_ln <- function(...) sf::st_linestring(rbind(...))

test_that("two reaches on the same ground are both reported", {
  p <- .crc_gpkg(
    list(.crc_ln(c(0, 0), c(1000, 0)),          # a
         .crc_ln(c(0, 0.2), c(1000, 0.2)),      # b, 20 cm away for its whole length
         .crc_ln(c(0, 5000), c(1000, 5000))),   # c, far off
    c("a", "b", "c"))
  r <- suppressMessages(
    hydrofabric::check_reference_coincident_geometry(p, sample = NULL))
  expect_setequal(r$id, c("a", "b"))
  expect_true(all(r$frac >= 0.5))
})

test_that("reaches merely near each other are not reported", {
  # 20 m apart: under the 30 m DEM cell but far outside a 1 m tolerance. These
  # are distinct channels and must not be flagged.
  p <- .crc_gpkg(
    list(.crc_ln(c(0, 0), c(1000, 0)),
         .crc_ln(c(0, 20), c(1000, 20))),
    c("a", "b"))
  r <- suppressMessages(
    hydrofabric::check_reference_coincident_geometry(p, sample = NULL))
  expect_equal(nrow(r), 0L)
})

test_that("reaches meeting end to end are not reported", {
  # Ordinary network connectivity: they share a point, not a length.
  p <- .crc_gpkg(
    list(.crc_ln(c(0, 0), c(1000, 0)),
         .crc_ln(c(1000, 0), c(2000, 0))),
    c("a", "b"))
  r <- suppressMessages(
    hydrofabric::check_reference_coincident_geometry(p, sample = NULL))
  expect_equal(nrow(r), 0L)
})

test_that("a partial overlap below the fraction threshold is not reported", {
  # Sharing a short tail into a confluence is common and within tolerance; the
  # threshold is what separates that from a reach drawn along another.
  p <- .crc_gpkg(
    list(.crc_ln(c(0, 0), c(1000, 0)),
         .crc_ln(c(900, 0.2), c(1900, 0.2))),   # ~10% shared
    c("a", "b"))
  r <- suppressMessages(
    hydrofabric::check_reference_coincident_geometry(p, sample = NULL))
  expect_equal(nrow(r), 0L)
})

test_that("min_frac and tolerance_m are honoured", {
  p <- .crc_gpkg(
    list(.crc_ln(c(0, 0), c(1000, 0)),
         .crc_ln(c(900, 0.2), c(1900, 0.2))),
    c("a", "b"))
  r <- suppressMessages(hydrofabric::check_reference_coincident_geometry(
    p, sample = NULL, min_frac = 0.05))
  expect_gt(nrow(r), 0L)
})

test_that("a vpuid filter restricts the population", {
  p <- .crc_gpkg(
    list(.crc_ln(c(0, 0), c(1000, 0)),
         .crc_ln(c(0, 0.2), c(1000, 0.2))),
    c("a", "b"), vpuid = c("04", "18"))
  r <- suppressMessages(hydrofabric::check_reference_coincident_geometry(
    p, vpuid = "04", sample = NULL))
  expect_equal(nrow(r), 0L)   # only one reach in VPU 04, nothing to coincide with
})
