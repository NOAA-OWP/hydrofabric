# check_reference_divide_pairing: find reference reaches whose downstream end
# lies outside their own divide. Most lost catchments are
# this: the split seed sits 1.1 to 3.4 km beyond the parent it is meant to cut,
# so the masked D8 flood has no cell to start from and the parent is abandoned.
# An outlet sitting marginally outside its boundary is ordinary, so the measure
# is distance, not containment.

.crp_gpkg <- function(lines, polys, line_ids, poly_ids, vpuid = NULL) {
  fl <- sf::st_sf(flowpath_id = line_ids,
                  geometry = sf::st_sfc(lines, crs = 5070))
  dv <- sf::st_sf(flowpath_id = poly_ids,
                  geometry = sf::st_sfc(polys, crs = 5070))
  if (!is.null(vpuid)) { fl$vpuid <- vpuid; dv$vpuid <- vpuid }
  p <- tempfile(fileext = ".gpkg")
  suppressWarnings(sf::st_write(fl, p, "flowpaths", quiet = TRUE))
  suppressWarnings(sf::st_write(dv, p, "divides", quiet = TRUE))
  p
}
.crp_ln <- function(...) sf::st_linestring(rbind(...))
.crp_box <- function(x0, y0, x1, y1)
  sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))

test_that("a reach ending well outside its own divide is reported", {
  # 'b' ends 2000 m east of its divide's edge.
  p <- .crp_gpkg(
    list(.crp_ln(c(100, 500), c(900, 500)),
         .crp_ln(c(100, 500), c(3000, 500))),
    list(.crp_box(0, 0, 1000, 1000), .crp_box(0, 0, 1000, 1000)),
    c("a", "b"), c("a", "b"))
  r <- suppressMessages(hydrofabric::check_reference_divide_pairing(p))
  expect_equal(r$id, "b")
  expect_equal(r$dist_m, 2000)
})

test_that("an outlet on the divide boundary is not reported", {
  # The normal case: a reach ends exactly at its catchment outlet, which lies on
  # the boundary. Distance 0, and it must not be flagged.
  p <- .crp_gpkg(
    list(.crp_ln(c(100, 500), c(1000, 500))),
    list(.crp_box(0, 0, 1000, 1000)), "a", "a")
  r <- suppressMessages(hydrofabric::check_reference_divide_pairing(p))
  expect_equal(nrow(r), 0L)
})

test_that("an endpoint just past the boundary is within tolerance", {
  # 100 m out is imprecision, not a broken pairing: the 250 m default is the
  # buffer the splitter masks with.
  p <- .crp_gpkg(
    list(.crp_ln(c(100, 500), c(1100, 500))),
    list(.crp_box(0, 0, 1000, 1000)), "a", "a")
  r <- suppressMessages(hydrofabric::check_reference_divide_pairing(p))
  expect_equal(nrow(r), 0L)
})

test_that("tolerance_m is honoured", {
  p <- .crp_gpkg(
    list(.crp_ln(c(100, 500), c(1100, 500))),
    list(.crp_box(0, 0, 1000, 1000)), "a", "a")
  r <- suppressMessages(
    hydrofabric::check_reference_divide_pairing(p, tolerance_m = 30))
  expect_equal(r$id, "a")
  expect_equal(r$dist_m, 100)
})

test_that("distance is measured pairwise by id, not against every divide", {
  # 'b' ends inside a DIFFERENT divide. Measuring against the nearest polygon
  # would score it 0; measuring against its own scores it 2000.
  p <- .crp_gpkg(
    list(.crp_ln(c(100, 500), c(900, 500)),
         .crp_ln(c(100, 500), c(3500, 500))),
    list(.crp_box(0, 0, 1000, 1000), .crp_box(0, 0, 1000, 1000),
         .crp_box(3000, 0, 4000, 1000)),
    c("a", "b"), c("a", "b", "c"))
  r <- suppressMessages(hydrofabric::check_reference_divide_pairing(p))
  expect_equal(r$id, "b")
  expect_equal(r$dist_m, 2500)
})

test_that("a reach with no divide of its own is skipped, not reported", {
  p <- .crp_gpkg(
    list(.crp_ln(c(100, 500), c(900, 500)),
         .crp_ln(c(5000, 500), c(9000, 500))),
    list(.crp_box(0, 0, 1000, 1000)), c("a", "b"), "a")
  r <- suppressMessages(hydrofabric::check_reference_divide_pairing(p))
  expect_equal(nrow(r), 0L)
})

test_that("the furthest reach is reported first", {
  p <- .crp_gpkg(
    list(.crp_ln(c(100, 500), c(2000, 500)),
         .crp_ln(c(100, 500), c(4000, 500))),
    list(.crp_box(0, 0, 1000, 1000), .crp_box(0, 0, 1000, 1000)),
    c("a", "b"), c("a", "b"))
  r <- suppressMessages(hydrofabric::check_reference_divide_pairing(p))
  expect_equal(r$id, c("b", "a"))
  expect_equal(r$dist_m, c(3000, 1000))
})

test_that("a vpuid filter restricts the population", {
  p <- .crp_gpkg(
    list(.crp_ln(c(100, 500), c(3000, 500)),
         .crp_ln(c(100, 500), c(3000, 500))),
    list(.crp_box(0, 0, 1000, 1000), .crp_box(0, 0, 1000, 1000)),
    c("a", "b"), c("a", "b"), vpuid = c("04", "18"))
  r <- suppressMessages(
    hydrofabric::check_reference_divide_pairing(p, vpuid = "04"))
  expect_equal(r$id, "a")
})
