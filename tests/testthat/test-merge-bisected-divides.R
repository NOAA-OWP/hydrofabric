# Divide parts split by geometry cleaning must be rejoined, not discarded.
suppressMessages({ library(sf) })

# Two disjoint squares sharing one divide_id: what a bisected divide looks like
# once cast to POLYGON. Areas are deliberately lopsided (a sliver plus the real
# catchment) because the discarding version could keep either one.
.bisected <- function(sliver_first = TRUE) {
  big    <- st_polygon(list(cbind(c(0, 2, 2, 0, 0), c(0, 0, 2, 2, 0))))     # 4
  sliver <- st_polygon(list(cbind(c(9, 9.01, 9.01, 9, 9), c(0, 0, 0.01, 0.01, 0))))
  g <- if (sliver_first) st_sfc(sliver, big, crs = 5070) else st_sfc(big, sliver, crs = 5070)
  st_sf(divide_id = c("cat-1", "cat-1"), areasqkm = as.numeric(st_area(g)) / 1e6,
        geom = g)
}

test_that("bisected parts are rejoined into one multipart divide", {
  d <- .bisected()
  r <- hydrofabric:::.merge_bisected_divides(d)
  expect_identical(r$n_merged, 1L)
  expect_identical(nrow(r$divides), 1L)
  expect_identical(as.character(st_geometry_type(r$divides)), "MULTIPOLYGON")
})

test_that("no area is lost, whichever part came first", {
  for (sf_first in c(TRUE, FALSE)) {
    d <- .bisected(sliver_first = sf_first)
    before <- as.numeric(sum(st_area(d)))
    r <- hydrofabric:::.merge_bisected_divides(d)
    expect_equal(as.numeric(sum(st_area(r$divides))), before, tolerance = 1e-6)
  }
})

test_that("the rejoined geometry is writable (sfc class not left stale)", {
  d <- .bisected()
  r <- hydrofabric:::.merge_bisected_divides(d)
  p <- tempfile(fileext = ".gpkg")
  expect_silent(st_write(r$divides, p, layer = "divides", quiet = TRUE))
  back <- st_read(p, "divides", quiet = TRUE)
  expect_identical(nrow(back), 1L)
  expect_equal(as.numeric(sum(st_area(back))), as.numeric(sum(st_area(d))),
               tolerance = 1e-6)
})

test_that("distinct divides are never merged together", {
  g <- st_sfc(
    st_polygon(list(cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)))),
    st_polygon(list(cbind(c(5, 6, 6, 5, 5), c(0, 0, 1, 1, 0)))), crs = 5070)
  d <- st_sf(divide_id = c("cat-1", "cat-2"), areasqkm = 1, geom = g)
  r <- hydrofabric:::.merge_bisected_divides(d)
  expect_identical(r$n_merged, 0L)
  expect_identical(nrow(r$divides), 2L)
  expect_identical(sort(r$divides$divide_id), c("cat-1", "cat-2"))
})

test_that("an already-clean multipart layer is reported unchanged", {
  d <- .bisected()
  r1 <- hydrofabric:::.merge_bisected_divides(d)
  r2 <- hydrofabric:::.merge_bisected_divides(r1$divides)
  expect_identical(r2$n_merged, 0L)
  expect_false(r2$changed)
})

test_that("an empty or column-less layer is a no-op", {
  expect_identical(hydrofabric:::.merge_bisected_divides(
    .bisected()[0, ])$n_merged, 0L)
  g <- st_sfc(st_polygon(list(cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)))), crs = 5070)
  expect_false(hydrofabric:::.merge_bisected_divides(st_sf(x = 1, geom = g))$changed)
})
