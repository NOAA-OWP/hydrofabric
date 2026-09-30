# An unowned divide must be absorbed by a neighbour, never deleted: dropping it
# removes that land from the layer, which a sum-based conservation check misses.
suppressMessages({ library(sf) })

.sq <- function(x0, y0, w = 1, h = 1)
  st_polygon(list(cbind(c(x0, x0 + w, x0 + w, x0, x0),
                        c(y0, y0, y0 + h, y0 + h, y0))))

# host A [0,1], orphan [1,2] (shares a full edge with A and B), host B [2,3]
.three <- function() {
  g <- st_sfc(.sq(0, 0), .sq(1, 0), .sq(2, 0), crs = 5070)
  st_sf(divide_id = c("a", "orphan", "b"),
        areasqkm = as.numeric(st_area(g)) / 1e6, geom = g)
}

test_that("an orphan is absorbed and total coverage is preserved", {
  d <- .three()
  before <- as.numeric(st_area(st_union(st_geometry(d))))
  r <- hydrofabric:::.absorb_orphan_divides(d, c(FALSE, TRUE, FALSE))
  expect_identical(r$n_absorbed, 1L)
  expect_identical(nrow(r$div), 2L)
  expect_false("orphan" %in% r$div$divide_id)
  expect_equal(as.numeric(st_area(st_union(st_geometry(r$div)))), before,
               tolerance = 1e-6)
})

test_that("areasqkm is refreshed on the absorbing divide", {
  d <- .three()
  r <- hydrofabric:::.absorb_orphan_divides(d, c(FALSE, TRUE, FALSE))
  host <- r$div[r$div$divide_id %in% c("a", "b"), ]
  expect_true(any(host$areasqkm > 1.5e-6))   # one host grew by the orphan's area
  expect_equal(sum(host$areasqkm),
               as.numeric(sum(st_area(host))) / 1e6, tolerance = 1e-9)
})

test_that("the neighbour sharing the most boundary wins", {
  # orphan shares a 1-unit edge with `big` and only 0.25 with `small`.
  g <- st_sfc(.sq(0, 0, 1, 1),          # big: full shared edge with orphan
              .sq(1, 0, 1, 1),          # orphan
              .sq(2, 0, 1, 0.25),       # small: partial shared edge
              crs = 5070)
  d <- st_sf(divide_id = c("big", "orphan", "small"),
             areasqkm = as.numeric(st_area(g)) / 1e6, geom = g)
  r <- hydrofabric:::.absorb_orphan_divides(d, c(FALSE, TRUE, FALSE))
  expect_identical(r$n_absorbed, 1L)
  big <- r$div[r$div$divide_id == "big", ]
  expect_equal(as.numeric(st_area(big)), 2, tolerance = 1e-6)  # 1 + orphan's 1
})

test_that("an orphan with no neighbour is KEPT, not dropped", {
  g <- st_sfc(.sq(0, 0), .sq(50, 50), crs = 5070)   # far apart
  d <- st_sf(divide_id = c("a", "orphan"),
             areasqkm = as.numeric(st_area(g)) / 1e6, geom = g)
  before <- as.numeric(sum(st_area(d)))
  r <- hydrofabric:::.absorb_orphan_divides(d, c(FALSE, TRUE))
  # Nearest-host fallback absorbs it; either way the area must survive.
  expect_equal(as.numeric(sum(st_area(r$div))), before, tolerance = 1e-6)
  expect_true(r$n_absorbed + r$n_kept == 1L)
})

test_that("all-orphan input keeps everything", {
  d <- .three()
  r <- hydrofabric:::.absorb_orphan_divides(d, c(TRUE, TRUE, TRUE))
  expect_identical(r$n_absorbed, 0L)
  expect_identical(nrow(r$div), 3L)
})

test_that("no orphans is a no-op", {
  d <- .three()
  r <- hydrofabric:::.absorb_orphan_divides(d, c(FALSE, FALSE, FALSE))
  expect_identical(r$n_absorbed, 0L)
  expect_identical(nrow(r$div), 3L)
})

test_that("the result is writable (sfc class not left stale)", {
  d <- .three()
  r <- hydrofabric:::.absorb_orphan_divides(d, c(FALSE, TRUE, FALSE))
  p <- tempfile(fileext = ".gpkg")
  expect_silent(st_write(r$div, p, layer = "divides", quiet = TRUE))
  expect_identical(nrow(st_read(p, "divides", quiet = TRUE)), 2L)
})

test_that("absorbing never adds a flowpath_id column that was absent", {
  # The coalesce caller passes a mask on a layer that has no flowpath_id at that
  # point in the pipeline. Creating one silently changes the schema, and an
  # all-NULL flowpath_id also makes .hf_network_divides treat every divide as
  # unowned, dropping the layer out of the area and coverage accounting.
  g <- st_sfc(.sq(0, 0), .sq(1, 0), crs = 5070)
  d <- st_sf(divide_id = c("a", "orphan"),
             areasqkm = as.numeric(st_area(g)) / 1e6, geom = g)
  expect_false("flowpath_id" %in% names(d))
  r <- hydrofabric:::.absorb_orphan_divides(d, c(FALSE, TRUE))
  expect_false("flowpath_id" %in% names(r$div))
})
