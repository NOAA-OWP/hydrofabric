# .absorb_divide_slivers: move a divide's sliver parts to the neighbour they
# share the most boundary with, conserving area.
#
# Fixtures are at km scale in EPSG:5070: the eligibility thresholds are absolute
# areas (km2), so a unit-square fixture would make every part look like a sliver.

.ads_sq <- function(x0, y0, w, h) sf::st_polygon(list(rbind(
  c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + h), c(x0, y0 + h), c(x0, y0))))

# A: 2 km x 2 km main body (4 km2) plus a detached 50 m x 50 m sliver (0.0025
# km2) sitting on B's edge. B: 2 km x 2 km, sharing A's right edge.
.ads_fixture <- function(sliver_w = 50, sliver_h = 50) {
  a_main   <- .ads_sq(0, 0, 2000, 2000)
  b        <- .ads_sq(2000, 0, 2000, 2000)
  a_sliver <- .ads_sq(2400, 2000, sliver_w, sliver_h)   # sits on B's top edge
  sf::st_sf(
    divide_id = c("cat-A", "cat-B"),
    geometry  = sf::st_sfc(
      sf::st_multipolygon(list(a_main, a_sliver)),
      sf::st_multipolygon(list(b)),
      crs = 5070))
}

.ads_area <- function(d) sum(as.numeric(sf::st_area(sf::st_geometry(d)))) / 1e6

test_that("a sliver moves to the neighbour it shares boundary with", {
  d <- .ads_fixture()
  before <- .ads_area(d)
  r <- hydrofabric:::.absorb_divide_slivers(d)

  expect_equal(r$n_absorbed, 1L)
  expect_equal(r$n_multi_before, 1L)
  expect_equal(r$n_multi_after, 0L)
  # cat-A is left single-part; cat-B gained the sliver
  np <- vapply(sf::st_geometry(r$div), length, integer(1))
  expect_equal(unname(np[r$div$divide_id == "cat-A"]), 1L)
  # AREA IS CONSERVED, the polygon changed owner, it was not dropped
  expect_equal(.ads_area(r$div), before, tolerance = 1e-9)
  expect_gt(as.numeric(sf::st_area(sf::st_geometry(r$div)[r$div$divide_id == "cat-B"])),
            as.numeric(sf::st_area(sf::st_geometry(d)[d$divide_id == "cat-B"])))
})

test_that("a part too large to be a sliver keeps its multi-part geometry", {
  # 1 km x 1 km second part = 1 km2, far above min_part_km2 and 20% of the
  # divide: a genuinely disjoint catchment, not a sliver.
  d <- .ads_fixture(sliver_w = 1000, sliver_h = 1000)
  before <- .ads_area(d)
  r <- hydrofabric:::.absorb_divide_slivers(d)
  expect_equal(r$n_absorbed, 0L)
  expect_equal(r$n_multi_after, 1L)
  expect_equal(.ads_area(r$div), before, tolerance = 1e-9)
})

test_that("the relative test also protects a small divide's own part", {
  # An absolutely tiny part that is nonetheless most of its divide must stay.
  tiny  <- .ads_sq(0, 0, 60, 60)          # 0.0036 km2
  tiny2 <- .ads_sq(2400, 2000, 50, 50)    # 0.0025 km2 -- 41% of the divide
  b     <- .ads_sq(2000, 0, 2000, 2000)
  d <- sf::st_sf(divide_id = c("cat-A", "cat-B"),
    geometry = sf::st_sfc(sf::st_multipolygon(list(tiny, tiny2)),
                          sf::st_multipolygon(list(b)), crs = 5070))
  r <- hydrofabric:::.absorb_divide_slivers(d)
  expect_equal(r$n_absorbed, 0L)   # under min_part_km2 but over max_part_frac
  expect_equal(r$n_multi_after, 1L)
})

test_that("a sliver with no adjacent divide is kept, never dropped", {
  a_main   <- .ads_sq(0, 0, 2000, 2000)
  a_island <- .ads_sq(9000, 9000, 50, 50)   # nothing near it
  b        <- .ads_sq(2000, 0, 2000, 2000)
  d <- sf::st_sf(divide_id = c("cat-A", "cat-B"),
    geometry = sf::st_sfc(sf::st_multipolygon(list(a_main, a_island)),
                          sf::st_multipolygon(list(b)), crs = 5070))
  before <- .ads_area(d)
  r <- hydrofabric:::.absorb_divide_slivers(d)
  expect_equal(r$n_absorbed, 0L)
  expect_equal(r$n_kept, 1L)                       # eligible but unplaceable
  expect_equal(r$n_multi_after, 1L)
  expect_equal(.ads_area(r$div), before, tolerance = 1e-9)
})

test_that("a single-part layer is returned untouched", {
  d <- sf::st_sf(divide_id = c("cat-A", "cat-B"),
    geometry = sf::st_sfc(sf::st_multipolygon(list(.ads_sq(0, 0, 2000, 2000))),
                          sf::st_multipolygon(list(.ads_sq(2000, 0, 2000, 2000))),
                          crs = 5070))
  r <- hydrofabric:::.absorb_divide_slivers(d)
  expect_equal(r$n_absorbed, 0L)
  expect_equal(r$n_multi_before, 0L)
  expect_identical(sf::st_geometry(r$div), sf::st_geometry(d))
})

test_that("the recipient is the neighbour sharing the MOST boundary", {
  # The sliver must not touch its own main body: parts that touch are one
  # polygon and the seam pass dissolves them before any absorption is considered.
  # It shares a 60 m edge with B and a 40 m edge with C, so B should win.
  a_main   <- .ads_sq(0, 0, 2000, 2000)
  a_sliver <- .ads_sq(2100, 2000, 60, 40)       # 100 m clear of a_main's corner
  b        <- .ads_sq(2000, 1000, 2000, 1000)   # top edge carries the sliver's 60 m base
  cc       <- .ads_sq(2160, 2000, 500, 500)     # shares only the 40 m right side
  d <- sf::st_sf(divide_id = c("cat-A", "cat-B", "cat-C"),
    geometry = sf::st_sfc(sf::st_multipolygon(list(a_main, a_sliver)),
                          sf::st_multipolygon(list(b)),
                          sf::st_multipolygon(list(cc)), crs = 5070))
  r <- hydrofabric:::.absorb_divide_slivers(d)
  expect_equal(r$n_absorbed, 1L)
  gain <- function(id) as.numeric(sf::st_area(sf::st_geometry(r$div)[r$div$divide_id == id])) -
                       as.numeric(sf::st_area(sf::st_geometry(d)[d$divide_id == id]))
  expect_gt(gain("cat-B"), 0)      # 60 m edge wins
  expect_equal(gain("cat-C"), 0, tolerance = 1e-6)
})

test_that("parts that merely touch are dissolved, not treated as separate", {
  # Two squares meeting at a single corner point are one polygon, but st_union
  # leaves them as two parts because they share no edge. Half the real
  # multi-part divides are exactly this. A 1 cm close merges them; area must not
  # move, and the part must NOT be handed to a neighbour.
  a1 <- .ads_sq(0, 0, 1000, 1000)
  a2 <- .ads_sq(1000, 1000, 400, 400)           # touches a1 at (1000,1000) only
  b  <- .ads_sq(2000, 0, 1000, 1000)            # far away
  d <- sf::st_sf(divide_id = c("cat-A", "cat-B"),
    geometry = sf::st_sfc(sf::st_multipolygon(list(a1, a2)),
                          sf::st_multipolygon(list(b)), crs = 5070))
  before <- .ads_area(d)
  r <- hydrofabric:::.absorb_divide_slivers(d)
  expect_equal(r$n_seam_dissolved, 1L)
  expect_equal(r$n_multi_after, 0L)
  expect_equal(r$n_absorbed, 0L)                # dissolved in place, not moved
  expect_equal(.ads_area(r$div), before, tolerance = 1e-6)
  expect_equal(as.numeric(sf::st_area(sf::st_geometry(r$div)[r$div$divide_id == "cat-B"])),
               1000 * 1000, tolerance = 1e-6)   # B unchanged
})

test_that("a real separation is never bridged by the seam pass", {
  # 30 m apart, one DEM cell. Must stay two parts.
  a1 <- .ads_sq(0, 0, 1000, 1000)
  a2 <- .ads_sq(1030, 0, 1000, 1000)
  d <- sf::st_sf(divide_id = "cat-A",
    geometry = sf::st_sfc(sf::st_multipolygon(list(a1, a2)), crs = 5070))
  r <- hydrofabric:::.absorb_divide_slivers(d)
  expect_equal(r$n_seam_dissolved, 0L)
  expect_equal(r$n_multi_after, 1L)
})

test_that("a part carrying its own flowpath's flowline is never absorbed", {
  # The AK gate failure. All five failing divides lost exactly the part their
  # own flowline ran through, leaving them with no member flowline and failing
  # divide_has_member_flowline. Such a part is the head of the catchment, not a
  # sliver, however small it is.
  d <- .ads_fixture()                       # cat-A: main body + detached sliver
  d$flowpath_id <- c("fp-A", "fp-B")
  # fp-A's flowline runs through the sliver, not the main body
  fl <- sf::st_sf(
    flowpath_id = "fp-A",
    geometry = sf::st_sfc(sf::st_linestring(rbind(c(2410, 2010), c(2440, 2040))),
                          crs = 5070))
  before <- .ads_area(d)

  r <- hydrofabric:::.absorb_divide_slivers(d, flowlines = fl)
  expect_equal(r$n_absorbed, 0L)            # refused
  expect_equal(.ads_area(r$div), before, tolerance = 1e-9)

  # and without the flowline argument the old behaviour is unchanged
  r2 <- hydrofabric:::.absorb_divide_slivers(d)
  expect_equal(r2$n_absorbed, 1L)
})

test_that("a sliver with no flowline of its own is still absorbed", {
  # The guard must not disable absorption generally: a part the reach never
  # touches is exactly what this pass exists to move.
  d <- .ads_fixture()
  d$flowpath_id <- c("fp-A", "fp-B")
  fl <- sf::st_sf(
    flowpath_id = "fp-A",
    geometry = sf::st_sfc(sf::st_linestring(rbind(c(100, 100), c(900, 900))),
                          crs = 5070))     # runs through the main body
  r <- hydrofabric:::.absorb_divide_slivers(d, flowlines = fl)
  expect_equal(r$n_absorbed, 1L)
})
