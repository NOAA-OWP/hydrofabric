# Late multi-part divide cleanup (build_nextgen Step 6n3).
#
# Step 6d5b absorbs slivers before the steps that reshape divides run, so
# anything those steps create survived to the product: on VPU 01 six divides
# reached the delivered layer multi-part, three with a real sliver and three
# with a part of no area at all.

test_that("a degenerate part carries no area, so dropping one conserves area", {
  sq <- function(x, y, s) sf::st_polygon(list(rbind(
    c(x, y), c(x + s, y), c(x + s, y + s), c(x, y + s), c(x, y))))
  main <- sq(0, 0, 1000)                    # 1 km2
  tiny <- sq(5000, 5000, 0.5)               # 0.25 m2, detached
  g <- sf::st_sfc(sf::st_multipolygon(list(main, tiny)), crs = 5070)
  ar <- vapply(g[[1]], function(pp)
    as.numeric(sf::st_area(sf::st_sfc(sf::st_polygon(pp), crs = 5070))), numeric(1))
  kept <- sf::st_sfc(sf::st_multipolygon(g[[1]][ar >= 1]), crs = 5070)
  expect_equal(length(kept[[1]]), 1L)
  expect_equal(as.numeric(sf::st_area(kept)), as.numeric(sf::st_area(g)),
               tolerance = 1e-6)
})

test_that("a sliver sharing a boundary is moved, not dropped", {
  # .absorb_divide_slivers is area-conserving by construction: the polygon
  # changes owner. That is the difference from dissolve-and-drop, which would
  # forfeit the land.
  sq <- function(x, y, w, h) sf::st_polygon(list(rbind(
    c(x, y), c(x + w, y), c(x + w, y + h), c(x, y + h), c(x, y))))
  a_main <- sq(0, 0, 1000, 1000)
  a_sliv <- sq(980, 1600, 20, 20)           # 400 m2, detached from a_main
  b_main <- sq(1000, 1500, 1000, 1000)      # abuts the sliver along x = 1000
  d <- sf::st_sf(divide_id = c("cat-a", "cat-b"),
    geometry = sf::st_sfc(sf::st_multipolygon(list(a_main, a_sliv)),
                          sf::st_multipolygon(list(b_main)), crs = 5070))
  before <- sum(as.numeric(sf::st_area(d)))
  r <- hydrofabric:::.absorb_divide_slivers(d)
  after <- sum(as.numeric(sf::st_area(r$div)))
  expect_equal(after, before, tolerance = 1e-6)
  expect_equal(r$n_multi_after, 0L)
  expect_equal(r$n_absorbed, 1L)
})

test_that("a sliver moves when the keeper part also holds a flowline", {
  # The guard exists so absorbing a part cannot leave a divide with no member
  # flowline. When the keeper holds one too, the divide keeps one either way.
  sq <- function(x, y, w, h) sf::st_polygon(list(rbind(
    c(x, y), c(x + w, y), c(x + w, y + h), c(x, y + h), c(x, y))))
  a_main <- sq(0, 0, 1000, 1000)
  a_sliv <- sq(980, 1600, 20, 20)
  b_main <- sq(1000, 1500, 1000, 1000)
  d <- sf::st_sf(divide_id = c("cat-a", "cat-b"),
    flowpath_id = c("fp-a", "fp-b"),
    geometry = sf::st_sfc(sf::st_multipolygon(list(a_main, a_sliv)),
                          sf::st_multipolygon(list(b_main)), crs = 5070))
  # one line through the keeper, one through the sliver
  fl <- sf::st_sf(flowpath_id = c("fp-a", "fp-a"),
    geometry = sf::st_sfc(
      sf::st_linestring(rbind(c(100, 100), c(900, 900))),      # keeper
      sf::st_linestring(rbind(c(985, 1605), c(995, 1615))),    # sliver
      crs = 5070))
  r <- hydrofabric:::.absorb_divide_slivers(d, flowlines = fl)
  expect_equal(r$n_absorbed, 1L)
  expect_equal(r$n_multi_after, 0L)
})

test_that("a sliver is protected when it is the only part holding a flowline", {
  # This is the Alaska case: absorbing the part would leave the divide with no
  # member flowline at all, so it stays multi-part instead.
  sq <- function(x, y, w, h) sf::st_polygon(list(rbind(
    c(x, y), c(x + w, y), c(x + w, y + h), c(x, y + h), c(x, y))))
  d <- sf::st_sf(divide_id = c("cat-a", "cat-b"),
    flowpath_id = c("fp-a", "fp-b"),
    geometry = sf::st_sfc(
      sf::st_multipolygon(list(sq(0, 0, 1000, 1000), sq(980, 1600, 20, 20))),
      sf::st_multipolygon(list(sq(1000, 1500, 1000, 1000))), crs = 5070))
  fl <- sf::st_sf(flowpath_id = "fp-a", geometry = sf::st_sfc(
    sf::st_linestring(rbind(c(985, 1605), c(995, 1615))), crs = 5070))  # sliver only
  r <- hydrofabric:::.absorb_divide_slivers(d, flowlines = fl)
  expect_equal(r$n_absorbed, 0L)
  expect_equal(r$n_multi_after, 1L)
})
