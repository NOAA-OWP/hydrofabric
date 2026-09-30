test_that("flowpath_single_part counts parts, not coordinates", {
  skip_if_not(requireNamespace("sf", quietly = TRUE))
  # The regression this guards: lengths() on an sfc_LINESTRING returns the
  # coordinate count, so a plain LINESTRING would read as multi-part. Both
  # geometry types must therefore appear in the fixture.
  ls1 <- sf::st_linestring(cbind(c(0, 1, 2, 3, 4), c(0, 0, 0, 0, 0)))
  ls2 <- sf::st_linestring(cbind(c(0, 1), c(1, 1)))
  ml  <- sf::st_multilinestring(list(
    cbind(c(0, 1), c(2, 2)), cbind(c(5, 6), c(2, 2))))

  fp <- sf::st_sf(
    flowpath_id = paste0("fp-", 1:3),
    geometry    = sf::st_sfc(ls1, ls2, ml, crs = 5070))

  res <- hydrofabric:::.hf_check_ngen(list(flowpaths = fp), FALSE, 0)
  e <- res$flowpath_single_part
  expect_identical(e$kind, "info")
  # Exactly one multi-part feature, despite ls1 having 5 coordinates.
  expect_match(e$msg, "^1 of 3 flowpath\\(s\\) are multi-part")
  expect_match(e$msg, "max 2 part\\(s\\)")

  # All single-part must report zero, which is the case that the lengths() bug
  # would have reported as 2 of 2.
  fp2 <- sf::st_sf(flowpath_id = c("fp-1", "fp-2"),
    geometry = sf::st_sfc(ls1, ls2, crs = 5070))
  e2 <- hydrofabric:::.hf_check_ngen(list(flowpaths = fp2), FALSE, 0)$flowpath_single_part
  expect_match(e2$msg, "^0 of 2 flowpath\\(s\\) are multi-part")
})

test_that(".divide_metrics separates overlap from a hole, and .report_divide_effect is quiet on a no-op", {
  skip_if_not(requireNamespace("sf", quietly = TRUE))
  # Kilometre-scale squares. The reporter's thresholds are in km2 with a 1 m2
  # floor, so a fixture of unit squares would fall under the floor and read as
  # "no effect", which is correct behaviour on real data, not a defect.
  sq <- function(x0, y0, w = 1000) sf::st_polygon(list(cbind(
    c(x0, x0 + w, x0 + w, x0, x0), c(y0, y0, y0 + w, y0 + w, y0))))

  # Two disjoint 1 km2 squares: sum == union, no overlap.
  d1 <- sf::st_sf(divide_id = c("cat-1", "cat-2"), flowpath_id = c("fp-1", "fp-2"),
    geometry = sf::st_sfc(sq(0, 0), sq(2000, 0), crs = 5070))
  m1 <- hydrofabric:::.divide_metrics(d1)
  expect_identical(m1$n, 2L)
  expect_equal(m1$sum_m2, m1$union_m2, tolerance = 1e-6)
  expect_equal(m1$sum_m2, 2e6, tolerance = 1e-6)

  # Half-overlapping pair: sum exceeds union by exactly the 0.5 km2 overlap, and
  # the footprint is smaller. This is the distinction a summed area cannot make.
  d2 <- sf::st_sf(divide_id = c("cat-1", "cat-2"), flowpath_id = c("fp-1", "fp-2"),
    geometry = sf::st_sfc(sq(0, 0), sq(500, 0), crs = 5070))
  m2 <- hydrofabric:::.divide_metrics(d2)
  expect_gt(m2$sum_m2, m2$union_m2)
  expect_equal(m2$sum_m2 - m2$union_m2, 5e5, tolerance = 1e-6)

  # A no-op must be stated, not left silent: "no line" and "no effect" must not
  # look the same to whoever reads the log.
  expect_message(
    hydrofabric:::.report_divide_effect("t", m1, m1, vpuid = "01"),
    "no effect")
  # A real change reports the deltas.
  expect_message(
    hydrofabric:::.report_divide_effect("t", m2, m1, vpuid = "01"),
    "overlap 0\\.5000 -> 0\\.0000 km2")
  # Unmeasurable input says so rather than inventing a delta.
  expect_message(
    hydrofabric:::.report_divide_effect("t", hydrofabric:::.divide_metrics(NULL),
      m1, vpuid = "01"),
    "effect not measurable")
})

test_that("network_divide_has_owner catches an unowned cat- divide that no flowpath references", {
  skip_if_not(requireNamespace("sf", quietly = TRUE))
  # The hole this closes: .hf_network_divides treats a divide with a NA
  # flowpath_id as a landscape divide whatever its prefix, so an unowned cat-
  # divide leaves every area and coverage figure without any check failing.
  # divide_flowpath_consistent cannot see this case, because it walks the
  # flowpaths and nothing points at the divide.
  sq <- function(x0, w = 1000) sf::st_polygon(list(cbind(
    c(x0, x0 + w, x0 + w, x0, x0), c(0, 0, w, w, 0))))

  mk <- function(owner2) sf::st_sf(
    divide_id   = c("cat-1", "cat-2", "lnd-9"),
    flowpath_id = c("fp-1", owner2, NA_character_),
    areasqkm    = c(1, 0.25, 5),
    geometry    = sf::st_sfc(sq(0), sq(2000), sq(4000), crs = 5070))

  fp <- sf::st_sf(
    flowpath_id = "fp-1", divide_id = "cat-1",
    geometry = sf::st_sfc(sf::st_linestring(cbind(c(100, 900), c(500, 500))), crs = 5070))

  ok <- hydrofabric:::.hf_check_ngen(list(flowpaths = fp, divides = mk("fp-2")), FALSE, 0)
  expect_true(ok$network_divide_has_owner$ok)
  expect_identical(ok$network_divide_has_owner$kind, "check")

  # cat-2 unowned, and no flowpath references it.
  bad <- hydrofabric:::.hf_check_ngen(list(flowpaths = fp, divides = mk(NA_character_)), FALSE, 0)
  expect_false(bad$network_divide_has_owner$ok)
  expect_match(bad$network_divide_has_owner$msg, "^1 cat- divide\\(s\\) have no owning flowpath")
  # The excluded area is named, since that is the quantity that goes missing.
  expect_match(bad$network_divide_has_owner$msg, "0\\.2500 km2 excluded")
  # An lnd- divide with no owner is legitimate and must not trip it.
  expect_match(bad$network_divide_has_owner$msg, "cat-2")

  # And divide_flowpath_consistent stays green here, which is precisely why the
  # separate gate is needed.
  expect_true(bad$divide_flowpath_consistent$ok)
})
