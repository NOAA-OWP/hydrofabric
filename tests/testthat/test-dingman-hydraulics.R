# Dingman's cross-section supplies the bf_area the trapezoid inversion needs.
# r = 1 is a triangle (A = WD/2, zero bottom width); large r approaches a
# rectangle (A = WD, bottom width -> top width). The RiverML path must be left
# untouched, since both sources feed the same formulas.

test_that("bankfull area follows the shape exponent to its limits", {
  d <- data.frame(flowpath_id = 1:3, owp_tw_bf_m = c(10, 10, 10),
                  owp_y_bf_m = c(2, 2, 2), owp_r_bf = c(1, 3, 1e6),
                  n_single = 0.05, slope_m_m = 0.01)
  out <- hydrofabric:::.as_canonical_hydraulics(d)
  W <- 10; D <- 2
  expect_equal(out$bf_area[1], W * D / 2)              # triangle
  expect_equal(out$bf_area[2], W * D * 3 / 4)          # between
  expect_equal(out$bf_area[3], W * D, tolerance = 1e-5) # rectangle in the limit
})

test_that("the derived bottom width spans zero to the full top width", {
  d <- data.frame(flowpath_id = 1:2, owp_tw_bf_m = 10, owp_y_bf_m = 2,
                  owp_r_bf = c(1, 1e6), n_single = 0.05, slope_m_m = 0.01)
  o <- hydrofabric:::.as_canonical_hydraulics(d)
  Bw <- (2 * o$bf_area / o$owp_y_bf) - o$owp_tw_bf
  expect_equal(Bw[1], 0)                        # triangle has no bed
  expect_equal(Bw[2], 10, tolerance = 1e-4)     # rectangle is all bed
})

test_that("columns are renamed to the canonical names", {
  d <- data.frame(flowpath_id = 1L, owp_tw_bf_m = 8, owp_y_bf_m = 1.5,
                  owp_r_bf = 3, n_single = 0.04, slope_m_m = 0.002)
  o <- hydrofabric:::.as_canonical_hydraulics(d)
  expect_true(all(c("feature_id","owp_tw_bf","owp_y_bf","bf_area",
                    "owp_roughness_bathy","slope","lengthkm") %in% names(o)))
  expect_identical(o$owp_roughness_bathy, 0.04)
  expect_identical(o$slope, 0.002)
})

test_that("supplied reach lengths become the weights", {
  d <- data.frame(flowpath_id = c(1L, 2L), owp_tw_bf_m = 8, owp_y_bf_m = 1.5,
                  owp_r_bf = 3, n_single = 0.04, slope_m_m = 0.002)
  len <- data.frame(reference_id = c(1L, 2L), lengthkm = c(3, 7))
  o <- hydrofabric:::.as_canonical_hydraulics(d, len)
  expect_identical(o$lengthkm, c(3, 7))
})

test_that("a RiverML table passes through untouched", {
  d <- data.frame(feature_id = 1L, owp_y_bf = 2, owp_tw_bf = 10, bf_area = 15,
                  owp_roughness_bathy = 0.05, slope = 0.01, lengthkm = 1)
  expect_false(hydrofabric:::.is_dingman_hydraulics(d))
  expect_identical(hydrofabric:::.as_canonical_hydraulics(d), d)
})

# Hawaii and Puerto Rico carry 14-digit NHD-HR reference ids. Matching the
# hydraulics table to the network on as.integer() ids turns every one of them
# into NA, and the domain silently takes the geomorphic proxy instead of the
# table describing it.

test_that("reaches above the 32-bit ceiling still find their hydraulics", {
  big <- c("80000100004003", "80000100003990")
  network <- data.frame(
    reference_id     = big,
    flowpath_id      = c("fp-1", "fp-1"),
    flowline_id      = c(1L, 2L),
    flowline_lengthkm = c(2, 3),
    lengthkm         = c(5, 5),
    hl_reference     = NA_character_,
    total_dasqkm     = c(10, 10),
    slope            = c(0.01, 0.01),
    stringsAsFactors = FALSE)
  hyd <- data.frame(flowpath_id = bit64::as.integer64(big),
                    owp_tw_bf_m = c(10, 20), owp_y_bf_m = c(2, 2),
                    owp_r_bf = 3, n_single = 0.05, slope_m_m = c(0.02, 0.04))

  out <- hydrofabric:::.calculate_flow_atts(
    network = network, hydraulics = hyd, hf_id = "flowpath_id")

  expect_identical(out$attr_source, "riverml")
  expect_false(any(is.na(out$BtmWdth)))
  expect_equal(out$So, 0.03)   # equal weights: mean of 0.02 and 0.04
})

test_that("reach lengths weight the reaches behind an aggregate", {
  big <- c("80000100004003", "80000100003990")
  network <- data.frame(
    reference_id     = big,
    flowpath_id      = c("fp-1", "fp-1"),
    flowline_lengthkm = c(1, 9),
    lengthkm         = c(10, 10),
    hl_reference     = NA_character_,
    total_dasqkm     = c(10, 10),
    slope            = c(0.01, 0.01),
    stringsAsFactors = FALSE)
  hyd <- data.frame(flowpath_id = bit64::as.integer64(big),
                    owp_tw_bf_m = c(10, 20), owp_y_bf_m = c(2, 2),
                    owp_r_bf = 3, n_single = 0.05, slope_m_m = c(0.02, 0.04))
  len <- data.frame(reference_id = big, lengthkm = c(1, 9))

  out <- hydrofabric:::.calculate_flow_atts(
    network = network, hydraulics = hyd, hf_id = "flowpath_id",
    reach_lengths = len)

  # 0.02 * 1/10 + 0.04 * 9/10
  expect_equal(out$So, 0.038)
})

# The OWP channel tables ship under two column spellings. The ML v2 tables drop
# the `owp_` prefix from the geometry columns, and a table in either spelling
# must canonicalize to the same thing rather than being passed through unread.

test_that("both column spellings of a Dingman table are recognized", {
  v1 <- data.frame(flowpath_id = 1:2, owp_tw_bf_m = c(10, 20), owp_y_bf_m = c(2, 2),
                   owp_r_bf = 3, n_single = 0.05, slope_m_m = c(0.02, 0.04))
  v2 <- data.frame(flowpath_id = 1:2, tw_bf_m = c(10, 20), y_bf_m = c(2, 2),
                   r_bf = 3, n_single = 0.05, slope_m_m = c(0.02, 0.04),
                   n_in_channel = 0.04, n_out_channel = 0.07,
                   n_single_confidence_score = 80)
  expect_true(hydrofabric:::.is_dingman_hydraulics(v1))
  expect_true(hydrofabric:::.is_dingman_hydraulics(v2))

  a <- hydrofabric:::.as_canonical_hydraulics(v1)
  b <- hydrofabric:::.as_canonical_hydraulics(v2)
  expect_equal(a$owp_tw_bf, b$owp_tw_bf)
  expect_equal(a$bf_area,   b$bf_area)
  expect_equal(a$slope,     b$slope)
  # The geometry canonicalizes identically whichever spelling it arrived in.
  expect_identical(setdiff(names(a), names(b)), character(0))
  # v2 carries an overbank roughness, which the v1 fixture has no column for;
  # the confidence scores are dropped either way.
  expect_identical(setdiff(names(b), names(a)), "owp_roughness_overbank")
  expect_false(any(grepl("confidence", names(b))))
})

test_that("a table in neither spelling is passed through untouched", {
  riverml <- data.frame(feature_id = 1:2, owp_y_bf = 1, owp_tw_bf = 5,
                        bf_area = 4, owp_roughness_bathy = 0.05, slope = 0.01,
                        lengthkm = 1)
  expect_false(hydrofabric:::.is_dingman_hydraulics(riverml))
  expect_identical(hydrofabric:::.as_canonical_hydraulics(riverml), riverml)
})

# Channel routing wants the in-channel roughness, and the compound-channel
# roughness is what the overbank estimate describes. Doubling the channel value
# is the stand-in for a source that supplies only one.

test_that("roughness comes from the in-channel column and nCC from the overbank", {
  network <- data.frame(
    reference_id = c("11", "12"), flowpath_id = c("fp-1", "fp-1"),
    flowline_lengthkm = c(1, 1), lengthkm = c(2, 2),
    hl_reference = NA_character_, total_dasqkm = c(5, 5), slope = c(0.01, 0.01),
    stringsAsFactors = FALSE)
  hyd <- data.frame(flowpath_id = c(11L, 12L), tw_bf_m = 10, y_bf_m = 2,
                    r_bf = 3, n_single = 0.20, n_in_channel = 0.05,
                    n_out_channel = 0.09, slope_m_m = 0.01)

  out <- hydrofabric:::.calculate_flow_atts(
    network = network, hydraulics = hyd, hf_id = "flowpath_id")

  expect_equal(out$n,   0.05)   # in-channel, not n_single
  expect_equal(out$nCC, 0.09)   # overbank, not 2 * n
})

test_that("a source with one roughness still gets the doubled nCC", {
  network <- data.frame(
    reference_id = "11", flowpath_id = "fp-1", flowline_lengthkm = 1,
    lengthkm = 1, hl_reference = NA_character_, total_dasqkm = 5, slope = 0.01,
    stringsAsFactors = FALSE)
  hyd <- data.frame(flowpath_id = 11L, tw_bf_m = 10, y_bf_m = 2, r_bf = 3,
                    n_single = 0.06, slope_m_m = 0.01)

  out <- hydrofabric:::.calculate_flow_atts(
    network = network, hydraulics = hyd, hf_id = "flowpath_id")

  expect_equal(out$n,   0.06)        # falls back to n_single
  expect_equal(out$nCC, 2 * 0.06)    # no overbank column, so the heuristic
})
