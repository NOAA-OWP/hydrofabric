# A source link is distrusted when the point sits further than max_link_m from
# the reach it names, and the nearest flowpath is used instead. Distance alone is
# the wrong test: a location named on a large river whose nearest feature drains a
# tiny fraction of it is a coordinate that missed the channel, not a bad link.
#
# Measured on the CONUS hydrolocation set, 1,800 rows have a coordinate and a link
# that disagree; 24 gages ended up on a reach draining under 1% of the one they
# named, 17 of those on rivers over 10,000 km2 while sitting under 1 km away.
#
# The rescue restores the link in exactly that case. It never moves the point.

decide <- function(named_da, near_da, ratio = 100, min_da = 100) {
  # TRUE = keep the source link despite the distance
  is.finite(named_da) && is.finite(near_da) &&
    named_da >= min_da && near_da < named_da / ratio
}

test_that("a large river is not abandoned for an adjacent minor feature", {
  # 05587450: named 444,031 km2, nearest drains 1.41 km2
  expect_true(decide(444030.77, 1.41))
  # 07138020: named 64,317 km2, nearest drains 0.12 km2
  expect_true(decide(64316.89, 0.12))
})

test_that("a genuinely mis-stated link on a small reach still falls back", {
  # named reach is small, so the comparison carries no information
  expect_false(decide(12.0, 0.4))
  expect_false(decide(99.9, 0.1))
})

test_that("a plausible nearest host still wins", {
  # nearest drains a comparable area: the coordinate is believable
  expect_false(decide(5000, 4800))
  expect_false(decide(5000, 60))      # 1.2%, above the 1% floor
})

test_that("the boundary is where the ratio says it is", {
  expect_false(decide(10000, 100))    # exactly 1% -- not below
  expect_true(decide(10000, 99.9))    # just under
})

test_that("missing drainage areas never trigger the rescue", {
  expect_false(decide(NA_real_, 1))
  expect_false(decide(5000, NA_real_))
  expect_false(decide(Inf, 1))
})

test_that("the thresholds are adjustable", {
  # a caller demanding a stricter ratio rescues fewer rows
  expect_true(decide(5000, 40, ratio = 100))
  expect_false(decide(5000, 40, ratio = 1000))
})

test_that("the rescue is wired into prepare_hydrolocations", {
  a <- names(formals(hydrofabric::prepare_hydrolocations))
  expect_true(all(c("link_keep_da_ratio", "link_keep_min_da_sqkm") %in% a))
  src <- paste(deparse(hydrofabric::prepare_hydrolocations), collapse = "\n")
  expect_true(grepl("link_keep_da_ratio", src, fixed = TRUE))
  # and it must adjust link_far, not the geometry
  expect_true(grepl("link_far[lf[keep]] <- FALSE", src, fixed = TRUE))
})
