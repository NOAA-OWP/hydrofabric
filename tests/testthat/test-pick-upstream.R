# A POI at an inlet belongs to the upstream reach continuing the same mainstem.
# Two ways the earlier selection got that wrong:
#
#   a candidate with no levelpathid made `== lp` yield NA, and NA-subscripting
#   injected empty rows into the match set;
#
#   where no upstream shared the mainstem, the choice fell to list order, which
#   at a confluence hands the POI a small tributary instead of the channel it is
#   on.

fp_tbl <- function(levelpath, da) {
  data.frame(flowpath_id = paste0("fp-", seq_along(levelpath)),
             levelpathid = levelpath,
             tot_drainage_area = da,
             stringsAsFactors = FALSE)
}

test_that("the upstream on the same mainstem wins, whatever the order", {
  # row 4 is downstream (levelpath 100); candidates are rows 1-3
  d <- fp_tbl(c(200, 100, 300, 100), c(500, 12, 900, 1000))
  # the big tributary (row 3, 900) must NOT win over the mainstem (row 2)
  expect_identical(hydrofabric:::.pick_upstream(c(1L, 3L, 2L), 4L, d), 2L)
})

test_that("a candidate with no levelpathid does not poison the match", {
  d <- fp_tbl(c(NA, 100, 300, 100), c(500, 12, 900, 1000))
  expect_identical(hydrofabric:::.pick_upstream(c(1L, 2L, 3L), 4L, d), 2L)
})

test_that("with no mainstem match, the largest drainage area wins, not the first", {
  # downstream is levelpath 100; no candidate shares it
  d <- fp_tbl(c(700, 800, 900, 100), c(4.3, 902.9, 15, 903))
  expect_identical(hydrofabric:::.pick_upstream(c(1L, 2L, 3L), 4L, d), 2L)
})

test_that("a single candidate is returned regardless", {
  d <- fp_tbl(c(700, 100), c(4, 900))
  expect_identical(hydrofabric:::.pick_upstream(1L, 2L, d), 1L)
})

test_that("no candidates yields NA rather than an error", {
  d <- fp_tbl(c(100, 100), c(1, 2))
  expect_true(is.na(hydrofabric:::.pick_upstream(integer(0), 2L, d)))
})

test_that("a downstream reach with no levelpathid falls back on drainage area", {
  d <- fp_tbl(c(700, 800, NA), c(4, 900, 1000))
  expect_identical(hydrofabric:::.pick_upstream(c(1L, 2L), 3L, d), 2L)
})

test_that("missing drainage areas do not select an unusable candidate", {
  d <- fp_tbl(c(700, 800, 100), c(NA, 50, 900))
  expect_identical(hydrofabric:::.pick_upstream(c(1L, 2L), 3L, d), 2L)
})
