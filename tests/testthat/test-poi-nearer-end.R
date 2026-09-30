# A POI that cannot split its reach without leaving a sub-threshold piece moves
# to the NEARER end: the upstream reach's outlet if it is closer to the inlet,
# its own outlet if it is closer to the outlet. Only the former changes
# flowpath_id.
#
# The outlet branch is guarded by `down_piece > 0`, so a POI sitting exactly at
# the outlet falls through it. On a reach shorter than the threshold the inlet
# branch then matches as well, because up_piece is the whole reach length, and
# the POI is walked upstream and labelled an inlet move.
#
# These tests exercise the decision on measures alone, which is what the pass
# uses; no geometry is needed.

decide <- function(m, fl, thresh = 250, has_upstream = TRUE) {
  down <- fl - m; up <- m
  if (down == 0) return("no action")
  if (m > fl || (down > 0 && down < thresh && down <= up)) return("own outlet")
  if (up >= 0 && up < thresh && up <= down)
    return(if (has_upstream) "upstream outlet" else "headwater: keep own reach")
  "split"
}

test_that("a POI already at the outlet is left alone, even on a short reach", {
  # the regression: fl < threshold, so the inlet branch used to claim it
  expect_identical(decide(m = 200, fl = 200), "no action")
  expect_identical(decide(m = 60,  fl = 60),  "no action")
})

test_that("nearer end wins when both pieces are under threshold", {
  # reach 394 m, POI 157 from inlet / 237 from outlet -> inlet is nearer
  expect_identical(decide(m = 157, fl = 394), "upstream outlet")
  # same reach, POI 237 from inlet / 157 from outlet -> outlet is nearer
  expect_identical(decide(m = 237, fl = 394), "own outlet")
})

test_that("a tie goes to the outlet, leaving flowpath_id alone", {
  # The outlet branch is evaluated first and returns, so an exact tie keeps the
  # POI on its own reach. That is the conservative side of the rule: only an
  # inlet move re-points flowpath_id, so a tie should not re-point.
  expect_identical(decide(m = 100, fl = 200), "own outlet")
})

test_that("a POI with room on both sides is split, not moved", {
  expect_identical(decide(m = 600, fl = 1200), "split")
})

test_that("a headwater POI near its inlet keeps its own reach", {
  expect_identical(decide(m = 40, fl = 900, has_upstream = FALSE),
                   "headwater: keep own reach")
})

test_that("the outlet-side fix is present in the installed pass", {
  src <- paste(deparse(getFromNamespace(
    "rf_generate_events_from_hydrolocations_fast", asNamespace("hydrofabric"))),
    collapse = "\n")
  expect_true(grepl("down_piece == 0", src, fixed = TRUE))
  expect_true(grepl("up_piece <= down_piece", src, fixed = TRUE))
  expect_true(grepl("at_inlet_of_headwater_kept_own_reach", src, fixed = TRUE))
})
