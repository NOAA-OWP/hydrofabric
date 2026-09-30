# flowpath_coverage measures an aggregated flowpath's line against its own
# catchment, which is not a defect when it falls short: the line is a mainstem
# walk and need not lie inside the catchment. It also only ever passed because
# Step 6j cleared the offenders and Step 6k deleted them. It is now INFO, and
# divide_has_member_flowline is the gate that carries meaning.
suppressMessages({ library(sf) })

# Kilometre scale, in metres (EPSG:5070). Two reasons it cannot be unit-scale:
# divide_has_member_flowline accepts a flowline within the fabric's 30 m tolerance,
# so a unit-square fixture puts a "far" line ~7 m away and inside tolerance; and
# the rows declare areasqkm = 1, which only matches a 1000 m square.
.KM <- 1000
.sq <- function(x0) { x <- x0 * .KM
  st_polygon(list(cbind(c(x, x + .KM, x + .KM, x, x),
                        c(0, 0, .KM, .KM, 0)))) }
.ln <- function(x0) { x <- x0 * .KM
  st_linestring(cbind(c(x + 200, x + 800), c(500, 500))) }

.fab <- function(fl_x) {
  n <- length(fl_x)
  ids <- paste0("fp-", seq_len(n))
  list(
    fp = st_sf(flowpath_id = ids, divide_id = paste0("cat-", seq_len(n)),
               geom = st_sfc(lapply(fl_x, .ln), crs = 5070)),
    dv = st_sf(divide_id = paste0("cat-", seq_len(n)), flowpath_id = ids,
               areasqkm = 1,
               geom = st_sfc(lapply(seq_len(n) - 1, .sq), crs = 5070)),
    fl = st_sf(flowpath_id = ids,
               geom = st_sfc(lapply(fl_x, .ln), crs = 5070)))
}

.inv <- function(f, ...) suppressMessages(
  check_invariants("ngen", flowpaths = f$fp, divides = f$dv, flowlines = f$fl,
                   strict = FALSE, ...))

test_that("flowpath_coverage is INFO, so a coarse aggregate cannot fail a build", {
  f <- .fab(c(0, 9))            # fp-2's line is far from cat-2
  r <- .inv(f, coverage = TRUE)
  expect_identical(r$checks$flowpath_coverage$kind, "info")
})

test_that("a divide containing none of its flowpath's flowlines is a failure", {
  f <- .fab(c(0, 9))
  r <- .inv(f)
  expect_false(r$checks$divide_has_member_flowline$ok)
  expect_match(r$checks$divide_has_member_flowline$msg, "cat-2")
})

test_that("a divide containing its flowpath's flowline passes", {
  f <- .fab(c(0, 1))            # both lines sit inside their own square
  r <- .inv(f)
  expect_true(r$checks$divide_has_member_flowline$ok)
})

test_that("a divide whose flowpath has no flowlines is not judged", {
  # Untestable rather than failing: absence of flowlines is a different concern,
  # already covered by flowpath_has_constituent.
  f <- .fab(c(0, 1))
  f$fl <- f$fl[f$fl$flowpath_id == "fp-1", ]
  r <- .inv(f)
  expect_true(r$checks$divide_has_member_flowline$ok)
  expect_match(r$checks$divide_has_member_flowline$msg, "of 1 network divide")
})

test_that("unowned and lnd- divides are excluded from the gate", {
  f <- .fab(c(0, 1))
  f$dv$divide_id[2]   <- "lnd-2"
  f$dv$flowpath_id[2] <- NA_character_
  r <- .inv(f)
  expect_true(r$checks$divide_has_member_flowline$ok)
  expect_match(r$checks$divide_has_member_flowline$msg, "of 1 network divide")
})
