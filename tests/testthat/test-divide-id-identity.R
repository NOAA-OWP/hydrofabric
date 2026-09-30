# The stage-4 precondition: a cat- divide is re-IDed by looking up ITS OWN
# integer in the VPU's flowpath pool, so cat-N must name a flowpath that exists.
# This is distinct from divide_flowpath_exists, which only checks the owner.
suppressMessages({ library(sf) })

.pt <- function(n = 1L) st_sfc(rep(list(st_point(c(0, 0))), n), crs = 5070)

.inv <- function(fp, dv)
  suppressMessages(check_invariants("ngen", flowpaths = fp, divides = dv,
                                    strict = FALSE))$checks

test_that("a correctly owned but misnamed divide is caught", {
  # The exact shape that broke the global merge twice: cat-78392 owned by
  # fp-78391. The owner exists, so divide_flowpath_exists passes, but stage 4
  # looks up integer 78392, which names no flowpath.
  fp <- st_sf(flowpath_id = "fp-78391", divide_id = "cat-78392", geom = .pt())
  dv <- st_sf(divide_id = "cat-78392", flowpath_id = "fp-78391", geom = .pt())
  ch <- .inv(fp, dv)
  expect_true(ch$divide_flowpath_exists$ok)          # passes, misleadingly
  expect_false(ch$divide_id_names_a_flowpath$ok)     # the real precondition
})

test_that("a divide named for its owner passes", {
  fp <- st_sf(flowpath_id = "fp-78391", divide_id = "cat-78391", geom = .pt())
  dv <- st_sf(divide_id = "cat-78391", flowpath_id = "fp-78391", geom = .pt())
  ch <- .inv(fp, dv)
  expect_true(ch$divide_flowpath_exists$ok)
  expect_true(ch$divide_id_names_a_flowpath$ok)
})

test_that("lnd- divides are exempt: they are not keyed to a flowpath", {
  # Stage 4 pools lnd- ids separately, so an lnd- id naming no flowpath is fine.
  fp <- st_sf(flowpath_id = "fp-1", divide_id = "cat-1", geom = .pt())
  dv <- st_sf(divide_id = c("cat-1", "lnd-99999"),
              flowpath_id = c("fp-1", NA_character_), geom = .pt(2))
  ch <- .inv(fp, dv)
  expect_true(ch$divide_id_names_a_flowpath$ok)
})

test_that("an unowned cat- divide is still caught if its integer is unmatched", {
  # Releasing a divide to unowned does NOT excuse it from the id space: it is
  # still cat-prefixed, so stage 4 will still try to map its integer.
  fp <- st_sf(flowpath_id = "fp-1", divide_id = "cat-1", geom = .pt())
  dv <- st_sf(divide_id = c("cat-1", "cat-4242"),
              flowpath_id = c("fp-1", NA_character_), geom = .pt(2))
  ch <- .inv(fp, dv)
  expect_false(ch$divide_id_names_a_flowpath$ok)
})

test_that("the offending divide_id is named in the message", {
  fp <- st_sf(flowpath_id = "fp-1", divide_id = "cat-1", geom = .pt())
  dv <- st_sf(divide_id = c("cat-1", "cat-777"),
              flowpath_id = c("fp-1", "fp-1"), geom = .pt(2))
  ch <- .inv(fp, dv)
  expect_false(ch$divide_id_names_a_flowpath$ok)
  expect_match(ch$divide_id_names_a_flowpath$msg, "cat-777")
})
