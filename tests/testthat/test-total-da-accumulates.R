# total_dasqkm is built by accumulating each flowline's own divide area down the
# flowline graph, so `total = incremental + sum(upstream totals)` holds exactly
# at every reach. That identity is what makes the figure checkable without an
# external baseline, the reference `flowline_total_dasqkm` cannot serve as
# one, being a reference attribute carried onto reconciled reaches, where a
# split fragment inherits its parent's value.

fl3 <- function(tot) data.frame(
  flowline_id   = c(1, 2, 3),
  flowline_toid = c(3, 3, 0),          # 1 and 2 drain into 3
  incremental_areasqkm = c(2, 4, 1),
  total_dasqkm = tot, stringsAsFactors = FALSE)

test_that("a consistent accumulation passes", {
  r <- hydrofabric:::.hf_check_ngen(list(flowlines = fl3(c(2, 4, 7))), FALSE, 0)
  expect_true(r$total_da_accumulates$ok)
})

test_that("a reach whose total omits an upstream contributor fails", {
  # 3 should be 1 + 2 + 4 = 7; 5 drops one tributary.
  r <- hydrofabric:::.hf_check_ngen(list(flowlines = fl3(c(2, 4, 5))), FALSE, 0)
  expect_false(r$total_da_accumulates$ok)
  expect_match(r$total_da_accumulates$msg, "1 flowline")
})

test_that("the check is skipped when the columns are absent", {
  r <- hydrofabric:::.hf_check_ngen(list(flowlines = data.frame(flowline_id = 1)),
                                    FALSE, 0)
  expect_null(r$total_da_accumulates)
})
