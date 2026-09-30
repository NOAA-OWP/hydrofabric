# aggregated_divide_link_resolves
#
# A flowpath claiming a divide must name one that exists. Where it does not,
# there is no polygon to place; the nextgen phase clears the unresolvable id and
# the flowpath ships as a real river owning no land. Asserted at aggregate grain because the finished product is
# always clean: Step 6d5 nulls the dangling id after the loss has happened.

.adl <- function(fp, div)
  hydrofabric:::.hf_check_aggregated(
    list(flowpaths = fp, divides = div), FALSE, 0)$aggregated_divide_link_resolves

.adl_fp <- function(divide_id, has_divide = TRUE)
  data.frame(flowpath_id = paste0("fp-", seq_along(divide_id)),
             divide_id = divide_id, has_divide = has_divide,
             stringsAsFactors = FALSE)

test_that("a flowpath naming a divide that does not exist is reported", {
  fp  <- .adl_fp(c("cat-1", "cat-2"))
  div <- data.frame(divide_id = "cat-1", stringsAsFactors = FALSE)
  chk <- .adl(fp, div)
  expect_match(chk$msg, "^1 of 2 flowpath")
  expect_match(chk$msg, "fp-2")
  expect_identical(chk$kind, "info")
})

test_that("all links resolving passes as a gate", {
  fp  <- .adl_fp(c("cat-1", "cat-2"))
  div <- data.frame(divide_id = c("cat-1", "cat-2"), stringsAsFactors = FALSE)
  chk <- .adl(fp, div)
  expect_true(chk$ok)
  expect_identical(chk$kind, "check")
  expect_match(chk$msg, "^0 of 2 flowpath")
})

test_that("a flowpath that claims no divide is not judged", {
  # Legitimately divide-less connectors must not be counted as dangling.
  fp  <- .adl_fp(c("cat-1", NA), has_divide = c(TRUE, FALSE))
  div <- data.frame(divide_id = "cat-1", stringsAsFactors = FALSE)
  chk <- .adl(fp, div)
  expect_match(chk$msg, "^0 of 1 flowpath")
  expect_true(chk$ok)
})

test_that("has_divide FALSE with a leftover divide_id is not judged", {
  # The flowpath is not claiming the divide, whatever the stale id says.
  fp  <- .adl_fp(c("cat-1", "cat-99"), has_divide = c(TRUE, FALSE))
  div <- data.frame(divide_id = "cat-1", stringsAsFactors = FALSE)
  chk <- .adl(fp, div)
  expect_match(chk$msg, "^0 of 1 flowpath")
})

test_that("the check is skipped when the columns are absent", {
  fp  <- data.frame(flowpath_id = "fp-1", stringsAsFactors = FALSE)
  div <- data.frame(divide_id = "cat-1", stringsAsFactors = FALSE)
  expect_null(.adl(fp, div))
})
