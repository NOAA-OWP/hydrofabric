# A build that finishes is not a build that passed. run_vpus once recorded a VPU
# as failed only when run_pipeline threw, so a product whose gates went red was
# reported "ok" and merged. These pin the verdict reader that closes that gap.

write_invariants <- function(dir, ok, checks) {
  work <- file.path(dir, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  writeLines(jsonlite::toJSON(list(ok = ok, stage = "ngen", checks = checks),
    auto_unbox = TRUE), file.path(work, "invariants.json"))
  file.path(dir, "ngen.gpkg")
}

test_that("a clean record reports ok with nothing named", {
  g <- write_invariants(file.path(tempdir(), "inv_clean"), TRUE,
    list(a = list(ok = TRUE, msg = "", kind = "check")))
  v <- hydrofabric:::.rv_invariant_verdict(g)
  expect_true(v$ok)
  expect_length(v$failed, 0L)
})

test_that("a failed gate is caught and named", {
  g <- write_invariants(file.path(tempdir(), "inv_red"), FALSE,
    list(good = list(ok = TRUE, msg = "", kind = "check"),
         poi_id_is_integer = list(ok = FALSE, msg = "2429 non-integer",
                                  kind = "check")))
  v <- hydrofabric:::.rv_invariant_verdict(g)
  expect_false(v$ok)
  expect_identical(v$failed, "poi_id_is_integer")
})

test_that("an informational failure does not fail the VPU", {
  # Info-level checks report populations that are not yet zero; gating on them
  # would fail builds that are correct.
  g <- write_invariants(file.path(tempdir(), "inv_info"), TRUE,
    list(connector_claim_valid = list(ok = FALSE, msg = "36 of 109",
                                      kind = "info")))
  v <- hydrofabric:::.rv_invariant_verdict(g)
  expect_true(v$ok)
  expect_length(v$failed, 0L)
})

test_that("an absent or unreadable record is unverified, not a pass", {
  expect_null(hydrofabric:::.rv_invariant_verdict(NA_character_))
  expect_null(hydrofabric:::.rv_invariant_verdict(
    file.path(tempdir(), "inv_missing", "ngen.gpkg")))
})

test_that("an entry with no kind is treated as a gate, not skipped", {
  # merge_vpus writes its structural checks as bare list(ok = ...). Filtering on
  # kind == "check" skipped 12 of 21 entries in a real merge record, is_dag among
  # them, so a red structural check would have reported green.
  g <- write_invariants(file.path(tempdir(), "inv_nokind"), FALSE,
    list(is_dag = list(ok = FALSE, n = 3)))
  v <- hydrofabric:::.rv_invariant_verdict(g)
  expect_false(v$ok)
  expect_identical(v$failed, "is_dag")
})
