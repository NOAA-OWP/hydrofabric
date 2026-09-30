test_that(".pipeline_fingerprint is deterministic and param-sensitive", {
  fp1 <- hydrofabric:::.pipeline_fingerprint(list(a = 1, b = "x"))
  fp2 <- hydrofabric:::.pipeline_fingerprint(list(a = 1, b = "x"))
  fp3 <- hydrofabric:::.pipeline_fingerprint(list(a = 2, b = "x"))
  expect_type(fp1, "character")
  expect_length(fp1, 1L)
  expect_identical(fp1, fp2)            # same build + params -> stable
  expect_false(identical(fp1, fp3))     # param change -> different fingerprint
  # a no-param fingerprint still carries the build identity
  expect_true(nzchar(hydrofabric:::.pipeline_fingerprint()))
})

test_that(".checkpoints_valid gates strictly on the stamped fingerprint", {
  d <- tempfile("ckpt_fp_"); dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  fp <- "1.0||built-stamp||params-sig"

  # nothing stamped yet -> not valid (must recompute)
  expect_false(hydrofabric:::.checkpoints_valid(d, fp))

  hydrofabric:::.stamp_fingerprint(d, fp)
  expect_true(file.exists(file.path(d, ".pipeline_fingerprint")))
  expect_true(hydrofabric:::.checkpoints_valid(d, fp))          # match -> reuse ok
  expect_false(hydrofabric:::.checkpoints_valid(d, "other-fp")) # mismatch -> invalid

  # re-stamping with a new fingerprint (simulating a reinstall) invalidates the old
  hydrofabric:::.stamp_fingerprint(d, "2.0||new-build||params-sig")
  expect_false(hydrofabric:::.checkpoints_valid(d, fp))
  expect_true(hydrofabric:::.checkpoints_valid(d, "2.0||new-build||params-sig"))
})
