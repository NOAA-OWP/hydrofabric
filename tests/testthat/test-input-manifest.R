test_that("a manifest round-trips and detects a changed or missing input", {
  a <- tempfile(); writeLines("alpha", a)
  b <- tempfile(); writeLines("beta",  b)
  js <- tempfile(fileext = ".json")

  m <- hf_input_manifest(reference = a, fdr = b, out = js)
  expect_equal(m$role, c("reference", "fdr"))
  expect_true(all(m$exists))
  expect_true(file.exists(js))

  expect_silent(suppressMessages(hf_verify_manifest(js)))

  # Same path, different content.
  writeLines("alpha changed", a)
  v <- suppressWarnings(suppressMessages(hf_verify_manifest(js)))
  expect_equal(v$status[v$role == "reference"], "changed")
  expect_equal(v$status[v$role == "fdr"], "ok")
  expect_error(suppressMessages(hf_verify_manifest(js, strict = TRUE)),
               "no longer match")

  unlink(b)
  v2 <- suppressWarnings(suppressMessages(hf_verify_manifest(js)))
  expect_equal(v2$status[v2$role == "fdr"], "missing")
})

test_that("inputs must be named, and NULLs are dropped", {
  f <- tempfile(); writeLines("x", f)
  expect_error(hf_input_manifest(f), "must be named")
  expect_error(hf_input_manifest(), "no inputs")
  expect_equal(nrow(hf_input_manifest(reference = f, lakes = NULL)), 1L)
})

test_that("a missing input is recorded rather than fatal", {
  m <- suppressWarnings(hf_input_manifest(reference = "/no/such/file.gpkg"))
  expect_false(m$exists)
  expect_true(is.na(m$md5))
})
