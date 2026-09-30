test_that("negative FAC readings become NA", {
  # -2147483648 is the INT32 sentinel observed in a source tile that declares
  # NoData as nan, so GDAL hands it back as an ordinary value.
  x <- hydrofabric:::.fac_drop_sentinels(c(10, -2147483648, 0, 3.5))
  expect_equal(x, c(10, NA, 0, 3.5))
})

test_that("zero is kept and existing NAs are preserved", {
  x <- hydrofabric:::.fac_drop_sentinels(c(0, NA_real_, 7))
  expect_equal(x, c(0, NA, 7))
})

test_that("a clean vector is returned unchanged and silently", {
  x <- c(1, 2, 3)
  expect_silent(out <- hydrofabric:::.fac_drop_sentinels(x))
  expect_equal(out, x)
})

test_that("an all-negative vector collapses to all NA", {
  expect_true(all(is.na(hydrofabric:::.fac_drop_sentinels(c(-1, -2147483648)))))
})

test_that("an empty vector is handled", {
  expect_length(hydrofabric:::.fac_drop_sentinels(numeric(0)), 0)
})
