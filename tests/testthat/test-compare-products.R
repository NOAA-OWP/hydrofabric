.mk_prod <- function(fp = c("fp-1", "fp-2", "fp-3"),
                     div_fp = c("fp-1", "fp-2", "fp-3"),
                     areas = c(10, 20, 30)) {
  f <- tempfile(fileext = ".gpkg")
  pts <- function(ids) sf::st_as_sf(
    data.frame(flowpath_id = ids, x = seq_along(ids), y = seq_along(ids)),
    coords = c("x", "y"), crs = 5070)
  sf::st_write(pts(fp), f, "flowpaths", quiet = TRUE)
  d <- pts(div_fp); d$areasqkm <- areas
  sf::st_write(d, f, "divides", quiet = TRUE)
  f
}

test_that("identical products report no regression", {
  a <- .mk_prod()
  r <- suppressMessages(hf_compare_products(a, a))
  expect_true(r$ok)
  expect_equal(nrow(r$lost_divide), 0L)
})

test_that("a flowpath losing its divide is named, and strict raises", {
  a <- .mk_prod()
  b <- .mk_prod(div_fp = c("fp-1", "fp-3"), areas = c(10, 30))
  r <- suppressWarnings(suppressMessages(hf_compare_products(a, b)))
  expect_false(r$ok)
  expect_equal(r$lost_divide$flowpath_id, "fp-2")
  expect_equal(r$counts$delta[r$counts$layer == "divides"], -1)
  expect_error(suppressMessages(hf_compare_products(a, b, strict = TRUE)),
               "lost a divide")
})

test_that("a gained divide is reported but is not a regression", {
  a <- .mk_prod(div_fp = c("fp-1", "fp-2"), areas = c(10, 20))
  b <- .mk_prod()
  r <- suppressMessages(hf_compare_products(a, b))
  expect_true(r$ok)
  expect_equal(r$gained_divide$flowpath_id, "fp-3")
})

test_that("an area change beyond tolerance is reported, within it is not", {
  a <- .mk_prod()
  b <- .mk_prod(areas = c(10, 26, 30))          # fp-2 +30%
  r <- suppressMessages(hf_compare_products(a, b, area_tol = 0.01))
  expect_true(r$ok)                              # nothing lost
  expect_equal(r$area_changed$flowpath_id, "fp-2")

  q <- suppressMessages(hf_compare_products(a, b, area_tol = 0.5))
  expect_equal(nrow(q$area_changed), 0L)
})
