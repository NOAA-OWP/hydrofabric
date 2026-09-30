# .clear_disjoint_fp_divides must not take a divide the gate would accept.
#
# The pass judged the aggregated flowpath line: zero overlap with its own
# polygon meant misassigned. But build_nextgen deliberately does not clip that
# line to its divide, because an aggregated mainstem legitimately runs outside
# the catchment where its member reaches do. The property that must hold is
# divide_has_member_flowline, a MEMBER flowline within 30 m. On VPU 09 all
# five divides this pass cleared satisfied that and were absorbed anyway.

.cdj_sq <- function(x0, y0, w, h) sf::st_polygon(list(rbind(
  c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + h), c(x0, y0 + h), c(x0, y0))))
.cdj_ln <- function(...) sf::st_linestring(rbind(...))

# A divide, a flowpath whose aggregated line runs entirely outside it, and a
# member flowline that sits inside it.
.cdj_case <- function(member_inside = TRUE) {
  div <- sf::st_sf(
    divide_id = "cat-1", flowpath_id = "fp-1",
    geometry = sf::st_sfc(.cdj_sq(0, 0, 1000, 1000), crs = 5070))
  fp <- sf::st_sf(
    flowpath_id = "fp-1", divide_id = "cat-1", has_divide = TRUE,
    geometry = sf::st_sfc(.cdj_ln(c(2000, 500), c(3000, 500)), crs = 5070))
  fl <- sf::st_sf(
    flowpath_id = "fp-1",
    geometry = sf::st_sfc(
      if (member_inside) .cdj_ln(c(100, 500), c(900, 500))     # inside the divide
      else               .cdj_ln(c(5000, 500), c(5900, 500)),  # far from it
      crs = 5070))
  list(fp = fp, div = div, fl = fl)
}

test_that("a divide whose member flowline is inside it is not cleared", {
  x <- .cdj_case(member_inside = TRUE)
  r <- hydrofabric:::.clear_disjoint_fp_divides(x$fp, x$div, flowlines = x$fl)
  expect_equal(nrow(r$report), 0L)
  expect_equal(as.character(r$fp$divide_id[1]), "cat-1")
  expect_true(r$fp$has_divide[1])
  expect_equal(nrow(r$div), 1L)          # the divide row survives
})

test_that("a divide with no member flowline anywhere near it is still cleared", {
  # The guard must not disable the repair: a genuine misassignment still goes.
  x <- .cdj_case(member_inside = FALSE)
  r <- hydrofabric:::.clear_disjoint_fp_divides(x$fp, x$div, flowlines = x$fl)
  expect_equal(nrow(r$report), 1L)
  expect_true(is.na(r$fp$divide_id[1]))
  expect_false(r$fp$has_divide[1])
})

test_that("behaviour is unchanged when no flowlines are supplied", {
  # Older callers, and any gpkg without the layer, must behave exactly as before.
  x <- .cdj_case(member_inside = TRUE)
  r <- hydrofabric:::.clear_disjoint_fp_divides(x$fp, x$div)
  expect_equal(nrow(r$report), 1L)
  expect_true(is.na(r$fp$divide_id[1]))
})

test_that("a member flowline just inside the 30 m tolerance protects the divide", {
  # The gate accepts within 30 m, so this pass must too, or the two disagree at
  # the boundary and divides are absorbed that the gate would pass.
  x <- .cdj_case()
  x$fl <- sf::st_sf(flowpath_id = "fp-1",
    geometry = sf::st_sfc(.cdj_ln(c(100, 1020), c(900, 1020)), crs = 5070))  # 20 m above
  r <- hydrofabric:::.clear_disjoint_fp_divides(x$fp, x$div, flowlines = x$fl)
  expect_equal(nrow(r$report), 0L)
})

test_that("a member flowline beyond the tolerance does not protect it", {
  x <- .cdj_case()
  x$fl <- sf::st_sf(flowpath_id = "fp-1",
    geometry = sf::st_sfc(.cdj_ln(c(100, 1100), c(900, 1100)), crs = 5070))  # 100 m above
  r <- hydrofabric:::.clear_disjoint_fp_divides(x$fp, x$div, flowlines = x$fl)
  expect_equal(nrow(r$report), 1L)
})
