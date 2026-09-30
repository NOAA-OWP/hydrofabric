test_that("the Group-B swap never takes a divide another flowpath already owns", {
  skip_if_not(requireNamespace("sf", quietly = TRUE))
  skip_if_not(requireNamespace("geos", quietly = TRUE))

  # Two catchments side by side. fp-1 is correctly placed inside cat-1. fp-2 is a
  # misplaced (disjoint) flowpath: it is assigned cat-2, which sits far away, but
  # its line actually runs through cat-1, the classic FDR headwater
  # misattribution the swap exists to fix.
  #
  # The swap must NOT hand cat-1 to fp-2, because fp-1 already owns it. Doing so
  # leaves two flowpaths on one divide, which is the state Step 5e was written to
  # undo. With no unowned overlapping divide available, the correct outcome is to
  # clear fp-2's assignment and let it become a connector.
  sq <- function(x0, w = 1000) sf::st_polygon(list(cbind(
    c(x0, x0 + w, x0 + w, x0, x0), c(0, 0, w, w, 0))))

  div <- sf::st_sf(
    divide_id = c("cat-1", "cat-2"),
    geometry  = sf::st_sfc(sq(0), sq(5000), crs = 5070))

  fp <- sf::st_sf(
    flowpath_id = c("fp-1", "fp-2"),
    divide_id   = c("cat-1", "cat-2"),
    has_divide  = c(TRUE, TRUE),
    geometry    = sf::st_sfc(
      sf::st_linestring(cbind(c(100, 900), c(500, 500))),   # inside cat-1
      sf::st_linestring(cbind(c(200, 800), c(300, 300))),   # also inside cat-1
      crs = 5070))

  res <- hydrofabric:::.clear_disjoint_fp_divides(fp, div, drop_orphan_cats = FALSE)

  # fp-2 was detected as disjoint from cat-2.
  expect_equal(nrow(res$report), 1L)
  expect_identical(res$report$flowpath_id, "fp-2")

  out <- res$fp
  did <- stats::setNames(as.character(out$divide_id), as.character(out$flowpath_id))

  # The owner keeps its divide.
  expect_identical(unname(did["fp-1"]), "cat-1")
  # The misplaced flowpath must not have been given the owned divide.
  expect_true(is.na(did[["fp-2"]]))
  if ("has_divide" %in% names(out))
    expect_false(out$has_divide[out$flowpath_id == "fp-2"])

  # The invariant this protects: no divide referenced by more than one flowpath.
  claimed <- did[!is.na(did)]
  expect_false(any(duplicated(claimed)))
})

test_that("the swap still reassigns when the overlapping divide is genuinely unowned", {
  skip_if_not(requireNamespace("sf", quietly = TRUE))
  skip_if_not(requireNamespace("geos", quietly = TRUE))

  # Same shape, except cat-1 has no owner. The swap is then the right move and
  # must still happen, the ownership guard must not disable the repair itself.
  sq <- function(x0, w = 1000) sf::st_polygon(list(cbind(
    c(x0, x0 + w, x0 + w, x0, x0), c(0, 0, w, w, 0))))

  div <- sf::st_sf(
    divide_id = c("cat-1", "cat-2"),
    geometry  = sf::st_sfc(sq(0), sq(5000), crs = 5070))

  fp <- sf::st_sf(
    flowpath_id = "fp-2",
    divide_id   = "cat-2",
    has_divide  = TRUE,
    geometry    = sf::st_sfc(
      sf::st_linestring(cbind(c(200, 800), c(300, 300))), crs = 5070))

  res <- hydrofabric:::.clear_disjoint_fp_divides(fp, div, drop_orphan_cats = FALSE)
  expect_identical(as.character(res$fp$divide_id), "cat-1")
  expect_true(res$fp$has_divide)
})
