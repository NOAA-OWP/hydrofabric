# Hermetic tests for the pure masked multi-source D8 divide partition
# (engine = "pure" in simple_split_divides_hfrunners). No real FDR/FAC needed:
# .label_parent_d8 works on a plain matrix and .labels_to_pieces on a synthetic
# label raster, so these are fast and deterministic.

test_that(".label_parent_d8 tiles a 1-D channel (first outlet downstream wins)", {
  # 5 cells flowing South (Esri code 4); outlets at row2 (label 1), row4 (label 2).
  fdr <- matrix(c(4, 4, 4, 4, 4), nrow = 5, ncol = 1)
  lab <- hydrofabric:::.label_parent_d8(fdr, matrix(c(2, 1, 4, 1), ncol = 2, byrow = TRUE),
    c(1L, 2L), terminal_label = 3L)
  # row1->outlet1, row2 outlet1, row3->outlet2, row4 outlet2, row5->exit->terminal
  expect_equal(as.vector(lab), c(1, 1, 2, 2, 3))
})

test_that(".label_parent_d8 handles nested outlets on one channel", {
  fdr <- matrix(rep(4, 6), nrow = 6, ncol = 1)             # all flow South
  lab <- hydrofabric:::.label_parent_d8(fdr, matrix(c(2, 1, 4, 1), ncol = 2, byrow = TRUE),
    c(1L, 2L), terminal_label = 3L)
  expect_equal(as.vector(lab), c(1, 1, 2, 2, 3, 3))
})

test_that(".label_parent_d8 routes a 2-D grid to its single outlet", {
  # cols 1-2 flow East (1); col 3 flows South (4); bottom of col 3 exits.
  fdr <- matrix(c(1, 1, 4,
    1, 1, 4,
    1, 1, 4), nrow = 3, byrow = TRUE)
  lab <- hydrofabric:::.label_parent_d8(fdr, matrix(c(2, 3), ncol = 2), 1L, terminal_label = 2L)
  expect_equal(lab[1, ], c(1, 1, 1))   # rows 1-2 drain through the (2,3) outlet
  expect_equal(lab[2, ], c(1, 1, 1))
  expect_equal(lab[3, ], c(2, 2, 2))   # row 3 drains past it -> terminal
})

test_that(".label_parent_d8 sends flats and cycles to terminal (no hang)", {
  # row2 has no valid flow code (0); a 2-cell mutual cycle elsewhere is terminal too.
  fdr <- matrix(c(4, 0, 4), nrow = 3, ncol = 1)            # row2 = flat
  lab <- hydrofabric:::.label_parent_d8(fdr, matrix(c(3, 1), ncol = 2), 1L, terminal_label = 9L)
  # row3 = outlet(1); row2 flat -> exit -> terminal(9); row1 -> row2(flat) -> terminal(9)
  expect_equal(as.vector(lab), c(9, 9, 1))

  cyc <- matrix(c(4, 64), nrow = 2, ncol = 1)              # row1->row2, row2->row1
  labc <- hydrofabric:::.label_parent_d8(cyc, matrix(integer(0), ncol = 2), integer(0), terminal_label = 7L)
  expect_equal(as.vector(labc), c(7, 7))
})

test_that(".label_parent_d8 keeps every non-NA cell labeled and NA cells at 0", {
  fdr <- matrix(c(4, NA, 4, 4), nrow = 4, ncol = 1)
  lab <- hydrofabric:::.label_parent_d8(fdr, matrix(c(4, 1), ncol = 2), 1L, terminal_label = 5L)
  expect_true(is.na(fdr[2, 1]) && lab[2, 1] == 0L)         # NA cell stays 0
  expect_true(all(lab[c(1, 3, 4), 1] != 0L))               # non-NA cells labeled
})

test_that(".labels_to_pieces yields disjoint single polygons tiling the parent", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  # 10x10 grid: left 5 cols = label 1, right 5 cols = label 2 (terminal).
  m <- matrix(c(rep(1L, 5), rep(2L, 5)), nrow = 10, ncol = 10, byrow = TRUE)
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 100,
    ymin = 0, ymax = 100, crs = "EPSG:5070")
  terra::values(r) <- as.vector(t(m))
  parent <- sf::st_sf(
    divide_id = "P", flowpath_id = "P",
    geometry  = sf::st_sfc(sf::st_polygon(list(rbind(
      c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)))), crs = 5070))

  pieces <- hydrofabric:::.labels_to_pieces(r, parent, label_to_id = c("1" = "A"),
    terminal_label = 2L, terminal_id = "B")

  expect_s3_class(pieces, "sf")
  expect_setequal(as.character(pieces$divide_id), c("A", "B"))
  # single POLYGON per piece
  expect_true(all(as.character(sf::st_geometry_type(pieces)) == "POLYGON"))
  # divide_id unique (flowpaths_no_shared_divide)
  expect_equal(anyDuplicated(pieces$divide_id), 0L)
  # disjoint (divides_non_overlapping)
  ov <- sf::st_overlaps(pieces, sparse = FALSE)
  expect_false(any(ov[upper.tri(ov)]))
  # exact tiling: union area == parent area, each half == 5000
  areas <- as.numeric(sf::st_area(pieces))
  expect_equal(sum(areas), 10000, tolerance = 1e-6)
  expect_true(all(abs(areas - 5000) < 1e-6))
})
