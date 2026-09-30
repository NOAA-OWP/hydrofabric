# The splitter states two properties in prose and relies on both. Neither was
# ever checked, and one of them was false.
#
#   1. members are ordered "upstream-most by point order"
#   2. the masked multi-source D8 partition is "non-overlap ... by construction"
#
# Both are decidable without rasters or geometry, so they are asserted here
# rather than asserted in a comment.

# ---- 1. topological order --------------------------------------------------

test_that("members are ordered most-upstream first, whatever their ids sort like", {
  # a -> b -> c leaving the group. Ids chosen so the lexicographic sort ("10"
  # before "2") disagrees with the chain, which is exactly the case that made
  # the previous id sort silently wrong.
  ids   <- c("2", "10", "7")
  dsids <- c("10", "7", "outside")     # 2 -> 10 -> 7 -> (sink)
  ord <- hydrofabric:::.topo_order_members(ids, dsids)
  expect_identical(ids[ord], c("2", "10", "7"))

  # same chain, rows presented in a different order: the answer must not change
  ord2 <- hydrofabric:::.topo_order_members(c("7", "2", "10"),
                                            c("outside", "10", "7"))
  expect_identical(c("7", "2", "10")[ord2], c("2", "10", "7"))
})

test_that("the order is a permutation, so no member is dropped or duplicated", {
  ids   <- c("5", "3", "40", "1")
  dsids <- c("3", "40", "1", "sink")   # 5 -> 3 -> 40 -> 1
  ord <- hydrofabric:::.topo_order_members(ids, dsids)
  expect_setequal(ord, seq_along(ids))
  expect_identical(ids[ord], c("5", "3", "40", "1"))
})

test_that("a chain that is not simple returns NULL rather than a false order", {
  # fork: two members share a successor
  expect_null(hydrofabric:::.topo_order_members(c("1", "2", "3"),
                                                c("3", "3", "out")))
  # cycle
  expect_null(hydrofabric:::.topo_order_members(c("1", "2"), c("2", "1")))
  # break: two disconnected pairs, so no single head
  expect_null(hydrofabric:::.topo_order_members(c("1", "2", "3", "4"),
                                                c("2", "out", "4", "out")))
  # duplicate ids
  expect_null(hydrofabric:::.topo_order_members(c("1", "1"), c("2", "2")))
})

test_that("a single member and an empty group are handled", {
  expect_identical(hydrofabric:::.topo_order_members("1", "out"), 1L)
  expect_identical(hydrofabric:::.topo_order_members(character(0), character(0)),
                   integer(0))
})

# ---- 2. the partition really partitions ------------------------------------

# A 5x5 grid draining east (FDR 1) into a column that drains south (FDR 4), so
# every cell has a path to the bottom-right, which then routes off the grid.
#
# The exit cell keeps a direction rather than NA: a cell with no flow direction
# is not a routable target, so an outlet seeded there collects nothing and its
# would-be drainage falls to the terminal label instead. `.snap_outlet_cell`
# only ever returns cells with valid FDR, so this mirrors production.
.fdr_grid <- function() {
  m <- matrix(1L, nrow = 5, ncol = 5)   # 1 = east
  m[, 5] <- 4L                          # last column drains south, off the grid
  m
}

test_that("every routable cell is claimed by exactly one outlet", {
  m <- .fdr_grid()
  outlets <- rbind(c(3L, 5L), c(5L, 5L))     # mid-column and the exit
  lab <- hydrofabric:::.label_parent_d8(m, outlets, c(1L, 2L), terminal_label = 3L)
  expect_false(is.null(lab))
  expect_identical(dim(lab), dim(m))

  # "non-overlap by construction" means each cell carries exactly one label,
  # a matrix cell cannot hold two, so the property to check is that no routable
  # cell is left unlabelled and every label is one of the ones we supplied.
  routable <- !is.na(m)
  expect_true(all(lab[routable] != 0L))
  expect_true(all(lab[routable] %in% c(1L, 2L, 3L)))
})

test_that("cells above a sub-outlet go to it, not to the one below", {
  m <- .fdr_grid()
  outlets <- rbind(c(3L, 5L), c(5L, 5L))
  lab <- hydrofabric:::.label_parent_d8(m, outlets, c(1L, 2L), terminal_label = 3L)
  # row 1 drains east then south through (3,5), so it belongs to outlet 1
  expect_identical(unique(lab[1, 1:4]), 1L)
  # row 4 is below (3,5) and reaches the grid exit at (5,5): outlet 2
  expect_identical(unique(lab[4, 1:4]), 2L)
})

test_that("the labelling is independent of the order outlets are supplied in", {
  m <- .fdr_grid()
  a <- hydrofabric:::.label_parent_d8(m, rbind(c(3L, 5L), c(5L, 5L)),
                                      c(1L, 2L), terminal_label = 3L)
  b <- hydrofabric:::.label_parent_d8(m, rbind(c(5L, 5L), c(3L, 5L)),
                                      c(2L, 1L), terminal_label = 3L)
  expect_identical(a, b)
})
