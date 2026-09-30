# .allocate_parent_parts: give a split parent's displaced land to the child
# reach nearest each part.
#
# The D8 partition can hand a split reach's fragment a sub-polygon that sits
# nowhere near it. What it replaced made two spatially arbitrary choices: the
# whole parent to the longest flowline, or the same leftover polygon to every
# displaced fragment. Fixtures are at km scale in EPSG:5070 so distances read
# as metres.

.app_sq <- function(x0, y0, w, h) sf::st_polygon(list(rbind(
  c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + h), c(x0, y0 + h), c(x0, y0))))

.app_line <- function(x0, y0, x1, y1)
  sf::st_sfc(sf::st_linestring(rbind(c(x0, y0), c(x1, y1))), crs = 5070)

.app_area <- function(g) sum(as.numeric(
  sf::st_area(sf::st_sfc(g, crs = 5070))))

test_that("two swapped partitions are exchanged", {
  # Parent split into a west and an east half. The partition handed the west
  # polygon to the east fragment and vice versa, so neither meets its own line.
  west <- .app_sq(0, 0, 1000, 1000)
  east <- .app_sq(1000, 0, 1000, 1000)
  lines <- list(.app_line(100, 500, 900, 500),     # fragment 1 lives west
                .app_line(1100, 500, 1900, 500))   # fragment 2 lives east
  geoms <- list(east, west)                        # swapped
  out <- hydrofabric:::.allocate_parent_parts(
    geoms, lines, disjoint = c(TRUE, TRUE), crs = sf::st_crs(5070))

  expect_false(is.null(out))
  expect_equal(.app_area(out), .app_area(geoms), tolerance = 1e-9)
  # each fragment now holds the half its own line runs through
  expect_true(sf::st_intersects(sf::st_sfc(out[[1]], crs = 5070),
                                lines[[1]], sparse = FALSE)[1, 1])
  expect_true(sf::st_intersects(sf::st_sfc(out[[2]], crs = 5070),
                                lines[[2]], sparse = FALSE)[1, 1])
})

test_that("a correctly placed partition is not disturbed", {
  # Fragment 1 is valid; only fragment 2 is displaced. The old code gave every
  # displaced row the parent-minus-valid remainder regardless of where it was.
  west <- .app_sq(0, 0, 1000, 1000)
  east <- .app_sq(1000, 0, 1000, 1000)
  lines <- list(.app_line(100, 500, 900, 500), .app_line(1100, 500, 1900, 500))
  out <- hydrofabric:::.allocate_parent_parts(
    list(west, west), lines, disjoint = c(FALSE, TRUE), crs = sf::st_crs(5070))
  expect_identical(out[[1]], west)   # valid row untouched
})

test_that("the whole parent is never handed to one fragment", {
  # The branch this replaces gave the union to whichever fragment had the
  # longest flowline, leaving the others displaced. Fragment 2's line is far
  # longer, but it must not end up owning the west half.
  west <- .app_sq(0, 0, 1000, 1000)
  east <- .app_sq(1000, 0, 1000, 1000)
  lines <- list(.app_line(400, 500, 600, 500),      # short
                .app_line(1010, 500, 1990, 500))    # long
  out <- hydrofabric:::.allocate_parent_parts(
    list(east, west), lines, disjoint = c(TRUE, TRUE), crs = sf::st_crs(5070))
  a1 <- .app_area(out[[1]]); a2 <- .app_area(out[[2]])
  expect_equal(a1, 1e6, tolerance = 1e-6)
  expect_equal(a2, 1e6, tolerance = 1e-6)
})

test_that("a multi-part displaced polygon is split across fragments", {
  # One fragment holds both halves as a MULTIPOLYGON. Parts are the unit of
  # reassignment, so they must separate.
  both <- sf::st_multipolygon(list(.app_sq(0, 0, 1000, 1000),
                                   .app_sq(1000, 0, 1000, 1000)))
  lines <- list(.app_line(100, 500, 900, 500), .app_line(1100, 500, 1900, 500))
  out <- hydrofabric:::.allocate_parent_parts(
    list(both, sf::st_polygon()), lines, disjoint = c(TRUE, TRUE),
    crs = sf::st_crs(5070))
  expect_equal(.app_area(out[[1]]), 1e6, tolerance = 1e-6)
  expect_equal(.app_area(out[[2]]), 1e6, tolerance = 1e-6)
})

test_that("area is conserved and no fragment claims the same ground", {
  parts <- list(.app_sq(0, 0, 1000, 1000), .app_sq(1000, 0, 1000, 1000),
                .app_sq(2000, 0, 1000, 1000))
  geoms <- list(parts[[3]], parts[[1]], parts[[2]])
  lines <- list(.app_line(100, 500, 900, 500),
                .app_line(1100, 500, 1900, 500),
                .app_line(2100, 500, 2900, 500))
  out <- hydrofabric:::.allocate_parent_parts(
    geoms, lines, disjoint = c(TRUE, TRUE, TRUE), crs = sf::st_crs(5070))
  expect_equal(.app_area(out), 3e6, tolerance = 1e-6)
  sfc <- sf::st_sfc(out, crs = 5070)
  ov  <- sf::st_intersection(sfc, sfc)
  expect_equal(sum(as.numeric(sf::st_area(ov))), 3e6, tolerance = 1e-6)
})

test_that("nothing is done when no partition is displaced", {
  west <- .app_sq(0, 0, 1000, 1000)
  expect_null(hydrofabric:::.allocate_parent_parts(
    list(west, west), list(.app_line(1, 1, 2, 2), .app_line(1, 1, 2, 2)),
    disjoint = c(FALSE, FALSE), crs = sf::st_crs(5070)))
})

test_that("a group with no locatable fragment is declined", {
  west <- .app_sq(0, 0, 1000, 1000)
  expect_null(hydrofabric:::.allocate_parent_parts(
    list(west, west), list(NULL, NULL), disjoint = c(TRUE, TRUE),
    crs = sf::st_crs(5070)))
})

test_that("an unlocatable fragment neither gives up land nor receives any", {
  west <- .app_sq(0, 0, 1000, 1000)
  east <- .app_sq(1000, 0, 1000, 1000)
  out <- hydrofabric:::.allocate_parent_parts(
    list(east, west), list(.app_line(100, 500, 900, 500), NULL),
    disjoint = c(TRUE, FALSE), crs = sf::st_crs(5070))
  expect_identical(out[[2]], west)                       # untouched
  expect_equal(.app_area(out), 2e6, tolerance = 1e-6)    # nothing lost
})

test_that("a single-fragment group is a no-op", {
  expect_null(hydrofabric:::.allocate_parent_parts(
    list(.app_sq(0, 0, 10, 10)), list(.app_line(1, 1, 2, 2)),
    disjoint = TRUE, crs = sf::st_crs(5070)))
})

test_that("a fragment that gave up land keeps a share of it", {
  # The starvation guard. Both parts lie at fragment 1's end, so nearest-wins
  # alone would give fragment 2 nothing and ship it with no divide. On a real AK
  # build that stripped 95 flowpaths of their catchment. Fragment 2 must get its
  # own nearest part back, and area must still be conserved.
  p1 <- .app_sq(0, 0, 500, 500)
  p2 <- .app_sq(500, 0, 500, 500)
  lines <- list(.app_line(100, 250, 400, 250),      # nearest to p1
                .app_line(700, 900, 900, 900))      # north of p2, far from p1
  out <- hydrofabric:::.allocate_parent_parts(
    list(sf::st_multipolygon(list(p1, p2)), sf::st_polygon()),
    lines, disjoint = c(TRUE, TRUE), crs = sf::st_crs(5070))

  expect_false(sf::st_is_empty(sf::st_sfc(out[[2]], crs = 5070)))
  expect_equal(.app_area(out), 5e5, tolerance = 1e-6)     # nothing created or lost
  expect_gt(.app_area(out[[1]]), 0)
})

test_that("a donor is never emptied to feed a starved fragment", {
  # Only fragments holding more than one part may donate, so the guard cannot
  # simply move the shortage around.
  p1 <- .app_sq(0, 0, 500, 500)
  lines <- list(.app_line(100, 250, 400, 250), .app_line(9000, 9000, 9400, 9000))
  out <- hydrofabric:::.allocate_parent_parts(
    list(p1, sf::st_polygon()), lines, disjoint = c(TRUE, TRUE),
    crs = sf::st_crs(5070))
  # one part, two fragments: the nearer one keeps it rather than both losing out
  expect_equal(.app_area(out), 2.5e5, tolerance = 1e-6)
  expect_false(sf::st_is_empty(sf::st_sfc(out[[1]], crs = 5070)))
})

test_that("a fragment whose parent all drains elsewhere is left with no land", {
  # Documented consequence of the rule, not an accident. Both parts of the
  # parent lie at fragment 1's end, so fragment 2: whose own partition was
  # displaced: receives nothing and ends up with an empty polygon rather than
  # the wrong catchment it held before. Area still lands entirely on fragment 1.
  #
  # The alternative would be to force one part onto every fragment, which puts
  # land under a reach it does not drain. If a build ever shows this producing
  # divide-less flowpaths at scale, that is the trade-off to revisit.
  # Only reachable when the fragment never held land to begin with: it is not
  # marked disjoint, so the guard has no claim to honour on its behalf.
  p1 <- .app_sq(0, 0, 500, 500)
  lines <- list(.app_line(100, 250, 400, 250), .app_line(5000, 250, 5800, 250))
  out <- hydrofabric:::.allocate_parent_parts(
    list(p1, sf::st_polygon()), lines, disjoint = c(TRUE, FALSE),
    crs = sf::st_crs(5070))

  expect_equal(.app_area(out[[1]]), 2.5e5, tolerance = 1e-6)
  expect_true(sf::st_is_empty(sf::st_sfc(out[[2]], crs = 5070)))
  expect_equal(.app_area(out), 2.5e5, tolerance = 1e-6)   # nothing lost
})

# .cut_parent_by_nearest_line: divide a split parent's catchment among the
# fragments of the reach, by which fragment's line each piece of land is
# nearest. Replaces repairing the D8 trace's output with recomputing it.

test_that("the parent is cut so each fragment gets the land nearest its line", {
  parent <- .app_sq(0, 0, 2000, 1000)
  lines  <- list(.app_line(100, 500, 900, 500),      # west half
                 .app_line(1100, 500, 1900, 500))    # east half
  out <- hydrofabric:::.cut_parent_by_nearest_line(parent, lines, sf::st_crs(5070))
  expect_length(out, 2)
  # area conserved exactly: the pieces tile the parent
  expect_equal(.app_area(out), 2e6, tolerance = 1e-6)
  # and split about evenly, since the lines are symmetric
  expect_equal(.app_area(out[[1]]), 1e6, tolerance = 1e4)
  expect_equal(.app_area(out[[2]]), 1e6, tolerance = 1e4)
  # each fragment's own line lies in its own piece
  for (i in 1:2)
    expect_true(sf::st_intersects(sf::st_sfc(out[[i]], crs = 5070),
                                  lines[[i]], sparse = FALSE)[1, 1])
})

test_that("no fragment is starved, even one the trace gave nothing", {
  # The failure the cut exists to remove: on AK this left 191 fragments with no
  # catchment at all, because whole-part reassignment cannot manufacture a piece
  # for a fragment the partition skipped.
  parent <- .app_sq(0, 0, 2000, 1000)
  lines  <- list(.app_line(100, 500, 1800, 500),   # long, near everything
                 .app_line(1900, 100, 1950, 150))  # short, tucked in a corner
  out <- hydrofabric:::.cut_parent_by_nearest_line(parent, lines, sf::st_crs(5070))
  expect_false(sf::st_is_empty(sf::st_sfc(out[[2]], crs = 5070)))
  expect_gt(.app_area(out[[2]]), 0)
  expect_equal(.app_area(out), 2e6, tolerance = 1e-6)
})

test_that("the pieces do not overlap", {
  parent <- .app_sq(0, 0, 3000, 1000)
  lines  <- list(.app_line(100, 500, 900, 500),
                 .app_line(1100, 500, 1900, 500),
                 .app_line(2100, 500, 2900, 500))
  out <- hydrofabric:::.cut_parent_by_nearest_line(parent, lines, sf::st_crs(5070))
  sfc <- sf::st_sfc(out, crs = 5070)
  # summed area equals the footprint: no ground is claimed twice
  expect_equal(sum(as.numeric(sf::st_area(sfc))),
               as.numeric(sf::st_area(sf::st_union(sfc))), tolerance = 1e-6)
  expect_equal(.app_area(out), 3e6, tolerance = 1e-6)
})

test_that("a single locatable fragment takes the whole parent", {
  parent <- .app_sq(0, 0, 1000, 1000)
  out <- hydrofabric:::.cut_parent_by_nearest_line(
    parent, list(.app_line(100, 500, 900, 500), NULL), sf::st_crs(5070))
  expect_equal(.app_area(out[[1]]), 1e6, tolerance = 1e-6)
  expect_true(sf::st_is_empty(sf::st_sfc(out[[2]], crs = 5070)))
})

test_that("a group with no locatable fragment is declined", {
  expect_null(hydrofabric:::.cut_parent_by_nearest_line(
    .app_sq(0, 0, 1000, 1000), list(NULL, NULL), sf::st_crs(5070)))
})
