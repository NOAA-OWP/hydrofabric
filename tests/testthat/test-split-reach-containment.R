# A split piece is a divide for a child only if it carries that child's reach.
# The D8 delineation can land the piece somewhere else entirely: measured on
# Alaska, 260 children sat a median 310 m from their own reach, and the nextgen phase
# then clears the pairing, so the child ships with no divide at all. These
# helpers decide that at split time instead.
suppressMessages({ library(sf) })

.sq <- function(x0, y0, w, h)
  st_polygon(list(cbind(c(x0, x0 + w, x0 + w, x0, x0),
                        c(y0, y0, y0 + h, y0 + h, y0))))
.ln <- function(x0, y0, x1, y1) st_linestring(cbind(c(x0, x1), c(y0, y1)))

test_that("a piece holding its reach passes, one that does not fails", {
  piece <- st_sfc(.sq(0, 0, 100, 100), crs = 5070)
  inside  <- st_sfc(.ln(10, 10, 90, 90),     crs = 5070)
  far     <- st_sfc(.ln(500, 500, 600, 600), crs = 5070)
  expect_true(hydrofabric:::.piece_holds_reach(piece, inside))
  expect_false(hydrofabric:::.piece_holds_reach(piece, far))
})

test_that("the 30 m tolerance is the boundary, matching the rest of the fabric", {
  piece <- st_sfc(.sq(0, 0, 100, 100), crs = 5070)
  near <- st_sfc(.ln(120, 10, 120, 90), crs = 5070)   # 20 m clear of the edge
  off  <- st_sfc(.ln(140, 10, 140, 90), crs = 5070)   # 40 m clear of the edge
  expect_true(hydrofabric:::.piece_holds_reach(piece, near))
  expect_false(hydrofabric:::.piece_holds_reach(piece, off))
})

test_that("absent reach geometry keeps the previous behaviour", {
  piece <- st_sfc(.sq(0, 0, 100, 100), crs = 5070)
  expect_true(hydrofabric:::.piece_holds_reach(piece, NULL))
  expect_true(hydrofabric:::.piece_holds_reach(NULL, NULL))
})

test_that("parts are ordered by the reach each one carries, not by row order", {
  # Two side-by-side halves. The upstream reach lies in the RIGHT half, so the
  # correct answer inverts the row order the caller would otherwise assume.
  parts <- st_sf(geometry = st_sfc(.sq(0, 0, 100, 100), .sq(100, 0, 100, 100),
                                   crs = 5070))
  up <- st_sfc(.ln(150, 10, 150, 90), crs = 5070)   # right half
  dn <- st_sfc(.ln(50, 10, 50, 90),   crs = 5070)   # left half
  expect_identical(hydrofabric:::.order_parts_by_reach(parts, up, dn), c(2L, 1L))
  # and the un-inverted case still resolves to row order
  expect_identical(hydrofabric:::.order_parts_by_reach(parts, dn, up), c(1L, 2L))
})

test_that("an ambiguous pairing still returns an order rather than nothing", {
  # Returning NULL here USED to mean "refuse the split". It no longer does:
  # the caller applies the order only when one is returned, so NULL now leaves
  # `lwgeom::st_split`'s part order in place, and that order is undefined.
  # Scoring by reach coverage and taking the best pairing is therefore strictly
  # better than declining, even when neither part is a clean match.
  parts <- st_sf(geometry = st_sfc(.sq(0, 0, 100, 100), .sq(100, 0, 100, 100),
                                   crs = 5070))
  one <- st_sfc(.ln(50, 10, 50, 90), crs = 5070)
  far <- st_sfc(.ln(9000, 9000, 9100, 9100), crs = 5070)
  expect_length(hydrofabric:::.order_parts_by_reach(parts, one, one), 2L)
  expect_length(hydrofabric:::.order_parts_by_reach(parts, one, far), 2L)
})

test_that("no order is possible without two parts or without both reaches", {
  parts <- st_sf(geometry = st_sfc(.sq(0, 0, 100, 100), .sq(100, 0, 100, 100),
                                   crs = 5070))
  one <- st_sfc(.ln(50, 10, 50, 90), crs = 5070)
  far <- st_sfc(.ln(9000, 9000, 9100, 9100), crs = 5070)
  expect_null(hydrofabric:::.order_parts_by_reach(parts, one, NULL))
  expect_null(hydrofabric:::.order_parts_by_reach(parts[1, ], one, far))
})
