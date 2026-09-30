# check_reference_severed_terminals, a reach with no downstream that ends on a
# channel carrying more water than it does. The attributes cannot answer this:
# terminalfl is unpopulated for whole sources and dnhydroseq is stored too
# narrowly to hold its own values. Geometry can. On the CONUS reference 446 of
# 15,877 terminals are severed this way, the largest a 324,094 km2 river.

.cst_gpkg <- function(ids, toids, areas, lines, vpuid = NULL) {
  d <- sf::st_sf(flowpath_id = as.character(ids),
                 flowpath_toid = as.character(toids),
                 totdasqkm = as.numeric(areas),
                 geometry = sf::st_sfc(lines, crs = 5070))
  if (!is.null(vpuid)) d$vpuid <- vpuid
  p <- tempfile(fileext = ".gpkg")
  suppressWarnings(sf::st_write(d, p, "flowpaths", quiet = TRUE))
  p
}
.cst_ln <- function(...) sf::st_linestring(rbind(...))

test_that("a terminal ending on a larger channel is reported", {
  p <- .cst_gpkg(c("main", "trib"), c("0", "0"), c(500, 20),
    list(.cst_ln(c(0, 0), c(10000, 0)), .cst_ln(c(5000, 3000), c(5000, 0))))
  r <- suppressMessages(hydrofabric::check_reference_severed_terminals(p))
  expect_equal(r$id, "trib")
  expect_equal(r$continues_into, "main")
})

test_that("a true outlet is not reported", {
  p <- .cst_gpkg(c("main", "trib"), c("0", "main"), c(500, 20),
    list(.cst_ln(c(0, 0), c(10000, 0)), .cst_ln(c(5000, 3000), c(5000, 0))))
  r <- suppressMessages(hydrofabric::check_reference_severed_terminals(p))
  expect_equal(nrow(r), 0L)
})

test_that("a neighbour draining less is not a receiver", {
  # The big river ends on a stub. Routing it into the stub would be backwards.
  p <- .cst_gpkg(c("big", "stub"), c("0", "0"), c(500, 2),
    list(.cst_ln(c(0, 0), c(5000, 0)), .cst_ln(c(5000, -200), c(5000, 200))))
  r <- suppressMessages(hydrofabric::check_reference_severed_terminals(p))
  expect_equal(nrow(r), 0L)
})

test_that("tolerance is honoured", {
  p <- .cst_gpkg(c("main", "trib"), c("0", "0"), c(500, 20),
    list(.cst_ln(c(0, 0), c(10000, 0)), .cst_ln(c(5000, 3000), c(5000, 100))))
  expect_equal(nrow(suppressMessages(
    hydrofabric::check_reference_severed_terminals(p))), 0L)
  expect_equal(suppressMessages(
    hydrofabric::check_reference_severed_terminals(p, tolerance_m = 150))$id,
    "trib")
})

test_that("the largest receiver wins when several are in reach", {
  # 'small' is itself terminal and also ends beside 'big', so it is reported
  # too: correctly. The assertion is about which receiver 'trib' is given.
  p <- .cst_gpkg(c("small", "big", "trib"), c("0", "0", "0"), c(50, 900, 20),
    list(.cst_ln(c(0, 10), c(10000, 10)), .cst_ln(c(0, -10), c(10000, -10)),
         .cst_ln(c(5000, 3000), c(5000, 0))))
  r <- suppressMessages(hydrofabric::check_reference_severed_terminals(p))
  expect_equal(r$continues_into[r$id == "trib"], "big")
  expect_setequal(r$id, c("small", "trib"))
})

test_that("results come back largest first", {
  p <- .cst_gpkg(c("main", "a", "b"), c("0", "0", "0"), c(900, 20, 100),
    list(.cst_ln(c(0, 0), c(10000, 0)), .cst_ln(c(3000, 3000), c(3000, 0)),
         .cst_ln(c(7000, 3000), c(7000, 0))))
  r <- suppressMessages(hydrofabric::check_reference_severed_terminals(p))
  expect_equal(r$id, c("b", "a"))
})

test_that("a vpuid filter restricts the population", {
  p <- .cst_gpkg(c("main", "trib"), c("0", "0"), c(500, 20),
    list(.cst_ln(c(0, 0), c(10000, 0)), .cst_ln(c(5000, 3000), c(5000, 0))),
    vpuid = c("04", "04"))
  expect_equal(suppressMessages(
    hydrofabric::check_reference_severed_terminals(p, vpuid = "04"))$id, "trib")
  expect_equal(nrow(suppressMessages(
    hydrofabric::check_reference_severed_terminals(p, vpuid = "18"))), 0L)
})

test_that("a layer with no terminals is skipped", {
  p <- .cst_gpkg(c("a", "b"), c("b", "a"), c(10, 20),
    list(.cst_ln(c(0, 0), c(100, 0)), .cst_ln(c(100, 0), c(200, 0))))
  expect_equal(nrow(suppressMessages(
    hydrofabric::check_reference_severed_terminals(p))), 0L)
})
