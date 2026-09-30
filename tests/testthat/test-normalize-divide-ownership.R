# A divide must carry one of its flowpath's flowlines, the rule
# `divide_has_member_flowline` states. This enforces it once geometry is final,
# for domains whose grid and vector network come from unrelated sources and
# cannot be split exactly.
#
# The property that matters most is the one asserted first: on a fabric that
# already conforms, this must do nothing at all.
suppressMessages({ library(sf) })

.sq <- function(x0, y0, w, h)
  st_polygon(list(cbind(c(x0, x0 + w, x0 + w, x0, x0),
                        c(y0, y0, y0 + h, y0 + h, y0))))
.ln <- function(x0, y0, x1, y1) st_linestring(cbind(c(x0, x1), c(y0, y1)))

# Two adjacent divides, two flowpaths, one flowline each. `swap` puts each
# flowline in the OTHER divide, which is the mispairing this repairs.
#
# The divides are 1 km across and the lines sit 100 m inside them, so a line in
# the wrong divide is well beyond the 30 m tolerance. At 100 m squares the
# "wrong" line lands 10 m from its neighbour and the gate correctly passes it,
# adjacency alone is not a mispairing.
.fab <- function(swap = FALSE, orphan = FALSE) {
  g <- st_sfc(.sq(0, 0, 1000, 1000), .sq(1000, 0, 1000, 1000), crs = 5070)
  d <- st_sf(divide_id = c("cat-1", "cat-2"), flowpath_id = c("fp-1", "fp-2"),
             type = c("network", "network"), areasqkm = c(1, 1), geom = g)
  # lines well inside cat-1 and cat-2 respectively
  l1 <- .ln(100, 100, 900, 900); l2 <- .ln(1100, 100, 1900, 900)
  owners <- if (swap) c("fp-2", "fp-1") else c("fp-1", "fp-2")
  fl <- st_sf(flowpath_id = owners, geom = st_sfc(l1, l2, crs = 5070))
  fp <- st_sf(flowpath_id = c("fp-1", "fp-2"), divide_id = c("cat-1", "cat-2"),
              has_divide = c(TRUE, TRUE), areasqkm = c(1, 1),
              geom = st_sfc(l1, l2, crs = 5070))

  if (orphan) {
    # A third divide next door whose owner DOES have a flowline, it just sits
    # 2 km away, inside cat-1. Nothing runs through cat-3, so there is no
    # flowpath to re-own it to and it has to be absorbed.
    #
    # The owner must have flowlines somewhere: the gate treats a flowpath with
    # none as untestable and skips it, so a divide owned by such a flowpath
    # never fails and is deliberately left alone.
    l3 <- .ln(200, 800, 800, 200)                    # inside cat-1
    d  <- rbind(d, st_sf(divide_id = "cat-3", flowpath_id = "fp-3",
                         type = "network", areasqkm = 1,
                         geom = st_sfc(.sq(2000, 0, 1000, 1000), crs = 5070)))
    fl <- rbind(fl, st_sf(flowpath_id = "fp-3", geom = st_sfc(l3, crs = 5070)))
    fp <- rbind(fp, st_sf(flowpath_id = "fp-3", divide_id = "cat-3",
                          has_divide = TRUE, areasqkm = 1,
                          geom = st_sfc(l3, crs = 5070)))
  }
  p <- tempfile(fileext = ".gpkg")
  st_write(d,  p, layer = "divides",   quiet = TRUE)
  st_write(fp, p, layer = "flowpaths", quiet = TRUE)
  st_write(fl, p, layer = "flowlines", quiet = TRUE)
  p
}

test_that("a conforming fabric is left completely untouched", {
  p <- .fab()
  before_d  <- st_read(p, "divides",   quiet = TRUE)
  before_fp <- st_read(p, "flowpaths", quiet = TRUE)
  r <- suppressMessages(hydrofabric:::normalize_divide_ownership(p))
  expect_identical(nrow(r), 0L)
  expect_identical(st_read(p, "divides",   quiet = TRUE)$flowpath_id,
                   before_d$flowpath_id)
  expect_identical(st_read(p, "flowpaths", quiet = TRUE)$divide_id,
                   before_fp$divide_id)
})

test_that("a swapped pairing is re-owned to the flowpath whose line it holds", {
  p <- .fab(swap = TRUE)
  r <- suppressMessages(hydrofabric:::normalize_divide_ownership(p))
  expect_identical(nrow(r), 2L)
  expect_true(all(r$action == "reowned"))
  d <- st_read(p, "divides", quiet = TRUE)
  # cat-1 holds fp-2's line after the swap, so cat-1 must now belong to fp-2
  expect_identical(as.character(d$flowpath_id[d$divide_id == "cat-1"]), "fp-2")
  expect_identical(as.character(d$flowpath_id[d$divide_id == "cat-2"]), "fp-1")
  # and the flowpaths agree
  fp <- st_read(p, "flowpaths", quiet = TRUE)
  expect_identical(as.character(fp$divide_id[fp$flowpath_id == "fp-2"]), "cat-1")
})

test_that("re-owning conserves ground: no divide geometry moves", {
  p <- .fab(swap = TRUE)
  before <- as.numeric(st_area(st_union(st_geometry(st_read(p, "divides", quiet = TRUE)))))
  suppressMessages(hydrofabric:::normalize_divide_ownership(p))
  after <- as.numeric(st_area(st_union(st_geometry(st_read(p, "divides", quiet = TRUE)))))
  expect_equal(after, before, tolerance = 1e-6)
})

test_that("a divide nobody's line touches is absorbed and its owner demoted", {
  p <- .fab(orphan = TRUE)   # cat-3 holds no line; its owner's line is 2 km away
  before <- as.numeric(st_area(st_union(st_geometry(st_read(p, "divides", quiet = TRUE)))))
  r <- suppressMessages(hydrofabric:::normalize_divide_ownership(p))
  expect_true("absorbed" %in% r$action)
  d  <- st_read(p, "divides",   quiet = TRUE)
  fp <- st_read(p, "flowpaths", quiet = TRUE)
  expect_false("cat-3" %in% as.character(d$divide_id))
  expect_true(is.na(fp$divide_id[fp$flowpath_id == "fp-3"]))
  expect_false(as.logical(fp$has_divide[fp$flowpath_id == "fp-3"]))
  # the ground stays in the layer
  after <- as.numeric(st_area(st_union(st_geometry(d))))
  expect_equal(after, before, tolerance = 1e-6)
})

test_that("a missing flowlines layer is a skip, not an error", {
  p <- .fab()
  con <- DBI::dbConnect(RSQLite::SQLite(), p)
  DBI::dbExecute(con, "DROP TABLE flowlines")
  DBI::dbExecute(con, "DELETE FROM gpkg_contents WHERE table_name='flowlines'")
  DBI::dbDisconnect(con)
  expect_silent(suppressMessages(hydrofabric:::normalize_divide_ownership(p)))
})
