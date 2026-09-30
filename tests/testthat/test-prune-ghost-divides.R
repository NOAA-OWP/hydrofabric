# A divide too small to be a catchment must be pruned wherever the prune runs,
# including after a step that has just rewritten geometry. Overlap removal
# resolves two overlapping divides by awarding the ground to one of them, which
# can leave the other as a scatter of near-zero parts; deciding from a stored
# area column instead of the polygon lets that residue reach the product with
# has_divide = TRUE, where it fails divide_has_member_flowline because its
# flowlines are kilometres from what is left of it.
suppressMessages({ library(sf) })

.sq <- function(x0, y0, w, h)
  st_polygon(list(cbind(c(x0, x0 + w, x0 + w, x0, x0),
                        c(y0, y0, y0 + h, y0 + h, y0))))

# host 100x100 m (0.01 km2, well above the floor) beside a 10x10 m residue
# (1e-4 km2, below it) that still CLAIMS 5 km2 in its areasqkm column.
.gpkg_with_ghost <- function(stated_ghost_area = 5) {
  g <- st_sfc(.sq(0, 0, 100, 100), .sq(100, 0, 10, 10), crs = 5070)
  d <- st_sf(divide_id   = c("cat-host", "cat-ghost"),
             flowpath_id = c("fp-host", "fp-ghost"),
             type        = c("network", "network"),
             areasqkm    = c(0.01, stated_ghost_area),
             geom        = g)
  fl <- st_sfc(st_linestring(cbind(c(10, 90), c(10, 90))),
               st_linestring(cbind(c(101, 109), c(1, 9))), crs = 5070)
  fp <- st_sf(flowpath_id = c("fp-host", "fp-ghost"),
              divide_id   = c("cat-host", "cat-ghost"),
              has_divide  = c(TRUE, TRUE),
              areasqkm    = c(0.01, stated_ghost_area),
              geom        = fl)
  p <- tempfile(fileext = ".gpkg")
  st_write(d,  p, layer = "divides",   quiet = TRUE)
  st_write(fp, p, layer = "flowpaths", quiet = TRUE)
  p
}

test_that("a collapsed divide is pruned even while its area column claims otherwise", {
  p <- .gpkg_with_ghost(stated_ghost_area = 5)
  pruned <- suppressMessages(hydrofabric:::prune_ghost_divides(p))
  expect_identical(as.character(pruned), "cat-ghost")

  d <- st_read(p, "divides", quiet = TRUE)
  expect_false("cat-ghost" %in% as.character(d$divide_id))
})

test_that("the pruned divide's flowpath becomes a connector", {
  p <- .gpkg_with_ghost()
  suppressMessages(hydrofabric:::prune_ghost_divides(p))
  fp <- st_read(p, "flowpaths", quiet = TRUE)
  ghost <- fp[fp$flowpath_id == "fp-ghost", ]
  expect_false(as.logical(ghost$has_divide))
  expect_true(is.na(ghost$divide_id))
  expect_equal(ghost$areasqkm, 0)
  # the untouched flowpath keeps its divide
  host <- fp[fp$flowpath_id == "fp-host", ]
  expect_identical(as.character(host$divide_id), "cat-host")
})

test_that("pruning conserves ground: the residue is absorbed, not deleted", {
  p <- .gpkg_with_ghost()
  before <- as.numeric(st_area(st_union(st_geometry(st_read(p, "divides", quiet = TRUE)))))
  suppressMessages(hydrofabric:::prune_ghost_divides(p))
  after <- as.numeric(st_area(st_union(st_geometry(st_read(p, "divides", quiet = TRUE)))))
  expect_equal(after, before, tolerance = 1e-6)
})

test_that("a divide above the floor is left alone even if its column understates it", {
  # The mirror of the first case: trusting the column would prune a real
  # catchment because some earlier step recorded it as empty.
  p <- .gpkg_with_ghost(stated_ghost_area = 5)
  d <- st_read(p, "divides", quiet = TRUE)
  d$areasqkm <- c(0, 0)                      # both columns claim nothing
  st_write(d, p, layer = "divides", delete_layer = TRUE, quiet = TRUE)
  pruned <- suppressMessages(hydrofabric:::prune_ghost_divides(p))
  expect_false("cat-host" %in% as.character(pruned))
  expect_true("cat-host" %in%
    as.character(st_read(p, "divides", quiet = TRUE)$divide_id))
})
