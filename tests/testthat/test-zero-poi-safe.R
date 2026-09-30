# A reference fabric can legitimately carry no POIs: most TDX-Hydro continental
# tiles do, and the pipeline then finds no `pois` layer where it expects one.
# Unguarded reads and joins against an absent or empty layer used to abort the
# run rather than skip the work.
#
# These tests cover the two shapes that matter and are easy to conflate: the
# layer being ABSENT, and the layer being PRESENT BUT EMPTY.

skip_if_not_installed("hfutils")

make_gpkg <- function(pois = c("absent", "empty", "populated"),
                      hydrolocations = TRUE) {
  pois <- match.arg(pois)
  f <- tempfile(fileext = ".gpkg")
  pt <- function(n) sf::st_sfc(lapply(seq_len(max(n, 1)), function(i) sf::st_point(c(i, 0))),
                               crs = 5070)
  fp <- sf::st_sf(flowpath_id = c("fp-1", "fp-2"), geom = pt(2))
  sf::st_write(fp, f, "flowpaths", quiet = TRUE)
  if (hydrolocations) {
    hl <- sf::st_sf(poi_id = c(1L, 2L), hl_reference = c("a", "b"),
                    hl_type = c("nwis", "wbout"), hl_class = NA_character_,
                    hl_link = c("x", "y"), geom = pt(2))
    sf::st_write(hl, f, "hydrolocations", quiet = TRUE)
  }
  if (pois != "absent") {
    p <- sf::st_sf(poi_id = c(1L, 2L), flowpath_id = c("fp-1", "fp-2"), geom = pt(2))
    if (pois == "empty") p <- p[0, ]
    sf::st_write(p, f, "pois", quiet = TRUE)
  }
  f
}

test_that("rollups skip a fabric with no pois layer at all", {
  f <- make_gpkg(pois = "absent")
  expect_null(hydrofabric:::.add_poi_rollups(f))
})

test_that("rollups skip a pois layer that exists but is empty", {
  # the case the layer_exists() guard alone does not catch
  f <- make_gpkg(pois = "empty")
  expect_null(hydrofabric:::.add_poi_rollups(f))
})

test_that("rollups skip when hydrolocations are absent", {
  f <- make_gpkg(pois = "populated", hydrolocations = FALSE)
  expect_null(hydrofabric:::.add_poi_rollups(f))
})

test_that("rollups still run when both layers are populated", {
  f <- make_gpkg(pois = "populated")
  expect_silent(hydrofabric:::.add_poi_rollups(f))
  p <- sf::st_read(f, "pois", quiet = TRUE)
  expect_true(any(grepl("^hl_", names(p))))
})

test_that("the zero-POI guards are present in the installed code", {
  agg <- paste(deparse(getFromNamespace("aggregate_fabric", asNamespace("hydrofabric"))),
               collapse = "\n")
  expect_true(grepl('layer_exists(s1_pois, "pois")', agg, fixed = TRUE))
  fin <- paste(deparse(getFromNamespace(".add_poi_rollups", asNamespace("hydrofabric"))),
               collapse = "\n")
  expect_true(grepl("nrow(pois) == 0", fin, fixed = TRUE))
  pro <- paste(deparse(getFromNamespace(".promote_hydrolocations", asNamespace("hydrofabric"))),
               collapse = "\n")
  expect_true(grepl('layer_exists(ngen_gpkg, "pois")', pro, fixed = TRUE))
})
