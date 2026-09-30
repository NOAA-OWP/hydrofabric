test_that("id_types_not_double catches floating-point identifiers and clears once retyped", {
  skip_if_not_installed("sf")
  skip_if_not_installed("RSQLite")

  tmp <- tempfile(fileext = ".gpkg")
  on.exit(unlink(tmp), add = TRUE)

  pts <- sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(1, 1)), crs = 4326)
  sf::st_write(sf::st_sf(
    flowline_id   = c(101, 102),        # REAL identifier -- the defect
    mainstem_id   = c(7, 8),            # REAL identifier -- the defect
    flowpath_id   = c("fp-1", "fp-2"),  # prefixed, legitimately TEXT
    lengthkm      = c(1.5, 2.25),       # a measure, not an id
    geom = pts), tmp, "flowlines", quiet = TRUE)

  chk <- hydrofabric:::.hf_id_types_not_double(tmp)
  expect_false(chk$ok)
  expect_match(chk$msg, "flowlines\\.flowline_id")
  expect_match(chk$msg, "flowlines\\.mainstem_id")

  enforce_integer_ids(tmp)

  chk2 <- hydrofabric:::.hf_id_types_not_double(tmp)
  expect_true(chk2$ok)
  expect_match(chk2$msg, "^0 id column")
})

test_that("id_types_not_double does not flag measures or prefixed string ids", {
  skip_if_not_installed("sf")
  skip_if_not_installed("RSQLite")

  tmp <- tempfile(fileext = ".gpkg")
  on.exit(unlink(tmp), add = TRUE)

  pts <- sf::st_sfc(sf::st_point(c(0, 0)), crs = 4326)
  sf::st_write(sf::st_sf(divide_id = "cat-1", nexus_id = "nex-1",
    areasqkm = 3.7, geom = pts), tmp, "divides", quiet = TRUE)

  chk <- hydrofabric:::.hf_id_types_not_double(tmp)
  expect_true(chk$ok)
})
