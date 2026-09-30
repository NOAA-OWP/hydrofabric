## End-to-end regression test for the raster-free aggregate phase (aggregate_fabric).
##
## Runs aggregate_fabric() on a committed minimal refactor-phase fixture (the smallest
## gage subset, tests/fixtures/stage1_min) and asserts the package's own
## aggregated-stage invariants pass plus a few stable structural facts. This
## exercises the real aggregation/divide-construction code path; assertions are
## contract-based (invariants + topology) rather than brittle geometry byte
## comparisons.
##
## Skipped where mapshaper is unavailable (the aggregate divide cleanup shells out
## to it via hfutils::clean_gpkg_layer).

test_that("aggregate_fabric produces an invariant-clean aggregated fabric", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.which("mapshaper")), "mapshaper not installed")

  stage1_dir <- testthat::test_path("..", "fixtures", "stage1_min", "vpu_01")
  skip_if_not(dir.exists(stage1_dir), "stage1_min fixture missing")

  out_dir <- withr::local_tempdir()
  res <- suppressWarnings(suppressMessages(
    aggregate_fabric(vpuid = "01", fabric_dir = stage1_dir,
      out_dir = out_dir, verbose = FALSE)))

  # The pipeline's own contract: aggregated-stage invariants must pass.
  expect_true(isTRUE(res$invariants$ok))

  gpkg <- file.path(out_dir, "vpu_01", "work", "out_vpu01.gpkg")
  expect_true(file.exists(gpkg))
  expect_setequal(sf::st_layers(gpkg)$name,
    c("divides", "flowpaths", "network", "pois", "small_segments"))

  fp  <- sf::read_sf(gpkg, "flowpaths", quiet = TRUE)
  dv  <- sf::read_sf(gpkg, "divides",   quiet = TRUE)
  net <- sf::read_sf(gpkg, "network",   quiet = TRUE)

  # Structural contracts
  expect_gt(nrow(fp), 0L)
  expect_equal(nrow(dv), nrow(fp))                       # one catchment per flowpath
  expect_true(hfutils::hf_network_is_dag(fp))                    # dendritic topology
  expect_false(anyNA(as.character(fp$flowpath_id)))

  # Regression snapshot for this fixture (current pipeline). A change here means
  # the aggregation produced a different network/divide count: investigate.
  #
  # 4 catchments, from 2 in the delivered build. The fixture holds 32.19 km2
  # and the count alone says nothing, the areas say why 4 is right:
  #
  #   delivered      5.35, 26.84                 one unit at 2.7x the target
  #   target only    2.32, 2.91, ... 9.51        nothing over 10, but two under
  #                                              the floor with nowhere to go
  #   three sizes    3.87, 5.35, 9.51, 13.46     floor cleared, ceiling kept
  #
  # Total area is 32.19 km2 in every case: the same land, in different units.
  # The 13.46 unit is what absorbing a sub-floor neighbour costs, and it is
  # allowed precisely because it lands under the 15 km2 ceiling.
  expect_equal(nrow(fp),  4L)
  expect_equal(nrow(dv),  4L)
  expect_equal(nrow(net), 13L)

  # The two properties the count is a proxy for, and neither may regress:
  # nothing is left below the floor, and nothing is built past the ceiling.
  ar <- suppressWarnings(as.numeric(dv$areasqkm))
  expect_true(all(ar >= 3))
  expect_true(all(ar <= 15))
  expect_equal(round(sum(ar), 2), 32.19)
})
