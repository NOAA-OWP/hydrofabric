test_that("construct_agg_divides emits one row per aggregate, even when non-contiguous", {
  skip_if_not(requireNamespace("sf", quietly = TRUE))
  # Regression: the function used to cast to POLYGON, which explodes a multi-part
  # divide into one row per part with the SAME divide_id. The layer then left the
  # function with more rows than aggregates (+23 on VPU 13), and Step 4d plus Step
  # 6d re-merged them: two repairs undoing one cast.
  #
  # A divide is legitimately multi-part when its members' catchments are not
  # contiguous. That must stay one row, keep both parts, and lose no area.
  sq <- function(x0, w = 1000) sf::st_polygon(list(cbind(
    c(x0, x0 + w, x0 + w, x0, x0), c(0, 0, w, w, 0))))

  rec_div <- sf::st_sf(
    flowpath_id = c("r1", "r2"),
    geometry    = sf::st_sfc(sq(0), sq(5000), crs = 5070))   # deliberately far apart
  network <- data.frame(flowpath_id = "A", rec_id_list = "r1,r2",
    stringsAsFactors = FALSE)
  fp <- sf::st_sf(flowpath_id = "A", flowpath_toid = "0",
    geometry = sf::st_sfc(
      sf::st_linestring(cbind(c(0, 900), c(500, 500))), crs = 5070))

  d <- suppressMessages(construct_agg_divides(network, rec_div, fp))

  expect_identical(nrow(d), 1L)
  expect_identical(length(unique(as.character(d$divide_id))), 1L)
  expect_true(all(as.character(sf::st_geometry_type(d)) == "MULTIPOLYGON"))
  # Both 1 km2 squares survive: the discard the cast enabled is gone.
  expect_equal(sum(d$areasqkm), 2, tolerance = 1e-6)
  expect_identical(length(sf::st_geometry(d)[[1]]), 2L)

  # Contiguous members must dissolve to a single part, not stay separate.
  rec_div2 <- sf::st_sf(
    flowpath_id = c("r1", "r2"),
    geometry    = sf::st_sfc(sq(0), sq(1000), crs = 5070))    # share an edge
  d2 <- suppressMessages(construct_agg_divides(network, rec_div2, fp))
  expect_identical(nrow(d2), 1L)
  expect_equal(sum(d2$areasqkm), 2, tolerance = 1e-6)
})

test_that("rec_id_partition and aggregated_toid_resolves replace the repairs they retired", {
  # Steps 4_dedup_recid, 4_toid and 4_selfloop were removed after reporting no-ops on
  # all 24 VPUs. Their properties are load-bearing: construct_agg_divides dissolves
  # rec_div by member, so a reach claimed twice double-counts its catchment, so they
  # are asserted rather than left unchecked.
  net <- data.frame(flowpath_id = c("1", "2"), rec_id_list = c("a,b", "c"),
    stringsAsFactors = FALSE)
  fp  <- data.frame(flowpath_id = c("1", "2"), flowpath_toid = c("2", "0"),
    stringsAsFactors = FALSE)

  ok <- hydrofabric:::.hf_check_aggregated(
    list(flowpaths = fp, network = net), FALSE, 0)
  expect_true(ok$rec_id_partition$ok)
  expect_true(ok$aggregated_toid_resolves$ok)
  expect_identical(ok$rec_id_partition$kind, "check")

  # A reach claimed by two aggregates must fail.
  dup <- hydrofabric:::.hf_check_aggregated(
    list(flowpaths = fp,
      network = data.frame(flowpath_id = c("1", "2"), rec_id_list = c("a,b", "b"),
        stringsAsFactors = FALSE)), FALSE, 0)
  expect_false(dup$rec_id_partition$ok)
  expect_match(dup$rec_id_partition$msg, "^1 reconciled reach")

  # A dangling toid and a self-loop must both fail.
  bad <- hydrofabric:::.hf_check_aggregated(
    list(flowpaths = data.frame(flowpath_id = c("1", "2"),
      flowpath_toid = c("9", "2"), stringsAsFactors = FALSE),
      network = net), FALSE, 0)
  expect_false(bad$aggregated_toid_resolves$ok)
  expect_match(bad$aggregated_toid_resolves$msg, "1 dangling.*1 self-loop")

  # A terminal toid (>= 1e9 nexus id or "0") is not dangling.
  term <- hydrofabric:::.hf_check_aggregated(
    list(flowpaths = data.frame(flowpath_id = c("1", "2"),
      flowpath_toid = c("2", "1000000001"), stringsAsFactors = FALSE),
      network = net), FALSE, 0)
  expect_true(term$aggregated_toid_resolves$ok)
})
