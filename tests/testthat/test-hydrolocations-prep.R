# Two vertical flowpaths at x=0 and x=1; points (0,0) and (1,1) snap to them.
.fp_fixture <- function() {
  a <- sf::st_linestring(rbind(c(0, -1), c(0, 2)))
  b <- sf::st_linestring(rbind(c(1, -1), c(1, 2)))
  sf::st_sf(flowpath_id = c("A", "B"), vpuid = c("01", "01"),
            geometry = sf::st_sfc(a, b, crs = 5070))
}

test_that("prepare_hydrolocations maps columns, snaps, and assigns ids", {
  skip_if_not_installed("sf")
  pts <- sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(1, 1)), crs = 5070)
  raw <- sf::st_sf(gage = c("01030500", "01030350"), net = c("usgs", "usgs"),
                   geometry = pts)
  hl <- prepare_hydrolocations(raw, flowpaths = .fp_fixture(),
    mapping = c(hl_reference = "gage", hl_source = "net", hl_class = "gage"))
  expect_true(all(c("flowpath_id", "poi_id", "hl_id", "hl_reference",
                    "pathmeas", "vpuid") %in% names(hl)))
  expect_setequal(hl$flowpath_id, c("A", "B"))          # snapped to nearest
  expect_setequal(hl$hl_reference, c("01030500", "01030350"))
  expect_equal(length(unique(hl$poi_id)), 2L)           # one poi_id per flowpath
  expect_false(anyNA(hl$hl_id))
})

test_that("prepare_hydrolocations drops rows without hl_reference (warns)", {
  skip_if_not_installed("sf")
  pts <- sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(1, 1)), crs = 5070)
  raw <- sf::st_sf(gage = c("01030500", NA), geometry = pts)
  expect_warning(hl <- prepare_hydrolocations(raw, flowpaths = .fp_fixture(),
                   mapping = c(hl_reference = "gage")), "dropping 1")
  expect_equal(nrow(hl), 1L)
})

test_that("prepare_hydrolocations drops POIs beyond max_dist (warns)", {
  skip_if_not_installed("sf")
  pts <- sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(5000, 5000)), crs = 5070)
  raw <- sf::st_sf(gage = c("near", "far"), geometry = pts)
  expect_warning(hl <- prepare_hydrolocations(raw, flowpaths = .fp_fixture(),
                   mapping = c(hl_reference = "gage"), max_dist = 1000),
                 "> 1000 m")
  expect_equal(hl$hl_reference, "near")
})

test_that("prepare_hydrolocations honours a credible source link over the nearest reach", {
  skip_if_not_installed("sf")
  # The point sits ON A, but its link names B 1 m away: well inside max_link_m,
  # so the link is authoritative about which reach it belongs to.
  raw <- sf::st_sf(gage = "01030500", flowpath_id = "B",
                   geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 5070))
  hl <- prepare_hydrolocations(raw, flowpaths = .fp_fixture(),
    mapping = c(hl_reference = "gage"))
  expect_equal(hl$flowpath_id, "B")
})

test_that("prepare_hydrolocations distrusts a source link that is too far away", {
  skip_if_not_installed("sf")
  # `flowpath_id` is not unique across reference domains, so a link can resolve
  # onto an unrelated reach by id collision. The separation between the raw point
  # and the reach its link names is what catches that, so the guard has to be
  # evaluated BEFORE the point is snapped onto that reach, or it reads zero by
  # construction and never fires.
  fp <- .fp_fixture()
  far <- sf::st_sf(flowpath_id = "C", vpuid = "01",
    geometry = sf::st_sfc(sf::st_linestring(rbind(c(10000, -1), c(10000, 2))), crs = 5070))
  fp <- rbind(fp, far)
  raw <- sf::st_sf(gage = "01030500", flowpath_id = "C",
                   geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 5070))
  expect_message(hl <- prepare_hydrolocations(raw, flowpaths = fp,
    mapping = c(hl_reference = "gage")), "more than 250m")
  expect_equal(hl$flowpath_id, "A")                       # fell back to nearest
  # and the geometry must not have been dragged 10 km onto the distrusted reach
  expect_lt(sf::st_coordinates(hl)[1, 1], 1)
})

test_that("prepare_hydrolocations errors on missing mapping / non-point / bad source", {
  skip_if_not_installed("sf")
  pts <- sf::st_sfc(sf::st_point(c(0, 0)), crs = 5070)
  raw <- sf::st_sf(gage = "01030500", geometry = pts)
  expect_error(prepare_hydrolocations(raw, flowpaths = .fp_fixture()), "hl_reference")
  expect_error(prepare_hydrolocations(raw, flowpaths = .fp_fixture(),
               mapping = c(hl_reference = "nope")), "not found")
  poly <- sf::st_sf(hl_reference = "x",
    geometry = sf::st_sfc(sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,0)))), crs = 5070))
  expect_error(prepare_hydrolocations(poly, flowpaths = .fp_fixture()), "POINT")
})

# Reaches carrying has_divide / totdasqkm, so the empty-nexus rule is testable.
# C has no local catchment but routes a large basin; D has neither.
.fp_da_fixture <- function() {
  mk <- function(x) sf::st_linestring(rbind(c(x, -1), c(x, 2)))
  sf::st_sf(flowpath_id = c("A", "C", "D"), vpuid = rep("01", 3),
            has_divide  = c(1L, 0L, 0L),
            totdasqkm   = c(12.5, 9178.96, 0),
            geometry = sf::st_sfc(mk(0), mk(1), mk(2), crs = 5070))
}

test_that("a POI that could only become an inflow-less nexus is dropped", {
  skip_if_not_installed("sf")
  pts <- sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(1, 0)),
                    sf::st_point(c(2, 0)), crs = 5070)
  raw <- sf::st_sf(gage = c("keep-catchment", "keep-routing", "drop-dead"),
                   type = c("nwis", "nwis", "wbout"), geometry = pts)
  hl <- suppressMessages(prepare_hydrolocations(raw,
    flowpaths = .fp_da_fixture(),
    mapping = c(hl_reference = "gage", hl_type = "type")))
  expect_equal(nrow(hl), 2L)
  expect_setequal(hl$hl_reference, c("keep-catchment", "keep-routing"))
  # the routing reach is kept precisely because has_divide == 0 alone is not
  # grounds to drop: it carries 9179 km2 from upstream
  expect_true("C" %in% hl$flowpath_id)
  expect_false("D" %in% hl$flowpath_id)
})

test_that("the rule is inert when the reference lacks the columns", {
  skip_if_not_installed("sf")
  pts <- sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(1, 1)), crs = 5070)
  raw <- sf::st_sf(gage = c("a", "b"), geometry = pts)
  hl <- prepare_hydrolocations(raw, flowpaths = .fp_fixture(),
    mapping = c(hl_reference = "gage"))
  expect_equal(nrow(hl), 2L)
})

test_that("the reach each POI is measured against survives the row sort", {
  skip_if_not_installed("sf")
  # hl_reference order is deliberately the reverse of flowpath_id order, so the
  # sort by flowpath_id displaces every row. If `idx` does not travel with the
  # rows, the along-reach position is read from another row's reach and the
  # collocation grouping below is computed from it.
  mk <- function(x) sf::st_linestring(rbind(c(x, -1), c(x, 2)))
  fp <- sf::st_sf(flowpath_id = c("A", "B", "C"), vpuid = rep("01", 3),
                  geometry = sf::st_sfc(mk(0), mk(10), mk(20), crs = 5070))
  pts <- sf::st_sfc(sf::st_point(c(20, 0)), sf::st_point(c(10, 0)),
                    sf::st_point(c(0, 0)), crs = 5070)
  raw <- sf::st_sf(gage = c("z", "y", "x"), geometry = pts)
  hl <- prepare_hydrolocations(raw, flowpaths = fp,
    mapping = c(hl_reference = "gage"))
  # each POI must end up on the reach nearest its own coordinate
  expect_identical(hl$flowpath_id[hl$hl_reference == "x"], "A")
  expect_identical(hl$flowpath_id[hl$hl_reference == "y"], "B")
  expect_identical(hl$flowpath_id[hl$hl_reference == "z"], "C")
  # three separate reaches -> three distinct POIs, never collapsed
  expect_equal(length(unique(hl$poi_id)), 3L)
})
