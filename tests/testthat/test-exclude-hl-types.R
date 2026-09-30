# Hydrolocation filtering at ingest.
#
# The OWP-specified set carries waterbody inlets, about half of it, which are
# not consumed as points of interest. Each one pins its own unit against merging
# and breaks the levelpath, stranding the unit below it too, so they dominate
# the sub-minimum population while naming nothing a consumer asks for.

hl <- function(poi, type) data.frame(
  poi_id = as.character(poi), hl_type = type, flowpath_id = paste0("fp-", poi),
  stringsAsFactors = FALSE)

test_that("rows of an excluded type are dropped", {
  l <- list(hydrolocations = hl(c(1, 2, 3), c("wbin", "nwis", "wbin")),
            pois = NULL)
  r <- hydrofabric:::.drop_hl_types(l, "wbin", "01")
  expect_equal(nrow(r$hydrolocations), 1L)
  expect_equal(r$hydrolocations$poi_id, "2")
})

test_that("a location keeps a wanted type recorded alongside an excluded one", {
  # Filtering is by row, not by location: poi 1 is recorded twice and survives
  # on its gage row.
  l <- list(hydrolocations = rbind(hl(1, "wbin"), hl(1, "nwis"), hl(2, "wbin")),
            pois = NULL)
  r <- hydrofabric:::.drop_hl_types(l, "wbin", "01")
  expect_equal(nrow(r$hydrolocations), 1L)
  expect_equal(r$hydrolocations$hl_type, "nwis")
})

test_that("an all-NA poi_id does not take the pois layer with it", {
  # poi_id is assigned by .normalize_poi_id straight after this, so the input
  # carries none. Matching the pois layer against NA would empty it.
  l <- list(hydrolocations = data.frame(poi_id = NA_character_,
              hl_type = c("wbin", "nwis"), stringsAsFactors = FALSE),
            pois = data.frame(poi_id = c("1", "2"), stringsAsFactors = FALSE))
  r <- hydrofabric:::.drop_hl_types(l, "wbin", "01")
  expect_equal(nrow(r$hydrolocations), 1L)
  expect_equal(nrow(r$pois), 2L)
})

test_that("the pois layer is filtered to the same surviving locations", {
  # The two layers are renumbered together straight after, so a poi_id left in
  # one but not the other would desync them.
  l <- list(hydrolocations = rbind(hl(1, "wbin"), hl(2, "nwis")),
            pois = data.frame(poi_id = c("1", "2"), stringsAsFactors = FALSE))
  r <- hydrofabric:::.drop_hl_types(l, "wbin", "01")
  expect_equal(r$pois$poi_id, "2")
})

test_that("nothing changes when no row matches, or the layer is untyped", {
  l <- list(hydrolocations = hl(c(1, 2), c("nwis", "rfc")), pois = NULL)
  expect_identical(hydrofabric:::.drop_hl_types(l, "wbin", "01"), l)
  u <- list(hydrolocations = data.frame(poi_id = "1"), pois = NULL)
  expect_identical(hydrofabric:::.drop_hl_types(u, "wbin", "01"), u)
})

test_that("run_pipeline drops waterbody inlets and nothing else by default", {
  expect_equal(eval(formals(hydrofabric::run_pipeline)$exclude_hl_types), "wbin")
})

test_that("inlets go; outlets, gages and coastal outlets stay", {
  # An outlet is where a lake discharges, so it is a real routing boundary and
  # earns its levelpath break. An inlet names no boundary a consumer asks for.
  l <- list(hydrolocations = rbind(hl(1, "wbin"), hl(2, "wbout"), hl(3, "nwis"),
                                   hl(4, "rfc"), hl(5, "owp_atl")),
            pois = NULL)
  r <- hydrofabric:::.drop_hl_types(l, "wbin", "01")
  expect_setequal(r$hydrolocations$hl_type, c("wbout", "nwis", "rfc", "owp_atl"))
})
