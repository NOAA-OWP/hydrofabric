.mk_hl <- function(layer = "hydrolocations", n = 3, vpuid = "01",
                   flowpath_id = c("1", "2", "3"), dup = FALSE) {
  f <- tempfile(fileext = ".gpkg")
  d <- sf::st_as_sf(data.frame(
    vpuid = vpuid, flowpath_id = flowpath_id,
    hl_reference = paste0("nwis-", seq_len(n)),
    hl_type = "nwis", hl_class = "gage", poi_id = NA_character_,
    x = seq_len(n), y = seq_len(n)), coords = c("x", "y"), crs = 5070)
  if (dup) d <- rbind(d, d[1, ])
  sf::st_write(d, f, layer, quiet = TRUE)
  f
}

test_that("a layer the build cannot read is rejected before anything else", {
  f <- .mk_hl(layer = "owp_hl")
  r <- suppressMessages(check_hydrolocation_input(f))
  expect_false(r$ok)
  expect_false(r$checks$hl_layer_readable$ok)
  # It stops there rather than reporting on rows the build will never see.
  expect_null(r$checks$hl_geometry_is_point)
  # Naming the layer explicitly gets past it.
  expect_true(suppressMessages(
    check_hydrolocation_input(f, layer = "owp_hl"))$checks$hl_layer_readable$ok)
})

test_that("a clean set passes, and strict raises on a failure", {
  f <- .mk_hl()
  r <- suppressMessages(check_hydrolocation_input(f))
  expect_true(r$ok)
  expect_error(suppressMessages(
    check_hydrolocation_input(.mk_hl(layer = "owp_hl"), strict = TRUE)),
    "NOT FOUND")
})

test_that("duplicate rows, unknown vpuid and mis-indexed ids are caught", {
  ref <- list(conus = list(ids = c("1", "2", "3"), vpuids = "01"))

  d <- suppressMessages(check_hydrolocation_input(.mk_hl(dup = TRUE), references = ref))
  expect_false(d$checks$hl_no_duplicate_rows$ok)

  v <- suppressMessages(check_hydrolocation_input(.mk_hl(vpuid = "03"), references = ref))
  expect_false(v$checks$hl_vpuid_known$ok)

  # Right vpuid, identifiers from another network: resolves 0%.
  m <- suppressMessages(check_hydrolocation_input(
    .mk_hl(flowpath_id = c("900", "901", "902")), references = ref))
  expect_false(m$checks$hl_indexed_conus$ok)
  expect_match(m$checks$hl_indexed_conus$msg, "0 of 3")

  ok <- suppressMessages(check_hydrolocation_input(.mk_hl(), references = ref))
  expect_true(ok$checks$hl_indexed_conus$ok)
})

test_that("hl_reference_index reads ids and vpuids without geometry", {
  f <- .mk_hl()
  # Reuse the hl fixture as a stand-in reference: same two columns.
  ix <- hl_reference_index(f, layer = "hydrolocations")
  expect_equal(ix$ids, c("1", "2", "3"))
  expect_equal(ix$vpuids, "01")
})
