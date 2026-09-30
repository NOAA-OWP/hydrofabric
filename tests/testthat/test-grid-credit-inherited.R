test_that("grid credits are carried forward from the merged VPU products", {
  mk <- function(attribution) {
    f <- tempfile(fileext = ".gpkg")
    con <- DBI::dbConnect(RSQLite::SQLite(), f)
    DBI::dbExecute(con, "CREATE TABLE gpkg_metadata (id INTEGER, md_standard_uri TEXT,
                         mime_type TEXT, metadata TEXT)")
    DBI::dbExecute(con, sprintf(
      "INSERT INTO gpkg_metadata VALUES (1, 'x', 'application/json', '%s')",
      jsonlite::toJSON(list(attribution = attribution), auto_unbox = TRUE)))
    DBI::dbDisconnect(con)
    f
  }
  hs <- hydrofabric:::.HF_ATTR_HYDROSHEDS
  gd <- hydrofabric:::.HF_ATTR_GEDTM

  expect_equal(hydrofabric:::.hf_grid_attribution_inherited(mk(paste("base.", hs))), hs)

  # A domain merges many VPUs; each distinct credit is carried once.
  two <- c(mk(paste("base.", hs)), mk(paste("base.", hs)), mk(paste("base.", gd)))
  expect_setequal(hydrofabric:::.hf_grid_attribution_inherited(two), c(hs, gd))

  # A product with no grid credit contributes none, and neither does a missing file.
  expect_length(hydrofabric:::.hf_grid_attribution_inherited(mk("base only.")), 0L)
  expect_length(hydrofabric:::.hf_grid_attribution_inherited("/no/such/file.gpkg"), 0L)
})
