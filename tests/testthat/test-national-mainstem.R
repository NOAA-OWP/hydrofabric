# Aggregation keys level paths by a number local to the VPU, with the reference
# level-path id carried alongside as orig_levelpathID. Finalize makes that id
# the mainstem_id and drops the helper column. A record that does not map
# one-to-one, or is missing, leaves the local ids.

.ms_fixture <- function(orig) {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  pt <- function(n) sf::st_sfc(lapply(seq_len(n), function(i) sf::st_point(c(i, i))), crs = 5070)
  g <- file.path(d, "ngen.gpkg")
  sf::st_write(sf::st_sf(flowpath_id = c("fp-1", "fp-2"), mainstem_id = c(1L, 2L),
    geometry = pt(2)), g, "flowpaths", quiet = TRUE)
  fl <- sf::st_sf(flowline_id = 1:3, mainstem_id = c(1L, 1L, 2L), geometry = pt(3))
  if (!is.null(orig)) fl$orig_levelpathID <- orig
  sf::st_write(fl, g, "flowlines", quiet = TRUE)
  con <- DBI::dbConnect(RSQLite::SQLite(), g)
  DBI::dbWriteTable(con, "network", data.frame(flowline_id = 1:3, mainstem_id = c(1L, 1L, 2L)))
  DBI::dbDisconnect(con)
  hydrofabric::enforce_integer_ids(g)
  g
}

.ms_col <- function(g, layer, col = "mainstem_id") {
  con <- DBI::dbConnect(RSQLite::SQLite(), g); on.exit(DBI::dbDisconnect(con))
  if (!col %in% DBI::dbListFields(con, layer)) return(NULL)
  DBI::dbGetQuery(con, sprintf('SELECT "%s" FROM %s ORDER BY rowid', col, layer))[[1]]
}

test_that("mainstem_id becomes the recorded reference level-path id", {
  g <- .ms_fixture(c(5e9, 5e9, 3e10))
  hydrofabric:::.national_mainstem_ids(g)
  expect_equal(.ms_col(g, "flowlines"), c(5e9, 5e9, 3e10))
  expect_equal(.ms_col(g, "flowpaths"), c(5e9, 3e10))
  expect_equal(.ms_col(g, "network"),   c(5e9, 5e9, 3e10))
  expect_null(.ms_col(g, "flowlines", "orig_levelpathID"))
})

test_that("a record that does not map one-to-one leaves the local ids", {
  g <- .ms_fixture(c(5e9, 7e9, 3e10))
  hydrofabric:::.national_mainstem_ids(g)
  expect_equal(.ms_col(g, "flowlines"), c(1, 1, 2))
})

test_that("fabrics without the record keep their local ids", {
  g <- .ms_fixture(NULL)
  hydrofabric:::.national_mainstem_ids(g)
  expect_equal(.ms_col(g, "flowpaths"), c(1, 2))
})
