# Extracted from test-national-mainstem.R:30

# prequel ----------------------------------------------------------------------
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
  g
}
.ms_col <- function(g, layer, col = "mainstem_id") {
  con <- DBI::dbConnect(RSQLite::SQLite(), g); on.exit(DBI::dbDisconnect(con))
  if (!col %in% DBI::dbListFields(con, layer)) return(NULL)
  DBI::dbGetQuery(con, sprintf('SELECT "%s" FROM %s ORDER BY rowid', col, layer))[[1]]
}

# test -------------------------------------------------------------------------
g <- .ms_fixture(c(5e9, 5e9, 3e10))
hydrofabric:::.national_mainstem_ids(g)
expect_equal(.ms_col(g, "flowlines"), c(5e9, 5e9, 3e10))
