test_that("enforce_integer_ids retypes integer-valued ids, guards strings/measures, indexes", {
  skip_if_not_installed("sf")
  skip_if_not_installed("RSQLite")

  tmp <- tempfile(fileext = ".gpkg")
  on.exit(unlink(tmp), add = TRUE)

  pts <- sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(1, 1)),
                    sf::st_point(c(2, 2)), crs = 4326)
  x <- sf::st_sf(
    flowpath_id   = c(101, 102, 103),           # REAL, integer-valued -> INTEGER
    flowpath_toid = c(102, 103, 0),             # REAL, integer-valued -> INTEGER
    divide_id     = c("cat-101", "cat-102", NA),# prefixed string       -> TEXT
    streamorder   = c(1, 2, 2),                 # integer-valued        -> INTEGER
    lengthkm      = c(1.5, 2.25, 0.3),          # real measure          -> REAL
    geom = pts)
  sf::st_write(x, tmp, "flowpaths", quiet = TRUE)

  coltype <- function(col) {
    con <- DBI::dbConnect(RSQLite::SQLite(), tmp); on.exit(DBI::dbDisconnect(con))
    ti <- DBI::dbGetQuery(con, 'PRAGMA table_info("flowpaths")')
    ti$type[ti$name == col]
  }

  enforce_integer_ids(tmp)

  expect_equal(coltype("flowpath_id"),   "INTEGER")
  expect_equal(coltype("flowpath_toid"), "INTEGER")
  expect_equal(coltype("streamorder"),   "INTEGER")
  expect_equal(coltype("divide_id"),     "TEXT")   # prefixed string untouched
  expect_equal(coltype("lengthkm"),      "REAL")   # measure untouched

  # values + string preserved and GDAL-readable
  y <- sf::st_read(tmp, "flowpaths", quiet = TRUE)
  expect_equal(as.numeric(y$flowpath_id), c(101, 102, 103))
  expect_equal(y$divide_id, c("cat-101", "cat-102", NA))

  # join keys indexed
  con <- DBI::dbConnect(RSQLite::SQLite(), tmp)
  idx <- DBI::dbGetQuery(con,
    "SELECT name FROM sqlite_master WHERE type='index' AND name LIKE 'idx_%'")$name
  DBI::dbDisconnect(con)
  expect_true(all(c("idx_flowpaths_flowpath_id", "idx_flowpaths_flowpath_toid",
                    "idx_flowpaths_divide_id") %in% idx))

  # idempotent: second run is a no-op and does not error
  expect_no_error(enforce_integer_ids(tmp))
  expect_equal(coltype("flowpath_id"), "INTEGER")
})
