# SQLite refuses to DROP a column an index references. Because the retype does
# all its columns inside one transaction, a single indexed column aborted the
# whole layer, which is how idx_flowpaths_upstream_id left hydroseq,
# mainstem_id, streamorder, has_divide, upstream_id and num_upstreams as REAL in
# every merged domain while the merge still reported success.

make_db <- function(indexed = TRUE) {
  f <- tempfile(fileext = ".gpkg")
  con <- DBI::dbConnect(RSQLite::SQLite(), f)
  DBI::dbExecute(con, "CREATE TABLE flowpaths (
     flowpath_id TEXT, hydroseq REAL, mainstem_id REAL, upstream_id REAL)")
  DBI::dbExecute(con, "INSERT INTO flowpaths VALUES ('fp-1', 1.0, 10.0, 1.0), ('fp-2', 2.0, 20.0, 2.0)")
  DBI::dbExecute(con, "CREATE INDEX idx_flowpaths_flowpath_id ON flowpaths(flowpath_id)")
  if (indexed) DBI::dbExecute(con, "CREATE INDEX idx_flowpaths_upstream_id ON flowpaths(upstream_id)")
  DBI::dbDisconnect(con)
  f
}
types <- function(f, cols) {
  con <- DBI::dbConnect(RSQLite::SQLite(), f); on.exit(DBI::dbDisconnect(con))
  ti <- DBI::dbGetQuery(con, "PRAGMA table_info(flowpaths)")
  toupper(stats::setNames(ti$type, ti$name)[cols])
}
idx_names <- function(f) {
  con <- DBI::dbConnect(RSQLite::SQLite(), f); on.exit(DBI::dbDisconnect(con))
  sort(DBI::dbGetQuery(con,
    "SELECT name FROM sqlite_master WHERE type='index' AND tbl_name='flowpaths' AND sql IS NOT NULL")$name)
}

test_that("an indexed column no longer blocks the retype", {
  f <- make_db(indexed = TRUE)
  con <- DBI::dbConnect(RSQLite::SQLite(), f)
  hydrofabric:::.retype_cols_integer(con, "flowpaths",
    c("hydroseq", "mainstem_id", "upstream_id"))
  DBI::dbDisconnect(con)
  expect_identical(unname(types(f, c("hydroseq", "mainstem_id", "upstream_id"))),
                   rep("INTEGER", 3))
})

test_that("the indexes are restored afterwards", {
  f <- make_db(indexed = TRUE)
  before <- idx_names(f)
  con <- DBI::dbConnect(RSQLite::SQLite(), f)
  hydrofabric:::.retype_cols_integer(con, "flowpaths", c("upstream_id"))
  DBI::dbDisconnect(con)
  expect_identical(idx_names(f), before)
})

test_that("values survive the retype unchanged", {
  f <- make_db(indexed = TRUE)
  con <- DBI::dbConnect(RSQLite::SQLite(), f)
  hydrofabric:::.retype_cols_integer(con, "flowpaths", c("hydroseq", "upstream_id"))
  d <- DBI::dbGetQuery(con, "SELECT flowpath_id, hydroseq, upstream_id FROM flowpaths ORDER BY flowpath_id")
  DBI::dbDisconnect(con)
  expect_identical(d$flowpath_id, c("fp-1", "fp-2"))
  expect_identical(d$hydroseq, c(1L, 2L))
  expect_identical(d$upstream_id, c(1L, 2L))
})

test_that("the un-indexed case still works (no regression)", {
  f <- make_db(indexed = FALSE)
  con <- DBI::dbConnect(RSQLite::SQLite(), f)
  hydrofabric:::.retype_cols_integer(con, "flowpaths", c("hydroseq", "mainstem_id"))
  DBI::dbDisconnect(con)
  expect_identical(unname(types(f, c("hydroseq", "mainstem_id"))), rep("INTEGER", 2))
})

test_that("enforce_integer_ids types the whole layer end to end", {
  f <- make_db(indexed = TRUE)
  # gpkg_contents is what the function reads to enumerate layers
  con <- DBI::dbConnect(RSQLite::SQLite(), f)
  DBI::dbExecute(con, "CREATE TABLE gpkg_contents (table_name TEXT)")
  DBI::dbExecute(con, "INSERT INTO gpkg_contents VALUES ('flowpaths')")
  DBI::dbDisconnect(con)
  hydrofabric::enforce_integer_ids(f)
  expect_identical(unname(types(f, c("hydroseq", "mainstem_id", "upstream_id"))),
                   rep("INTEGER", 3))
})
