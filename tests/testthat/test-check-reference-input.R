# Reference-input validation. Fixtures are built in-line so the test needs no
# data: a healthy network, a 32-bit-wrapped toid column, and an unresolvable one.
suppressMessages({ library(sf); library(DBI); library(RSQLite) })

# Minimal reference: 4 reaches draining 1 -> 2 -> 3 -> 4(terminal).
.make_ref <- function(path, ids, toids, toid_type = "INTEGER") {
  g <- st_sfc(lapply(seq_along(ids), function(i)
    st_linestring(cbind(c(i, i + 1), c(0, 1)))), crs = 5070)
  fp <- st_sf(flowpath_id = ids, flowpath_toid = toids, lengthkm = 1,
              geom = g)
  if (file.exists(path)) unlink(path)
  st_write(fp, path, layer = "flowpaths", quiet = TRUE)
  # Force the declared column type when the test needs a narrow one.
  if (!identical(toid_type, "INTEGER")) {
    con <- dbConnect(SQLite(), path); on.exit(dbDisconnect(con))
    trg <- dbGetQuery(con, "SELECT name, sql FROM sqlite_master
                            WHERE type='trigger' AND tbl_name='flowpaths' AND sql IS NOT NULL")
    for (nm in trg$name) dbExecute(con, sprintf('DROP TRIGGER "%s"', nm))
    dbExecute(con, sprintf('ALTER TABLE flowpaths RENAME COLUMN flowpath_toid TO _t'))
    dbExecute(con, sprintf('ALTER TABLE flowpaths ADD COLUMN flowpath_toid %s', toid_type))
    dbExecute(con, 'UPDATE flowpaths SET flowpath_toid = _t')
    for (s in trg$sql) dbExecute(con, s)
  }
  path
}

test_that("a healthy reference passes", {
  p <- .make_ref(tempfile(fileext = ".gpkg"),
                 ids = c(101, 102, 103, 104), toids = c(102, 103, 104, 0))
  r <- suppressMessages(check_reference_input(p, strict = TRUE))
  expect_true(r$ok)
  expect_true(r$toid_resolves$ok)
  expect_true(r$toid_not_truncated$ok)
})

test_that("terminals are excluded from the resolution denominator", {
  # 3 of 4 reaches are terminal; the one non-terminal link resolves -> 100%.
  p <- .make_ref(tempfile(fileext = ".gpkg"),
                 ids = c(1, 2, 3, 4), toids = c(2, 0, 0, 0))
  r <- suppressMessages(check_reference_input(p, strict = TRUE))
  expect_true(r$ok)
})

test_that("a negative toid is reported as a 32-bit wrap and is fatal", {
  p <- .make_ref(tempfile(fileext = ".gpkg"),
                 ids = c(80000100000001, 80000100000002),
                 toids = c(-2055822587, 0))
  expect_error(suppressMessages(check_reference_input(p, strict = TRUE)),
               "toid_not_truncated")
  r <- suppressWarnings(suppressMessages(
    check_reference_input(p, strict = FALSE)))
  expect_false(r$ok)
  expect_false(r$toid_not_truncated$ok)
})

test_that("unresolvable topology is fatal even when values look plausible", {
  # Positive, in-range, but naming ids that do not exist.
  p <- .make_ref(tempfile(fileext = ".gpkg"),
                 ids = c(1, 2, 3, 4), toids = c(9001, 9002, 9003, 0))
  expect_error(suppressMessages(check_reference_input(p, strict = TRUE)),
               "toid_resolves")
})

test_that("a narrow column holding large ids warns but does not block", {
  # Intact data (resolves, no negatives) in a 32-bit column near the ceiling:
  # a latent hazard, not present damage, so the build must still proceed.
  p <- .make_ref(tempfile(fileext = ".gpkg"),
                 ids = c(2000565230, 2000565231), toids = c(2000565231, 0),
                 toid_type = "MEDIUMINT")
  expect_warning(r <- suppressMessages(check_reference_input(p, strict = TRUE)),
                 "at risk")
  expect_true(r$ok)                      # not blocked
  expect_false(r$id_width_safe$ok)       # but flagged
  expect_identical(r$id_width_safe$severity, "risk")
})

test_that("a missing layer or column is an explicit error", {
  p <- .make_ref(tempfile(fileext = ".gpkg"), ids = c(1, 2), toids = c(2, 0))
  expect_error(check_reference_input(p, layer = "nope"), "not in")
  expect_error(check_reference_input(p, toid = "nope"), "not in layer")
})
