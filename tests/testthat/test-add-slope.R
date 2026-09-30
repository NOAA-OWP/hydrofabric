# Flowpath slope is a length-weighted mean over member reaches. Its source can
# be a slope table (`id`, `slope_m_per_m`) or a channel table carrying
# `slope_m_m`; the same values must give the same slope whichever shape holds
# them, so a build needs no separate slope input when the channel table has it.

skip_if_not_installed("arrow")

.slope_fixture <- function(dir) {
  fp <- sf::st_sf(
    flowpath_id  = c("fp-1", "fp-2"),
    member_comid = c("101,102", "103"),
    vpuid        = "01",
    geometry     = sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(1, 1)),
                              crs = 5070))
  ngen <- file.path(dir, "ngen.gpkg")
  sf::st_write(fp, ngen, "flowpaths", quiet = TRUE)

  ref <- file.path(dir, "ref.gpkg")
  con <- DBI::dbConnect(RSQLite::SQLite(), ref)
  DBI::dbWriteTable(con, "flowpaths", data.frame(
    flowpath_id = c(101L, 102L, 103L), lengthkm = c(1, 3, 2), vpuid = "01"))
  DBI::dbDisconnect(con)
  list(ngen = ngen, ref = ref)
}

.slope_after <- function(src) {
  d <- withr::local_tempdir()
  f <- .slope_fixture(d)
  hydrofabric:::.add_slope(f$ngen, f$ref, "01", slope_parquet = src)
  sf::st_read(f$ngen, "flowpaths", quiet = TRUE)$slope
}

test_that("a channel table and a slope table give the same slope", {
  d <- withr::local_tempdir()
  slope_tbl <- file.path(d, "slope.parquet")
  arrow::write_parquet(data.frame(id = c("101", "102", "103"),
    slope_m_per_m = c(0.01, 0.03, 0.002)), slope_tbl)
  chan_tbl <- file.path(d, "channel.parquet")
  arrow::write_parquet(data.frame(flowpath_id = c(101, 102, 103),
    tw_bf_m = 5, y_bf_m = 1, r_bf = 2, slope_m_m = c(0.01, 0.03, 0.002)),
    chan_tbl)

  from_slope <- .slope_after(slope_tbl)
  expect_equal(from_slope, c((0.01 * 1 + 0.03 * 3) / 4, 0.002))
  expect_identical(.slope_after(chan_tbl), from_slope)
})

test_that("a table with no slope column is not used as the source", {
  d <- withr::local_tempdir()
  no_slope <- file.path(d, "channel.parquet")
  arrow::write_parquet(data.frame(flowpath_id = 101:103, tw_bf_m = 5), no_slope)
  expect_message(s <- .slope_after(no_slope), "carries no slope column")
  expect_true(all(is.na(s)))
})

test_that("slope table columns are recognised in either shape", {
  cols <- hydrofabric:::.slope_table_cols
  expect_identical(cols(c("id", "slope_m_per_m")),
                   c(id = "id", slope = "slope_m_per_m"))
  expect_identical(cols(c("flowpath_id", "tw_bf_m", "slope_m_m")),
                   c(id = "flowpath_id", slope = "slope_m_m"))
  expect_null(cols(c("flowpath_id", "tw_bf_m")))
})
