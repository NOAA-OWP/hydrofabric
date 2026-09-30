test_that("run_vpus merges shared and per-VPU arguments, and reports per VPU", {
  seen <- list()
  local_mocked_bindings(
    run_pipeline = function(...) {
      a <- list(...); seen[[a$vpuid]] <<- a
      if (a$vpuid == "02") stop("boom")
      list(vpuid = a$vpuid, status = "ok",
           ngen_gpkg = file.path(a$out_dir, a$vpuid, "ngen.gpkg"), skipped = FALSE)
    },
    .package = "hydrofabric")
  res <- run_vpus(
    inputs  = list("01" = list(), "02" = list(),
                   "20" = list(fdr_path = "hi/fdr.vrt")),
    out_dir = "out", fdr_path = "conus/fdr.vrt", fac_path = "conus/fac.vrt")

  expect_equal(res$status, c("ok", "error", "ok"))
  expect_equal(res$message[2], "boom")
  # Shared arguments reach every VPU; a per-VPU value overrides the shared one.
  expect_equal(seen[["01"]]$fdr_path, "conus/fdr.vrt")
  expect_equal(seen[["20"]]$fdr_path, "hi/fdr.vrt")
  expect_equal(seen[["20"]]$fac_path, "conus/fac.vrt")
})

test_that("a character vector of vpuids is accepted, and errors can stop the run", {
  local_mocked_bindings(
    run_pipeline = function(...) {
      a <- list(...)
      if (a$vpuid == "02") stop("boom")
      list(vpuid = a$vpuid, status = "ok", ngen_gpkg = "x.gpkg", skipped = FALSE)
    },
    .package = "hydrofabric")
  expect_error(run_vpus(c("01", "02"), out_dir = "out", continue_on_error = FALSE),
               "02 failed")
})

test_that("a partial set is not merged", {
  local_mocked_bindings(
    run_pipeline = function(...) {
      a <- list(...)
      if (a$vpuid == "02") stop("boom")
      list(vpuid = a$vpuid, status = "ok", ngen_gpkg = "x.gpkg", skipped = FALSE)
    },
    merge_vpus = function(...) stop("merge_vpus should not be reached"),
    .package = "hydrofabric")
  expect_error(run_vpus(c("01", "02"), out_dir = "out", merge = TRUE),
               "refusing to merge")
})
