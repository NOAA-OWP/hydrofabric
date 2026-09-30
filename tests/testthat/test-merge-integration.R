## Stage-4 (global merge) integration test.
##
## Runs merge_vpus() on two tiny committed per-VPU ngen fixtures (closed
## sub-basins from a real VPU, relabelled VPU 01 + 02; VPU 01 also contains one
## isolated single-reach catchment) and asserts the merged product is clean,
## mainstem-preserving, and fully hydroseq-ordered. This seam shipped two silent
## regressions with zero end-to-end coverage: mainstem_id remapped through the
## flowpath-id table (-> 78% NULL), and the hydroseq recompute omitting isolated
## nodes as graph vertices (-> 31-45% NA on island domains). The isolated reach
## in VPU 01 exercises that second path.

test_that("merge_vpus yields a clean, mainstem- and hydroseq-preserving merge", {
  skip_on_cran()
  fdir <- testthat::test_path("..", "fixtures", "merge_min")
  f1 <- file.path(fdir, "ngen_vpu01.gpkg")
  f2 <- file.path(fdir, "ngen_vpu02.gpkg")
  skip_if_not(file.exists(f1) && file.exists(f2), "merge_min fixtures missing")

  out <- withr::local_tempdir()
  res <- suppressWarnings(suppressMessages(
    merge_vpus(gpkg_map = list("01" = f1, "02" = f2),
               out_dir = out, overwrite = TRUE, strict = FALSE)))

  sc <- file.path(out, "superconus_ngen.gpkg")
  expect_true(file.exists(sc))

  # The merge's own contract passes, including the mainstem + hydroseq guards.
  expect_true(isTRUE(res$invariants$superconus$ok))
  expect_true(res$invariants$superconus$checks$mainstem_id_populated$ok)
  expect_true(res$invariants$superconus$checks$hydroseq_valid$ok)

  fp <- sf::read_sf(sc, "flowpaths", quiet = TRUE)
  expect_equal(nrow(fp), 21L)                               # 15 + 6 conserved
  expect_true(all(startsWith(as.character(fp$flowpath_id), "fp-")))
  expect_false(any(duplicated(fp$flowpath_id)))             # globally unique ids
  expect_equal(sum(is.na(fp$mainstem_id)), 0L)              # mainstem carried through
  expect_equal(sum(is.na(fp$hydroseq)), 0L)                 # incl. the isolated reach

  # Every non-terminal flowpath_toid resolves to a real nexus (no disconnect).
  nx   <- sf::read_sf(sc, "nexus", quiet = TRUE)
  toid <- as.character(fp$flowpath_toid)
  term <- grepl("^tnx-|^ter", toid)
  expect_true(all(toid[!term] %in% as.character(nx$nexus_id)))
})
