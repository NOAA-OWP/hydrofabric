test_that("lineage_build normalizes split-suffixed and unsplit refactored IDs", {
  refactored <- data.frame(
    flowpath_id = c("100", "200.1", "200.2", "300"),
    stringsAsFactors = FALSE
  )
  refactored <- sf::st_as_sf(refactored, geometry = sf::st_sfc(
    sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)),
    sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)),
    sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)),
    sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE))))

  reconciled <- data.frame(
    ID = c(1L, 2L, 3L),
    member_flowpath_id = c("100", "200.1,200.2,200", "300"),
    stringsAsFactors = FALSE)
  reconciled <- sf::st_as_sf(reconciled, geometry = sf::st_sfc(
    sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)),
    sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)),
    sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE))))

  lg <- lineage_build(refactored, reconciled)

  expect_s3_class(lg, "data.table")
  expect_equal(nrow(lg), 4L)
  expect_setequal(lg$refactored_id, c("100", "200.1", "200.2", "300"))
  expect_setequal(lg$nhd_comid, c("100", "200", "300"))
  expect_equal(lg[refactored_id == "200.1", split_part], 1L)
  expect_equal(lg[refactored_id == "100",   split_part], 0L)
  # split children of 200 share a reconciled group
  expect_equal(length(unique(lg[nhd_comid == "200", reconciled_id])), 1L)
})

test_that("lineage_attach_ngen applies the fp-/cat- prefix", {
  lineage <- data.table::data.table(
    nhd_comid = "100", split_part = 0L, refactored_id = "100",
    reconciled_id = 1L, aggregated_id = "1",
    final_divide_id = NA_character_, final_flowpath_id = NA_character_)
  lg <- lineage_attach_ngen(lineage)
  expect_equal(lg$final_flowpath_id, "fp-1")
  expect_equal(lg$final_divide_id,   "cat-1")
})

test_that("lineage_member_csv reconstructs member lists by group", {
  lineage <- data.table::data.table(
    nhd_comid = c("100", "200", "200", "300"),
    split_part = c(0L, 1L, 2L, 0L),
    refactored_id = c("100", "200.1", "200.2", "300"),
    reconciled_id = c(1L, 2L, 2L, 3L),
    aggregated_id = c("A", "A", "A", "B"),
    final_divide_id = NA_character_, final_flowpath_id = NA_character_)
  csv <- lineage_member_csv(lineage, "aggregated_id")
  a <- csv[aggregated_id == "A", member_comid]
  expect_equal(sort(strsplit(a, ",", fixed = TRUE)[[1]]),
    c("100", "200.1", "200.2"))
})

test_that("check_invariants refactored stage passes for clean inputs", {
  refactored <- sf::st_as_sf(
    data.frame(flowpath_id = c("A", "B", "C"),
      flowpath_toid = c("B", "C", "0")),
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)),
      sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)),
      sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE))))
  reconciled <- sf::st_as_sf(
    data.frame(ID = c(1L, 2L), member_flowpath_id = c("A,B", "C")),
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)),
      sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE))))
  out <- suppressMessages(
    check_invariants("refactored",
      refactored = refactored, reconciled = reconciled))
  expect_true(out$ok)
})

test_that("check_invariants refactored catches dangling toid", {
  refactored <- sf::st_as_sf(
    data.frame(flowpath_id = c("A", "B"),
      flowpath_toid = c("B", "GHOST")),
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)),
      sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE))))
  expect_error(
    suppressMessages(check_invariants("refactored", refactored = refactored)),
    "refactored_topology")
})

test_that("check_invariants ngen flags duplicate flowline_id in the flowlines layer", {
  net <- data.frame(flowpath_id = "fp-1", flowline_id = c(1L, 2L, 3L),
    flowline_toid = c(2L, 3L, 0L), reference_id = c(10L, 11L, 12L),
    vpuid = "01", stringsAsFactors = FALSE)
  # flowlines layer with a duplicated flowline_id (the regression signature)
  fl_dup <- data.frame(flowline_id = c(1L, 1L, 2L, 3L),
    reconciled_id = c(5L, 5L, 6L, 7L), stringsAsFactors = FALSE)
  res <- suppressMessages(check_invariants("ngen",
    network = net, flowlines = fl_dup, strict = FALSE))
  expect_false(res$checks$flowlines_id_unique$ok)

  # unique per-segment flowline_id passes
  fl_ok <- data.frame(flowline_id = c(1L, 2L, 3L),
    reconciled_id = c(5L, 5L, 6L), stringsAsFactors = FALSE)
  res2 <- suppressMessages(check_invariants("ngen",
    network = net, flowlines = fl_ok, strict = FALSE))
  expect_true(res2$checks$flowlines_id_unique$ok)
})

test_that("check_invariants ngen flags implausible reference-COMID attribution on short flowlines", {
  fl <- data.frame(flowline_id = c(1L, 2L, 3L), reconciled_id = c(5L, 6L, 7L),
    lengthkm = c(0.5, 0.5, 0.5), stringsAsFactors = FALSE)
  # levelpath-inflation bug: two 0.5 km siblings each INHERIT the same 15-COMID
  # levelpath set (members shared across flowlines, not owned), the sharing is
  # the signature the invariant flags, not the raw count.
  net_bad <- data.frame(flowpath_id = "fp-1",
    flowline_id = c(rep(1L, 15), rep(2L, 15), 3L), flowline_toid = 0L,
    reference_id = c(101:115, 101:115, 201L), vpuid = "01", stringsAsFactors = FALSE)
  res <- suppressMessages(check_invariants("ngen",
    network = net_bad, flowlines = fl, strict = FALSE))
  expect_false(res$checks$flowline_reference_plausible$ok)

  # plausible A: each short flowline maps to only a few reference COMIDs
  net_ok <- data.frame(flowpath_id = "fp-1", flowline_id = c(1L, 1L, 2L),
    flowline_toid = 0L, reference_id = c(101L, 102L, 201L),
    vpuid = "01", stringsAsFactors = FALSE)
  res2 <- suppressMessages(check_invariants("ngen",
    network = net_ok, flowlines = fl, strict = FALSE))
  expect_true(res2$checks$flowline_reference_plausible$ok)

  # plausible B: a 0.5 km flowline with 15 *private* COMIDs (dense reference
  # tiling, each member owned by exactly one flowline) must NOT be flagged,
  # raw count alone is legitimate; only shared inheritance is the bug.
  net_private <- data.frame(flowpath_id = "fp-1",
    flowline_id = c(rep(1L, 15), 2L), flowline_toid = 0L,
    reference_id = c(101:115, 201L), vpuid = "01", stringsAsFactors = FALSE)
  res3 <- suppressMessages(check_invariants("ngen",
    network = net_private, flowlines = fl, strict = FALSE))
  expect_true(res3$checks$flowline_reference_plausible$ok)
})
