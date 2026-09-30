# split_siblings_share_catchment
#
# Splitting a reference reach at a POI splits its catchment: each fragment
# drains the land beside it. A 7.6 km river once shipped with a NULL divide; its member was 1764202.1, and sibling 1764202.2 held the whole
# 24.05 km2 parent catchment. This states the property that was violated.

.ssc_args <- function(fp) list(
  flowpaths = fp,
  nexus = data.frame(nexus_id = "tnx-1", nexus_toid = "", stringsAsFactors = FALSE))

.ssc <- function(fp) hydrofabric:::.hf_check_ngen(.ssc_args(fp), FALSE, 0)

test_that("a starved split sibling is reported", {
  fp <- data.frame(
    flowpath_id   = c("fp-1", "fp-2"),
    flowpath_toid = c("tnx-1", "tnx-1"),
    member_comid  = c("1764202.2", "1764202.1"),
    areasqkm      = c(24.05243, 0),
    has_divide    = c(TRUE, FALSE),
    stringsAsFactors = FALSE)
  chk <- .ssc(fp)$split_siblings_share_catchment
  expect_false(is.null(chk))
  expect_match(chk$msg, "^1 of 1 split reference reach")
  expect_match(chk$msg, "fp-2")
  expect_identical(chk$kind, "info")
})

test_that("siblings that both carry catchment pass as a gate", {
  fp <- data.frame(
    flowpath_id   = c("fp-1", "fp-2"),
    flowpath_toid = c("tnx-1", "tnx-1"),
    member_comid  = c("1764202.2", "1764202.1"),
    areasqkm      = c(16.7, 7.35),
    has_divide    = c(TRUE, TRUE),
    stringsAsFactors = FALSE)
  chk <- .ssc(fp)$split_siblings_share_catchment
  expect_true(chk$ok)
  expect_identical(chk$kind, "check")
  expect_match(chk$msg, "^0 of 1 split reference reach")
})

test_that("fragments all in one flowpath are not judged", {
  # Not a split at product grain: nothing was divided between owners.
  fp <- data.frame(
    flowpath_id   = "fp-1",
    flowpath_toid = "tnx-1",
    member_comid  = "1764202.1,1764202.2",
    areasqkm      = 0,
    has_divide    = FALSE,
    stringsAsFactors = FALSE)
  chk <- .ssc(fp)$split_siblings_share_catchment
  expect_match(chk$msg, "^0 of 0 split reference reach")
  expect_true(chk$ok)
})

test_that("unsplit reaches are outside its scope", {
  # A divide-less unsplit reach is a different mechanism (its divide is cleared
  # as disjoint and absorbed), so this check must not claim it.
  fp <- data.frame(
    flowpath_id   = c("fp-1", "fp-2"),
    flowpath_toid = c("tnx-1", "tnx-1"),
    member_comid  = c("41011712", "41151412"),
    areasqkm      = c(1.674, 0),
    has_divide    = c(TRUE, FALSE),
    stringsAsFactors = FALSE)
  chk <- .ssc(fp)$split_siblings_share_catchment
  expect_match(chk$msg, "^0 of 0 split reference reach")
})

test_that("all siblings starved is not reported, only a mixed split", {
  # If every fragment lacks catchment the reference reach had none to give;
  # that is a divide-less reach carried faithfully, not a failed split.
  fp <- data.frame(
    flowpath_id   = c("fp-1", "fp-2"),
    flowpath_toid = c("tnx-1", "tnx-1"),
    member_comid  = c("41151412.1", "41151412.2"),
    areasqkm      = c(0, 0),
    has_divide    = c(FALSE, FALSE),
    stringsAsFactors = FALSE)
  chk <- .ssc(fp)$split_siblings_share_catchment
  expect_match(chk$msg, "^0 of 1 split reference reach")
  expect_true(chk$ok)
})

test_that("the check is skipped when member_comid is absent", {
  fp <- data.frame(flowpath_id = "fp-1", flowpath_toid = "tnx-1",
                   areasqkm = 1, has_divide = TRUE, stringsAsFactors = FALSE)
  expect_null(.ssc(fp)$split_siblings_share_catchment)
})

test_that("a starved pair is found when it is not the first group", {
  # The single-group fixtures above cannot catch an indexing slip: with one
  # group the group's indices span the whole vector, so a subset and the whole
  # read alike. Three groups with the starved pair LAST is what distinguishes
  # them, an earlier version scored every group from the first group's areas
  # and reported 0 of 707 on a VPU whose true count was 29.
  fp <- data.frame(
    flowpath_id   = c("fpA1","fpA2","fpB1","fpB2","fpC1","fpC2"),
    flowpath_toid = rep("tnx-1", 6),
    member_comid  = c("100.1","100.2","200.1","200.2","300.1","300.2"),
    areasqkm      = c(10, 5,  7, 3,  24.05, 0),
    has_divide    = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE)
  chk <- .ssc(fp)$split_siblings_share_catchment
  expect_match(chk$msg, "^1 of 3 split reference reach")
  expect_match(chk$msg, "fpC2")
  expect_identical(chk$kind, "info")
})

test_that("every group is scored on its own areas", {
  # Same shape, starved pair FIRST. Both orderings must give the same count.
  fp <- data.frame(
    flowpath_id   = c("fpA1","fpA2","fpB1","fpB2","fpC1","fpC2"),
    flowpath_toid = rep("tnx-1", 6),
    member_comid  = c("100.1","100.2","200.1","200.2","300.1","300.2"),
    areasqkm      = c(24.05, 0,  7, 3,  10, 5),
    has_divide    = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE)
  chk <- .ssc(fp)$split_siblings_share_catchment
  expect_match(chk$msg, "^1 of 3 split reference reach")
  expect_match(chk$msg, "fpA2")
})

test_that("several starved groups are all counted", {
  fp <- data.frame(
    flowpath_id   = c("fpA1","fpA2","fpB1","fpB2","fpC1","fpC2"),
    flowpath_toid = rep("tnx-1", 6),
    member_comid  = c("100.1","100.2","200.1","200.2","300.1","300.2"),
    areasqkm      = c(10, 0,  7, 3,  24.05, 0),
    has_divide    = c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE)
  chk <- .ssc(fp)$split_siblings_share_catchment
  expect_match(chk$msg, "^2 of 3 split reference reach")
})
