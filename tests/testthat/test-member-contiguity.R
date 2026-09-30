# Flowpath member contiguity and cross-grain topology agreement
#
# An aggregated flowpath contracts a set of reconciled reaches into one
# feature. Two properties have to hold and neither was asserted anywhere:
#
#   1. the members form ONE connected chain, a flowpath with two members that
#      both leave it has two outlets;
#   2. the flowpath-grain topology (flowpath_toid through the nexus layer) names
#      the same downstream flowpath the flowline-grain topology reaches.
#
# The broken fixture is the real VPU 18 interleave, reduced: fp-A owns reaches
# on both sides of fp-B, so walking flowline_toid gives fp-A -> fp-B -> fp-A
# while fp-A declares fp-C. The declared flowpath graph is still a DAG in that
# state, which is why every pre-existing gate passed.
#
#   flowline  1(fp-A,spine) -> 2(fp-B) -> 3(fp-A,off-spine) -> 4(fp-C) -> 0
#   declared  fp-A -> nex-C -> fp-C          (but the spine reaches fp-B)

.mc_broken <- function() {
  list(
    flowlines = data.frame(
      flowline_id          = 1:4,
      flowline_toid        = c(2L, 3L, 4L, 0L),
      flowpath_id          = c("fp-A", "fp-B", "fp-A", "fp-C"),
      flowpath_constituent = c(TRUE, TRUE, FALSE, TRUE),
      lengthkm             = c(1.0, 0.5, 0.7, 2.0),
      stringsAsFactors     = FALSE),
    flowpaths = data.frame(
      flowpath_id   = c("fp-A", "fp-B", "fp-C"),
      flowpath_toid = c("nex-C", "nex-A", "tnx-C"),
      lengthkm      = c(1.0, 0.5, 2.0),
      # carried because other checks in .hf_check_ngen read it
      has_divide    = c(TRUE, TRUE, TRUE),
      stringsAsFactors = FALSE),
    # Terminal nexuses are real rows in the product: "tnx-N" with an empty
    # nexus_toid. Omitting them leaves flowpath_toid unresolvable and other
    # checks in the same function build a frame with NA row names.
    nexus = data.frame(
      nexus_id   = c("nex-A", "nex-B", "nex-C", "tnx-C"),
      nexus_toid = c("fp-A", "fp-B", "fp-C", ""),
      stringsAsFactors = FALSE))
}

# fp-A owns 1 and 2 as one chain; its single exit reaches fp-B, which is what it
# declares.
.mc_clean <- function() {
  list(
    flowlines = data.frame(
      flowline_id          = 1:3,
      flowline_toid        = c(2L, 3L, 0L),
      flowpath_id          = c("fp-A", "fp-A", "fp-B"),
      flowpath_constituent = c(TRUE, TRUE, TRUE),
      lengthkm             = c(1.0, 0.5, 2.0),
      stringsAsFactors     = FALSE),
    flowpaths = data.frame(
      flowpath_id   = c("fp-A", "fp-B"),
      flowpath_toid = c("nex-B", "tnx-B"),
      lengthkm      = c(1.5, 2.0),
      has_divide    = c(TRUE, TRUE),
      stringsAsFactors = FALSE),
    nexus = data.frame(
      nexus_id   = c("nex-A", "nex-B", "tnx-B"),
      nexus_toid = c("fp-A", "fp-B", ""),
      stringsAsFactors = FALSE))
}

.mc <- function(args) hydrofabric:::.hf_check_ngen(args, FALSE, 0)

test_that("an interleaved flowpath is reported as discontiguous", {
  r <- .mc(.mc_broken())
  chk <- r$flowpath_members_contiguous
  expect_false(is.null(chk))
  expect_match(chk$msg, "^1 of 3 flowpath")
  expect_match(chk$msg, "fp-A")
  # INFO while the population is non-zero, per the suite's one-directional
  # promotion policy, it must not fail a build that already ships this.
  expect_identical(chk$kind, "info")
})

test_that("a contiguous fabric passes contiguity as a gate", {
  r <- .mc(.mc_clean())
  chk <- r$flowpath_members_contiguous
  expect_true(chk$ok)
  expect_identical(chk$kind, "check")
  expect_match(chk$msg, "^0 of 2 flowpath")
})

test_that("a declared downstream the flowlines do not reach is reported", {
  r <- .mc(.mc_broken())
  chk <- r$flowpath_topology_matches_flowlines
  expect_false(is.null(chk))
  # fp-A declares fp-C; its routing spine reaches fp-B. fp-C is terminal (tnx).
  expect_match(chk$msg, "^1 of 2 non-terminal flowpath")
  expect_match(chk$msg, "fp-A")
  expect_identical(chk$kind, "info")
})

test_that("agreeing grains pass topology as a gate", {
  r <- .mc(.mc_clean())
  chk <- r$flowpath_topology_matches_flowlines
  expect_true(chk$ok)
  expect_identical(chk$kind, "check")
})

test_that("length is compared against constituent members only", {
  # The off-spine member (0.7 km) must NOT count: a tributary inside the
  # catchment contributes area but no length, so fp-A's 1.0 km is correct even
  # though its members total 1.7 km. Getting this wrong flags most of the fabric.
  r <- .mc(.mc_broken())
  chk <- r$flowpath_length_matches_members
  expect_true(chk$ok)
  expect_identical(chk$kind, "check")
  expect_match(chk$msg, "^0 of 3 flowpath")
})

test_that("a genuine length disagreement is reported", {
  a <- .mc_broken()
  a$flowpaths$lengthkm[a$flowpaths$flowpath_id == "fp-A"] <- 5.0
  chk <- .mc(a)$flowpath_length_matches_members
  # Non-zero reports as INFO, and .hf_info() fixes ok = TRUE on purpose:
  # .hf_report_checks folds every entry's `ok` into the stage verdict, so an
  # info carrying FALSE would fail the build it is meant not to block. The
  # count in the message is the signal; the kind is the contract.
  expect_identical(chk$kind, "info")
  expect_true(chk$ok)
  expect_match(chk$msg, "^1 of 3 flowpath")
})

test_that("the checks are skipped, not failed, when flowlines are absent", {
  a <- .mc_clean(); a$flowlines <- NULL
  r <- .mc(a)
  expect_null(r$flowpath_members_contiguous)
  expect_null(r$flowpath_topology_matches_flowlines)
  expect_null(r$flowpath_length_matches_members)
})

test_that("a terminal member is a chain end, not a second outlet", {
  # fp-A owns 1 and 2 as one chain (1 -> 2 -> fp-B) and also owns 3, a dangling
  # reference reach with no downstream at all. That is one way out, not two.
  # Counting "no downstream" as "leaves the flowpath" reported a second outlet
  # for every flowpath containing a terminal or dangling member, on VPU 01 it
  # turned 0 real cases into 116, and it hid the real result behind noise.
  a <- .mc_clean()
  a$flowlines <- rbind(a$flowlines, data.frame(
    flowline_id = 4L, flowline_toid = 0L, flowpath_id = "fp-A",
    flowpath_constituent = FALSE, lengthkm = 0.3, stringsAsFactors = FALSE))
  chk <- .mc(a)$flowpath_members_contiguous
  expect_true(chk$ok)
  expect_identical(chk$kind, "check")
  expect_match(chk$msg, "^0 of 2 flowpath")
})

test_that("an NA downstream is treated the same as a zero one", {
  a <- .mc_clean()
  a$flowlines <- rbind(a$flowlines, data.frame(
    flowline_id = 4L, flowline_toid = NA_integer_, flowpath_id = "fp-A",
    flowpath_constituent = FALSE, lengthkm = 0.3, stringsAsFactors = FALSE))
  expect_true(.mc(a)$flowpath_members_contiguous$ok)
})

test_that("two real exits are still reported when a terminal is also present", {
  # The terminal must not MASK a genuine second outlet either.
  a <- .mc_broken()
  a$flowlines <- rbind(a$flowlines, data.frame(
    flowline_id = 5L, flowline_toid = 0L, flowpath_id = "fp-A",
    flowpath_constituent = FALSE, lengthkm = 0.2, stringsAsFactors = FALSE))
  chk <- .mc(a)$flowpath_members_contiguous
  expect_match(chk$msg, "^1 of 3 flowpath")
  expect_match(chk$msg, "fp-A")
})

test_that("contiguity is measured without flowpath_constituent", {
  a <- .mc_broken()
  a$flowlines$flowpath_constituent <- NULL
  chk <- .mc(a)$flowpath_members_contiguous
  expect_match(chk$msg, "^1 of 3 flowpath")
})
