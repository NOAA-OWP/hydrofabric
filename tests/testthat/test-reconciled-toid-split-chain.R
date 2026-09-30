# .fix_reconciled_toid_from_reference must honour the split chain
#
# split_flowlines cuts a reference reach into ".1" (upstream) ... ".N"
# (downstream) and chains .k -> .(k+1), only .N inheriting the parent's toid.
# Every fragment inherits the parent's identical hydroseq, so a pass that
# resolves a fragment through its BASE COMID sees one reach where there are
# several: each fragment looks up the parent's single reference toCOMID and the
# upstream fragment is pointed PAST its own sibling, straddling whichever group
# holds the reference downstream.
#
# The fixture is the real VPU 18 case, reduced. Reference reach 10 (split into
# 10.1 / 10.2 at a POI) drains to reach 20, which drains to 30:
#
#   reference   10 -> 20 -> 30
#   post-split  10.1 -> 10.2 -> 20
#   reconciled  g4{10.1}  g3{10.2}  g2{20}  g1{30}
#
# so g4's downstream is g3, NOT g2. Before the fix both g3 and g4 resolved to g2.

.rt_ref <- function() {
  list(
    # base COMID -> reference toCOMID
    toc = c("10" = "20", "20" = "30", "30" = "0"),
    # base COMID -> reference hydroseq (downstream = smaller)
    hs  = c("10" = 100, "20" = 90, "30" = 80))
}

.rt_rec <- function(members, toID) {
  x <- data.frame(ID = seq_along(members), toID = toID)
  x$member_flowpath_id <- members
  x
}

test_that("the upstream fragment of a split reach chains to its sibling", {
  r <- .rt_ref()
  rec <- .rt_rec(
    members = list("30", "20", "10.2", "10.1"),
    #                g1    g2      g3      g4
    toID    = c(0L, 1L, 2L, 2L))   # g4 -> g2 is the defect
  out <- hydrofabric:::.fix_reconciled_toid_from_reference(rec, r$toc, r$hs)

  expect_equal(out$toID[out$ID == 4L], 3)   # g4 -> g3, its own sibling
  expect_equal(out$toID[out$ID == 3L], 2)   # g3 -> g2, unchanged (already right)
  expect_equal(out$toID[out$ID == 2L], 1)
})

test_that("an unsplit reach is left on the reference toCOMID path", {
  r <- .rt_ref()
  rec <- .rt_rec(members = list("30", "20", "10"), toID = c(0L, 1L, 2L))
  out <- hydrofabric:::.fix_reconciled_toid_from_reference(rec, r$toc, r$hs)
  # 10 -> 20 is group 2; nothing about the reference chain changes
  expect_equal(out$toID[out$ID == 3L], 2)
  expect_equal(out$toID[out$ID == 2L], 1)
})

test_that("fragments that collapsed into ONE group stay internal", {
  # Both fragments in a single group: the sibling is intra-group, so the group
  # takes the reference toCOMID of its most-downstream member.
  r <- .rt_ref()
  rec <- .rt_rec(members = list("30", "20", c("10.1", "10.2")),
                 toID = c(0L, 1L, 2L))
  out <- hydrofabric:::.fix_reconciled_toid_from_reference(rec, r$toc, r$hs)
  expect_equal(out$toID[out$ID == 3L], 2)
})

test_that("a three-way split chains .1 -> .2 -> .3", {
  r <- .rt_ref()
  rec <- .rt_rec(
    members = list("30", "20", "10.3", "10.2", "10.1"),
    #                g1    g2      g3      g4      g5
    toID    = c(0L, 1L, 2L, 2L, 2L))
  out <- hydrofabric:::.fix_reconciled_toid_from_reference(rec, r$toc, r$hs)
  expect_equal(out$toID[out$ID == 5L], 4)   # .1 -> .2
  expect_equal(out$toID[out$ID == 4L], 3)   # .2 -> .3
  expect_equal(out$toID[out$ID == 3L], 2)   # .3 -> reference downstream
})

test_that("the outlet member is the largest part, not an arbitrary hydroseq tie", {
  # A group holding .1 and .2 plus a separate group holding .3. Both fragments
  # in the first group carry the parent's identical hydroseq, so without the
  # part tie-break the outlet could be .1 and the group would point at .2,
  # which is inside itself.
  r <- .rt_ref()
  rec <- .rt_rec(members = list("30", "20", "10.3", c("10.1", "10.2")),
                 toID = c(0L, 1L, 2L, 2L))
  out <- hydrofabric:::.fix_reconciled_toid_from_reference(rec, r$toc, r$hs)
  expect_equal(out$toID[out$ID == 4L], 3)   # outlet is .2, sibling .3 is group 3
})

test_that("a downstream reach that is itself split is entered at .1", {
  # Reach 10 (unsplit) drains to reach 20, which was split into 20.1 / 20.2
  # landing in different groups. Flow enters 20 at its UPSTREAM fragment, so 10
  # must resolve to 20.1's group. A base-COMID map cannot express that, it has
  # one slot for "20" and lands on whichever fragment it saw last. Groups are
  # ordered here so that "last seen" is the WRONG fragment (20.2).
  r <- .rt_ref()
  rec <- .rt_rec(
    members = list("30", "20.1", "20.2", "10"),
    #                g1      g2      g3     g4
    toID    = c(0L, 3L, 1L, 3L))
  out <- hydrofabric:::.fix_reconciled_toid_from_reference(rec, r$toc, r$hs)
  expect_equal(out$toID[out$ID == 4L], 2)   # 10 -> 20.1's group, not 20.2's
  expect_equal(out$toID[out$ID == 2L], 3)   # 20.1 -> 20.2
  expect_equal(out$toID[out$ID == 3L], 1)   # 20.2 -> 30
})

test_that("multi-digit split parts chain correctly (.9 -> .10, not .9 -> .1)", {
  # A reach split into ten or more fragments. String-sorting the suffix would put
  # ".10" before ".9"; the part must be compared as an integer, which is why the
  # package's own parse_member_part does the parsing.
  r <- .rt_ref()
  rec <- .rt_rec(
    members = list("30", "20", "10.10", "10.9"),
    #                g1    g2       g3      g4
    toID    = c(0L, 1L, 2L, 2L))
  out <- hydrofabric:::.fix_reconciled_toid_from_reference(rec, r$toc, r$hs)
  expect_equal(out$toID[out$ID == 4L], 3)   # .9 -> .10
  expect_equal(out$toID[out$ID == 3L], 2)   # .10 is the last part -> reference downstream
})

test_that("the entry fragment is the smallest SURVIVING part", {
  # ".1" did not survive the collapse; flow enters reach 20 at ".2".
  r <- .rt_ref()
  rec <- .rt_rec(members = list("30", "20.2", "20.3", "10"),
                 #               g1     g2      g3     g4
                 toID    = c(0L, 3L, 1L, 2L))
  out <- hydrofabric:::.fix_reconciled_toid_from_reference(rec, r$toc, r$hs)
  expect_equal(out$toID[out$ID == 4L], 2)   # 10 -> 20.2's group
})

test_that("a missing member_flowpath_id column is a no-op", {
  rec <- data.frame(ID = 1:2, toID = c(0L, 1L))
  r <- .rt_ref()
  expect_identical(
    hydrofabric:::.fix_reconciled_toid_from_reference(rec, r$toc, r$hs), rec)
})
