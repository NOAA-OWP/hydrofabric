# .build_hw_collapse_fp candidate selection in ideal mode.
#
# The collapse exists to merge across levelpaths, which .cs_group cannot do: its
# repair helpers return early on a single-group chunk, so a whole chunk under the
# floor survives. It has to serve two jobs at once, and the OG aggregation keeps
# them separate within a levelpath:
#
#   .assign_id                      never build a group past the target
#   .pinch_sides/.middle_massage    repair what falls short, overshooting only
#                                   as much as the floor demands
#
# Selecting on the donor alone collapsed both into one rule and let two 9.9 km2
# units merge into 19.8.

.ci_net <- function(ids, toids, members, poi = NA_character_, lp = "1") {
  list(flowpaths = sf::st_sf(
    flowpath_id   = as.character(ids),
    flowpath_toid = as.character(toids),
    flowline_id   = as.character(members),
    member_comid  = as.character(members),
    levelpathid   = rep_len(as.character(lp), length(ids)),
    poi_id        = rep_len(poi, length(ids)),
    areasqkm      = NA_real_,
    lengthkm      = 1,
    geometry = sf::st_sfc(lapply(seq_along(ids), function(i)
      sf::st_linestring(rbind(c(i, 0), c(i + 1, 0)))), crs = 5070)))
}
# member -> rec_div area, the areas the collapse actually sums
.ci_area <- function(...) c(...)

cand <- function(net, area, ideal = 10, amin = 3, amax = 15)
  hydrofabric:::.build_hw_collapse_fp(net, min_area_sqkm = amin,
    min_length_km = 0.25, rec_touch = NULL, rec_area = area,
    ideal_size_sqkm = ideal, connector_ids = NULL, max_area_sqkm = amax)

test_that("two units that would exceed the target together do not merge", {
  # 9.9 + 9.9 = 19.8. Both are under the target on their own, which is what the
  # donor-only test looked at, and the result is nearly double it.
  net <- .ci_net(c("a", "b"), c("b", "0"), c("ma", "mb"))
  r <- cand(net, .ci_area(ma = 9.9, mb = 9.9))
  expect_equal(nrow(r), 0L)
})

test_that("two units that fit under the target do merge", {
  net <- .ci_net(c("a", "b"), c("b", "0"), c("ma", "mb"))
  r <- cand(net, .ci_area(ma = 4, mb = 5))   # 9 <= 10
  expect_equal(r$flowpath_id, "a")
  expect_equal(r$becomes, "b")
})

test_that("a merge landing exactly on the target is allowed", {
  net <- .ci_net(c("a", "b"), c("b", "0"), c("ma", "mb"))
  expect_equal(nrow(cand(net, .ci_area(ma = 4, mb = 6))), 1L)   # exactly 10
  expect_equal(nrow(cand(net, .ci_area(ma = 4.1, mb = 6))), 0L) # 10.1
})

test_that("a sub-floor unit merges even when the receiver is far past target", {
  # The floor is a mandate, not a preference: a unit under min_area_sqkm is not a
  # valid computational unit, while an oversized one is a valid unit of the wrong
  # size. There is one downstream receiver, so refusing is permanent.
  net <- .ci_net(c("a", "b"), c("b", "0"), c("ma", "mb"))
  r <- cand(net, .ci_area(ma = 1, mb = 40))
  expect_equal(r$flowpath_id, "a")
})

test_that("a sub-floor unit does merge when the receiver has room for it", {
  net <- .ci_net(c("a", "b"), c("b", "0"), c("ma", "mb"))
  r <- cand(net, .ci_area(ma = 1, mb = 6))   # 7 <= 10
  expect_equal(r$flowpath_id, "a")
})

test_that("a unit already past the target is never a donor", {
  net <- .ci_net(c("a", "b"), c("b", "0"), c("ma", "mb"))
  expect_equal(nrow(cand(net, .ci_area(ma = 12, mb = 1))), 0L)
})

test_that("a POI outlet is never collapsed away, at any size", {
  net <- .ci_net(c("a", "b"), c("b", "0"), c("ma", "mb"),
                 poi = c("gage-1", NA_character_))
  expect_equal(nrow(cand(net, .ci_area(ma = 1, mb = 2))), 0L)
})

test_that("area is summed over members, not read off the flowpath", {
  # 'a' carries two member reaches; together they are past the target, so no
  # merge, even though either alone would be under it.
  net <- .ci_net(c("a", "b"), c("b", "0"), c("m1,m2", "mb"))
  expect_equal(nrow(cand(net, .ci_area(m1 = 6, m2 = 6, mb = 1))), 0L)
})

test_that("a confluence is not a donor, so junctions survive", {
  # x and y both drain to z: z has two inflows and must not collapse into w.
  net <- .ci_net(c("x", "y", "z", "w"), c("z", "z", "w", "0"),
                 c("mx", "my", "mz", "mw"))
  r <- cand(net, .ci_area(mx = 1, my = 1, mz = 1, mw = 1))
  expect_false("z" %in% r$flowpath_id)
})

test_that("the floor still reaches a sliver whose neighbour is under target", {
  net <- .ci_net(c("a", "b"), c("b", "0"), c("ma", "mb"))
  r <- cand(net, .ci_area(ma = 0.5, mb = 2))
  expect_equal(r$flowpath_id, "a")
})

test_that("a receiver's headroom is spent once, not once per donor", {
  # Both donors drain into 'b'. Each pairwise sum is 9 and passes on its own;
  # applied together in the same round they would leave the receiver at 13.
  # Every accepted candidate IS applied in the same round, so the budget has to
  # be per receiver rather than per merge.
  net <- .ci_net(c("a", "c", "b", "d"), c("b", "b", "d", "0"),
                 c("ma", "mc", "mb", "md"), lp = c("1", "2", "1", "1"))
  ar  <- .ci_area(ma = 4, mc = 4, mb = 5, md = 1)
  ids <- names(ar)
  touch <- stats::setNames(lapply(ids, function(i) ids), ids)
  r <- hydrofabric:::.build_hw_collapse_fp(net, min_area_sqkm = 3,
    min_length_km = 0.25, rec_touch = touch, rec_area = ar,
    ideal_size_sqkm = 10, connector_ids = NULL, max_area_sqkm = 15)
  expect_equal(nrow(r), 1L)
  expect_equal(max(table(r$becomes)), 1L)
})

test_that("a sub-floor donor is still taken when the receiver has no headroom", {
  # The floor outranks the target, so a unit below min_area_sqkm merges even
  # though the receiver is already full.
  net <- .ci_net(c("a", "c", "b", "d"), c("b", "b", "d", "0"),
                 c("ma", "mc", "mb", "md"), lp = c("1", "2", "1", "1"))
  ar  <- .ci_area(ma = 1, mc = 4, mb = 9, md = 1)
  ids <- names(ar)
  touch <- stats::setNames(lapply(ids, function(i) ids), ids)
  r <- hydrofabric:::.build_hw_collapse_fp(net, min_area_sqkm = 3,
    min_length_km = 0.25, rec_touch = touch, rec_area = ar,
    ideal_size_sqkm = 10, connector_ids = NULL, max_area_sqkm = 15)
  expect_true("a" %in% r$flowpath_id)   # 1 km2, below the floor
  expect_false("c" %in% r$flowpath_id)  # 4 km2, would take the receiver to 13
})

test_that("the ceiling bounds ordinary merges and yields to the floor", {
  # 10 is what an ordinary merge aims at and stops at; 15 is the size we prefer
  # not to exceed, which a floor rescue is allowed to ignore because the
  # alternative is a unit that should not exist at all.
  net <- .ci_net(c("a", "b"), c("b", "0"), c("ma", "mb"))
  expect_equal(nrow(cand(net, .ci_area(ma = 2, mb = 9))),  1L)  # floor donor, 11
  expect_equal(nrow(cand(net, .ci_area(ma = 2, mb = 14))), 1L)  # floor donor, 16
  expect_equal(nrow(cand(net, .ci_area(ma = 4, mb = 9))),  0L)  # ordinary, 13
})

test_that("an ordinary donor gets no access to the band above the target", {
  # Only clearing the floor justifies passing 10; an ordinary merge stops there.
  net <- .ci_net(c("a", "b"), c("b", "0"), c("ma", "mb"))
  expect_equal(nrow(cand(net, .ci_area(ma = 4, mb = 9))), 0L)   # 13, not sub-floor
})

test_that("the upstream absorb is never blocked by the ceiling", {
  # A sub-floor POI outlet grows up its own mainstem, but not past the ceiling:
  # rescuing one criterion by missing the other is the trade the ceiling refuses.
  net <- .ci_net(c("u", "s"), c("s", "0"), c("mu", "ms"),
                 poi = c(NA_character_, "gage-9"), lp = c("1", "1"))
  ids   <- c("mu", "ms")
  touch <- stats::setNames(lapply(ids, function(i) ids), ids)
  small <- function(ar, amax) hydrofabric:::.build_absorb_small_fp(net,
    min_area_sqkm = 3, rec_touch = touch, rec_area = ar,
    connector_ids = NULL, max_area_sqkm = amax)
  # s is 2 km2 and a POI outlet; u is its mainstem upstream
  # This is the last path a fixed outlet has, so the ceiling reports and never
  # refuses: blocking it would strand the unit below the floor permanently.
  expect_equal(nrow(small(c(mu = 6,  ms = 2), 15)), 1L)   # 8, under the ceiling
  expect_equal(nrow(small(c(mu = 40, ms = 2), 15)), 1L)   # 42, taken and reported
})
