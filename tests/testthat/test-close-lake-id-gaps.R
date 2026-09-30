# t-route collapses every flowpath sharing a lake_id into one node, so those
# flowpaths must form a connected subgraph. Membership is a length-weighted
# majority of member reaches and cannot see the network, so a reach mostly
# outside the lake is rejected even when both neighbours are in it.
#
# Step 6h closes those gaps, but 6k later deletes divide-less leaf flowpaths and
# rewires flowpath_toid, which can open new ones. These
# tests cover the closure running against settled topology.

make_gpkg <- function(lake_ids, toids, dir = tempfile()) {
  dir.create(dir, showWarnings = FALSE)
  g <- file.path(dir, "t.gpkg")
  n <- length(lake_ids)
  ids <- paste0("fp-", seq_len(n))
  pt <- function(i) sf::st_sfc(sf::st_linestring(rbind(c(i, 0), c(i + 1, 0))), crs = 5070)
  fp <- sf::st_sf(
    flowpath_id   = ids,
    flowpath_toid = toids,
    lake_id       = lake_ids,
    geom = do.call(c, lapply(seq_len(n), pt))
  )
  # one nexus per flowpath, nex-i sits between fp-i and its successor
  nex_ids <- unique(toids[!is.na(toids)])
  nx <- sf::st_sf(
    nexus_id   = nex_ids,
    nexus_toid = sub("^nex-", "fp-", nex_ids),
    geom = do.call(c, lapply(seq_along(nex_ids), function(i) sf::st_sfc(sf::st_point(c(i, 0)), crs = 5070)))
  )
  sf::st_write(fp, g, "flowpaths", quiet = TRUE)
  sf::st_write(nx, g, "nexus", quiet = TRUE)
  g
}

read_lakes <- function(g) as.character(sf::st_read(g, "flowpaths", quiet = TRUE)$lake_id)

test_that("a single non-member between two members of one lake is closed", {
  # fp-1(L) -> fp-2(NA) -> fp-3(L)
  g <- make_gpkg(c("L", NA, "L"), c("nex-2", "nex-3", NA))
  out <- hydrofabric:::close_lake_id_gaps(g)
  expect_identical(out, "fp-2")
  expect_identical(read_lakes(g), c("L", "L", "L"))
})

test_that("a fabric with no gaps is left unchanged", {
  g <- make_gpkg(c("L", "L", NA), c("nex-2", "nex-3", NA))
  before <- read_lakes(g)
  out <- hydrofabric:::close_lake_id_gaps(g)
  expect_length(out, 0)
  expect_identical(read_lakes(g), before)
})

test_that("a gap between two DIFFERENT lakes is not closed", {
  # fp-1(A) -> fp-2(NA) -> fp-3(B): not one waterbody, must stay a gap
  g <- make_gpkg(c("A", NA, "B"), c("nex-2", "nex-3", NA))
  out <- hydrofabric:::close_lake_id_gaps(g)
  expect_length(out, 0)
  expect_identical(read_lakes(g)[2], NA_character_)
})

test_that("a run longer than max_step is reported rather than absorbed", {
  # four non-members between two members: a genuine river between two arms
  g <- make_gpkg(c("L", NA, NA, NA, NA, "L"),
                 c("nex-2", "nex-3", "nex-4", "nex-5", "nex-6", NA))
  expect_length(hydrofabric:::close_lake_id_gaps(g, max_step = 3L), 0)
  # with the bound raised, the same gap closes
  expect_length(hydrofabric:::close_lake_id_gaps(g, max_step = 6L), 4)
})

test_that("flowlines and network inherit the closed stamp", {
  g <- make_gpkg(c("L", NA, "L"), c("nex-2", "nex-3", NA))
  fl <- sf::st_sf(flowline_id = c("fl-1", "fl-2"), flowpath_id = c("fp-2", "fp-3"),
                  lake_id = c(NA_character_, "L"),
                  geom = sf::st_sfc(sf::st_linestring(rbind(c(0,0), c(1,0))),
                                    sf::st_linestring(rbind(c(1,0), c(2,0))), crs = 5070))
  sf::st_write(fl, g, "flowlines", quiet = TRUE)
  net <- data.frame(flowpath_id = c("fp-1","fp-2","fp-3"), lake_id = c("L", NA, "L"))
  sf::st_write(net, g, "network", quiet = TRUE)

  hydrofabric:::close_lake_id_gaps(g)
  expect_identical(as.character(sf::st_read(g, "flowlines", quiet = TRUE)$lake_id), c("L","L"))
  expect_identical(as.character(sf::st_read(g, "network", quiet = TRUE)$lake_id), c("L","L","L"))
})

test_that("the closure satisfies the lake_members_contiguous rule it exists for", {
  # Assert the gate's own condition directly: no non-member may have a member
  # upstream and a member downstream of the same lake.
  g  <- make_gpkg(c("L", NA, "L"), c("nex-2", "nex-3", NA))
  hydrofabric:::close_lake_id_gaps(g)
  fp <- sf::st_read(g, "flowpaths", quiet = TRUE)
  nx <- sf::st_read(g, "nexus", quiet = TRUE)
  fid <- as.character(fp$flowpath_id)
  lk  <- as.character(fp$lake_id); lk[!nzchar(lk)] <- NA_character_
  n2f <- stats::setNames(as.character(nx$nexus_toid), as.character(nx$nexus_id))
  nxt <- match(unname(n2f[as.character(fp$flowpath_toid)]), fid)
  ups <- split(seq_along(fid), nxt)
  bad <- character(0)
  for (i in which(is.na(lk))) {
    dn <- nxt[i]
    if (is.na(dn) || is.na(lk[dn])) next
    u <- ups[[as.character(i)]]
    if (is.null(u)) next
    if (any(!is.na(lk[u]) & lk[u] == lk[dn])) bad <- c(bad, fid[i])
  }
  expect_length(bad, 0)
})
