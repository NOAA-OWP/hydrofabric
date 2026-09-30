## Unit tests for the nested-set upstream index (add_upstream_index).
##
## The contract: for every flowpath, the range query
##   upstream_id in (upstream_id[X], upstream_id[X] + num_upstreams[X]]
## must return exactly the set of flowpaths reachable by walking flowpath_toid
## upstream, and num_upstreams must equal that set's size. These tests build a
## small synthetic fp -> nexus -> fp tree and cross-check the index against a
## brute-force upstream walk, then confirm the divergence tripwire fires.

# Brute-force upstream set by walking the fp -> nexus -> fp graph.
.walk_upstream <- function(fp, nex, x) {
  id <- as.character(fp$flowpath_id)
  nexmap <- stats::setNames(as.character(nex$nexus_toid), as.character(nex$nexus_id))
  ds <- unname(nexmap[as.character(fp$flowpath_toid)]); ds[!(ds %in% id)] <- NA
  children <- split(id, ds)
  acc <- character(0); stk <- children[[x]]
  while (length(stk)) { h <- stk[1]; stk <- stk[-1]; acc <- c(acc, h); stk <- c(stk, children[[h]]) }
  acc
}

test_that("range query equals the upstream walk on a branching tree", {
  #        7   (outlet)
  #        |
  #        6
  #       / \
  #      5   2
  #     /|   |
  #    4 3   1
  fp <- data.frame(
    flowpath_id  = as.character(1:7),
    flowpath_toid = paste0("nex-", c(2, 6, 5, 5, 6, 7, "T")),
    stringsAsFactors = FALSE)
  nex <- data.frame(
    nexus_id  = paste0("nex-", 1:7),
    nexus_toid = as.character(1:7),
    stringsAsFactors = FALSE)

  idx <- hfutils::hf_upstream_index(fp, nex)

  expect_equal(attr(idx, "n_bad"), 0L)
  expect_equal(attr(idx, "n_outlets"), 1L)
  expect_equal(nrow(idx), 7L)
  expect_equal(length(unique(idx$upstream_id)), 7L)   # a valid pre-order labelling

  u  <- stats::setNames(idx$upstream_id, idx$flowpath_id)
  nu <- stats::setNames(idx$num_upstreams, idx$flowpath_id)
  for (x in idx$flowpath_id) {
    lo <- u[[x]]
    rng <- idx$flowpath_id[idx$upstream_id > lo & idx$upstream_id <= lo + nu[[x]]]
    walk <- .walk_upstream(fp, nex, x)
    expect_setequal(rng, walk)
    expect_equal(nu[[x]], length(walk))
  }

  # spot-check the known counts: 6 has {5,4,3,2,1}, 7 has all 6 others
  expect_equal(nu[["7"]], 6L)
  expect_equal(nu[["6"]], 5L)
  expect_equal(nu[["1"]], 0L)   # headwater
})

test_that("a divergent nexus (two downstreams) trips the tree guard", {
  # Nexus 9 has two distinct downstreams (flowpath 2 and flowpath 3): a
  # divergence the fp -> nexus -> fp resolution cannot represent as a tree.
  fp <- data.frame(
    flowpath_id   = as.character(1:3),
    flowpath_toid = paste0("nex-", c(9, "T", "T")),
    stringsAsFactors = FALSE)
  nex <- data.frame(
    nexus_id   = c("nex-9", "nex-9"),
    nexus_toid = c("2", "3"),
    stringsAsFactors = FALSE)
  expect_warning(hfutils::hf_upstream_index(fp, nex), "divergence")
  idx <- suppressWarnings(hfutils::hf_upstream_index(fp, nex))
  expect_gt(attr(idx, "n_divergences"), 0L)
})

test_that("a duplicated flowpath_id trips the tree guard", {
  fp <- data.frame(
    flowpath_id   = c("1", "1", "2"),      # 1 duplicated
    flowpath_toid = paste0("nex-", c(2, 2, "T")),
    stringsAsFactors = FALSE)
  nex <- data.frame(
    nexus_id   = paste0("nex-", 1:2),
    nexus_toid = as.character(1:2),
    stringsAsFactors = FALSE)
  expect_warning(hfutils::hf_upstream_index(fp, nex), "not a tree")
})

test_that("upstream_index_monotone invariant passes on a clean indexed fabric", {
  skip_if_not(requireNamespace("sf", quietly = TRUE))
  fx <- testthat::test_path("..", "fixtures", "merge_min", "ngen_vpu01.gpkg")
  skip_if_not(file.exists(fx), "merge_min fixture missing")
  g <- withr::local_tempfile(fileext = ".gpkg")
  file.copy(fx, g, overwrite = TRUE)
  suppressMessages(add_upstream_index(g))
  fp <- sf::st_read(g, "flowpaths", quiet = TRUE)
  nx <- sf::st_read(g, "nexus", quiet = TRUE)

  res <- hydrofabric:::.hf_check_ngen(list(flowpaths = fp, nexus = nx), FALSE, 0)
  expect_true(res$upstream_index_monotone$ok)

  # corrupt one edge: set a downstream flowpath's num_upstreams below its
  # upstream neighbour -> the invariant must catch it
  fpid <- as.character(fp$flowpath_id)
  nexmap <- stats::setNames(as.character(nx$nexus_toid), as.character(nx$nexus_id))
  di <- match(unname(nexmap[as.character(fp$flowpath_toid)]), fpid)
  k <- which(!is.na(di))[1]
  fp$num_upstreams[di[k]] <- fp$num_upstreams[k] - 1L
  res2 <- hydrofabric:::.hf_check_ngen(list(flowpaths = fp, nexus = nx), FALSE, 0)
  expect_false(res2$upstream_index_monotone$ok)
})
