# check_reference_downstream_resolves finds a reach that names a downstream the
# file does not contain. A true outlet names none; a truncated reach still
# carries the pointer to the neighbour that was removed, and becomes a terminal
# nexus whose flow is never routed on. The rate is domain-dependent.

.crd_gpkg <- function(d) {
  p <- tempfile(fileext = ".gpkg")
  d$geom <- sf::st_sfc(lapply(seq_len(nrow(d)), function(i)
    sf::st_linestring(rbind(c(i, 0), c(i + 1, 0)))), crs = 5070)
  suppressWarnings(sf::st_write(sf::st_as_sf(d), p, "flowpaths", quiet = TRUE))
  p
}

test_that("a toid naming an absent reach is reported", {
  p <- .crd_gpkg(data.frame(flowpath_id = c(1, 2), flowpath_toid = c(2, 99),
    dnhydroseq = c(0, 0), hydroseq = c(10, 20)))
  r <- suppressMessages(hydrofabric::check_reference_downstream_resolves(p))
  expect_equal(r$id, "2")
  expect_equal(r$names_id, "99")
  expect_equal(r$dangling, "flowpath_toid")
})

test_that("a true outlet is not reported", {
  # toid 0 and no dnhydroseq: the network legitimately stops here.
  p <- .crd_gpkg(data.frame(flowpath_id = c(1, 2), flowpath_toid = c(2, 0),
    dnhydroseq = c(0, 0), hydroseq = c(10, 20)))
  r <- suppressMessages(hydrofabric::check_reference_downstream_resolves(p))
  expect_equal(nrow(r), 0L)
})

test_that("a truncated reach is caught by its dangling dnhydroseq", {
  # toid says outlet, but dnhydroseq still names a reach that is gone: the
  # signature of a network cut at the edge of its coverage.
  p <- .crd_gpkg(data.frame(flowpath_id = c(1, 2), flowpath_toid = c(2, 0),
    dnhydroseq = c(20, 777), hydroseq = c(10, 20)))
  r <- suppressMessages(hydrofabric::check_reference_downstream_resolves(p))
  expect_equal(r$id, "2")
  expect_equal(r$dangling, "dnhydroseq")
})

test_that("a dnhydroseq that does resolve is not reported", {
  p <- .crd_gpkg(data.frame(flowpath_id = c(1, 2), flowpath_toid = c(2, 0),
    dnhydroseq = c(20, 10), hydroseq = c(10, 20)))
  r <- suppressMessages(hydrofabric::check_reference_downstream_resolves(p))
  expect_equal(nrow(r), 0L)
})

test_that("a live toid suppresses the dnhydroseq test for that reach", {
  # The reach drains somewhere real; dnhydroseq disagreeing is not a truncation.
  p <- .crd_gpkg(data.frame(flowpath_id = c(1, 2), flowpath_toid = c(2, 1),
    dnhydroseq = c(555, 10), hydroseq = c(10, 20)))
  r <- suppressMessages(hydrofabric::check_reference_downstream_resolves(p))
  expect_equal(nrow(r), 0L)
})

test_that("the downstream search spans the whole file, not just the vpu", {
  # Reach 2 is in VPU 04 and drains into reach 3 in VPU 05. Restricting the
  # search to VPU 04 would call that a truncation; it is ordinary.
  p <- .crd_gpkg(data.frame(flowpath_id = c(1, 2, 3),
    flowpath_toid = c(2, 3, 0), dnhydroseq = c(0, 0, 0),
    hydroseq = c(10, 20, 30), vpuid = c("04", "04", "05")))
  r <- suppressMessages(
    hydrofabric::check_reference_downstream_resolves(p, vpuid = "04"))
  expect_equal(nrow(r), 0L)
})

test_that("a vpuid filter restricts which reaches are reported", {
  p <- .crd_gpkg(data.frame(flowpath_id = c(1, 2), flowpath_toid = c(88, 99),
    dnhydroseq = c(0, 0), hydroseq = c(10, 20), vpuid = c("04", "18")))
  r <- suppressMessages(
    hydrofabric::check_reference_downstream_resolves(p, vpuid = "04"))
  expect_equal(r$id, "1")
})

test_that("a layer without the hydroseq columns still checks toid", {
  p <- .crd_gpkg(data.frame(flowpath_id = c(1, 2), flowpath_toid = c(2, 99)))
  r <- suppressMessages(hydrofabric::check_reference_downstream_resolves(p))
  expect_equal(r$id, "2")
})

test_that("an overflowed dnhydroseq is refused, not counted", {
  # The CONUS reference stores dnhydroseq as a 32-bit integer while hydroseq runs
  # past 13 billion, so most pointers have wrapped negative. Counting those as
  # unresolved would report an overflow as a truncated network.
  p <- .crd_gpkg(data.frame(flowpath_id = c(1, 2), flowpath_toid = c(2, 0),
    dnhydroseq = c(-116718348, -1154826488), hydroseq = c(1e10, 1.3e10)))
  r <- suppressMessages(hydrofabric::check_reference_downstream_resolves(p))
  expect_equal(nrow(r), 0L)
})

test_that("a sound dnhydroseq is still used", {
  p <- .crd_gpkg(data.frame(flowpath_id = c(1, 2), flowpath_toid = c(2, 0),
    dnhydroseq = c(20, 777), hydroseq = c(10, 20)))
  r <- suppressMessages(hydrofabric::check_reference_downstream_resolves(p))
  expect_equal(r$dangling, "dnhydroseq")
})
