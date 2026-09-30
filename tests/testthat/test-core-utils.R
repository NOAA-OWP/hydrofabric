## Unit tests for foundational pure helpers: lineage CSV, DAG/hydroseq topology
## utilities, and deterministic parameter hashing.
## (id formatting moved to hfutils::hf_fmt_id: tested in hfutils/test-fmt-id.R;
##  gpkg write-primitives moved to hfutils::gpkg_update_col/gpkg_update_geom/gpkg_exec.)

test_that("lineage member-CSV helpers round-trip and normalize", {
  expect_identical(hydrofabric:::expand_member_csv("a, b ,c"), c("a", "b", "c"))   # trims
  expect_identical(hydrofabric:::expand_member_csv(c("a,b", "b,c"), unique = TRUE), c("a", "b", "c"))
  expect_identical(hydrofabric:::expand_member_csv(c(NA, "", "x")), "x")           # drops NA/empty
  expect_identical(hydrofabric:::expand_member_csv(character(0)), character(0))

  expect_identical(hydrofabric:::compose_member_csv(c("a", "b", "b", "")), "a,b")  # dedupe + drop empty
  expect_identical(hydrofabric:::compose_member_csv(character(0)), NA_character_)

  # round-trip
  expect_identical(hydrofabric:::compose_member_csv(hydrofabric:::expand_member_csv("a,b,c")), "a,b,c")
})

test_that("member-part parse/strip split id and .part", {
  p <- hydrofabric:::parse_member_part(c("123.2", "456", "789.10"))
  expect_identical(p$hf_id, c("123", "456", "789"))
  expect_identical(p$hf_id_part, c(2L, 1L, 10L))       # no suffix -> part 1
})

test_that("hfutils::hf_network_is_dag detects cycles; terminals are not edges", {
  chain <- data.frame(flowpath_id = c("1", "2", "3"),
    flowpath_toid = c("2", "3", "0"), stringsAsFactors = FALSE)
  expect_true(hfutils::hf_network_is_dag(chain))
  # terminal sentinels (0 / NA / >=1e9 / not-in-ids) must not create self-loops
  term <- data.frame(flowpath_id = c("1", "2"),
    flowpath_toid = c("2", "1000000001"), stringsAsFactors = FALSE)
  expect_true(hfutils::hf_network_is_dag(term))
  cyc <- data.frame(flowpath_id = c("1", "2"),
    flowpath_toid = c("2", "1"), stringsAsFactors = FALSE)
  expect_false(hfutils::hf_network_is_dag(cyc))
})

test_that("hfutils::hf_network_is_dag accepts a nondendritic (divergent) network", {
  # Nondendritic: flowpath 1 splits to two downstream paths (distributary /
  # braid / delta). Out-degree > 1 is acyclic but not a simple tree, it must
  # be accepted as a valid DAG, not mistaken for a cycle or error. (T-Route
  # requires a DAG; the AC names acyclic + nondendritic tests explicitly.)
  divergent <- data.frame(
    flowpath_id   = c("1", "1", "2", "3"),
    flowpath_toid = c("2", "3", "0", "0"), stringsAsFactors = FALSE)
  expect_true(hfutils::hf_network_is_dag(divergent))

  # A divergence that also reconverges downstream stays acyclic (1->2, 1->3,
  # 2->4, 3->4): diamond, nondendritic, still a DAG.
  diamond <- data.frame(
    flowpath_id   = c("1", "1", "2", "3", "4"),
    flowpath_toid = c("2", "3", "4", "4", "0"), stringsAsFactors = FALSE)
  expect_true(hfutils::hf_network_is_dag(diamond))
})

test_that("hfutils::hf_recompute_hydroseq puts outlets lowest (NHD convention)", {
  chain <- data.frame(flowpath_id = c("1", "2", "3"),
    flowpath_toid = c("2", "3", "0"), stringsAsFactors = FALSE)
  hs <- hfutils::hf_recompute_hydroseq(chain)$hydroseq
  expect_length(hs, 3L)
  expect_false(anyNA(hs))
  expect_equal(hs[3], min(hs))   # fp 3 is the outlet -> lowest hydroseq
  expect_equal(hs[1], max(hs))   # fp 1 is the headwater -> highest
})

test_that("hfutils::hf_break_cycles makes a cyclic network acyclic", {
  cyc <- data.frame(flowpath_id = c("1", "2"),
    flowpath_toid = c("2", "1"), stringsAsFactors = FALSE)
  expect_false(hfutils::hf_network_is_dag(cyc))
  broken <- suppressMessages(hfutils::hf_break_cycles(cyc))
  expect_true(hfutils::hf_network_is_dag(broken))
  expect_equal(nrow(broken), 2L)   # severs an edge, never drops a node
})

# parameter hashing (canonical_params_json / params_sha) moved out of this package.

test_that(".hf_attribution credits TDX and HydroSHEDS except on HI and PRVI", {
  old <- Sys.getenv("HF_DATA_ATTRIBUTION", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("HF_DATA_ATTRIBUTION")
  else Sys.setenv(HF_DATA_ATTRIBUTION = old), add = TRUE)
  Sys.unsetenv("HF_DATA_ATTRIBUTION")

  # Every domain but Hawaii and PRVI touches TDX-Hydro and HydroSHEDS in some
  # capacity, so both credits are owed alongside CC BY-SA 4.0.
  expect_match(hydrofabric:::.hf_attribution("19"), "NGA", fixed = TRUE)
  expect_match(hydrofabric:::.hf_attribution("19"), "CC BY-SA 4.0", fixed = TRUE)
  expect_match(hydrofabric:::.hf_attribution("09"), "NGA", fixed = TRUE)
  expect_match(hydrofabric:::.hf_attribution("09"), "HydroSHEDS", fixed = TRUE)

  # Hawaii and PRVI are NHDPlus HR throughout, a US Government work in the
  # public domain, so they carry the Lynker notice alone.
  expect_identical(hydrofabric:::.hf_attribution("20"), hydrofabric:::.HF_ATTR_LS)
  expect_identical(hydrofabric:::.hf_attribution("21"), hydrofabric:::.HF_ATTR_LS)

  # An unknown vpuid credits both rather than falling through uncredited.
  expect_match(hydrofabric:::.hf_attribution(NA_character_), "NGA", fixed = TRUE)

  # An explicit operator override always wins over the source-derived default.
  Sys.setenv(HF_DATA_ATTRIBUTION = "Operator Override")
  expect_identical(hydrofabric:::.hf_attribution("19"), "Operator Override")
})
