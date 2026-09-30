test_that("divide_flowpath_consistent gates all three copies of the link", {
  # The flowpath<->divide link is stored redundantly in flowpaths.divide_id,
  # divides.flowpath_id and network.divide_id. This check replaced a finalize
  # repair step (.reconcile_divide_ids), so it must fail on a disagreement in
  # any of the three, not just the two spatial layers.
  fp <- data.frame(
    flowpath_id = paste0("fp-", 1:3),
    divide_id   = paste0("cat-", 1:3),
    stringsAsFactors = FALSE)
  dv <- data.frame(
    divide_id   = paste0("cat-", 1:3),
    flowpath_id = paste0("fp-", 1:3),
    stringsAsFactors = FALSE)
  # flowline_id / flowline_toid are part of the real network schema and are read
  # by other checks in the same function, so the fixture carries them.
  net <- data.frame(
    flowpath_id   = paste0("fp-", c(1, 1, 2, 3)),
    divide_id     = paste0("cat-", c(1, 1, 2, 3)),
    flowline_id   = as.character(1:4),
    flowline_toid = c("2", "3", "3", "0"),
    stringsAsFactors = FALSE)

  res <- hydrofabric:::.hf_check_ngen(
    list(flowpaths = fp, divides = dv, network = net), FALSE, 0)
  expect_true(res$divide_flowpath_consistent$ok)
  expect_identical(res$divide_flowpath_consistent$kind, "check")

  # (1) divides.flowpath_id points at a different flowpath
  dv_bad <- dv
  dv_bad$flowpath_id[2] <- "fp-3"
  expect_false(hydrofabric:::.hf_check_ngen(
    list(flowpaths = fp, divides = dv_bad, network = net), FALSE, 0
  )$divide_flowpath_consistent$ok)

  # (2) a network row disagrees with its own flowpath's divide_id
  net_bad <- net
  net_bad$divide_id[2] <- "cat-2"
  expect_false(hydrofabric:::.hf_check_ngen(
    list(flowpaths = fp, divides = dv, network = net_bad), FALSE, 0
  )$divide_flowpath_consistent$ok)

  # (3) network divide_id null where the flowpath has one
  net_na <- net
  net_na$divide_id[3] <- NA_character_
  expect_false(hydrofabric:::.hf_check_ngen(
    list(flowpaths = fp, divides = dv, network = net_na), FALSE, 0
  )$divide_flowpath_consistent$ok)

  # Absent network layer must not fail the check, it is optional.
  expect_true(hydrofabric:::.hf_check_ngen(
    list(flowpaths = fp, divides = dv), FALSE, 0
  )$divide_flowpath_consistent$ok)
})
