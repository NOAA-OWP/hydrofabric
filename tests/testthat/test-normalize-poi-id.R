# A poi_id is a single integer by contract, and the same value keys the join from
# a POI event to its reconciled reach. The AK OWP hydrolocation set ships a text
# handle instead. Renumbering has to be atomic across every layer that carries
# the key: an earlier fix rewrote one table and left the other, the join then
# matched nothing, and every POI was dropped from the product.

test_that("text ids are renumbered identically across both layers", {
  lyrs <- list(
    hydrolocations = data.frame(poi_id = c("wbin-2-409089", "wbin-3-113398",
                                           "wbin-2-409089", NA)),
    pois           = data.frame(poi_id = c("wbin-3-113398", "wbin-2-409089")))
  out <- hydrofabric:::.normalize_poi_id(lyrs)
  h <- out$hydrolocations$poi_id
  p <- out$pois$poi_id

  expect_false(anyNA(suppressWarnings(as.integer(h[!is.na(h)]))))
  # the same handle must land on the same integer in both layers
  expect_identical(h[1L], p[2L])
  expect_identical(h[2L], p[1L])
  # a repeated handle keeps one id, and NA stays NA
  expect_identical(h[1L], h[3L])
  expect_true(is.na(h[4L]))
})

test_that("renumbering cannot collide with an integer already in use", {
  lyrs <- list(
    hydrolocations = data.frame(poi_id = c("7", "12", "wbin-1-5")),
    pois           = data.frame(poi_id = "wbin-1-5"))
  out <- hydrofabric:::.normalize_poi_id(lyrs)
  v <- as.integer(out$hydrolocations$poi_id)
  expect_identical(v[1:2], c(7L, 12L))     # real ids untouched
  expect_gt(v[3L], 12L)                     # remapped above the highest in use
  expect_identical(out$pois$poi_id[1L], out$hydrolocations$poi_id[3L])
})

test_that("an all-integer or absent column is left alone", {
  lyrs <- list(hydrolocations = data.frame(poi_id = c("1", "2")),
               pois           = data.frame(poi_id = c("1", "2")))
  expect_identical(hydrofabric:::.normalize_poi_id(lyrs), lyrs)

  bare <- list(hydrolocations = data.frame(x = 1), pois = NULL)
  expect_identical(hydrofabric:::.normalize_poi_id(bare), bare)
})
