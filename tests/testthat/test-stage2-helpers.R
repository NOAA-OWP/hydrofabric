## Stage-2 helper unit tests.
## .small_segments_layer must reproduce the inline Step 4e filter.


suppressMessages(library(sf))

# Horizontal linestring from (x0,y) to (x1,y): length = x1 - x0.
.hline <- function(x0, x1, y) st_linestring(matrix(c(x0, x1, y, y), 2, 2))
# Axis-aligned rectangle polygon over [xmin,xmax] x [ymin,ymax].
.rect  <- function(xmin, xmax, ymin, ymax)
  st_polygon(list(matrix(c(xmin, xmax, xmax, xmin, xmin,
    ymin, ymin, ymax, ymax, ymin), 5, 2)))

# Reference implementation: the ORIGINAL Step 4a0 per-group selection loop
# (sf::st_intersects per candidate). The vectorized helper must match this.
.old_pick <- function(group_vec, lens_all, flat_sfc, div_sfc_by_group, has_div) {
  vapply(seq_len(length(has_div)), function(k) {
    grp <- which(group_vec == k)
    if (length(grp) == 1L) return(grp[1L])
    longest <- grp[which.max(lens_all[grp])]
    ll <- max(lens_all[grp])
    if (isTRUE(has_div[k]) && !is.null(div_sfc_by_group[[k]])) {
      nl   <- grp[grp != longest]
      cand <- nl[lens_all[nl] >= 0.5 * ll]
      if (length(cand)) {
        hits <- vapply(cand, function(j)
          isTRUE(st_intersects(flat_sfc[j], div_sfc_by_group[[k]],
            sparse = FALSE)[1, 1]), logical(1L))
        if (any(hits)) {
          hc <- cand[hits]
          return(hc[which.max(lens_all[hc])])
        }
      }
    }
    longest
  }, integer(1L))
}



test_that("Step 4e: small_segments filter + classification", {
  fp <- st_sf(
    flowpath_id   = as.character(1:5),
    flowpath_toid = c("2", "3", "4", "1000000001", "0"),
    areasqkm      = c(1, 5, 10, 2, 4),     # rows 1,4 small by area
    lengthkm      = c(2, 0.1, 3, 1, 0.2),  # rows 2,5 small by length
    poi_id        = c(NA, NA, NA, NA, "p9"),
    geometry      = st_sfc(lapply(1:5, function(i) .hline(0, 1, i)))
  )
  small <- hydrofabric:::.small_segments_layer(fp)
  # rows kept: 1 (area<3), 2 (len<0.25), 4 (area<3), 5 (len<0.25)
  expect_setequal(small$flowpath_id, c("1", "2", "4", "5"))
  st <- setNames(small$small_type, small$flowpath_id)
  expect_equal(unname(st["5"]), "poi_protected")  # poi_id present wins
  expect_equal(unname(st["4"]), "terminal")        # toid >= 1e9
  expect_equal(unname(st["1"]), "other")
  expect_equal(unname(st["2"]), "other")
})
