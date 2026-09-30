test_that(".separate_longer explodes a delimited column, recycling other cols", {
  df <- data.frame(id = c(1L, 2L, 3L),
    x  = c("a,b,c", "d", "e,f"),
    stringsAsFactors = FALSE)

  out <- hydrofabric:::.separate_longer(df, "x", ",")
  expect_equal(nrow(out), 6L)
  expect_equal(out$x,  c("a", "b", "c", "d", "e", "f"))   # tokens in order
  expect_equal(out$id, c(1L, 1L, 1L, 2L, 3L, 3L))         # other cols recycled
  expect_equal(rownames(out), as.character(seq_len(6L)))  # row names reset

  # row count == total token count
  expect_equal(nrow(out),
    sum(lengths(strsplit(df$x, ",", fixed = TRUE))))

  # single token / no delimiter is a no-op on row count
  expect_equal(hydrofabric:::.separate_longer(data.frame(x = "a"), "x", ",")$x, "a")
  expect_equal(nrow(hydrofabric:::.separate_longer(data.frame(x = c("a", "b")), "x", ",")), 2L)
})
