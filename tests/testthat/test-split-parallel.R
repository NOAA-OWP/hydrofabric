# Regression guard for the ffa1a97 (2026-06-15) `use_parallel <- FALSE`
# diagnostic that silently serialized reconcile_divides in every build.
# The divide splitter must derive parallelism from its `parallel` argument,
# never hard-disable it.

test_that("split_divides_d8 does not hard-disable parallelism", {
  fn  <- getFromNamespace("split_divides_d8", "hydrofabric")
  src <- paste(deparse(body(fn)), collapse = "\n")

  # No unconditional force-sequential override (the regression).
  expect_false(
    grepl("use_parallel\\s*<-\\s*FALSE", src),
    label = "split_divides_d8 contains a hard `use_parallel <- FALSE` override"
  )

  # `use_parallel` must be derived from the `parallel` argument.
  expect_true(
    grepl("use_parallel\\s*<-[^\\n]*\\bparallel\\b", src),
    label = "`use_parallel` should be derived from the `parallel` argument"
  )
})
