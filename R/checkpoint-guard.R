#' Version guard for stage checkpoint reuse
#'
#' `run_pipeline(use_checkpoints = TRUE)` resumes a build from existing stage
#' checkpoints (e.g. after an OOM/crash) instead of recomputing from scratch.
#' That is only safe when the checkpoints were produced by the *same* pipeline
#' build and params: otherwise a code fix (reinstall) or a param change would
#' silently reuse a stale/buggy intermediate. These helpers fingerprint the
#' current build+params and gate reuse on a match.
#'
#' The fingerprint is intentionally coarse: any mismatch (a reinstall flips the
#' `Built` stamp; a param edit flips the params component) forces a clean
#' recompute, which is always correct. Over-invalidation only costs time.
#'
#' @param params Optional list of stage params to fold into the fingerprint.
#' @return A single character string identifying the build + params.
#' @keywords internal
.pipeline_fingerprint <- function(params = NULL) {
  d <- tryCatch(
    utils::packageDescription("hydrofabric", fields = c("Version", "Built")),
    error = function(e) list(Version = "dev", Built = "dev"))
  parts <- c(as.character(d[["Version"]]), as.character(d[["Built"]]))
  if (!is.null(params))
    parts <- c(parts, tryCatch(
      paste(utils::capture.output(utils::str(params)), collapse = "|"),
      error = function(e) ""))
  paste(parts, collapse = "||")
}

#' @rdname dot-pipeline_fingerprint
#' @param work_dir Stage work directory holding the `.pipeline_fingerprint` stamp.
#' @param fp Current fingerprint from [.pipeline_fingerprint()].
#' @return `.checkpoints_valid`: `TRUE` iff `work_dir` was stamped with `fp`.
#' @keywords internal
.checkpoints_valid <- function(work_dir, fp) {
  f <- file.path(work_dir, ".pipeline_fingerprint")
  file.exists(f) && identical(readLines(f, warn = FALSE)[1L], fp)
}

#' @rdname dot-pipeline_fingerprint
#' @return `.stamp_fingerprint`: invisibly `NULL`; writes the stamp file.
#' @keywords internal
.stamp_fingerprint <- function(work_dir, fp) {
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(fp, file.path(work_dir, ".pipeline_fingerprint"))
  invisible(NULL)
}
