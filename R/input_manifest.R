# Recording which files a build actually read.
#
# The three input gates catch damaged or mismatched inputs, but none of them can
# tell you whether the file in front of you is the file a result was produced
# from. Inputs get regenerated in place under the same name, the hydrolocation
# set and the CONUS reference both did during this delivery, and nothing in a
# GeoPackage says which revision it was. A manifest closes that: it records what
# was read, so a later run can be checked against it rather than assumed
# equivalent.
#
# The hash is MD5 via tools::md5sum, which is base R and needs no dependency.
# This detects a different or changed file, which is the whole job; it is not a
# defence against a deliberately forged one.

.file_fingerprint <- function(path) {
  ok <- !is.na(path) & nzchar(path) & file.exists(path)
  data.frame(
    path  = path,
    exists = ok,
    bytes = ifelse(ok, file.size(path), NA_real_),
    mtime = ifelse(ok, format(file.mtime(path), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      NA_character_),
    md5   = ifelse(ok, unname(tools::md5sum(path)), NA_character_),
    stringsAsFactors = FALSE)
}

#' Record the inputs a build read
#'
#' @description
#' Fingerprints each supplied file, size, modification time and MD5, and
#' optionally writes the result as JSON. Pair it with [hf_verify_manifest()] to
#' confirm a later run is reading the same inputs.
#'
#' Hashing streams the whole file, so a multi-gigabyte reference takes a little
#' time. It is worth it once per build rather than per VPU.
#'
#' @param ... Named file paths, e.g. `reference = "...", fdr = "...", fac = "..."`.
#'   Names become the manifest's roles. `NULL` entries are dropped.
#' @param out Optional path to write the manifest as JSON. `NULL` (default)
#'   returns it without writing.
#' @returns A data frame with one row per input: `role`, `path`, `exists`,
#'   `bytes`, `mtime`, `md5`.
#' @examples
#' f <- tempfile(); writeLines("x", f)
#' hf_input_manifest(reference = f)
#' @export
hf_input_manifest <- function(..., out = NULL) {
  args <- list(...)
  args <- args[!vapply(args, is.null, logical(1L))]
  if (!length(args)) stop("hf_input_manifest: no inputs supplied", call. = FALSE)
  if (is.null(names(args)) || any(!nzchar(names(args))))
    stop("hf_input_manifest: every input must be named", call. = FALSE)

  m <- do.call(rbind, lapply(names(args), function(nm) {
    fp <- .file_fingerprint(as.character(args[[nm]]))
    cbind(role = rep(nm, nrow(fp)), fp, stringsAsFactors = FALSE)
  }))
  rownames(m) <- NULL

  missing <- m$role[!m$exists]
  if (length(missing))
    .hf_log_warn("hf_input_manifest: not found -- ", paste(missing, collapse = ", "))

  if (!is.null(out)) {
    jsonlite::write_json(list(
      software    = paste("hydrofabric", utils::packageVersion("hydrofabric")),
      recorded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      inputs      = m), out, auto_unbox = TRUE, pretty = TRUE)
    message("[manifest] wrote ", out)
  }
  m
}

#' Verify inputs against a recorded manifest
#'
#' @description
#' Re-fingerprints the files a manifest names and reports each as `ok`,
#' `changed` (the hash differs) or `missing`. Use it before a rebuild to
#' establish that the inputs are the ones a previous result came from.
#'
#' @param manifest Path to a manifest written by [hf_input_manifest()], or the
#'   data frame it returned.
#' @param strict If `TRUE`, any file that changed or went missing raises an
#'   error. Default `FALSE` reports and returns.
#' @returns A data frame with `role`, `path`, `status`, and the recorded and
#'   current hashes. Returned invisibly when everything matches.
#' @examples
#' f <- tempfile(); writeLines("x", f)
#' m <- hf_input_manifest(reference = f)
#' hf_verify_manifest(m)
#' @export
hf_verify_manifest <- function(manifest, strict = FALSE) {
  rec <- if (is.character(manifest) && length(manifest) == 1L) {
    j <- jsonlite::fromJSON(manifest)
    j$inputs
  } else manifest
  stopifnot(is.data.frame(rec), all(c("role", "path", "md5") %in% names(rec)))

  now <- .file_fingerprint(rec$path)
  status <- ifelse(!now$exists, "missing",
    ifelse(is.na(rec$md5) | now$md5 == rec$md5, "ok", "changed"))
  out <- data.frame(role = rec$role, path = rec$path, status = status,
    recorded_md5 = rec$md5, current_md5 = now$md5, stringsAsFactors = FALSE)

  for (i in seq_len(nrow(out)))
    message(sprintf("[manifest] %-8s %-16s %s", out$status[i], out$role[i],
      basename(out$path[i])))

  bad <- out$status != "ok"
  if (any(bad)) {
    msg <- sprintf("hf_verify_manifest: %d input(s) no longer match: %s",
      sum(bad), paste(out$role[bad], collapse = ", "))
    if (strict) stop(msg, call. = FALSE) else .hf_log_warn(msg)
    return(out)
  }
  message(sprintf("[manifest] all %d input(s) match", nrow(out)))
  invisible(out)
}
