#' Hydrofabric lineage CSV utilities
#'
#' Throughout the refactor and aggregate phases several columns store multi-valued lineage as a
#' comma-separated string: \code{member_comid} (CSV of original NHDPlus COMIDs,
#' possibly carrying a \code{.N} split suffix per child), and the
#' (mis-named) \code{flowline_id} column on Stage-2 \code{network} / \code{flowpaths}
#' (CSV of Stage-1 \code{rec_id}s).  These helpers provide one canonical
#' parser/composer so the dozens of inline \code{strsplit()} + \code{trimws()}
#' + \code{nzchar()} sites collapse to a single implementation.
#'
#' \itemize{
#'   \item \code{expand_member_csv(x)}: vector or single string of CSV -> flat
#'         character vector with whitespace stripped and empty tokens dropped.
#'         Use \code{unique=TRUE} to deduplicate.
#'   \item \code{compose_member_csv(v)}: character vector -> single CSV string
#'         with deterministic ordering (input order preserved, duplicates
#'         dropped on first occurrence).
#'   \item \code{parse_member_part(x)}: split a COMID lineage token of the form
#'         \code{"COMID"} or \code{"COMID.N"} into base \code{hf_id} (numeric
#'         COMID as character) and \code{hf_id_part} (integer N, defaulting to
#'         1L when absent).  Vectorised.
#' }
#'
#' These functions are intentionally pure: no I/O, no spatial side effects.
#'
#' @name lineage_helpers
#' @keywords internal
NULL


#' @rdname lineage_helpers
#' @param x A character vector of CSV strings, or a single CSV string.
#' @param unique If TRUE, deduplicate the expanded vector (input order preserved).
#' @return Character vector of tokens.
#' @keywords internal
expand_member_csv <- function(x, unique = FALSE) {
  if (length(x) == 0L) return(character(0L))
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) return(character(0L))
  toks <- trimws(unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE))
  toks <- toks[nzchar(toks)]
  if (isTRUE(unique)) toks <- toks[!duplicated(toks)]
  toks
}


#' @rdname lineage_helpers
#' @param v Character vector of tokens.
#' @return Single CSV string (length-1 character).
#' @keywords internal
compose_member_csv <- function(v) {
  if (length(v) == 0L) return(NA_character_)
  v <- as.character(v)
  v <- v[!is.na(v) & nzchar(v)]
  if (length(v) == 0L) return(NA_character_)
  paste(v[!duplicated(v)], collapse = ",")
}


#' @rdname lineage_helpers
#' @return \code{parse_member_part}: data.frame with columns \code{hf_id}
#'   (character, base COMID) and \code{hf_id_part} (integer, 1L when no
#'   suffix).
#' @keywords internal
parse_member_part <- function(x) {
  x <- as.character(x)
  bits <- strsplit(x, ".", fixed = TRUE)
  hf_id <- vapply(bits, function(b) if (length(b) >= 1L) b[1L] else NA_character_,
    character(1L))
  hf_id_part <- vapply(bits, function(b) {
    if (length(b) < 2L) return(1L)
    n <- suppressWarnings(as.integer(b[2L]))
    if (is.na(n)) 1L else n
  }, integer(1L))
  data.frame(hf_id = hf_id, hf_id_part = hf_id_part,
    stringsAsFactors = FALSE)
}

