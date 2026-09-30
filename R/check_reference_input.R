#' Validate a reference fabric before building from it
#'
#' A build inherits whatever its reference asserts. Two classes of input damage
#' are silent at read time yet destroy the result:
#'
#' * **Truncated identifiers.** GeoPackage `MEDIUMINT` columns are 32-bit
#'   (`OFTInteger`). A 14-digit NHD-HR identifier written into one wraps modulo
#'   2^32, so `flowpath_toid` becomes a negative or unrelated number while
#'   `flowpath_id` (64-bit `INTEGER`) survives. Nothing errors; the topology is
#'   simply gone.
#' * **Unresolvable topology.** When `flowpath_toid` names identifiers that are
#'   not present in `flowpath_id`, aggregation finds no connected reaches and
#'   degenerates to one flowpath per reach: the fabric is produced, passes its
#'   invariants, and is silently unaggregated.
#'
#' Both are cheap to detect and expensive to discover afterwards, so this runs
#' before a build rather than after.
#'
#' @param reference_gpkg Path to a reference fabric GeoPackage.
#' @param layer Flowpaths layer name. Default `"flowpaths"`.
#' @param id,toid Identifier and downstream-identifier column names.
#' @param min_resolved Minimum share of non-terminal `toid` values that must
#'   name an existing `id`. Terminals (`toid` 0 or `NA`) are excluded from the
#'   denominator, since they legitimately point nowhere.
#' @param strict If `TRUE` (default) a failure raises an error; if `FALSE` it
#'   warns and returns the findings.
#' @returns Invisibly, a list with `ok` and one entry per check.
#' @examples
#' \dontrun{
#' check_reference_input("conus_reference.gpkg")
#' }
#' @export
check_reference_input <- function(reference_gpkg,
                                  layer = "flowpaths",
                                  id    = "flowpath_id",
                                  toid  = "flowpath_toid",
                                  min_resolved = 0.90,
                                  strict = TRUE) {
  stopifnot(is.character(reference_gpkg), length(reference_gpkg) == 1L,
    file.exists(reference_gpkg))

  con <- DBI::dbConnect(RSQLite::SQLite(), reference_gpkg)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  tabs <- DBI::dbListTables(con)
  if (!layer %in% tabs)
    stop(sprintf("layer '%s' not in %s", layer, basename(reference_gpkg)),
      call. = FALSE)
  cols <- DBI::dbGetQuery(con, sprintf('PRAGMA table_info("%s")', layer))
  for (cl in c(id, toid))
    if (!cl %in% cols$name)
      stop(sprintf("column '%s' not in layer '%s'", cl, layer), call. = FALSE)

  INT32_MAX <- 2147483647

  q <- DBI::dbGetQuery(con, sprintf('
    SELECT COUNT(*) n,
           SUM(CASE WHEN "%s" IS NULL OR "%s" = 0 THEN 1 ELSE 0 END) n_terminal,
           SUM(CASE WHEN "%s" < 0 THEN 1 ELSE 0 END) n_negative,
           MIN("%s") min_toid, MAX("%s") max_toid, MAX("%s") max_id
    FROM "%s"', toid, toid, toid, toid, toid, id, layer))

  n_res <- DBI::dbGetQuery(con, sprintf('
    SELECT COUNT(*) n FROM "%s"
    WHERE "%s" IS NOT NULL AND "%s" <> 0
      AND "%s" IN (SELECT "%s" FROM "%s")',
    layer, toid, toid, toid, id, layer))$n

  n_nonterm <- q$n - q$n_terminal
  share     <- if (n_nonterm > 0L) n_res / n_nonterm else NA_real_
  declared  <- toupper(cols$type[match(c(id, toid), cols$name)])
  names(declared) <- c(id, toid)

  # Two severities. `fatal` checks describe damage that is already present and
  # would silently produce a wrong fabric, so they stop the build. `risk` checks
  # describe a latent hazard in an input that is currently intact: reporting
  # them must not block a valid build.
  checks <- list()
  .chk <- function(ok, msg, severity = "fatal")
    list(ok = isTRUE(ok), msg = msg, severity = severity)

  # A negative toid cannot arise from real data: it is a wrapped 64-bit value.
  checks$toid_not_truncated <- .chk(
    q$n_negative == 0L,
    sprintf("%d %s value(s) are negative (32-bit wrap of a 64-bit id); min = %s",
      q$n_negative, toid, format(q$min_toid, scientific = FALSE)))

  # A 32-bit column holding ids near or past the limit is corruption waiting to
  # happen even when the current values survived.
  narrow <- declared[declared %in% c("MEDIUMINT", "SMALLINT", "TINYINT", "INT")]
  at_risk <- length(narrow) > 0L && q$max_id > INT32_MAX * 0.9
  checks$id_width_safe <- .chk(
    !at_risk,
    if (at_risk)
      sprintf(paste0("max %s = %s exceeds 90%% of the 32-bit limit while %s ",
        "declared narrow (%s); widen to INTEGER before ids grow"),
        id, format(q$max_id, scientific = FALSE),
        paste(names(narrow), collapse = "/"), paste(narrow, collapse = "/"))
    else sprintf("max %s = %s fits its declared type(s) (%s)",
      id, format(q$max_id, scientific = FALSE),
      paste(sprintf("%s=%s", names(declared), declared), collapse = " ")),
    severity = "risk")

  checks$toid_resolves <- .chk(
    !is.na(share) && share >= min_resolved,
    sprintf("%d of %d non-terminal %s resolve to an existing %s (%.1f%%, floor %.0f%%)",
      n_res, n_nonterm, toid, id, 100 * share, 100 * min_resolved))

  fatal <- names(Filter(function(c) !c$ok && c$severity == "fatal", checks))
  risk  <- names(Filter(function(c) !c$ok && c$severity == "risk",  checks))
  ok    <- length(fatal) == 0L

  for (nm in names(checks))
    message(sprintf("[reference:%s] %s %-20s %s", basename(reference_gpkg),
      if (checks[[nm]]$ok) "ok" else if (checks[[nm]]$severity == "risk") "!" else "x",
      nm, checks[[nm]]$msg))

  if (length(risk))
    warning(sprintf("reference input at risk (%s): %s",
      basename(reference_gpkg), paste(risk, collapse = ", ")), call. = FALSE)
  if (!ok) {
    msg <- sprintf("reference input unusable (%s): %s",
      basename(reference_gpkg), paste(fatal, collapse = ", "))
    if (isTRUE(strict)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  invisible(c(list(ok = ok), checks))
}
