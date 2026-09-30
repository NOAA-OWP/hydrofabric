#' Enforce integer typing and indexes on hydrofabric id columns
#'
#' Normalizes a GeoPackage so identifier and topology columns are stored as
#' `INTEGER` (not floating point) and the join keys are indexed. This keeps ids
#' from serializing to scientific notation, shrinks them on disk, and speeds the
#' `toid`-based joins the pipeline and downstream consumers rely on.
#'
#' The cast is **value-guarded**: a column is retyped only when every non-null
#' value is integer-valued. Prefixed string ids (`fp-`, `cat-`, `nex-`, ...) and
#' genuine real-valued measures are therefore left untouched; only columns that
#' really hold integers (reference COMIDs, hydrosequence, level paths, stream
#' order, flags, ...) are converted.
#'
#' Works entirely through SQL (no geometry is read into memory), so it is safe
#' on the largest domains where materializing geometry would exhaust RAM.
#'
#' @param gpkg Path to a GeoPackage.
#' @param layers Character vector of layers to process. `NULL` (default)
#'   processes every layer/table registered in `gpkg_contents`.
#' @param cast Logical; retype integer-valued id columns to `INTEGER`. Default
#'   `TRUE`.
#' @param index Logical; create indexes on id / toid join keys. Default `TRUE`.
#' @returns `gpkg`, invisibly.
#' @examples
#' # ids arriving as doubles are stored as INTEGER; prefixed string ids are not
#' f <- tempfile(fileext = ".gpkg")
#' con <- DBI::dbConnect(RSQLite::SQLite(), f)
#' DBI::dbExecute(con, "CREATE TABLE flowpaths (flowpath_id TEXT, hydroseq REAL)")
#' DBI::dbExecute(con, "INSERT INTO flowpaths VALUES ('fp-1', 42.0)")
#' DBI::dbExecute(con, "CREATE TABLE gpkg_contents (table_name TEXT)")
#' DBI::dbExecute(con, "INSERT INTO gpkg_contents VALUES ('flowpaths')")
#' DBI::dbDisconnect(con)
#'
#' enforce_integer_ids(f)
#'
#' con <- DBI::dbConnect(RSQLite::SQLite(), f)
#' DBI::dbGetQuery(con, "PRAGMA table_info(flowpaths)")[, c("name", "type")]
#' DBI::dbDisconnect(con)
#' @export
enforce_integer_ids <- function(gpkg, layers = NULL, cast = TRUE, index = TRUE) {
  stopifnot(is.character(gpkg), length(gpkg) == 1L, file.exists(gpkg))

  # id / topology columns that should be INTEGER when integer-valued
  int_candidates <- c(
    "flowpath_id", "flowpath_toid", "divide_id", "flowline_id", "flowline_toid",
    "nexus_id", "nexus_toid", "poi_id", "mainstemlp", "dnlevelpat", "terminalpa",
    "hydroseq", "dnhydroseq", "flowpath_hydroseq", "mainstem_id", "streamorder",
    "terminalfl", "startflag", "has_divide", "has_flowpath",
    # Nested-set upstream index. The geometry layers get these typed INTEGER when
    # add_upstream_index writes them, but `network` is rebuilt by a join that
    # promotes them to double, so they reached the product as REAL.
    "upstream_id", "num_upstreams")
  # join keys worth indexing (string or integer)
  index_candidates <- c(
    "flowpath_id", "flowpath_toid", "divide_id", "flowline_id", "flowline_toid",
    "nexus_id", "nexus_toid", "poi_id")

  con <- DBI::dbConnect(RSQLite::SQLite(), gpkg)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  contents <- tryCatch(
    DBI::dbGetQuery(con, "SELECT table_name FROM gpkg_contents")$table_name,
    error = function(e) DBI::dbListTables(con))
  layers <- if (is.null(layers)) contents else intersect(layers, contents)
  tables <- DBI::dbListTables(con)

  is_int_valued <- function(layer, col) {
    q <- DBI::dbGetQuery(con, sprintf(
      'SELECT COUNT(*) bad FROM "%s" WHERE "%s" IS NOT NULL AND "%s" <> CAST("%s" AS INTEGER)',
      layer, col, col, col))
    isTRUE(q$bad == 0)
  }

  for (layer in layers) {
    if (!layer %in% tables) next
    ti      <- DBI::dbGetQuery(con, sprintf('PRAGMA table_info("%s")', layer))
    present <- ti$name

    if (isTRUE(cast)) {
      todo <- intersect(int_candidates, present)
      todo <- todo[toupper(ti$type[match(todo, ti$name)]) != "INTEGER"]
      todo <- todo[vapply(todo, function(c) is_int_valued(layer, c), logical(1))]
      if (length(todo)) .retype_cols_integer(con, layer, todo)
    }

    if (isTRUE(index)) {
      for (ic in intersect(index_candidates, present)) {
        DBI::dbExecute(con, sprintf(
          'CREATE INDEX IF NOT EXISTS "idx_%s_%s" ON "%s"("%s")', layer, ic, layer, ic))
      }
    }
  }
  invisible(gpkg)
}

# Retype the given columns of `layer` to INTEGER in place, without touching
# geometry. GeoPackage rtree triggers on the geometry column call Spatialite
# functions (ST_IsEmpty, ...) that plain RSQLite lacks, and SQLite compiles a
# table's triggers when preparing any UPDATE on it, so the triggers are dropped
# for the duration and restored verbatim from sqlite_master afterwards.
#
# Indexes get the same treatment, and for a harder reason: SQLite refuses to
# DROP a column that an index references, so a single indexed column aborts the
# whole retype and the layer ships with every candidate still untyped. That is
# how `idx_flowpaths_upstream_id` (created by add_upstream_index) left
# hydroseq, mainstem_id, streamorder, has_divide, upstream_id and num_upstreams
# as REAL in all four merged domains. Auto-indexes carry no SQL and cannot be
# recreated, so they are excluded rather than dropped.
.retype_cols_integer <- function(con, layer, cols) {
  trg <- DBI::dbGetQuery(con, sprintf(
    "SELECT name, sql FROM sqlite_master WHERE type='trigger' AND tbl_name='%s' AND sql IS NOT NULL",
    layer))
  idx <- DBI::dbGetQuery(con, sprintf(
    "SELECT name, sql FROM sqlite_master WHERE type='index' AND tbl_name='%s' AND sql IS NOT NULL",
    layer))
  DBI::dbBegin(con)
  ok <- tryCatch({
    for (nm in trg$name) DBI::dbExecute(con, sprintf('DROP TRIGGER "%s"', nm))
    for (nm in idx$name) DBI::dbExecute(con, sprintf('DROP INDEX "%s"', nm))
    for (col in cols) {
      tmp <- paste0("__int_", col)
      DBI::dbExecute(con, sprintf('ALTER TABLE "%s" ADD COLUMN "%s" INTEGER', layer, tmp))
      DBI::dbExecute(con, sprintf('UPDATE "%s" SET "%s" = CAST("%s" AS INTEGER)', layer, tmp, col))
      DBI::dbExecute(con, sprintf('ALTER TABLE "%s" DROP COLUMN "%s"', layer, col))
      DBI::dbExecute(con, sprintf('ALTER TABLE "%s" RENAME COLUMN "%s" TO "%s"', layer, tmp, col))
    }
    for (s in idx$sql) DBI::dbExecute(con, s)
    for (s in trg$sql) DBI::dbExecute(con, s)
    DBI::dbCommit(con)
    TRUE
  }, error = function(e) {
    DBI::dbRollback(con)
    .hf_log_warn(sprintf("enforce_integer_ids: could not retype %s (%s): %s",
      layer, paste(cols, collapse = ", "), conditionMessage(e)))
    FALSE
  })
  invisible(ok)
}
