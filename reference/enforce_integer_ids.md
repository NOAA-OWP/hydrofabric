# Enforce integer typing and indexes on hydrofabric id columns

Normalizes a GeoPackage so identifier and topology columns are stored as
`INTEGER` (not floating point) and the join keys are indexed. This keeps
ids from serializing to scientific notation, shrinks them on disk, and
speeds the `toid`-based joins the pipeline and downstream consumers rely
on.

## Usage

``` r
enforce_integer_ids(gpkg, layers = NULL, cast = TRUE, index = TRUE)
```

## Arguments

- gpkg:

  Path to a GeoPackage.

- layers:

  Character vector of layers to process. `NULL` (default) processes
  every layer/table registered in `gpkg_contents`.

- cast:

  Logical; retype integer-valued id columns to `INTEGER`. Default
  `TRUE`.

- index:

  Logical; create indexes on id / toid join keys. Default `TRUE`.

## Value

`gpkg`, invisibly.

## Details

The cast is **value-guarded**: a column is retyped only when every
non-null value is integer-valued. Prefixed string ids (`fp-`, `cat-`,
`nex-`, ...) and genuine real-valued measures are therefore left
untouched; only columns that really hold integers (reference COMIDs,
hydrosequence, level paths, stream order, flags, ...) are converted.

Works entirely through SQL (no geometry is read into memory), so it is
safe on the largest domains where materializing geometry would exhaust
RAM.

## Examples

``` r
# ids arriving as doubles are stored as INTEGER; prefixed string ids are not
f <- tempfile(fileext = ".gpkg")
con <- DBI::dbConnect(RSQLite::SQLite(), f)
DBI::dbExecute(con, "CREATE TABLE flowpaths (flowpath_id TEXT, hydroseq REAL)")
#> [1] 0
DBI::dbExecute(con, "INSERT INTO flowpaths VALUES ('fp-1', 42.0)")
#> [1] 1
DBI::dbExecute(con, "CREATE TABLE gpkg_contents (table_name TEXT)")
#> [1] 0
DBI::dbExecute(con, "INSERT INTO gpkg_contents VALUES ('flowpaths')")
#> [1] 1
DBI::dbDisconnect(con)

enforce_integer_ids(f)

con <- DBI::dbConnect(RSQLite::SQLite(), f)
DBI::dbGetQuery(con, "PRAGMA table_info(flowpaths)")[, c("name", "type")]
#>          name    type
#> 1 flowpath_id    TEXT
#> 2    hydroseq INTEGER
DBI::dbDisconnect(con)
```
