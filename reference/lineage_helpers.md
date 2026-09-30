# Hydrofabric lineage CSV utilities

Throughout the refactor and aggregate phases several columns store
multi-valued lineage as a comma-separated string: `member_comid` (CSV of
original NHDPlus COMIDs, possibly carrying a `.N` split suffix per
child), and the (mis-named) `flowline_id` column on Stage-2 `network` /
`flowpaths` (CSV of Stage-1 `rec_id`s). These helpers provide one
canonical parser/composer so the dozens of inline
[`strsplit()`](https://rdrr.io/r/base/strsplit.html) +
[`trimws()`](https://rdrr.io/r/base/trimws.html)

- [`nzchar()`](https://rdrr.io/r/base/nchar.html) sites collapse to a
  single implementation.

## Usage

``` r
expand_member_csv(x, unique = FALSE)

compose_member_csv(v)

parse_member_part(x)
```

## Arguments

- x:

  A character vector of CSV strings, or a single CSV string.

- unique:

  If TRUE, deduplicate the expanded vector (input order preserved).

- v:

  Character vector of tokens.

## Value

Character vector of tokens.

Single CSV string (length-1 character).

`parse_member_part`: data.frame with columns `hf_id` (character, base
COMID) and `hf_id_part` (integer, 1L when no suffix).

## Details

- `expand_member_csv(x)`: vector or single string of CSV -\> flat
  character vector with whitespace stripped and empty tokens dropped.
  Use `unique=TRUE` to deduplicate.

- `compose_member_csv(v)`: character vector -\> single CSV string with
  deterministic ordering (input order preserved, duplicates dropped on
  first occurrence).

- `parse_member_part(x)`: split a COMID lineage token of the form
  `"COMID"` or `"COMID.N"` into base `hf_id` (numeric COMID as
  character) and `hf_id_part` (integer N, defaulting to 1L when absent).
  Vectorised.

These functions are intentionally pure: no I/O, no spatial side effects.
