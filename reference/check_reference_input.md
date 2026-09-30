# Validate a reference fabric before building from it

A build inherits whatever its reference asserts. Two classes of input
damage are silent at read time yet destroy the result:

## Usage

``` r
check_reference_input(
  reference_gpkg,
  layer = "flowpaths",
  id = "flowpath_id",
  toid = "flowpath_toid",
  min_resolved = 0.9,
  strict = TRUE
)
```

## Arguments

- reference_gpkg:

  Path to a reference fabric GeoPackage.

- layer:

  Flowpaths layer name. Default `"flowpaths"`.

- id, toid:

  Identifier and downstream-identifier column names.

- min_resolved:

  Minimum share of non-terminal `toid` values that must name an existing
  `id`. Terminals (`toid` 0 or `NA`) are excluded from the denominator,
  since they legitimately point nowhere.

- strict:

  If `TRUE` (default) a failure raises an error; if `FALSE` it warns and
  returns the findings.

## Value

Invisibly, a list with `ok` and one entry per check.

## Details

- **Truncated identifiers.** GeoPackage `MEDIUMINT` columns are 32-bit
  (`OFTInteger`). A 14-digit NHD-HR identifier written into one wraps
  modulo 2^32, so `flowpath_toid` becomes a negative or unrelated number
  while `flowpath_id` (64-bit `INTEGER`) survives. Nothing errors; the
  topology is simply gone.

- **Unresolvable topology.** When `flowpath_toid` names identifiers that
  are not present in `flowpath_id`, aggregation finds no connected
  reaches and degenerates to one flowpath per reach: the fabric is
  produced, passes its invariants, and is silently unaggregated.

Both are cheap to detect and expensive to discover afterwards, so this
runs before a build rather than after.

## Examples

``` r
if (FALSE) { # \dontrun{
check_reference_input("conus_reference.gpkg")
} # }
```
