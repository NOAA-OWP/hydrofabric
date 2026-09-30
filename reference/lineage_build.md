# Build the initial lineage table from refactor_flowlines output

Build the initial lineage table from refactor_flowlines output

## Usage

``` r
lineage_build(refactored, reconciled)
```

## Arguments

- refactored:

  An sf (or data frame) of refactored flowlines with a `flowpath_id`
  column containing base COMIDs and `.N` split suffixes.

- reconciled:

  An sf with columns `ID` (integer) and `member_flowpath_id`
  (comma-separated refactored IDs).

## Value

A data.table with one row per refactored flowline, sorted by
`reconciled_id` then `refactored_id`.

## Examples

``` r
if (FALSE) { # \dontrun{
refactored <- sf::st_read("refactored.gpkg", "refactored")
reconciled <- sf::st_read("reconciled.gpkg", "reconciled")
lin <- lineage_build(refactored, reconciled)
} # }
```
