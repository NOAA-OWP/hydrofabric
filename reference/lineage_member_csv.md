# Recompute `member_comid` CSV strings from a lineage table

Utility to validate that the lineage table carries the same information
as the existing `member_comid` strings.

## Usage

``` r
lineage_member_csv(lineage, group_by = "aggregated_id")
```

## Arguments

- lineage:

  A lineage data.table.

- group_by:

  Column to group by. Default `"aggregated_id"` (NextGen /
  aggregated-flowpath member lists).

## Value

data.table with columns `group_by` and `member_comid` (CSV).

## Examples

``` r
lin <- data.table::data.table(
  refactored_id = c("101", "102", "103"),
  aggregated_id = c("fp-1", "fp-1", "fp-2"))
lineage_member_csv(lin, group_by = "aggregated_id")
#>    aggregated_id member_comid
#>           <char>       <char>
#> 1:          fp-1      101,102
#> 2:          fp-2          103
```
