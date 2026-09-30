# Attach aggregated flowpath IDs to a lineage table

Attach aggregated flowpath IDs to a lineage table

## Usage

``` r
lineage_attach_aggregated(lineage, network)
```

## Arguments

- lineage:

  A lineage `data.table` from
  [`lineage_build()`](https://lynker-spatial.github.io/hydrofabric/reference/lineage_build.md).

- network:

  A data.frame or sf object with columns `flowpath_id` (aggregated,
  integer or character) and `rec_id_list` (comma-separated reconciled
  IDs): typically `sf::st_read(out_gpkg, "network")`.

## Value

The lineage table with `aggregated_id` populated (in place).

## Examples

``` r
if (FALSE) { # \dontrun{
lin     <- lineage_build(refactored, reconciled)
network <- sf::st_read("out_vpu01.gpkg", "network")
lin <- lineage_attach_aggregated(lin, network)
} # }
```
