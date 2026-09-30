# Prepare Hydrologic Network

Adds hydrosequence (`hydroseq`), geometric measures (area/length), and
total contributing drainage area (`tot_drainage_area`) to the
`flowpaths` element of a network list. Also standardizes geometry column
names, normalizes names to lowercase, and drops duplicates.

## Usage

``` r
prepare_network(network_list)
```

## Arguments

- network_list:

  A list containing `sf` objects:

  - `flowpaths`: flowpath/flowline features with at least `flowpath_id`,
    `flowpath_toid`.

  - `divides`: catchment/divide features (ideally with `divide_id`).

## Value

The input list with updated `flowpaths` and `divides` `sf` objects.

## Details

- `tot_drainage_area` is computed by accumulating per-flowpath
  `areasqkm` downstream; `NA` areas are treated as zero to enable
  accumulation.

- The function expects a list with `flowpaths` and `divides` `sf`
  objects.

- If duplicate rows or duplicate IDs are found, they are removed with a
  warning.

## Examples

``` r
if (FALSE) { # \dontrun{
nl <- list(flowpaths = sf::st_read("rec_fp_agg.gpkg", "flowpaths"),
           divides   = sf::st_read("rec_div.gpkg", "divides"))
nl <- prepare_network(nl)
} # }
```
