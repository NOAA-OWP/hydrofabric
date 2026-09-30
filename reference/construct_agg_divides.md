# Construct aggregated divides from Stage-1 rec_div via rec_id_list

Given the network table from
[`aggregate_network()`](https://lynker-spatial.github.io/hydrofabric/reference/aggregate_network.md)
(which carries `flowpath_id` -\> `rec_id_list`) and the Stage-1
`rec_div` polygons, unions the `rec_div` members for each aggregated
flowpath.

No raster, no snapping, no mapshaper. Stage-3 `fill_small_holes` and
`clean_gpkg_layer` handle any residual topology.

## Usage

``` r
construct_agg_divides(network, rec_div, flowpaths)
```

## Arguments

- network:

  Network `data.frame` from
  [`aggregate_network()`](https://lynker-spatial.github.io/hydrofabric/reference/aggregate_network.md)
  with `flowpath_id` and `rec_id_list`.

- rec_div:

  Stage-1 reconciled divide `sf` with `flowpath_id` (= rec_id).

- flowpaths:

  Aggregated flowpath `sf` from
  [`aggregate_network()`](https://lynker-spatial.github.io/hydrofabric/reference/aggregate_network.md).

## Value

`sf` with `divide_id`, `flowpath_toid`, `areasqkm`, `has_flowline`,
`ds_id`, `type`, geometry (POLYGON).

## Examples

``` r
if (FALSE) { # \dontrun{
agg     <- aggregate_network(flowpath = s1_fp, divide = rec_div, pois = pois)
rec_div <- sf::st_read("rec_div.gpkg", "divides")
divides <- construct_agg_divides(network   = agg$network,
                                 rec_div   = rec_div,
                                 flowpaths = agg$flowpaths)
} # }
```
