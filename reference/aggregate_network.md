# Aggregate flowpath network topology without divide geometry

Flowpath-only replacement for the aggregate-to-distribution step. Runs
mainstem aggregation and headwater collapse using purely graph/tabular
operations, no `union_polygons`,
[`terra::makeValid`](https://rspatial.github.io/terra/reference/is.valid.html),
`st_snap`, or
[`terra::relate`](https://rspatial.github.io/terra/reference/relate.html).
Returns aggregated flowpaths, network table, and poi layer.

The returned `divides` element is a geometry-free stub (`data.frame`)
carrying only `divide_id` and `flowpath_toid`. Pass it with `rec_div` to
[`construct_agg_divides()`](https://NOAA-OWP.github.io/hydrofabric/reference/construct_agg_divides.md)
to build real divide geometry.

## Usage

``` r
aggregate_network(
  flowpath,
  divide,
  pois = NULL,
  vpu = NULL,
  max_junction = 2,
  ideal_size_sqkm = 10,
  min_area_sqkm = 3,
  max_area_sqkm = 15,
  min_length_km = 0.25,
  connector_ids = NULL
)
```

## Arguments

- flowpath:

  Reconciled flowpath `sf` from Stage-1 `rec_fp_agg`.

- divide:

  Reconciled divide `sf` from Stage-1 `rec_div`. Used only by
  [`prepare_network()`](https://NOAA-OWP.github.io/hydrofabric/reference/prepare_network.md)
  for initial topology setup.

- pois:

  POI `sf` from the refactor phase.

- vpu:

  VPU identifier written to network table.

- max_junction, ideal_size_sqkm, min_area_sqkm, max_area_sqkm,
  min_length_km:

  Aggregation parameters.

- connector_ids:

  Optional character vector of catchment-less connector reach ids (e.g.
  divide-less reaches over waterbodies). Such reaches are held out of
  the small-catchment collapse/absorb steps so they are preserved rather
  than merged away. `NULL` (default) treats none as connectors.

## Value

Named list: `flowpaths` (sf), `divides` (stub data.frame), `network`
(data.frame), `pois` (sf or NULL).

## Examples

``` r
if (FALSE) { # \dontrun{
s1_fp   <- sf::st_read("rec_fp_agg.gpkg", "flowpaths")
rec_div <- sf::st_read("rec_div.gpkg", "divides")
pois    <- sf::st_read("pois.gpkg", "pois")
agg <- aggregate_network(flowpath = s1_fp, divide = rec_div,
                         pois = pois, vpu = "01")
} # }
```
