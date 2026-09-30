# Extend aggregated flowpath geometry along its main-stem sub-flowlines

`aggregate_to_distribution` often returns a flowpath whose geometry only
covers the outlet portion of its catchment. Sub-flowlines attributed to
the flowpath (via `add_flowlines`) may extend upstream along eligible
branches to the divide headwater but are not part of the aggregated
flowpath line. This function rebuilds each flowpath's geometry as the
main-stem path through its sub-flowlines: from the outlet, walking
upstream at each confluence along the branch with the largest
`tot_drainage_areasqkm` until reaching a headwater.

## Usage

``` r
extend_flowpaths_to_mainstem(
  input_gpkg,
  flowpath_layer = "flowpaths",
  flowlines_layer = "flowlines",
  divides_layer = "divides"
)
```

## Arguments

- input_gpkg:

  Path to the intermediate GeoPackage (must contain `flowpaths` and
  `flowlines` layers).

- flowpath_layer:

  Flowpath layer name. Default `"flowpaths"`.

- flowlines_layer:

  Sub-flowline layer name. Default `"flowlines"`.

- divides_layer:

  Retained for backward compatibility; no longer read.

## Value

Invisibly, the updated flowpaths sf object. Writes back to `input_gpkg`
(overwrites flowpath_layer).

## Examples

``` r
if (FALSE) { # \dontrun{
extend_flowpaths_to_mainstem("out_vpu01.gpkg")
} # }
```
