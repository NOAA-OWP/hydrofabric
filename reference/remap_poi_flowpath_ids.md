# Remap POI flowpath_ids to aggregated flowpath_ids

The aggregate-to-distribution step aggregates reconciled flowlines into
aggregated flowpaths but does not update each POI's `flowpath_id` from
the pre-aggregation reconciled flowline ID to the post-aggregation
aggregated flowpath ID. The result: the `pois` layer's `flowpath_id`
points to a reconciled ID (e.g. `6603`) that doesn't exist in the output
flowpaths layer; the POI should point to the aggregated flowpath (e.g.
`4325`) whose `poi_id` column correctly identifies it.

This helper reads the `pois` and `flowpaths` layers from `input_gpkg`
and rewrites `pois$flowpath_id` using the aggregated flowpath's `poi_id`
column as the authoritative mapping.

## Usage

``` r
remap_poi_flowpath_ids(
  input_gpkg,
  pois_layer = "pois",
  flowpath_layer = "flowpaths"
)
```

## Arguments

- input_gpkg:

  Path to the GeoPackage (must have `pois` and `flowpaths` layers).

- pois_layer:

  Name of the pois layer. Default `"pois"`.

- flowpath_layer:

  Name of the flowpaths layer. Default `"flowpaths"`.

## Value

Invisibly, integer count of remapped rows.
