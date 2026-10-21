# Partition and Reconcile Divides with Raster-Based Splitting

Performs strict raster-based divide partitioning of reconciled divides.
Consumes GeoPackage outputs or in-memory flowpath/divide objects (the
in-memory entry point
[`partition_reconcile_divides_from_objects()`](https://NOAA-OWP.github.io/hydrofabric/reference/partition_reconcile_divides_from_objects.md)
is what
[`reconcile_divides()`](https://NOAA-OWP.github.io/hydrofabric/reference/reconcile_divides.md)
calls during the refactor phase).

For each parent divide linked to multiple final flowpaths, the function:

1.  Builds split points from downstream endpoints of member flowpaths.

2.  Runs raster-based splitting via `split_divides_d8()`.

3.  Assigns one output polygon per member flowpath.

4.  Enforces hard-fail validation on coverage, overlap, and area
    conservation.

## Usage

``` r
partition_reconcile_divides(
  input_gpkg,
  output_gpkg,
  flowline_layer = "flowline",
  divides_layer = "divides",
  events_layer = "events",
  split_events_layer = "split_events",
  lookup_table_name = "lookup_table",
  fdr_path,
  fac_path,
  write_outputs = TRUE
)
```

## Arguments

- input_gpkg:

  Character path to input GeoPackage.

- output_gpkg:

  Character path to output GeoPackage.

- flowline_layer:

  Flowline layer name. Default `"flowline"`.

- divides_layer:

  Divides layer name. Default `"divides"`.

- events_layer:

  Events layer name. Default `"events"`.

- split_events_layer:

  Split events layer name. Default `"split_events"`.

- lookup_table_name:

  Lookup table name. Default `"lookup_table"`.

- fdr_path:

  Path to D8 flow-direction raster. Required.

- fac_path:

  Path to flow-accumulation raster. Required.

- write_outputs:

  Logical; write to `output_gpkg`. Default `TRUE`.

## Value

Invisible list with updated `flowline`, `events`, `split_events`,
`divides`, and `lookup_table`.

## Examples

``` r
if (FALSE) { # \dontrun{
partition_reconcile_divides(
  input_gpkg  = "refactored.gpkg",
  output_gpkg = "rec_div.gpkg",
  fdr_path    = "fdr.tif",
  fac_path    = "fac.tif")
} # }
```
