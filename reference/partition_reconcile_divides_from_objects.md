# Partition and Reconcile Divides from In-Memory Objects

Partition and Reconcile Divides from In-Memory Objects

## Usage

``` r
partition_reconcile_divides_from_objects(
  flowline,
  divides,
  pois = NULL,
  lookup_table = NULL,
  events = NULL,
  split_events = NULL,
  output_gpkg = NULL,
  fdr_path,
  fac_path,
  fallback_fdr_path = NULL,
  fallback_fac_path = NULL,
  write_outputs = TRUE,
  num_workers = NULL
)
```

## Arguments

- flowline:

  Final flowline `sf`.

- divides:

  Final divides `sf` from collapse stage.

- pois:

  sf of points of interest (POIs); their locations seed divide split
  points.

- lookup_table:

  Lookup table `data.frame`.

- events:

  Optional events `sf` (passed through).

- split_events:

  Optional split-events `sf` (passed through).

- output_gpkg:

  Output GeoPackage path when `write_outputs = TRUE`.

- fdr_path:

  Path to D8 flow-direction raster. Required.

- fac_path:

  Path to flow-accumulation raster. Required.

- fallback_fdr_path, fallback_fac_path:

  Optional paths to fallback FDR / FAC rasters, retried when a split
  returns NULL on the primary rasters (split point on NoData). `NULL`
  disables the retry (default).

- write_outputs:

  Logical; write outputs to `output_gpkg`.

- num_workers:

  Integer parallel worker count for the raster splitter. Default 4.

## Value

Invisible list with updated layers.

## Examples

``` r
if (FALSE) { # \dontrun{
partition_reconcile_divides_from_objects(
  flowline    = out$refactored,
  divides     = divides_sf,
  pois        = pois_sf,
  lookup_table = out$reconciled,
  output_gpkg = "rec_div.gpkg",
  fdr_path    = "fdr.tif",
  fac_path    = "fac.tif")
} # }
```
