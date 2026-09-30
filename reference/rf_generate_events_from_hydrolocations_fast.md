# Generate events from hydrolocations (data.table + sf nearest-feature, GEOS)

Reads only the necessary layers and uses vectorized nearest-feature
queries to snap hydrolocations onto reference flowpaths and emit split
events. This is the sole event-generation entry point used by the
refactor phase.

## Usage

``` r
rf_generate_events_from_hydrolocations_fast(
  hydrolocations_gpkg,
  flowpaths_gpkg = NULL,
  hydroloc_layer = NULL,
  flowpath_layer = "reference_flowpaths",
  out_gpkg = NULL,
  force_lr = TRUE,
  tolerance_m = 50,
  collapse_threshold_m = 250,
  min_area_m = 0,
  snap_threshold_m = 250
)
```

## Arguments

- hydrolocations_gpkg:

  Path to the GeoPackage of hydrolocations (POIs) to snap onto
  flowpaths.

- flowpaths_gpkg:

  Optional path to the reference flowpaths GeoPackage; defaults to
  `hydrolocations_gpkg` when `NULL`.

- hydroloc_layer:

  Hydrolocations layer name; auto-detected when `NULL`.

- flowpath_layer:

  Reference flowpaths layer name. Default `"reference_flowpaths"`.

- out_gpkg:

  Optional output path; when `NULL` the events `sf` is returned instead
  of written.

- force_lr:

  Logical; request left/right bank attribution where available. Default
  `TRUE`.

- tolerance_m:

  numeric (meters). Maximum distance from a flowpath for a hydrolocation
  to match; farther candidates are dropped. Default 50.

- collapse_threshold_m:

  numeric (meters). An event whose split would leave a piece shorter
  than this collapses onto the nearest flowpath endpoint. Default 250.

- min_area_m:

  numeric. Minimum drainage-area filter for events; `0` disables it
  (default).

- snap_threshold_m:

  numeric (meters). Distance within which an event snaps to a flowpath
  outlet/inlet node. Default 250.

## Value

An `sf` of events; when `out_gpkg` is supplied it is written there and
returned invisibly.
