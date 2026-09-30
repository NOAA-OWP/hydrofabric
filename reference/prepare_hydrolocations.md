# Prepare an externally-supplied hydrolocation set for embedding

Maps a raw hydrolocation point dataset (e.g. the gage / assimilation
points supplied by OWP) to the canonical `hydrolocations` schema that
[`run_pipeline()`](https://lynker-spatial.github.io/hydrofabric/reference/run_pipeline.md)
embeds, and optionally writes it as a GeoPackage suitable for the
`hydrolocations_gpkg` argument.

## Usage

``` r
prepare_hydrolocations(
  x,
  flowpaths,
  mapping = character(0),
  max_dist = 1000,
  max_link_m = 250,
  link_keep_da_ratio = 100,
  link_keep_min_da_sqkm = 100,
  out_gpkg = NULL,
  layer = NULL
)
```

## Arguments

- x:

  An `sf` of hydrolocation POINTs, or a path to a point GeoPackage.

- flowpaths:

  Reference flowpaths to snap onto: an `sf` or a path to the reference
  GeoPackage (its `flowpaths` layer is read). Required.

- mapping:

  Named character vector mapping *canonical* column -\> *source* column,
  e.g. `c(hl_reference = "gage_id", hl_type = "type")`. Canonical
  columns: `hl_reference` (REQUIRED), `hl_type`, `hl_link`, `hl_class`,
  `hl_source`. Anything unmapped/absent is filled `NA`.

- max_dist:

  Max snap distance (m); POIs farther than this are dropped.

- max_link_m:

  Max distance a point may be from the reach its source `flowpath_id`
  prescribes before that link is distrusted and the nearest flowpath
  used instead. Reported when it triggers. Default 1000.

- link_keep_da_ratio:

  A distrusted link is restored when the nearest flowpath drains less
  than this fraction of the reach the source names (default 100, i.e.
  under 1%). Guards against relocating a large-river location onto an
  adjacent minor feature.

- link_keep_min_da_sqkm:

  Only apply that rescue when the named reach drains at least this much,
  so the comparison is meaningful. Default 100.

- out_gpkg:

  Optional output path. When given, writes the `hydrolocations` layer
  there and returns the path; otherwise returns the prepared `sf`.

- layer:

  Layer to read when `x` is a GeoPackage path (default: first).

## Value

The prepared hydrolocations `sf`, or `out_gpkg` when written.

## Details

The POIs are **snapped onto `flowpaths`** and assigned `flowpath_id`,
`poi_id` (one per flowpath, co-located POIs bundled), and
linear-reference `pathmeas`/`reachmeas`. This pre-assignment is
required: the refactor phase's spatial re-snap fallback is lossy (it
drops many POIs when the input lacks `flowpath_id`/`poi_id`), so raw
points MUST be snapped here first. Rows without `hl_reference`, or
farther than `max_dist` from any flowpath, are dropped with a warning.

## Examples

``` r
if (FALSE) { # \dontrun{
# OWP delivers gage points with columns `siteid` and `network`:
prepare_hydrolocations(
  "owp_gages.gpkg", flowpaths = "reference.gpkg",
  mapping  = c(hl_reference = "siteid", hl_source = "network", hl_class = "gage"),
  out_gpkg = "hydrolocations.gpkg")
# -> run_pipeline(reference_gpkg = ref, hydrolocations_gpkg = "hydrolocations.gpkg", ...)
} # }
```
