# Refactor a reference hydrofabric (split, collapse, reconcile; refactor phase)

Consumes a reference fabric + FDR/FAC rasters; produces refactored
flowpaths, reconciled flowpaths + divides (raster-partitioned), and a
COMID-\>reconciled_id lineage table. This is the raster-heavy stage;
aggregate phase has no raster dependency.

## Usage

``` r
refactor_fabric(
  vpuid,
  reference_gpkg,
  fdr_path,
  fac_path,
  fallback_fdr_path = NULL,
  fallback_fac_path = NULL,
  hydrolocations_gpkg = NULL,
  exclude_hl_types = NULL,
  out_dir,
  tolerance_m = 500,
  snap_threshold_m = 250,
  split_flines_meters = 1e+09,
  collapse_flines_meters = 250,
  force_lr = TRUE,
  n_cores = NULL,
  use_checkpoints = FALSE
)
```

## Arguments

- vpuid:

  VPU identifier (e.g. `"01"`, `"21"`).

- reference_gpkg:

  Path to the reference fabric. May be multi-VPU (auto-sliced) or a
  single-VPU fabric.

- fdr_path:

  Path (or `/vsis3/...`) to the D8 flow-direction raster.

- fac_path:

  Path (or `/vsis3/...`) to the flow-accumulation raster.

- fallback_fdr_path, fallback_fac_path:

  Optional paths to fallback FDR / FAC rasters used when the primary
  rasters return NoData. `NULL` disables the fallback (default); ignored
  when equal to the primary paths.

- hydrolocations_gpkg:

  Optional path to a separate hydrolocations GeoPackage (POIs to embed).
  When `NULL`, hydrolocations are read from `reference_gpkg` if present.

- exclude_hl_types:

  Hydrolocation `hl_type` values dropped at ingest; see
  [`run_pipeline()`](https://NOAA-OWP.github.io/hydrofabric/reference/run_pipeline.md).

- out_dir:

  Parent output directory; `vpu_{vpuid}/` is created.

- tolerance_m:

  POI snap tolerance (m), used in `generate_events`.

- snap_threshold_m:

  Warn when a POI snap distance exceeds this (m).

- split_flines_meters:

  Split flowlines longer than this (m); the default (1e9) effectively
  disables Stage-1 splitting (routing-unit splitting is done later by
  [`add_flowlines()`](https://NOAA-OWP.github.io/hydrofabric/reference/add_flowlines.md)).

- collapse_flines_meters:

  Minimum inter-confluence segment length (m).

- force_lr:

  Force left/right event assignment.

- n_cores:

  Worker count for `reconcile_divides`.

- use_checkpoints:

  Logical; if `TRUE`, reuse existing stage-1 checkpoints instead of
  recomputing.

## Value

Invisibly, a list with `files` (named character, output paths),
`timings` (named list of elapsed seconds), `invariants` (named list).

## Details

Writes to `{out_dir}/vpu_{vpuid}/work/`:

- `refactored.gpkg` (single layer: flowpaths)

- `reconciled.gpkg` (layers: flowpaths, divides, events)

- `pois.gpkg` (event-bearing POIs, deduplicated)

- `hydrolocations.gpkg` (copy of the reference layer for the nextgen
  phase)

- `lineage.parquet` (COMID -\> reconciled_id, list cols coerced to CSV)

- `timings.json`, `invariants.json`

- `provenance.json` (written by the caller: `refactor_fabric()` only
  returns the record; writing happens once the caller knows the image
  digest + upstream ids).

## Examples

``` r
if (FALSE) { # \dontrun{
refactor_fabric(vpuid          = "01",
                reference_gpkg = "reference_01.gpkg",
                fdr_path       = "fdr.tif",
                fac_path       = "fac.tif",
                out_dir        = "out",
                collapse_flines_meters = 300)
} # }
```
