# Run the end-to-end hydrofabric pipeline for a single VPU

Local convenience wrapper that runs the refactor, aggregate and nextgen
phases sequentially on the same filesystem. Each stage is also exported
and can be called directly.

## Usage

``` r
run_pipeline(
  vpuid,
  reference_gpkg,
  fdr_path,
  fac_path,
  fallback_fdr_path = NULL,
  fallback_fac_path = NULL,
  hydrolocations_gpkg = NULL,
  exclude_hl_types = "wbin",
  out_dir = "out",
  from = "refactor",
  through = "nextgen",
  stages = NULL,
  n_cores = NULL,
  skip_if_exists = FALSE,
  verbose = FALSE,
  hydraulics_path = NULL,
  lakes_path = NULL,
  lake_mask_path = NULL,
  slope_parquet = NULL,
  use_checkpoints = FALSE,
  check_input = TRUE,
  ...
)
```

## Arguments

- vpuid:

  VPU identifier (character, e.g. `"01"`, `"03N"`, `"21"`).

- reference_gpkg:

  Path to the reference fabric (multi-VPU or single-VPU; extraction is
  auto-detected).

- fdr_path, fac_path:

  Paths to the FDR / FAC rasters.

- fallback_fdr_path, fallback_fac_path:

  Optional paths to a fallback FDR/FAC pair. `NULL` disables fallback
  (default).

- hydrolocations_gpkg:

  Optional path to a separate GeoPackage holding the `hydrolocations`
  (POIs) + `pois` layers. `NULL` (default) reads them from
  `reference_gpkg` (all-in-one). When supplied, the reference need only
  carry `flowpaths` + `divides`, and POIs are sourced from this file.

- exclude_hl_types:

  Hydrolocation `hl_type` values to drop before the build sees them. The
  default removes waterbody INLETS (`wbin`), 33,926 of the 69,036
  supplied locations. An inlet names no computational boundary a
  consumer asks for, and each one costs twice: it pins its own unit
  against merging and breaks the levelpath, stranding the unit below it
  too, which is why they dominate the sub-minimum population.

  Waterbody OUTLETS (`wbout`) are kept. An outlet is a real routing
  boundary , it is where a lake discharges, so it earns its break in a
  way an inlet does not. Gages (`nwis`, `rfc`) and terminal/coastal
  outlets (the `owp_*` types) are kept for the same reason. Lake
  handling is independent of this either way: lakes come from
  `lakes_path`, a separate input.

  `NULL` or `character(0)` keeps every type. Filtering is by row, so a
  location carrying a wanted type as well survives on that type.

- out_dir:

  Directory for outputs. A `vpu_{vpuid}/` subfolder is created; the
  stage outputs all land there so this call is backwards-compatible with
  the monolithic pre-split runs.

- from, through:

  Stage range to run, by name (`"refactor"`, `"aggregate"`, `"nextgen"`,
  in order). Defaults run the full pipeline. Use `through = "aggregate"`
  to stop before the ngen output, or `from = "aggregate"` to resume from
  an existing refactored fabric.

- stages:

  Optional integer stage numbers (e.g. `c(2, 3)`) as an alias for
  `from`/`through`: when supplied, `from`/`through` are set to the
  min/max stage. `NULL` (default) uses `from`/`through` directly.

- n_cores:

  Worker count for parallel steps. Default `detectCores() - 2`.

- skip_if_exists:

  If `TRUE`, skip a VPU whose target output (the ngen gpkg when
  `through = "nextgen"`, otherwise the aggregate output gpkg) already
  exists. Default `FALSE` (always rebuild). This is a plain
  already-built check, not content-hash reuse.

- verbose:

  Pass through to individual steps for extra logging.

- hydraulics_path:

  Optional path to a RiverML channel table parquet, either its v2 output
  (`tw_bf_m`, `y_bf_m`, `r_bf`, roughness, `slope_m_m`) or the earlier
  channel-geometry table. When supplied, routing attributes
  (`flowpath-attributes`, `flowline-attributes`) are computed in Step
  6g, and a `slope_m_m` column also supplies flowpath slope when
  `slope_parquet` is `NULL`. `NULL` skips Step 6g (default).

- lakes_path:

  Optional path to the NWM lakes GeoPackage; forwarded to
  [`nextgen_fabric()`](https://NOAA-OWP.github.io/hydrofabric/reference/nextgen_fabric.md)
  (Step 6h). `NULL` skips lakes (default).

- lake_mask_path:

  Optional path to the lake-mask parquet (reach -\> waterbody
  membership); forwarded to
  [`nextgen_fabric()`](https://NOAA-OWP.github.io/hydrofabric/reference/nextgen_fabric.md).
  Required when `lakes_path` is supplied.

- slope_parquet:

  Optional path to a per-reference-id slope parquet; forwarded to
  [`nextgen_fabric()`](https://NOAA-OWP.github.io/hydrofabric/reference/nextgen_fabric.md).
  `NULL` (default) takes slope from `hydraulics_path` when that table
  carries a `slope_m_m` column, and otherwise leaves slope `NA`.

- use_checkpoints:

  If `TRUE`, allow stage checkpoint reuse where supported (currently the
  refactor phase). Default `FALSE` recomputes from scratch.

- check_input:

  If `TRUE` (default), validate the inputs before building:
  [`check_reference_input()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_reference_input.md)
  aborts on damaged identifiers or unresolvable topology, and
  [`check_grid_coverage()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_grid_coverage.md)
  aborts when a primary grid does not cover the domain, which is what a
  configuration inheriting another domain's grid looks like.

- ...:

  Stage tuning parameters, forwarded to whichever stage declares each
  (e.g. `min_area_sqkm = 5`, `coverage_min = 0.95`). See
  [`refactor_fabric()`](https://NOAA-OWP.github.io/hydrofabric/reference/refactor_fabric.md),
  [`aggregate_fabric()`](https://NOAA-OWP.github.io/hydrofabric/reference/aggregate_fabric.md),
  [`nextgen_fabric()`](https://NOAA-OWP.github.io/hydrofabric/reference/nextgen_fabric.md)
  for the full set and defaults.

## Value

Invisibly, a list with `vpuid`, `status`, `ngen_gpkg`, `out_gpkg`,
`timings`, `invariants`, `skipped`.

## Examples

``` r
if (FALSE) { # \dontrun{
run_pipeline(vpuid          = "01",
             reference_gpkg = "reference_01.gpkg",
             fdr_path       = "fdr.tif",
             fac_path       = "fac.tif",
             out_dir        = "out",
             through        = "nextgen")
} # }
```
