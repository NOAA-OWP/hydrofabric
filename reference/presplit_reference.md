# Pre-split a multi-VPU reference fabric into per-VPU slices

Reads each layer of `reference_gpkg` **once**, groups rows by `vpuid`,
and writes one `ref_vpu{id}.gpkg` per VPU into
`<out_dir>/vpu_{id}/work/`. When these files exist, refactor Step 0 is
an instant skip: eliminating the per-VPU random read against the large
source file.

## Usage

``` r
presplit_reference(
  reference_gpkg,
  out_dir,
  vpuids = NULL,
  layers = c("flowpaths", "divides", "hydrolocations", "pois"),
  overwrite = FALSE
)
```

## Arguments

- reference_gpkg:

  Path to the multi-VPU reference GeoPackage.

- out_dir:

  Parent output directory (same as `run_pipeline`).

- vpuids:

  Character vector of VPU IDs to split out. `NULL` (default) splits all
  VPUs found in the `flowpaths` layer.

- layers:

  Layers to extract. Default matches refactor Step 0.

- overwrite:

  If `TRUE`, overwrite existing per-VPU slices. Default `FALSE`.

## Value

Invisibly, a named character vector `vpuid -> path` for every file
written.

## Examples

``` r
if (FALSE) { # \dontrun{
presplit_reference(reference_gpkg = "reference_conus.gpkg",
                   out_dir        = "out",
                   vpuids         = c("01", "02"))
} # }
```
