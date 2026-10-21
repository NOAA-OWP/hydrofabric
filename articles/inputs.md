# Pipeline inputs

A build reads the inputs below. Everything else in the product is
derived.

| Input | Contents | Required |
|----|----|----|
| Reference fabric | GeoPackage with `flowpaths` and `divides` | yes |
| Hydrolocations | GeoPackage of points with `hl_reference` | yes, or a layer in the reference |
| Grids | D8 flow-direction and flow-accumulation rasters | yes |
| Channel table | RiverML parquet (`hydraulics_path`) | for routing attributes and reach slope |
| Lakes and lake mask | lake set and reach-to-lake parquet (`lakes_path`, `lake_mask_path`) | for the `lakes` layer |

## Reference fabric

| Layer | Columns |
|----|----|
| `flowpaths` | `flowpath_id`, `flowpath_toid`, `vpuid`, `lengthkm`, `hydroseq`, `source`, `has_divide`, `divide_id`, and a level-path id (`mainstemlp`) |
| `divides` | `divide_id`, `flowpath_id`, `vpuid`, `areasqkm`, `has_flowpath` |

The level-path id becomes the product’s `mainstem_id`. A national
reference is sliced by `vpuid`, and a row whose `vpuid` names no VPU
being built is not used.

[`run_pipeline()`](https://NOAA-OWP.github.io/hydrofabric/reference/run_pipeline.md)
checks the reference before building with
[`check_reference_input()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_reference_input.md),
which catches truncated identifiers and topology that does not resolve:

``` r

check_reference_input("reference.gpkg")
```

## Hydrolocations

Hydrolocations are the gages, assimilation and routing points that
become POIs. The minimum is `POINT` geometry and an `hl_reference` such
as `nwis-07187000`.
[`prepare_hydrolocations()`](https://NOAA-OWP.github.io/hydrofabric/reference/prepare_hydrolocations.md)
snaps raw points to the reference and assigns `flowpath_id`, `poi_id`
and linear-reference measures:

``` r

library(hydrofabric)
prepare_hydrolocations(
  "gages.gpkg",
  flowpaths = "reference.gpkg",
  mapping   = c(hl_reference = "siteid"),
  out_gpkg  = "hydrolocations.gpkg"
)
```

Points more than `max_dist` (1000 m) from a flowpath are dropped with a
warning. The build excludes `hl_type` values named in `exclude_hl_types`
(default `"wbin"`) and bridge hydrolocations (`hl_type = "nbi"`).
[`check_hydrolocation_input()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_hydrolocation_input.md)
validates a set against the references before a build.

## Grids

The grids are supplied per domain as file paths, a VRT mosaic or
GeoTIFF, read through GDAL, and need not share the fabric’s CRS.

| Domain     | Grids                                        | Resolution   |
|------------|----------------------------------------------|--------------|
| SuperCONUS | NHDPlus V2.1, with HydroSHEDS v2 as fallback | 30 m         |
| Alaska     | HydroSHEDS v2                                | 1 arc-second |
| Hawaii     | NHDPlus HR                                   | 10 m         |
| PRVI       | NHDPlus HR                                   | 10 m         |

`fallback_fdr_path` and `fallback_fac_path` name an optional second
pair, used where the primary grid cannot delineate a divide. Values
taken from it are flagged in `fac_source`.

## Channel table and slope

The RiverML channel table supplies the routing attributes and, through
its `slope_m_m` column, the flowpath slope. A separate `slope_parquet`
(`id`, `slope_m_per_m`) takes precedence when given.

## Recording the inputs

[`hf_input_manifest()`](https://NOAA-OWP.github.io/hydrofabric/reference/hf_input_manifest.md)
records each input with its size and MD5, and
[`hf_verify_manifest()`](https://NOAA-OWP.github.io/hydrofabric/reference/hf_verify_manifest.md)
checks a later run against it:

``` r

hf_input_manifest(
  reference_conus = "conus_reference.gpkg",
  fdr_conus       = "conus/fdr.vrt",
  fac_conus       = "conus/fac.vrt",
  out = "inputs_manifest.json")

hf_verify_manifest("inputs_manifest.json", strict = TRUE)
```

A VRT names tiles stored beside it, so record or check the tiles too.
[`hf_output_manifest()`](https://NOAA-OWP.github.io/hydrofabric/reference/hf_output_manifest.md)
and
[`hf_verify_output_manifest()`](https://NOAA-OWP.github.io/hydrofabric/reference/hf_verify_output_manifest.md)
do the same for the products.
