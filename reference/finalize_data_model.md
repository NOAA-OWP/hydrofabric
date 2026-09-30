# Finalize NextGen output to the Lynker hydrofabric data model spec

Applies the final column renames, adds VPU/streamorder/slope fields,
rolls up hydrolocation attributes on POIs, and promotes a dedicated
`hydrolocations` layer. Mutates the provided ngen GeoPackage in place.

See `docs/data_model.md` for the full mapping against the spec.

## Usage

``` r
finalize_data_model(
  ngen_gpkg,
  reference_gpkg = NULL,
  vpuid,
  events_layer = "pois",
  events_gpkg = NULL,
  slope_parquet = NULL
)
```

## Arguments

- ngen_gpkg:

  Path to the NextGen GeoPackage produced by `nextgen_from_gpkg` +
  cleanup steps.

- reference_gpkg:

  Path to the reference fabric (used to source `slope` per member
  reference COMID). If `NULL`, slope stays `NA`.

- vpuid:

  Single VPU identifier to stamp on every layer (character).

- events_layer:

  Source layer for hydrolocation attributes. Default `"pois"`, the
  pre-ngen pois layer carries the full hl\_\* set (hl_class, hl_type,
  hl_source, hl_uri, etc). The `events` layer in the intermediate gpkg
  only has a subset, so prefer `pois`.

- events_gpkg:

  Optional: if the hydrolocation source layer isn't in `ngen_gpkg`, read
  it from this gpkg (typically the pre-ngen intermediate).

- slope_parquet:

  Optional path to a parquet of per-reference-id slope used to populate
  `slope`: either a slope table (`id`, `slope_m_per_m`) or a channel
  table (`flowpath_id`, `slope_m_m`). `NULL` leaves `slope` `NA`
  (default).

## Value

Invisibly, the input `ngen_gpkg` path.

## Examples

``` r
if (FALSE) { # \dontrun{
finalize_data_model(ngen_gpkg      = "ngen_vpu01.gpkg",
                    reference_gpkg = "reference_01.gpkg",
                    vpuid          = "01",
                    events_gpkg    = "out_vpu01.gpkg")
} # }
```
