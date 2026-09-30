# Build the NextGen data model from refactored + aggregated layers (nextgen phase)

Consumes stage-1 (refactored + reconciled + hydrolocations) AND stage-2
(`out_vpu{id}.gpkg`) artifacts; produces the ngen-ready GeoPackage with
model attributes and aspatial network.

## Usage

``` r
nextgen_fabric(
  vpuid,
  fabric_dir,
  fdr_path,
  fac_path,
  fallback_fdr_path = NULL,
  fallback_fac_path = NULL,
  out_dir,
  coverage_min = 0.9,
  threshold_m2 = 1000,
  threshold_sqkm = 0.001,
  gap_fill_area = 1000,
  overlap_rule = "max-area",
  n_cores = NULL,
  hydraulics_path = NULL,
  lakes_path = NULL,
  lake_mask_path = NULL,
  slope_parquet = NULL,
  reference_gpkg = NULL,
  verbose = FALSE
)
```

## Arguments

- vpuid:

  VPU identifier.

- fabric_dir:

  Path to the vpu working directory holding the refactored + aggregated
  fabric.

- fdr_path, fac_path:

  Raster paths (local or `/vsis3/...`).

- fallback_fdr_path, fallback_fac_path:

  Optional paths to a fallback FDR/FAC pair used when the primary
  returns NA. `NULL` disables the fallback (default). Ignored when equal
  to the primary paths.

- out_dir:

  Parent output directory; `vpu_{vpuid}/` is created.

- coverage_min:

  Divide-coverage invariant + coalesce/clip floor (0-1).

- threshold_m2:

  `fill_small_holes` area threshold on ngen divides (m^2).

- threshold_sqkm:

  `prune_ghost_divides` area threshold (km^2).

- gap_fill_area:

  Gap-fill area for
  [`hfutils::clean_gpkg_layer`](https://lynker-spatial.github.io/hfutils/reference/clean_gpkg_layer.html)
  (ngen phase).

- overlap_rule:

  Divide overlap resolution rule.

- n_cores:

  Worker count for `repair_split_divides_at_pois`.

- hydraulics_path:

  Optional path to a RiverML channel table parquet, either its v2 output
  (e.g. `superconus_ml_attributes_v2.parquet`) or the earlier
  channel-geometry table. When supplied, Step 6g computes flow
  attributes for both `flowpath_id` and `flowline_id` and writes
  `flowpath-attributes` and `flowline-attributes` layers to the ngen
  GeoPackage. Reaches with no match in the table fall back through a
  geomorphic proxy. `NULL` skips Step 6g (default).

- lakes_path:

  Optional path to the NWM lakes GeoPackage (e.g. `nwm_lakes.gpkg`).
  When supplied, Step 6h identifies the outlet flowline of each lake
  (minimum `flowline_hydroseq` among COMIDs belonging to the lake),
  joins NWM lake attributes, and writes a `lakes` layer to the ngen
  GeoPackage. The reach -\> `lake_id` membership used for outlet
  fallback and member stamping comes from `lake_mask_path`. `NULL` skips
  Step 6h (default).

- lake_mask_path:

  Optional path to a lake-mask parquet: reach -\> waterbody membership
  with integer columns `comid` (reference flowline id) and `wbareacomi`
  (waterbody / `lake_id`). Required when `lakes_path` is supplied. The
  CONUS mask is materialized from NHDPlusV2 wbareacomi; oCONUS /
  superconus domains supply an equivalent with the same schema.

- slope_parquet:

  Optional path to a per-reference-id slope parquet, passed to
  [`finalize_data_model()`](https://lynker-spatial.github.io/hydrofabric/reference/finalize_data_model.md).
  `NULL` (default) takes slope from `hydraulics_path` when that table
  carries a `slope_m_m` column, and otherwise leaves `slope` `NA`.

- reference_gpkg:

  Optional path to the reference fabric GeoPackage (source of member
  slope / divide geometry). `NULL` disables features needing it.

- verbose:

  Propagated to internal helpers.

## Value

Invisibly, a list with `files`, `timings`, `invariants`.

## Details

`ngen_vpu{id}.gpkg` writes to `{out_dir}/vpu_{vpuid}/`; all other files
write to `{out_dir}/vpu_{vpuid}/work/`:

- `ngen_vpu{id}.gpkg`: layers: `flowpaths`, `divides`, `nexus`, plus
  aspatial network and model-attribute layers

- `timings.json`, `invariants.json`

## Examples

``` r
if (FALSE) { # \dontrun{
nextgen_fabric(vpuid      = "01",
               fabric_dir = "out/vpu_01",
               fdr_path   = "fdr.tif",
               fac_path   = "fac.tif",
               out_dir    = "out",
               coverage_min = 0.95)
} # }
```
