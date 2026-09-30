# Merge per-VPU fabrics into domain products with globally unique IDs

Assigns globally unique IDs across every VPU in `gpkg_map`, then writes
one GeoPackage per domain:

- `superconus_ngen.gpkg`: VPUs 01-18

- `ak_ngen.gpkg`: VPU 19

- `hi_ngen.gpkg`: VPU 20

- `prvi_ngen.gpkg`: VPU 21

- `guam_ngen.gpkg`: VPU 22, only when it is in `gpkg_map` (not part of
  the NWM v4 delivery)

## Usage

``` r
merge_vpus(
  gpkg_map,
  out_dir,
  mods_csv = system.file("extdata", "rpu_vpu_out.csv", package = "hydrofabric"),
  overwrite = FALSE,
  strict = TRUE,
  verbose = FALSE,
  grouping = c("vpu", "tdx_continent"),
  target_crs = "EPSG:4326"
)
```

## Arguments

- gpkg_map:

  Named list mapping vpuid -\> absolute path to `ngen_vpu*.gpkg`. All
  VPUs that should participate in the global ID pool must be present.
  Example:
  `list("01" = "/path/ngen_vpu01.gpkg", "19" = "/path/ngen_vpu19.gpkg")`

- out_dir:

  Directory where output gpkgs are written.

- mods_csv:

  Path to the cross-VPU COMID modifications CSV (CONUS only). Defaults
  to the vendored copy in `inst/extdata/rpu_vpu_out.csv`.

- overwrite:

  If `TRUE`, delete existing output files before writing.

- strict:

  If `TRUE` (default), the merge is fail-fast: unresolved cross-VPU
  topology and invariant failures abort the run.

- verbose:

  Emit per-layer progress messages.

- grouping:

  How to partition VPUs into output domains: `"vpu"` (default) groups by
  the leading integer of `vpuid` (country/VPU domain); `"tdx_continent"`
  groups by the leading digit of the TDX-Hydro HydroBASINS level-2 id
  (continent).

- target_crs:

  Output coordinate reference system for the merged domain products
  (default `"EPSG:4326"`).

## Value

Invisibly, a list with `outfiles`, `fp_map`, `poi_map`, `cross_vpu_fix`,
`invariants` (one per domain product).

## Details

IDs are assigned in contiguous blocks in VPU order, so a flowpath id
such as `fp-450003` belongs to exactly one VPU in any domain.

## Examples

``` r
if (FALSE) { # \dontrun{
gpkg_map <- list("01" = "out/vpu_01/ngen_vpu01.gpkg",
                 "19" = "out/vpu_19/ngen_vpu19.gpkg")
merge_vpus(gpkg_map, out_dir = "out/global", overwrite = TRUE)
} # }
```
