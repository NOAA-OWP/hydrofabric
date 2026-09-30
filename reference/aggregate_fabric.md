# Aggregate a refactored fabric to the distribution network (aggregate phase)

Consumes stage-1 artifacts; produces the aggregated network + divides

- remapped POIs in a single `out_vpu{id}.gpkg`. No raster access.

## Usage

``` r
aggregate_fabric(
  vpuid,
  fabric_dir,
  out_dir,
  min_fragment_m = 5,
  threshold_m2 = 1000,
  gap_fill_area = 0,
  overlap_rule = "max-area",
  snap_interval = 1e-08,
  min_length_km = 0.25,
  min_area_sqkm = 3,
  ideal_size_sqkm = 10,
  max_area_sqkm = 15,
  verbose = FALSE
)
```

## Arguments

- vpuid:

  VPU identifier.

- fabric_dir:

  Path to the vpu working directory holding the refactored fabric.

- out_dir:

  Parent output directory; `vpu_{vpuid}/` is created.

- min_fragment_m:

  Sub-threshold flowpath collapse length (m), Step 3c.

- threshold_m2:

  `fill_small_holes` area threshold on divides (m^2).

- gap_fill_area:

  Gap-fill area for
  [`hfutils::clean_gpkg_layer`](https://lynker-spatial.github.io/hfutils/reference/clean_gpkg_layer.html).

- overlap_rule:

  Divide overlap resolution rule.

- snap_interval:

  mapshaper snap interval (degrees); collapses degenerate edges.

- min_length_km:

  Aggregate-to-distribution flowpath length threshold (km).

- min_area_sqkm:

  Aggregate-to-distribution divide area threshold (km^2). The minimum a
  unit must reach, applied ahead of `max_area_sqkm`.

- ideal_size_sqkm:

  Divide area the collapse aims a merged unit at (km^2).

- max_area_sqkm:

  Divide area a merged unit is not carried past (km^2), except where no
  receiver has room for a unit below `min_area_sqkm`.

- verbose:

  Propagated to internal helpers.

## Value

Invisibly, a list with `files`, `timings`.

## Details

Reads stage-1 outputs from disk, splits the `rec_div` layer into
landscape vs network on read, and exposes the stage's tunable thresholds
as documented arguments.

Writes to `{out_dir}/vpu_{vpuid}/work/`:

- `out_vpu{id}.gpkg`: layers: `flowpaths` (aggregated), `divides`,
  `network` (the aspatial topology table), `pois`, `small_segments`, and
  optional `landscape_divides`

- `timings.json`

## Examples

``` r
if (FALSE) { # \dontrun{
aggregate_fabric(vpuid      = "01",
                 fabric_dir = "out/vpu_01",
                 out_dir    = "out",
                 min_area_sqkm = 5)
} # }
```
