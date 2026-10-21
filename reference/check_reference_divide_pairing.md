# Find reference reaches that end outside their own divide

A reach and its divide are a pair: the reach drains the divide, and its
downstream end is the divide's outlet. Where the pair is wrong the two
describe different places, and nothing downstream can reconcile them.

The splitter is where it shows. A parent divide is partitioned by
seeding a masked D8 flood at each member's downstream endpoint, and the
mask is the parent buffered by 250 m. An endpoint beyond that buffer has
no grid cell to seed from, the partition is abandoned, and the whole
parent passes to the downstream-most child while its siblings ship
owning no land. Widening the mask is not the cure: it would let every
parent's flood escape its own catchment.

Most endpoints sit marginally outside their divide because an outlet
lies on the boundary, so the measure is distance rather than
containment. Observed on the CONUS reference for VPU 04: of 150,744
paired reaches, 104,830 endpoints are within 30 m of their divide, and
147 (0.098%) are beyond 250 m, the worst 10.4 km away.

Report-only, and deliberately not wired into
[`run_pipeline()`](https://NOAA-OWP.github.io/hydrofabric/reference/run_pipeline.md).
The pairing belongs to the reference and is not ours to change, and this
asks a question about a reference release rather than about a VPU.

## Usage

``` r
check_reference_divide_pairing(
  reference_gpkg,
  flowpath_layer = "flowpaths",
  divide_layer = "divides",
  id = "flowpath_id",
  vpuid = NA_character_,
  tolerance_m = 250,
  sample = NULL,
  strict = FALSE
)
```

## Arguments

- reference_gpkg:

  Path to the reference GeoPackage.

- flowpath_layer, divide_layer:

  Layer names.

- id:

  Identifier column, present in both layers.

- vpuid:

  Optional VPU to restrict to, when the layers carry `vpuid`.

- tolerance_m:

  Distance beyond which a pairing is reported. The default is the 250 m
  buffer the splitter masks with, so it is the distance at which the
  pairing stops being usable rather than merely imprecise.

- sample:

  Number of pairs to test, or `NULL` for all.

- strict:

  If `TRUE`, warn when any are found.

## Value

Invisibly, a data frame of `id`, `dist_m` for the reaches found,
furthest first.
