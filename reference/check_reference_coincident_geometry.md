# Find reference reaches drawn on top of one another

Two reaches that occupy the same ground are indistinguishable to anyone
reading the fabric: one click returns both, and they appear to be one
channel counted twice. The pipeline carries such geometry through
faithfully, so the only place to catch it is the input.

The failure is not exact duplication, which is rare and easy to test
for. It is near-coincidence: reaches running within a metre of each
other for kilometres without sharing a vertex. A strict topological
overlap test sees almost nothing there, while a mouse click cannot tell
the difference. So this measures the fraction of each reach that lies
within `tolerance_m` of a different reach.

Observed on the CONUS reference at a confluence in VPU 04: four reaches
terminating at one point with the same bearing, sharing 2.3 to 5.5 km
with each other, three of them from a preliminary source and the fourth
a 522,054 km2 river drawn as a single straight 9.5 km segment.
Coincident geometry also misleads aggregation, which reads spatial
adjacency: those four tributaries were contracted into one flowpath
carrying four outlets.

Report-only, and deliberately not wired into
[`run_pipeline()`](https://lynker-spatial.github.io/hydrofabric/reference/run_pipeline.md).
The geometry belongs to the reference and is not ours to change, and
this asks a question about a reference release rather than about a VPU,
so it runs once per input rather than once per build.

## Usage

``` r
check_reference_coincident_geometry(
  reference_gpkg,
  layer = "flowpaths",
  id = "flowpath_id",
  vpuid = NA_character_,
  tolerance_m = 1,
  min_frac = 0.5,
  sample = 5000L,
  strict = FALSE
)
```

## Arguments

- reference_gpkg:

  Path to the reference GeoPackage.

- layer:

  Flowline layer name.

- id:

  Identifier column.

- vpuid:

  Optional VPU to restrict to, when the layer carries `vpuid`.

- tolerance_m:

  Distance within which two reaches count as coincident. 1 m is three
  orders of magnitude below the 30 m DEM cell, so it cannot catch
  channels that are merely close.

- min_frac:

  Fraction of a reach that must lie within `tolerance_m` of another
  before it is reported.

- sample:

  Number of reaches to test, or `NULL` for all. The default keeps the
  cost bounded on a national reference; the reported rate is an estimate
  and the message says so.

- strict:

  If `TRUE`, warn when any are found.

## Value

Invisibly, a data frame of `id`, `frac`, `length_m` for the reaches
found, most coincident first.
