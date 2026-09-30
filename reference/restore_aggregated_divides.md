# Restore aggregated divide areas dropped by PRD / aggregation

Two kinds of area-loss can affect an aggregated divide between the
reference NHD input and the final output:

1.  **Aggregation-stage loss**: the aggregate-to-distribution step
    internally runs
    [`hfutils::clean_geometry`](https://lynker-spatial.github.io/hfutils/reference/clean_geometry.html)
    on the aggregated divides, whose small-part glue logic can crush a
    reconciled-split-sibling polygon (e.g. the `.2` half of an reference
    COMID) to a ~100 m2 residual. Detected by comparing each aggregated
    cat's area to the **sum of its member reconciled-divide areas**
    (from `rec_div`).

2.  **PRD-stage loss**: the refactor-split may drop one of the split
    children entirely (typically a narrow sliver that falls below a
    cleanup threshold inside `reconcile_divides`), so its area never
    appears in `rec_div` at all. The sibling's divide survives; the lost
    half vanishes into no catchment. Detected only when
    `reference_divides` is supplied, by comparing the aggregated cat's
    area to the **sum of reference NHD divide areas for its member
    COMIDs**.

Whichever reference is larger (rec_div-based or reference-based) wins.
Shortfalls beyond `threshold` trigger a rebuild from the union of:

- rec_div geometries for members whose reconciled divide exists, plus

- reference-divide geometries for members missing from rec_div (the PRD
  dropouts). Divide polygons are not otherwise modified.

## Usage

``` r
restore_aggregated_divides(
  input_gpkg,
  rec_div,
  reference_divides = NULL,
  threshold = 0.9,
  abs_shortfall_m2 = 10000,
  split_gap_m2 = 5e+05,
  divide_layer = "divides",
  network_layer = "network",
  verbose = FALSE
)
```

## Arguments

- input_gpkg:

  Path to the intermediate GeoPackage (contains `divides` and `network`
  layers).

- rec_div:

  sf object of reconciled divides returned by
  [`reconcile_divides()`](https://lynker-spatial.github.io/hydrofabric/reference/reconcile_divides.md)
  (before aggregation). Must contain `flowpath_id` (= reconciled integer
  ID) and geometry.

- reference_divides:

  Optional sf of reference NHD divides (the pipeline's pre-refactor
  input). Must contain `divide_id` (the reference COMID as character)
  and geometry. When supplied, enables detection of PRD-stage drops
  where a split child was never produced.

- threshold:

  Fraction. Aggregated divides with area \<
  `threshold * max( expected_rec, expected_ref)` are restored. Default
  `0.95`.

- abs_shortfall_m2:

  numeric (m^2). A divide is also restored when its absolute area
  shortfall vs expected exceeds this, catching small-percentage gaps on
  large divides. Default `1e4`.

- split_gap_m2:

  numeric (m^2). Minimum reference-divide gap area that qualifies as a
  genuine PRD split dropout worth filling. Default `5e5`.

- divide_layer:

  Name of the divides layer. Default `"divides"`.

- network_layer:

  Name of the network layer. Default `"network"`.

- verbose:

  If TRUE, emit a per-divide breakdown of area shortfalls and recovery
  deltas, plus a source-breakdown (rec_div / ref / both) per fix.

## Value

Invisibly, character vector of restored aggregated divide ids.
