# Detect and clear flowpath-divide misassignments

Finds flowpaths whose line geometry is spatially disjoint from the
divide they're assigned to (coverage = 0, distance above
`distance_threshold_m`) and clears the bad assignment. The flowpath
remains in the network (topology is preserved via attribute IDs) but
loses its catchment.

Root cause is an aggregation join artifact that appears at VPU scale:
reconciled flowlines with missing/default levelpathid (=1) can be
grouped into aggregate flowpaths whose member subs are geographically
unrelated, leading to a bogus cat pairing. The proper fix is to call
this check before aggregation (aggregate step 3c) so the bad pair never
enters aggregate_to_distribution.

## Usage

``` r
resolve_disjoint_flowpath_divides(
  input_gpkg,
  flowpath_layer = "flowpaths",
  divides_layer = "divides",
  distance_threshold_m = 100,
  drop_orphan_cats = TRUE
)
```

## Arguments

- input_gpkg:

  Path to the GeoPackage (mutated in place).

- flowpath_layer, divides_layer:

  Layer names.

- distance_threshold_m:

  Kept for back-compat; no longer used as a gate (coverage = 0 is the
  sole criterion). Reported in output for triage.

- drop_orphan_cats:

  If TRUE, remove the now-orphan divide polygons that no flowpath
  references. Default TRUE.

## Value

Invisibly, a data.frame of cleared
`{flowpath_id, old_divide_id, distance_m}` rows.
