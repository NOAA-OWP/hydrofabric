# Coalesce fragmented flowpaths with their downstream same-mainstem neighbor

After aggregation + PRD + restore, some flowpaths end up with severely
fragmented divides (sliver-dominated MULTIPOLYGONs with tiny total area)
while their line geometry is valid. The missing territory was absorbed
by the same-mainstem downstream neighbor, so the line overshoots its
divide. This pass merges such flowpaths into their immediate downstream
neighbor when the merge is safe.

Merge rule: fp_up is absorbed into fp_down when ALL hold:

1.  fp_up.flowpath_toid == fp_down.flowpath_id (directly adjacent)

2.  fp_up.mainstem == fp_down.mainstem (same stream)

3.  No POI attached to fp_up (preserves POI=outlet semantics)

4.  fp_up's current line coverage by its own divide is below
    `coverage_threshold`

Chains cascade: if fp_up -\> fp_mid -\> fp_down and both fp_up and
fp_mid are absorbable, all three merge into fp_down.

## Usage

``` r
coalesce_fragmented_flowpaths(
  input_gpkg,
  flowpath_layer = "flowpaths",
  divides_layer = "divides",
  flowlines_layer = "flowlines",
  pois_layer = "pois",
  network_layer = "network",
  coverage_threshold = 0.9
)
```

## Arguments

- input_gpkg:

  Path to the intermediate GeoPackage.

- flowpath_layer, divides_layer, flowlines_layer, pois_layer,
  network_layer:

  Layer names. Defaults match the pre-ngen pipeline.

- coverage_threshold:

  fp_up coverage below this triggers merge candidacy (default 0.9).

## Value

Invisibly, a data.frame of `{up, down}` merge pairs.
