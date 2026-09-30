# Split absorbed chains that receive routed flow into their own flowpaths

Split absorbed chains that receive routed flow into their own flowpaths

## Usage

``` r
.split_junction_absorbed(gpkg, s1_fp, vpuid = NA_character_)
```

## Arguments

- gpkg:

  Stage-2 GeoPackage holding `flowpaths` and `network`, before
  `construct_agg_divides` has run.

- s1_fp:

  Reconciled reaches (`sf`) carrying `flowpath_id`, `flowpath_toid`,
  `levelpathid` and geometry.

- vpuid:

  VPU label for messages.

## Value

Invisibly, the number of chains split out.
