# Re-merge divide parts split by geometry cleaning

`mapshaper -clean`, reached through
[`hfutils::clean_gpkg_layer()`](https://lynker-spatial.github.io/hfutils/reference/clean_gpkg_layer.html),
can bisect a narrow divide into disjoint parts. Casting that result to
`POLYGON` turns the parts into separate rows sharing one `divide_id`,
and keeping only one of them discards the rest. That silently destroys
catchment area: across a national build it cost 886 km2, with every VPU
still inside the 0.5% conservation tolerance, so nothing flagged it. It
also strands the owning flowpath outside what remains of its catchment,
which later reads as low coverage and gets the flowpath cleared and then
deleted.

## Usage

``` r
.merge_bisected_divides(d)
```

## Arguments

- d:

  Divides `sf`.

## Value

A list with `divides` (the repaired layer), `n_merged` (divides whose
parts were rejoined) and `changed` (whether anything was altered).

## Details

A bisected divide is a cleaning artifact rather than two catchments: the
parts were one divide before cleaning and belong to one flowpath after
it. So normalize to `MULTIPOLYGON`, which is lossless and already the
layer's type elsewhere, and union any rows sharing a `divide_id` back
into one multipart feature.
