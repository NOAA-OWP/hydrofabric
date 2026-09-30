# Clear disjoint flowpath-divide assignments (in-memory)

Core detection and clearing logic for flowpath-divide misassignments.
Operates on in-memory sf objects; used by
`resolve_disjoint_flowpath_divides` and by the aggregate phase before
aggregation.

## Usage

``` r
.clear_disjoint_fp_divides(fp, div, drop_orphan_cats = TRUE, flowlines = NULL)
```

## Arguments

- fp:

  sf with flowpath_id, divide_id, and line geometry.

- div:

  sf with divide_id and polygon geometry.

- drop_orphan_cats:

  If TRUE, remove divides no longer referenced.

## Value

list(fp = sf, div = sf, report = data.frame of cleared pairs)
