# Dissolve multi-part divide polygons to their largest component

Lightweight replacement for
[`hfutils::clean_geometry`](https://lynker-spatial.github.io/hfutils/reference/clean_geometry.html)'s
intra-divide fragmentation handling. For each divide whose geometry is a
MULTIPOLYGON with disjoint components, keeps the LARGEST component as
the sole geometry for that ID and drops the smaller components.

Rationale: `clean_geometry`'s "glue small orphan parts to nearest
main_part by shared-boundary length" logic is buggy in the presence of
single-part sibling divides, it can crush a single-part polygon to a
~100 m2 residual. For our pipeline, inter-divide topology (overlaps,
gaps between IDs) is handled separately by
[`hfutils::clean_gpkg_layer`](https://lynker-spatial.github.io/hfutils/reference/clean_gpkg_layer.html)
(Step 4b + 6d), so we only need the intra-divide "pick largest piece"
behavior here.

Orphan fragments that are dropped become geometry gaps, which the
downstream mapshaper pass absorbs into neighboring divides by
`gap_fill_area`.

## Usage

``` r
dissolve_multipart_divides(divides, id_col = "flowpath_id")
```

## Arguments

- divides:

  An `sf` object of divides with a `flowpath_id` (or equivalent ID)
  column and polygon geometry.

- id_col:

  Name of the ID column. Default `"flowpath_id"`.

## Value

Input `sf` with MULTIPOLYGON divides reduced to their largest POLYGON
component. Single-part polygons are untouched.
