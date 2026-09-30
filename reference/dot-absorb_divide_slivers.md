# Absorb a divide's sliver parts into the neighbour they share the most boundary with

`construct_agg_divides` unions each aggregate's member catchments, so an
aggregate whose members are not spatially contiguous gets a genuinely
multi-part divide. Most of that multi-part-ness is not a real disjoint
catchment: on CONUS, 1,192 of 1,300 multi-part divides have a second
part under 1% of the divide (median 900 m2). Those are slivers left by
clipping, snapping and union arithmetic.

Nothing downstream reduces them. `dissolve_multipart_divides` runs in
`reconcile_divides`, which is BEFORE the aggregation that creates the
fragmentation, and `hfutils::clean_gpkg_layer(gap_fill_area = )` fills
holes BETWEEN divides, a divide's own second part is an owned polygon,
not a gap, so it is never a gap-fill candidate.

This pass moves each qualifying sliver into the adjacent divide it
shares the most boundary with, which is the divide that land most
plausibly drains to. It is area-conserving by construction: the polygon
changes owner, it is never dropped. That distinction matters, the
dissolve-and-drop approach used at the refactor phase would forfeit 285
km2 nationally, and for a genuinely disjoint divide it would discard a
real catchment.

Deliberately conservative. A part must be BOTH absolutely small
(`min_part_km2`) and small relative to its own divide (`max_part_frac`)
before it is eligible, so the ~53 divides whose second part is a
substantial fraction of a real catchment keep their multi-part geometry.
A sliver with no adjacent divide is kept as well: carrying a multi-part
divide is better than losing the land.

## Usage

``` r
.absorb_divide_slivers(
  div,
  min_part_km2 = 0.01,
  max_part_frac = 0.01,
  min_shared_m = 1,
  flowlines = NULL
)
```

## Arguments

- div:

  Divides `sf`.

- min_part_km2:

  Absolute ceiling on an absorbable part, km2. Default 0.01 (1 hectare),
  which covers 1,194 of the 1,300 CONUS cases.

- max_part_frac:

  Ceiling on an absorbable part as a fraction of its own divide's total
  area. Default 0.01.

- min_shared_m:

  Minimum shared boundary with the recipient, metres. A part touching
  only at a point would union into a still-disjoint polygon.

## Value

A list with `div` (repaired layer), `n_absorbed`, `n_kept` (eligible
slivers with no usable neighbour), `area_moved_km2`, and
`n_multi_before` / `n_multi_after`.
