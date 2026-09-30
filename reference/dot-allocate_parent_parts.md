# Allocate a split parent's displaced land to the child reach nearest each part

When a reference catchment is split among the fragments of a split
reach, the D8 partition can return sub-polygons that sit nowhere near
the fragment they were assigned to. The land itself is correct, it is
the parent catchment and no area is missing, but the assignment of land
to fragment is not.

This resolves the assignment the way the partition is itself trying to:
every polygon part of the displaced land goes to the fragment whose own
flowline is nearest it. Fragments whose partition already meets their
flowline keep what they hold and may gain parts, so only displaced land
moves, the parent's area is conserved exactly, and no two fragments can
claim the same ground.

The alternatives it replaces were both spatially arbitrary: handing the
whole parent to whichever fragment had the longest flowline, and handing
every displaced fragment the same leftover polygon.

A fragment with no flowline cannot receive land, since there is nothing
to measure against; a group in which no fragment has one is left
untouched.

## Usage

``` r
.allocate_parent_parts(geoms, lines, disjoint, crs)
```

## Arguments

- geoms:

  List of `sfg`, the current partition polygon of each fragment.

- lines:

  List of `sfc` or `NULL`, each fragment's member flowlines.

- disjoint:

  Logical marking fragments whose polygon does not meet their own
  flowline.

- crs:

  CRS shared by `geoms` and `lines`.

## Value

A list of `sfg` the same length as `geoms`, or `NULL` when there is
nothing to reallocate.
