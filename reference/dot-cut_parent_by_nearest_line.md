# Cut a split parent's catchment among its fragments by nearest reach

The partition of a split reference catchment among the fragments of the
reach is a nearest-channel assignment: land drains to the piece of
channel closest to it. The D8 trace computes that by following flow
directions, which needs the vector reach and the flow grid to agree on
where the channel is. Where they do not, Alaska routes TDX-Hydro reaches
over HydroSHEDS grids, the trace lands elsewhere and returns
sub-polygons displaced from the fragment they were assigned to, or none
at all.

This computes the same assignment geometrically. Each fragment's line is
sampled at `step` metres, the samples generate a nearest-point
tessellation, cells are dissolved by the fragment that generated them,
and the result is clipped to the parent. Three properties hold by
construction rather than by repair: the pieces tile the parent exactly,
so area is conserved; no two fragments can claim the same ground; and
every fragment whose line meets the parent receives a share, so none is
left without a catchment.

This is a partition of one real catchment among its own children,
clipped to that catchment. It creates no land and fills no gap between
divides.

## Usage

``` r
.cut_parent_by_nearest_line(parent, lines, crs, step = 25)
```

## Arguments

- parent:

  `sfg`, the parent catchment to divide.

- lines:

  List of `sfc` or `NULL`, each fragment's member flowlines.

- crs:

  CRS shared by `parent` and `lines`.

- step:

  Spacing in CRS units at which each line is sampled. 25 m sits below
  the 30 m DEM cell, so the cut resolves detail the grid cannot.

## Value

A list of `sfg` the same length as `lines`, or `NULL` if the cut cannot
be computed and the caller should leave the group alone.
