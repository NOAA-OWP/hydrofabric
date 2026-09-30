# Flag whether each flowline is a direct constituent of its flowpath's line

Adds a logical `flowpath_constituent` attribute to the `flowlines`
layer. A flowline is a *constituent* (`TRUE`) when its geometry lies on
the line that its aggregate flowpath traces from outlet to head: i.e. it
is directly part of `flowpath_id`. Member flowlines that enter the
flowpath as tributaries (a branch that drains into the stem but whose
channel the flowpath line does not trace) are `FALSE`.

The test is geometric and matches the flowpath line by construction:
[`extend_flowpaths_to_mainstem()`](https://lynker-spatial.github.io/hydrofabric/reference/extend_flowpaths_to_mainstem.md)
builds each flowpath line by merging its stem member reaches, so a
constituent's geometry is a sub-segment of the line (coverage ~1) while
a tributary branches away (coverage low). A flowline is a constituent
when at least `cov_threshold` of its length lies within `tol_m` of its
flowpath's line. Every flowpath is guaranteed at least one constituent:
if no member clears the threshold (a flowpath still carrying the
aggregated original, non-member line), its single best-covered member is
flagged.

## Usage

``` r
flag_flowpath_constituents(
  input_gpkg,
  flowpath_layer = "flowpaths",
  flowlines_layer = "flowlines",
  tol_m = 5,
  cov_threshold = 0.5
)
```

## Arguments

- input_gpkg:

  Path to a GeoPackage containing `flowpath_layer` and
  `flowlines_layer`.

- flowpath_layer, flowlines_layer:

  Layer names. Defaults `"flowpaths"` / `"flowlines"`.

- tol_m:

  Buffer (m) around the flowpath line within which a flowline segment
  counts as "on the line". Default 5.

- cov_threshold:

  Fraction of a flowline's length that must lie on the line for it to be
  a constituent. Default 0.5.

## Value

Invisibly, the number of constituent flowlines. Adds
`flowpath_constituent` (logical) to `flowlines_layer`, overwriting it in
place.

## Routing contract

OWP routes on the flowlines coincident with flowpaths
(`flowpath_constituent == TRUE`), that subset IS the routing network,
not a descriptive overlay. It must be closed under `flowline_toid`: a
coincident flowline's downstream must resolve to another coincident
flowline. Do not fix a routing gap by filtering the mask; dropping a
non-coincident connector reach that carries a coincident reach's
downstream link would leave the routing set disconnected. A coincident
flowline whose `flowline_toid` points to a non-coincident reach is a
routing-set gap (detected by the `coincident_flowline_toid_closed`
invariant); the remedy is to re-point the routing toid across the
intervening non-coincident reaches to the next coincident flowline.
`flowpath_constituent == FALSE` reaches are off-network tributary /
connector geometry; a flowline's `flowline_toid` may legitimately name a
reach in a different flowpath (a tributary confluence).

## Examples

``` r
if (FALSE) { # \dontrun{
flag_flowpath_constituents("ngen_vpu01.gpkg")
} # }
```
