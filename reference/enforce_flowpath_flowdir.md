# Enforce upstream-\>downstream orientation on aggregated flowpath geometry

[`extend_flowpaths_to_mainstem()`](https://lynker-spatial.github.io/hydrofabric/reference/extend_flowpaths_to_mainstem.md)
rebuilds each aggregate flowpath line with `geos_line_merge`, which does
not preserve flow direction, and (unlike the reconciled `flowlines`,
guarded by the `reconciled_flow_direction` invariant) the aggregated
`flowpaths` layer has no orientation check, so a fraction of flowpaths
ship reversed (their geometry runs downstream-\>upstream). This pass
re-orients every flowpath so its OUTLET is the LAST vertex.

The true outlet is located from the reliably-oriented `flowlines` layer:
each flowpath's outlet sub-flowline is the member whose `flowline_toid`
leaves the flowpath, and its downstream (end) vertex is the flowpath
outlet. A flowpath whose START is nearer that outlet than its END is
reversed and is flipped with
[`sf::st_reverse()`](https://r-spatial.github.io/sf/reference/geos_unary.html).
Genuinely branching multi-part geometry (no single meaningful direction)
is left untouched.

## Usage

``` r
enforce_flowpath_flowdir(
  gpkg,
  flowpath_layer = "flowpaths",
  flowlines_layer = "flowlines"
)
```

## Arguments

- gpkg:

  Path to a GeoPackage containing `flowpath_layer` and
  `flowlines_layer`.

- flowpath_layer, flowlines_layer:

  Layer names. Defaults `"flowpaths"` / `"flowlines"`.

## Value

Invisibly, the number of flowpaths reversed. Overwrites `flowpath_layer`
in place.

## Examples

``` r
if (FALSE) { # \dontrun{
enforce_flowpath_flowdir("ngen_vpu01.gpkg")
} # }
```
