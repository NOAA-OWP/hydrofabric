# Flag the routeable flowlines (the routing network)

Adds a logical `routeable` attribute to the `flowlines` layer marking
the reaches in the routing network: the flowpath-constituent
(coincident) flowlines plus the incremental reaches that bridge a
coincident flowline to the next coincident flowline downstream. Routing
follows the native `flowline_toid`; the `routeable == TRUE` subgraph is
closed under it (every non-terminal routeable flowline drains to a
routeable flowline), so a consumer can route on the routeable set
without a coincident reach draining into a non-coincident connector that
is filtered out ("flying water").

A bridge reach is a real channel carrying flow between mainstems; it is
non-coincident only because aggregation traced a different mainstem.
When a coincident flowline never rejoins a coincident flowline
downstream, the bridge is extended to the terminal instead, so the
reach's flow exits the domain rather than draining into a filtered-out
connector.

Requires `flowpath_constituent` on the flowlines layer (see
[`flag_flowpath_constituents()`](https://lynker-spatial.github.io/hydrofabric/reference/flag_flowpath_constituents.md),
which must run first).

## Usage

``` r
flag_routeable_flowlines(input_gpkg, flowlines_layer = "flowlines")
```

## Arguments

- input_gpkg:

  Path to a GeoPackage containing `flowlines_layer`.

- flowlines_layer:

  Layer name. Default `"flowlines"`.

## Value

Invisibly, the number of routeable flowlines. Adds `routeable` (logical)
to `flowlines_layer`, overwriting it in place.

## Examples

``` r
if (FALSE) { # \dontrun{
flag_routeable_flowlines("ngen_vpu01.gpkg")
} # }
```
