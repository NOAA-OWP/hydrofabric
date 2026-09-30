# Make each lake's flowpaths a contiguous subgraph

t-route collapses every flowpath sharing a `lake_id` into one waterbody
node. That is only valid if those flowpaths are connected: a non-member
sitting between two members turns `lake -> gap -> lake` into
`node_L -> gap -> node_L`, a cycle, and routing with waterbodies enabled
fails on it.

Membership is a length-weighted majority of a flowpath's member reaches,
a rule that decides each flowpath alone and cannot see the network: a
reach mostly outside the lake is rejected even when both its neighbours
are members. An earlier closure runs before divide-less leaf flowpaths
are demoted, which rewires `flowpath_toid` and can open a gap that was
not there before, so this runs once the topology is final.

It only ever *adds* a `lake_id` to a flowpath already bracketed by one
lake, so it cannot change a fabric with no gaps and cannot move a
flowpath between lakes.

The stamp is propagated to `flowlines` and `network`, which inherit
`lake_id` from their flowpath and would otherwise disagree with it.

## Usage

``` r
close_lake_id_gaps(
  input_gpkg,
  max_step = 5L,
  flowpath_layer = "flowpaths",
  nexus_layer = "nexus",
  flowline_layer = "flowlines",
  network_layer = "network"
)
```

## Arguments

- input_gpkg:

  Path to the NextGen GeoPackage.

- max_step:

  Longest run of non-members to bridge. Beyond this the gap is a genuine
  river between two arms of one waterbody and is left to be reported
  rather than silently absorbed.

- flowpath_layer, nexus_layer, flowline_layer, network_layer:

  Layer names.

## Value

Invisibly, a character vector of the `flowpath_id`s stamped.
