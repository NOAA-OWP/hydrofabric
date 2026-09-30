# The invariant model

The invariant suite tests the structural properties a consumer depends
on. Every phase records its results;
[`nextgen_fabric()`](https://lynker-spatial.github.io/hydrofabric/reference/nextgen_fabric.md)
stops the build when a gating check fails, and
[`merge_vpus()`](https://lynker-spatial.github.io/hydrofabric/reference/merge_vpus.md)
runs its own suite on each merged domain, adding cross-VPU conservation
and completeness checks.

## Gate or information

A **gate** is a property the build guarantees; its failure stops the
build. An **informational** check reports a desirable property without
stopping the build. A check becomes a gate once it holds everywhere and
a step in the build keeps it that way.

## What each phase checks

| Phase | What is checked |
|----|----|
| Refactor | unique ids, every reference reach reconciled, reconciled reaches with members, flow direction, a divide for every reconciled reach, landscape divides kept, valid divide polygons, POIs at reach outlets |
| Aggregate | unique flowpaths and divides, resolvable `toid`s and divide links, single-part flowpaths each lying in its divide, reconciled reaches assigned to exactly one flowpath, divide area conserved, no overlap, no shared divide |
| Nextgen | `fp-` and `cat-` ids, acyclic flowpath-nexus and flowline graphs, resolvable `toid`s, at most one divide per flowpath and an owner for every network divide, divide coverage and area conserved with no overlaps or gaps, flowpath members contiguous and agreeing with the flowlines, a closed `routeable` network, every reference id traceable, drainage area accumulating downstream, `mainstem_id` present, positive slope, integer ids, lake membership, the POI contract, the upstream index |
| Merge | id uniqueness and prefixes across the domain, resolvable `toid`s, divide ownership, an acyclic graph, every VPU present, feature counts, divide area and drainage area conserved, consistent CRS, and the nextgen network checks on the merged network |

The POI contract is described in
[`vignette("poi-contract")`](https://lynker-spatial.github.io/hydrofabric/articles/poi-contract.md).

## Running the suite

The same suite runs on any GeoPackage:

``` r

check_invariants(
  "ngen",
  flowpaths = sf::read_sf("ngen_vpu01.gpkg", "flowpaths"),
  divides   = sf::read_sf("ngen_vpu01.gpkg", "divides"),
  nexus     = sf::read_sf("ngen_vpu01.gpkg", "nexus"),
  flowlines = sf::read_sf("ngen_vpu01.gpkg", "flowlines"),
  strict    = TRUE
)
```

`strict = FALSE` collects every result instead of stopping at the first
failure.

## The record

Each build writes `invariants.json` with every check’s name, level and
result: in `work/` for a VPU, and beside each merged domain product. A
check whose input was not supplied leaves no entry, so absence means
*not evaluated*, not *passed*; the record’s `inputs` block lists what
the suite was given.

To compare a rebuild against the product it replaces, use
[`hf_compare_products()`](https://lynker-spatial.github.io/hydrofabric/reference/hf_compare_products.md),
which names any flowpath that lost its divide.
