# Test suite

The package ships with unit and functional tests in `tests/testthat/`.
They build their inputs in-line or from the small fixtures in
`tests/fixtures/`, so they need no network and no external data.

## Running the tests

``` r

library(testthat)
library(hydrofabric)
test_dir("tests/testthat")
```

Two end-to-end tests, a full aggregation of a fixture basin and a merge
of per-VPU products, run only when `NOT_CRAN=true` is set; the
aggregation test also needs mapshaper on `PATH`. Without them they are
skipped, not passed.

## Network topology

The acyclic and non-dendritic network tests are in `test-core-utils.R`:
a cycle is detected, a non-dendritic (divergent) network is accepted as
a valid graph, and a cyclic network is made acyclic. The same property
is checked on every product by the invariant suite.

## Checking a product with hfutils

A finished product can be checked with `hfutils` alone, without this
package or its inputs.
[`hf_check_invariants()`](https://lynker-spatial.github.io/hfutils/reference/hf_check_invariants.html)
runs the structural checks on one GeoPackage, and
[`hf_check_merge_invariants()`](https://lynker-spatial.github.io/hfutils/reference/hf_check_merge_invariants.html)
adds the conservation and completeness checks for a merged domain:

``` r

library(hfutils)

gpkg <- "prvi_ngen.gpkg"
fp   <- sf::read_sf(gpkg, "flowpaths")
dv   <- sf::read_sf(gpkg, "divides")
nx   <- sf::read_sf(gpkg, "nexus")

hf_network_is_dag(fp)                          # TRUE when the network is acyclic
hf_check_invariants("ngen", flowpaths = fp, divides = dv, nexus = nx,
                    strict = FALSE)
hf_check_merge_invariants(list(flowpaths = fp, divides = dv, nexus = nx),
                          gpkg = gpkg)
```

Each returns a list with `ok` and one entry per check. `strict = TRUE`
stops at the first failure.
[`hydrofabric::check_invariants()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_invariants.md)
runs the full pipeline suite, which adds the checks that need the
build’s intermediates.

## What is covered

| Area | What the tests assert | Files |
|----|----|----|
| Network topology | acyclic graphs, hydrosequence, split-reach chains, member contiguity, drainage area accumulating downstream, the upstream index | `core-utils`, `reconciled-toid-split-chain`, `member-contiguity`, `total-da-accumulates`, `upstream-index` |
| Input checks | reference ids and topology, coincident or mispaired reaches, severed terminals, grid coverage, input manifests, hydrolocation sets | `check-reference-*`, `grid-coverage`, `input-manifest`, `hydrolocation-input` |
| Divides | divides partitioned from the grids, one divide per flowpath at most, slivers and orphans absorbed with area conserved, split pieces kept with their reach | `pure-partition`, `partition-contract`, `absorb-*`, `allocate-parent-parts`, `normalize-divide-ownership`, `split-*`, `merge-bisected-divides`, `prune-ghost-divides`, `agg-divides-single-row`, `aggregated-divide-link`, `clear-disjoint-member-test`, `disjoint-swap-ownership`, `divide-flowpath-consistent`, `divide-member-flowline`, `final-multipart-cleanup`, `report-layer-effect` |
| Aggregation | collapse toward the target divide size, and an end-to-end aggregation that passes its invariants | `collapse-ideal-size`, `stage2-regression`, `stage2-helpers` |
| Hydrolocations and POIs | snapping, source-link plausibility, type filtering, placement at reach outlets, integer `poi_id`s, drainage-area reporting, fabrics with no POIs | `hydrolocations-prep`, `exclude-hl-types`, `pick-upstream`, `poi-nearer-end`, `normalize-poi-id`, `poi-drainage-area-report`, `zero-poi-safe`, `link-plausibility` |
| Channel attributes | RiverML channel geometry, flowpath slope, FAC sentinels, lake ids | `dingman-hydraulics`, `add-slope`, `fac-sentinels`, `close-lake-id-gaps` |
| Identifiers and merge | integer id columns, divide and flowpath id agreement, national `mainstem_id`, the domain merge, VPU runs that record their invariants | `enforce-integer-ids*`, `id-types-not-double`, `divide-id-identity`, `national-mainstem`, `merge-integration`, `run-vpus*`, `invariants-and-lineage` |
| Build tooling | checkpoint reuse, product comparison, installation from public packages alone, table helpers | `checkpoint-fingerprint`, `compare-products`, `standalone-dependencies`, `separate-longer` |
| Licensing | the attribution stamped into each product, including grid and lake credits | `attribution-policy`, `grid-credit-inherited`, `lakes-credit` |

File names omit the `test-` prefix.
