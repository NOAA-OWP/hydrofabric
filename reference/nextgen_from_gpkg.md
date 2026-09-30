# Build a NextGen hydrofabric from a single GeoPackage

A self-contained function that reads all inputs–flowpaths, divides,
network, pois, flowlines, and events–from a single GeoPackage produced
by the `add_flowlines` step of the pipeline. No separate `flowline_gpkg`
is required because `add_flowlines` writes full-resolution split
segments (with `id`, `toid`, `flowpath_id`, `areasqkm`, `fac_value`,
etc.) directly into the primary GeoPackage.

## Usage

``` r
nextgen_from_gpkg(
  flowpath_gpkg,
  collapse_method = "upstream",
  nexus_prefix = "nex-",
  terminal_nexus_prefix = "tnx-",
  divide_prefix = "cat-",
  flowpath_prefix = "fp-",
  term_add = 1e+09,
  outfile = NULL
)
```

## Arguments

- flowpath_gpkg:

  Path to a GeoPackage containing at minimum: `flowpaths`, `divides`,
  `network`. Should also contain `pois`, `flowlines` (from
  `add_flowlines`), and `events`.

- collapse_method:

  Character, `"upstream"` (default) or `"downstream"`. Controls nexus ID
  keying at confluences. `"upstream"` keys the nexus on the
  **receiving** flowpath (`nex-X` where X = toid) so all donors to fp X
  share one nexus: fewer nexuses, shared at confluences. `"downstream"`
  keys the nexus on the **donating** flowpath (`nex-X` where X =
  flowpath_id) so each donor gets its own outlet nexus, one nexus per
  fp. Terminal nexuses are unaffected by this parameter.

- nexus_prefix:

  Character prefix for nexus IDs. Default `"nex-"`.

- terminal_nexus_prefix:

  Prefix for terminal nexus IDs. Default `"tnx-"`.

- divide_prefix:

  Prefix for divide IDs. Default `"cat-"`.

- flowpath_prefix:

  Prefix for flowpath IDs. Default `"fp-"`.

- term_add:

  Numeric offset added to terminal `flowpath_toid` values. Default
  `1e9`.

- outfile:

  Optional path for output GeoPackage. If `NULL` (default) returns the
  list invisibly instead of writing.

## Value

If `outfile` is specified, writes the result to that GeoPackage and
returns `NULL` invisibly. Otherwise returns a named list with elements
`nexus`, `flowpaths`, `divides`, `flowlines`, `incremental_areas`,
`pois` (when present), and `network`.

## Details

Topology steps:

1.  Build nexus geometries from flowpath endpoints.

2.  Identify terminal and headwater nexus nodes.

3.  Assign prefixed IDs to flowpaths, divides, and nexuses.

4.  Join the sub-flowline segments (from the `flowlines` layer) to their
    parent aggregated flowpaths via `flowlines$flowpath_id`.

5.  Compute `flowline_incr` (FAC-derived incremental drainage area in
    km2 per sub-flowline:
    `total_dasqkm[i] - max(total_dasqkm[upstream donors])`).

6.  Build `incremental_areas` from divides joined to flowlines.

7.  Join hydrolocations from the `events` layer into the network table.

## Examples

``` r
if (FALSE) { # \dontrun{
nextgen_from_gpkg("out_vpu01.gpkg", outfile = "ngen_vpu01.gpkg")
ng <- nextgen_from_gpkg("out_vpu01.gpkg")
} # }
```
