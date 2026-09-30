# Add a nested-set upstream index to a hydrofabric GeoPackage

Stamps `upstream_id` (depth-first pre-order ordinal from each outlet)
and `num_upstreams` (strict upstream count) on the flowpaths, flowlines,
divides, and network layers, so the flowpaths upstream of X are those
with `upstream_id` in
`(upstream_id[X], upstream_id[X] + num_upstreams[X]]`, an O(1) range
filter instead of a traversal.

## Usage

``` r
add_upstream_index(
  ngen_gpkg,
  flowpaths_layer = "flowpaths",
  nexus_layer = "nexus",
  propagate_to = c("flowlines", "divides", "network")
)
```

## Arguments

- ngen_gpkg:

  Path to the NextGen GeoPackage.

- flowpaths_layer:

  Flowpath layer name. Default `"flowpaths"`.

- nexus_layer:

  Nexus layer name. Default `"nexus"`.

- propagate_to:

  Layers to stamp by `flowpath_id`. Default
  `c("flowlines", "divides", "network")`.

## Value

Invisibly, a data.frame of `flowpath_id`, `upstream_id`, and
`num_upstreams`. Adds those integer columns to the layers in place.

## Details

The values are specific to this file's network (a per-VPU, merged, or
subset gpkg each gets a correct-scope index) and are not a persistent
key. The compute (fp -\> nexus -\> fp resolution + nested set) is
[`hfutils::hf_upstream_index()`](https://lynker-spatial.github.io/hfutils/reference/hf_upstream_index.html);
reads and writes are attribute-only via DBI, so it is safe at
merged/national scale.

## Examples

``` r
if (FALSE) { # \dontrun{
add_upstream_index("ngen_vpu01.gpkg")
# flowpaths strictly upstream of fp-8:
library(DBI)
con <- dbConnect(RSQLite::SQLite(), "ngen_vpu01.gpkg")
anchor <- dbGetQuery(con, "SELECT upstream_id u, num_upstreams n
                           FROM flowpaths WHERE flowpath_id = '8'")
dbGetQuery(con, sprintf(
  "SELECT flowpath_id FROM flowpaths WHERE upstream_id > %d AND upstream_id <= %d",
  anchor$u, anchor$u + anchor$n))
} # }
```
