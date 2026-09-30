# Prune ghost divides (areas below a threshold)

Removes network divides whose `areasqkm` is below a threshold. These are
usually produced by imbalanced PRD watershed splits where two sibling
reconciled flowlines share a parent reference COMID and one side
receives essentially no area (often \<500 m2). The corresponding
aggregated flowpath is a through-reach with no meaningful local
drainage; keeping a microscopic "ghost" polygon makes the visualization
misleading and the schema inconsistent (has_divide = TRUE with a 100 m2
polygon).

For each ghost divide: the divide row is dropped and its flowpath is
updated with `has_divide = FALSE`, `divide_id = NA`, `areasqkm = 0`. A
subsequent `hfutils::clean_gpkg_layer(..., gap_fill_area = > threshold)`
call will absorb the resulting tiny gaps into neighbor divides.

## Usage

``` r
prune_ghost_divides(
  input_gpkg,
  threshold_sqkm = 0.001,
  divide_layer = "divides",
  flowpath_layer = "flowpaths"
)
```

## Arguments

- input_gpkg:

  Path to the GeoPackage (must contain `divides` and `flowpaths`
  layers).

- threshold_sqkm:

  Divides with `areasqkm` strictly less than this are pruned. Default
  `0.001` (1000 m2).

- divide_layer:

  Name of the divides layer. Default `"divides"`.

- flowpath_layer:

  Name of the flowpaths layer. Default `"flowpaths"`.

## Value

Invisibly, a list of pruned `divide_id` values.
