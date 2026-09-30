# Check that the supplied grids cover the domain

Compares each raster's extent against the reference fabric's extent. A
primary grid that does not cover the domain is an error: the most common
cause is a per-domain configuration inheriting a grid meant for
somewhere else, which otherwise builds without complaint. A fallback
that does not reach the domain is reported and left alone, since it
simply never fires.

## Usage

``` r
check_grid_coverage(
  reference_gpkg,
  fdr_path,
  fac_path,
  fallback_fdr_path = NULL,
  fallback_fac_path = NULL,
  vpuid = NA_character_,
  min_coverage = 0.99,
  strict = TRUE
)
```

## Arguments

- reference_gpkg:

  Path to the reference fabric.

- fdr_path, fac_path:

  Primary grid paths.

- fallback_fdr_path, fallback_fac_path:

  Optional fallback grid paths.

- vpuid:

  Optional VPU being built. When supplied and the reference carries a
  `vpuid` column, coverage is measured against that VPU's extent rather
  than the whole reference, so a per-VPU grid tile is judged on the
  ground it is meant to cover.

- min_coverage:

  Minimum share of the reference extent a primary grid must cover.
  Default `0.99`.

- strict:

  If `TRUE` (default) an uncovered primary grid raises an error.

## Value

Invisibly, a named numeric of coverage fractions.
