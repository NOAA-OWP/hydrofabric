# Compare a rebuilt product against the one it replaces

Diffs two NextGen GeoPackages on the things an intra-build invariant
cannot see: whether features were lost, and whether a flowpath that
owned a divide still owns one. A flowpath silently reclassified as a
connector passes every check in the suite, so a rebuild that drops
divides can otherwise look clean.

## Usage

``` r
hf_compare_products(old_gpkg, new_gpkg, area_tol = 0.01, strict = FALSE)
```

## Arguments

- old_gpkg, new_gpkg:

  Paths to the previous and rebuilt products.

- area_tol:

  Relative change in a flowpath's divide area to report. Default `0.01`
  (1%).

- strict:

  If `TRUE`, raise an error when a flowpath loses its divide or a layer
  loses features. Default `FALSE` reports and returns.

## Value

Invisibly, a list with `ok`, `counts`, `lost_divide`, `gained_divide`
and `area_changed`.

## What this can and cannot decide

A difference is only evidence of a regression when the inputs were the
same. Hold them constant and any loss is attributable to the code: that
is the case this was written for, and the one where `strict = TRUE`
belongs.

Across an input change it reports rather than decides. Hydrolocations
force reach splits, so a different hydrolocation set legitimately adds,
removes and resizes divides; a divide whose POI is no longer in the set
is *supposed* to disappear. Used as a gate there it fails a correct
build, so read the output and account for the differences instead.

With both changed at once nothing is attributable. When a new input set
and new code both need to land, build one unit on the new inputs with
the current code first and compare that, one VPU is a couple of hours
against a national run's tens of hours, and it buys an interpretable
answer.

## Examples

``` r
if (FALSE) { # \dontrun{
# Same inputs, changed code: a loss is a regression, so gate on it.
hf_compare_products("before/ngen_vpu04.gpkg", "after/ngen_vpu04.gpkg",
                    strict = TRUE)

# New hydrolocation set: report and account for the differences.
hf_compare_products("baseline/ngen_vpu04.gpkg", "rebuild/ngen_vpu04.gpkg")
} # }
```
