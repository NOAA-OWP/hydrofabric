# Pipeline invariant checks

Pipeline invariant checks

## Usage

``` r
check_invariants(
  stage,
  ...,
  strict = TRUE,
  coverage = FALSE,
  coverage_min = 0.9,
  area = FALSE,
  area_baseline_m2 = NULL,
  coverage_baseline_m2 = NULL,
  area_tol_frac = 0.005
)
```

## Arguments

- stage:

  Character. One of `"refactored"`, `"reconciled"`, `"aggregated"`,
  `"ngen"`.

- ...:

  Stage-specific named arguments. For `"ngen"`: `flowpaths`, `divides`,
  `nexus`, and optionally `network` (the aspatial network table, used to
  check flowline_id global uniqueness, flowline_toid referential
  integrity, and flowline graph acyclicity) and `lakes` (the waterbody
  layer, used to check that lake_id stamps intersect their lake).

- strict:

  Logical. If `TRUE` (default) any failed check throws an error.

- coverage:

  Logical. Run expensive flowpath-in-catchment coverage check.

- coverage_min:

  Minimum fraction of a flowpath inside its catchment.

- area:

  Logical. Run the geometric non-overlap check (computes the union of
  all network divides: moderately expensive on large VPUs). The cheap
  total-area, area-conservation, and shared-divide checks always run for
  stages that carry a `divides` layer, regardless of this flag.

- area_baseline_m2:

  Numeric. Total area (m2) of the upstream divide set this stage must
  conserve (typically the sum of the input divides' areas). When
  supplied, `divides_area_conserved` asserts the output total stays
  within `area_tol_frac` of it. `NULL` (default) reports the total as
  info only.

- coverage_baseline_m2:

  Numeric. Union area (m2) of the upstream divide set, i.e. the ground
  it covers rather than the sum of its parts. When supplied,
  `divides_coverage_conserved` reports the change in footprint. This is
  what detects lost ground: removing an overlap and punching a hole both
  reduce the summed area, but only a hole reduces the union. Requires
  `area = TRUE`, which is when the union is computed.

- area_tol_frac:

  Numeric. Fractional tolerance for the area-conservation and
  non-overlap checks. Default `0.005` (0.5%).

## Value

Invisibly a list with `ok`, `stage`, and `checks`.

## Examples

``` r
if (FALSE) { # \dontrun{
check_invariants("ngen",
                 flowpaths = sf::st_read("ngen_vpu01.gpkg", "flowpaths"),
                 divides   = sf::st_read("ngen_vpu01.gpkg", "divides"),
                 nexus     = sf::st_read("ngen_vpu01.gpkg", "nexus"),
                 strict = FALSE, area = TRUE)
} # }
```
