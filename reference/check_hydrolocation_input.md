# Validate a hydrolocation set before building from it

Checks a hydrolocation GeoPackage against what the pipeline requires of
it. Gated checks are properties a build depends on; informational ones
report shape without judging it.

Gated: the layer is named so the pipeline can read it; geometry is POINT
and non-empty; the CRS is defined; `hl_reference` and `flowpath_id` are
populated; no row is an exact duplicate of another; every `vpuid` names
a unit some supplied reference carries; and every `flowpath_id` resolves
to a reach in its domain's reference.

Reported: how many rows share an `hl_reference` (collocation is by
design, since several hydrolocations may collapse onto one POI), how
many arrive with a `poi_id`, and the mix of `hl_type`.

## Usage

``` r
check_hydrolocation_input(
  hl_gpkg,
  layer = "hydrolocations",
  references = list(),
  check_distance = TRUE,
  max_dist_m = 200,
  strict = FALSE
)
```

## Arguments

- hl_gpkg:

  Path to the hydrolocation GeoPackage.

- layer:

  Layer to read. Default `"hydrolocations"`, which is the name the
  pipeline reads; a set whose layer is called anything else is rejected
  unless this is set explicitly, because the build cannot find it.

- references:

  Named list of references to resolve against, each a list of `ids`
  (character `flowpath_id`) and `vpuids`. Built by
  [`hl_reference_index()`](https://lynker-spatial.github.io/hydrofabric/reference/hl_reference_index.md)
  when not supplied directly.

- check_distance:

  If `TRUE` (default), also measure how far each point lies from the
  reach it names. Resolution is not placement: an identifier can exist
  and still name a reach on the other side of the domain, which is what
  a statewide point set indexed against a clipped reference produces.

- max_dist_m:

  Reporting threshold for that distance, in metres. Default `200`. The
  result is informational, not a gate: across the NWM v4 inputs the
  population is nowhere zero, and a threshold tight enough to catch
  genuine mis-attribution also flags several thousand forecast points
  that sit a few hundred metres off the reference line and snap cleanly.
  Read the spread the check reports, a median near zero with a long tail
  is snapping noise, a median in the kilometres is a point set indexed
  against the wrong domain.

- strict:

  If `TRUE`, a failed gate raises an error. Default `FALSE` reports
  every check and returns.

## Value

Invisibly, a list with `ok`, `stage` and one entry per check.

## Examples

``` r
if (FALSE) { # \dontrun{
refs <- list(conus = hl_reference_index("conus_reference.gpkg"))
check_hydrolocation_input("merged_hydrolocations.gpkg", references = refs)
} # }
```
