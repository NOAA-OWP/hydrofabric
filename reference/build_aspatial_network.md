# Build an aspatial long-form network crosswalk

Replaces the spatial `network` layer in a hydrofabric GeoPackage with an
aspatial attribute table keyed on `(flowline_id, reference_id)`, one row
per sub-flowline x reference-COMID member. The table is rebuilt from the
flowlines layer and the reconciled-group memberships so it always
matches the trimmed topology produced by `extend_flowpaths_to_mainstem`:
flowlines on a flowpath's mainstem levelpathid keep the mainstem's
`flowpath_id`; flowlines on absorbed tributaries (different levelpathid)
report a synthetic `fp-t<levelpathid>` identity and inherit the
mainstem's downstream nexus as their `flowpath_toid`.

Output columns:

- `flowpath_id`, `flowpath_toid` : aggregated flowpath topology

- `mainstem_id` : levelpathid of the flowpath's mainstem; group all
  flowlines on the same mainstem channel

- `is_mainstem` : TRUE when this flowline is on the flowpath's mainstem
  (levelpathid == mainstem_id); FALSE for absorbed tributaries

- `flowline_id`, `flowline_toid` : split-segment topology

- `reference_id` : originating reference COMID

- `reference_percent` : share of reference NHD length represented by
  this flowline's reconciled group

- `poi_id` : aggregation POI (nullable)

- `hydrolocation_id` : `hl_id` from hydrolocations

- `divide_id`, `divide_areasqkm` : divide this flowpath drains (all
  flowlines in a flowpath, mainstem and tributary: carry the same
  `divide_id`)

- `flowline_area` : drainage area at flowline outlet (from FAC
  extraction)

## Usage

``` r
build_aspatial_network(
  gpkg,
  refactored = NULL,
  reconciled = NULL,
  hydrolocations = NULL,
  verbose = FALSE
)
```

## Arguments

- gpkg:

  Path to the ngen GeoPackage produced by `nextgen_from_gpkg`.

- refactored:

  `sf` of refactored flowlines from `refactor_flowlines`
  (`out$refactored`). Used to compute `reference_percent` per reference
  COMID.

- reconciled:

  `sf` of reconciled flowlines (`out$reconciled`): must contain `ID` and
  `member_flowpath_id`. Used to map each flowline's levelpathid back to
  the reference COMIDs it represents.

- hydrolocations:

  `sf` of hydrolocations with `poi_id` and `hl_id`, or path to a gpkg
  with a `hydrolocations` layer. If `NULL`, `hydrolocation_id` is
  sourced from the `events` layer if present.

- verbose:

  If TRUE, emit a rich breakdown: row class counts, tributary-flowline
  share, reference_percent distribution, hydrolocation match rate,
  divide-area / flowline-area summaries. Default FALSE.

## Value

Invisibly, the path to `gpkg`.

## Examples

``` r
if (FALSE) { # \dontrun{
build_aspatial_network(
  gpkg           = "ngen_vpu01.gpkg",
  refactored     = sf::st_read("refactored.gpkg", "refactored"),
  reconciled     = sf::st_read("reconciled.gpkg", "reconciled"),
  hydrolocations = "hydrolocations.gpkg")
} # }
```
