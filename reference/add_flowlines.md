# Add Sub-Flowline Segments to Aggregated Hydrofabric

Takes the full-resolution reconciled flowlines (as passed to
`aggregate_to_distribution`), splits any that exceed
`split_flines_meters`, assigns each segment its parent aggregated
`flowpath_id` and `divide_id` via the `network` lookup table in
`input_gpkg`, remaps topology so the resulting segments form a valid
directed graph, and attributes each segment with total and incremental
contributing area derived from a flow-accumulation raster (FAC).

If `events` is supplied (the full event set from
`rf_generate_events_from_hydrolocations_fast`), each event is assigned
the aggregated `flowpath_id` and `divide_id` of the reconciled flowline
it belongs to, and written as an `events` layer to `input_gpkg`.

## Usage

``` r
add_flowlines(
  flowlines,
  input_gpkg,
  split_flines_meters = 2000,
  fac,
  fallback_fac = NULL,
  poi_events = NULL,
  junction_events = NULL
)
```

## Arguments

- flowlines:

  An `sf` object of full-resolution reconciled flowlines (the data
  passed to `aggregate_to_distribution`). Must have integer columns `ID`
  and `toID`, and a character column `member_flowpath_id`
  (comma-separated NHD IDs per reconciled flowline).

- input_gpkg:

  Path to a GeoPackage produced by `aggregate_to_distribution` (must
  contain a `network` table). The `flowlines` (and optionally `events`)
  layer will be written here.

- split_flines_meters:

  Numeric. Maximum individual segment length in metres. Default `2000`.

- fac:

  Character path to a flow-accumulation GeoTIFF, or a
  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  already in memory.

- fallback_fac:

  Optional character path (or `NULL`) to a fallback FAC raster used when
  `fac` returns NA after the standard 2-cell buffer fill. Values filled
  from the fallback are tagged `fac_source` naming that grid in the
  output. When `fallback_fac` resolves to the same file as `fac` (path
  equality, i.e. a domain configured with the fallback grid as its
  primary), the fallback lookup is skipped automatically. Default
  `NULL`.

- poi_events:

  Optional `sf` object of snapped POI events returned by
  `refactor_flowlines` (i.e. `out$events`). Must contain a
  `reconciled_id` column (integer reconciled flowline ID) and `poi_id`.
  These events carry geometries repositioned onto the reconciled
  network.

- junction_events:

  Optional `sf` object of junction events from
  `rf_generate_events_from_hydrolocations_fast` (i.e.
  `filter(events, event_kind == "junction")`). Must contain a
  `flowpath_id` column with original NHD flowline IDs.

## Value

Invisibly, the `sf` data frame of flowline segments. The layer
`flowlines` is written/overwritten in `input_gpkg`. If `events` is
provided, an `events` layer is also written/overwritten.

## Details

Workflow:

1.  Read the `network` table from `input_gpkg` and expand its
    comma-separated `flowline_id` column to build a per-ID lookup:
    reconciled flowline ID -\> aggregated `flowpath_id`,
    `flowpath_toid`, `divide_id`.

2.  Split `flowlines` segments that exceed `split_flines_meters` using
    `split_flowlines`.

3.  Assign `flowpath_id` / `divide_id` to each segment via the lookup.

4.  Remap topology: segments within the same parent flowpath chain
    upstream-to-downstream; the terminal segment of each chain points to
    the first segment of the downstream chain (by `flowpath_toid`).

5.  Extract the FAC value at the downstream endpoint of each segment to
    derive `tot_drainage_areasqkm` (FAC \* cell area in km2).

6.  Compute `areasqkm` (incremental) as `tot_drainage_areasqkm` of this
    segment minus the maximum `tot_drainage_areasqkm` of immediate
    upstream segments. Head segments with no upstream neighbour use
    their own total.

7.  If `poi_events` and/or `junction_events` are provided, each is
    attributed with the aggregated `flowpath_id` and `divide_id`, then
    the combined set is written as an `events` layer in `input_gpkg`.
    POI events are mapped via their `reconciled_id` (integer reconciled
    flowline ID from `refactor_flowlines`) directly through the network
    lookup. Junction events are mapped via their NHD `flowpath_id`
    column through `flowlines$member_flowpath_id` to the reconciled
    integer ID.

8.  Write flowline segments as a `flowlines` layer to `input_gpkg` and
    return the `sf` object invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
reconciled <- sf::st_read("reconciled.gpkg", "reconciled")
add_flowlines(flowlines = reconciled,
              input_gpkg = "out_vpu01.gpkg",
              split_flines_meters = 2000,
              fac = "fac.tif")
} # }
```
