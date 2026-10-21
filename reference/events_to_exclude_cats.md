# Build exclude_cats from event/split_event data

Returns a character vector of segment IDs that should be protected from
collapse because they host (or are immediately downstream of) an event.
Pass the result directly as the `exclude_cats` argument to
[`refactor_flowlines()`](https://NOAA-OWP.github.io/hydrofabric/reference/refactor_flowlines.md).

## Usage

``` r
events_to_exclude_cats(
  events = NULL,
  split_events = NULL,
  flowlines = NULL,
  protect_downstream = TRUE,
  only_split_candidates = FALSE
)
```

## Arguments

- events:

  sf or data.frame of events. Should have a `flowpath_id` column
  identifying the segment each event sits on.

- split_events:

  sf or data.frame of split events. Also expected to have `flowpath_id`.
  When both `events` and `split_events` are provided, IDs from both are
  used.

- flowlines:

  Optional sf or data.frame with flowpath topology columns
  (`flowpath_id`/`flowpath_id` and `flowpath_toid`/`flowpath_toid`).
  When supplied and `protect_downstream = TRUE`, the immediate
  downstream neighbor of each event-hosting segment is also included:
  matching the reference refactor-nhdplus behaviour which adds both the
  event flowpath_id and its `flowpath_toid` to `exclude_cats`.

- protect_downstream:

  Logical (default `TRUE`). Include the downstream (`toid`) neighbor of
  each event segment when `flowlines` is provided.

- only_split_candidates:

  Logical (default `FALSE`). When `TRUE` and `split_events` has a
  `split_candidate` column, restrict to rows where
  `split_candidate == TRUE`. Useful if only confirmed split-interior
  events should force protection.

## Value

Named character vector of unique segment IDs to exclude from collapse.

## Examples

``` r
if (FALSE) { # \dontrun{
excl <- events_to_exclude_cats(events    = poi_events,
                               flowlines = flines,
                               protect_downstream = TRUE)
refactor_flowlines(flines, split_flines_meters = 2000,
                   collapse_flines_meters = 1000,
                   out_refactored = "refactored.gpkg",
                   out_reconciled = "reconciled.gpkg",
                   exclude_cats = excl)
} # }
```
