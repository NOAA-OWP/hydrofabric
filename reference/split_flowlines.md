# Split flowlines (sequential, flowpath_id-based)

Split flowlines by events and/or a max segment length. Returns semantic
IDs with .part suffix ordered upstream (.1) -\> downstream (.N).

## Usage

``` r
split_flowlines(
  flines,
  max_length = NULL,
  events = NULL,
  avoid = NA,
  min_split_m = 10
)
```

## Arguments

- flines:

  sf LINESTRINGs with columns:

  - flowpath_id (coerced to character)

  - one downstream column: flowpath_toid \| toCOMID \| toid (kept as-is)

  - LENGTHKM (will be recomputed from geometry) CRS should be projected
    in meters for length-based splitting.

- max_length:

  numeric (meters). If NULL, only event splits apply.

- events:

  data.frame of event matches with **per-event measures**: required:
  flowpath_id, REACHCODE, FromMeas, ToMeas, REACH_meas optional:
  identifier (passed through to output)

- avoid:

  vector of flowpath_id to exclude from *length* splitting (events still
  apply).

- min_split_m:

  numeric (meters). Minimum piece length: an event split that would
  create a fragment shorter than this is dropped. Default 10.

## Value

sf with split flowlines; includes event\_\* fields and recomputed
LENGTHKM.

## Examples

``` r
if (FALSE) { # \dontrun{
flines <- sf::st_read("refactored.gpkg", "refactored")
segs <- split_flowlines(flines, max_length = 2000)
} # }
```
