# Refactor reference flowlines into reconciled reaches

Refactors an NHDPlus flowline network by splitting long segments,
collapsing short inter-confluence segments, and reconciling collapsed
groups into a clean topology. Combines the reference accuracy guards
with data.table O(n) performance. Two corrections over
fast_collapse_flowlines: (1) a problem-headwater guard that marks
confluence-adjacent headwaters as uncollapsible ("-9999") rather than
forcing them downstream; and (2) a Stage 4/5 threshold split that
applies mainstem_thresh to single-upstream (mainstem-run) segments and
thresh to multi-upstream (confluence) segments.

## Usage

``` r
refactor_flowlines(
  flines,
  split_flines_meters,
  collapse_flines_meters,
  collapse_flines_main_meters = NULL,
  out_refactored,
  out_reconciled,
  three_pass = TRUE,
  purge_non_dendritic = TRUE,
  exclude_cats = NULL,
  events = NULL,
  cores = 1L,
  warn = FALSE
)
```

## Arguments

- flines:

  sf data.frame with columns flowpath_id, flowpath_toid, lengthkm,
  hydroseq, levelpathid. totdasqkm optional.

- split_flines_meters:

  numeric max segment length (meters).

- collapse_flines_meters:

  numeric minimum inter-confluence segment length (meters).

- collapse_flines_main_meters:

  numeric minimum mainstem segment length (meters). Defaults to
  collapse_flines_meters.

- out_refactored:

  character path for split+collapsed gpkg output.

- out_reconciled:

  character path for reconciled gpkg output.

- three_pass:

  logical run collapse at 10%, 50%, 100% of threshold.

- purge_non_dendritic:

  logical (unused, kept for API compatibility).

- exclude_cats:

  integer/character IDs to protect from collapse.

- events:

  sf point event locations (passed to split_flowlines). Their geometry
  also determines which reach outlets a collapse may not merge away.

- cores:

  integer cores for parallel split (uses future/furrr when \>1).

- warn:

  logical print progress messages.

## Value

invisibly returns list(refactored, reconciled) as sf data.frames.

## Examples

``` r
if (FALSE) { # \dontrun{
flines <- sf::st_read("reference_01.gpkg", "flowpaths")
out <- refactor_flowlines(flines,
                          split_flines_meters    = 2000,
                          collapse_flines_meters = 1000,
                          out_refactored = "refactored.gpkg",
                          out_reconciled = "reconciled.gpkg")
} # }
```
