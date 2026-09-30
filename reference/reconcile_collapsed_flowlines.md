# Reconcile collapsed flowlines

Equivalent to the reference reconcile-collapsed-flowlines step but
data.table-based for speed. Returns sf or data.frame depending on
whether `geom` is supplied.

## Usage

``` r
reconcile_collapsed_flowlines(flines, geom = NULL, id = "flowpath_id")
```

## Arguments

- flines:

  data.frame output of the collapse step (`.reactor_collapse`).

- geom:

  sf data.frame with geometry column and id column for line union.

- id:

  character name of the ID column in geom (default "flowpath_id").

## Value

sf data.frame (if geom supplied) or data.frame with columns: ID, toID,
lengthkm, member_flowpath_id (list), levelpathid, hydroseq.

## Examples

``` r
if (FALSE) { # \dontrun{
collapsed <- .reactor_collapse(flines)
rec <- reconcile_collapsed_flowlines(collapsed, geom = flowline_sf,
                                     id = "flowpath_id")
} # }
```
