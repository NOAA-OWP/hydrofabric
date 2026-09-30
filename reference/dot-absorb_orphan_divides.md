# Absorb an unowned divide into its neighbour instead of deleting it

When a divide's flowpath is swapped away or cleared, that divide is left
unowned. Deleting it removes its polygon from the layer, and with it
that land's coverage: on a single CONUS VPU this cost 97 km2, invisible
because `divides_area_conserved` compares sums at a 0.5% tolerance and a
hole and an overlap removal look alike in a sum.

## Usage

``` r
.absorb_orphan_divides(div, orphans_mask)
```

## Arguments

- div:

  Divides `sf`.

- orphans_mask:

  Logical vector over `div` rows marking unowned divides.

## Value

A list with `div` (the repaired layer), `n_absorbed` and `n_kept`
(orphans that found no neighbour and were retained).

## Details

Unowned is not the same as empty. The land still drains somewhere, so
the polygon is merged into the neighbouring divide it shares the most
boundary with, which conserves coverage while keeping one divide per
flowpath. A divide with no neighbour at all is kept rather than dropped,
since losing area is worse than carrying an unowned polygon forward.
