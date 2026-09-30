# Attach NextGen IDs (`fp-N`, `cat-N`) to a lineage table

`nextgen_from_gpkg` prefixes aggregated IDs with `"fp-"` and `"cat-"`
deterministically, so this attach is a trivial string paste, no join to
the network table required.

## Usage

``` r
lineage_attach_ngen(lineage, flowpath_prefix = "fp-", divide_prefix = "cat-")
```

## Arguments

- lineage:

  A lineage `data.table` that already has `aggregated_id`.

- flowpath_prefix:

  Prefix for final flowpath IDs. Default `"fp-"`.

- divide_prefix:

  Prefix for final divide IDs. Default `"cat-"`.

## Value

Lineage with `final_flowpath_id` and `final_divide_id` populated.

## Examples

``` r
if (FALSE) { # \dontrun{
lin <- lineage_attach_aggregated(lin, network)
lin <- lineage_attach_ngen(lin)
} # }
```
