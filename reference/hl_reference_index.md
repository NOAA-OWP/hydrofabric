# Index a reference for hydrolocation checking

Reads the `flowpath_id` and `vpuid` values a reference carries, in the
shape
[`check_hydrolocation_input()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_hydrolocation_input.md)
expects. Attribute-only, so it does not read geometry.

## Usage

``` r
hl_reference_index(reference_gpkg, layer = "flowpaths")
```

## Arguments

- reference_gpkg:

  Path to a reference fabric.

- layer:

  Flowpaths layer name. Default `"flowpaths"`.

## Value

A list with `ids` and `vpuids`.

## Examples

``` r
if (FALSE) { # \dontrun{
hl_reference_index("conus_reference.gpkg")
} # }
```
