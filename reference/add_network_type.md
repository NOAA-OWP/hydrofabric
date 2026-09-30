# Annotate Flowpaths and Divides with Type Info

Adds boolean flags `has_divide` to flowpaths and `has_flowpath` to
divides.

## Usage

``` r
add_network_type(network_list)
```

## Arguments

- network_list:

  List of `flowpaths` and `divides` (sf objects).

## Value

Updated `network_list` with type info columns.

## Examples

``` r
if (FALSE) { # \dontrun{
nl <- list(flowpaths = sf::st_read("out_vpu01.gpkg", "flowpaths"),
           divides   = sf::st_read("out_vpu01.gpkg", "divides"))
nl <- add_network_type(nl)
} # }
```
