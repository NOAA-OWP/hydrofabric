# Verify products against an output manifest

Verify products against an output manifest

## Usage

``` r
hf_verify_output_manifest(manifest, strict = FALSE)
```

## Arguments

- manifest:

  Path to a manifest written by
  [`hf_output_manifest()`](https://lynker-spatial.github.io/hydrofabric/reference/hf_output_manifest.md).

- strict:

  If `TRUE`, raise an error when a product is missing or altered.

## Value

Invisibly, a data frame with `name`, `exists`, `md5_ok`.

## Examples

``` r
if (FALSE) { # \dontrun{
hf_verify_output_manifest("global/output_manifest.json", strict = TRUE)
} # }
```
