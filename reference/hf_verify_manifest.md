# Verify inputs against a recorded manifest

Re-fingerprints the files a manifest names and reports each as `ok`,
`changed` (the hash differs) or `missing`. Use it before a rebuild to
establish that the inputs are the ones a previous result came from.

## Usage

``` r
hf_verify_manifest(manifest, strict = FALSE)
```

## Arguments

- manifest:

  Path to a manifest written by
  [`hf_input_manifest()`](https://lynker-spatial.github.io/hydrofabric/reference/hf_input_manifest.md),
  or the data frame it returned.

- strict:

  If `TRUE`, any file that changed or went missing raises an error.
  Default `FALSE` reports and returns.

## Value

A data frame with `role`, `path`, `status`, and the recorded and current
hashes. Returned invisibly when everything matches.

## Examples

``` r
f <- tempfile(); writeLines("x", f)
m <- hf_input_manifest(reference = f)
hf_verify_manifest(m)
#> [manifest] ok       reference        file1f92498166d
#> [manifest] all 1 input(s) match
```
