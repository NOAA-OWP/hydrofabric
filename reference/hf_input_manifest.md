# Record the inputs a build read

Fingerprints each supplied file, size, modification time and MD5, and
optionally writes the result as JSON. Pair it with
[`hf_verify_manifest()`](https://NOAA-OWP.github.io/hydrofabric/reference/hf_verify_manifest.md)
to confirm a later run is reading the same inputs.

Hashing streams the whole file, so a multi-gigabyte reference takes a
little time. It is worth it once per build rather than per VPU.

## Usage

``` r
hf_input_manifest(..., out = NULL)
```

## Arguments

- ...:

  Named file paths, e.g. `reference = "...", fdr = "...", fac = "..."`.
  Names become the manifest's roles. `NULL` entries are dropped.

- out:

  Optional path to write the manifest as JSON. `NULL` (default) returns
  it without writing.

## Value

A data frame with one row per input: `role`, `path`, `exists`, `bytes`,
`mtime`, `md5`.

## Examples

``` r
f <- tempfile(); writeLines("x", f)
hf_input_manifest(reference = f)
#>        role                             path exists bytes                mtime
#> 1 reference /tmp/RtmpXwoem0/file1fed4dd6525a   TRUE     2 2026-11-12T13:14:57Z
#>                                md5
#> 1 401b30e3b8b5d629635a5c613cdb7919
```
