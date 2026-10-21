# Record what a build produced

Writes a manifest of finished products: the md5, size, layer counts and
invariant verdict of each, together with the commit that built them and
the digest of the input manifest they were built from. Verifiable with
[`hf_verify_manifest()`](https://NOAA-OWP.github.io/hydrofabric/reference/hf_verify_manifest.md),
the same way inputs are.

## Usage

``` r
hf_output_manifest(
  gpkgs,
  out,
  input_manifest = NULL,
  git_sha = Sys.getenv("HF_GIT_SHA", "")
)
```

## Arguments

- gpkgs:

  Paths to the finished products.

- out:

  Path to write the manifest JSON to.

- input_manifest:

  Optional path to the input manifest used for the build; its md5 is
  recorded so the two can be tied together.

- git_sha:

  Commit that produced the products. Defaults to `HF_GIT_SHA`.

## Value

Invisibly, the manifest as a list. Written to `out` as JSON.

## Details

This is the other half of
[`hf_input_manifest()`](https://NOAA-OWP.github.io/hydrofabric/reference/hf_input_manifest.md).
Without it a product cannot be checked against what was shipped, and the
tie between a delivered file, the code that made it and the inputs it
consumed lives only in a build log.

The input manifest is recorded by digest rather than copied, so one line
establishes which input set was used without duplicating it.

## Examples

``` r
if (FALSE) { # \dontrun{
hf_output_manifest(
  Sys.glob("global/*_ngen.gpkg"),
  out = "global/output_manifest.json",
  input_manifest = "owp-inputs/national_inputs_manifest.json")
} # }
```
