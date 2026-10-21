# Build a set of VPUs, then merge them

A loop over
[`run_pipeline()`](https://NOAA-OWP.github.io/hydrofabric/reference/run_pipeline.md)
with per-VPU inputs, so a national build is a list of paths rather than
a script. One VPU failing does not stop the rest by default; the return
value says which succeeded, and
[`merge_vpus()`](https://NOAA-OWP.github.io/hydrofabric/reference/merge_vpus.md)
can be run from it once they all have.

`inputs` is a named list, one entry per `vpuid`, each holding the
arguments that differ for that VPU. Anything passed in `...` applies to
every VPU, so shared paths are given once.

## Usage

``` r
run_vpus(
  inputs,
  out_dir,
  ...,
  continue_on_error = TRUE,
  merge = FALSE,
  merge_dir = NULL
)
```

## Arguments

- inputs:

  Named list of per-VPU argument lists, keyed by `vpuid`.

- out_dir:

  Parent output directory; each VPU writes `vpu_{vpuid}/`.

- ...:

  Arguments common to every VPU, forwarded to
  [`run_pipeline()`](https://NOAA-OWP.github.io/hydrofabric/reference/run_pipeline.md).

- continue_on_error:

  If `TRUE` (default), record a failure and carry on. `FALSE` stops at
  the first one.

- merge:

  If `TRUE`, merge the VPUs that succeeded once the loop finishes.
  Default `FALSE`.

- merge_dir:

  Directory for the merged domain products when `merge = TRUE`. Defaults
  to `file.path(out_dir, "global")`.

## Value

A data frame with one row per VPU: `vpuid`, `status`
(`"ok"`/`"error"`/`"invariant"`/`"skipped"`), `ngen_gpkg`, `minutes`,
and `message`. `"invariant"` means the build finished but its gates went
red; it counts as a failure and blocks the merge, because a product that
does not pass its own suite is not one to publish. Returned invisibly
when every VPU succeeded.

## Examples

``` r
if (FALSE) { # \dontrun{
conus <- c("01","02","03N","03S","03W","04","05","06","07","08","09",
           "10L","10U","11","12","13","14","15","16","17","18")

# Shared inputs once, per-VPU overrides where they differ.
inputs <- stats::setNames(vector("list", length(conus)), conus)
inputs[["20"]] <- list(reference_gpkg = "hi_reference.gpkg",
                       fdr_path = "hi/hi_fdr.vrt", fac_path = "hi/hi_fac.vrt")

res <- run_vpus(
  inputs,
  out_dir        = "out",
  reference_gpkg = "conus_reference.gpkg",
  fdr_path       = "conus/fdr.vrt",
  fac_path       = "conus/fac.vrt",
  skip_if_exists = TRUE,
  merge          = TRUE)

subset(res, status != "ok")
} # }
```
