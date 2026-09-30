# Version guard for stage checkpoint reuse

`run_pipeline(use_checkpoints = TRUE)` resumes a build from existing
stage checkpoints (e.g. after an OOM/crash) instead of recomputing from
scratch. That is only safe when the checkpoints were produced by the
*same* pipeline build and params: otherwise a code fix (reinstall) or a
param change would silently reuse a stale/buggy intermediate. These
helpers fingerprint the current build+params and gate reuse on a match.

## Usage

``` r
.pipeline_fingerprint(params = NULL)

.checkpoints_valid(work_dir, fp)

.stamp_fingerprint(work_dir, fp)
```

## Arguments

- params:

  Optional list of stage params to fold into the fingerprint.

- work_dir:

  Stage work directory holding the `.pipeline_fingerprint` stamp.

- fp:

  Current fingerprint from `.pipeline_fingerprint()`.

## Value

A single character string identifying the build + params.

`.checkpoints_valid`: `TRUE` iff `work_dir` was stamped with `fp`.

`.stamp_fingerprint`: invisibly `NULL`; writes the stamp file.

## Details

The fingerprint is intentionally coarse: any mismatch (a reinstall flips
the `Built` stamp; a param edit flips the params component) forces a
clean recompute, which is always correct. Over-invalidation only costs
time.
