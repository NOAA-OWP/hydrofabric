# Reconcile divides to refactored flowlines (CRS-aware, lr-capable)

Raster-aware divide reconciliation with two key properties:

\(1\) CRS from rasters: all metric spatial operations (gap fill, area
calculations, MULTIPOLYGON cleanup) use the native CRS of the FDR raster
instead of a hardcoded EPSG:5070. Required for datasets outside CONUS
(e.g. PRVI, Alaska).

\(2\) lr flag: when lr = TRUE, each reconciled divide polygon is further
split left/right along its associated flowline using lwgeom::st_split,
matching the reference split_catchment_divide(lr = TRUE) approach.
Output flowpath_id values are suffixed "\_L" and "\_R".

## Usage

``` r
reconcile_divides(
  flowline,
  divides,
  fdr_path,
  fac_path,
  fallback_fdr_path = NULL,
  fallback_fac_path = NULL,
  lookup_table = NULL,
  lr = FALSE,
  lr_densify_m = 100,
  lr_snap_thresh_m = 100,
  parallel = TRUE,
  num_workers = NULL
)
```

## Arguments

- flowline:

  sf LINESTRING. Refactored flowlines with flowpath_id and flowpath_toid
  columns. Typically out\$refactored from refactor_flowlines().

- divides:

  sf POLYGON. Input divide polygons with flowpath_id.

- fdr_path:

  Character. Path to D8 flow-direction raster (required).

- fac_path:

  Character. Path to flow-accumulation raster (required).

- fallback_fdr_path, fallback_fac_path:

  Optional paths to fallback FDR / FAC rasters, used when the primary
  rasters return NoData at a split point. `NULL` disables the fallback
  (default); ignored when equal to the primary paths.

- lookup_table:

  The out\$reconciled sf from refactor_flowlines(), or NULL. When
  provided (must have ID and member_flowpath_id columns) it is used for
  reconciled-aware deduplication before splitting.

- lr:

  Logical. If TRUE, apply left/right bank split along the flowline after
  D8 reconciliation. Requires lwgeom. Default FALSE.

- lr_densify_m:

  Numeric. dfMaxLength for st_segmentize when densifying divide boundary
  before lr split. Default 100. Units match raster CRS.

- lr_snap_thresh_m:

  Numeric. Distance threshold (raster CRS units) for snapping flowline
  endpoints to divide boundary. Default 100.

- parallel:

  Logical. Enable parallel multi-parent divide splitting. Default TRUE.

- num_workers:

  Integer or NULL. Parallel workers. NULL uses 4.

## Value

sf POLYGON with columns divide_id, flowpath_id, and geometry. When lr =
TRUE, flowpath_id values are suffixed "\_L" / "\_R".

## See also

refactor_flowlines()

## Examples

``` r
if (FALSE) { # \dontrun{
out <- refactor_flowlines(flines, split_flines_meters = 2000,
                          collapse_flines_meters = 1000,
                          out_refactored = "refactored.gpkg",
                          out_reconciled = "reconciled.gpkg")
rec_div <- reconcile_divides(flowline      = out$refactored,
                             divides       = divides_sf,
                             fdr_path      = "fdr.tif",
                             fac_path      = "fac.tif",
                             lookup_table  = out$reconciled)
} # }
```
