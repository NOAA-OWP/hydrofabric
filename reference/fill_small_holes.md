# Fill tiny interior holes in polygons

Removes interior rings (holes) whose area is below a threshold. Common
fix for raster-to-vector artifacts where polygon unioning leaves
sub-pixel slivers as holes.

## Usage

``` r
fill_small_holes(x, threshold_m2 = 1000, layer = NULL, inplace = TRUE)
```

## Arguments

- x:

  An `sf` object with POLYGON or MULTIPOLYGON geometry, or a file path
  to a GeoPackage layer.

- threshold_m2:

  Numeric. Interior rings with area below this threshold (in square
  meters, using the object's projected CRS) are removed. Default `1000`
  (0.001 sqkm).

- layer:

  If `x` is a file path, the layer name to read/write.

- inplace:

  If `x` is a file path, whether to write the filled result back to the
  same layer. Default `TRUE`.

## Value

An `sf` object with tiny holes removed. If `inplace = TRUE` and `x` is a
file path, writes to disk and returns invisibly.
