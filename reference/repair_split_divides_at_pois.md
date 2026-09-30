# Repair fp\<-\>cat joins where PRD missed a POI-triggered divide split

After `resolve_disjoint_flowpath_divides` clears bogus fp\<-\>cat joins,
a small class of legitimate flowpaths is left with `has_divide = FALSE`
even though their LINE sits entirely inside the downstream neighbour's
cat on the same NHD base. Pattern: upstream half (`.1`) of a split NHD,
with a POI at its outlet: PRD should have carved a new divide boundary
at that POI but the full divide ended up with the downstream half.

This pass re-runs the raster-based splitter (`split_divides_d8`) for
just those cases, using the POI location (= fp_up outlet = fp_down
inlet) as the split point. Upstream piece becomes fp_up's new cat;
downstream piece replaces fp_down's existing cat.

## Usage

``` r
repair_split_divides_at_pois(
  input_gpkg,
  flowpath_layer = "flowpaths",
  divides_layer = "divides",
  pois_layer = "pois",
  network_layer = "network",
  fdr_path,
  fac_path,
  fallback_fdr_path = NULL,
  fallback_fac_path = NULL,
  coverage_threshold = 0.9,
  num_workers = 4L
)
```

## Arguments

- input_gpkg:

  Path to the pre-ngen GeoPackage (layers: flowpaths, divides, pois;
  optional: network).

- flowpath_layer:

  Name of the flowpaths layer in `input_gpkg`. Default `"flowpaths"`.

- divides_layer:

  Name of the divides layer in `input_gpkg`. Default `"divides"`.

- pois_layer:

  Name of the POIs layer in `input_gpkg`. Default `"pois"`.

- network_layer:

  Name of the network layer in `input_gpkg`. Default `"network"`.

- fdr_path, fac_path:

  Paths to FDR / FAC rasters (same files PRD used earlier in the
  pipeline).

- fallback_fdr_path, fallback_fac_path:

  Optional paths to a fallback FDR/FAC pair. When a split attempt
  returns NULL with the primary rasters (split point lands on NoData),
  the splitter is retried with the fallback. `NULL` disables the retry
  (default). Ignored when equal to the primary paths (path equality
  check, same as `add_flowlines`).

- coverage_threshold:

  Minimum fraction of fp_up's line that must sit inside fp_down's cat
  for the pair to qualify. Default 0.9.

- num_workers:

  Parallel workers for the splitter. Default 4.

## Value

Invisibly, a data.frame of repaired
`{fp_up, fp_down, new_cat, area_km2}` rows.
