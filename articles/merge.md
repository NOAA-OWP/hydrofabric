# Merge

[`merge_vpus()`](https://NOAA-OWP.github.io/hydrofabric/reference/merge_vpus.md)
is the last phase: it combines the per-VPU NextGen products into domain
products with globally unique identifiers, and resolves the topology
that crosses VPU boundaries. It reads only finished `ngen_vpu*.gpkg`
files, so it needs neither the reference nor the rasters.

Every VPU must have passed its invariants first. Merge all VPUs in one
call: the merge assigns one id pool across them, so merging a subset
assigns different ids.

## Inputs and outputs

**Input**: `gpkg_map`, a named list mapping each `vpuid` to its
`ngen_vpu*.gpkg`.

**Output**, written to `out_dir`:

- `superconus_ngen.gpkg` (VPUs 01-18)
- `ak_ngen.gpkg` (VPU 19)
- `hi_ngen.gpkg` (VPU 20)
- `prvi_ngen.gpkg` (VPU 21)
- `guam_ngen.gpkg`, only when VPU 22 is supplied

## What it does

1.  **Assign global ids.** Contiguous integers across every domain,
    carrying the `fp-`, `cat-`/`lnd-`, `nex-` and `tnx-` prefixes. The
    pool is assigned in VPU order, so no two domains can collide.
2.  **Fix cross-VPU topology.** A flowpath whose downstream lies in
    another VPU is rewired to the global id of that reach.
3.  **Recompute domain-scale attributes.** Hydrosequence and total
    drainage area are accumulated over the merged network rather than
    carried from the per-VPU inputs.
4.  **Add the upstream index** to each merged product, so `upstream_id`
    and `num_upstreams` describe the domain network rather than a single
    VPU.

## Running it

``` r

library(hydrofabric)

vpus <- c("01","02","03N","03S","03W","04","05","06","07","08","09",
          "10L","10U","11","12","13","14","15","16","17","18","19","20","21")

merge_vpus(
  gpkg_map = as.list(stats::setNames(
    file.path("out", paste0("vpu_", vpus), sprintf("ngen_vpu%s.gpkg", vpus)), vpus)),
  out_dir  = "out/global"
)
```

## Parameters

| Parameter | Default | Controls |
|----|----|----|
| `mods_csv` | bundled `extdata/rpu_vpu_out.csv` | Cross-VPU connectivity edits applied during the merge |
| `overwrite` | `FALSE` | Replace existing domain products |
| `strict` | `TRUE` | Abort a domain whose checks fail, rather than writing it |
| `grouping` | `c("vpu", "tdx_continent")` | How VPUs are grouped into domains |
| `target_crs` | `"EPSG:4326"` | CRS of the merged products |

The merged products are written in `target_crs` (EPSG:4326 by default),
not the per-VPU working CRS.

## Validation

Each merged domain is checked before it is written: structural checks on
the merged tables (id uniqueness, prefix integrity, `toid` resolution,
divide ownership, an acyclic graph), the conservation and completeness
checks in
[`hfutils::hf_check_merge_invariants()`](https://lynker-spatial.github.io/hfutils/reference/hf_check_merge_invariants.html)
(every VPU present, feature counts, divide area and drainage area
conserved from the per-VPU inputs), and the ngen suite run over the
re-identified network (flowline id uniqueness, `flowline_toid`
resolution, the flowline DAG, and the POI contract).

A domain that fails under `strict = TRUE` is not written.

## Next

The field-level schema is in `docs/data_model.md` in the source
repository, and [the upstream
index](https://NOAA-OWP.github.io/hydrofabric/articles/upstream-index.md)
covers the range queries each domain supports.
