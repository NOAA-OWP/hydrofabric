# Find reference reaches whose downstream neighbour is missing

A reach with no downstream is an outlet: it reaches the sea, or a closed
basin, and the network legitimately stops. A reach that *names* a
downstream the file does not contain is something else. The network was
truncated, and the reach now reads as an outlet to anything that
consumes it.

Routing is where it shows. Every reach with no downstream becomes a
terminal nexus, and flow arriving there is not routed on, so a truncated
network underreports flow at each cut. The pipeline cannot repair it:
there is no downstream reach to route into.

The reference carries the evidence itself, in the value-added attributes
it inherits. `flowpath_toid` says where a reach drains and `dnhydroseq`
names the same neighbour by its `hydroseq`. A true outlet has neither. A
truncated reach still carries `dnhydroseq`, pointing at a reach that was
removed, and that dangling pointer is what this reports.

Observed on the CONUS reference at the Canadian border: of VPU 04's
2,999 terminal nexuses, 2,948 trace to a reach naming a downstream by
`dnhydroseq` that exists nowhere in the file. VPU 01, whose terminals
are mostly the New England coastline, has 9. The `terminalfl` column
cannot be used to tell the two apart, because it is unpopulated for
`NWM_v3_0_CONUS_Prelim_ls_mod`: 14 reaches of that source's 109,850
carry the flag.

Report-only, and deliberately not wired into
[`run_pipeline()`](https://lynker-spatial.github.io/hydrofabric/reference/run_pipeline.md).
The coverage belongs to the reference and is not ours to extend, and
this asks a question about a reference release rather than about a VPU.

## Usage

``` r
check_reference_downstream_resolves(
  reference_gpkg,
  layer = "flowpaths",
  id = "flowpath_id",
  toid = "flowpath_toid",
  dn = "dnhydroseq",
  seq = "hydroseq",
  vpuid = NA_character_,
  strict = FALSE
)
```

## Arguments

- reference_gpkg:

  Path to the reference GeoPackage.

- layer:

  Flowline layer name.

- id, toid, dn, seq:

  Column names for the identifier, the downstream identifier, the
  downstream hydrologic sequence, and the hydrologic sequence.

- vpuid:

  Optional VPU to restrict the *reported* reaches to. The search for the
  named downstream always spans the whole file, since a reach may
  legitimately drain into the next VPU.

- strict:

  If `TRUE`, warn when any are found.

## Value

Invisibly, a data frame of `id`, `dangling` (which pointer dangles) and
`names_id`, the neighbour that is missing.
