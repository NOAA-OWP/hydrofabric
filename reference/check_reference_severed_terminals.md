# Find reference reaches that terminate on a channel that carries on

A reach with no downstream becomes a terminal nexus, and flow arriving
there is not routed on. That is right for an outlet: a coast, a closed
basin. It is wrong where the channel plainly continues, and a routing
consumer sees the flow simply stop.

The test is geometric because the attributes cannot answer it.
`terminalfl` is unpopulated for whole sources, and `dnhydroseq` is
stored too narrowly to hold the values it points at. What does
distinguish the two is what lies at the reach's downstream end: nothing,
or another reach carrying more water than the one that just stopped.

Observed on the CONUS reference: of 15,877 terminal reaches, 822 end
within a DEM cell of another reach, and 446 of those end on one that
drains more, including a 324,094 km2 river in VPU 13 that stops 30 m
short of its own continuation and the Saint John at Fort Kent, where
12,334 km2 stops 13 m from a reach carrying 12,363 km2. They are spread
nationwide rather than gathered at the borders: VPU 17 has 76, VPU 18
69, VPU 16 68.

Report-only, and deliberately not wired into
[`run_pipeline()`](https://lynker-spatial.github.io/hydrofabric/reference/run_pipeline.md).
The topology belongs to the reference, and repairing it means writing
edges into an input rather than working around it in a build. This
reports what a staged reference still carries.

## Usage

``` r
check_reference_severed_terminals(
  reference_gpkg,
  layer = "flowpaths",
  id = "flowpath_id",
  toid = "flowpath_toid",
  area = "totdasqkm",
  vpuid = NA_character_,
  tolerance_m = 30,
  strict = FALSE
)
```

## Arguments

- reference_gpkg:

  Path to the reference GeoPackage.

- layer:

  Flowline layer name.

- id, toid:

  Identifier and downstream-identifier columns.

- area:

  Column holding total upstream drainage area, used to tell a receiver
  from a neighbouring headwater.

- vpuid:

  Optional VPU to restrict to, when the layer carries `vpuid`.

- tolerance_m:

  Distance within which a reach counts as ending on another. The default
  is the 30 m DEM cell. A wider radius reports unrelated channel heads:
  at 100 m VPU 01 scores 17 and at 5 m it scores none.

- strict:

  If `TRUE`, warn when any are found.

## Value

Invisibly, a data frame of `id`, `area`, `continues_into` and
`into_area`, largest first.
