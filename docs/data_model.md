# hydrofabric Output Data Model

> **Scope:** The NextGen GeoPackage schema written by `nextgen_fabric()`. Mapped against the upstream
> [Lynker Hydrofabric Data Model](https://github.com/lynker-spatial/hfdatamodel/blob/main/docs/data_model.md).

This document enumerates every column each layer emits and how it maps to the
upstream spec.

Terminology:
- **Spec**: the Lynker Hydrofabric Data Model linked above.
- **Emitted**: the column name in the GeoPackage the pipeline produces.
- **Status**: ✓ matches, ◆ renamed, ✖ not emitted, + emitted beyond the spec.

---

## Relationship to HY_Features (WaterML2 Part 3)

The Lynker Hydrofabric Data Model is an **extension of the HY_Features (OGC
WaterML2 Part 3) conceptual model**, specialized with the concrete feature
types and attributes hydrologic and hydraulic modelers need. HY_Features
supplies the abstract catchment / hydro-network concepts; the model below
realizes them as routing-ready features and adds the finer flowline grain and
model attributes that NextGen and T-Route consume.

| HY_Features concept | This data model |
|---|---|
| `HY_Catchment` | **divides**: drainage-area polygons (`divide_id`) |
| `HY_FlowPath` | **flowpaths**: aggregated routing linestrings (`flowpath_id`) |
| `HY_HydroNexus` | **nexus**: flow-handoff nodes between flowpaths (`nexus_id`) |
| catchment realization, finer grain | **flowlines**: reconciled reaches carried for area-weighted post-processing (an extension beyond core HY_Features) |
| `HY_HydroLocation` | **hydrolocations** / **pois**: gages, assimilation and routing points (`hl_reference`, `poi_id`) |
| referenced waterbody | **lakes / waterbodies**: lake membership flagged on flowpaths |

The network is a directed acyclic graph of `flowpath to nexus to flowpath` edges,
isolated only at known terminal sinks (ocean / endorheic basins).

## Design targets

The build targets the geometric and attribute characteristics needed for
continental hydrologic and hydraulic modeling:

- **Divide area**: divides are aggregated toward a working band with a 3 km²
  floor (`min_area_sqkm`, `aggregate_fabric()`). The floor is the enforced
  parameter; the upper end of the band is what the collapse works toward, not a
  hard bound.
- **Routing-unit (flowline) length**: flowlines carry a 250 m minimum-length
  floor (`collapse_flines_meters`) with breaks forced at hydrologic points of
  interest and at junctions; longer segments are subdivided by
  `add_flowlines(split_flines_meters = …)` (default 2000 m; lower it to tighten
  the routing grain). Where a junction or POI prevents a target length, flow is
  routed to the nearest upstream/downstream node.
- **Gage points**: gage assimilation points are represented as hydrologic
  points of interest unless within the 250 m snap threshold of an existing POI.
- **Topology**: a directed acyclic graph suitable for hydraulic routing with
  T-Route, isolated only at known sinks (enforced by the invariant suite).
- **Attribute sufficiency**: emitted attributes support hydraulic routing with
  T-Route and area-weighted post-processing that derives streamflow at flowline
  outlets (drainage area, length, slope, channel geometry; see the layer tables
  below).

---

## Layers

| Spec entity | Emitted layer | Notes |
|---|---|---|
| Flowpaths | `flowpaths` | ✓ |
| Divides | `divides` | Includes both network and landscape divides (`has_flowline` flag) |
| (none) | `landscape_divides` | Pre-ngen only; merged into `divides` post-ngen |
| POIs | `pois` | ✓ |
| Hydrolocations | `hydrolocations` | ✓ Dedicated layer emitted by `finalize_data_model`; see below |
| Network | `network` | Aspatial attribute table (one row per flowline × reference member). Schema diverges from spec; see below. |
| Flowlines | `flowlines` | Reconciled reaches |
| Incremental Area | *(none)* | ✖ (only `flowlines.incremental_area_id` pointer exists) |
| Nexus | `nexus` | ✓ |
| Referenced waterbody | `lakes` | + Waterbody / reservoir parameters (points in CONUS, polygons elsewhere); see below |
| (none) | `flowpath-attributes` | + Routing channel geometry, roughness and Muskingum params, keyed on `flowpath_id` (aspatial); see below |
| (none) | `flowline-attributes` | + Same routing params keyed on `flowline_id` (aspatial); see below |

---

## Flowpaths

| Spec column | Emitted | Status | Notes |
|---|---|---|---|
| `flowpath_id` | `flowpath_id` | ✓ | Prefixed `fp-` |
| `flowpath_toid` | `flowpath_toid` | ✓ | Points to a `nexus_id`, not a flowpath; the spec allows a self-reference, this model routes via nexus |
| `vpuid` | `vpuid` | ✓ | |
| `poi_id` | `poi_id` | ✓ | |
| `mainstem_id` | `mainstem_id` | ✓ | The reference's level-path id, the same in every VPU |
| `areasqkm` | `areasqkm` | ✓ | |
| `lengthkm` | `lengthkm` | ✓ | |
| `total_dasqkm` | `total_dasqkm` | ✓ | |
| `hydroseq` | `hydroseq` | ✓ | |
| `streamorder` | `streamorder` | ✓ | Strahler order via `hfutils::get_streamorder` on the flattened fp to nexus to fp graph |
| `slope` | `slope` | ✓ | Length-weighted mean across member reference ids; sentinel 1e-5 excluded when ≥1 non-sentinel member survives; NA when no valid member exists (or reference input lacks slope) |
| (none) | `rec_id_list` | + | Comma-separated list of member reconciled flowlines |
| (none) | `member_comid` | + | Comma-separated Lynker Reference Fabric ids (source lineage), matched to NHD COMIDs where a correspondence exists |
| (none) | `has_divide` | + | Boolean: TRUE for network, FALSE for connector/artificial paths |
| (none) | `divide_id` | + | FK back to the paired divide; not in the spec's Flowpath table |
| (none) | `lake_id` | + | FK to `lakes.lake_id` when the flowpath is a lake member; `NA` otherwise |
| (none) | `upstream_id` | + | Nested-set pre-order index; upstream-of-X is `upstream_id` in `(u, u + n]`. See `vignette("upstream-index")` |
| (none) | `num_upstreams` | + | Size of that flowpath's upstream subtree |
| (none) | `ibt` | + | 1 where the merge connected the flowpath to a downstream flowpath in another VPU, else 0 |

---

## Divides

| Spec column | Emitted | Status | Notes |
|---|---|---|---|
| `divide_id` | `divide_id` | ✓ | Prefixed `cat-` (network) or `lnd-` (landscape) |
| `flowpath_id` | `flowpath_id` | ✓ | |
| `flowpath_toid` | `flowpath_toid` | ✓ | Points to the downstream `nexus_id`; denormalized from the paired flowpath |
| `vpuid` | `vpuid` | ✓ | |
| `areasqkm` | `areasqkm` | ✓ | |
| `type` | `type` | ✓ | Values: `network`, `landscape`. This binary split reflects the three-way feature model (network with divide, network without divide, landscape divide with no network) rather than the spec's `network`/`coastal`/`internal` |
| (none) | `has_divide` | + | Boolean: TRUE when the divide has a paired network flowpath |
| (none) | `has_flowline` | + | Boolean mirror of `type == "network"`; redundant but useful |
| (none) | `ds_id` | + | Downstream divide id (used by some consumers); not in spec |
| (none) | `small_reason` | + | Diagnostic string flagging why a divide stayed below the target area band (`NA` when not undersized) |
| (none) | `upstream_id` | + | Nested-set pre-order index, denormalized from the owning flowpath by `flowpath_id`; see `vignette("upstream-index")` |
| (none) | `num_upstreams` | + | Size of that flowpath's upstream subtree, denormalized alongside `upstream_id` |
| (none) | `orig_divide_id` | + | Not populated in this release |

---

## POIs

| Spec column | Emitted | Status | Notes |
|---|---|---|---|
| `poi_id` | `poi_id` | ✓ | |
| `flowpath_id` | `flowpath_id` | ✓ | Outlet convention: POI is on the flowpath whose OUTLET is the POI |
| `vpuid` | `vpuid` | ✓ | |
| `hl_classes` | `hl_classes` | ✓ | Comma-separated distinct `hl_class` values aggregated per POI |
| `hl_types` | `hl_types` | ✓ | Comma-separated distinct `hl_type` values aggregated per POI |
| `hl_count` | `hl_count` | ✓ | Count of hydrolocations under this POI |
| (none) | `upstream_id` | + | Nested-set pre-order index, denormalized from the owning flowpath by `flowpath_id`; see `vignette("upstream-index")` |
| (none) | `num_upstreams` | + | Size of that flowpath's upstream subtree, denormalized alongside `upstream_id` |

Rollups are populated by `finalize_data_model` via the dedicated
`hydrolocations` layer (see below).

---

## Hydrolocations

Dedicated `hydrolocations` layer emitted by `finalize_data_model`, sourced
from the pre-ngen `pois` layer (which carries the full hl_* attribute set;
the pre-ngen `events` layer has only a subset).

| Spec column | Emitted | Status | Notes |
|---|---|---|---|
| `poi_id` | `poi_id` | ✓ | FK to `pois.poi_id` |
| `flowpath_id` | `flowpath_id` | ✓ | Re-resolved from ngen `pois.flowpath_id` to pick up the `fp-` prefix |
| `vpuid` | `vpuid` | ✓ | |
| `hl_type` | `hl_type` | ✓ | |
| `hl_link` | `hl_link` | ✓ | May be pipe-separated when POIs were collapsed within 250 m |
| `hl_class` | `hl_class` | ✓ | |
| `hl_source` | `hl_source` | ✓ | |
| `hl_reference` | `hl_reference` | ✓ | Same pipe-separated alias carry-through for collapsed POIs |
| `hl_uri` | `hl_uri` | ✓ | |
| (none) | `hl_id` | + | Original hydrolocation primary key, carried through when the input supplies one (the OWP set does not) |
| (none) | `hl_flowpath_id_source` | + | The reference reach the source assigned, carried through when the input supplies it (the OWP set does not) |
| (none) | `hl_reach_change_reason` | + | Why the location moved, or `NA` if it did not. See `vignette("poi-contract")` |
| (none) | `poi_da_ref_sqkm`, `poi_da_ngen_sqkm` | + | Gage drainage area from the reference and from this fabric |
| (none) | `poi_da_ratio` | + | Relative difference, reported without a threshold. `NA` where the reference carries no usable area, i.e. unchecked rather than clean |
| (none) | `poi_drop_reason` | + | Why a hydrolocation carries no POI. Present only in a product where at least one hydrolocation could not be placed |

---

## Network

Spec lists a single wide denormalized table keyed on `flowpath_id` with a
mirror of flowpath attributes plus crosswalks to divide / POI / hydrolocation
/ external hydrofabric and per-flowline breakdowns.

We emit the network as an **aspatial attribute table** with one row per
`(flowline_id, reference_id)`, which is finer than the spec's one-row-per-fp.
A flowpath that aggregates N flowlines × M reference members yields N×M rows.

| Spec column | Emitted | Status | Notes |
|---|---|---|---|
| `flowpath_id` | `flowpath_id` | ✓ | `NA` for tributary flowlines (not a flowpath member) |
| `flowpath_toid` | `flowpath_toid` | ✓ | |
| `mainstem_id` | `mainstem_id` | ✓ | Denormalized from flowpaths by `finalize_data_model`; `NA` on tributary rows |
| `hydroseq` | `hydroseq` | ✓ | Denormalized from flowpaths; `NA` on tributary rows |
| `areasqkm` | `areasqkm` | ✓ | Flowpath-scale (not divide), denormalized from flowpaths |
| `lengthkm` | `lengthkm` | ✓ | Flowpath-scale, denormalized from flowpaths |
| `total_dasqkm` | `total_dasqkm` | ✓ | Flowpath-scale, denormalized from flowpaths |
| `vpuid` | `vpuid` | ✓ | |
| `divide_id` | `divide_id` | ✓ | |
| `poi_id` | `poi_id` | ✓ | |
| `hl_reference` | `hl_reference` | ✓ | Joined from hydrolocations on poi_id; comma-separated when a POI has multiple linked hydrolocations |
| `hf_id` | `reference_id` | ◆ | Emitted as `reference_id`. NHD/Lynker reference id (base, without `.N` split suffix) |
| `hf_toid` | `hf_toid` | ✓ | Joined from reference fabric on `reference_id` |
| `hf_source` | `hf_source` | ✓ | Per-row `source` from the reference, or `"nhdplusv2"` when the reference carries none |
| `hf_hydroseq` | `hf_hydroseq` | ✓ | Joined from reference fabric on `reference_id` |
| `flowline_id` | `flowline_id` | ✓ | |
| `flowline_toid` | `flowline_toid` | ✓ | |
| `incremental_area_id` | (none) | ✖ | Not emitted on `network` |
| `flowline_areasqkm` | `flowline_areasqkm` | ✓ | Maps to the flowline-scale incremental catchment (`incremental_areasqkm` on flowlines) |
| `flowline_hydroseq` | `flowline_hydroseq` | ✓ | |
| `flowline_lengthkm` | `flowline_lengthkm` | ✓ | |
| `flowline_total_dasqkm` | `flowline_total_dasqkm` | ✓ | |
| (none) | `is_mainstem` | + | Boolean: TRUE when this flowline lies on the flowpath's traced mainstem line |
| (none) | `slope` | + | Flowpath-scale slope, denormalized from flowpaths |
| (none) | `flowline_slope` | + | Flowline-scale slope for this member reach |
| (none) | `lake_id` | + | FK to `lakes.lake_id` for lake-member rows; `NA` otherwise |
| (none) | `reference_percent` | + | Share of the source reference reach length represented by this flowline's reconciled group; useful provenance, not in spec |
| (none) | `divide_areasqkm` | + | Denormalized from divides; useful |
| (none) | `hydrolocation_id` | + | `hl_id` from hydrolocations; useful |
| (none) | `upstream_id` | + | Nested-set pre-order index, denormalized from the owning flowpath by `flowpath_id`; see `vignette("upstream-index")` |
| (none) | `num_upstreams` | + | Size of that flowpath's upstream subtree, denormalized alongside `upstream_id` |

**Grain note:** the spec is a wide denormalized mirror of flowpath
attributes (one row per flowpath). We emit the network at a finer grain,
one row per `(flowline_id, reference_id)`, and denormalize the flowpath-scale
columns on top so every spec column is present. Tributary rows (flowlines
not on any flowpath's mainstem) have `NA` in the flowpath-scale columns.
All `hf_*` crosswalks are joined from the reference fabric on `reference_id`.

---

## Flowlines

| Spec column | Emitted | Status | Notes |
|---|---|---|---|
| `flowline_id` | `flowline_id` | ✓ | |
| `flowline_toid` | `flowline_toid` | ✓ | |
| `flowpath_id` | `flowpath_id` | ✓ | `NA` for flowlines attributed to tributaries (not on any flowpath's mainstem) |
| `flowpath_toid` | `flowpath_toid` | ✓ | Downstream `nexus_id` of the parent flowpath; `NA` on tributary rows |
| `vpuid` | `vpuid` | ✓ | |
| `mainstem_id` | `mainstem_id` | ✓ | The reference's level-path id, the same in every VPU |
| `flowline_hydroseq` | `flowline_hydroseq` | ✓ | |
| `flowpath_hydroseq` | `flowpath_hydroseq` | ✓ | Denormalized from `flowpaths.hydroseq`; `NA` on tributary rows |
| `incremental_area_id` | `incremental_area_id` | ✓ | |
| `incremental_areasqkm` | `incremental_areasqkm` | ✓ | |
| `lengthkm` | `lengthkm` | ✓ | |
| `total_dasqkm` | `total_dasqkm` | ✓ | Authoritative drainage area: divide `areasqkm` accumulated downstream over the topology |
| `streamorder` | `streamorder` | ✓ | Denormalized/propagated to fine-scale flowlines |
| `slope` | `slope` | ✓ | Flowline-scale slope (floored to 1e-5 for routing validity) |
| (none) | `reconciled_id` | + | Reconciled-flowline id this segment was split from (source lineage) |
| (none) | `reconciled_toid` | + | Downstream reconciled-flowline id |
| (none) | `fac_value` | + | FAC raster sample at the outlet cell (provenance) |
| (none) | `fac_source` | + | Which raster the FAC sample came from: `primary`, or a label naming the fallback grid |
| (none) | `lake_id` | + | FK to `lakes.lake_id` when the flowline is a lake member; `NA` otherwise |
| (none) | `flowpath_constituent` | + | Logical: TRUE marks flowlines on their flowpath's traced line; FALSE are off-network tributary/connector reaches |
| (none) | `routeable` | + | Logical: TRUE marks the routing network: the coincident (flowpath_constituent) reaches plus the incremental reaches that bridge each coincident reach to the next coincident reach downstream, so the routeable subgraph is closed under flowline_toid |
| (none) | `upstream_id` | + | Nested-set pre-order index, denormalized from the owning flowpath by `flowpath_id`; see `vignette("upstream-index")` |
| (none) | `num_upstreams` | + | Size of that flowpath's upstream subtree, denormalized alongside `upstream_id` |

---

## Nexus

| Spec column | Emitted | Status | Notes |
|---|---|---|---|
| `nexus_id` | `nexus_id` | ✓ | |
| `nexus_toid` | `nexus_toid` | ✓ | |
| `vpuid` | `vpuid` | ✓ | |

---

## Lakes

Referenced waterbodies carrying waterbody / reservoir parameters, one feature
per waterbody. The geometry follows the lake set a domain is given: points in
CONUS, from the NWM lake parameters, and polygons in Alaska (HydroLAKES), Hawaii
and PRVI (NHDPlus HR). Lake membership comes from the reference network's
waterbody attribution rather than from lake geometry. Maps to the spec's referenced-waterbody concept; the
full parameter set is a NextGen/NWM extension not enumerated in the core spec,
so these columns are marked `+`.

| Spec column | Emitted | Status | Notes |
|---|---|---|---|
| referenced waterbody | `lake_id` | + | Waterbody / reservoir primary key; FK target for `lake_id` on flowpaths, flowlines, network |
| (none) | `flowline_id` | + | Member flowline id |
| (none) | `flowpath_id` | + | Member flowpath id |
| (none) | `reference_id` | + | Reference-fabric id of the associated reach |
| (none) | `LkArea` | + | NWM lake surface area |
| (none) | `LkMxE` | + | NWM lake maximum elevation |
| (none) | `WeirC` | + | Weir coefficient |
| (none) | `WeirL` | + | Weir length |
| (none) | `OrificeC` | + | Orifice coefficient |
| (none) | `OrificeA` | + | Orifice area |
| (none) | `OrificeE` | + | Orifice elevation |
| (none) | `lat` | + | Latitude |
| (none) | `lon` | + | Longitude |
| (none) | `time` | + | NWM parameter timestamp |
| (none) | `WeirE` | + | Weir elevation |
| (none) | `ascendingIndex` | + | NWM ascending index |
| (none) | `ifd` | + | Initial fractional depth |
| (none) | `Dam_Length` | + | Dam length |
| (none) | `crs` | + | Coordinate reference system string |
| (none) | `domain` | + | Modeling domain label |
| (none) | `x` | + | Projected x coordinate |
| (none) | `y` | + | Projected y coordinate |
| (none) | `poi_id` | + | Associated POI id |
| (none) | `comid` | + | Associated NHD/Lynker reference id |
| (none) | `reservoir_index_AnA` | + | Reservoir index for the AnA configuration |
| (none) | `reservoir_index_Extended_AnA` | + | Reservoir index for Extended AnA |
| (none) | `reservoir_index_GDL_AK` | + | Reservoir index for GDL AK |
| (none) | `reservoir_index_Medium_Range` | + | Reservoir index for Medium Range |
| (none) | `reservoir_index_Short_Range` | + | Reservoir index for Short Range |
| (none) | `res_id` | + | Reservoir id |
| (none) | `vpuid` | + | VPU of the feature |

---

## Flowpath-attributes

Aspatial table of routing channel and Muskingum-Cunge parameters, keyed on
`flowpath_id`, computed from the RiverML channel table the build is given (its v2 output in
this delivery). Not part of the core spec entity set (a routing
extension), so columns are marked `+`.

| Spec column | Emitted | Status | Notes |
|---|---|---|---|
| (none) | `flowpath_id` | + | FK to `flowpaths.flowpath_id` |
| (none) | `lengthm` | + | Reach length (m) |
| (none) | `BtmWdth` | + | Channel bottom width |
| (none) | `TopWdth` | + | Channel top width |
| (none) | `ChSlp` | + | Channel side slope |
| (none) | `n` | + | Manning's n (main channel) |
| (none) | `nCC` | + | Manning's n (compound / floodplain) |
| (none) | `So` | + | Channel bed slope |
| (none) | `TopWdthCC` | + | Compound-channel top width |
| (none) | `hl_reference` | + | Hydrolocation reference for the reach |
| (none) | `MusX` | + | Muskingum X parameter |
| (none) | `MusK` | + | Muskingum K parameter |
| (none) | `attr_source` | + | Where the values came from: `riverml` = the RiverML channel table the build was given (its v2 output in this delivery), `routelink` = a routelink fallback, `geomorphic_proxy` = estimated from drainage area |
| (none) | `vpuid` | + | VPU of the feature |

---

## Flowline-attributes

Aspatial table with the same routing channel and Muskingum-Cunge parameter set as `flowpath-attributes`, keyed on `flowline_id`.

| Spec column | Emitted | Status | Notes |
|---|---|---|---|
| (none) | `flowline_id` | + | FK to `flowlines.flowline_id` |
| (none) | `lengthm` | + | Reach length (m) |
| (none) | `BtmWdth` | + | Channel bottom width |
| (none) | `TopWdth` | + | Channel top width |
| (none) | `ChSlp` | + | Channel side slope |
| (none) | `n` | + | Manning's n (main channel) |
| (none) | `nCC` | + | Manning's n (compound / floodplain) |
| (none) | `So` | + | Channel bed slope |
| (none) | `TopWdthCC` | + | Compound-channel top width |
| (none) | `hl_reference` | + | Hydrolocation reference for the reach |
| (none) | `MusX` | + | Muskingum X parameter |
| (none) | `MusK` | + | Muskingum K parameter |
| (none) | `attr_source` | + | Where the values came from: `riverml` = the RiverML channel table the build was given (its v2 output in this delivery), `routelink` = a routelink fallback, `geomorphic_proxy` = estimated from drainage area |
| (none) | `vpuid` | + | VPU of the feature |

---

## Differences from the spec

- **Network grain**: the spec is flowpath-scale wide; `network` is
  flowline × reference-member long-form.
- **Divide `type`**: a binary `network` / `landscape` classification, matching
  the three-way feature model (network with divide, network without divide, and
  landscape divide with no network), rather than the spec's `network` /
  `coastal` / `internal`.

### Columns emitted beyond the spec
- `flowpaths.rec_id_list`: comma-separated member reconciled-flowline list
- `flowpaths.member_comid`: reference-id lineage
- `flowpaths.has_divide`, `divides.has_flowline`, `divides.has_divide`: the three-way split
  (network / landscape / connector). This is what lets the pipeline
  handle all reference features without dropping any
- `divides.small_reason`: diagnostic for undersized divides
- `flowlines.reconciled_id` / `reconciled_toid`, `flowlines.fac_value` / `fac_source`: provenance
- `flowlines.flowpath_constituent`: on-flowpath-line flag (TRUE on the traced flowpath line)
- `flowlines.routeable`: routing-network membership (coincident reaches + bridge reaches; closed under flowline_toid)
- `flowpaths.lake_id`, `flowlines.lake_id`, `network.lake_id`: FK to the `lakes` layer
- `network.reference_percent`, `network.hydrolocation_id`,
  `network.divide_areasqkm`, `network.is_mainstem`, `network.slope`,
  `network.flowline_slope`: denormalization / provenance convenience fields
- `lakes`, `flowpath-attributes`, `flowline-attributes`: NWM waterbody and
  channel and Muskingum routing-parameter layers (extensions beyond the core spec)

---

## Reading the product: three rules that are easy to get wrong

Three properties of the model surprise first-time users. Each is intentional.

**1. `flowpath_toid` names a nexus, not a flowpath.** Every flowpath drains to a
nexus, and `nex-N` is the *inlet* nexus of `fp-N`. To get the downstream
flowpath, resolve through the `nexus` layer:

```sql
SELECT f.flowpath_id, n.nexus_toid AS downstream_flowpath
FROM flowpaths f JOIN nexus n ON f.flowpath_toid = n.nexus_id;
```

Terminal flowpaths point at a `tnx-N` nexus, which is a real row whose
`nexus_toid` is empty.

**2. Aggregate flowline output only over `routeable` (or
`flowpath_constituent`).** The `flowlines` layer keeps every reconciled reach,
including channels the reference fabric deliberately severed from the network:
minor divergence paths are drawn but receive no flow, carry no catchment
and have zero drainage area. Because a severed channel still *starts* at the
trunk's node, its line stacks on the trunk there and appears to share the
trunk's downstream reach. It is not a duplicate and it does not double-count:
`total_dasqkm` satisfies `DA = Σ(inflow DA) + own increment` exactly at both
grains. `routeable` marks the closed routing network; sum over that and the
stubs drop out.

**3. `has_divide = FALSE` flowpaths have no catchment, by design.** The
flowpath-to-divide relationship is 1:1 or 1:0, never 2:1, and no divide is
ever referenced by two flowpaths. A divide-less flowpath is a connector: it
routes water through terrain a neighbouring catchment already accounts for, so
`areasqkm` is 0 and its line necessarily runs inside a neighbour's divide. Code
that assumes a catchment for every flowpath must skip these; they contribute no
area.

Related: `lengthkm` on a flowpath is the length of its drawn mainstem stem, which
equals the sum of its `flowpath_constituent` member reaches, not the sum of
*all* its members. Non-constituent members are tributaries inside the catchment;
they contribute area but no length.

---

## Invariants enforced on the output

Checked by `check_invariants("ngen", ...)` at the end of each pipeline run, which
writes `invariants.json` to the VPU's `work/` directory. Representative checks:

- `ngen_fp_prefix` / `ngen_cat_prefix`: identifiers carry their `fp-` and
  `cat-`/`lnd-` prefixes
- `network_is_dag`: the fp to nexus to fp graph has no cycles (connectors
  included)
- `flowpath_coverage`: ≥ 0.90 of each flowpath's line length sits inside its
  catchment
- `divide_flowpath_exists`, `divide_id_names_a_flowpath`,
  `flowpath_divide_assigned`, `network_divide_has_owner`: referential integrity
  between flowpaths and divides in both directions, with the landscape and
  connector carve-outs accounted for
- `network_covers_all_reference_ids`: every reference id the build consumed is
  still traceable in the network table

See `vignette("invariants")`.

---

## Licensing

The layers described here are data, licensed separately from the software that
writes them. See the License section of the [README](../README.md).
