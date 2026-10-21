# Package index

## Pipeline

End-to-end pipeline entry points.

- [`run_pipeline()`](https://NOAA-OWP.github.io/hydrofabric/reference/run_pipeline.md)
  : Run the end-to-end hydrofabric pipeline for a single VPU
- [`run_vpus()`](https://NOAA-OWP.github.io/hydrofabric/reference/run_vpus.md)
  : Build a set of VPUs, then merge them
- [`presplit_reference()`](https://NOAA-OWP.github.io/hydrofabric/reference/presplit_reference.md)
  : Pre-split a multi-VPU reference fabric into per-VPU slices
- [`prepare_hydrolocations()`](https://NOAA-OWP.github.io/hydrofabric/reference/prepare_hydrolocations.md)
  : Prepare an externally-supplied hydrolocation set for embedding
- [`check_reference_input()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_reference_input.md)
  : Validate a reference fabric before building from it
- [`check_grid_coverage()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_grid_coverage.md)
  : Check that the supplied grids cover the domain
- [`check_hydrolocation_input()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_hydrolocation_input.md)
  : Validate a hydrolocation set before building from it
- [`check_reference_coincident_geometry()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_reference_coincident_geometry.md)
  : Find reference reaches drawn on top of one another
- [`check_reference_divide_pairing()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_reference_divide_pairing.md)
  : Find reference reaches that end outside their own divide
- [`check_reference_downstream_resolves()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_reference_downstream_resolves.md)
  : Find reference reaches whose downstream neighbour is missing
- [`check_reference_severed_terminals()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_reference_severed_terminals.md)
  : Find reference reaches that terminate on a channel that carries on
- [`hf_compare_products()`](https://NOAA-OWP.github.io/hydrofabric/reference/hf_compare_products.md)
  : Compare a rebuilt product against the one it replaces
- [`hl_reference_index()`](https://NOAA-OWP.github.io/hydrofabric/reference/hl_reference_index.md)
  : Index a reference for hydrolocation checking
- [`hf_input_manifest()`](https://NOAA-OWP.github.io/hydrofabric/reference/hf_input_manifest.md)
  : Record the inputs a build read
- [`hf_verify_manifest()`](https://NOAA-OWP.github.io/hydrofabric/reference/hf_verify_manifest.md)
  : Verify inputs against a recorded manifest
- [`hf_output_manifest()`](https://NOAA-OWP.github.io/hydrofabric/reference/hf_output_manifest.md)
  : Record what a build produced
- [`hf_verify_output_manifest()`](https://NOAA-OWP.github.io/hydrofabric/reference/hf_verify_output_manifest.md)
  : Verify products against an output manifest

## Refactor

Split/collapse flowlines and raster-partition divides.

- [`refactor_fabric()`](https://NOAA-OWP.github.io/hydrofabric/reference/refactor_fabric.md)
  : Refactor a reference hydrofabric (split, collapse, reconcile;
  refactor phase)
- [`refactor_flowlines()`](https://NOAA-OWP.github.io/hydrofabric/reference/refactor_flowlines.md)
  : Refactor reference flowlines into reconciled reaches
- [`reconcile_collapsed_flowlines()`](https://NOAA-OWP.github.io/hydrofabric/reference/reconcile_collapsed_flowlines.md)
  : Reconcile collapsed flowlines
- [`split_flowlines()`](https://NOAA-OWP.github.io/hydrofabric/reference/split_flowlines.md)
  : Split flowlines (sequential, flowpath_id-based)
- [`events_to_exclude_cats()`](https://NOAA-OWP.github.io/hydrofabric/reference/events_to_exclude_cats.md)
  : Build exclude_cats from event/split_event data
- [`reconcile_divides()`](https://NOAA-OWP.github.io/hydrofabric/reference/reconcile_divides.md)
  : Reconcile divides to refactored flowlines (CRS-aware, lr-capable)
- [`partition_reconcile_divides()`](https://NOAA-OWP.github.io/hydrofabric/reference/partition_reconcile_divides.md)
  : Partition and Reconcile Divides with Raster-Based Splitting
- [`partition_reconcile_divides_from_objects()`](https://NOAA-OWP.github.io/hydrofabric/reference/partition_reconcile_divides_from_objects.md)
  : Partition and Reconcile Divides from In-Memory Objects

## Aggregate

Aggregate the network and divides to distribution scale.

- [`aggregate_fabric()`](https://NOAA-OWP.github.io/hydrofabric/reference/aggregate_fabric.md)
  : Aggregate a refactored fabric to the distribution network (aggregate
  phase)
- [`aggregate_network()`](https://NOAA-OWP.github.io/hydrofabric/reference/aggregate_network.md)
  : Aggregate flowpath network topology without divide geometry
- [`construct_agg_divides()`](https://NOAA-OWP.github.io/hydrofabric/reference/construct_agg_divides.md)
  : Construct aggregated divides from Stage-1 rec_div via rec_id_list
- [`prepare_network()`](https://NOAA-OWP.github.io/hydrofabric/reference/prepare_network.md)
  : Prepare Hydrologic Network
- [`add_network_type()`](https://NOAA-OWP.github.io/hydrofabric/reference/add_network_type.md)
  : Annotate Flowpaths and Divides with Type Info
- [`extend_flowpaths_to_mainstem()`](https://NOAA-OWP.github.io/hydrofabric/reference/extend_flowpaths_to_mainstem.md)
  : Extend aggregated flowpath geometry along its main-stem
  sub-flowlines

## NextGen

Build the NextGen product and finalize the data model.

- [`nextgen_fabric()`](https://NOAA-OWP.github.io/hydrofabric/reference/nextgen_fabric.md)
  : Build the NextGen data model from refactored + aggregated layers
  (nextgen phase)
- [`nextgen_from_gpkg()`](https://NOAA-OWP.github.io/hydrofabric/reference/nextgen_from_gpkg.md)
  : Build a NextGen hydrofabric from a single GeoPackage
- [`add_flowlines()`](https://NOAA-OWP.github.io/hydrofabric/reference/add_flowlines.md)
  : Add Sub-Flowline Segments to Aggregated Hydrofabric
- [`build_aspatial_network()`](https://NOAA-OWP.github.io/hydrofabric/reference/build_aspatial_network.md)
  : Build an aspatial long-form network crosswalk
- [`enforce_flowpath_flowdir()`](https://NOAA-OWP.github.io/hydrofabric/reference/enforce_flowpath_flowdir.md)
  : Enforce upstream-\>downstream orientation on aggregated flowpath
  geometry
- [`flag_flowpath_constituents()`](https://NOAA-OWP.github.io/hydrofabric/reference/flag_flowpath_constituents.md)
  : Flag whether each flowline is a direct constituent of its flowpath's
  line
- [`flag_routeable_flowlines()`](https://NOAA-OWP.github.io/hydrofabric/reference/flag_routeable_flowlines.md)
  : Flag the routeable flowlines (the routing network)
- [`add_upstream_index()`](https://NOAA-OWP.github.io/hydrofabric/reference/add_upstream_index.md)
  : Add a nested-set upstream index to a hydrofabric GeoPackage
- [`finalize_data_model()`](https://NOAA-OWP.github.io/hydrofabric/reference/finalize_data_model.md)
  : Finalize NextGen output to the Lynker hydrofabric data model spec

## Global merge

Merge per-VPU outputs into domain products with global ids.

- [`merge_vpus()`](https://NOAA-OWP.github.io/hydrofabric/reference/merge_vpus.md)
  : Merge per-VPU fabrics into domain products with globally unique IDs

## QA, lineage & invariants

The correctness contract and reference-id traceability.

- [`check_invariants()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_invariants.md)
  : Pipeline invariant checks

- [`enforce_integer_ids()`](https://NOAA-OWP.github.io/hydrofabric/reference/enforce_integer_ids.md)
  : Enforce integer typing and indexes on hydrofabric id columns

- [`lineage_build()`](https://NOAA-OWP.github.io/hydrofabric/reference/lineage_build.md)
  : Build the initial lineage table from refactor_flowlines output

- [`lineage_attach_aggregated()`](https://NOAA-OWP.github.io/hydrofabric/reference/lineage_attach_aggregated.md)
  : Attach aggregated flowpath IDs to a lineage table

- [`lineage_attach_ngen()`](https://NOAA-OWP.github.io/hydrofabric/reference/lineage_attach_ngen.md)
  :

  Attach NextGen IDs (`fp-N`, `cat-N`) to a lineage table

- [`lineage_member_csv()`](https://NOAA-OWP.github.io/hydrofabric/reference/lineage_member_csv.md)
  :

  Recompute `member_comid` CSV strings from a lineage table
