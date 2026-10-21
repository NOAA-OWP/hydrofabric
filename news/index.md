# Changelog

## hydrofabric 0.1.0

First public release. Builds the NextGen hydrofabric for the National
Water Model from the Lynker Reference Fabric and D8 flow-direction and
flow-accumulation grids.

### Pipeline

- [`refactor_fabric()`](https://NOAA-OWP.github.io/hydrofabric/reference/refactor_fabric.md)
  splits and collapses reference reaches and reconciles them and their
  divides against the grids.
- [`aggregate_fabric()`](https://NOAA-OWP.github.io/hydrofabric/reference/aggregate_fabric.md)
  aggregates the reconciled network to distribution scale, conserving
  divide area.
- [`nextgen_fabric()`](https://NOAA-OWP.github.io/hydrofabric/reference/nextgen_fabric.md)
  builds the NextGen data model: flowpaths, divides, nexuses, flowlines,
  network, hydrolocations, POIs, lakes, and RiverML channel routing
  attributes.
- [`merge_vpus()`](https://NOAA-OWP.github.io/hydrofabric/reference/merge_vpus.md)
  merges per-VPU products into the SuperCONUS, Alaska, Hawaii and PRVI
  domain products with globally unique identifiers, and records an
  upstream index for range queries
  ([`add_upstream_index()`](https://NOAA-OWP.github.io/hydrofabric/reference/add_upstream_index.md)).
- [`run_pipeline()`](https://NOAA-OWP.github.io/hydrofabric/reference/run_pipeline.md)
  runs the phases for one VPU;
  [`run_vpus()`](https://NOAA-OWP.github.io/hydrofabric/reference/run_vpus.md)
  runs a set.

### Validation

- Every phase is checked by
  [`check_invariants()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_invariants.md):
  topology (acyclic, resolvable `toid`s), divide coverage and area
  conservation, one divide per flowpath at most, and the POI placement
  contract. A VPU that fails does not ship.
- Inputs are checked before a build:
  [`check_reference_input()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_reference_input.md),
  [`check_reference_severed_terminals()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_reference_severed_terminals.md)
  and
  [`check_grid_coverage()`](https://NOAA-OWP.github.io/hydrofabric/reference/check_grid_coverage.md).

### Provenance

- [`hf_input_manifest()`](https://NOAA-OWP.github.io/hydrofabric/reference/hf_input_manifest.md)
  and
  [`hf_output_manifest()`](https://NOAA-OWP.github.io/hydrofabric/reference/hf_output_manifest.md)
  record every input and product with its checksum;
  [`hf_verify_manifest()`](https://NOAA-OWP.github.io/hydrofabric/reference/hf_verify_manifest.md)
  and
  [`hf_verify_output_manifest()`](https://NOAA-OWP.github.io/hydrofabric/reference/hf_verify_output_manifest.md)
  check them.
- Every product carries its license, attribution and build metadata in
  its `gpkg_metadata` table.
- [`prepare_hydrolocations()`](https://NOAA-OWP.github.io/hydrofabric/reference/prepare_hydrolocations.md)
  snaps and assigns a raw hydrolocation set to a reference before a
  build.

### License

Code is Apache-2.0 and data CC BY-SA 4.0, both © Lynker Corporation. See
`LICENSE.md`, `NOTICE` and `THIRD_PARTY_NOTICES.md`.
