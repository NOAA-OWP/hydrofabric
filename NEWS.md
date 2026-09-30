# hydrofabric 0.1.0

First public release. Builds the NextGen hydrofabric for the National Water
Model from the Lynker Reference Fabric and D8 flow-direction and
flow-accumulation grids.

## Pipeline

* `refactor_fabric()` splits and collapses reference reaches and reconciles them
  and their divides against the grids.
* `aggregate_fabric()` aggregates the reconciled network to distribution scale,
  conserving divide area.
* `nextgen_fabric()` builds the NextGen data model: flowpaths, divides, nexuses,
  flowlines, network, hydrolocations, POIs, lakes, and RiverML channel routing
  attributes.
* `merge_vpus()` merges per-VPU products into the SuperCONUS, Alaska, Hawaii and
  PRVI domain products with globally unique identifiers, and records an
  upstream index for range queries (`add_upstream_index()`).
* `run_pipeline()` runs the phases for one VPU; `run_vpus()` runs a set.

## Validation

* Every phase is checked by `check_invariants()`: topology (acyclic, resolvable
  `toid`s), divide coverage and area conservation, one divide per flowpath at
  most, and the POI placement contract. A VPU that fails does not ship.
* Inputs are checked before a build: `check_reference_input()`,
  `check_reference_severed_terminals()` and `check_grid_coverage()`.

## Provenance

* `hf_input_manifest()` and `hf_output_manifest()` record every input and
  product with its checksum; `hf_verify_manifest()` and
  `hf_verify_output_manifest()` check them.
* Every product carries its license, attribution and build metadata in its
  `gpkg_metadata` table.
* `prepare_hydrolocations()` snaps and assigns a raw hydrolocation set to a
  reference before a build.

## License

Code is Apache-2.0 and data CC BY-SA 4.0, both © Lynker Corporation. See
`LICENSE.md`, `NOTICE` and `THIRD_PARTY_NOTICES.md`.
