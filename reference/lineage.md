# Pipeline ID lineage table

A normalized lineage table: one row per leaf NHD flowline (including
post-refactor split children), with columns tracking the ID progression
through each pipeline stage. Replaces comma-separated `member_comid`
string parsing with data.table joins.

Schema:

- `nhd_comid` (character): original NHD COMID, no suffix

- `split_part` (integer) : 0 if unsplit; 1, 2, ... for post-split
  children

- `refactored_id` (character): e.g. `"85000100004325.1"` or
  `"85000100004325"`

- `reconciled_id` (integer) : reconciled group ID after
  refactor_flowlines

- `aggregated_id` (character): set by
  [`lineage_attach_aggregated()`](https://lynker-spatial.github.io/hydrofabric/reference/lineage_attach_aggregated.md)

- `final_divide_id`(character): set by
  [`lineage_attach_ngen()`](https://lynker-spatial.github.io/hydrofabric/reference/lineage_attach_ngen.md)

- `final_flowpath_id`(character): set by
  [`lineage_attach_ngen()`](https://lynker-spatial.github.io/hydrofabric/reference/lineage_attach_ngen.md)
