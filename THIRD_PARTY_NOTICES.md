# Third-party notices

`hydrofabric` is licensed under **Apache-2.0** (© 2025-2026 Lynker
Corporation). That license applies to the copyright holder’s own code in
this repository.

This package **depends on** third-party R packages that carry their own
licenses. Those packages are **not redistributed** with this source;
they are installed separately from CRAN (or their upstream source) by
the end user, so their licenses govern the user’s use of those packages,
independent of the Apache-2.0 grant here. Several are copyleft
(GPL/LGPL/MPL); using `hydrofabric` requires installing them and
complying with their terms.

This notice is provided for transparency and to answer reuse questions;
it adds no constraint to the Apache-2.0 grant on this repository’s code.

## Runtime dependencies (hydrofabric + hfutils)

| License family | Packages |
|----|----|
| **GPL-2 / GPL-3 (≥)** | `igraph` (GPL ≥2), `terra` (GPL ≥3), `lwgeom` (GPL-2), `units` (GPL-2) |
| **GPL-2 \| MIT (dual)** | `sf`, may be used under MIT |
| **LGPL (≥2.1)** | `DBI`, `future`, `RSQLite` |
| **MPL-2.0 (dual w/ file LICENSE)** | `data.table` |
| **MIT** | `cli`, `dbplyr`, `dplyr`, `duckdb`, `furrr`, `geos`, `glue`, `httr2`, `jsonlite`, `logger`, `rlang`, `rmapshaper`, `yyjsonr` |
| **Apache-2.0** | `arrow`, `hfutils` |
| **Part of R (base)** | `methods`, `parallel`, `stats`, `utils` |

## Notes for redistributors

- The copyleft licenses above (GPL/LGPL/MPL) impose obligations on
  **distribution of those packages**. Distributing `hydrofabric`’s own
  source does not redistribute them. If you bundle these dependencies
  into a combined distribution (e.g. a Docker image or a vendored
  archive), comply with each package’s license for that distribution.
- `sf` and `data.table` offer permissive options (`sf`: MIT;
  `data.table`: MPL-2.0, which is file-level copyleft only).
- This list reflects runtime dependencies. Build and test tooling
  (`knitr`, `rmarkdown`, `testthat`, `withr`) is declared in `Suggests`
  and is not loaded by the pipeline, so it does not travel with a
  distribution of this code and is out of scope here.
- Transitive dependencies carry their own licenses; standard CRAN
  tooling resolves them at install time.

## Data

The data this pipeline produces is licensed separately, under CC BY-SA
4.0. The attribution each domain requires is in `NOTICE`, and every
product carries it in its own metadata.
