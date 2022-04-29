
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hydrofabric

<!-- badges: start -->

[![R-CMD-check](https://github.com/NOAA-OWP/hydrofabric/workflows/R-CMD-check/badge.svg)](https://github.com/NOAA-OWP/hydrofabric/actions)
<!-- badges: end -->

### Overview

There are three major types of network refactoring needed to meet a
broad set of needs:

1.  One based on a flowline length criteria (routing)
2.  One that aims for a uniform catchment size (rainfall-runoff)
3.  A POI version that forces things down to a set critical locations
    (PRMS).

Additionally many instances of hydrofabric creation requires the
creation of model specific attributes, data files, and formates
(releases).

`hydrofabric` is a set of packages that work in harmony to meet these
needs. The package is designed to make it easy to install and load core
packages across users and organizations in a single command.

## Installation

``` r
# Install the development version from GitHub
# install.packages("remotes")
remotes::install_github("NOAA-OWP/hydrofabric")
```

## Usage

``` r
library(hydrofabric)
#> ── Attaching packages ────────────────────────────────────── hydrofabric0.0.2 ──
#> ✔ nhdplusTools    0.5.2          ✔ zonal           0.0.1     
#> ✔ hyRefactor      0.4.6.9011     ✔ opendap.catalog 0.0.0.9000
#> ✔ hyRelease       0.0.0.9000     ✔ eHydRo          0.0.0.9000
#> ✔ hyAggregate     0.0.1
#> ── Conflicts ──────────────────────────────────────── hydrofabric_conflicts() ──
#> ✖ hyAggregate::flowpaths_to_linestrings() masks hyRefactor::flowpaths_to_linestrings()
#> ✖ opendap.catalog::search()               masks base::search()
#> ✖ opendap.catalog::weighting_grid()       masks zonal::weighting_grid()
```

`library(hydrofabric)` will load the core packages:

-   [nhdplusTools](https://github.com/usgs-r/nhdplusTools/) for network
    manipulation

-   [hyRefactor](https://github.com/dblodgett-usgs/hyRefactor) for
    network factoring

-   [hyAggregate](https://github.com/mikejohnson51/hyAggregate) for
    network aggregation

-   [hyRelease](https://github.com/mikejohnson51/hyRelease) for data
    releases running elected subroutines

-   [opendap.catalog](https://github.com/mikejohnson51/opendap.catalog)
    for accessing remote data resources for parameter and attributes
    estimation

-   [zonal](https://github.com/mikejohnson51/zonal) for catchment
    parameter estimation

-   [eHydRo](https://github.com/mikejohnson51/eHydRo) for accessing Army
    Corp bathymetry data and embedding it in the a DEM.

## Code of Conduct

Please note that the project is released with a [Contributor Code of
Conduct](). By contributing to this project, you agree to abide by its
terms.
