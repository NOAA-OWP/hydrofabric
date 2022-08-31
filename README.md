
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hydrofabric

<!-- badges: start -->

[![R-CMD-check](https://github.com/NOAA-OWP/hydrofabric/workflows/R-CMD-check/badge.svg)](https://github.com/NOAA-OWP/hydrofabric/actions)
<!-- badges: end -->

### Overview

This repository serves two purpose. (1) it provides a dedicated landing
page to access the community nextgen artifacts and (2) it provides a
single source download collection of R packages designed for data
science.

## Data

The data and documentation for the Nextgen community resources can be
found [here](https://noaa-owp.github.io/hydrofabric/) with detailed
documentation of the geopackage contents
[here](https://noaa-owp.github.io/hydrofabric/schema.html)

## Installation

``` r
# Install the development version from GitHub
# install.packages("remotes")
remotes::install_github("NOAA-OWP/hydrofabric")
```

## Usage

``` r
library(hydrofabric)
#> ── Attaching packages ────────────────────────────────────── hydrofabric0.0.3 ──
#> ✔ nhdplusTools    0.5.7          ✔ zonal           0.0.1     
#> ✔ hydrofab        0.4.7          ✔ opendap.catalog 0.0.0.9000
#> ✔ ngen.hydrofab   0.0.3
#> ── Conflicts ──────────────────────────────────────── hydrofabric_conflicts() ──
#> ✖ opendap.catalog::search() masks base::search()
```

`library(hydrofabric)` will load the core packages:

-   [nhdplusTools](https://github.com/usgs-r/nhdplusTools/) for network
    manipulation
-   [hydrofab](https://github.com/mikejohnson51/hydrofab) a toolset for
    “fabricating” multiscale hydrofabrics
-   [ngen.hydrofab](https://github.com/mikejohnson51/ngen.hydrofab)
    Nextgen extensions for hydrofab
-   [opendap.catalog](https://github.com/mikejohnson51/opendap.catalog)
    for accessing remote data resources for parameter and attributes
    estimation
-   [zonal](https://github.com/mikejohnson51/zonal) for catchment
    parameter estimation

## Code of Conduct

Please note that the project is released with a [Contributor Code of
Conduct](). By contributing to this project, you agree to abide by its
terms.
