
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
[here](https://noaa-owp.github.io/hydrofabric/schema.html).

## Installation

``` r
# Install the development version from GitHub
# install.packages("remotes")
remotes::install_github("NOAA-OWP/hydrofabric")
```

## Usage

``` r
library(hydrofabric)
#> ── Attaching packages ────────────────────────────────────── hydrofabric0.0.6 ──
#> ✔ dplyr         1.1.1      ✔ nhdplusTools  0.6.2 
#> ✔ terra         1.7.21     ✔ hydrofab      0.5.0 
#> ✔ ngen.hydrofab 0.0.3      ✔ zonal         0.0.2 
#> ✔ climateR      0.3.0      ✔ glue          1.6.2 
#> ✔ sf            1.0.12
#> ── Conflicts ──────────────────────────────────────── hydrofabric_conflicts() ──
#> ✖ terra::intersect() masks dplyr::intersect()
#> ✖ glue::trim()       masks terra::trim()
#> ✖ terra::union()     masks dplyr::union()
```

`library(hydrofabric)` will load the core packages:

- [nhdplusTools](https://github.com/usgs-r/nhdplusTools/) for network
  manipulation
- [hydrofab](https://github.com/mikejohnson51/hydrofab) a toolset for
  “fabricating” multiscale hydrofabrics
- [ngen.hydrofab](https://github.com/mikejohnson51/ngen.hydrofab)
  Nextgen extensions for hydrofab
- [climateR](https://github.com/mikejohnson51/climateR) for accessing
  remote data resources for parameter and attributes estimation
- [zonal](https://github.com/mikejohnson51/zonal) for catchment
  parameter estimation

Additionally it will load key r-spatial and data science libraries
including `terra`, `sf`, `dplyr` and `glue`

## Code of Conduct

Please note that the project is released with a [Contributor Code of
Conduct](). By contributing to this project, you agree to abide by its
terms.
