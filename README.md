
<!-- README.md is generated from README.Rmd. Please edit that file -->

<br>

## NOAA Next Generation Modeling Framework Hydrofabric

<!-- badges: start -->

[![R-CMD-check](https://github.com/NOAA-OWP/hydrofabric/workflows/R-CMD-check/badge.svg)](https://github.com/NOAA-OWP/hydrofabric/actions)
<!-- badges: end -->

<br>

``` r
Johnson, J. M. (2022). National Hydrologic Geospatial Fabric (hydrofabric)
for the Next Generation (NextGen) Hydrologic Modeling Framework,
HydroShare, http://www.hydroshare.org/resource/129787b468aa4d55ace7b124ed27dbde
```

### Overview

This repository serves two purpose. (1) It provides a dedicated landing
page to access the Next Generation Modeling Framework (NextGen)
hydrofabric artifacts and (2) a meta package of R packages that are
designed for hydroscience.

NextGen artifacts are distributed by *NHDPlusV2* **V**ector
**P**rocessing **U**nits and are generated from a set of national
reference datasets built in collaboration between NOAA, the USGS, and
Lynker for federal water modeling efforts. These artifacts are designed
to be easily updated, manipulated, and quality controlled to meet the
needs of a wide range of modeling tasks while leveraging the best
possible input data.

## How do I get it?

NextGen artifacts are publicly available through a partnership with
Lynker and the NOAA OWP. For each VPU a geopackage that contains all
tables, spatial data, and lookups relevant to a hydrofabric data model

### [Lynker’s NextGen Data Artifacts](https://lynker-spatial.s3.amazonaws.com/index.html)

<img src="man/figures/lynker-spatial.png" width="1529" style="display: block; margin: auto;" />

## Package Installation

``` r
# install.packages("remotes")
remotes::install_github("NOAA-OWP/hydrofabric")
```

## Usage

``` r
library(hydrofabric)
```

    ## The legacy packages maptools, rgdal, and rgeos, underpinning the sp package,
    ## which was just loaded, were retired in October 2023.
    ## Please refer to R-spatial evolution reports for details, especially
    ## https://r-spatial.org/r/2023/05/15/evolution4.html.
    ## It may be desirable to make the sf package available;
    ## package maintainers should consider adding sf to Suggests:.

    ## ── Attaching packages ────────────────────────────────────────────────────────── hydrofabric0.0.6 ──

    ## ✔ dplyr         1.1.3       ✔ nhdplusTools  1.0.1  
    ## ✔ terra         1.7.46      ✔ hydrofab      0.5.0  
    ## ✔ ngen.hydrofab 0.0.3       ✔ zonal         0.0.2  
    ## ✔ climateR      0.3.1.4     ✔ glue          1.6.2  
    ## ✔ sf            1.0.14

    ## ── Conflicts ──────────────────────────────────────────────────────────── hydrofabric_conflicts() ──
    ## ✖ terra::intersect() masks dplyr::intersect()
    ## ✖ glue::trim()       masks terra::trim()
    ## ✖ terra::union()     masks dplyr::union()

`library(hydrofabric)` will load the core packages:

- [nhdplusTools](https://github.com/usgs-r/nhdplusTools/) for network
  manipulation
- [hydrofab](https://github.com/mikejohnson51/hydrofab) a tool set for
  “fabricating” multiscale hydrofabrics
- [ngen.hydrofab](https://github.com/mikejohnson51/ngen.hydrofab)
  NextGen extensions for hydrofab
- [climateR](https://github.com/mikejohnson51/climateR) for accessing
  remote data resources for parameter and attributes estimation
- [zonal](https://github.com/mikejohnson51/zonal) for catchment
  parameter estimation

Additionally it will load key spatial data science libraries: `arrow`,
`terra`, `sf`, `dplyr` and `glue`

# Background

The NextGen artifacts are a *model application* dataset built to meet
the aims of [NextGen](https://github.com/NOAA-OWP/ngen). By design,
these artifacts are derived from a set of general authoritative data
products outlined in figure 1 that have been built in close
collaboration with the USGS.

<div class="figure" style="text-align: center">

<img src="man/figures/roadmap.png" alt="Figure 1" width="1433" />
<p class="caption">
Figure 1
</p>

</div>

These include a set of base data that improves the network topology and
geometry validity while defining a set of community hydrolocations
(POIs). These 4 data products are used to build an intermediate
refactored network from which one hydrofabric network has been
aggregated to a set of community hydrolocations (minimal network), and
one has been aggregated to a more consistent size (3-10 sqkm) with
enforced POI locations (target distribution). NextGen specifically is
derived from the target size aggregated product while the upcoming
developments on the [National Hydrologic Model
(NHM)](https://www.usgs.gov/mission-areas/water-resources/science/national-hydrologic-model-infrastructure)
will be built from the community minimal network.

While these two aggregations serve a wide range of federal modeling
needs, our focus on open source software development and workflows allow
interested parties to build there own networks starting with either the
4 reference datasets, or the refactored network!

# Resources

- The hydrofabric builds on the OGC [HY_Features conceptual
  model](https://docs.opengeospatial.org/is/14-111r6/14-111r6.html), the
  [Hydrofabric Logical model](https://docs.ogc.org/per/22-040.html), and
  the proposed [Hydrofabric Data
  Model](https://noaa-owp.github.io/hydrofabric/articles/hf_dm.html).

- The reference, refactor, minimal, and target hydrofabrics can all be
  accessed
  [here](https://www.sciencebase.gov/catalog/item/60be0e53d34e86b93891012b).
  A high level introduction to these resources can be found on the [USGS
  Water Data blog](https://waterdata.usgs.gov/blog/hydrofabric/).

# Hydrofabric Subsetter

We have created a NextGen hydrofabric subsetter. GO binaries of these
can be installed at the [release
page](https://github.com/LynkerIntel/hfsubset/releases).

<img src="man/figures/logos.png" width="1796" style="display: block; margin: auto;" />

## Questions:

<a href = "mailto:jjohnson@lynker.com?subject=Nexgen Hydrofabric Questions">
Mike Johnson</a> (Hydrofabric Lead),
<a href = "mailto:trey.flowers@noaa.gov?subject=Nexgen Hydrofabric Questions">
Trey Flowers </a> (Director, OWP Analysis and Prediction Division)

**Disclaimer**: These data are preliminary or provisional and are
subject to revision. They are being provided to meet the need for timely
best science. The data have not received final approval by the National
Oceanic and Atmospheric Administration (NOAA) or the U.S. Geological
Survey (USGS) and are provided on the condition that neither NOAA, the
USGS, nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the data.
