# NOAA OWP Next Generation Modeling Framework Hydrofabric


## What is it?

Next Generation Modeling Framework (NextGen) hydrofabric artifacts are distributed by _NHDPlusV2_ **V**ector **P**rocessing **U**nits. They are generated from a set of national reference datasets built in collaboration between the USGS and NOAA for federal water modeling efforts. These artifacts are designed to be easily updated, manipulated, and quality controlled to meet the needs of a wide range of modeling tasks while leveraging the best possible input data.

## How do I get it?

NextGen artifacts are publicly available through a partnership with Lynker and the NOAA OWP. For each VPU two artifacts are available:
 
  - a geopackage that contains all tables, spatial data, and lookups relvant to a hydrofabric logical model
  - a zip folder containing the files needed to run the [Next Generation (NextGen) Water Modelling Framework](https://github.com/NOAA-OWP/ngen).

These can be programatically accessed using the respective URL patterns:

### s3

```
s3://nextgen-hydrofabric/{version}/nextgen_{VPU}.gpkg
```

```
s3://nextgen-hydrofabric/{version}/nextgen_{VPU}.zip
```

### https

```
https://nextgen-hydrofabric.s3.amazonaws.com/{version}/nextgen_{VPU}.gpkg
```

```
https://nextgen-hydrofabric.s3.amazonaws.com/{version}/nextgen_{VPU}.zip
```

### Versions

`v1.0` was released September 01, 2022 as an intial beta product.

### Disclaimer:

These data are preliminary or provisional and are subject to revision. They are being provided to meet the need for timely best science. The data have not received final approval by the National Oceanic and Atmospheric Administration (NOAA) or the U.S. Geological Survey (USGS) and are provided on the condition that neither NOAA, the USGS, nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the data.



