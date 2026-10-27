# Running the workflow

> Source: consolidated v3 technical documentation, draft v0.3. Preserve **TO CONFIRM** items pending code and dataset verification.


## 9. User Execution Template


### 9.1 Software

- R with sf and dplyr.

- GeoPackage-compatible GIS for QA (QGIS recommended for visual inspection).

- Projected working CRS appropriate to the domain; Alaska production uses EPSG:3338.


### 9.2 Configuration block to expose

- Input GeoPackage/path

- Waterbody layer

- Flowpath layer

- Output GeoPackage/path

- Working CRS

- Field mappings / accepted aliases

- Alaska transition tolerances: SIDE_SAMPLE_M, INSIDE_LENGTH_TOL_M, MAX_DOWNSTREAM_STEPS

- Parameter overwrite behavior for hydraulic enrichment


### 9.3 Recommended run order


| Step | Action | Result |
| --- | --- | --- |
| A | Run domain hydrolocation selection script | Creates selected inlet/outlet points; Alaska also creates detailed QA layers. |
| B | Review QA layers in GIS | Confirm inlet/outlet placement, multi-crossing behavior, and final outlet choice. |
| C | Run schema/parameter script | Build common Hydrofabric fields; attach/estimate outlet hydraulic controls where applicable. |
| D | Validate schema and counts | Check IDs, nulls, coordinates, duplicate points, and expected one-outlet-per-waterbody behavior. |
| E | Publish products | Stage versioned GeoPackages and documentation to GitHub/AWS. |

