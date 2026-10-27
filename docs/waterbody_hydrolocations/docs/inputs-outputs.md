# Input data and output schemas

> Source: consolidated v3 technical documentation, draft v0.3. Preserve **TO CONFIRM** items pending code and dataset verification.


## 6. Input Data Requirements


| Dataset / field | Status | Use |
| --- | --- | --- |
| Waterbody polygons | Required | Polygon/MultiPolygon waterbody geometry; unique waterbody identifier. |
| Flowpaths | Required | Hydrofabric flowpath geometry with stable flowpath ID. |
| Downstream ID | Required/preferred | flowpath_toid / toid / equivalent for downstream connectivity. |
| totdasqkm | Preferred | Total contributing drainage area used in outlet ranking. |
| hydroseq | Preferred | Hydrologic sequence used as a ranking fallback. |
| lengthkm | Preferred | Flowpath length used as an additional ranking fallback. |
| Waterbody attributes | Optional for location; required for enrichment | e.g., nidid, waterbody_type, H_m, Dam_Length, LkMxE. |
| CRS | Required | Waterbody and flowpath layers must be compatible; processing should use a projected CRS for distance tests. |



### 6.1 Field-name normalization

Production scripts should use a field-selection helper so common variants can be accepted without rewriting the core algorithm (for example flowpath_id / ID / id / LINKNO and flowpath_toid / toid / DSLINKNO). Repository documentation should list accepted aliases and identify which are mandatory for each script.

