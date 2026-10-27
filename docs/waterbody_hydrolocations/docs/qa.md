# Quality assurance and validation

> Source: consolidated v3 technical documentation, draft v0.3. Preserve **TO CONFIRM** items pending code and dataset verification.


### 7.7 Validation and known limitations

- Verify H > 0 where used in R = L/H, nonnegative estimated lengths/areas, and valid coefficient units.

- Check that LkMxE, WeirE and OrificeE use the same vertical datum and that the estimated outlet elevations are physically consistent.

- Compare estimates with observed spillway lengths, crest elevations and outlet dimensions where available; report MAE, RMSE, median bias and missingness by region and structure type.

- Do not interpret generalized routing controls as surveyed or calibrated hydraulic design values. Replace estimates with authoritative engineering records or calibrated parameters when available.


## 10. QA/QC Checklist

- Every target waterbody has the expected number of selected outlets (normally one).

- Inlets correspond to true entry transitions and are not created from reaches that begin inside the waterbody.

- Temporary exits that re-enter downstream are not selected as final outlets.

- Same-reach inlet/outlet cases are retained when physically/topologically valid.

- No inlet duplicates the final outlet at the same coordinate.

- flowpath_id and waterbody IDs resolve to source features.

- Coordinates hl_x / hl_y are valid longitude/latitude values.

- hl_type, hl_class, hl_reference, and poi_id follow the agreed schema.

- Outlet hydraulic fields are populated only where appropriate and source/estimated values are traceable.

- Counts and QA diagnostics are recorded with each published version.


### 10.1 Known validation benchmark


| Alaska final schema work documented 618 inlet records and 312 outlet records. These counts are useful as a regression benchmark for the specific validated Alaska input version, not as a universal expected count. |
| --- |

