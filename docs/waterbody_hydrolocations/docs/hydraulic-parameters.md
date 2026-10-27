# Hydraulic parameter estimation

> Source: consolidated v3 technical documentation, draft v0.3. Preserve **TO CONFIRM** items pending code and dataset verification.


### 7.1 Alaska hydraulic enrichment

The schema/parameter script reads validated inlet/outlet layers, joins flowpath and waterbody attributes, transforms output coordinates to EPSG:4326, and populates outlet hydraulic controls while keeping inlet records lean.



| Group | Fields / behavior |
| --- | --- |
| Source / joined | wb_id, nidid, H_m, LkMxE, Dam_Length (where available) |
| Weir | WeirC, WeirL, WeirE |
| Orifice | OrificeC, OrificeA, OrificeE |
| Behavior | Fill missing parameter values by default; optionally overwrite existing values when explicitly enabled |



### 7.2 Empirical and engineering basis for hydraulic parameter estimation

The PR/USVI reservoir inventory supplied observed or compiled dam height, dam length, spillway information and maximum lake/pool elevations. Relationships between these structural attributes and required routing controls were examined empirically. Standard weir and orifice equations and published coefficient ranges then informed physically plausible parameter forms and initial values. The resulting rules were simplified for data-limited oCONUS applications. These are empirically informed, first-order routing estimates—not site-specific engineering design parameters.


- Empirical step: compare available PR/USVI structural and spillway attributes to estimate missing routing fields.

- Engineering step: use governing discharge equations, coefficient ranges and dam/spillway guidance to bound initial estimates.

- Regionalization step: retain observed attributes when available; otherwise estimate from H, L and LkMxE, with minimum-data fallbacks.


### 7.3 Definitions, equations, and formula hierarchy

H = dam height (H_m), in metres; L = dam length (Dam_Length), in metres; R = L/H, where H > 0; LkMxE = maximum lake/pool elevation, in metres. Elevation-based outputs require a consistent vertical datum. The following is the formulation explicitly documented in the accompanying Hydraulic Parameter Estimation Methods and References note; any differing historical script variants must be versioned separately.



| Parameter | Formula or rule | Interpretation |
| --- | --- | --- |
| WeirE | LkMxE − 0.15H | Estimated crest elevation, 15% of dam height below maximum pool elevation. |
| WeirL — detailed, controlled | R ≥ 6: 0.24L; R < 6: 0.50L | Uses dam aspect ratio and spillway/control class. |
| WeirL — detailed, other/unknown | R ≥ 6: 0.32L; R < 6: 0.27L | Alternative class-dependent rule. |
| WeirL — generalized | 0.247L if L is available; 1.40H if only H is available | Minimum-data regionalization. |
| WeirC | 1.7 if H < 15 m; 2.0 if 15 ≤ H < 50 m; 2.1 if H ≥ 50 m | Representative weir discharge coefficient in the SI-form equation. |
| OrificeA | D = 0.05H; A = π(D/2)² | Equivalent circular outlet area. |
| OrificeE | WeirE − 0.20H | Representative low-level outlet elevation. |
| OrificeC | 0.65 | Representative orifice discharge coefficient. |


The governing hydraulic forms are Q = Cw × WeirL × h^(3/2) for weir discharge and Q = Co × OrificeA × √(2gh) for orifice discharge, where h is the applicable hydraulic head and g is gravitational acceleration. Coefficient values depend on the equation and units used; they should not be transferred into another discharge formulation without checking conventions.



### 7.4 Minimum-data selection and legacy variants

Preserve authoritative observed or engineering-derived parameters where available. Where values are absent, use the best supported regional method; otherwise use the documented generalized formula and finally the height-only fallback. Record the chosen method for each output parameter. This is the documented estimation intent, not a verified implementation order for every regional script.


Historical project analyses also considered WeirL = 0.19L + 0.39H and relative crest-height relationships such as 0.846H. Earlier notes used a large-dam WeirC value of 2.2 and a relative OrificeE formulation. These are distinct experimental or legacy variants, NOT substitutions for the attached methods note. Confirm which formulas are actually present in each production script before release.



### 7.5 Engineering references and evidence limitations


| Reference | Relevance |
| --- | --- |
| U.S. Bureau of Reclamation (1987), Design of Small Dams, 3rd ed. | Dam/spillway and outlet works design, discharge coefficients. |
| Brater, E. F. & King, H. W. (1976), Handbook of Hydraulics for the Solution of Hydraulic Engineering Problems, 6th ed., McGraw-Hill. | Classical weir, orifice, gate and coefficient relationships. |
| U.S. Bureau of Reclamation, Water Measurement Manual / Water Measurement with Flumes and Weirs. | Supporting guidance on weir measurement and coefficient behavior. |


Reference provenance qualification: the retained notes do not preserve a complete record of all sources actually consulted during initial PR/USVI development. These references are relevant to and support the equations and coefficient choices, but the earliest source consultation history cannot be established from the available files.



### 7.6 Parameter provenance and proposed method registry


| Recommended metadata | Purpose |
| --- | --- |
| parameter_method_id / version | Distinct identifier for each rule, e.g., weirl_controlled_Rge6_v1. |
| parameter_source / source_id | Input dataset and dam/waterbody identifier. |
| parameter_quality | Observed, regionally estimated, generalized, or fallback. |
| parameter_units / vertical_datum | Physical units and vertical reference, particularly for WeirE and OrificeE. |
| parameter_qa_flag | Missing inputs, out-of-range values, physical plausibility and datum checks. |


These are proposed metadata extensions for v4 traceability. Maintain in a companion table if the existing v3 GeoPackage schema cannot accommodate them.


## Equation quick reference

The source guide gives the following governing forms:

\[
Q_{weir}=C_w\,L_{weir}\,h^{3/2}
\]

\[
Q_{orifice}=C_o\,A_{orifice}\,\sqrt{2gh}
\]

These equations are reproduced for readability; see the formula table above for the empirical parameterization and method-specific caveats.
