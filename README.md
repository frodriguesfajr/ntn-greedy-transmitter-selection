# NTN Greedy Transmitter Selection

Reproducibility code for the manuscript:

**On Positioning in Non-Terrestrial Networks: Greedy Algorithms for Optimizing Transmitter Allocation**

This repository contains the MATLAB implementation used to reproduce the nominal NTN transmitter-selection experiments, HAPS-availability analysis, LEO-compensation analysis, fixed-size selection-necessity experiment, robustness experiments, and paper figures.

## Main scenario

The nominal candidate pool contains 26 transmitters:

- 4 synthetic HAPS
- 8 Starlink LEO satellites
- 7 O3b/O3b mPOWER MEO satellites
- 7 GEO satellites

The receiver is fixed in Rio de Janeiro, Brazil:

- Latitude: `-22.8596582 deg`
- Longitude: `-43.2303236 deg`
- Height: `10 m`

Reference epoch:

- `01-Aug-2026 12:00:00 UTC`

Elevation masks:

- HAPS: `15 deg`
- LEO/MEO/GEO: `5 deg`

The nominal transmitter-selection target is a positional bound of:

```text
0.6 m
```

## Reproducibility design

Orbital propagation is performed from the archived TLEs using the included Vallado SGP4 implementation with WGS-72 constants.

The main reproduction workflow uses the archived orbital data, fixed experiment data, and the MATLAB scripts provided in this repository.

The main orbital chain is:

```text
TLE
  -> Vallado SGP4
  -> TEME-to-ECEF conversion
  -> orbital catalog
  -> candidate pool
  -> link budget
  -> transmitter-selection experiments
```

The fixed Earth-orientation parameters used by the paper experiment are defined in:

```text
src/paper_epoch_eop.m
```

The SGP4 wrapper used by the repository is:

```text
src/propagate_tle_vallado.m
```

The validated Vallado routines used by the wrapper are stored in:

```text
third_party/vallado_sgp4/
```

## Atmospheric-loss reference

The link-budget model includes clear-sky gaseous attenuation based on ITU-R P.618/P.676.

For the fixed scenario evaluated in the manuscript, the corresponding atmospheric gaseous-loss values are provided in:

```text
data/atmospheric_gas_loss_reference.csv
```

These values are treated as fixed reference inputs for the experiment, while all remaining link-budget quantities are recomputed by `build_link_budget.m`.

This preserves the atmospheric-loss assumptions used in the original experiment while keeping the reproduction workflow transparent and deterministic.

## Repository structure

```text
ntn-greedy-transmitter-selection/
|
|-- build_orbital_catalog.m
|-- build_candidate_pool.m
|-- build_link_budget.m
|
|-- run_FA.m
|-- run_BE.m
|-- run_HAPS_restrictions.m
|-- run_LEO_availability.m
|-- run_MC_outliers.m
|-- run_selection_necessity.m
|-- generate_paper_figures.m
|
|-- setup_paths.m
|
|-- data/
|-- results/
|-- src/
|-- tests/
|-- third_party/
`-- tle/
```

### `tle/`

Archived TLEs used in the experiments, organized by architecture and NORAD catalog identifier.

### `data/`

Fixed input data, including the HAPS definition, nominal pool definition, and atmospheric gaseous-loss reference.

### `src/`

Core propagation and support functions.

### `tests/`

Validation scripts used to compare the Vallado implementation with independent reference data and with the original paper baseline.

### `results/`

Numerical outputs and paper figures retained for reproducibility and inspection. Regenerable intermediate MATLAB files are not required to be versioned.

## Running the reproduction

Clone the repository and open MATLAB in the repository root.

First configure the MATLAB path:

```matlab
clear
clc

setup_paths
```

Then execute the main pipeline in this order:

```matlab
build_orbital_catalog
build_candidate_pool
build_link_budget

run_FA
run_BE

run_HAPS_restrictions
run_LEO_availability

run_MC_outliers
run_selection_necessity

generate_paper_figures
```

The scripts write their outputs to `results/`.

## Nominal Forward Adding result

For the nominal 26-transmitter pool and a `0.6 m` target, Forward Adding returns:

```text
Final K            : 12
Final bound        : 0.583634461214 m
PDOP               : 1.391808654754
Subset evaluations : 15098
```

Selection order:

```text
1 8 9 18 2 17 19 11 5 3 10 14
```

Final composition:

```text
3 HAPS + 5 LEO + 4 MEO + 0 GEO
```

## Nominal Backward Elimination result

Backward Elimination returns the same final 12-transmitter subset:

```text
Final K            : 12
Final bound        : 0.583634461214 m
PDOP               : 1.391808654754
Subset evaluations : 285
```

Selected pool indices:

```text
1 2 3 5 8 9 10 11 14 17 18 19
```

## HAPS-availability experiment

The HAPS-restriction experiment considers maximum HAPS availability of 4, 2, and 0.

Representative results are:

| Maximum HAPS | Final K | Positional bound [m] | Target reached |
|---:|---:|---:|:---:|
| 4 | 12 | 0.583634 | Yes |
| 2 | 12 | 0.592053 | Yes |
| 0 | 22 | 0.620169 | No |

## LEO-compensation experiment

With no HAPS, the experiment increases the number of available LEO candidates while retaining 7 MEO and 7 GEO candidates.

| Available LEO | Final K | Final bound [m] | Target reached |
|---:|---:|---:|:---:|
| 8  | 22 | 0.620169 | No |
| 12 | 20 | 0.595092 | Yes |
| 16 | 20 | 0.595092 | Yes |
| 20 | 20 | 0.595092 | Yes |

## Fixed-size selection-necessity experiment

`run_selection_necessity.m` evaluates whether the transmitter identities matter when the subset size is fixed.

The experiment uses `K=12`, matching the subset size returned by FA and BE, and compares:

- FA-selected subset
- BE-selected subset
- the 12 transmitters with the largest `C/N0`
- the 12 transmitters with the largest elevation
- 50,000 uniformly sampled random subsets of the same size

The random experiment uses:

```text
Random subsets : 50000
Subset size    : K = 12
Target         : 0.6 m
Random generator: rng(2,'twister')
```

Representative results are:

| Method | K | Positional bound [m] | Target reached |
|---|---:|---:|:---:|
| FA selected | 12 | 0.583634 | Yes |
| BE selected | 12 | 0.583634 | Yes |
| Highest-C/N0 subset | 12 | 0.726835 | No |
| Highest-elevation subset | 12 | 1.559634 | No |

For the 50,000 random subsets:

```text
Random median bound                : 0.830169 m
Random 5th / 95th percentile      : 0.667815 / 1.201139 m
P(random bound <= 0.6 m)          : 0.0060 %
P(random bound <= FA/BE bound)     : 0.0000 %
FA/BE bound reduction vs. median  : 29.6969 %
```

The experiment therefore evaluates selection at a fixed number of measurements: it does not compare the 12-transmitter subset against the full 26-transmitter pool as if both represented the same resource constraint.

The generated numerical outputs are stored under:

```text
results/selection_necessity/
```

## Robustness experiment

`run_MC_outliers.m` evaluates the nominal geometry under Gaussian-impulsive pseudorange errors.

The experiment uses:

```text
Monte Carlo trials : 10000
Outliers B         : 0, 2, 4, 6
Outlier magnitude  : 20 to 100 sigma
Random generator   : rng(2,'twister')
```

The generated summary and detailed trial data are written under:

```text
results/MC_outliers/
```

## Paper figures

`generate_paper_figures.m` generates the paper figures as **separate files** so that their final placement can be controlled directly in LaTeX.

The main exported PDFs are:

```text
selection_visibility.pdf
time_sensitivity.pdf
haps_availability.pdf
leo_compensation.pdf
robustness.pdf
selection_necessity.pdf
```

PNG versions are also generated.

The figures are stored under:

```text
results/paper_figures/
```

## Validation

The `tests/` directory contains the main propagation-validation scripts.

Examples include:

```text
compare_old_vs_vallado_catalog.m
validate_all_orbital_tles.m
validate_vallado_vs_paper_baseline.m
```

The Vallado SGP4 implementation was also checked against independently generated STK orbital states before being adopted as the propagation method used by this repository.

The validation scripts are intentionally kept separate from the main reproduction pipeline.

## Notes on generated MAT files

Some scripts save `.mat` files as intermediate or detailed result files.

These files are outputs rather than hidden external dependencies. The orbital catalog, candidate pool, link budget, and experiment-specific MAT files can be regenerated by executing the corresponding scripts.

## Citation

If this repository is used in academic work, please cite the associated manuscript when its final bibliographic information becomes available.

The orbital propagation implementation is based on the Vallado SGP4 software accompanying:

> D. A. Vallado, P. Crawford, R. Hujsak, and T. S. Kelso,  
> “Revisiting Spacetrack Report #3,” AIAA/AAS Astrodynamics Specialist Conference, 2006.

## Status

This repository contains the reproducibility code corresponding to the current manuscript version. Numerical outputs should be regenerated after modifying TLEs, propagation settings, link-budget parameters, selection thresholds, or candidate-pool definitions.
