# qra-domino-fuel-storage

MATLAB implementation of the probabilistic risk assessment pipeline described in

> Y. Cuba Arana, D. Avila, J. L. Orozco, R. Quiza, G. N. Marichal,
> **Towards a Digital Twin for the Quantitative Risk Assessment of Domino Effects in Fuel Storage Facilities**,
> submitted to *Fire* (MDPI), 2026.

The pipeline simulates primary tank fires and their fire-driven escalation (domino effect) in a fuel storage facility by Monte Carlo sampling, condenses the resulting population of thermal-radiation fields into a position-dependent Gamma distribution, encodes the Gamma parameter maps in a Generalized Regression Neural Network (GRNN) surrogate, and derives from it receptor-level risk indicators: annual exceedance frequencies of harm thresholds, expected number of fatalities (sheltered and outdoor population) and Expected Annual Loss of every asset of the complex.

The repository contains the code and the input data of the five-tank case study reported in the paper, so that every figure and table can be regenerated.

## Requirements

- MATLAB (developed and tested with R2020a
- Statistics and Machine Learning Toolbox (`mle`, `gamcdf`, `gampdf`, `normcdf`)
- Deep Learning Toolbox (`newgrnn`)

No other dependencies. All scripts are plain MATLAB files; no compilation is needed.

## Repository layout

```
.
├── m01_generate_scenarios.m        Monte Carlo generation of accident scenarios and escalation
├── m02_compute_radiation_areas.m   Radiation field of each fire configuration on the computation grid
├── m03_consolidate_data.m          (auxiliary) Flattened, frequency-weighted dataset (x, y, I)
├── m04_organize_data_dist.m        Scenario-wise radiation sample for the distribution fitting
├── m05_fit_gammas.m                Point-wise Gamma fit of the radiation intensity (maps alpha, beta)
├── m06_fir_gnn_model.m             GRNN surrogates of the Gamma parameter maps
├── m07_forescart_radiation_level.m Exceedance-frequency contour maps (Fig. 5)
├── m08_evaluate_point.m            (auxiliary) Exceedance frequencies at a single receptor
├── m09_people_affected.m           Expected number of fatalities per building (Table 7, Fig. 6)
├── m10_economic_risk.m             Expected Annual Loss of every asset (Table 8, Fig. 7)
│
├── poolfire.m            Pool fire model: burning rate and radiated power
├── radiation.m           Point-source radiation received by a tank (escalation analysis)
├── rad2point.m           Point-source radiation received at an arbitrary point of the terrain
├── escalation.m          Probit escalation model (time to failure of atmospheric tanks)
├── people_affected.m     Vulnerability models for outdoor and sheltered population
├── regression_metrics.m  R2, RMSE and maximum error of the surrogate models
├── grnn2pseudocode.m     Writes the trained GRNN as platform-independent pseudocode
│
├── data/                 Intermediate and final results (.mat), created by the scripts
│   └── scenarios/        One .mat file per fire configuration (from m02)
└── figures/              Figures (.fig, .png), created by the scripts
    ├── scenarios/
    └── radlevels/
```

The folders `data/`, `data/scenarios/`, `figures/`, `figures/scenarios/` and `figures/radlevels/` must exist before running the scripts (create them, or keep the empty `.gitkeep` files of the repository).

## Running the pipeline

Run the numbered scripts in order from the repository root:

```
m01_generate_scenarios      % ~10^8 Monte Carlo trials; the longest step
m02_compute_radiation_areas % 32 radiation fields on a 101 x 101 grid
m04_organize_data_dist
m05_fit_gammas
m06_fir_gnn_model
m07_forescart_radiation_level   % run once per threshold (edit rad_level: 2, 5, 12.5, 37.5)
m09_people_affected
m10_economic_risk
```

`m03_consolidate_data` and `m08_evaluate_point` are auxiliary and are not required by the main chain.

Every script reads its inputs from `data/` and writes its outputs to `data/` and `figures/`, so each step can be re-run independently once the previous ones have been executed.

## Case study and parameters

The layout of the facility (five vertical fixed-roof atmospheric tanks storing light naphtha, jet fuel and diesel), the industrial machinery, the buildings of the surrounding complex with their occupancy and replacement values, and the thermophysical properties of the fuels are defined at the top of the scripts that use them and correspond to Tables 1–4 of the paper.

The main model parameters are:

| Parameter | Value | Where |
|---|---|---|
| Primary accident probability per tank-year, `ACCIDENT_LH` | 2 × 10⁻⁵ | `m01` |
| Monte Carlo trials, `SCENARIOS_COUNT` | 10⁸ | `m01` |
| Burning-rate constant, `k` | 8.76 × 10⁻⁴ kg m⁻² s⁻¹ | `poolfire.m` |
| Combustion efficiency, `lambda_c` | 0.92 | `poolfire.m` |
| Radiative fraction, `chi_r` | 0.30 (constant) | `poolfire.m` |
| Pool enlargement beyond the shell, `delta` | 5 m | `poolfire.m` |
| Atmospheric transmissivity | piecewise `a·R^b` | `radiation.m`, `rad2point.m` |
| GRNN spread | 2 m | `m06` |
| Harm thresholds | 2, 5, 12.5, 37.5 kW m⁻² | `m02`, `m07`, `m08` |
| Outdoor fraction / exposure time | 0.05 / 20 s | `m09` |
| Vulnerability bands per asset class | see header of `m10` | `m10` |

## Reproducibility notes

- `m01` does not fix a random seed. Runs are statistically equivalent but not bit-identical; add `rng(<seed>)` at the top of `m01` to reproduce a run exactly.
- The number of retained scenarios `n` (numerator of the annual accident probability `n/N`) is written by `m01` as `size(id, 1)` in `data/scenarios.mat`. Scripts `m07`–`m10` currently declare `n` explicitly at their top; after running `m01`, set it to the value of your run (n = 10 048 for the run reported in the paper).
- Escalation trials at successive steps are independent (see `escalation.m` and Section 2.4 of the paper).
- The exact configuration mixture (Eq. for the exact distribution in the paper) can be evaluated at any receptor from the files in `data/scenarios/` (fields `Z` and counts `C`) and `rad2point.m`; it is the reference against which the Gamma condensation should be checked.

## License

The code is released under the MIT License (see `LICENSE`). The input data of the case study are synthetic and are released under the same terms.

## Funding

This work was funded by the Matanzas Territorial Delegation of the Cuban Ministry of Science, Technology and Environment, grant PT211MT002-022.

## Contact

Ramón Quiza [email: <ramon.quiza@xymbot.com>] — Centre for Advanced and Sustainable Manufacturing Studies, University of Matanzas / Xymbot Digital Solutions.
