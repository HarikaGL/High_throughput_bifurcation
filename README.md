# Exact ensemble bifurcation analysis reveals hidden bistability and makes signalling dynamics predictable from kinetic parameters

Code for the article by **Guturu L Harika and K. Sriram** (Department of Computational Biology, IIIT-Delhi).

The code maps the kinetic parameter space of the ordered, distributive dual phosphorylation–dephosphorylation cycle onto ten codimension-1 bifurcation classes. It uses an exact parametrisation of every steady state by the free-enzyme ratio u = K/P, which finds all branches, including detached ones. On top of that labelling it adds:

- validation against MATCONT;
- normal-form coefficients;
- machine-learning classifiers;
- an analysis of the transitions between neighbouring classes.

> **Data and results** (labelled ensemble of 40,000 parameter sets, complete steady-state results, machine-learning and transition results) are archived in Zenodo: **https://doi.org/10.5281/zenodo.23609482**.
> This repository holds the code only. To reproduce the figures without rerunning the labelling, download the Zenodo archive and copy the result folders into the matching folders here (see *Data* below).

## Repository layout

| Folder | Content | Supplementary File |
|---|---|---|
| `A_Ensemble_labelling_and_ML/` | exact ten-class labelling (`Ensemble_LHS_v3.m`), figure scripts (`Fig_inverse_method_v3.m`, `Fig_criterion_v3.m`, ...), machine learning (`ML_v3_10class_parallel.m`, `ML_v3_feature_engineering.m`, `ML_v3_inference_time.m`, `ML_v3_SHAP_per_class.m`, misclassification scripts), class transitions (`Class_transitions/`), helper functions (`support/`) | 2, 3 |
| `B_Single_start_comparison/` | independent re-implementation of the exact analysis, verification tests and comparison with single-start continuation (`make_comparison_figures.m`) | 4 |
| `C_Mathematica_derivations/` | Mathematica 15 notebooks deriving Eqs (1)–(14) and checking them on ensemble sets | 6 |
| `S1_CRNToolbox/` | Chemical Reaction Network Toolbox input (`MAPKFULL.NET`) and reports | 1 |

`README.txt` (Supplementary File 5) gives the full run order, inputs, outputs, settings and run times. Supplementary Table 17 of the article lists the script behind every figure and table.

## Requirements

- MATLAB R2020b or newer (runs were made with R2026a). Toolboxes:
  - Statistics and Machine Learning Toolbox;
  - Parallel Computing Toolbox (optional);
  - Signal Processing Toolbox (only for the single-start MATCONT check).
- [MATCONT 7p6](https://sourceforge.net/projects/matcont/), needed only for the MATCONT validation and the single-start labels.
- Python 3 via MATLAB's `pyenv`, with `numpy`, `lightgbm`, `scikit-learn`, `xgboost` and `catboost`. Needed for the machine-learning scripts only; the models are called from MATLAB, and the repository contains MATLAB code only (no stand-alone Python scripts).
- Wolfram Mathematica 15, needed for folder C only.
- The class-transition scripts need base MATLAB only and also run in GNU Octave 6 or newer.

## Quick start

```matlab
% 1. exact labelling (place A_Ensemble_labelling_and_ML inside the MATCONT 7p6 folder)
cd A_Ensemble_labelling_and_ML
ITERATION = 1;  NumWorkers = 0;   % edit at the top of the script, then run
Ensemble_LHS_v3                   % repeat for ITERATION = 2, 3, 4

% 2. machine learning (Python lightgbm etc. through pyenv)
ML_v3_10class_parallel            % cfg.nWorkers = 0 runs serially
ML_v3_feature_engineering         % feature engineering and JSD selection (base MATLAB)
ML_v3_inference_time              % time per set of the trained classifier
ML_v3_misclassified_diagrams; ML_v3_misclassified_figure17
ML_v3_SHAP_per_class

% 3. ready reckoner (Fig 7) and class transitions (base MATLAB or Octave)
cd Class_transitions
TR_00_run_all

% 4. figures of the main text from the labelled ensemble (base MATLAB or Octave)
cd ..
Fig_inverse_method_v3             % Fig 3, inverse method for one parameter set
Fig_criterion_v3                  % Fig 4, criterion (10)
```

The exact labeller for a single parameter set is:

```matlab
addpath('A_Ensemble_labelling_and_ML/support')
p = [KTOT PTOT k1 k11 k2 k3 k33 k4 k5 k55 k6 k7 k77 k8];          % nM, nM^-1 s^-1, s^-1
S = ss_structure([1 p], struct('ATOT_min', 0, 'ATOT_max', 1e6, 'N', 10000));
L = exact_labels(S);   L.class                                     % class 1..10
```

## Units

- **Concentrations** (K_TOT, P_TOT, A_TOT and all species): nM.
- **Binding constants** k1, k3, k5, k7: nM⁻¹ s⁻¹.
- **Unbinding and catalytic constants** k11, k2, k33, k4, k55, k6, k77, k8: s⁻¹.

## Data

The Zenodo archive contains the full supporting package, with the same folder names as this repository plus all results:

- `Ensemble_iter<i>_v3_1e6/`, holding for each iteration i:
  - Excel/CSV labels;
  - normal forms;
  - single-start labels;
  - MAT files with every branch.
- `ML_v3_results_parallel/` and `ML_v3_results_SHAP/`;
- `Class_transitions/results/` and `Figures_combined_v3/`;
- `B_Single_start_comparison/3_results/`.

Copy them into `A_Ensemble_labelling_and_ML/` (or `B_Single_start_comparison/`) to rerun any figure script without recomputing the labels.

## Citation

Please cite the article and the Zenodo archive (see `CITATION.cff`).

## Licence

Code: MIT (see `LICENSE`). Data in the Zenodo archive: CC BY 4.0.
