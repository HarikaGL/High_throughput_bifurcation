===============================================================================
FOLDER B - SINGLE-START COMPARISON (Supplementary File 4)
High-throughput bifurcation analysis across kinetic parameter ensembles in
systems biology models
Harika G. L., K. Sriram - Department of Computational Biology, IIIT-Delhi
===============================================================================

This folder contains the MATLAB code and results used to compare labels
obtained by single-start continuation (one ode15s start at ATOT = 0.1 followed
by one MATCONT continuation) with the exact steady-state structure (inverse
method).

The single-start labels used in the paper were computed on the SAME 40,000
parameter sets as the exact labels, with SingleStart_v3.m (folder A; outputs in
A/Ensemble_iter<i>_v3_1e6/SingleStart/: ML_Ready_Dataset_singlestart.csv,
SingleStart_vs_Exact.csv, SingleStart_summary.txt with core-hours). The run
completed for 39,638 sets (362 ended on negative concentrations). With
SINGLE_START = 'v3' (default) make_comparison_figures.m reads these files.
The labelling and machine learning of the manuscript are in folder A.

-------------------------------------------------------------------------------
1. LAYOUT
-------------------------------------------------------------------------------
2_codes/
  a_comparison_figures/make_comparison_figures.m
        Figs 8, 9 and 15B and Supplementary Figs 3-6 of the paper (output files
        EX_fig*, BP_fig*, SS_fig*, ML_fig4_analytic_conditions), and every number quoted
        in their legends (written to figure_numbers.txt). Contains an independent
        MATLAB implementation of the exact analysis (u = K/P parametrisation,
        folds, branch points, stability, normal-form coefficients).
  b_data_generation_v2/Ensemble_LHS_v2.m
        exact labelling on the single-start ensemble (8-class scheme, W = 1e6)
  b_data_generation_v2/Compare_SingleStart_vs_Exact.m
        single-start and exact labelling of the same LHS sets
  d_ML_label_comparison/ML_RC_based2.m
        ML analyses (baseline, confound audit, feature-set ablation, per
        iteration, Class 7 vs 8, pooled confusion) with the exact 8-class and
        the single-start 7-class labels
3_results/
  a_comparison_figures/        output of make_comparison_figures.m (.pdf/.png/.fig,
                               figure_numbers.txt, cache/)
  b_data_generation_v2/        output of Ensemble_LHS_v2.m
  c_ML_label_comparison/       output of ML_RC_based2.m

-------------------------------------------------------------------------------
2. RUNNING
-------------------------------------------------------------------------------
make_comparison_figures.m expects the project layout of the authors' computer;
set ROOT at the top of the script (data folders: Ensemble_iterN_v2_1e6,
data_gen/MatCont7p6/Ensemble_v3/Ensemble_iterN_v3_1e6/SingleStart/ and the
output of the ML scripts).
Base MATLAB R2020a or newer; 20-40 min the first time, ~1 min with the cache.
Ensemble_LHS_v2.m needs MATCONT 7p6 (run from inside the
MATCONT folder); parfor needs the Parallel Computing Toolbox.
The ML scripts call LightGBM, RandomForest and NaiveBayes through MATLAB's
Python interface (Python with lightgbm and scikit-learn).
