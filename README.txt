===============================================================================
SUPPLEMENTARY CODE AND DATA  (Supplementary Files 2-6)
Hidden bistability and the geometry of kinetic parameter space in dual
phosphorylation
Guturu L Harika, K. Sriram - Department of Computational Biology, IIIT-Delhi
Code:    https://github.com/HarikaGL/High_throughput_bifurcation
Archive: https://doi.org/10.5281/zenodo.23609482
===============================================================================

This package contains all code, data and results needed to reproduce every
figure, table and number of the manuscript and its Supplementary Information.
Supplementary Table 15 lists the code behind each figure and table; the same
information is given at the end of every figure legend.

Units: concentrations (K_TOT, P_TOT, A_TOT and all species) in nM; binding
constants k1, k3, k5, k7 in nM^-1 s^-1; unbinding and catalytic constants
k11, k2, k33, k4, k55, k6, k77, k8 in s^-1.

-------------------------------------------------------------------------------
1. FOLDER LAYOUT
-------------------------------------------------------------------------------
README.txt                      this file (Supplementary File 5)
A_Ensemble_labelling_and_ML/    exact ten-class labelling, machine learning,
                                class transitions, figure scripts and all their
                                results (Supplementary Files 2 and 3)
B_Single_start_comparison/      independent re-implementation of the exact
                                analysis and comparison with single-start
                                continuation (Supplementary File 4)
C_Mathematica_derivations/      Mathematica 15 notebooks of the analytical
                                results, Eqs (1)-(14) (Supplementary File 6)
S1_CRNToolbox/                  Chemical Reaction Network Toolbox input and
                                reports (Supplementary File 1)

Inside A_Ensemble_labelling_and_ML/:
  Ensemble_LHS_v3.m                 exact labelling of one sampling iteration
  Manuscript_numbers_v3.m           numbers of Tables 1-3 and 6 and of the criteria
  Ensemble_LHS_tiles_v3.m, Ensemble_class_table_v3.m, Ensemble_class_wheel_v3.m,
  Ensemble_class_wheel_single_v3.m, Ensemble_iteration_panels_v3.m,
  Fig_normal_forms_v3.m, Ensemble_bif_regenerate_v3.m       figure scripts
  SingleStart_v3.m                  single-start MATCONT labels of the same sets
  ML_v3_10class_parallel.m          machine learning (all models, all analyses)
  ML_v3_SHAP_per_class.m            per-class TreeSHAP contributions
  ML_v3_misclassified_diagrams.m, ML_v3_misclassified_figure17.m   Fig 17
  Class_transitions/                Fig 11 and the class-transition analysis
                                    (TR_00_run_all.m runs TR_01 ... TR_05;
                                    helper functions in tr_lib/)
  support/                          exact labeller (ss_structure.m, exact_labels.m,
                                    classify_branches.m), model files (MAPKfull.m,
                                    rhs6.m, jac6.m), normal forms, MATCONT
                                    interface, plotting, feature colours
  Ensemble_iter<i>_v3_1e6/          results of iteration i (i = 1..4)
  ML_v3_results_parallel/, ML_v3_results_SHAP/, Figures_combined_v3/,
  Class_transitions/results/        results of the scripts above

-------------------------------------------------------------------------------
2. SOFTWARE
-------------------------------------------------------------------------------
  MATLAB R2026a (R2020b or newer should work), with
    Statistics and Machine Learning Toolbox  lhsdesign, fitcsvm, tsne, fitcknn,
                                             perfcurve
    Parallel Computing Toolbox               optional (parfor / parfeval);
                                             without it set NumWorkers = 0 or
                                             cfg.nWorkers = 0
    Signal Processing Toolbox                findpeaks; only for the single-start
                                             run of the MATCONT check
  MATCONT 7p6 (CL_MatCont)                   https://sourceforge.net/projects/matcont/
                                             (labelling validation and single-start
                                             labels only)
  Python 3 through MATLAB's pyenv, with numpy, lightgbm, scikit-learn, xgboost,
    catboost                                 machine-learning scripts only (the
                                             models are called from MATLAB; there
                                             are no stand-alone Python scripts)
  Wolfram Mathematica 15                     folder C only
  The class-transition scripts (A/Class_transitions) need base MATLAB only and
  also run in GNU Octave 6 or newer.
  Hardware of the reported runs: AMD EPYC 9334 (32 cores), 251 GB RAM, Linux,
  64 MATLAB workers.

-------------------------------------------------------------------------------
3. HOW TO RUN
-------------------------------------------------------------------------------
3.1  Placement
  Ensemble_LHS_v3.m and SingleStart_v3.m call MATCONT. Copy the folder
  A_Ensemble_labelling_and_ML into the MATCONT 7p6 folder:
      MatCont7p6/
        init.m, Continuer/, Equilibrium/, LimitPoint/, BranchPoint/, ...  (MATCONT)
        A_Ensemble_labelling_and_ML/
  Every other script only needs to be run from its own folder; they find
  support/ and the results by relative paths.

3.2  Order of the runs (all in A_Ensemble_labelling_and_ML)
  Step 1  Ensemble_LHS_v3.m       exact ten-class labelling
          Set ITERATION = 1, 2, 3, 4 (one run each) and NumWorkers.
          Seed 20260930; 10,000 parameter sets per iteration; ATOT window
          [0, 1e6] nM. For every set: all steady states from the u = K/P
          parametrisation, folds, branch points, stability, ten-class label,
          flag "LP coincides with BP", normal-form coefficients; MATCONT seeded
          on every branch for a random subsample of 1,000 sets (validation).
          Output: Ensemble_iter<i>_v3_1e6/
             Excel_Files/Bifurcation_Analysis_Results_v3.xlsx
             Excel_Files/ML_Ready_Dataset_v3_iter<i>.csv   (input of the ML)
             Excel_Files/Normal_Forms_v3_iter<i>.csv       (one row per LP / BP)
             MAT_Files/Complete_Bifurcation_Results_v3.mat (all branches)
             Figures/ (Supplementary Figs 1, 2), Core_Hours_iter<i>.txt
          Run time: 8-18 min per iteration on 64 workers.
          SingleStart_v3.m (after step 1) computes single-start MATCONT labels of
          the same sets -> Ensemble_iter<i>_v3_1e6/SingleStart/.

  Step 2  Figure scripts and numbers (read the output of step 1; base MATLAB)
          Ensemble_LHS_tiles_v3.m          Fig 1 panel [II]
          Ensemble_class_table_v3.m        Fig 1 panel [III]
          Ensemble_class_wheel_v3.m        Fig 1 panel [IV]
          Ensemble_iteration_panels_v3.m   Figs 2-5
          Fig_normal_forms_v3.m            Fig 10
          Ensemble_class_wheel_single_v3.m class wheel of Fig 12
          Ensemble_bif_regenerate_v3.m     per-iteration class figures
          Manuscript_numbers_v3.m          Tables 1, 2, 3, 6 and the criteria
          Output: Figures_combined_v3/, Manuscript_numbers_v3.txt.
          Fig 1 and Fig 12 were assembled from these panels (Fig 1 panel [I]
          drawn by hand).

  Step 3  ML_v3_10class_parallel.m        machine learning
          Reads Ensemble_iter<i>_v3_1e6/Excel_Files/ML_Ready_Dataset_v3_iter<i>.csv;
          removes the 32 flagged sets; feature set A = the 14 sampled parameters
          and their 91 log-ratios (other sets for the ablation). Stratified
          70/10/20 split (seed 42 + repeat), 3 repeats; 5-fold CV inside the
          training set; LightGBM, XGBoost, CatBoost, random forest, gradient
          boosting, naive Bayes, most-frequent baseline. cfg.nWorkers = 64 was
          used; cfg.nWorkers = 0 runs serially. Feature colours and the colour
          scale of the confusion matrices: support/feature_colors.m, cm_graded.m.
          Output: ML_v3_results_parallel/ -> Figs 13, 14, 15A, 16;
          Supplementary Figs 7-16; Tables 7-9; Supplementary Tables 5-12.

  Step 4  ML_v3_misclassified_diagrams.m, ML_v3_misclassified_figure17.m   Fig 17
          Refit the pooled LightGBM model (same data, split and seed), list every
          misclassified test set (misclassified/Misclassified_test.csv) and draw
          exact diagrams of representative sets.

  Step 5  ML_v3_SHAP_per_class.m          Table 10; Supplementary Figs 18, 19;
          Supplementary Table 16 (Note 16). Same data, splits and LightGBM
          settings as step 3; TreeSHAP (pred_contrib) on held-out test rows.
          About 2-5 min. Output: ML_v3_results_SHAP/.

  Step 6  Class_transitions/TR_00_run_all.m
          TR_01_ready_reckoner.m    Fig 11 (ready reckoner, open circles, panel B)
          TR_02_scan.m              single-parameter scan of 1,560 sets
          TR_03_two_class_sweeps.m  two-class sweeps through each boundary
          TR_04_figures.m           Supplementary Figs 20-22, Tables 17, 18
          TR_05_window_sets.m       Supplementary Figs 23, 24, Table 19 (sets
                                    whose detached fold lies above 1e6; reads
                                    the misclassified test sets of step 4)
          Output: Class_transitions/results/, Figures_combined_v3/.
          Run time: about 50 min serial in GNU Octave; minutes with parfor.


3.3  Folder B (Supplementary File 4)
  2_codes/a_comparison_figures/make_comparison_figures.m draws Figs 7-9 and 15B
  and Supplementary Figs 3-6, contains an independent implementation of the exact
  analysis and runs the verification tests of the Methods; every number quoted
  in their legends is written to 3_results/a_comparison_figures/figure_numbers.txt.
  Set ROOT at the top of the script to the folder that contains A and B.
  2_codes/b_data_generation_v2/, c_single_start_pipeline/ and
  d_ML_label_comparison/ hold the earlier single-start pipeline and its
  comparison, kept for reference. Details: B_Single_start_comparison/README_B.txt.

3.4  Folder C (Supplementary File 6)
  Analytical_Derivations_StepByStep_final.nb   derivations of Eqs (1)-(14) and of
                                               the normal-form coefficients
  Excel_Parameter_Examples_Condition_Checks_final.nb   the conditions evaluated
                                               on one ensemble set of every class
  Open in Mathematica 15 and evaluate all cells. Details: README.txt inside.

-------------------------------------------------------------------------------
4. DATA FILES (main entry points)
-------------------------------------------------------------------------------
  A/Ensemble_iter<i>_v3_1e6/Excel_Files/ML_Ready_Dataset_v3_iter<i>.csv
      one row per parameter set: 14 parameters, Exact_Class_v3 (1-10),
      Use_For_ML (0 for the 32 flagged sets), detached branch, singularities
  A/Ensemble_iter<i>_v3_1e6/Excel_Files/Bifurcation_Analysis_Results_v3.xlsx
      complete dataset, class summaries, normal forms, MATCONT validation
  A/Ensemble_iter<i>_v3_1e6/Excel_Files/Normal_Forms_v3_iter<i>.csv
  A/Ensemble_iter<i>_v3_1e6/SingleStart/SingleStart_vs_Exact.csv
  A/Ensemble_iter<i>_v3_1e6/MAT_Files/Complete_Bifurcation_Results_v3.mat
      every branch of every set (MATLAB)
  Ten classes (support/class_info_v3.m):
     1 LP_S_ONLY          5 LP_DETACHED_ONLY    9 NO_SING_biphasic
     2 LP_S_DETACHED      6 BP_ONLY            10 NO_SING_monotone
     3 LP_Z_ONLY          7 LP_BP_LPgreaterBP  11 OTHER (no set)
     4 LP_Z_DETACHED      8 LP_BP_LPlessBP

-------------------------------------------------------------------------------
5. LICENCE AND CITATION
-------------------------------------------------------------------------------
  Code: MIT licence. Data and results: CC BY 4.0.
  Please cite the article and the Zenodo archive (DOI above).
