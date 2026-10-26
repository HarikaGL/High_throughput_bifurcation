A_Ensemble_labelling_and_ML  -  ten-class bifurcation labelling of the dual futile cycle (MAPKfull)
(See ../README.txt for the full run order, including machine learning and class transitions.)
================================================================================

SET-UP (Windows or Linux)
  Copy the whole A_Ensemble_labelling_and_ML folder INTO the MatCont7p6 folder:

    MatCont7p6\
      init.m, Continuer\, Equilibrium\, ...       (MATCONT 7p6, unchanged)
      A_Ensemble_labelling_and_ML\
        Ensemble_LHS_v3.m             main script: run this
        Ensemble_class_wheel_v3.m     class wheel: all 4 iterations in one circle, wedge
                                      = share of the ensemble, discs = real diagrams
        Ensemble_class_table_v3.m     table: parameter sets | labels | one diagram per
                                      class (editable text in the .fig)
        Ensemble_iteration_panels_v3.m  one table-row figure per iteration:
                                      constraint/bounds | class bar plot | overlaid
                                      representative diagrams per class
        Ensemble_LHS_tiles_v3.m       LHS tile view of real parameter sets (scaled
                                      position in each sampling range)
        SingleStart_v3.m              single-start MATCONT (ode15s at ATOT = 0.1, one
                                      forward continuation, branch switching at BPs)
                                      on the same parameter sets, parfor; writes
                                      Ensemble_iter<i>_v3_1e6\SingleStart\ (labels,
                                      comparison with the exact labels, core-hours).
                                      Run after Ensemble_LHS_v3.m.
        README.txt
        support\                      all helper functions (added to the path
                                      automatically; nothing to set by hand)

  Open Ensemble_LHS_v3.m in MATLAB, set ITERATION (1-4) and NumWorkers, press Run.
  The script finds support\ and MATCONT (the parent folder) by itself and writes
  its results to A_Ensemble_labelling_and_ML\Ensemble_iter<ITERATION>_v3_1e6\.
  Requirements: MATLAB with the Statistics and Machine Learning Toolbox (lhsdesign),
  the Parallel Computing Toolbox (optional, for parfor), and the Signal Processing
  Toolbox (findpeaks, only for the legacy single-start run in the MATCONT check).

SAME PARAMETER SETS AS v2
  Seed 20260930, same bounds/constraints per ITERATION, numSamples = 10000, same
  order of random-number calls. The exact analysis (ss_structure) is unchanged, so
  every parameter set gets the same singularities as in v2; only the labels are
  refined. The MATCONT check uses its own random stream (cfg.check_seed).

TEN CLASSES                                   (support\class_info_v3.m)
   1 LP_S_ONLY          5 LP_DETACHED_ONLY    9 NO_SING_biphasic
   2 LP_S_DETACHED      6 BP_ONLY            10 NO_SING_monotone
   3 LP_Z_ONLY          7 LP_BP_LPgreaterBP  11 OTHER   (0 = FAILED)
   4 LP_Z_DETACHED      8 LP_BP_LPlessBP
  LP+BP sets whose LP coincides with the BP are kept, reported and flagged
  (LP_Coincides_BP = 1, Use_For_ML = 0).

SETTINGS YOU MAY CHANGE (top of Ensemble_LHS_v3.m)
  ITERATION        1..4
  NumWorkers       parallel workers (0 = serial)
  cfg.matcont_n    samples validated by MATCONT seeded on every branch (default 1000)
  cfg.check_n      samples in the visual MATCONT check (default 10; e.g. 100)
  cfg.save_png_if  per-sample diagrams, e.g. @(R) R.multistationary

OUTPUT  (Ensemble_iter<i>_v3_1e6\)
  Excel_Files\Bifurcation_Analysis_Results_v3.xlsx
      Complete_Dataset, Class_Summary_v3, Class_Summary_v2, Normal_Forms,
      Validation, Confusion_Exact_vs_MATCONT, MATCONT_check, Settings
  Excel_Files\Normal_Forms_v3_iter<i>.csv    one row per fold / branch point
  Excel_Files\ML_Ready_Dataset_v3_iter<i>.csv (+ _for_ML.csv without flagged sets)
  Core_Hours_iter<i>.txt   core-hours of this run: exact labelling, normal forms and
      MATCONT validation (sum of per-sample times measured on the workers = serial
      CPU time), wall clock, allocated core-hours (wall x workers), CPU, RAM
  ..\Core_Hours_all_iterations.txt   one line per run (all iterations, in this folder)
  MAT_Files\Complete_Bifurcation_Results_v3.mat
  Figures\  (.png and .fig)
      Class_Distribution_v3, Branches_by_Class_v3, Confusion_Exact_vs_MATCONT,
      Normal_Forms_v3, Live_by_class
      MATCONT_check\MATCONT_check_WITH_inverse_method_pXX
          MATCONT seeded on every exact branch, grey over the exact branches
      MATCONT_check\MATCONT_check_WITHOUT_inverse_method_pXX
          original single start (ode15s at ATOT = 0.1, one continuation), grey over
          the exact branches: shows what it misses
      MATCONT_check\MATCONT_check_summary.csv, MATCONT_check_results.mat

NORMAL FORMS                                   (support\normal_forms.m)
  Computed for every fold (LP) and branch point (BP) of every sample at the exact
  location from the inverse method.
    LP: a = 1/2 p'B(q,q) (same definition and normalisation as MATCONT's s.data.a:
        ||q|| = 1, p'q = 1; sign fixed by q(2) > 0), beta = p' dF/dATOT,
        kappa = -a/beta (ATOT - ATOT_LP ~ kappa xi^2; > 0 for a fold that opens to
        larger ATOT), kappa_rel = kappa*ATOT_LP.
    BP: a, b, c of the branching equation a xi^2 + 2 b xi eta + c eta^2 = 0,
        disc = b^2 - ac, slopes of the two crossing branches, pitchfork flag (a = 0).
  MATCONT computes a only for folds and only on the branches it continues, and gives
  no BP coefficient, so it cannot cover all cases. For the MATCONT-validated samples
  its value is stored (MC_LP_nf_a, a_MATCONT) and compared with |a| (MC_LP_a_reldiff).

support\ FILES
  MAPKfull.m             model for MATCONT (unchanged)
  mc_init.m              MATCONT initialisation on client and workers
  ss_structure.m         exact steady-state structure (inverse method)  - unchanged
  ss_states_at.m, ss_state_from_uA.m, model_constants.m, jac6.m, rhs6.m  - unchanged
  analyze_sample.m       one sample: exact labels, normal forms, seeded MATCONT
  exact_labels.m         folds/BPs from the exact structure + normal forms
  classify_branches.m    ten-class rule (also returns the v2 and 7-class labels)
  class_info_v3.m        class names, plot labels, colours
  load_ensemble_v3.m     reads the 40,000 labelled sets (feature and Fig 7 scripts)
  normal_forms.m         normal-form coefficients (LP and BP)
  old_analyze_sample.m   original single-start MATCONT procedure (for the check)
  matcont_check.m        the two MATCONT check figures
  diagram_window.m       axis limits per class (linear or log-log) for the figures
  draw_branches.m, plot_sample.m, save_sample_png.m, auto_xlim.m, y_range.m  plotting
  feature_colors.m (+ .csv)  one fixed colour per ML feature (all importance figures)
  cm_graded.m            colour scale of the confusion matrices
  ml_v3_core.m           shared ML routines (data, splits, models)
