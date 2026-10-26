% ======================= TR_00_run_all.m ===============================================
% Class transitions of the dual PdP cycle: runs everything in order.
%   TR_01_ready_reckoner     Fig 11 (main text): the ten classes, their regions and diagrams
%   TR_02_scan               single-parameter scan + exact-condition moves + symmetry breaking
%   TR_03_two_class_sweeps   two-class switches (sweeps of one quantity through a boundary)
%   TR_04_figures            Supplementary figures (transition map, two-class switches,
%                            symmetry breaking) and the LaTeX rows of the transition table
%   TR_05_window_sets        Supplementary Figs 27, 28 and Table 21: sets whose detached fold
%                            lies above ATOT = 1e6 (open circles of Fig 11) and the classifier
%                            errors they cause (reads Misclassified_test.csv)
% Place this folder (Class_transitions) inside A_Ensemble_labelling_and_ML: it uses the
% labeller in ../support and the ensemble files ../Ensemble_iter<i>_v3_1e6/Excel_Files/.
% Base MATLAB only; TR_02 and TR_03 use parfor (open a pool with parpool to run them in
% parallel). Run times on one core of GNU Octave: TR_01 1 min, TR_02 25 min,
% TR_03 25 min, TR_04 2 min, TR_05 2 min (MATLAB is faster; ~ 92,000 + 56,000 exact labels).
% Results: Class_transitions/results/ (results/TR_window_sets.mat caches the open circles;
% delete it to recompute); figures: ../Figures_combined_v3/Fig11_ready_reckoner
% and ../Figures_combined_v3/Class_transitions.
% ======================================================================================
clc;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
for s = {'TR_01_ready_reckoner.m', 'TR_02_scan.m', 'TR_03_two_class_sweeps.m', 'TR_04_figures.m', ...
         'TR_05_window_sets.m'}
    run(fullfile(HERE, s{1}));
    clearvars -except HERE s      % each script starts from its own default settings
end
