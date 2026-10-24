function [cls, lab, flag] = classify_v3(R, source)
%CLASSIFY_V3  Ten-class bifurcation label of one sample (Ensemble_LHS_v3).
%
%   [cls, lab, flag] = classify_v3(R)            exact label (inverse method)
%   [cls, lab, flag] = classify_v3(R, 'matcont') MATCONT validation label
%
%   R is the result struct of analyze_sample.m (one sample). For results of the v3
%   analyze_sample (R.class_version == 3) the stored ten-class label is returned
%   directly; for Ensemble_LHS_v2 results it is derived from the v2 label as below
%   (identical rule to classify_branches.m). The ten classes are
%   built from the nine v2 classes by splitting the main-branch hysteresis classes
%   according to whether the sample ALSO has a fold on a detached branch inside the
%   ATOT window, and by grouping the classes as limit points / branch points /
%   no singularity:
%
%     #   label                    v2 class   rule
%     1   LP_S_ONLY                1          main-branch fold, KTOT > PTOT, no detached fold
%     2   LP_S_DETACHED            1          as 1, plus a fold on a detached branch
%     3   LP_Z_ONLY                2          main-branch fold, KTOT < PTOT, no detached fold
%     4   LP_Z_DETACHED            2          as 3, plus a fold on a detached branch
%     5   LP_DETACHED_ONLY         7          fold only on a detached branch
%     6   BP_ONLY                  3          branch point, no fold (pitchfork)
%     7   LP_BP_LPgreaterBP        4          branch point + fold, x2(LP) > x2(BP)
%     8   LP_BP_LPlessBP           5          branch point + fold, x2(LP) < x2(BP)
%     9   NO_SING_biphasic         6          no singularity, interior maximum of x2
%    10   NO_SING_monotone         8          no singularity, monotone
%    11   OTHER                    9          anything else
%     0   FAILED                   -          computation failed
%
%   "Fold on a detached branch" = a limit point inside the window whose branch type
%   is 'detached' (R.n_LP_detached > 0 for the exact label; 'detached' in
%   R.mc_LP_type for the MATCONT label).
%
%   flag.lp_coincides_bp : LP+BP sample whose LP coincides with the BP to 1e-9
%                          (relative). These are the turning point of the crossing
%                          branch at a pitchfork, not a separate fold (31 sets in
%                          Iteration 3 and 1 in Iteration 4 of the v2 runs). They keep
%                          their label and are reported, but are excluded from the ML.
%   flag.detached_in_window : a detached-branch fold lies inside the window.

if nargin < 2, source = 'exact'; end
CI = class_info_v3(); names = CI.names;
flag = struct('lp_coincides_bp', false, 'detached_in_window', false);
isV3 = isfield(R, 'class_version') && isequal(R.class_version, 3);

switch lower(source)
    case 'exact'
        if ~R.exact_ok, cls = 0; lab = 'FAILED'; return; end
        if isV3
            cls = R.exact_class; lab = names{cls};
            flag.lp_coincides_bp = R.lp_coincides_bp; flag.detached_in_window = R.detached_in_window;
            return
        end
        c2   = R.exact_class;
        hasD = isfield(R, 'n_LP_detached') && ~isempty(R.n_LP_detached) && R.n_LP_detached > 0;
        LP = R.exact_LP_ATOT; BP = R.exact_BP_ATOT;
    case 'matcont'
        if ~isfield(R, 'mc_ok') || ~R.mc_ok, cls = 0; lab = 'N/A'; return; end
        if isV3
            cls = R.mc_class; lab = names{cls};
            t = R.mc_LP_type; flag.detached_in_window = iscell(t) && any(strcmp(t, 'detached'));
            return
        end
        c2   = R.mc_class;
        t    = R.mc_LP_type; if ischar(t), t = {t}; end
        hasD = iscell(t) && any(strcmp(t, 'detached'));
        LP = R.mc_LP_ATOT; BP = R.mc_BP_ATOT;
    otherwise
        error('classify_v3: unknown source %s', source);
end
flag.detached_in_window = hasD;

map = [1 3 6 7 8 9 5 10 11];      % v2 class -> v3 class (before the detached split)
if c2 < 1 || c2 > 9, cls = 11; lab = names{cls}; return; end
cls = map(c2);
if cls == 1 && hasD, cls = 2; end
if cls == 3 && hasD, cls = 4; end
if any(cls == [7 8]) && ~isempty(LP) && ~isempty(BP)
    flag.lp_coincides_bp = min(abs(LP(:) - BP(1))) <= 1e-9 * abs(BP(1));
end
lab = names{cls};
end
