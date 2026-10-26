% ======================= Ensemble_LHS_v3.m =========================================
% Bifurcation-class labelling of the dual futile cycle (MAPKfull.m) over an LHS
% ensemble, TEN classes, with normal-form coefficients and a visual MATCONT check.
%
% HOW TO RUN
%   Put the whole folder A_Ensemble_labelling_and_ML inside the MatCont7p6 folder:
%       MatCont7p6/
%         A_Ensemble_labelling_and_ML/
%           Ensemble_LHS_v3.m              <- this script (run it)
%           README.txt
%           support/                       <- all helper functions (added automatically)
%   Open this script, set ITERATION (1..4) below and press Run. The paths to the helper
%   functions and to MATCONT are set automatically; results are written to
%   Ensemble_v3/Ensemble_iter<ITERATION>_v3_1e6/.
%
% THE TEN CLASSES (support/class_info_v3.m, support/classify_branches.m)
%    limit points    1 LP_S_ONLY          main-branch fold, S-shaped, no detached fold
%                    2 LP_S_DETACHED      S-shaped AND a fold on a detached branch
%                    3 LP_Z_ONLY          main-branch fold, Z-shaped, no detached fold
%                    4 LP_Z_DETACHED      Z-shaped AND a fold on a detached branch
%                    5 LP_DETACHED_ONLY   fold only on a detached branch
%    branch points   6 BP_ONLY            branch point, no fold (pitchfork)
%                    7 LP_BP_LPgreaterBP  branch point + fold, x2(LP) > x2(BP)
%                    8 LP_BP_LPlessBP     branch point + fold, x2(LP) < x2(BP)
%    no singularity  9 NO_SING_biphasic   interior maximum of x2
%                   10 NO_SING_monotone   monotone
%                   11 OTHER              anything else;   0 = FAILED
%   "Fold on a detached branch" = inside the ATOT window. LP+BP sets whose LP coincides
%   with the BP (turning point of the crossing branch at a pitchfork; all LP+BP sets of
%   Iteration 3 and one of Iteration 4 in the v2 runs) keep their label, are reported
%   everywhere and flagged (LP_Coincides_BP = 1, Use_For_ML = 0): too few for
%   iteration-wise learning, so they are left out of the ML data.
%
% SAME PARAMETER SETS AS v2
%   Same seed (20260930), same bounds and constraints per ITERATION, same number of
%   samples and the same order of random-number calls (lhsdesign, then the randperm that
%   picks the MATCONT-validated samples). The MATCONT check below draws its samples
%   from a SEPARATE random stream, so it does not change anything else.
%
% NORMAL FORMS (support/normal_forms.m)
%   For EVERY fold and branch point of EVERY sample, at the exact location given by the
%   inverse method:  folds: a (= MATCONT's coefficient, same normalisation), beta,
%   kappa = -a/beta;  branch points: a, b, c of the branching equation, discriminant,
%   slopes of the crossing branches, pitchfork flag. MATCONT reports a only for folds
%   and only on the continued branches, so it cannot cover all cases; its value
%   (s.data.a) is stored for the MATCONT-validated samples and compared (|a|).
%
% MATCONT CHECK (support/matcont_check.m), on cfg.check_n random samples:
%   Figures/MATCONT_check/MATCONT_check_WITH_inverse_method_pXX     (.png, .fig)
%       MATCONT seeded on every exact branch, grey over the exact branches
%   Figures/MATCONT_check/MATCONT_check_WITHOUT_inverse_method_pXX  (.png, .fig)
%       the original single start (ode15s at ATOT = 0.1, one continuation), grey over
%       the exact branches: shows what it misses (detached branches, ...)
% ===================================================================================
clear; close all; clc

%% ================= PATHS (automatic) ================================================
HERE = fileparts(mfilename('fullpath'));          % the Ensemble_v3 folder
if isempty(HERE), HERE = pwd; end
SUPPORT_DIR = fullfile(HERE, 'support');          % helper functions + MAPKfull.m
addpath(SUPPORT_DIR);
cd(HERE);                                         % so no older copy of a helper shadows ours
if ~check_worker_path(SUPPORT_DIR)
    error('An older analyze_sample.m shadows support/analyze_sample.m: %s', which('analyze_sample'));
end
MATCONT_DIR = fileparts(HERE);                    % Ensemble_v3 sits inside MatCont7p6
if ~exist(fullfile(MATCONT_DIR, 'init.m'), 'file') || ~exist(fullfile(MATCONT_DIR, 'Continuer'), 'dir')
    error(['MATCONT not found in %s.\nPut the Ensemble_v3 folder inside the MatCont7p6 ' ...
           'folder (or set MATCONT_DIR by hand).'], MATCONT_DIR);
end
PROJECT_DIR = SUPPORT_DIR;                        % added to the path of every worker

%% ================= USER SETTINGS ====================================================
ITERATION    = 1;           % 1..4: presets identical to the v2 runs (see below)
numSamples   = 10000;       % as in the v2 runs (changing it changes the LHS design)
rngSeed      = 20260930;    % SAME seed as v2 -> same parameter sets
NumWorkers   = 64;          % 0 = run serially

% ---- iteration presets (read from the Settings sheets of the v2 output) ----------
switch ITERATION
    case 1, type = 1; tot_list = {};            bounds_set = 'S_dominant';  % kinase-dominant
    case 2, type = 1; tot_list = {};            bounds_set = 'Z_dominant';  % phosphatase-dominant
    case 3, type = 2; tot_list = {};            bounds_set = 'wide';        % full symmetry
    case 4, type = 2; tot_list = {'k5','k55'};  bounds_set = 'wide';        % first-binding asymmetry
    otherwise, error('ITERATION must be 1, 2, 3 or 4');
end
outputDir = fullfile(HERE, sprintf('Ensemble_iter%d_v3_1e6', ITERATION));

cfg = struct();
%cfg.ATOT_window    = [0 1e4];  % SAME window for every sample; labels refer to it
cfg.ATOT_window    = [0 1e6];  % SAME window for every sample; labels refer to it

cfg.ATOT_seed      = 0.1;      % where the main branch is seeded (the old base value)
%cfg.exact_N        = 6000;     % grid points per branch in ss_structure
cfg.exact_N        = 10000;     % grid points per branch in ss_structure

cfg.run_matcont    = true;     % independent MATCONT check (validation statistic)
cfg.matcont_n      = 1000;     % random subsample that gets the MATCONT check
                               % (Inf = all samples; ~1-2 s per sample per core)
cfg.keep_mc_curves = true;
cfg.neg_tol        = 1e-8;     % relative tolerance for "negative concentration"
cfg.curve_points   = 600;      % points stored per branch (for plots)

% ---- visual checking while the ensemble runs ----
cfg.live_plot          = true;                   % live figures, updated as samples finish
cfg.live_max_per_class = 40;                     % curves overlaid per class panel
cfg.live_flag          = @(R) R.multistationary; % samples shown in the "latest flagged" grid
                                                 % e.g. @(R) R.mushroom_candidate || R.exact_class == 9
cfg.save_png_if        = @(R) false;             % per-sample PNG + .fig into Figures/Samples,
                                                 % e.g. @(R) R.multistationary  or  @(R) true
cfg.fig_per_class      = 6;                      % samples per panel in Branches_by_Class_v3

% ---- MATCONT performance check (grey MATCONT curves over the exact branches) -------
cfg.check_n    = 10;            % number of random samples to check (e.g. 10 or 100; 0 = off)
cfg.check_seed = rngSeed + 1;   % separate random stream: does not change the samples

%% ================= OUTPUT FOLDERS ===================================================
figuresDir = fullfile(outputDir,'Figures');
excelDir   = fullfile(outputDir,'Excel_Files');
matDir     = fullfile(outputDir,'MAT_Files');
sampleFigDir = fullfile(figuresDir,'Samples');
for d = {outputDir, figuresDir, excelDir, matDir}
    if ~exist(d{1},'dir'), mkdir(d{1}); end
end
total_start_time = tic;

%% ================= PARAMETERS AND LHS ===============================================
param_names = {'ATOT','KTOT','PTOT','k1','k11','k2','k3','k33','k4','k5','k55','k6','k7','k77','k8'};
Param_Base_Vec = [0.1,100,100,0.02,1,0.01,0.0024,1,0.21,0.198,1,0.01,0.0024,1,0.2];
numParams = numel(Param_Base_Vec);
vary_list = param_names(2:end);

switch bounds_set
    case 'wide'
        for p = vary_list, bounds.(p{1}) = [1e-3 1e1]; end
        bounds.KTOT = [10 500]; bounds.PTOT = [10 500];
    case 'Z_dominant'
        bounds.KTOT=[10 100];   bounds.PTOT=[100 500]; bounds.k1=[1e-1 1e1];  bounds.k11=[1e-3 1e-1];
        bounds.k2=[1e-3 1e-1];  bounds.k3=[1e-1 1e1];  bounds.k33=[1e-3 1e-1]; bounds.k4=[1e-1 1e1];
        bounds.k5=[1e-3 1e-2];  bounds.k55=[1e-3 1e1]; bounds.k6=[1e-3 1e-1];  bounds.k7=[1e-3 1e1];
        bounds.k77=[1e-3 1e1];  bounds.k8=[1e-1 1e1];
    case 'S_dominant'
        bounds.KTOT=[100 500];  bounds.PTOT=[10 100];  bounds.k1=[1e-3 1e-2];  bounds.k11=[1e-3 1e1];
        bounds.k2=[1e-3 1e-1];  bounds.k3=[1e-3 1];    bounds.k33=[1e-3 1e1];  bounds.k4=[1e-1 1e1];
        bounds.k5=[1e-1 1e1];   bounds.k55=[1e-3 1e-1];bounds.k6=[1e-3 1e-1];  bounds.k7=[1e-1 1e1];
        bounds.k77=[1e-3 1e-1]; bounds.k8=[1e-3 1e1];
    otherwise
        error('Unknown bounds_set %s', bounds_set);
end

rng(rngSeed, 'twister');
lhsNorm = lhsdesign(numSamples, numParams, 'criterion','maximin','iterations',50);
paramSamples = repmat(Param_Base_Vec, numSamples, 1);   % ATOT stays at its base value
logsample = @(name, col) 10.^(log10(bounds.(name)(1)) + ...
                              lhsNorm(:,col)*(log10(bounds.(name)(2)) - log10(bounds.(name)(1))));
idx = @(name) find(strcmp(param_names, name));
if type == 1
    for p = 2:numParams
        paramSamples(:,p) = logsample(param_names{p}, p);
    end
else
    % kinase/phosphatase symmetry: PTOT=KTOT, k5=k1, k55=k11, k6=k2, k7=k3, k77=k33, k8=k4
    pairs = {'KTOT','PTOT'; 'k1','k5'; 'k11','k55'; 'k2','k6'; 'k3','k7'; 'k33','k77'; 'k4','k8'};
    for i = 1:size(pairs,1)
        v = logsample(pairs{i,1}, idx(pairs{i,1}));
        paramSamples(:, idx(pairs{i,1})) = v;
        paramSamples(:, idx(pairs{i,2})) = v;
    end
    for i = 1:numel(tot_list)                    % symmetry-breaking parameters
        paramSamples(:, idx(tot_list{i})) = logsample(tot_list{i}, idx(tot_list{i}));
    end
end
fprintf('Generated %d parameter sets (type %d, %s bounds)\n', numSamples, type, bounds_set);

%% ================= MATCONT OPTIONS ==================================================
needMC = cfg.run_matcont || cfg.check_n > 0;
if needMC
    mc_init(MATCONT_DIR, PROJECT_DIR);
    W = cfg.ATOT_window(2);
    opt = contset;
    opt = contset(opt,'InitStepSize', 1e-4);
    opt = contset(opt,'MinStepSize', 1e-10);
    opt = contset(opt,'MaxStepSize', W/1e4);    % scaled to the window, not fixed at 1e-3
    opt = contset(opt,'MaxNewtoniters', 10);
    opt = contset(opt,'MaxCorriters', 15);
    opt = contset(opt,'MaxTestiters', 20);
    opt = contset(opt,'VarTolerance', 1e-8);
    opt = contset(opt,'FunTolerance', 1e-8);
    opt = contset(opt,'TestTolerance', 1e-7);
    opt = contset(opt,'Singularities', 1);
    opt = contset(opt,'Locators', [1 1 1]);
    opt = contset(opt,'Adapt', 3);
    opt = contset(opt,'MaxNumPoints', 20000);   % enough arclength to cross the window
    opt = contset(opt,'CheckClosed', 50);
    cfg.contset = opt;

    % options of the ORIGINAL single-start procedure (used only by the MATCONT check)
    oldOpt = contset;
    oldOpt = contset(oldOpt,'InitStepSize',1);
    oldOpt = contset(oldOpt,'MaxStepSize',1e-3);
    oldOpt = contset(oldOpt,'MinStepSize',1e-5);
    oldOpt = contset(oldOpt,'MaxNewtoniters',10);
    oldOpt = contset(oldOpt,'MaxCorriters',15);
    oldOpt = contset(oldOpt,'MaxTestiters',10);
    oldOpt = contset(oldOpt,'VarTolerance',1e-7);
    oldOpt = contset(oldOpt,'FunTolerance',1e-7);
    oldOpt = contset(oldOpt,'TestTolerance',1e-6);
    oldOpt = contset(oldOpt,'TestFunctions',1);
    oldOpt = contset(oldOpt,'Locators',[1 1 1]);
    oldOpt = contset(oldOpt,'Singularities',1);
    oldOpt = contset(oldOpt,'Adapt',5);
    oldOpt = contset(oldOpt,'MaxNumPoints',1500);
    oldOpt = contset(oldOpt,'CheckClosed',50);
    oldOpt = contset(oldOpt,'Backward',0);
end
addpath(SUPPORT_DIR);   % keep the helper functions first on the path

%% ================= PARALLEL POOL ====================================================
useParallel = NumWorkers > 0 && license('test','Distrib_Computing_Toolbox');
if useParallel
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool('local', NumWorkers); end
    % A pool left open from an earlier run keeps that run's current folder and path on
    % the workers, where an older analyze_sample.m (e.g. the v2 one in data_gen) would
    % shadow support/analyze_sample.m. Move every worker to this folder and check.
    wait(parfevalOnAll(pool, @cd, 0, HERE));
    wait(parfevalOnAll(pool, @addpath, 0, PROJECT_DIR));
    if needMC
        wait(parfevalOnAll(pool, @mc_init, 0, MATCONT_DIR, PROJECT_DIR));
    end
    wait(parfevalOnAll(pool, @addpath, 0, PROJECT_DIR));   % helpers first on the path
    okW = fetchOutputs(parfevalOnAll(pool, @check_worker_path, 1, SUPPORT_DIR));
    if ~all(okW)
        error(['%d of %d workers do not use support/analyze_sample.m (an older copy shadows ' ...
               'it). Close the pool with  delete(gcp(''nocreate''))  and run again.'], ...
               sum(~okW), numel(okW));
    end
    nw = pool.NumWorkers;
else
    nw = 0;
end

%% ================= MAIN LOOP (parallel) =============================================
fprintf('\n=== Analysing %d samples on %d workers ===\n', numSamples, nw);
analysis_start_time = tic;
resC = cell(numSamples,1);
doMC = false(numSamples,1);
if cfg.run_matcont
    doMC(randperm(numSamples, min(numSamples, cfg.matcont_n))) = true;   % uses rngSeed stream
end
on_result([], cfg, numSamples);                   % reset progress counter and live figures
if useParallel
    % Workers cannot draw on the client's screen, so each finished sample is sent back
    % through a DataQueue and plotted here, on the client, as soon as it arrives.
    dq = parallel.pool.DataQueue;
    afterEach(dq, @(d) on_result(d, cfg, numSamples));
    parfor (k = 1:numSamples, nw)
        R = analyze_sample(paramSamples(k,:), cfg, doMC(k));
        if cfg.save_png_if(R), save_sample_png(R, paramSamples(k,:), k, sampleFigDir); end
        resC{k} = R;
        send(dq, struct('k', k, 'R', R, 'P', paramSamples(k,:)));
    end
else
    for k = 1:numSamples
        R = analyze_sample(paramSamples(k,:), cfg, doMC(k));
        if cfg.save_png_if(R), save_sample_png(R, paramSamples(k,:), k, sampleFigDir); end
        resC{k} = R;
        on_result(struct('k', k, 'R', R, 'P', paramSamples(k,:)), cfg, numSamples);
    end
end
hLive = findobj('Type','figure','Name','Live: all samples by class');
if cfg.live_plot && ~isempty(hLive)
    saveas(hLive(1), fullfile(figuresDir, 'Live_by_class.png'));
    savefig(hLive(1), fullfile(figuresDir, 'Live_by_class.fig'));   % editable MATLAB figure
end
res = [resC{:}]';
if ~isfield(res, 'class_version')
    error(['The results were made by an older analyze_sample.m (not support/). Close the ' ...
           'pool with  delete(gcp(''nocreate''))  and run again.']);
end
analysis_time = toc(analysis_start_time);
fprintf('Main loop finished in %.1f s\n', analysis_time);

%% ================= COLLECT RESULTS ==================================================
CI = class_info_v3();
class_names_v3 = CI.names; class_names_v2 = CI.names_v2; class_names_7 = CI.names_7;
nC3 = CI.n;

exact_ok   = [res.exact_ok]';
num_v3     = [res.exact_class]';      num_v3(~exact_ok) = 0;      % 0 = failed
lab_v3     = {res.exact_label}';      lab_v3(~exact_ok) = {'FAILED'};
num_v2     = [res.exact_class_v2]';   num_v2(~exact_ok) = 0;
lab_v2     = {res.exact_label_v2}';   lab_v2(~exact_ok) = {'FAILED'};
num_7      = [res.exact_class_old7]'; num_7(~exact_ok) = 8;       % 8 = failed (as before)
lab_7      = repmat({'FAILED'}, numSamples, 1);
lab_7(exact_ok) = class_names_7(num_7(exact_ok));
lp_eq_bp   = [res.lp_coincides_bp]';
det_in_win = [res.detached_in_window]';
use_ml     = exact_ok & ~lp_eq_bp & num_v3 ~= 11;
fprintf('LP+BP samples whose LP coincides with the BP (flagged, excluded from ML): %d\n', sum(lp_eq_bp));

mc_ok  = [res.mc_ok]';
mc_v3  = [res.mc_class]';    mc_v3(~mc_ok) = 0;
mc_lab_v3 = {res.mc_label}'; mc_lab_v3(~mc_ok) = {'N/A'};
mc_v2  = [res.mc_class_v2]'; mc_v2(~mc_ok) = 0;
mc_lab_v2 = {res.mc_label_v2}'; mc_lab_v2(~mc_ok) = {'N/A'};

vec2str = @(v) mat2str(v(:)', 6);
T = array2table(paramSamples, 'VariableNames', param_names);
T.Exact_Class_v3        = num_v3;
T.Exact_Label_v3        = lab_v3;
T.LP_Coincides_BP       = lp_eq_bp;
T.Use_For_ML            = use_ml;
T.Detached_Fold_In_Window = det_in_win;
T.Exact_Class_v2        = num_v2;
T.Exact_Label_v2        = lab_v2;
T.Multistationary       = [res.multistationary]';
T.Bistable              = [res.bistable]';
T.Max_Steady_States     = [res.max_ss]';
T.Has_Detached_Branch   = [res.has_detached]';
T.Detached_Fold_ATOT    = [res.detached_fold_ATOT]';
T.Symmetric_Degenerate  = [res.degenerate]';
T.LP_ATOT_values        = arrayfun(@(r) vec2str(r.exact_LP_ATOT), res, 'uni', 0);
T.LP_x2_values          = arrayfun(@(r) vec2str(r.exact_LP_x2),   res, 'uni', 0);
T.LP_branch             = arrayfun(@(r) strjoin(r.exact_LP_type, ';'), res, 'uni', 0);
T.LP_nf_a               = arrayfun(@(r) vec2str(r.exact_LP_a),    res, 'uni', 0);
T.LP_nf_beta            = arrayfun(@(r) vec2str(r.exact_LP_beta), res, 'uni', 0);
T.LP_nf_kappa           = arrayfun(@(r) vec2str(r.exact_LP_kappa), res, 'uni', 0);
T.LP_nf_kappa_rel       = arrayfun(@(r) vec2str(r.exact_LP_kappa_rel), res, 'uni', 0);
T.BP_ATOT_values        = arrayfun(@(r) vec2str(r.exact_BP_ATOT), res, 'uni', 0);
T.BP_x2_values          = arrayfun(@(r) vec2str(r.exact_BP_x2),   res, 'uni', 0);
T.BP_nf_a               = arrayfun(@(r) vec2str(r.exact_BP_a),    res, 'uni', 0);
T.BP_nf_b               = arrayfun(@(r) vec2str(r.exact_BP_b),    res, 'uni', 0);
T.BP_nf_c               = arrayfun(@(r) vec2str(r.exact_BP_c),    res, 'uni', 0);
T.BP_nf_disc            = arrayfun(@(r) vec2str(r.exact_BP_disc), res, 'uni', 0);
T.BP_pitchfork          = arrayfun(@(r) vec2str(double(r.exact_BP_pitchfork)), res, 'uni', 0);
T.Bistable_ATOT_ranges  = arrayfun(@(r) mat2str(r.bistable_ATOT, 6), res, 'uni', 0);
T.Multiclass_Label_7    = lab_7;      % backward-compatible 7-class label (exact)
T.Class_Number_7        = num_7;
if cfg.run_matcont
    T.MC_Checked         = doMC;
    T.MC_Class_v3        = mc_v3;
    T.MC_Label_v3        = mc_lab_v3;
    T.MC_Class_v2        = mc_v2;
    T.MC_Label_v2        = mc_lab_v2;
    T.MC_LP_ATOT_values  = arrayfun(@(r) vec2str(r.mc_LP_ATOT), res, 'uni', 0);
    T.MC_LP_nf_a         = arrayfun(@(r) vec2str(r.mc_LP_a),    res, 'uni', 0);
    T.MC_LP_a_reldiff    = arrayfun(@(r) vec2str(r.mc_LP_a_reldiff), res, 'uni', 0);
    T.MC_BP_ATOT_values  = arrayfun(@(r) vec2str(r.mc_BP_ATOT), res, 'uni', 0);
    T.MC_H_ATOT_values   = arrayfun(@(r) vec2str(r.mc_H_ATOT),  res, 'uni', 0);
    T.MC_ATOT_reached    = arrayfun(@(r) vec2str(r.mc_ATOT_reached), res, 'uni', 0);
    T.MC_Negative        = [res.mc_negative]';
    T.MC_Agrees_Exact    = mc_ok & exact_ok & (mc_v3 == num_v3);
end
T.Message = {res.message}';

% one row per singularity, with its normal-form coefficients
NFrows = cell(0, 18);
for k = 1:numSamples
    r = res(k);
    if ~r.exact_ok, continue; end
    for i = 1:numel(r.exact_LP_ATOT)
        aMC = NaN;
        if r.mc_ok && ~isempty(r.mc_LP_ATOT)
            [d, j] = min(abs(r.mc_LP_ATOT - r.exact_LP_ATOT(i)));
            if d <= 1e-3*max(1, r.exact_LP_ATOT(i)), aMC = r.mc_LP_a(j); end
        end
        NFrows(end+1,:) = {k, num_v3(k), lab_v3{k}, 'LP', r.exact_LP_type{i}, r.exact_LP_ATOT(i), ...
            r.exact_LP_x2(i), r.exact_LP_a(i), r.exact_LP_beta(i), r.exact_LP_kappa(i), ...
            r.exact_LP_kappa_rel(i), NaN, NaN, NaN, NaN, NaN, false, aMC}; %#ok<SAGROW>
    end
    for i = 1:numel(r.exact_BP_ATOT)
        NFrows(end+1,:) = {k, num_v3(k), lab_v3{k}, 'BP', 'sym', r.exact_BP_ATOT(i), ...
            r.exact_BP_x2(i), r.exact_BP_a(i), NaN, NaN, NaN, r.exact_BP_b(i), r.exact_BP_c(i), ...
            r.exact_BP_disc(i), r.exact_BP_slopes(i,1), r.exact_BP_slopes(i,2), ...
            r.exact_BP_pitchfork(i), NaN}; %#ok<SAGROW>
    end
end
TNF = cell2table(NFrows, 'VariableNames', {'Sample','Exact_Class_v3','Exact_Label_v3','Kind', ...
    'Branch','ATOT','x2','a','beta','kappa','kappa_rel','b','c','disc','slope1','slope2', ...
    'pitchfork','a_MATCONT'});
fprintf('Normal forms: %d folds, %d branch points\n', sum(strcmp(TNF.Kind,'LP')), sum(strcmp(TNF.Kind,'BP')));
rd = [res.mc_LP_a_reldiff]; rd = rd(isfinite(rd));
if ~isempty(rd)
    fprintf('Fold coefficient |a|, MATCONT vs computed: median rel. difference %.2g (n = %d)\n', median(rd), numel(rd));
end

%% ================= VALIDATION STATISTIC ============================================
V = struct();
fprintf('\n=== Class distribution (exact labels, window ATOT in [%g, %g]) ===\n', cfg.ATOT_window);
counts_v3 = histc(num_v3(exact_ok), 1:nC3); counts_v3 = counts_v3(:);
flag_v3   = histc(num_v3(exact_ok & lp_eq_bp), 1:nC3); flag_v3 = flag_v3(:);
for i = 1:nC3
    fprintf('  %2d %-20s %6d (%5.1f%%)', i, class_names_v3{i}, counts_v3(i), 100*counts_v3(i)/numSamples);
    if flag_v3(i) > 0, fprintf('   [%d with LP = BP, excluded from ML]', flag_v3(i)); end
    fprintf('\n');
end
counts_v2 = histc(num_v2(exact_ok), 1:9); counts_v2 = counts_v2(:);
fprintf('  exact computation failed: %d\n', sum(~exact_ok));
[p_ms, lo_ms, hi_ms] = wilson(sum([res.multistationary]), numSamples);
[p_bs, lo_bs, hi_bs] = wilson(sum([res.bistable]), numSamples);
fprintf('Multistationary: %.1f%% (95%% CI %.1f-%.1f); bistable: %.1f%% (%.1f-%.1f)\n', ...
    100*p_ms, 100*lo_ms, 100*hi_ms, 100*p_bs, 100*lo_bs, 100*hi_bs);
V.class_names = class_names_v3; V.counts_exact = counts_v3; V.counts_lp_eq_bp = flag_v3;
V.class_names_v2 = class_names_v2; V.counts_exact_v2 = counts_v2; V.n_lp_coincides_bp = sum(lp_eq_bp);
V.multistationary = [p_ms lo_ms hi_ms]; V.bistable = [p_bs lo_bs hi_bs];
V.mc_a_median_reldiff = median(rd);

if cfg.run_matcont
    both = exact_ok & mc_ok;
    nb = sum(both);
    [p_ag, lo_ag, hi_ag] = wilson(sum(mc_v3(both) == num_v3(both)), nb);
    ms_exact = [res.multistationary]'; ms_exact = ms_exact(both);
    ms_mc = arrayfun(@(r) ~isempty(r.mc_LP_ATOT) || ~isempty(r.mc_BP_ATOT), res); ms_mc = ms_mc(both);
    [sens, sl, sh] = wilson(sum(ms_mc & ms_exact), sum(ms_exact));
    [spec, pl, ph] = wilson(sum(~ms_mc & ~ms_exact), sum(~ms_exact));
    both = both & num_v3 > 0 & mc_v3 > 0;
    CM = accumarray([num_v3(both), mc_v3(both)], 1, [nC3 nC3]);
    fprintf('\n=== Validation: MATCONT (all branches seeded) vs exact ===\n');
    fprintf('Class agreement: %.2f%% (95%% CI %.2f-%.2f), n = %d\n', 100*p_ag, 100*lo_ag, 100*hi_ag, nb);
    fprintf('Detection of any fold/BP: sensitivity %.2f%% (%.2f-%.2f), specificity %.2f%% (%.2f-%.2f)\n', ...
        100*sens, 100*sl, 100*sh, 100*spec, 100*pl, 100*ph);
    fprintf('MATCONT runs failed: %d of %d checked samples; negative concentrations inside window: %d\n', ...
        sum(~mc_ok & exact_ok & doMC), sum(doMC), sum([res.mc_negative]));
    V.agreement = [p_ag lo_ag hi_ag nb]; V.sensitivity = [sens sl sh]; V.specificity = [spec pl ph];
    V.confusion_rows_exact_cols_matcont = CM;
    Tcm = array2table(CM, 'VariableNames', matlab.lang.makeValidName(strcat('MC_', class_names_v3)), ...
                          'RowNames', strcat('Exact_', class_names_v3));
end

%% ================= SAVE (Excel + CSV + MAT) ========================================
excelFilePath = fullfile(excelDir, 'Bifurcation_Analysis_Results_v3.xlsx');
writetable(T, excelFilePath, 'Sheet', 'Complete_Dataset');
Tsum = table((1:nC3)', class_names_v3', counts_v3, 100*counts_v3/numSamples, flag_v3, counts_v3 - flag_v3, ...
    'VariableNames', {'Class_Number','Class_Label','Count','Percentage','LP_Coincides_BP','Count_For_ML'});
writetable(Tsum, excelFilePath, 'Sheet', 'Class_Summary_v3');
Tsum2 = table((1:9)', class_names_v2', counts_v2, 100*counts_v2/numSamples, ...
    'VariableNames', {'Class_Number','Class_Label','Count','Percentage'});
writetable(Tsum2, excelFilePath, 'Sheet', 'Class_Summary_v2');
writetable(TNF, excelFilePath, 'Sheet', 'Normal_Forms');
writetable(TNF, fullfile(excelDir, sprintf('Normal_Forms_v3_iter%d.csv', ITERATION)));
Tv = table({'Multistationary';'Bistable'}, [p_ms;p_bs]*100, [lo_ms;lo_bs]*100, [hi_ms;hi_bs]*100, ...
    'VariableNames', {'Quantity','Percent','CI95_low','CI95_high'});
if cfg.run_matcont
    Tv = [Tv; table({'MATCONT_vs_exact_class_agreement';'Fold_detection_sensitivity';'Fold_detection_specificity'}, ...
        100*[p_ag;sens;spec], 100*[lo_ag;sl;pl], 100*[hi_ag;sh;ph], ...
        'VariableNames', {'Quantity','Percent','CI95_low','CI95_high'})];
    writetable(Tcm, excelFilePath, 'Sheet', 'Confusion_Exact_vs_MATCONT', 'WriteRowNames', true);
end
writetable(Tv, excelFilePath, 'Sheet', 'Validation');
Tcfg = table({'ITERATION';'numSamples';'ATOT_window';'ATOT_seed';'exact_N';'type';'tot_list'; ...
              'bounds_set';'rngSeed';'matcont_n';'check_n';'check_seed'}, ...
    {ITERATION; numSamples; mat2str(cfg.ATOT_window); cfg.ATOT_seed; cfg.exact_N; type; ...
     strjoin(tot_list,','); bounds_set; rngSeed; cfg.matcont_n; cfg.check_n; cfg.check_seed}, ...
    'VariableNames', {'Setting','Value'});
writetable(Tcfg, excelFilePath, 'Sheet', 'Settings');

% ML-ready dataset: features = sampled parameters, labels = exact classes.
% Use_For_ML = 0 for the flagged LP = BP sets (reported, but left out of the ML).
ML = T(:, [param_names(2:end), {'LP_ATOT_values','LP_x2_values','BP_ATOT_values','BP_x2_values', ...
    'Multistationary','Bistable','Has_Detached_Branch','Detached_Fold_In_Window','LP_Coincides_BP', ...
    'Use_For_ML','Multiclass_Label_7','Class_Number_7','Exact_Label_v2','Exact_Class_v2', ...
    'Exact_Label_v3','Exact_Class_v3'}]);
ML.Iteration = repmat(ITERATION, height(ML), 1);
writetable(ML, fullfile(excelDir, sprintf('ML_Ready_Dataset_v3_iter%d.xlsx', ITERATION)), 'Sheet', 'ML_Ready_Data');
writetable(ML, fullfile(excelDir, sprintf('ML_Ready_Dataset_v3_iter%d.csv', ITERATION)));
writetable(ML(ML.Use_For_ML, :), fullfile(excelDir, sprintf('ML_Ready_Dataset_v3_iter%d_for_ML.csv', ITERATION)));

exact_class_v3 = num_v3; exact_label_v3 = lab_v3; mc_class_v3 = mc_v3; %#ok<NASGU>
save(fullfile(matDir, 'Complete_Bifurcation_Results_v3.mat'), 'res', 'paramSamples', 'param_names', ...
    'bounds', 'cfg', 'type', 'tot_list', 'bounds_set', 'rngSeed', 'ITERATION', 'V', 'T', 'TNF', ...
    'analysis_time', 'exact_class_v3', 'exact_label_v3', 'mc_class_v3', 'class_names_v3', ...
    'lp_eq_bp', 'det_in_win', 'use_ml', 'doMC', '-v7.3');
save(fullfile(matDir, 'Parameter_Samples.mat'), 'paramSamples', 'param_names', 'bounds', 'numSamples');

%% ================= FIGURES (linear axes; .png and .fig) ============================
% 1. class distribution, exact vs MATCONT
fig1 = figure('Color','w','Position',[100 100 1100 450]);
if cfg.run_matcont
    counts_mc = histc(mc_v3(mc_v3 > 0), 1:nC3); counts_mc = counts_mc(:);
    bar([counts_v3 counts_mc]); legend({'Exact','MATCONT (all branches seeded)'}, 'Location','northwest');
else
    bar(counts_v3);
end
set(gca,'XTick',1:nC3,'XTickLabel',class_names_v3,'XTickLabelRotation',40,'FontSize',9,'TickLabelInterpreter','none');
ylabel('Samples'); title(sprintf('Iteration %d: class distribution, ATOT window [%g, %g]', ITERATION, cfg.ATOT_window)); grid on;
saveas(fig1, fullfile(figuresDir, 'Class_Distribution_v3.png'));
savefig(fig1, fullfile(figuresDir, 'Class_Distribution_v3.fig'));

% 2. example branch diagrams per class (solid = stable, dashed = unstable,
%    blue = branch from ATOT = 0, red = detached branch, purple = off-symmetric branch)
fig2 = figure('Color','w','Position',[50 50 1700 950]);
for c = 1:nC3
    ax = subplot(3,4,c); hold(ax,'on'); box(ax,'on');
    ks = find(num_v3 == c & ~lp_eq_bp);
    ks = ks(1:min(cfg.fig_per_class, numel(ks)));
    if ~isempty(ks)
        xl = auto_xlim(res(ks)); ymax = 0;
        for k = ks'
            draw_branches(ax, res(k).curves, 0.9);
            if ~isempty(res(k).exact_LP_ATOT)
                plot(ax, res(k).exact_LP_ATOT, res(k).exact_LP_x2, 'ro', 'MarkerFaceColor','r', 'MarkerSize',4);
            end
            if ~isempty(res(k).exact_BP_ATOT)
                plot(ax, res(k).exact_BP_ATOT, res(k).exact_BP_x2, 'gs', 'MarkerFaceColor','g', 'MarkerSize',5);
            end
            yl = y_range(res(k).curves, xl); ymax = max(ymax, yl(2));
        end
        xlim(ax, xl); ylim(ax, [0 ymax]);
    end
    nf = sum(num_v3 == c & lp_eq_bp);
    t = sprintf('%d %s (n=%d)', c, class_names_v3{c}, sum(num_v3==c));
    if nf > 0, t = sprintf('%s, %d LP=BP not drawn', t, nf); end
    title(ax, t, 'Interpreter','none', 'FontSize', 9);
    xlabel(ax, 'A_{TOT}'); ylabel(ax, 'x_2 (App)');
end
saveas(fig2, fullfile(figuresDir, 'Branches_by_Class_v3.png'));
savefig(fig2, fullfile(figuresDir, 'Branches_by_Class_v3.fig'));

% 3. confusion matrix exact vs MATCONT
if cfg.run_matcont
    fig3 = figure('Color','w','Position',[100 100 700 600]);
    imagesc(log10(1 + CM)); colorbar; axis square;
    set(gca,'XTick',1:nC3,'XTickLabel',class_names_v3,'YTick',1:nC3,'YTickLabel',class_names_v3, ...
        'XTickLabelRotation',40,'FontSize',8,'TickLabelInterpreter','none');
    xlabel('MATCONT class'); ylabel('Exact class');
    for i = 1:nC3, for j = 1:nC3
        if CM(i,j) > 0, text(j, i, num2str(CM(i,j)), 'HorizontalAlignment','center', 'FontSize',8); end
    end, end
    title(sprintf('Agreement %.1f%% (n=%d)', 100*p_ag, nb));
    saveas(fig3, fullfile(figuresDir, 'Confusion_Exact_vs_MATCONT.png'));
    savefig(fig3, fullfile(figuresDir, 'Confusion_Exact_vs_MATCONT.fig'));
end

% 4. normal forms: fold curvature and branch-point coefficients
if height(TNF) > 0
    fig4 = figure('Color','w','Position',[100 100 1200 420]);
    isLP = strcmp(TNF.Kind,'LP'); isBP = ~isLP;
    subplot(1,2,1);
    if any(isLP)
        histogram(sign(TNF.kappa_rel(isLP)).*log10(abs(TNF.kappa_rel(isLP)) + 1), 60);
        xlabel('sign(\kappa_{rel}) log_{10}(1 + |\kappa_{rel}|)'); ylabel('folds');
        title(sprintf('Folds: \\kappa = -a/\\beta, n = %d', sum(isLP)));
    end
    subplot(1,2,2);
    if any(isBP)
        scatter(TNF.b(isBP), TNF.a(isBP), 10, 'filled'); grid on;
        xlabel('b'); ylabel('a');
        title(sprintf('Branch points: n = %d, pitchfork (a = 0): %d', sum(isBP), sum(TNF.pitchfork(isBP))));
    end
    saveas(fig4, fullfile(figuresDir, 'Normal_Forms_v3.png'));
    savefig(fig4, fullfile(figuresDir, 'Normal_Forms_v3.fig'));
end

%% ================= MATCONT CHECK (grey MATCONT over the exact branches) ============
check_time = 0;
if cfg.check_n > 0
    tCheck = tic;
    Ccheck = matcont_check(paramSamples, res, cfg, oldOpt, cfg.check_n, cfg.check_seed, ...
        fullfile(figuresDir, 'MATCONT_check'), useParallel);
    writetable(Ccheck, excelFilePath, 'Sheet', 'MATCONT_check');
    check_time = toc(tCheck);
end

%% ================= CORE HOURS (txt) =================================================
write_core_hours(fullfile(outputDir, sprintf('Core_Hours_iter%d.txt', ITERATION)), ...
    fullfile(HERE, 'Core_Hours_all_iterations.txt'), res, doMC, ITERATION, numSamples, ...
    nw, analysis_time, check_time, cfg.check_n, toc(total_start_time));

fprintf('\nDone in %.1f min. Results in %s\n', toc(total_start_time)/60, outputDir);

%% ================= LOCAL FUNCTIONS =================================================
function [p, lo, hi] = wilson(k, n)
% Wilson score interval, 95%
if n == 0, p = NaN; lo = NaN; hi = NaN; return; end
z = 1.959964; p = k/n;
c = (p + z^2/(2*n)) / (1 + z^2/n);
h = z*sqrt(p*(1-p)/n + z^2/(4*n^2)) / (1 + z^2/n);
lo = max(0, c - h); hi = min(1, c + h);
end

function write_core_hours(fname, fall, res, doMC, ITERATION, numSamples, nw, tLoop, tCheck, nCheck, tTotal)
% Core-hours of one iteration. Per-sample times are measured inside analyze_sample
% (tic/toc on the worker, one core per worker), so their sum is the serial CPU time
% in core-hours and does not depend on the number of workers. The allocated
% core-hours (wall clock x workers) are given as well.
ok   = [res.exact_ok];
tEx  = [res.t_exact];  tNF = [res.t_nf];  tMC = [res.t_mc];
cEx  = [res.cpu_exact]; cMC = [res.cpu_mc];
mc   = doMC(:)' & [res.mc_ran];
h    = @(t) sum(t(isfinite(t)))/3600;
cores = max(1, nw);
[cpuName, ramGB] = machine_info();
L = {};
L{end+1} = sprintf('Ensemble_LHS_v3 - core hours, ITERATION %d', ITERATION);
L{end+1} = sprintf('Date: %s', datestr(now, 'yyyy-mm-dd HH:MM'));
L{end+1} = sprintf('Machine: %s | CPU: %s | RAM: %s GB | MATLAB %s', ...
    strtrim(computer), cpuName, ramGB, version);
try, nc = feature('numcores'); catch, nc = NaN; end
L{end+1} = sprintf('Workers: %d (physical cores visible to MATLAB: %d)', nw, nc);
L{end+1} = '';
L{end+1} = sprintf('Samples: %d (exact OK: %d), MATCONT-validated: %d', numSamples, sum(ok), sum(mc));
L{end+1} = '';
L{end+1} = 'CORE-HOURS (sum of per-sample times = serial CPU time)';
L{end+1} = sprintf('  Exact labelling (inverse method, all samples)   : %10.4f core-h  (%.4f s/sample mean, %.4f median)', ...
    h(tEx), mean(tEx(isfinite(tEx))), median(tEx(isfinite(tEx))));
L{end+1} = sprintf('    of which labels + normal forms               : %10.6f core-h  (%.2f ms/sample mean)', ...
    h(tNF), 1e3*mean(tNF(isfinite(tNF))));
if any(mc)
    L{end+1} = sprintf('  MATCONT validation (seeded on every branch)     : %10.4f core-h  (%.3f s/validated sample mean, %.3f median)', ...
        h(tMC), mean(tMC(mc)), median(tMC(mc)));
end
L{end+1} = sprintf('  TOTAL main loop                                 : %10.4f core-h', h(tEx) + h(tMC));
L{end+1} = sprintf('  (CPU time reported by cputime: exact %.4f, MATCONT %.4f core-h)', h(cEx), h(cMC));
L{end+1} = '';
L{end+1} = 'WALL CLOCK';
L{end+1} = sprintf('  Main loop        : %10.1f s = %.3f h  -> allocated %.4f core-h (%d workers)', tLoop, tLoop/3600, tLoop*cores/3600, cores);
L{end+1} = sprintf('  MATCONT check    : %10.1f s (%d samples, seeded + single-start runs) -> allocated %.4f core-h', tCheck, nCheck, tCheck*cores/3600);
L{end+1} = sprintf('  Whole script     : %10.1f s = %.3f h  -> allocated %.4f core-h', tTotal, tTotal/3600, tTotal*cores/3600);
txt = strjoin(L, newline);
fid = fopen(fname, 'w'); fprintf(fid, '%s\n', txt); fclose(fid);
fprintf('\n%s\nSaved %s\n', txt, fname);
% one line per run in a combined file (all iterations)
newFile = ~exist(fall, 'file');
fid = fopen(fall, 'a');
if newFile
    fprintf(fid, 'date\titeration\tsamples\tworkers\texact_core_h\tnf_core_h\tmatcont_core_h\tmatcont_n\ttotal_core_h\tloop_wall_h\tscript_wall_h\tCPU\n');
end
fprintf(fid, '%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%d\t%.4f\t%.4f\t%.4f\t%s\n', datestr(now,'yyyy-mm-dd HH:MM'), ...
    ITERATION, numSamples, nw, h(tEx), h(tNF), h(tMC), sum(mc), h(tEx)+h(tMC), tLoop/3600, tTotal/3600, cpuName);
fclose(fid);
end

function [cpuName, ramGB] = machine_info()
cpuName = 'unknown'; ramGB = 'unknown';
try
    if ispc
        cpuName = strtrim(getenv('PROCESSOR_IDENTIFIER'));
        [~, out] = system('wmic cpu get name /value');
        t = regexp(out, 'Name=([^\r\n]+)', 'tokens', 'once'); if ~isempty(t), cpuName = strtrim(t{1}); end
        [~, m] = memory; ramGB = sprintf('%.0f', m.PhysicalMemory.Total/2^30);
    elseif isunix && ~ismac
        [~, out] = system('grep -m1 "model name" /proc/cpuinfo');
        t = regexp(out, ':\s*(.+)', 'tokens', 'once'); if ~isempty(t), cpuName = strtrim(t{1}); end
        [~, out] = system('grep MemTotal /proc/meminfo');
        kb = sscanf(regexprep(out, '[^0-9]', ' '), '%f', 1); if ~isempty(kb), ramGB = sprintf('%.0f', kb/2^20); end
    elseif ismac
        [~, out] = system('sysctl -n machdep.cpu.brand_string'); cpuName = strtrim(out);
        [~, out] = system('sysctl -n hw.memsize'); ramGB = sprintf('%.0f', str2double(out)/2^30);
    end
catch
end
end

function on_result(d, cfg, N)
% Runs on the client for every finished sample: progress + live figures (linear axes).
persistent n t0 axC cnt figF axF nextF names xmax ymax
if isempty(d)                                   % reset
    CI = class_info_v3(); names = CI.names;
    n = 0; t0 = tic; cnt = zeros(1,numel(names)); nextF = 1;
    xmax = zeros(1,numel(names)); ymax = zeros(1,numel(names));
    if cfg.live_plot
        f1 = figure('Name','Live: all samples by class','Color','w','Position',[20 60 1500 850]);
        axC = gobjects(1,numel(names));
        for c = 1:numel(names)
            axC(c) = subplot(3,4,c,'Parent',f1); hold(axC(c),'on'); box(axC(c),'on');
            title(axC(c), sprintf('%s (n=0)', names{c}), 'Interpreter','none', 'FontSize', 9);
            xlabel(axC(c),'A_{TOT}'); ylabel(axC(c),'x_2');
        end
        figF = figure('Name','Live: latest flagged samples','Color','w','Position',[80 40 1500 850]);
        axF = gobjects(1,12);
        for j = 1:12, axF(j) = subplot(3,4,j,'Parent',figF); end
    end
    return
end
n = n + 1;
R = d.R;
if mod(n, max(1, round(N/100))) == 0 || n == N
    fprintf('  %d / %d samples (%.0f s)\n', n, N, toc(t0));
end
if ~cfg.live_plot || ~R.exact_ok, return; end
c = R.exact_class;
if ~(c >= 1 && c <= numel(names)), return; end
cnt(c) = cnt(c) + 1;
if isvalid(axC(c))
    if cnt(c) <= cfg.live_max_per_class && ~R.lp_coincides_bp
        draw_branches(axC(c), R.curves, 0.7);
        if ~isempty(R.exact_LP_ATOT)
            plot(axC(c), R.exact_LP_ATOT, R.exact_LP_x2, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 3);
        end
        if ~isempty(R.exact_BP_ATOT)
            plot(axC(c), R.exact_BP_ATOT, R.exact_BP_x2, 'gs', 'MarkerFaceColor','g', 'MarkerSize', 4);
        end
        xl = auto_xlim(R); yl = y_range(R.curves, xl);
        xmax(c) = max(xmax(c), xl(2)); ymax(c) = max(ymax(c), yl(2));
        xlim(axC(c), [0 xmax(c)]); ylim(axC(c), [0 ymax(c)]);
    end
    title(axC(c), sprintf('%s (n=%d)', names{c}, cnt(c)), 'Interpreter','none', 'FontSize', 9);
end
if cfg.live_flag(R) && isvalid(figF)
    plot_sample(R, d.P, d.k, axF(nextF));
    nextF = mod(nextF, 12) + 1;
end
drawnow limitrate
end
