% ======================= Ensemble_LHS_v2.m =========================================
% Bifurcation-class labelling of the dual futile cycle (MAPKfull.m) over an LHS
% ensemble, relabelled with the exact (inverse-method) structure.
%
% What changed relative to Ensemble_LHS.m
%  1. Labels come from the EXACT steady-state structure (ss_structure.m): every
%     steady state for every ATOT in a fixed window, computed from an explicit 1-D
%     parametrisation (Wang & Sontag 2008). No initial condition, no ode15s, no grid
%     in ATOT; folds and branch points are located exactly; stability from eigenvalues.
%  2. MATCONT is kept as an independent check, but it is now seeded on EVERY branch of
%     steady states (including detached branches that cannot be reached from low
%     ATOT), run in both directions, and trimmed to the same ATOT window.
%     The agreement between the two labellings is reported as the validation statistic.
%  3. The ATOT window is the same for all samples (it was set implicitly by
%     MaxStepSize*MaxNumPoints = 1.5 units of arclength from ATOT = 0.1 before).
%  4. The old Class 7 ("OTHER") is split: 7 LP_DETACHED (new), 8 NO_SING_monotone,
%     9 OTHER. Class_Number_7 is still written (7/8/9 -> 7) for backward comparison.
%  5. Bugs fixed: FAILED samples wrote their label into the wrong row
%     (category{current_category_index}); any value < 0 (even round-off) sent a whole
%     sample to class 7; Hopf was tested as 'HB' (MATCONT writes 'H '); the ML column
%     'ATOT_LP' held concentrations, not ATOT values.
%  6. The per-sample work runs in a parfor loop; a fixed RNG seed makes it reproducible.
% ===================================================================================
clear; close all; clc

%% ================= USER SETTINGS ====================================================
MATCONT_DIR  = '';          % folder with MATCONT's init.m ('' = already on the path)
PROJECT_DIR  = pwd;         % folder with MAPKfull.m and the helper .m files
outputDir    = 'Ensemble_open_v2_1e6';
numSamples   = 100000;
rngSeed      = 20260930;    % reproducible LHS design
NumWorkers   = 64;          % 0 = run serially

type     = 1;               % 1 = free sampling, 2 = kinase/phosphatase symmetry imposed
tot_list = {};              % type 2 only: parameters sampled independently to break
                            % symmetry, e.g. {'k5','k55'}; {} = fully symmetric
% bounds_set = 'S_dominant';  % 'S_dominant' | 'Z_dominant' | 'wide'
bounds_set = 'wide';  % 'S_dominant' | 'Z_dominant' | 'wide'

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
cfg.save_png_if        = @(R) false;             % per-sample PNG into Figures/Samples, e.g.
                                                 % @(R) R.multistationary  or  @(R) true

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
if cfg.run_matcont
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
end

%% ================= PARALLEL POOL ====================================================
useParallel = NumWorkers > 0 && license('test','Distrib_Computing_Toolbox');
if useParallel
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool('local', NumWorkers); end
    wait(parfevalOnAll(pool, @addpath, 0, PROJECT_DIR));
    if cfg.run_matcont
        wait(parfevalOnAll(pool, @mc_init, 0, MATCONT_DIR, PROJECT_DIR));
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
end
res = [resC{:}]';
analysis_time = toc(analysis_start_time);
fprintf('Main loop finished in %.1f s\n', analysis_time);

%% ================= COLLECT RESULTS ==================================================
class_names_v2 = {'LP_ONLY_S-shaped','LP_ONLY_Z-shaped','BP_ONLY','LP_BP_LPgreaterBP', ...
    'LP_BP_LPlessBP','NO_SING_biphasic','LP_DETACHED','NO_SING_monotone','OTHER'};
class_names_7 = [class_names_v2(1:6), {'OTHER'}];

exact_ok   = [res.exact_ok]';
exact_cls  = [res.exact_class]';
exact_old7 = [res.exact_class_old7]';
mc_ok      = [res.mc_ok]';
mc_cls     = [res.mc_class]';

lab_v2 = repmat({'FAILED'}, numSamples, 1);
lab_v2(exact_ok) = class_names_v2(exact_cls(exact_ok));
lab_7 = repmat({'FAILED'}, numSamples, 1);
lab_7(exact_ok) = class_names_7(exact_old7(exact_ok));
num_v2 = exact_cls;  num_v2(~exact_ok) = 0;        % 0 = failed
num_7  = exact_old7; num_7(~exact_ok) = 8;         % 8 = failed (as before)

mc_lab = repmat({'N/A'}, numSamples, 1);
mc_lab(mc_ok) = class_names_v2(mc_cls(mc_ok));

vec2str = @(v) mat2str(v(:)', 6);
T = array2table(paramSamples, 'VariableNames', param_names);
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
T.BP_ATOT_values        = arrayfun(@(r) vec2str(r.exact_BP_ATOT), res, 'uni', 0);
T.BP_x2_values          = arrayfun(@(r) vec2str(r.exact_BP_x2),   res, 'uni', 0);
T.Bistable_ATOT_ranges  = arrayfun(@(r) mat2str(r.bistable_ATOT, 6), res, 'uni', 0);
T.Multiclass_Label_7    = lab_7;      % backward-compatible 7-class label (exact)
T.Class_Number_7        = num_7;
if cfg.run_matcont
    T.MC_Class_v2        = mc_cls;
    T.MC_Label_v2        = mc_lab;
    T.MC_LP_ATOT_values  = arrayfun(@(r) vec2str(r.mc_LP_ATOT), res, 'uni', 0);
    T.MC_BP_ATOT_values  = arrayfun(@(r) vec2str(r.mc_BP_ATOT), res, 'uni', 0);
    T.MC_H_ATOT_values   = arrayfun(@(r) vec2str(r.mc_H_ATOT),  res, 'uni', 0);
    T.MC_ATOT_reached    = arrayfun(@(r) vec2str(r.mc_ATOT_reached), res, 'uni', 0);
    T.MC_Negative        = [res.mc_negative]';
    T.MC_Agrees_Exact    = mc_ok & exact_ok & (mc_cls == exact_cls);
end
T.Message = {res.message}';

%% ================= VALIDATION STATISTIC ============================================
V = struct();
fprintf('\n=== Class distribution (exact labels, window ATOT in [%g, %g]) ===\n', cfg.ATOT_window);
counts_v2 = histc(num_v2(exact_ok), 1:9); counts_v2 = counts_v2(:);
for i = 1:9
    fprintf('  %d %-20s %6d (%5.1f%%)\n', i, class_names_v2{i}, counts_v2(i), 100*counts_v2(i)/numSamples);
end
fprintf('  exact computation failed: %d\n', sum(~exact_ok));
[p_ms, lo_ms, hi_ms] = wilson(sum([res.multistationary]), numSamples);
[p_bs, lo_bs, hi_bs] = wilson(sum([res.bistable]), numSamples);
fprintf('Multistationary: %.1f%% (95%% CI %.1f-%.1f); bistable: %.1f%% (%.1f-%.1f)\n', ...
    100*p_ms, 100*lo_ms, 100*hi_ms, 100*p_bs, 100*lo_bs, 100*hi_bs);
V.class_names = class_names_v2; V.counts_exact = counts_v2;
V.multistationary = [p_ms lo_ms hi_ms]; V.bistable = [p_bs lo_bs hi_bs];

if cfg.run_matcont
    both = exact_ok & mc_ok;
    nb = sum(both);
    [p_ag, lo_ag, hi_ag] = wilson(sum(mc_cls(both) == exact_cls(both)), nb);
    ms_exact = [res.multistationary]'; ms_exact = ms_exact(both);
    ms_mc = arrayfun(@(r) ~isempty(r.mc_LP_ATOT) || ~isempty(r.mc_BP_ATOT), res); ms_mc = ms_mc(both);
    [sens, sl, sh] = wilson(sum(ms_mc & ms_exact), sum(ms_exact));
    [spec, pl, ph] = wilson(sum(~ms_mc & ~ms_exact), sum(~ms_exact));
    CM = accumarray([exact_cls(both), mc_cls(both)], 1, [9 9]);
    fprintf('\n=== Validation: MATCONT (all branches seeded) vs exact ===\n');
    fprintf('Class agreement: %.2f%% (95%% CI %.2f-%.2f), n = %d\n', 100*p_ag, 100*lo_ag, 100*hi_ag, nb);
    fprintf('Detection of any fold/BP: sensitivity %.2f%% (%.2f-%.2f), specificity %.2f%% (%.2f-%.2f)\n', ...
        100*sens, 100*sl, 100*sh, 100*spec, 100*pl, 100*ph);
    fprintf('MATCONT runs failed: %d of %d checked samples; negative concentrations inside window: %d\n', ...
        sum(~mc_ok & exact_ok & doMC), sum(doMC), sum([res.mc_negative]));
    V.agreement = [p_ag lo_ag hi_ag nb]; V.sensitivity = [sens sl sh]; V.specificity = [spec pl ph];
    V.confusion_rows_exact_cols_matcont = CM;
    Tcm = array2table(CM, 'VariableNames', matlab.lang.makeValidName(strcat('MC_', class_names_v2)), ...
                          'RowNames', strcat('Exact_', class_names_v2));
end

%% ================= SAVE (Excel + CSV + MAT) ========================================
excelFilePath = fullfile(excelDir, 'Bifurcation_Analysis_Results.xlsx');
writetable(T, excelFilePath, 'Sheet', 'Complete_Dataset');
Tsum = table((1:9)', class_names_v2', counts_v2, 100*counts_v2/numSamples, ...
    'VariableNames', {'Class_Number','Class_Label','Count','Percentage'});
writetable(Tsum, excelFilePath, 'Sheet', 'Class_Summary_v2');
Tv = table({'Multistationary';'Bistable'}, [p_ms;p_bs]*100, [lo_ms;lo_bs]*100, [hi_ms;hi_bs]*100, ...
    'VariableNames', {'Quantity','Percent','CI95_low','CI95_high'});
if cfg.run_matcont
    Tv = [Tv; table({'MATCONT_vs_exact_class_agreement';'Fold_detection_sensitivity';'Fold_detection_specificity'}, ...
        100*[p_ag;sens;spec], 100*[lo_ag;sl;pl], 100*[hi_ag;sh;ph], ...
        'VariableNames', {'Quantity','Percent','CI95_low','CI95_high'})];
    writetable(Tcm, excelFilePath, 'Sheet', 'Confusion_Exact_vs_MATCONT', 'WriteRowNames', true);
end
writetable(Tv, excelFilePath, 'Sheet', 'Validation');
Tcfg = table({'ATOT_window';'ATOT_seed';'exact_N';'type';'tot_list';'bounds_set';'rngSeed'}, ...
    {mat2str(cfg.ATOT_window); cfg.ATOT_seed; cfg.exact_N; type; strjoin(tot_list,','); bounds_set; rngSeed}, ...
    'VariableNames', {'Setting','Value'});
writetable(Tcfg, excelFilePath, 'Sheet', 'Settings');

% ML-ready dataset: features = sampled parameters, labels = exact classes
ML = T(:, [param_names(2:end), {'LP_ATOT_values','LP_x2_values','BP_ATOT_values','BP_x2_values', ...
    'Multistationary','Bistable','Has_Detached_Branch','Multiclass_Label_7','Class_Number_7', ...
    'Exact_Label_v2','Exact_Class_v2'}]);
writetable(ML, fullfile(excelDir, 'ML_Ready_Dataset.xlsx'), 'Sheet', 'ML_Ready_Data');
writetable(ML, fullfile(excelDir, 'ML_Ready_Dataset.csv'));

save(fullfile(matDir, 'Complete_Bifurcation_Results_v2.mat'), 'res', 'paramSamples', 'param_names', ...
    'bounds', 'cfg', 'type', 'tot_list', 'bounds_set', 'rngSeed', 'V', 'T', 'analysis_time', '-v7.3');
save(fullfile(matDir, 'Parameter_Samples.mat'), 'paramSamples', 'param_names', 'bounds', 'numSamples');

%% ================= FIGURES =========================================================
% 1. class distribution, exact vs MATCONT
fig1 = figure('Color','w','Position',[100 100 1100 450]);
if cfg.run_matcont
    counts_mc = histc(mc_cls(mc_ok), 1:9); counts_mc = counts_mc(:);
    bar([counts_v2 counts_mc]); legend({'Exact','MATCONT (all branches)'}, 'Location','northwest');
else
    bar(counts_v2);
end
set(gca,'XTick',1:9,'XTickLabel',class_names_v2,'XTickLabelRotation',40,'FontSize',9);
ylabel('Samples'); title(sprintf('Class distribution, ATOT window [%g, %g]', cfg.ATOT_window)); grid on;
saveas(fig1, fullfile(figuresDir, 'Class_Distribution_v2.png'));

% 2. example branch diagrams per class (solid = stable, dashed = unstable,
%    blue = branch from ATOT = 0, red = detached branch, purple = off-symmetric branch)
fig2 = figure('Color','w','Position',[50 50 1500 850]);
for c = 1:9
    subplot(3,3,c); hold on; box on;
    ks = find(num_v2 == c);
    ks = ks(1:min(25, numel(ks)));
    for k = ks'
        draw_branches(gca, res(k).curves, 0.8);
        if ~isempty(res(k).exact_LP_ATOT)
            plot(res(k).exact_LP_ATOT, res(k).exact_LP_x2, 'ro', 'MarkerFaceColor','r', 'MarkerSize',4);
        end
        if ~isempty(res(k).exact_BP_ATOT)
            plot(res(k).exact_BP_ATOT, res(k).exact_BP_x2, 'gs', 'MarkerFaceColor','g', 'MarkerSize',5);
        end
    end
    set(gca,'XScale','log','YScale','log'); xlim(max(cfg.ATOT_window, [1e-3 0]));
    title(sprintf('%s (n=%d)', class_names_v2{c}, sum(num_v2==c)), 'Interpreter','none');
    xlabel('ATOT'); ylabel('x_2 (App)');
end
saveas(fig2, fullfile(figuresDir, 'Branches_by_Class_v2.png'));

% 3. confusion matrix exact vs MATCONT
if cfg.run_matcont
    fig3 = figure('Color','w','Position',[100 100 700 600]);
    imagesc(log10(1 + CM)); colorbar; axis square;
    set(gca,'XTick',1:9,'XTickLabel',class_names_v2,'YTick',1:9,'YTickLabel',class_names_v2, ...
        'XTickLabelRotation',40,'FontSize',8,'TickLabelInterpreter','none');
    xlabel('MATCONT class'); ylabel('Exact class');
    for i = 1:9, for j = 1:9
        if CM(i,j) > 0, text(j, i, num2str(CM(i,j)), 'HorizontalAlignment','center', 'FontSize',8); end
    end, end
    title(sprintf('Agreement %.1f%% (n=%d)', 100*p_ag, nb));
    saveas(fig3, fullfile(figuresDir, 'Confusion_Exact_vs_MATCONT.png'));
end

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

function on_result(d, cfg, N)
% Runs on the client for every finished sample: progress + live figures.
persistent n t0 axC cnt figF axF nextF names
if isempty(d)                                   % reset
    n = 0; t0 = tic; cnt = zeros(1,9); nextF = 1;
    names = {'LP_ONLY_S-shaped','LP_ONLY_Z-shaped','BP_ONLY','LP_BP_LPgreaterBP', ...
             'LP_BP_LPlessBP','NO_SING_biphasic','LP_DETACHED','NO_SING_monotone','OTHER'};
    if cfg.live_plot
        f1 = figure('Name','Live: all samples by class','Color','w','Position',[20 60 1250 800]);
        axC = gobjects(1,9);
        for c = 1:9
            axC(c) = subplot(3,3,c,'Parent',f1); hold(axC(c),'on'); box(axC(c),'on');
            set(axC(c),'XScale','log','YScale','log');
            title(axC(c), sprintf('%s (n=0)', names{c}), 'Interpreter','none', 'FontSize', 9);
            xlabel(axC(c),'ATOT'); ylabel(axC(c),'x_2');
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
cnt(c) = cnt(c) + 1;
if isvalid(axC(c))
    if cnt(c) <= cfg.live_max_per_class
        draw_branches(axC(c), R.curves, 0.7);
        if ~isempty(R.exact_LP_ATOT)
            plot(axC(c), R.exact_LP_ATOT, R.exact_LP_x2, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 3);
        end
        if ~isempty(R.exact_BP_ATOT)
            plot(axC(c), R.exact_BP_ATOT, R.exact_BP_x2, 'gs', 'MarkerFaceColor','g', 'MarkerSize', 4);
        end
    end
    title(axC(c), sprintf('%s (n=%d)', names{c}, cnt(c)), 'Interpreter','none', 'FontSize', 9);
end
if cfg.live_flag(R) && isvalid(figF)
    plot_sample(R, d.P, d.k, axF(nextF));
    nextF = mod(nextF, 12) + 1;
end
drawnow limitrate
end
