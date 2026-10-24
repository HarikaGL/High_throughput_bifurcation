% ======================= Compare_SingleStart_vs_Exact.m ==============================
% Runs BOTH labelling methods on the SAME LHS parameter sets (same seed) and compares.
%
%   OLD : one steady state (ode15s) -> single MATCONT continuation (+ BP extensions)
%         -> 7 classes                       (old_analyze_sample.m)
%   NEW : exact steady-state structure (ss_structure via analyze_sample.m) -> 9 classes
%
% The two methods do not look at the same ATOT range: the old continuation only covers
% what 1500 points x MaxStepSize 1e-3 reach from ATOT = 0.1, the new one covers
% cfg.ATOT_window. So the new method is run twice:
%   NEW (full window)    : exactly as in Ensemble_LHS_v2 (what you will publish)
%   NEW (matched window) : per sample, window = [0, largest ATOT the old method reached]
%                          -> isolates the METHOD difference (single start vs all branches)
% Old vs NEW(matched) = what the single-start strategy misses (detached branches etc.)
% Old vs NEW(full)    = the above plus the effect of the wider ATOT window.
%
% Needs on the path: MAPKfull.m, mc_init.m, analyze_sample.m, ss_structure.m,
% draw_branches.m (all from Ensemble_LHS_v2) and old_analyze_sample.m (this folder).
% The MATCONT-vs-exact validation inside v2 is switched off here (it is not the point).
% =====================================================================================
clear; close all; clc

%% ================= USER SETTINGS ====================================================
MATCONT_DIR = '';                 % folder with MATCONT's init.m ('' = already on path)
PROJECT_DIR = pwd;
outputDir   = 'Compare_Old_vs_New';
numSamples  = 10000;               % old method is slow (several s/sample/core): start small
rngSeed     = 20260930;           % SAME seed/design for both methods
NumWorkers  = 64;                 % 0 = serial
type        = 1;                  % 1 free, 2 kinase/phosphatase symmetry
tot_list    = {};                 % type 2 only
bounds_set  = 'S_dominant';             % 'S_dominant' | 'Z_dominant' | 'wide'
run_new_matched = true;           % also run the new method in the old method's ATOT range
maxPerClass = 40;                 % curves overlaid per class panel

cfg = struct();                   % settings of the NEW method (as in Ensemble_LHS_v2)
cfg.ATOT_window    = [0 1e6];     % [0 1e6] to match the v2 production run
cfg.ATOT_seed      = 0.1;
cfg.exact_N        = 6000;
cfg.run_matcont    = false;
cfg.matcont_n      = 0;
cfg.keep_mc_curves = false;
cfg.neg_tol        = 1e-8;
cfg.curve_points   = 600;
cfg.live_plot = false; cfg.live_max_per_class = 40;
cfg.live_flag = @(R) false; cfg.save_png_if = @(R) false;

xl_plot = [1e-3 cfg.ATOT_window(2)];   % x-limits of the ensemble figures (log axis)

%% ================= OUTPUT FOLDERS ===================================================
figuresDir = fullfile(outputDir,'Figures');
excelDir   = fullfile(outputDir,'Excel_Files');
matDir     = fullfile(outputDir,'MAT_Files');
for d = {outputDir, figuresDir, excelDir, matDir}
    if ~exist(d{1},'dir'), mkdir(d{1}); end
end
total_start = tic;

%% ================= PARAMETERS AND LHS (identical to Ensemble_LHS_v2) ================
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
paramSamples = repmat(Param_Base_Vec, numSamples, 1);
logsample = @(name, col) 10.^(log10(bounds.(name)(1)) + ...
    lhsNorm(:,col)*(log10(bounds.(name)(2)) - log10(bounds.(name)(1))));
idx = @(name) find(strcmp(param_names, name));
if type == 1
    for p = 2:numParams
        paramSamples(:,p) = logsample(param_names{p}, p);
    end
else
    pairs = {'KTOT','PTOT'; 'k1','k5'; 'k11','k55'; 'k2','k6'; 'k3','k7'; 'k33','k77'; 'k4','k8'};
    for i = 1:size(pairs,1)
        v = logsample(pairs{i,1}, idx(pairs{i,1}));
        paramSamples(:, idx(pairs{i,1})) = v;
        paramSamples(:, idx(pairs{i,2})) = v;
    end
    for i = 1:numel(tot_list)
        paramSamples(:, idx(tot_list{i})) = logsample(tot_list{i}, idx(tot_list{i}));
    end
end
fprintf('Generated %d parameter sets (type %d, %s bounds, seed %d) - shared by both methods\n', ...
    numSamples, type, bounds_set, rngSeed);

%% ================= MATCONT SETUP AND POOL ===========================================
mc_init(MATCONT_DIR, PROJECT_DIR);
oldOpt = contset;                                   % legacy options, unchanged
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

useParallel = NumWorkers > 0 && license('test','Distrib_Computing_Toolbox');
if useParallel
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool('local', NumWorkers); end
    wait(parfevalOnAll(pool, @addpath, 0, PROJECT_DIR));
    wait(parfevalOnAll(pool, @mc_init, 0, MATCONT_DIR, PROJECT_DIR));
    nw = pool.NumWorkers;
else
    nw = 0;
end

%% ================= RUN 1: OLD METHOD ================================================
fprintf('\n=== OLD method: %d samples ===\n', numSamples);
t0 = tic;
resOldC = run_parallel(@(k) old_analyze_sample(paramSamples(k,:), oldOpt, cfg.curve_points), ...
    numSamples, nw, 'old');
tOld = toc(t0);
resOld = [resOldC{:}]';  clear resOldC

%% ================= RUN 2: NEW METHOD, FULL WINDOW ===================================
fprintf('\n=== NEW method, window [%g %g] ===\n', cfg.ATOT_window);
t0 = tic;
resNewC = run_parallel(@(k) analyze_sample(paramSamples(k,:), cfg, false), numSamples, nw, 'new full');
tNew = toc(t0);
resNew = [resNewC{:}]';  clear resNewC

%% ================= RUN 3: NEW METHOD, WINDOW MATCHED TO THE OLD RUN =================
wMatch = arrayfun(@(r) r.ATOT_reached, resOld);
wMatch(isnan(wMatch)) = median(wMatch(~isnan(wMatch)));
wMatch = max(wMatch, 1.01*cfg.ATOT_seed);
tMatch = NaN;
if run_new_matched
    fprintf('\n=== NEW method, per-sample window = ATOT range reached by OLD ===\n');
    cfgM = cell(numSamples,1);
    for k = 1:numSamples
        c = cfg; c.ATOT_window = [0 wMatch(k)]; cfgM{k} = c;
    end
    t0 = tic;
    resMC = run_parallel(@(k) analyze_sample(paramSamples(k,:), cfgM{k}, false), numSamples, nw, 'new matched');
    tMatch = toc(t0);
    resMatch = [resMC{:}]';  clear resMC
end

%% ================= COLLECT LABELS ===================================================
names7 = {'LP_ONLY_S-shaped','LP_ONLY_Z-shaped','BP_ONLY','LP_BP_LPgreaterBP', ...
          'LP_BP_LPlessBP','NO_SING_biphasic','OTHER','FAILED/NEG'};
names9 = {'LP_ONLY_S-shaped','LP_ONLY_Z-shaped','BP_ONLY','LP_BP_LPgreaterBP', ...
          'LP_BP_LPlessBP','NO_SING_biphasic','LP_DETACHED','NO_SING_monotone','OTHER'};

old_cls    = [resOld.class7]';
old_status = {resOld.status}';
old_neg    = [resOld.negative]';
oldFold    = arrayfun(@(r) ~isempty(r.LP_ATOT), resOld);
oldBP      = arrayfun(@(r) ~isempty(r.BP_ATOT), resOld);

[new7F, new9F, okF] = new_labels(resNew);
if run_new_matched, [new7M, ~, okM] = new_labels(resMatch); else, new7M = new7F; okM = okF; end %#ok<NASGU>

msF   = [resNew.multistationary]';
detF  = [resNew.has_detached]';
newFoldF = arrayfun(@(r) ~isempty(r.exact_LP_ATOT), resNew);

cOld  = accumarray(old_cls, 1, [8 1]);
cNewF = accumarray(new7F,   1, [8 1]);
cNewM = accumarray(new7M,   1, [8 1]);
c9    = accumarray(max(new9F,1), double(new9F>0), [9 1]);
nFail9 = sum(new9F==0);

%% ================= STATISTICS =======================================================
CM_F = accumarray([old_cls, new7F], 1, [8 8]);     % rows OLD, cols NEW
CM_M = accumarray([old_cls, new7M], 1, [8 8]);
[agF, agF_lo, agF_hi] = wilson(sum(old_cls == new7F), numSamples);
[agM, agM_lo, agM_hi] = wilson(sum(old_cls == new7M), numSamples);

oldAny = oldFold | oldBP;
[sens, sens_lo, sens_hi] = wilson(sum(oldAny & msF), sum(msF));
[spec, spec_lo, spec_hi] = wilson(sum(~oldAny & ~msF), sum(~msF));
[detAg, ~, ~] = wilson(sum(old_cls(detF) == new7F(detF)), max(1,sum(detF)));

fprintf('\n==================== RESULTS (n = %d, seed %d) ====================\n', numSamples, rngSeed);
fprintf('Run time: old %.1f s | new full %.1f s | new matched %.1f s  (%d workers)\n', tOld, tNew, tMatch, nw);
fprintf('\n%-3s %-20s %8s %8s %8s\n', '#', 'class (7+1 bins)', 'OLD', 'NEWmatch', 'NEWfull');
for i = 1:8
    fprintf('%-3d %-20s %8d %8d %8d\n', i, names7{i}, cOld(i), cNewM(i), cNewF(i));
end
fprintf('\nOLD: failed %d, negative-concentration %d (counted in bin 8 here; the legacy 7-class count put them in OTHER)\n', ...
    sum(strcmp(old_status,'FAILED')), sum(old_neg));
fprintf('NEW full: failed %d | multistationary %d (%.1f%%) | detached branch %d (%.1f%%) | LP_DETACHED class %d\n', ...
    nFail9, sum(msF), 100*mean(msF), sum(detF), 100*mean(detF), c9(7));
fprintf('\nClass agreement OLD vs NEW matched window : %.2f%% (95%% CI %.2f-%.2f)\n', 100*agM, 100*agM_lo, 100*agM_hi);
fprintf('Class agreement OLD vs NEW full window    : %.2f%% (95%% CI %.2f-%.2f)\n', 100*agF, 100*agF_lo, 100*agF_hi);
fprintf('OLD detects a fold/BP in multistationary samples (NEW full as truth): sensitivity %.2f%% (%.2f-%.2f), specificity %.2f%% (%.2f-%.2f)\n', ...
    100*sens, 100*sens_lo, 100*sens_hi, 100*spec, 100*spec_lo, 100*spec_hi);
fprintf('Samples with a detached branch: %d; OLD gives the same class in %.1f%% of them\n', sum(detF), 100*detAg);

off = CM_M - diag(diag(CM_M));
[vals, ord] = sort(off(:), 'descend');
fprintf('\nLargest OLD -> NEW(matched) reassignments:\n');
for q = 1:min(6, sum(vals>0))
    [i,j] = ind2sub(size(off), ord(q));
    fprintf('  %-20s -> %-20s : %d\n', names7{i}, names7{j}, vals(q));
end

%% ================= FIGURES: ENSEMBLE DIAGRAMS =======================================
% OLD: grey/black = MATCONT main branch, dashed = BP-emanating branches, red = LP, green = BP
plot_old_ensemble(resOld, old_cls, names7, maxPerClass, xl_plot, ...
    fullfile(figuresDir,'OLD_Branches_by_Class.png'));

% NEW: exactly the v2 figure (all branches, detached branch in red)
fig2 = figure('Color','w','Position',[50 50 1500 850]);
for c = 1:9
    ax = subplot(3,3,c); hold(ax,'on'); box(ax,'on');
    ks = find(new9F == c); ks = ks(1:min(maxPerClass, numel(ks)));
    for k = ks'
        draw_branches(ax, resNew(k).curves, 0.8);
        if ~isempty(resNew(k).exact_LP_ATOT)
            plot(ax, resNew(k).exact_LP_ATOT, resNew(k).exact_LP_x2, 'ro', 'MarkerFaceColor','r','MarkerSize',4);
        end
        if ~isempty(resNew(k).exact_BP_ATOT)
            plot(ax, resNew(k).exact_BP_ATOT, resNew(k).exact_BP_x2, 'gs', 'MarkerFaceColor','g','MarkerSize',5);
        end
    end
    set(ax,'XScale','log','YScale','log'); xlim(ax, xl_plot);
    title(ax, sprintf('%s (n=%d)', names9{c}, sum(new9F==c)), 'Interpreter','none');
    xlabel(ax,'ATOT'); ylabel(ax,'x_2 (App)');
end
saveas(fig2, fullfile(figuresDir,'NEW_Branches_by_Class.png'));

%% ================= FIGURES: BAR PLOTS ===============================================
f = figure('Color','w','Position',[100 100 1000 450]);
bar_pct(gca, cOld, names7, sprintf('OLD method, 7+1 classes (n=%d)', numSamples), numSamples);
saveas(f, fullfile(figuresDir,'OLD_Class_Distribution.png'));

f = figure('Color','w','Position',[100 100 1100 450]);
bar_pct(gca, c9, names9, sprintf('NEW method, 9 classes, ATOT window [%g %g] (n=%d, failed %d)', ...
    cfg.ATOT_window, numSamples, nFail9), numSamples);
saveas(f, fullfile(figuresDir,'NEW_Class_Distribution.png'));

%% ================= FIGURES: COMPARISON ==============================================
f = figure('Color','w','Position',[100 100 1300 480]);
subplot(1,2,1);
if run_new_matched, bar(1:8, [cOld cNewM cNewF]); legend({'OLD','NEW matched window','NEW full window'}, 'Location','northeast');
else, bar(1:8, [cOld cNewF]); legend({'OLD','NEW full window'}, 'Location','northeast'); end
set(gca,'XTick',1:8,'XTickLabel',names7,'XTickLabelRotation',40,'TickLabelInterpreter','none','FontSize',9);
ylabel('Samples'); title('Class counts, same 7+1 bins'); grid on;
subplot(1,2,2);
dd = [cNewM - cOld, cNewF - cOld]/numSamples*100;
bar(1:8, dd); legend({'NEW matched - OLD','NEW full - OLD'}, 'Location','best');
set(gca,'XTick',1:8,'XTickLabel',names7,'XTickLabelRotation',40,'TickLabelInterpreter','none','FontSize',9);
ylabel('Difference (% of samples)'); title('What the old method over/under-counts'); grid on;
saveas(f, fullfile(figuresDir,'Compare_Class_Counts.png'));

f = figure('Color','w','Position',[50 50 1500 650]);
conf_plot(subplot(1,2,1), CM_M, names7, sprintf('OLD vs NEW matched window (agreement %.1f%%)', 100*agM));
conf_plot(subplot(1,2,2), CM_F, names7, sprintf('OLD vs NEW full window (agreement %.1f%%)', 100*agF));
saveas(f, fullfile(figuresDir,'Compare_Confusion.png'));

% examples: old curve grey underneath, new curves on top
dis = find(old_cls ~= new7M);
if run_new_matched
    plot_examples(dis, resOld, resMatch, old_cls, new7M, names7, wMatch, ...
        fullfile(figuresDir,'Examples_Disagreement_matched.png'));
end
det = find(detF);
plot_examples(det, resOld, resNew, old_cls, new7F, names7, repmat(cfg.ATOT_window(2),numSamples,1), ...
    fullfile(figuresDir,'Examples_Detached_full.png'));

%% ================= SAVE =============================================================
T = array2table(paramSamples, 'VariableNames', param_names);
T.Old_Status        = old_status;
T.Old_Class7        = old_cls;
T.Old_Label         = reshape(names7(old_cls), [], 1);
T.Old_ATOT_reached  = [resOld.ATOT_reached]';
T.Old_Min_State     = [resOld.min_state]';
T.Old_LP_ATOT       = arrayfun(@(r) mat2str(r.LP_ATOT,6), resOld, 'uni', 0);
T.Old_BP_ATOT       = arrayfun(@(r) mat2str(r.BP_ATOT,6), resOld, 'uni', 0);
T.New_Class9_full   = new9F;
T.New_Label9_full   = reshape(names9(max(new9F,1)), [], 1);
T.New_Class7_full   = new7F;
T.New_Class7_matched = new7M;
T.Matched_window_max = wMatch;
T.Multistationary   = msF;
T.Has_Detached      = detF;
T.New_LP_ATOT       = arrayfun(@(r) mat2str(r.exact_LP_ATOT,6), resNew, 'uni', 0);
T.New_BP_ATOT       = arrayfun(@(r) mat2str(r.exact_BP_ATOT,6), resNew, 'uni', 0);
T.Agree_matched     = old_cls == new7M;
T.Agree_full        = old_cls == new7F;

xf = fullfile(excelDir,'Old_vs_New_Comparison.xlsx');
writetable(T, xf, 'Sheet','Per_Sample');
Tc = table((1:8)', names7', cOld, cNewM, cNewF, ...
    'VariableNames', {'Bin','Class','OLD','NEW_matched_window','NEW_full_window'});
writetable(Tc, xf, 'Sheet','Class_Counts');
Ts = table({'Agreement_OLD_vs_NEWmatched';'Agreement_OLD_vs_NEWfull';'OLD_fold_or_BP_sensitivity';'OLD_fold_or_BP_specificity'}, ...
    100*[agM;agF;sens;spec], 100*[agM_lo;agF_lo;sens_lo;spec_lo], 100*[agM_hi;agF_hi;sens_hi;spec_hi], ...
    'VariableNames', {'Quantity','Percent','CI95_low','CI95_high'});
writetable(Ts, xf, 'Sheet','Agreement');
writetable(array2table(CM_M,'VariableNames',matlab.lang.makeValidName(strcat('NEW_',names7)), ...
    'RowNames',strcat('OLD_',names7)), xf, 'Sheet','Confusion_matched','WriteRowNames',true);
writetable(array2table(CM_F,'VariableNames',matlab.lang.makeValidName(strcat('NEW_',names7)), ...
    'RowNames',strcat('OLD_',names7)), xf, 'Sheet','Confusion_full','WriteRowNames',true);
Tt = table({'OLD';'NEW_full';'NEW_matched'}, [tOld;tNew;tMatch], 'VariableNames', {'Method','Seconds'});
writetable(Tt, xf, 'Sheet','Runtime');

save(fullfile(matDir,'Old_vs_New_Results.mat'), 'resOld','resNew','paramSamples','param_names', ...
    'bounds','cfg','rngSeed','wMatch','old_cls','new7F','new9F','new7M','CM_M','CM_F', ...
    'tOld','tNew','tMatch','-v7.3');
if run_new_matched, save(fullfile(matDir,'Old_vs_New_Results.mat'),'resMatch','-append'); end
fprintf('\nDone in %.1f min. Results in %s\n', toc(total_start)/60, outputDir);

%% ================= LOCAL FUNCTIONS (client side only) ===============================
function out = run_parallel(fun, N, nw, label)
out = cell(N,1);
progress_tick('reset');
if nw > 0
    dq = parallel.pool.DataQueue;
    afterEach(dq, @(~) progress_tick(label, N));
    parfor (k = 1:N, nw)
        out{k} = feval(fun, k);
        send(dq, 1);
    end
else
    for k = 1:N
        out{k} = fun(k);
        progress_tick(label, N);
    end
end
end

function progress_tick(label, N)
persistent n
if nargin == 1 && strcmp(label,'reset'), n = 0; return; end
n = n + 1;
if mod(n, max(1, round(N/20))) == 0 || n == N
    fprintf('  [%s] %d / %d\n', label, n, N);
end
end

function [c7, c9, ok] = new_labels(res)
ok = [res.exact_ok]';
c9 = [res.exact_class]';      c9(~ok) = 0;
c7 = [res.exact_class_old7]'; c7(~ok) = 8;
end

function [p, lo, hi] = wilson(k, n)
if n == 0, p = NaN; lo = NaN; hi = NaN; return; end
z = 1.959964; p = k/n;
c = (p + z^2/(2*n)) / (1 + z^2/n);
h = z*sqrt(p*(1-p)/n + z^2/(4*n^2)) / (1 + z^2/n);
lo = max(0, c - h); hi = min(1, c + h);
end

function bar_pct(ax, counts, labels, ttl, N)
bar(ax, counts, 'FaceColor', [0.25 0.5 0.75]); hold(ax,'on');
set(ax,'XTick',1:numel(counts),'XTickLabel',labels,'XTickLabelRotation',40, ...
    'TickLabelInterpreter','none','FontSize',9);
for i = 1:numel(counts)
    text(ax, i, counts(i), sprintf('%d\n(%.1f%%)', counts(i), 100*counts(i)/N), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
end
ylabel(ax,'Samples'); title(ax, ttl, 'Interpreter','none'); grid(ax,'on');
ylim(ax, [0 max(counts)*1.15 + 1]);
end

function conf_plot(ax, CM, names, ttl)
axes(ax); %#ok<LAXES>
imagesc(log10(1 + CM)); colorbar; axis square;
set(ax,'XTick',1:numel(names),'XTickLabel',names,'YTick',1:numel(names),'YTickLabel',names, ...
    'XTickLabelRotation',40,'FontSize',8,'TickLabelInterpreter','none');
xlabel('NEW class'); ylabel('OLD class'); title(ttl,'Interpreter','none');
for i = 1:size(CM,1)
    for j = 1:size(CM,2)
        if CM(i,j) > 0, text(j, i, num2str(CM(i,j)), 'HorizontalAlignment','center','FontSize',8,'Color','w'); end
    end
end
end

function plot_old_curve(ax, r, col, lw)
if ~isempty(r.main), plot(ax, r.main(1,:), r.main(2,:), '-', 'Color', col, 'LineWidth', lw); end
for j = 1:numel(r.ext)
    e = r.ext{j};
    plot(ax, e(1,:), e(2,:), '--', 'Color', col, 'LineWidth', lw);
end
end

function plot_old_ensemble(res, cls, names, maxPer, xl, fname)
f = figure('Color','w','Position',[50 50 1500 750]);
for c = 1:7
    ax = subplot(2,4,c); hold(ax,'on'); box(ax,'on');
    ks = find(cls == c); ks = ks(1:min(maxPer, numel(ks)));
    for k = ks'
        plot_old_curve(ax, res(k), [0.15 0.15 0.15], 0.8);
        if ~isempty(res(k).LP_ATOT)
            plot(ax, res(k).LP_ATOT, res(k).LP_x2, 'ro', 'MarkerFaceColor','r','MarkerSize',4);
        end
        if ~isempty(res(k).BP_ATOT)
            plot(ax, res(k).BP_ATOT, res(k).BP_x2, 'gs', 'MarkerFaceColor','g','MarkerSize',5);
        end
    end
    set(ax,'XScale','log','YScale','log'); xlim(ax, xl);
    title(ax, sprintf('%s (n=%d)', names{c}, sum(cls==c)), 'Interpreter','none');
    xlabel(ax,'ATOT'); ylabel(ax,'x_2 (App)');
end
ax = subplot(2,4,8); axis(ax,'off');
text(ax, 0, 0.6, sprintf('Failed: %d\nNegative concentration: %d\n(not drawn)', ...
    sum(~[res.ok]), sum([res.negative])), 'FontSize', 11);
saveas(f, fname);
end

function plot_examples(idx, resOld, resNew, oldCls, newCls, names, xmax, fname)
if isempty(idx), return; end
idx = idx(1:min(6, numel(idx)));
f = figure('Color','w','Position',[50 50 1500 800]);
for j = 1:numel(idx)
    k = idx(j);
    ax = subplot(2,3,j); hold(ax,'on'); box(ax,'on');
    plot_old_curve(ax, resOld(k), [0.7 0.7 0.7], 3.5);          % OLD = thick grey
    draw_branches(ax, resNew(k).curves, 1.2);                   % NEW on top
    if ~isempty(resNew(k).exact_LP_ATOT)
        plot(ax, resNew(k).exact_LP_ATOT, resNew(k).exact_LP_x2, 'ro','MarkerFaceColor','r','MarkerSize',4);
    end
    set(ax,'XScale','log','YScale','log'); xlim(ax, [1e-3 xmax(k)]);
    title(ax, sprintf('#%d  old: %s | new: %s', k, names{oldCls(k)}, names{newCls(k)}), ...
        'Interpreter','none','FontSize',8);
    xlabel(ax,'ATOT'); ylabel(ax,'x_2 (App)');
end
saveas(f, fname);
end