% =================== ML_v3_misclassified_diagrams.m ==================================
% Which parameter sets does the classifier get wrong, and what do their bifurcation
% diagrams look like?
%
% The result files of ML_v3_10class_parallel.m store the confusion COUNTS, not the
% individual predictions. This script therefore re-creates the pooled LightGBM model
% of the confusion matrix EXACTLY as the ML script does it:
%   * same input files (ML_Ready_Dataset_v3_iter<i>.csv), same row filtering
%     (Use_For_ML, classes < 50 sets per iteration), same 14 features -> feature set A
%     (14 log-parameters + 91 pairwise log-ratios),
%   * same stratified train / validation / test split (key = iteration x class,
%     seed 43 = cfg.seed + 1, 70 / 10 / 20 %), same model settings and seed,
%   * the model is fitted on the TRAINING rows only (support/ml_v3_core.m),
% and predicts the held-out TEST rows (as in Fig "Pooled confusion matrix") and/or the
% VALIDATION rows. The accuracy is printed next to the value in ML_v3_summary.txt so
% you can check that the run is reproduced.
%
% For every confusion  true class -> predicted class  (largest first) it draws
% nPerPair representative misclassified parameter sets: their EXACT bifurcation
% diagram (x2 = APP vs ATOT; solid = stable, dashed = unstable; blue = main branch,
% red = detached branch, purple = crossing branch; red circles = LP, green squares =
% BP), titled with the true and the predicted class, plus (optional) one correctly
% classified set of the predicted class for comparison. A short text in each panel
% gives the singularities (number of main / detached folds, loop width, detached-fold
% ATOT, BP), i.e. why the set belongs to its true class.
%
% Run from the Ensemble_v3 folder. Needs Python with lightgbm (as for the ML script)
% and Ensemble_iter<i>_v3_1e6/MAT_Files/Complete_Bifurcation_Results_v3.mat for the
% diagrams (each file is loaded once; only the chosen sets are kept).
% Output: ML_v3_results_parallel/labels_Exact_Class_v3/misclassified/
%   Misclassified_<set>.csv            every misclassified set (iteration, sample,
%                                      true / predicted class, singularities, params)
%   Confusion_pairs_<set>.csv          count of every confusion
%   Misclassified_overview_<set>.png/.fig   top pairs, one row per pair
%   Pair_<true>_to_<pred>_<set>.png/.fig    one figure per confusion pair
% ===================================================================================
clear; close all; clc;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
addpath(fullfile(HERE, 'support')); cd(HERE);

%% ---------------- USER SETTINGS ----------------------------------------------------
evalSet       = 'test';      % 'test' (= pooled confusion matrix of the paper),
                             % 'validation', or 'both'
model         = 'LightGBM';  % any model of ml_v3_core: LightGBM, XGBoost, RandomForest, ...
nPerPair      = 3;           % misclassified examples drawn per confusion pair
maxPairs      = 12;          % confusion pairs drawn (largest counts first)
minPairCount  = 1;           % ignore confusions with fewer sets
nOverview     = 6;           % pairs (rows) in the overview figure
showReference = true;        % add one correctly classified set of the predicted class
pickSeed      = 11;          % which misclassified sets are drawn (random, seeded)
pickMode      = 'spread';    % 'spread' : evenly spread over the misclassified sets of the pair
                             % 'random' : random
logZoom       = 30;          % log-log window: min(singular ATOT)/logZoom .. max*logZoom
LW = 2.5;                    % line width of the branches
FS = 11;                     % font size
reportedAcc   = 92.61;       % pooled LightGBM test accuracy in ML_v3_summary.txt (check)
resultsRoot   = fullfile(HERE, 'ML_v3_results_parallel', 'labels_Exact_Class_v3');

% ---- settings of ML_v3_10class_parallel.m (do not change: needed to reproduce it) --
cfg = struct();
cfg.iterIds   = [1 2 3 4];
cfg.featNames = {'KTOT','PTOT','k1','k11','k2','k3','k33','k4','k5','k55','k6','k7','k77','k8'};
cfg.labelCol  = 'Exact_Class_v3';  cfg.nameCol = 'Exact_Label_v3';
cfg.useForMLCol = 'Use_For_ML';    cfg.minClassPerIter = 50;
cfg.seed = 42; cfg.testFrac = 0.20; cfg.valFrac = 0.10; cfg.nBins = 30; cfg.tolEq = 1e-9;
cfg.useClassWeights = true; cfg.nThreads = -1;
cfg.lgbRounds = 300; cfg.lgbLeaves = 63; cfg.rfTrees = 300; cfg.xgbRounds = 300; cfg.xgbDepth = 6;
cfg.catIter = 300; cfg.catDepth = 6; cfg.gbRounds = 200; cfg.gbDepth = 5; cfg.gbFast = false;

%% ---------------- 1. data, exactly as the ML script ---------------------------------
CI = class_info_v3();
X = []; It = []; Smp = []; lab = strings(0,1); use = false(0,1);
for it = cfg.iterIds
    f = fullfile(HERE, sprintf('Ensemble_iter%d_v3_1e6', it), 'Excel_Files', ...
        sprintf('ML_Ready_Dataset_v3_iter%d.csv', it));
    if ~isfile(f), error('Missing %s', f); end
    T = readtable(f);
    Xk = zeros(height(T), numel(cfg.featNames));
    for j = 1:numel(cfg.featNames)
        c = T.(cfg.featNames{j}); if ~isnumeric(c), c = str2double(strtrim(string(c))); end
        Xk(:,j) = c;
    end
    u = T.(cfg.useForMLCol);
    if ~islogical(u) && ~isnumeric(u), u = strcmpi(strtrim(string(u)),'true') | strtrim(string(u)) == "1"; end
    X = [X; Xk]; It = [It; repmat(it, height(T), 1)]; Smp = [Smp; (1:height(T))']; %#ok<AGROW>
    lab = [lab; strtrim(string(T.(cfg.labelCol)))]; use = [use; logical(u(:))]; %#ok<AGROW>
    fprintf('loaded iteration %d: %d rows\n', it, height(T));
end
bad = any(isnan(X), 2);
X = X(~bad,:); It = It(~bad); Smp = Smp(~bad); lab = lab(~bad); use = use(~bad);
keep = ~ismissing(lab) & lab ~= "" & lab ~= "0" & use;
for a = unique(It)'
    la = lab(keep & It == a); [uq, ~, g] = unique(la); cnt = accumarray(g, 1);
    for q = find(cnt < cfg.minClassPerIter)', keep = keep & ~(It == a & lab == uq(q)); end
end
X = X(keep,:); It = It(keep); Smp = Smp(keep); lab = lab(keep);
lblNum = str2double(lab);
catStr = string(sort(unique(lblNum)))';
[~, y] = ismember(string(lblNum), catStr);
clsNum = double(catStr);                         % class number of every code
fprintf('ML data: %d rows, %d classes\n', numel(y), numel(catStr));

% feature set A: 14 log-parameters + 91 pairwise log-ratios
L = log10(X); p = size(L,2);
pairs = nchoosek(1:p, 2);
Rlog = L(:,pairs(:,1)) - L(:,pairs(:,2));
names = [cfg.featNames, arrayfun(@(q) sprintf('%s/%s', cfg.featNames{pairs(q,1)}, ...
    cfg.featNames{pairs(q,2)}), 1:size(pairs,1), 'UniformOutput', false)];
fsA = struct('name','A: 14 params + 91 log-ratios', 'X', [L, Rlog], 'names', {names}, 'jsd', false);

% the pooled split of the ML script (first repeat) and the model seed of its confusion
sp = stratSplit(It*1000 + y, cfg.testFrac, cfg.seed + 1, cfg.valFrac);
switch lower(evalSet)
    case 'test',       ev = sp.te; setName = 'test';
    case 'validation', ev = sp.va; setName = 'validation';
    case 'both',       ev = [sp.te; sp.va]; setName = 'test+validation';
    otherwise, error('evalSet must be test, validation or both');
end
fprintf('train %d | validation %d | test %d  -> predicting the %s rows with %s\n', ...
    numel(sp.tr), numel(sp.va), numel(sp.te), setName, model);

%% ---------------- 2. fit on the training rows, predict -------------------------------
res  = ml_v3_core('evalFS', fsA, y, sp.tr, ev, cfg, model, cfg.seed + 1);
pred = res.pred(:); yt = y(ev);
acc  = 100*mean(pred == yt);
fprintf('%s accuracy on the %s rows: %.2f%%', model, setName, acc);
if strcmpi(evalSet, 'test') && strcmp(model, 'LightGBM')
    fprintf('   (ML_v3_summary.txt: %.2f%%)', reportedAcc);
end
fprintf('\n');
if strcmpi(evalSet, 'both')
    nT = numel(sp.te);
    fprintf('   test %.2f%% | validation %.2f%%\n', 100*mean(pred(1:nT) == yt(1:nT)), ...
        100*mean(pred(nT+1:end) == yt(nT+1:end)));
end

%% ---------------- 3. tables of the misclassified sets ---------------------------------
outDir = fullfile(resultsRoot, 'misclassified');
if ~exist(outDir, 'dir'), mkdir(outDir); end
wrong = find(pred ~= yt);
rows  = ev(wrong);
tTrue = clsNum(yt(wrong))'; tPred = clsNum(pred(wrong))';
Tm = table(It(rows), Smp(rows), tTrue(:), string(CI.names(tTrue(:)))', tPred(:), ...
    string(CI.names(tPred(:)))', 'VariableNames', ...
    {'Iteration','Sample','TrueClass','TrueLabel','PredictedClass','PredictedLabel'});
Tm = [Tm, array2table(X(rows,:), 'VariableNames', cfg.featNames)];
[pp, ~, gp] = unique([tTrue(:) tPred(:)], 'rows');
cnt = accumarray(gp, 1);
nTrueClass = arrayfun(@(c) sum(clsNum(yt) == c), pp(:,1));
Tp = table(pp(:,1), string(CI.names(pp(:,1)))', pp(:,2), string(CI.names(pp(:,2)))', cnt, ...
    100*cnt./nTrueClass, 'VariableNames', {'TrueClass','TrueLabel','PredictedClass', ...
    'PredictedLabel','Count','PctOfTrueClass'});
Tp = sortrows(Tp, 'Count', 'descend');
Tp = Tp(Tp.Count >= minPairCount, :);
tag = strrep(setName, '+', '_');
writetable(Tp, fullfile(outDir, sprintf('Confusion_pairs_%s.csv', tag)));
fprintf('\n%d misclassified sets; confusion pairs (largest first):\n', numel(wrong));
disp(Tp(1:min(height(Tp), maxPairs), :));

%% ---------------- 4. choose the sets to draw ----------------------------------------
nPairs = min(maxPairs, height(Tp));
rsP = RandStream('mt19937ar', 'Seed', pickSeed);
pick = cell(nPairs, 1); ref = zeros(nPairs, 1);
for q = 1:nPairs
    in = find(Tm.TrueClass == Tp.TrueClass(q) & Tm.PredictedClass == Tp.PredictedClass(q));
    k = min(nPerPair, numel(in));
    if strcmpi(pickMode, 'random')
        sel = in(randperm(rsP, numel(in), k));
    else
        sel = in(round(linspace(1, numel(in), k)));          % ordered later by ATOT
    end
    pick{q} = sel;
    if showReference                                          % correct set of the predicted class
        ok = ev(pred == yt & clsNum(yt)' == Tp.PredictedClass(q));
        if ~isempty(ok), ref(q) = ok(randperm(rsP, numel(ok), 1)); end
    end
end
need = unique([cell2mat(cellfun(@(s) rows(s), pick, 'UniformOutput', false)); ref(ref > 0)]);
% load the exact results of the chosen sets
Rsel = containers.Map('KeyType','double','ValueType','any');
for it = cfg.iterIds
    nn = need(It(need) == it);
    if isempty(nn), continue; end
    f = fullfile(HERE, sprintf('Ensemble_iter%d_v3_1e6', it), 'MAT_Files', 'Complete_Bifurcation_Results_v3.mat');
    fprintf('Loading %s (%d sets) ...\n', f, numel(nn));
    S = load(f, 'res');
    for r = nn(:)', Rsel(r) = S.res(Smp(r)); end
    clear S
end

% add singularity summaries to the table of all misclassified sets that were loaded
Tm.Description = repmat("", height(Tm), 1);
for i = 1:height(Tm)
    r = rows(i);
    if isKey(Rsel, r), Tm.Description(i) = string(strjoin(describe(Rsel(r)), '; ')); end
end
writetable(Tm, fullfile(outDir, sprintf('Misclassified_%s.csv', tag)));

%% ---------------- 5. one figure per confusion pair ---------------------------------
nCol = nPerPair + showReference;
for q = 1:nPairs
    fig = figure('Color','w', 'Units','pixels', 'Position',[40 80 430*nCol 470], 'Visible','on');
    tl = tiledlayout(fig, 1, nCol, 'TileSpacing','compact', 'Padding','compact');
    title(tl, sprintf('True class %d (%s) predicted as class %d (%s): %d of %d sets (%.1f%%), %s set, %s', ...
        Tp.TrueClass(q), CI.labels{Tp.TrueClass(q)}, Tp.PredictedClass(q), CI.labels{Tp.PredictedClass(q)}, ...
        Tp.Count(q), nTrueClass(pp(:,1) == Tp.TrueClass(q) & pp(:,2) == Tp.PredictedClass(q)), ...
        Tp.PctOfTrueClass(q), setName, model), 'FontSize', FS+2, 'FontWeight','bold', 'Interpreter','tex');
    for i = 1:numel(pick{q})
        r = rows(pick{q}(i));
        drawDiagram(nexttile(tl), Rsel(r), It(r), Smp(r), Tp.TrueClass(q), Tp.PredictedClass(q), CI, false, logZoom, LW, FS);
    end
    for i = numel(pick{q})+1:nPerPair, nexttile(tl); axis off; end
    if showReference && ref(q) > 0
        r = ref(q);
        drawDiagram(nexttile(tl), Rsel(r), It(r), Smp(r), Tp.PredictedClass(q), Tp.PredictedClass(q), CI, true, logZoom, LW, FS);
    end
    base = fullfile(outDir, sprintf('Pair_%d_to_%d_%s', Tp.TrueClass(q), Tp.PredictedClass(q), tag));
    saveBoth(fig, base);
    if q > nOverview, close(fig); end
end

%% ---------------- 6. overview: top pairs, one row each ------------------------------
nO = min(nOverview, nPairs);
fig = figure('Color','w', 'Units','pixels', 'Position',[20 20 380*nCol 360*nO]);
tl = tiledlayout(fig, nO, nCol, 'TileSpacing','compact', 'Padding','compact');
title(tl, sprintf('Misclassified parameter sets (%s, %s set, accuracy %.2f%%): true class -> predicted class', ...
    model, setName, acc), 'FontSize', FS+3, 'FontWeight','bold');
for q = 1:nO
    for i = 1:nPerPair
        ax = nexttile(tl);
        if i <= numel(pick{q})
            r = rows(pick{q}(i));
            drawDiagram(ax, Rsel(r), It(r), Smp(r), Tp.TrueClass(q), Tp.PredictedClass(q), CI, false, logZoom, LW, FS-1);
        else
            axis(ax, 'off');
        end
    end
    if showReference
        ax = nexttile(tl);
        if ref(q) > 0
            r = ref(q);
            drawDiagram(ax, Rsel(r), It(r), Smp(r), Tp.PredictedClass(q), Tp.PredictedClass(q), CI, true, logZoom, LW, FS-1);
        else
            axis(ax, 'off');
        end
    end
end
saveBoth(fig, fullfile(outDir, sprintf('Misclassified_overview_%s', tag)));
fprintf('\nAll outputs in %s\n', outDir);

%% ============== LOCAL FUNCTIONS =====================================================
function drawDiagram(ax, R, it, smp, cTrue, cPred, CI, isRef, logZoom, LW, FS)
% exact bifurcation diagram of one set; log-log when a detached branch or a branch
% point is present (branches on very different scales), linear zoom otherwise
hold(ax, 'on'); box(ax, 'on');
hasLog = R.n_LP_detached > 0 || ~isempty(R.exact_BP_ATOT) || any(cTrue == [2 4 5 6 7 8]) || any(cPred == [2 4 5 6 7 8]);
if hasLog && ~isempty([R.exact_LP_ATOT(:); R.exact_BP_ATOT(:)])
    [xl, yl] = diagram_window(R, 'log', logZoom); isLog = true;
else
    [xl, yl] = diagram_window(R, 'lin', logZoom); isLog = false;
end
draw_branches(ax, R.curves, LW);
if ~isempty(R.exact_LP_ATOT)
    plot(ax, R.exact_LP_ATOT, R.exact_LP_x2, 'o', 'MarkerFaceColor', [0.9 0.05 0.05], ...
        'MarkerEdgeColor','k', 'MarkerSize', 8, 'LineWidth', 1);
end
if ~isempty(R.exact_BP_ATOT)
    plot(ax, R.exact_BP_ATOT, R.exact_BP_x2, 's', 'MarkerFaceColor', [0 0.62 0.2], ...
        'MarkerEdgeColor','k', 'MarkerSize', 9, 'LineWidth', 1);
end
if isLog, set(ax, 'XScale','log', 'YScale','log'); end
if all(isfinite(xl)) && xl(2) > xl(1), xlim(ax, xl); end
if all(isfinite(yl)) && yl(2) > yl(1), ylim(ax, yl); end
set(ax, 'FontSize', FS, 'FontWeight','bold', 'LineWidth', 1.5);
xlabel(ax, 'A_{TOT}'); ylabel(ax, 'Steady state (APP)');
if isRef
    title(ax, {sprintf('REFERENCE: class %d, correctly predicted', cTrue), ...
        sprintf('%s | iter %d, sample %d', CI.labels{cTrue}, it, smp)}, ...
        'FontSize', FS, 'Color', [0.35 0.35 0.35], 'Interpreter','tex');
else
    title(ax, {sprintf('True %d: %s', cTrue, CI.labels{cTrue}), ...
        sprintf('\\color[rgb]{0.8,0,0}Predicted %d: %s', cPred, CI.labels{cPred}), ...
        sprintf('\\color{black}iter %d, sample %d', it, smp)}, 'FontSize', FS, 'Interpreter','tex');
end
d = describe(R);
text(ax, 0.03, 0.97, d, 'Units','normalized', 'VerticalAlignment','top', 'FontSize', FS-2, ...
    'FontWeight','normal', 'BackgroundColor', [1 1 1], 'Interpreter','none', 'Margin', 1);
end

function d = describe(R)
% singularities of one set in words (why it belongs to its exact class)
d = {};
nM = R.n_LP_main; nD = R.n_LP_detached; nB = numel(R.exact_BP_ATOT);
types = R.exact_LP_type; if ischar(types), types = {types}; end
Am = R.exact_LP_ATOT(strcmp(types, 'main'));
if isempty(Am) && nM > 0, Am = R.exact_LP_ATOT(~strcmp(types, 'detached')); end
if nM >= 2
    w = (max(Am) - min(Am))/mean(Am);
    d{end+1} = sprintf('main-branch folds: %d at A = %.3g-%.3g (loop width %.2g%%)', nM, min(Am), max(Am), 100*w);
elseif nM == 1
    d{end+1} = sprintf('main-branch fold: 1 at A = %.3g', Am);
else
    d{end+1} = 'main branch: no fold';
end
if nD > 0
    d{end+1} = sprintf('detached fold at A = %.3g', min(R.detached_fold_ATOT));
elseif R.has_detached
    d{end+1} = 'detached branch: fold above window';
end
Ao = R.exact_LP_ATOT(strcmp(types, 'offsym'));
if ~isempty(Ao), d{end+1} = sprintf('crossing-branch fold at A = %.3g', Ao(1)); end
if nB > 0, d{end+1} = sprintf('branch point at A = %.3g', R.exact_BP_ATOT(1)); end
d{end+1} = sprintf('max. %d steady states', R.max_ss);
end

function sp = stratSplit(key, testFrac, seed, valFrac)
% identical to ML_v3_10class_parallel.m
    if nargin < 4; valFrac = 0; end
    rng(seed);
    n = numel(key); part = zeros(n,1);
    u = unique(key);
    for i = 1:numel(u)
        idx = find(key == u(i)); m = numel(idx);
        if m < 2; continue; end
        perm = idx(randperm(m));
        nTe = min(max(1, round(testFrac*m)), m-1);
        nVa = min(round(valFrac*m), m - 1 - nTe);
        part(perm(1:nTe)) = 2;
        part(perm(nTe+1:nTe+nVa)) = 1;
    end
    sp.tr = find(part == 0); sp.va = find(part == 1); sp.te = find(part == 2);
end

function saveBoth(fig, base)
savefig(fig, [base '.fig']);
try
    exportgraphics(fig, [base '.png'], 'Resolution', 300);
catch
    print(fig, base, '-dpng', '-r300');
end
fprintf('saved %s (.png, .fig)\n', base);
end
