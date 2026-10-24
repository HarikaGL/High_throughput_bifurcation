%% ============================================================
%  PER-CLASS FEATURE CONTRIBUTIONS (TreeSHAP) WITHIN EACH SAMPLING ITERATION
%                                               ML_v3_SHAP_per_class.m
%  Replaces the split-gain importance of the manuscript's within-iteration table
%  ("Within-iteration feature importance recovers the analytical conditions") by
%  per-CLASS contributions. Split gain is summed over all trees of a multiclass model and
%  cannot say WHICH class a feature informs; TreeSHAP (exact for tree models, built into
%  LightGBM as pred_contrib) gives the contribution of every feature to every class's
%  raw score, for every parameter set, with a sign.
%
%  Same data, features, model and splits as section ML-d of ML_v3_10class_parallel.m:
%    * rows: Use_For_ML = 1; classes with fewer than cfg.minClassPerIter rows in an
%      iteration are removed (as in the main pipeline)
%    * features: set A = 14 log10 parameters + 91 pairwise log10 ratios; constant and
%      duplicated columns are dropped using the TRAINING rows only
%    * model: LightGBM, 300 rounds, 63 leaves, learning rate 0.1, min_child_samples 3,
%      inverse-frequency sample weights
%    * splits: stratified 70 / 10 / 20 % (train / validation / test), seeds
%      cfg.seed + 100*a + r (a = iteration index, r = repeat), as in ML-d
%  SHAP values are computed on the held-out TEST rows of each iteration and averaged over
%  cfg.nRepeats repeats. Analyses are run within each iteration (never pooled), because
%  pooled features carry the sampling design (confound audit of the main pipeline).
%
%  OUTPUT (folder ML_v3_results_SHAP/):
%   SHAP_iter<i>_meanabs.csv        class x feature: mean |SHAP| over all test rows
%   SHAP_iter<i>_meanabs_own.csv    class x feature: mean |SHAP| over test rows OF that class
%   SHAP_iter<i>_direction.csv      class x feature: correlation between the feature value
%                                   and its SHAP value for that class (+: a larger value
%                                   raises the class score; -: lowers it)
%   SHAP_iter<i>_parameters.csv     class x 14 parameters: mean |attribution| when the SHAP
%                                   value of each log-ratio log(a/b) is split equally between
%                                   a and b (SHAP values are additive, so sums are legitimate)
%   SHAP_top_features_per_class.csv iteration, class, rank, feature, mean |SHAP|, share of the
%                                   class total (%), direction, test accuracy
%   SHAP_table_rows.tex             LaTeX rows for the manuscript table
%   figures/SHAP_per_class_heatmap.png/.pdf    four panels (one per iteration): classes x
%                                   union of their top features; colour = mean |SHAP|
%                                   (normalised per class), arrow = direction
%   figures/SHAP_per_class_parameters.png/.pdf  per-class attribution to the 14 parameters
%   figures/SHAP_dependence_top.png/.pdf        SHAP value against feature value for the
%                                   top feature of each class (one panel per class)
%   SHAP_summary.txt                every number printed by the run
%
%  HOW TO RUN: put this file next to ML_v3_10class_parallel.m (folder A, which contains
%  Ensemble_iter<i>_v3_1e6/ and support/) and press Run. Needs Python with numpy and
%  lightgbm through pyenv (the same set-up as ML_v3_10class_parallel.m). No parallel
%  toolbox needed; run time ~ 2-5 min.
% =============================================================
clear; clc; close all;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE); HERE = pwd; end
cd(HERE);
if exist(fullfile(HERE,'support'),'dir'); addpath(fullfile(HERE,'support')); end

cfg = struct();
cfg.iterIds    = [1 2 3 4];
cfg.iterFiles  = arrayfun(@(k) fullfile(sprintf('Ensemble_iter%d_v3_1e6',k), 'Excel_Files', ...
                          sprintf('ML_Ready_Dataset_v3_iter%d.csv',k)), cfg.iterIds, 'UniformOutput', false);
cfg.labelCol   = 'Exact_Class_v3';
cfg.useForMLCol= 'Use_For_ML';
cfg.minClassPerIter = 50;
cfg.featNames  = {'KTOT','PTOT','k1','k11','k2','k3','k33','k4','k5','k55','k6','k7','k77','k8'};
cfg.nRepeats   = 3;            % as in ML-d of the main pipeline
cfg.testFrac   = 0.20;  cfg.valFrac = 0.10;  cfg.seed = 42;
cfg.lgbRounds  = 300;   cfg.lgbLeaves = 63;  cfg.nThreads = -1;
cfg.topK       = 5;            % features listed per class in the table
cfg.heatK      = 4;            % features per class shown in the heatmap (union per iteration)
cfg.FS = 12; cfg.LW = 1.5;

outDir = fullfile(HERE, 'ML_v3_results_SHAP'); figDir = fullfile(outDir, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end
logFile = fopen(fullfile(outDir, 'SHAP_summary.txt'), 'w');
say = @(varargin) sayBoth(logFile, varargin{:});

%% ---- Python ---------------------------------------------------------------------
np  = py.importlib.import_module('numpy');
lgb = py.importlib.import_module('lightgbm');
say('LightGBM %s (Python %s)\n', char(py.getattr(lgb, '__version__')), char(pyenv().Version));

%% ---- feature names of set A ----------------------------------------------------------
fn = cfg.featNames; p = numel(fn);
pairs = nchoosek(1:p, 2);
names = [fn, arrayfun(@(q) sprintf('%s/%s', fn{pairs(q,1)}, fn{pairs(q,2)}), 1:size(pairs,1), 'UniformOutput', false)];
% parameter membership of every feature (for the per-parameter attribution)
memb = zeros(numel(names), p);
memb(sub2ind(size(memb), 1:p, 1:p)) = 1;
for q = 1:size(pairs,1)
    memb(p+q, pairs(q,1)) = 0.5; memb(p+q, pairs(q,2)) = 0.5;
end

TOP = {};  R = struct();
for a = 1:numel(cfg.iterIds)
    it = cfg.iterIds(a);
    %% ---- data of this iteration ------------------------------------------------------
    T = readtable(fullfile(HERE, cfg.iterFiles{a}));
    use = logical(double(T.(cfg.useForMLCol)));
    lab = double(T.(cfg.labelCol));
    keep = use & isfinite(lab);
    [u, ~, g] = unique(lab(keep)); cnt = accumarray(g, 1);
    small = u(cnt < cfg.minClassPerIter);
    keep = keep & ~ismember(lab, small);
    X14 = zeros(height(T), p);
    for j = 1:p, X14(:, j) = double(T.(fn{j})); end
    X14 = X14(keep, :); y = lab(keep);
    L = log10(X14);
    X = [L, L(:, pairs(:,1)) - L(:, pairs(:,2))];
    cls = unique(y); K = numel(cls); [~, yc] = ismember(y, cls);
    say('\n=== Iteration %d: %d rows, classes %s', it, numel(y), mat2str(cls'));
    if ~isempty(small), say(' (removed small classes %s)', mat2str(small')); end
    say('\n');

    nF = numel(names);
    A_all = zeros(K, nF); A_own = zeros(K, nF); Dir = zeros(K, nF); Par = zeros(K, p);
    nUsed = zeros(1, nF); acc = zeros(1, cfg.nRepeats);
    for r = 1:cfg.nRepeats
        sp = stratSplit(yc, cfg.testFrac, cfg.seed + 100*a + r, cfg.valFrac);
        Xtr = X(sp.tr, :); Xte = X(sp.te, :); ytr = yc(sp.tr); yte = yc(sp.te);
        % columns: non-constant and de-duplicated on the training rows (as ml_v3_core)
        cols = find(std(Xtr, 0, 1) > 0);
        [~, ia] = unique(round(Xtr(:, cols), 9)', 'rows', 'stable'); cols = cols(sort(ia));
        F = numel(cols);
        % inverse-frequency weights
        cn = accumarray(ytr, 1, [K 1]); sw = numel(ytr) ./ (K .* cn); sw = sw(ytr);
        % same conversion as ml_v3_core (fitPredictPy)
        Xtr_np = np.array(double(Xtr(:, cols))); Xtr_np = Xtr_np.reshape(int32(numel(ytr)), int32(F));
        Xte_np = np.array(double(Xte(:, cols))); Xte_np = Xte_np.reshape(int32(numel(yte)), int32(F));
        y_np = np.array(double(ytr(:) - 1)); y_np = y_np.ravel();
        sw_np = np.array(double(sw(:)));     sw_np = sw_np.ravel();
        clf = lgb.LGBMClassifier(pyargs('n_estimators', int32(cfg.lgbRounds), ...
            'num_leaves', int32(cfg.lgbLeaves), 'learning_rate', 0.1, 'min_child_samples', int32(3), ...
            'importance_type', 'gain', 'random_state', int32(cfg.seed + r), 'verbose', int32(-1), ...
            'n_jobs', int32(cfg.nThreads)));
        clf.fit(Xtr_np, y_np, pyargs('sample_weight', sw_np));
        pred = double(py.array.array('d', np.ravel(clf.predict(Xte_np)))) + 1;
        acc(r) = 100*mean(pred(:) == yte(:));
        % TreeSHAP: n x K(F+1) (multiclass, class blocks, last column of each block = bias)
        %           n x (F+1)  (binary: contribution to the log-odds of the second class)
        C = clf.booster_.predict(Xte_np, pyargs('pred_contrib', true));
        c = double(py.array.array('d', np.ravel(C)));
        n = numel(yte);
        if K > 2
            c = reshape(c, K*(F+1), n)';                      % numpy is row-major
            S = zeros(n, K, F);
            for k = 1:K, S(:, k, :) = reshape(c(:, (k-1)*(F+1) + (1:F)), n, 1, F); end
        else
            c = reshape(c, F+1, n)';
            S = zeros(n, 2, F); S(:, 2, :) = reshape(c(:, 1:F), n, 1, F); S(:, 1, :) = -S(:, 2, :);
        end
        for k = 1:K
            Sk = squeeze(S(:, k, :)); if F == 1, Sk = Sk(:); end
            A_all(k, cols) = A_all(k, cols) + mean(abs(Sk), 1);
            A_own(k, cols) = A_own(k, cols) + mean(abs(Sk(yte == k, :)), 1);
            cc = zeros(1, F);
            for j = 1:F
                v = corrcoef(Xte(:, cols(j)), Sk(:, j));
                if numel(v) > 1 && isfinite(v(1,2)), cc(j) = v(1,2); end
            end
            Dir(k, cols) = Dir(k, cols) + cc;
            Par(k, :) = Par(k, :) + mean(abs(Sk * memb(cols, :)), 1);
        end
        nUsed(cols) = nUsed(cols) + 1;
    end
    A_all = A_all / cfg.nRepeats; A_own = A_own / cfg.nRepeats;
    Dir = Dir ./ max(nUsed, 1); Par = Par / cfg.nRepeats;
    act = nUsed > 0;
    say('LightGBM test accuracy %.2f +/- %.2f %% (%d repeats), %d features used\n', ...
        mean(acc), std(acc), cfg.nRepeats, sum(act));

    %% ---- tables -------------------------------------------------------------------
    rowN = arrayfun(@(c) sprintf('Class_%d', c), cls, 'UniformOutput', false);
    vn = matlab.lang.makeValidName(names(act));
    writetable(array2table(A_all(:, act), 'VariableNames', vn, 'RowNames', rowN), ...
        fullfile(outDir, sprintf('SHAP_iter%d_meanabs.csv', it)), 'WriteRowNames', true);
    writetable(array2table(A_own(:, act), 'VariableNames', vn, 'RowNames', rowN), ...
        fullfile(outDir, sprintf('SHAP_iter%d_meanabs_own.csv', it)), 'WriteRowNames', true);
    writetable(array2table(Dir(:, act), 'VariableNames', vn, 'RowNames', rowN), ...
        fullfile(outDir, sprintf('SHAP_iter%d_direction.csv', it)), 'WriteRowNames', true);
    writetable(array2table(Par, 'VariableNames', fn, 'RowNames', rowN), ...
        fullfile(outDir, sprintf('SHAP_iter%d_parameters.csv', it)), 'WriteRowNames', true);
    for k = 1:K
        [v, o] = sort(A_all(k, :), 'descend'); tot = sum(A_all(k, :));
        line = sprintf('  Class %2d:', cls(k));
        for q = 1:cfg.topK
            j = o(q); d = '+'; if Dir(k, j) < 0, d = '-'; end
            TOP(end+1, :) = {it, cls(k), q, names{j}, v(q), 100*v(q)/tot, d, mean(acc)}; %#ok<SAGROW>
            line = [line sprintf(' %s (%.3f, %s)', names{j}, v(q), d)]; %#ok<AGROW>
        end
        say('%s\n', line);
    end
    R(a).it = it; R(a).cls = cls; R(a).A = A_all; R(a).Dir = Dir; R(a).Par = Par;
    R(a).names = names; R(a).acc = mean(acc); R(a).X = X; R(a).y = yc;
    % keep the SHAP values of the last repeat for the dependence plots
    R(a).S = S; R(a).Xte = Xte; R(a).cols = cols;
end
TT = cell2table(TOP, 'VariableNames', {'Iteration','Class','Rank','Feature','MeanAbsSHAP', ...
                                       'SharePct','Direction','TestAccuracy'});
writetable(TT, fullfile(outDir, 'SHAP_top_features_per_class.csv'));

%% ---- LaTeX rows for the manuscript table ----------------------------------------------
fid = fopen(fullfile(outDir, 'SHAP_table_rows.tex'), 'w');
for a = 1:numel(R)
    for k = 1:numel(R(a).cls)
        sel = TT(TT.Iteration == R(a).it & TT.Class == R(a).cls(k) & TT.Rank <= 3, :);
        f = cellfun(@(s, d) sprintf('%s (%s)', texName(s), arrow(d)), sel.Feature, sel.Direction, 'UniformOutput', false);
        fprintf(fid, '%d & %d & %s \\\\\n', R(a).it, R(a).cls(k), strjoin(f, ', '));
    end
end
fclose(fid);

%% ---- figure 1: class x feature heatmaps --------------------------------------------
fig = figure('Color','w','Units','centimeters','Position',[1 1 42 12]);
tl = tiledlayout(fig, 1, numel(R), 'TileSpacing','compact','Padding','compact');
cmap = [linspace(1,0.12,256)', linspace(1,0.35,256)', linspace(1,0.65,256)'];
for a = 1:numel(R)
    ax = nexttile(tl); A = R(a).A; K = numel(R(a).cls);
    sel = [];
    for k = 1:K
        [~, o] = sort(A(k, :), 'descend'); sel = [sel, o(1:min(cfg.heatK, end))]; %#ok<AGROW>
    end
    sel = unique(sel, 'stable');
    M = A(:, sel) ./ max(A(:, sel), [], 2);
    imagesc(ax, M'); colormap(ax, cmap); caxis(ax, [0 1]);
    for k = 1:K
        for j = 1:numel(sel)
            if M(k, j) >= 0.15
                if R(a).Dir(k, sel(j)) >= 0, s = char(8593); else, s = char(8595); end
                if M(k, j) > 0.6, tc = 'w'; else, tc = 'k'; end
                text(ax, k, j, s, 'HorizontalAlignment','center', 'FontSize', cfg.FS, 'FontWeight','bold', 'Color', tc);
            end
        end
    end
    set(ax, 'XTick', 1:K, 'XTickLabel', arrayfun(@(c) sprintf('%d', c), R(a).cls, 'UniformOutput', false), ...
        'YTick', 1:numel(sel), 'YTickLabel', colLabel(names(sel)), ...
        'TickLabelInterpreter', 'tex', 'FontSize', cfg.FS-2, 'TickLength', [0 0], 'LineWidth', cfg.LW);
    xlabel(ax, 'Class'); title(ax, sprintf('Iteration %d (test acc. %.1f%%)', R(a).it, R(a).acc), 'FontSize', cfg.FS);
end
cb = colorbar(ax); cb.Layout.Tile = 'east'; cb.Label.String = 'mean |SHAP| (relative to the class maximum)';
title(tl, 'Per-class feature contributions (TreeSHAP, LightGBM, held-out test rows); arrow: a larger value raises (\uparrow) or lowers (\downarrow) the class score', 'FontSize', cfg.FS);
saveBoth(fig, figDir, 'SHAP_per_class_heatmap');

%% ---- figure 2: attribution to the 14 parameters -------------------------------------
fig = figure('Color','w','Units','centimeters','Position',[1 1 42 10]);
tl = tiledlayout(fig, 1, numel(R), 'TileSpacing','compact','Padding','compact');
for a = 1:numel(R)
    ax = nexttile(tl); P = R(a).Par ./ max(R(a).Par, [], 2);
    imagesc(ax, P); colormap(ax, cmap); caxis(ax, [0 1]);
    set(ax, 'YTick', 1:numel(R(a).cls), 'YTickLabel', arrayfun(@(c) sprintf('Class %d', c), R(a).cls, 'UniformOutput', false), ...
        'XTick', 1:p, 'XTickLabel', colLabel(fn), 'TickLabelInterpreter','tex', ...
        'FontSize', cfg.FS-2, 'TickLength', [0 0]);
    xtickangle(ax, 90); title(ax, sprintf('Iteration %d', R(a).it), 'FontSize', cfg.FS);
end
cb = colorbar(ax); cb.Layout.Tile = 'east'; cb.Label.String = 'mean |attribution| (relative)';
title(tl, 'Per-class attribution to the 14 parameters (SHAP of each log-ratio split equally between its two parameters)', 'FontSize', cfg.FS);
saveBoth(fig, figDir, 'SHAP_per_class_parameters');

%% ---- figure 3: dependence of the top feature of every class ---------------------------
nP = sum(arrayfun(@(r) numel(r.cls), R));
fig = figure('Color','w','Units','centimeters','Position',[1 1 44 18]);
tl = tiledlayout(fig, 2, ceil(nP/2), 'TileSpacing','compact','Padding','compact');
for a = 1:numel(R)
    for k = 1:numel(R(a).cls)
        ax = nexttile(tl); hold(ax, 'on'); box(ax, 'on');
        [~, j] = max(R(a).A(k, :)); jj = find(R(a).cols == j, 1);
        if isempty(jj), axis(ax, 'off'); continue; end
        sk = squeeze(R(a).S(:, k, jj));
        scatter(ax, R(a).Xte(:, j), sk, 8, feature_colors(names(j)), 'filled', 'MarkerFaceAlpha', 0.5);
        yline(ax, 0, 'k:');
        xlabel(ax, ['log_{10}(' texLabel(names{j}) ')']); ylabel(ax, 'SHAP');
        title(ax, sprintf('Iter %d, Class %d', R(a).it, R(a).cls(k)), 'FontSize', cfg.FS-1);
        set(ax, 'FontSize', cfg.FS-3);
    end
end
title(tl, 'Top feature of every class: SHAP value against the feature value (held-out test rows)', 'FontSize', cfg.FS);
saveBoth(fig, figDir, 'SHAP_dependence_top');
fclose(logFile);
fprintf('\nResults in %s\n', outDir);

%% =====================================================================================
function sayBoth(fid, varargin)
% print to the command window and to SHAP_summary.txt
    fprintf(1, varargin{:}); fprintf(fid, varargin{:});
end

function sp = stratSplit(key, testFrac, seed, valFrac)
% identical to stratSplit of ML_v3_10class_parallel.m
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

function s = texLabel(name)
% 'k11/k55' -> 'k_{11}/k_{55}', 'KTOT' -> 'K_{TOT}'
    s = regexprep(name, 'KTOT', 'K_{TOT}'); s = regexprep(s, 'PTOT', 'P_{TOT}');
    s = regexprep(s, 'k(\d+)', 'k_{$1}');
end

function L = colLabel(nm)
% tick labels in the colour of each feature (feature_colors.m), as in every bar chart
    C = feature_colors(nm); L = cell(size(nm));
    for i = 1:numel(nm)
        L{i} = sprintf('\\color[rgb]{%.3f,%.3f,%.3f}\\bf %s', C(i,1), C(i,2), C(i,3), texLabel(nm{i}));
    end
end

function s = texName(name)
    s = regexprep(name, 'KTOT', '\\KTOT'); s = regexprep(s, 'PTOT', '\\PTOT');
    s = regexprep(s, 'k(\d+)', 'k_{$1}'); s = ['$' s '$'];
end

function s = arrow(d)
    if d == '+', s = '$\uparrow$'; else, s = '$\downarrow$'; end
end

function saveBoth(fig, figDir, name)
    try
        exportgraphics(fig, fullfile(figDir, [name '.png']), 'Resolution', 300);
        exportgraphics(fig, fullfile(figDir, [name '.pdf']), 'ContentType', 'vector');
    catch
        print(fig, fullfile(figDir, [name '.png']), '-dpng', '-r300');
        print(fig, fullfile(figDir, [name '.pdf']), '-dpdf', '-painters');
    end
    savefig(fig, fullfile(figDir, [name '.fig']));
end
