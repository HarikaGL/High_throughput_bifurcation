% =================== ML_v3_misclassified_figure17.m ===================================
% Fig 17 of the manuscript as ONE figure: representative misclassified parameter sets
% for the six most frequent confusions of the pooled LightGBM classifier, panels (A)-(F),
% one row per confusion (three misclassified sets + one correctly classified set of
% the predicted class).
%
% The selection of the confusions and of the drawn sets is IDENTICAL to
% ML_v3_misclassified_diagrams.m (same data, split, model, seed, pickSeed, pickMode,
% nPerPair), so the panels show the same parameter sets as before. Only the layout
% is new: the figure is drawn at its printed size (the text width of the manuscript,
% 16.6 cm), so the font sizes below are the sizes on paper.
%
%   * one heading per row:  (A)  Class 4 (...) predicted as class 5 (...): n of N test sets (%)
%   * one short line above each panel: iteration and sample (grey for the reference set)
%   * y label only in the first column, x label only in the bottom row,
%     at most 4 labelled decades per log axis, one legend for the whole figure
%   * no text inside the axes (nothing can overlap a curve); the singularities of every
%     drawn set are written to Fig17_panels.csv instead
%
% Run from the Ensemble_v3 folder (needs Python with lightgbm, as the ML script, and
% Ensemble_iter<i>_v3_1e6/MAT_Files/Complete_Bifurcation_Results_v3.mat).
% The first run fits the model and saves everything the figure needs to
% Fig17_cache.mat; later runs (e.g. after changing fonts) load the cache and take
% seconds. Set useCache = false to recompute.
% Output (ML_v3_results_parallel/labels_Exact_Class_v3/misclassified/):
%   Fig17.pdf (vector, 16.6 x 19.5 cm)  Fig17.png (600 dpi)  Fig17.fig  Fig17_panels.csv
% and, if the folder exists, a copy of Fig17.pdf / Fig17.png in
%   ../../../Manuscript_Final/npj_v2/1_Manuscript/figures/
% ===================================================================================
clear; close all; clc;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
addpath(fullfile(HERE, 'support')); cd(HERE);

%% ---------------- USER SETTINGS ----------------------------------------------------
useCache      = true;        % reuse Fig17_cache.mat if it exists (no model refit)
% selection (keep as in ML_v3_misclassified_diagrams.m to draw the same sets)
evalSet       = 'test';
model         = 'LightGBM';
nPerPair      = 3;           % misclassified sets per row
nRows         = 6;           % confusions shown, largest first -> (A)...(F)
showReference = true;        % 4th column: correctly classified set of the predicted class
pickSeed      = 11;
pickMode      = 'spread';
logZoom       = 30;          % log-log window: min(singular ATOT)/logZoom .. max*logZoom
% layout: sizes on paper
FIG_W   = 16.6;              % cm, = \textwidth of the manuscript (472 pt)
FIG_H   = 19.5;              % cm, leaves room for the caption on the same page
FONT    = 'Helvetica';
FS_LAB  = 10;                % (A)...(F)
FS_HEAD = 8;                 % row heading
FS_SUB  = 7;                 % iteration / sample above each panel
FS_AX   = 7.5;               % tick labels
FS_AXL  = 8;                 % axis labels
FS_LEG  = 8;                 % legend
LW      = 1.3;               % branch line width (pt)
MS_LP   = 4.5;  MS_BP = 5;   % marker sizes (pt)
maxDec  = 4;                 % at most this many labelled decades per log axis
resultsRoot = fullfile(HERE, 'ML_v3_results_parallel', 'labels_Exact_Class_v3');
outDir      = fullfile(resultsRoot, 'misclassified');
cacheFile   = fullfile(outDir, 'Fig17_cache.mat');
msFigDir    = fullfile(HERE, '..', '..', '..', 'Manuscript_Final', 'npj_v2', '1_Manuscript', 'figures');

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

if ~exist(outDir, 'dir'), mkdir(outDir); end
CI = class_info_v3();

if useCache && isfile(cacheFile)
    fprintf('Loading %s (set useCache = false to refit the model)\n', cacheFile);
    load(cacheFile, 'P', 'acc', 'setName');
else
    %% ---------------- 1. data, exactly as the ML script -----------------------------
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
    clsNum = double(catStr);
    L = log10(X); p = size(L,2);
    pairs = nchoosek(1:p, 2);
    Rlog = L(:,pairs(:,1)) - L(:,pairs(:,2));
    names = [cfg.featNames, arrayfun(@(q) sprintf('%s/%s', cfg.featNames{pairs(q,1)}, ...
        cfg.featNames{pairs(q,2)}), 1:size(pairs,1), 'UniformOutput', false)];
    fsA = struct('name','A: 14 params + 91 log-ratios', 'X', [L, Rlog], 'names', {names}, 'jsd', false);
    sp = stratSplit(It*1000 + y, cfg.testFrac, cfg.seed + 1, cfg.valFrac);
    switch lower(evalSet)
        case 'test',       ev = sp.te; setName = 'test';
        case 'validation', ev = sp.va; setName = 'validation';
        case 'both',       ev = [sp.te; sp.va]; setName = 'test+validation';
        otherwise, error('evalSet must be test, validation or both');
    end

    %% ---------------- 2. fit on the training rows, predict ---------------------------
    res  = ml_v3_core('evalFS', fsA, y, sp.tr, ev, cfg, model, cfg.seed + 1);
    pred = res.pred(:); yt = y(ev);
    acc  = 100*mean(pred == yt);
    fprintf('%s accuracy on the %s rows: %.2f%% (ML_v3_summary.txt: 92.61%%)\n', model, setName, acc);

    %% ---------------- 3. confusion pairs ---------------------------------------------
    wrong = find(pred ~= yt); rows = ev(wrong);
    tTrue = clsNum(yt(wrong))'; tPred = clsNum(pred(wrong))';
    [pp, ~, gp] = unique([tTrue(:) tPred(:)], 'rows');
    cnt = accumarray(gp, 1);
    nTrueClass = arrayfun(@(c) sum(clsNum(yt) == c), pp(:,1));
    Tp = table(pp(:,1), pp(:,2), cnt, nTrueClass, 100*cnt./nTrueClass, 'VariableNames', ...
        {'TrueClass','PredictedClass','Count','NTrue','Pct'});
    Tp = sortrows(Tp, 'Count', 'descend');
    nWrong = numel(wrong);

    %% ---------------- 4. choose the sets (as ML_v3_misclassified_diagrams.m) ---------
    maxPairs = 12; nPairs = min(maxPairs, height(Tp));   % same RandStream sequence
    rsP = RandStream('mt19937ar', 'Seed', pickSeed);
    pick = cell(nPairs, 1); ref = zeros(nPairs, 1);
    for q = 1:nPairs
        in = find(tTrue(:) == Tp.TrueClass(q) & tPred(:) == Tp.PredictedClass(q));
        k = min(nPerPair, numel(in));
        if strcmpi(pickMode, 'random'), sel = in(randperm(rsP, numel(in), k));
        else, sel = in(round(linspace(1, numel(in), k))); end
        pick{q} = sel;
        if showReference
            ok = ev(pred == yt & clsNum(yt)' == Tp.PredictedClass(q));
            if ~isempty(ok), ref(q) = ok(randperm(rsP, numel(ok), 1)); end
        end
    end
    nR = min(nRows, nPairs);
    need = unique([cell2mat(cellfun(@(s) rows(s), pick(1:nR), 'UniformOutput', false)); ref(ref(1:nR) > 0)]);
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
    % everything the figure needs, in plain structs
    P = struct('trueClass',{},'predClass',{},'count',{},'nTrue',{},'pct',{},'sets',{});
    for q = 1:nR
        s = struct('R',{},'iter',{},'sample',{},'isRef',{});
        for i = 1:numel(pick{q})
            r = rows(pick{q}(i));
            s(end+1) = struct('R', Rsel(r), 'iter', It(r), 'sample', Smp(r), 'isRef', false); %#ok<AGROW>
        end
        if showReference && ref(q) > 0
            r = ref(q);
            s(end+1) = struct('R', Rsel(r), 'iter', It(r), 'sample', Smp(r), 'isRef', true); %#ok<AGROW>
        end
        P(q) = struct('trueClass', Tp.TrueClass(q), 'predClass', Tp.PredictedClass(q), ...
            'count', Tp.Count(q), 'nTrue', Tp.NTrue(q), 'pct', Tp.Pct(q), 'sets', s);
    end
    fprintf('%d misclassified %s sets; the %d rows shown account for %d of them\n', ...
        nWrong, setName, nR, sum([P.count]));
    save(cacheFile, 'P', 'acc', 'setName', '-v7.3');
    fprintf('saved %s\n', cacheFile);
end

%% ---------------- 5. the figure ------------------------------------------------------
nR = numel(P); nCol = nPerPair + showReference;
% vertical budget of one row (cm): heading | panel line | axes | tick labels | gap
H_HEAD = 0.46; H_SUB = 0.34; H_TICK = 0.40; H_GAP = 0.06;
TOP = 0.05; BOT = 1.05;                                  % bottom: x label + legend
rowH = (FIG_H - TOP - BOT)/nR;
AX_H = rowH - H_HEAD - H_SUB - H_TICK - H_GAP;
% horizontal: y label + tick labels left of every axes
LEFT0 = 1.25; GAP_X = 0.95; RIGHT = 0.35;
AX_W = (FIG_W - LEFT0 - RIGHT - (nCol-1)*GAP_X)/nCol;
fprintf('panel size %.2f x %.2f cm\n', AX_W, AX_H);

fig = figure('Color','w', 'Units','centimeters', 'Position',[2 1 FIG_W FIG_H], ...
    'PaperUnits','centimeters', 'PaperSize',[FIG_W FIG_H], 'PaperPosition',[0 0 FIG_W FIG_H], ...
    'Renderer','painters', 'InvertHardcopy','off');
% overlay for all text outside the axes; data units = cm from the bottom-left corner
ov = axes(fig, 'Units','centimeters', 'Position',[0 0 FIG_W FIG_H], 'XLim',[0 FIG_W], ...
    'YLim',[0 FIG_H], 'Visible','off', 'HitTest','off'); hold(ov, 'on');

info = {};                                                % rows of Fig17_panels.csv
for q = 1:nR
    yTop = FIG_H - TOP - (q-1)*rowH;                     % top of this row block
    % --- row heading
    cT = P(q).trueClass; cP = P(q).predClass;
    head = sprintf(['\\fontsize{%g}\\bf(%c)\\rm\\fontsize{%g}   Class %d (%s) ', ...
        'predicted as class %d (%s):  %s of %s test sets (%.1f%%)'], ...
        FS_LAB, 'A'+q-1, FS_HEAD, cT, CI.labels{cT}, cP, CI.labels{cP}, ...
        commas(P(q).count), commas(P(q).nTrue), P(q).pct);
    ht = text(ov, 0.12, yTop - H_HEAD/2, head, 'FontName',FONT, 'Interpreter','tex', ...
        'VerticalAlignment','middle', 'HorizontalAlignment','left', 'Color',[0 0 0]);
    shrinkToWidth(ht, FIG_W - 0.2);
    % --- panels
    yAx = yTop - H_HEAD - H_SUB - AX_H;                  % bottom of the axes
    for i = 1:nCol
        xAx = LEFT0 + (i-1)*(AX_W + GAP_X);
        if i > numel(P(q).sets), continue; end
        s = P(q).sets(i);
        ax = axes(fig, 'Units','centimeters', 'Position',[xAx yAx AX_W AX_H]);
        drawPanel(ax, s.R, cT, cP, logZoom, LW, MS_LP, MS_BP);
        styleAxes(ax, FONT, FS_AX, maxDec);
        if i == 1, ylabel(ax, 'APP', 'FontSize', FS_AXL); end
        if q == nR, xlabel(ax, 'A_{TOT}', 'FontSize', FS_AXL); end
        % panel line above the axes
        if s.isRef
            sub = sprintf('reference: iter. %d, sample %d', s.iter, s.sample);
            col = [0.42 0.42 0.42]; ang = 'italic';
        else
            sub = sprintf('iter. %d, sample %d', s.iter, s.sample);
            col = [0.15 0.15 0.15]; ang = 'normal';
        end
        hs = text(ov, xAx + AX_W/2, yAx + AX_H + H_SUB*0.45, sub, 'FontName',FONT, ...
            'FontSize',FS_SUB, 'FontAngle',ang, 'Color',col, 'Interpreter','none', ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle');
        shrinkToWidth(hs, AX_W + GAP_X - 0.15);
        info(end+1,:) = {char('A'+q-1), i, cT, cP, s.isRef, s.iter, s.sample, ...
            strjoin(describe(s.R), '; ')}; %#ok<AGROW>
    end
end

% --- one legend for the whole figure (bottom)
lg = axes(fig, 'Units','centimeters', 'Position',[0 0.05 FIG_W 0.45], 'Visible','off'); hold(lg,'on');
h = [plot(lg, nan, nan, '-', 'Color',[0 0.35 0.75], 'LineWidth',LW), ...
     plot(lg, nan, nan, '-', 'Color',[0.85 0.2 0.1], 'LineWidth',LW), ...
     plot(lg, nan, nan, '-', 'Color',[0.25 0.25 0.25], 'LineWidth',LW), ...
     plot(lg, nan, nan, '--','Color',[0.25 0.25 0.25], 'LineWidth',max(0.6, 0.7*LW)), ...
     plot(lg, nan, nan, 'o', 'MarkerFaceColor',[0.9 0.05 0.05], 'MarkerEdgeColor','k', 'MarkerSize',MS_LP+0.5), ...
     plot(lg, nan, nan, 's', 'MarkerFaceColor',[0 0.62 0.2], 'MarkerEdgeColor','k', 'MarkerSize',MS_BP+0.5)];
legend(lg, h, {'main branch','detached branch','stable','unstable','limit point','branch point'}, ...
    'Orientation','horizontal', 'Box','off', 'Location','north', 'FontSize',FS_LEG, 'FontName',FONT);

%% ---------------- 6. save -------------------------------------------------------------
base = fullfile(outDir, 'Fig17');
savefig(fig, [base '.fig']);
try, print(fig, [base '.pdf'], '-dpdf', '-vector'); catch, print(fig, [base '.pdf'], '-dpdf', '-painters'); end
print(fig, [base '.png'], '-dpng', '-r600');
Ti = cell2table(info, 'VariableNames', {'Panel','Column','TrueClass','PredictedClass', ...
    'IsReference','Iteration','Sample','Singularities'});
writetable(Ti, fullfile(outDir, 'Fig17_panels.csv'));
fprintf('saved %s.pdf / .png / .fig and Fig17_panels.csv\n', base);
if isfolder(msFigDir)
    copyfile([base '.pdf'], fullfile(msFigDir, 'Fig17.pdf'));
    copyfile([base '.png'], fullfile(msFigDir, 'Fig17.png'));
    fprintf('copied Fig17.pdf / .png to %s\n', msFigDir);
end

%% ============== LOCAL FUNCTIONS =====================================================
function drawPanel(ax, R, cTrue, cPred, logZoom, LW, msLP, msBP)
% exact bifurcation diagram of one set (as ML_v3_misclassified_diagrams.m): log-log
% when a detached branch or a branch point is present, linear zoom otherwise
hold(ax, 'on'); box(ax, 'on');
hasLog = R.n_LP_detached > 0 || ~isempty(R.exact_BP_ATOT) || any(cTrue == [2 4 5 6 7 8]) || any(cPred == [2 4 5 6 7 8]);
if hasLog && ~isempty([R.exact_LP_ATOT(:); R.exact_BP_ATOT(:)])
    [xl, yl] = diagram_window(R, 'log', logZoom); isLog = true;
else
    [xl, yl] = diagram_window(R, 'lin', logZoom); isLog = false;
end
draw_branches(ax, R.curves, LW);
if ~isempty(R.exact_LP_ATOT)
    plot(ax, R.exact_LP_ATOT, R.exact_LP_x2, 'o', 'MarkerFaceColor',[0.9 0.05 0.05], ...
        'MarkerEdgeColor','k', 'MarkerSize',msLP, 'LineWidth',0.5);
end
if ~isempty(R.exact_BP_ATOT)
    plot(ax, R.exact_BP_ATOT, R.exact_BP_x2, 's', 'MarkerFaceColor',[0 0.62 0.2], ...
        'MarkerEdgeColor','k', 'MarkerSize',msBP, 'LineWidth',0.5);
end
if isLog, set(ax, 'XScale','log', 'YScale','log'); end
if all(isfinite(xl)) && xl(2) > xl(1), xlim(ax, xl); end
if all(isfinite(yl)) && yl(2) > yl(1), ylim(ax, yl); end
end

function styleAxes(ax, font, fs, maxDec)
set(ax, 'FontName',font, 'FontSize',fs, 'FontWeight','normal', 'LineWidth',0.6, ...
    'TickLength',[0.025 0.025], 'TickDir','in', 'XMinorTick','off', 'YMinorTick','off', ...
    'Layer','top', 'XColor',[0.1 0.1 0.1], 'YColor',[0.1 0.1 0.1]);
for xy = 'XY'
    if strcmp(get(ax, [xy 'Scale']), 'log')
        thinLogTicks(ax, xy, maxDec);
    else
        set(ax, [xy 'Tick'], niceTicks(get(ax, [xy 'Lim']), 4));
        try, ax.([xy 'Axis']).Exponent = 0; catch, end   % no "x10^n" above the axes
    end
end
end

function thinLogTicks(ax, xy, maxN)
% label at most maxN decades, evenly spaced, so tick labels never collide
lim = log10(get(ax, [xy 'Lim']));
d = ceil(lim(1) - 1e-9):floor(lim(2) + 1e-9);
if isempty(d), return; end
step = ceil(numel(d)/maxN);
d = d(1:step:end);
set(ax, [xy 'Tick'], 10.^d, [xy 'TickLabel'], arrayfun(@(e) sprintf('10^{%d}', e), d, 'UniformOutput', false));
try, set(ax, 'TickLabelInterpreter', 'tex'); catch, end
end

function t = niceTicks(lim, maxN)
% at most maxN round tick values inside lim
r = diff(lim); if r <= 0 || ~isfinite(r), t = lim; return; end
raw = r/maxN; m = 10^floor(log10(raw)); f = raw/m;
if f <= 1, s = m; elseif f <= 2, s = 2*m; elseif f <= 2.5, s = 2.5*m; elseif f <= 5, s = 5*m; else, s = 10*m; end
t = (ceil(lim(1)/s)*s):s:(lim(2) + 1e-12*r);
end

function shrinkToWidth(h, wMax)
% reduce the font size of text h (cm data units) until it is at most wMax wide
for k = 1:12
    e = get(h, 'Extent');
    if e(3) <= wMax, return; end
    set(h, 'FontSize', get(h, 'FontSize') - 0.25);
    s = get(h, 'String');
    if ischar(s) && ~isempty(strfind(s, '\fontsize{'))            % scale explicit sizes too
        s = regexprep(s, '\\fontsize\{([\d.]+)\}', '\\fontsize{${num2str(str2double($1)*0.97)}}');
        set(h, 'String', s);
    end
end
end

function d = describe(R)
% singularities of one set in words (as ML_v3_misclassified_diagrams.m)
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

function s = commas(n)
s = sprintf('%d', round(n)); k = numel(s) - 3;
while k > 0, s = [s(1:k) ',' s(k+1:end)]; k = k - 3; end
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
