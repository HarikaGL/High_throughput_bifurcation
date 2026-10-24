% =========================================================================
%  PUBLICATION-READY ENSEMBLE BIFURCATION DIAGRAMS - ALL 10 CLASSES (v3)
%  For the output of Ensemble_LHS_v3.m (exact labels, inverse method).
%  Run it from the Ensemble_v3 folder (paths are set automatically).
%  Also works on Ensemble_LHS_v2 output: the v2 labels are converted to the ten
%  classes on the fly with classify_v3.m (identical to what Ensemble_LHS_v3 writes).
%
%  Reads  <outputDir>/MAT_Files/Complete_Bifurcation_Results_v3.mat (or ..._v2.mat)
%  and draws, for every class present, representative exact bifurcation
%  diagrams x2 (= APP) versus ATOT:
%    * every branch of steady states stored in res(k).curves
%        - 'main' / 'sym'        : branch connected to low ATOT  (class colour)
%        - 'detached' / 'offsym' : detached branch / crossing branch (darker shade)
%    * solid = stable, dashed = unstable (eigenvalues of the analytic Jacobian)
%    * linear axes; each panel's ATOT range is chosen from its singularities
%    * red circles = limit points (LP), green squares = branch points (BP),
%      magenta triangle = maximum of a biphasic branch
%    * optional: MATCONT validation curves (res(k).mc_curves), thin black dotted
%
%  Set outputDirs to one iteration folder (one figure) or to several
%  (one figure per folder, plus one combined figure with all classes).
% =========================================================================
clear; close all; clc;

%% ---------- PATHS (automatic: this script sits in Ensemble_v3, helpers in support/) ----
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
addpath(fullfile(HERE, 'support')); cd(HERE);

%% ---------- USER SETTINGS ----------
% result folders (relative to Ensemble_v3, or full paths)
outputDirs = {'Ensemble_iter1_v3_1e6', 'Ensemble_iter2_v3_1e6', ...
              'Ensemble_iter3_v3_1e6', 'Ensemble_iter4_v3_1e6'};
% (to use the existing v2 runs instead, give their full paths, e.g.
%  fullfile(fileparts(fileparts(HERE)), 'Ensemble_iter1_v2_1e6') )
makeCombined     = true;   % one extra figure with all classes from all folders
maxCurves        = 3;      % representative parameter sets per class panel
useLogX          = false;  % linear axes (set true for log axes)
useLogY          = false;
xLimits          = [];     % [] = automatic per panel (see below); or fixed, e.g. [0 500]
xAutoFactor      = 2.5;    % automatic x-range: [0, xAutoFactor * largest LP/BP (or peak) ATOT]
xDefaultMax      = 200;    % automatic x-range for panels without any singularity
showMatcont      = false;  % overlay MATCONT validation curves when available
skipArtefacts    = true;   % do not draw LP+BP sets whose LP coincides with the BP (flagged;
                           % they are still counted in the panel titles)
rngSeed          = 42;     % reproducible choice of representatives
outFigBase       = 'Ensemble_Bifurcation_All_Classes_v3';
figResolution    = 300;

%% ---------- CLASSES (ten classes, support/class_info_v3.m) ----------
CI = class_info_v3();
classTypeStrings = CI.names(:); classLabels = CI.labels(:); classColors = CI.colors;
nC = numel(classTypeStrings);

%% ---------- LOAD AND PLOT EACH FOLDER ----------
resComb = []; selComb = cell(nC,1); countComb = zeros(nC,1); win = [0 Inf];
for d = 1:numel(outputDirs)
    outputDir  = outputDirs{d};
    if ~exist(outputDir, 'dir'), outputDir = fullfile(HERE, outputDir); end
    [~, dirName] = fileparts(outputDir);
    resultFile = fullfile(outputDir, 'MAT_Files', 'Complete_Bifurcation_Results_v3.mat');
    if ~exist(resultFile, 'file')
        resultFile = fullfile(outputDir, 'MAT_Files', 'Complete_Bifurcation_Results_v2.mat');
    end
    if ~exist(resultFile, 'file')
        warning('No results file in %s, skipped', outputDir); continue;
    end
    fprintf('Loading %s ...\n', resultFile);
    S = load(resultFile, 'res', 'cfg');
    res = S.res(:);
    if isfield(S, 'cfg') && isfield(S.cfg, 'ATOT_window'), win = S.cfg.ATOT_window; end

    classIndices = classIndicesOf(res, classTypeStrings, skipArtefacts);
    classAll     = classIndicesOf(res, classTypeStrings, false);
    for c = 1:nC
        fprintf('  %d %-22s %6d\n', c, classTypeStrings{c}, numel(classIndices{c}));
    end

    rng(rngSeed);
    sel = cell(nC,1);
    for c = 1:nC, sel{c} = pickRepresentatives(classIndices{c}, maxCurves); end

    counts = cellfun(@numel, classAll);
    fig = drawClassPanels(res, sel, counts, classLabels, classColors, ...
        useLogX, useLogY, xLimits, [xAutoFactor xDefaultMax], showMatcont, ...
        sprintf('Ensemble bifurcation diagrams - %s', strrep(dirName, '_', '\_')));
    saveFigure(fig, fullfile(outputDir, 'Figures'), sprintf('%s_%s', outFigBase, dirName), figResolution);

    % keep only the selected parameter sets (not the whole file) for the combined figure
    if makeCombined
        for c = 1:nC
            countComb(c) = countComb(c) + counts(c);
            if isempty(selComb{c}) && ~isempty(sel{c})     % first folder that has this class
                resComb = [resComb; res(sel{c})]; %#ok<AGROW>
                selComb{c} = (numel(resComb)-numel(sel{c})+1 : numel(resComb))';
            end
        end
    end
    clear S res
end

%% ---------- COMBINED FIGURE (all classes, all folders) ----------
if makeCombined && numel(outputDirs) > 1 && ~isempty(resComb)
    fig = drawClassPanels(resComb, selComb, countComb, classLabels, classColors, useLogX, useLogY, ...
        xLimits, [xAutoFactor xDefaultMax], showMatcont, 'Ensemble bifurcation diagrams - all iterations');
    saveFigure(fig, fullfile(HERE, 'Figures_combined_v3'), outFigBase, figResolution);
end
fprintf('Done.\n');

%% ============== LOCAL FUNCTIONS ==============
function idx = classIndicesOf(res, classTypeStrings, skipArtefacts)
% ten-class labels from classify_v3 (same rule as Ensemble_LHS_v3), computed once
n = numel(res); lab = cell(n,1); art = false(n,1);
for i = 1:n
    [~, lab{i}, fl] = classify_v3(res(i), 'exact');
    art(i) = fl.lp_coincides_bp;
end
idx = cell(numel(classTypeStrings),1);
for c = 1:numel(classTypeStrings)
    k = find(strcmp(lab, classTypeStrings{c}));
    if skipArtefacts, k = k(~art(k)); end
    idx{c} = k;
end
end

function idxOut = pickRepresentatives(idxIn, maxN)
if isempty(idxIn), idxOut = []; return; end
if numel(idxIn) <= maxN, idxOut = idxIn(:); return; end
p = randperm(numel(idxIn));
idxOut = idxIn(p(1:maxN)); idxOut = idxOut(:);
end

function xl = panelXLim(rs, xLimits, xAuto, useLogX)
% x-range of one panel: fixed if given; otherwise from the singularities of the
% parameter sets drawn (LP, BP, or the peak of a biphasic branch)
if ~isempty(xLimits), xl = xLimits; return; end
m = 0;
for r = rs(:)'
    m = max([m; r.exact_LP_ATOT(:); r.exact_BP_ATOT(:)]);
    if strcmp(char(r.exact_label), 'NO_SING_biphasic') && ~isempty(r.curves)
        types = arrayfun(@(cb) char(cb.type), r.curves, 'UniformOutput', false);
        b = find(strcmp(types, 'main'), 1); if isempty(b), b = 1; end
        [~, i] = max(r.curves(b).x2); m = max(m, r.curves(b).ATOT(i));
    end
end
if m > 0, xl = [0 xAuto(1)*m]; else, xl = [0 xAuto(2)]; end
if useLogX, xl(1) = xl(2)*1e-4; end
end

function fig = drawClassPanels(res, sel, counts, classLabels, classColors, ...
        useLogX, useLogY, xLimits, xAuto, showMatcont, figTitle)
active = find(~cellfun(@isempty, sel));
nPlots = numel(active);
if nPlots == 0, fig = figure; return; end
nCols = min(4, nPlots); nRows = ceil(nPlots / nCols);
fig = figure('Color','w', 'Position',[60 60 420*nCols 360*nRows], ...
    'Name', figTitle, 'NumberTitle','off');
for i = 1:nPlots
    c = active(i);
    ax = subplot(nRows, nCols, i); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    xl = panelXLim(res(sel{c}), xLimits, xAuto, useLogX);
    for k = sel{c}(:)'
        r = res(k);
        if showMatcont, plotMatcont(ax, r, xl, useLogX || useLogY); end
        plotBranches(ax, r, classColors(c,:), 2, xl, useLogX || useLogY);
        addMarkers(ax, r, c == 9);   % class 9 = biphasic
    end
    if useLogX, set(ax, 'XScale', 'log'); end
    if useLogY, set(ax, 'YScale', 'log'); end
    xlim(ax, xl);
    title(ax, sprintf('%d  %s  (n = %d)', c, classLabels{c}, counts(c)), 'FontSize',11, 'FontWeight','bold');
    xlabel(ax, 'A_{TOT}', 'FontSize',10);
    ylabel(ax, 'x_2 = [APP]', 'FontSize',10);
end
% legend (once, in the first panel)
ax1 = subplot(nRows, nCols, 1); hold(ax1,'on');
h = [plot(ax1, nan, nan, 'k-', 'LineWidth',2), plot(ax1, nan, nan, 'k--', 'LineWidth',2), ...
     plot(ax1, nan, nan, 'o', 'MarkerFaceColor','r', 'MarkerEdgeColor','k'), ...
     plot(ax1, nan, nan, 's', 'MarkerFaceColor','g', 'MarkerEdgeColor','k')];
legend(ax1, h, {'stable','unstable','LP','BP'}, 'Location','best', 'FontSize',8, 'Box','off');
try
    sgtitle(figTitle, 'FontSize',15, 'FontWeight','bold');
catch
end
drawnow;
end

function plotBranches(ax, r, col, lw, xl, logAxes)
% all branches of steady states; solid = stable, dashed = unstable
if ~isfield(r, 'curves') || isempty(r.curves), return; end
dark = 0.55*col;                                       % detached / crossing branch
for b = 1:numel(r.curves)
    cb = r.curves(b);
    A = cb.ATOT(:); x = cb.x2(:);
    if isempty(A), continue; end
    if isfield(cb, 'stable') && ~isempty(cb.stable), st = logical(cb.stable(:)); else, st = true(size(A)); end
    t = char(cb.type);
    if any(strcmp(t, {'detached','offsym'})), c = dark; else, c = col; end
    ok = isfinite(A) & isfinite(x) & A >= 0 & x >= 0 & A <= xl(2);   % clip to the panel range
    if logAxes, ok = ok & A > 0 & x > 0; end
    plotStability(ax, A(ok), x(ok), st(ok), c, lw);
end
end

function plotStability(ax, A, x, st, col, lw)
% stable and unstable parts as separate polylines; neighbouring points are
% shared so that the two styles join without gaps
ys = x; ys(~st) = NaN;
yu = x; yu(st)  = NaN;
sw = find(diff(st) ~= 0);                 % switches between stable and unstable
for j = sw(:)'
    yu(j) = x(j); yu(j+1) = x(j+1);
end
plot(ax, A, ys, '-',  'Color', col, 'LineWidth', lw);
plot(ax, A, yu, '--', 'Color', col, 'LineWidth', lw);
end

function plotMatcont(ax, r, xl, logAxes)
if ~isfield(r, 'mc_ran') || ~r.mc_ran || ~isfield(r, 'mc_curves') || isempty(r.mc_curves), return; end
for b = 1:numel(r.mc_curves)
    A = r.mc_curves(b).ATOT(:); x = r.mc_curves(b).x2(:);
    ok = A >= 0 & x >= 0 & A <= xl(2);
    if logAxes, ok = ok & A > 0 & x > 0; end
    plot(ax, A(ok), x(ok), ':', 'Color', [0 0 0], 'LineWidth', 1);
end
end

function addMarkers(ax, r, isBiphasic)
if ~isempty(r.exact_LP_ATOT)
    plot(ax, r.exact_LP_ATOT, r.exact_LP_x2, 'o', 'MarkerSize',7, ...
        'MarkerFaceColor','r', 'MarkerEdgeColor','k');
end
if ~isempty(r.exact_BP_ATOT)
    plot(ax, r.exact_BP_ATOT, r.exact_BP_x2, 's', 'MarkerSize',7, ...
        'MarkerFaceColor','g', 'MarkerEdgeColor','k');
end
if isBiphasic && ~isempty(r.curves)
    % maximum of x2 along the main branch (exact branch: no findpeaks needed)
    types = arrayfun(@(cb) char(cb.type), r.curves, 'UniformOutput', false);
    b = find(strcmp(types, 'main'), 1); if isempty(b), b = 1; end
    [pk, i] = max(r.curves(b).x2);
    plot(ax, r.curves(b).ATOT(i), pk, '^', 'MarkerSize',9, ...
        'MarkerFaceColor','m', 'MarkerEdgeColor','k');
end
end

function saveFigure(fig, saveDir, name, res)
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
p = fullfile(saveDir, name);
fprintf('Saving %s (.png, .pdf, .fig)\n', p);
try, savefig(fig, [p '.fig']); catch, end
try
    exportgraphics(fig, [p '.png'], 'Resolution', res);
    exportgraphics(fig, [p '.pdf'], 'ContentType', 'vector');
catch
    print(fig, p, '-dpng', sprintf('-r%d', res));
end
end
