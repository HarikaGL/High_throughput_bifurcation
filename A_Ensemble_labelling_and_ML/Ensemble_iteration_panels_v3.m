% ======================= Ensemble_iteration_panels_v3.m =============================
% One figure per iteration (four figures), laid out as a table row:
%   | ITERATION i: Constraint | Bounds | DISTRIBUTION | REPRESENTATIVE ENSEMBLE BIFURCATIONS |
% Column 1: coloured iteration header + constraint and bounds (editable text)
% Column 2: heading "DISTRIBUTION" inside the cell; number of parameter sets per class
%           (ten classes), count and percentage above every bar, numbered key, total N
% Column 3: heading "REPRESENTATIVE ENSEMBLE BIFURCATIONS" inside the cell; up to
%           nCurves randomly chosen parameter sets of every class that occurs in the
%           iteration, overlaid (class colour; solid = stable, dashed = unstable;
%           red circle = LP, green square = BP, magenta triangle = maximum of a
%           biphasic response); linear or log-log axes per class
% All four figures have exactly the same size (fixed figure size, fixed layout in
% pixels, printed at a fixed resolution); no separate grey header row.
% (The previous layout is kept in Ensemble_iteration_panels_v3_old.m.)
% Column 1 is narrow (smaller font), so the bar plot gets more width; the count and
% percentage labels above the bars are shrunk automatically until no two overlap.
%
% Run from the Ensemble_v3 folder (paths set automatically). Reads
%   Ensemble_iter<i>_v3_1e6/MAT_Files/Complete_Bifurcation_Results_v3.mat
% Saves Figures_combined_v3/Iteration_<i>_panel_v3.png/.pdf/.fig
% ===================================================================================
clear; close all; clc;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
addpath(fullfile(HERE, 'support')); cd(HERE);

%% ---------------- USER SETTINGS ----------------------------------------------------
iterationsToDraw = 1:4;
nCurves      = 5;         % parameter sets overlaid per class panel
rngSeed      = 11;         % which sets are drawn
maxPanels    = 4;          % class panels in column 3
skipLPeqBP   = true;       % do not draw the flagged LP = BP sets (they stay in the bars)
W = 2200; H = 370;         % figure size in pixels (identical for all iterations)
iterColors = [0.66 0.82 0.52;   % iteration 1 (green)
              0.55 0.69 0.89;   % iteration 2 (blue)
              0.84 0.60 0.84;   % iteration 3 (violet)
              0.96 0.71 0.53];  % iteration 4 (orange)
constraintText = {'Unconstrained', 'Unconstrained', ...
                  {'Symmetric','(Pitchfork','bifurcation)'}, ...
                  {'Asymmetric','(Symmetry','breaking)'}};
boundsText = {{'Set 1','(S - dominant)','High K_{TOT}','Low P_{TOT}'}, ...
              {'Set 2','(Z - dominant)','Low K_{TOT}','High P_{TOT}'}, ...
              {'Set 3 (General)','Wide ranges','for all','parameters'}, ...
              {'Set 3 (General)','Wide ranges','for all','parameters'}};
keyText = {'LP S-shaped','LP S + detached','LP Z-shaped','LP Z + detached', ...
           'LP detached only','BP only','LP+BP (LP > BP)','LP+BP (LP < BP)', ...
           'Biphasic','Monotone','Other'};
panelTitles = {'LP S-shaped','LP S + detached','LP Z-shaped','LP Z + detached', ...
               'LP detached only','BP only','LP+BP (LP>BP)','LP+BP (LP<BP)', ...
               'Biphasic','Monotone','Other'};
iconScale = {'lin','log','lin','log','log','log','log','log','lin','lin','lin'};
logZoom   = 30;
LW.table = 2.5; LW.axes = 2.0; LW.curve = 1.8; LW.bar = 1.5;
FS.head = 15; FS.cell = 12;                  % column 1 only (narrow column)
FS.bar = 12; FS.barLbl = 9; FS.barLblMin = 6; FS.key = 10; FS.panel = 12; FS.top = 16;
% column edges (normalised): col 1 = Constraint | Bounds, col 2 = bars, col 3 = panels
xE = [0.004 0.078 0.148 0.445 0.996];
% layout in pixels
PX.margin  = 6;     % outer margin of the table
PX.head    = 40;    % height of the in-cell headings of columns 2 and 3
PX.h1      = 58;    % ITERATION i band
PX.h2      = 48;    % Constraint | Bounds band
PX.barBot  = 34;    % space below the bar axes (tick labels)
PX.panBot  = 58;    % space below the class panels (tick labels + x label)
PX.title   = 24;    % space above the class panels (panel title)
outName = 'Iteration_%d_panel_v3';

CI = class_info_v3(); nC = 10;     % classes 1..10 in the bars
nx = @(p) p/W; ny = @(p) p/H;      % pixels -> normalised figure units

for ITER = iterationsToDraw
    %% ---------------- data ----------------------------------------------------------
    f = fullfile(HERE, sprintf('Ensemble_iter%d_v3_1e6', ITER), 'MAT_Files', ...
                 'Complete_Bifurcation_Results_v3.mat');
    if ~exist(f, 'file'), warning('Missing %s, skipped', f); continue; end
    fprintf('Iteration %d: loading %s ...\n', ITER, f);
    S = load(f, 'res', 'exact_class_v3', 'lp_eq_bp');
    cls = S.exact_class_v3(:); flg = S.lp_eq_bp(:); N = numel(cls);
    counts = histc(cls(cls > 0), 1:nC); counts = counts(:)';
    present = find(counts > 0);
    present = present(1:min(maxPanels, end));

    %% ---------------- figure and table frame -----------------------------------------
    fig = figure('Color','w', 'Units','pixels', 'Position', [30 80 W H], ...
        'Name', sprintf('Iteration %d', ITER), 'NumberTitle','off', ...
        'PaperPositionMode','auto', 'InvertHardcopy','off');
    ax0 = axes('Parent', fig, 'Position', [0 0 1 1], 'XLim', [0 1], 'YLim', [0 1], 'Visible', 'off');
    hold(ax0, 'on');
    yT = 1 - ny(PX.margin); yB = ny(PX.margin);
    h1 = ny(PX.h1); h2 = ny(PX.h2); hH = ny(PX.head);
    col = iterColors(ITER,:);
    patch(ax0, xE([1 3 3 1]), [yT yT yT-h1 yT-h1], col, 'EdgeColor', 'none');
    patch(ax0, xE([1 3 3 1]), [yT-h1 yT-h1 yT-h1-h2 yT-h1-h2], col, 'EdgeColor', 'none');
    text(ax0, mean(xE([1 3])), yT - h1/2, sprintf('ITERATION %d', ITER), 'FontSize', FS.head, ...
        'FontWeight','bold', 'HorizontalAlignment','center');
    text(ax0, mean(xE([1 2])), yT - h1 - h2/2, 'Constraint', 'FontSize', FS.head, 'FontWeight','bold', ...
        'HorizontalAlignment','center');
    text(ax0, mean(xE([2 3])), yT - h1 - h2/2, 'Bounds', 'FontSize', FS.head, 'FontWeight','bold', ...
        'HorizontalAlignment','center');
    yBody = yT - h1 - h2;
    yTxt  = (yBody + yB)/2;                       % body text centred in the remaining cell
    text(ax0, mean(xE([1 2])), yTxt, constraintText{ITER}, 'FontSize', FS.cell, ...
        'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle');
    text(ax0, mean(xE([2 3])), yTxt, boundsText{ITER}, 'FontSize', FS.cell, ...
        'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle');
    % headings of columns 2 and 3, inside their cells
    text(ax0, mean(xE([3 4])), yT - hH/2, 'DISTRIBUTION', 'FontSize', FS.top, ...
        'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle');
    text(ax0, mean(xE([4 5])), yT - hH/2, 'REPRESENTATIVE ENSEMBLE BIFURCATIONS', 'FontSize', FS.top, ...
        'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle');
    k = [0 0 0];
    plot(ax0, xE([1 5 5 1 1]), [yT yT yB yB yT], '-', 'Color', k, 'LineWidth', LW.table);
    plot(ax0, xE([1 3]), [1 1]*(yT - h1), '-', 'Color', k, 'LineWidth', LW.table);
    plot(ax0, xE([1 3]), [1 1]*yBody, '-', 'Color', k, 'LineWidth', LW.table);
    plot(ax0, [1 1]*xE(2), [yT - h1, yB], '-', 'Color', k, 'LineWidth', LW.table);
    for j = 3:4, plot(ax0, [1 1]*xE(j), [yT yB], '-', 'Color', k, 'LineWidth', LW.table); end

    %% ---------------- column 2: bar plot ---------------------------------------------
    bBot = yB + ny(PX.barBot); bTop = yT - hH - ny(4);
    axB = axes('Parent', fig, 'Position', [xE(3)+0.046, bBot, xE(4)-xE(3)-0.054, bTop-bBot]);
    hold(axB, 'on'); box(axB, 'on'); grid(axB, 'on');
    for c = 1:nC
        bar(axB, c, counts(c), 0.75, 'FaceColor', CI.colors(c,:), 'EdgeColor', 'k', 'LineWidth', LW.bar);
    end
    yMax = max(counts); yl = [0 max(10, 1.9*yMax)];
    set(axB, 'XLim', [0.4 nC+0.6], 'YLim', yl);
    hL = gobjects(nC, 1);
    for c = 1:nC
        if counts(c) == 0, lbl = {'0'}; else, lbl = {sprintf('(%.1f%%)', 100*counts(c)/N), sprintf('%d', counts(c))}; end
        if any(flg(cls == c)), lbl{2} = sprintf('%d*', counts(c)); end
        hL(c) = text(axB, c, counts(c) + 0.02*yl(2), lbl, 'FontSize', FS.barLbl, 'FontWeight','bold', ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom');
    end
    fit_bar_labels(hL, FS.barLbl, FS.barLblMin);
    key = arrayfun(@(c) sprintf('%d - %s', c, keyText{c}), 1:nC, 'UniformOutput', false);
    text(axB, 0.02, 0.98, key, 'Units','normalized', 'FontSize', FS.key - 1, 'FontWeight','bold', ...
        'VerticalAlignment','top', 'HorizontalAlignment','left');
    text(axB, 0.98, 0.98, sprintf('Total N = %d', N), 'Units','normalized', 'FontSize', FS.key + 1, ...
        'FontAngle','italic', 'FontWeight','bold', 'HorizontalAlignment','right', 'VerticalAlignment','top');
    if any(flg)
        text(axB, 0.98, 0.88, sprintf('* incl. %d LP = BP', sum(flg)), 'Units','normalized', ...
            'FontSize', FS.key, 'HorizontalAlignment','right', 'VerticalAlignment','top');
    end
    set(axB, 'XLim', [0.4 nC+0.6], 'YLim', yl, 'XTick', 1:nC, 'FontSize', FS.bar, 'FontWeight','bold', ...
        'LineWidth', LW.axes, 'Layer', 'top');
    ylabel(axB, 'Number of Parameter Sets', 'FontSize', FS.bar + 1, 'FontWeight','bold');

    %% ---------------- column 3: representative bifurcations --------------------------
    rs = RandStream('mt19937ar', 'Seed', rngSeed + ITER);
    nP = numel(present);
    w = (xE(5) - xE(4))/maxPanels;
    pBot = yB + ny(PX.panBot);
    pTopMax = yT - hH - ny(PX.title);
    pw = w - 0.042;                               % panel width (normalised)
    ph = min(pTopMax - pBot, pw*W/H);             % square panels where possible
    pw = ph*H/W;                                  % keep them square if height limits
    for p = 1:nP
        c = present(p);
        cand = find(cls == c & ~(skipLPeqBP & flg));
        allFlagged = isempty(cand);
        if allFlagged, cand = find(cls == c); end
        sel = cand(randperm(rs, numel(cand), min(nCurves, numel(cand))));
        Rs = S.res(sel);
        x0 = xE(4) + (p-1)*w + (maxPanels - nP)*w/2;
        axP = axes('Parent', fig, 'Position', [x0 + 0.032 + (w - 0.042 - pw)/2, ...
                   pBot + (pTopMax - pBot - ph)/2, pw, ph]);
        hold(axP, 'on'); box(axP, 'on');
        [xl, ylim2, isLog] = class_window(Rs, iconScale{c}, logZoom);
        lc = 0.8*CI.colors(c,:);                          % slightly darker class colour
        for R = Rs(:)'
            draw_overlay(axP, R, lc, LW.curve, isLog, c == 9);
        end
        if isLog, set(axP, 'XScale','log', 'YScale','log'); end
        xlim(axP, xl); ylim(axP, ylim2);
        set(axP, 'LineWidth', LW.axes, 'FontSize', FS.panel - 2, 'FontWeight','bold', 'Layer','top');
        if isLog
            set(axP, 'XTick', 10.^(ceil(log10(xl(1))):floor(log10(xl(2)))), ...
                     'YTick', 10.^(ceil(log10(ylim2(1))):2:floor(log10(ylim2(2)))));
        else
            xl = [xl(1) nice_up(xl(2))]; ylim2 = [ylim2(1) nice_up(ylim2(2))];
            xlim(axP, xl); ylim(axP, ylim2);
            set(axP, 'XTick', xl, 'YTick', ylim2);
        end
        t = sprintf('%s (n=%d)', panelTitles{c}, numel(sel));
        if allFlagged, t = [t ', LP = BP']; end
        title(axP, t, 'FontSize', FS.panel, 'FontWeight','bold');
        xlabel(axP, 'A_{TOT}', 'FontSize', FS.panel, 'FontWeight','bold');
        ylabel(axP, 'APP_{steady state}', 'FontSize', FS.panel - 1, 'FontWeight','bold');
    end

    %% ---------------- save ------------------------------------------------------------
    outDir = fullfile(HERE, 'Figures_combined_v3');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    base = fullfile(outDir, sprintf(outName, ITER));
    savefig(fig, [base '.fig']);
    drawnow;
    print(fig, [base '.png'], '-dpng', '-r300');      % fixed size: identical for all iterations
    try
        exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
    catch
        print(fig, [base '.pdf'], '-dpdf', '-painters');
    end
    fprintf('Saved %s (.png, .pdf, .fig)\n', base);
    clear S
end

%% ============== LOCAL FUNCTIONS =====================================================
function fit_bar_labels(h, fs, fsMin)
% shrink the labels above the bars until neighbouring labels are separated by a gap
% (labels are centred on x = 1, 2, ...; Extent is in data units)
gap = 0.06;
while true
    drawnow;
    E = get(h, 'Extent'); E = vertcat(E{:});
    hw = E(:,3)/2;                                % half-widths
    if all(hw(1:end-1) + hw(2:end) + gap <= 1) || fs <= fsMin, break; end
    fs = fs - 0.5; set(h, 'FontSize', fs);
end
end

function [xl, yl, isLog] = class_window(Rs, scaleMode, logZoom)
% common window of several samples: median of the per-sample windows (diagram_window),
% geometric median for log axes
L = zeros(numel(Rs), 4); lg = false(numel(Rs), 1);
for i = 1:numel(Rs)
    [a, b, lg(i)] = diagram_window(Rs(i), scaleMode, logZoom);
    L(i,:) = [a b];
end
isLog = any(lg) && strcmpi(scaleMode, 'log');
if isLog
    L = L(lg,:);
    xl = 10.^[median(log10(L(:,1))), median(log10(L(:,2)))];
    yl = 10.^[median(log10(L(:,3))), median(log10(L(:,4)))];
else
    xl = [min(L(:,1)), median(L(:,2))];
    yl = [min(L(:,3)), median(L(:,4))];
end
if ~(xl(2) > xl(1)), xl(2) = xl(1)*10 + 1; end
if ~(yl(2) > yl(1)), yl(2) = yl(1)*10 + 1; end
end

function v = nice_up(v)
% round up to 1, 2, 2.5 or 5 x 10^k (clean axis limits)
if ~(v > 0), return; end
e = 10^floor(log10(v)); m = [1 2 2.5 5 10];
v = e*m(find(m*e >= v*(1 - 1e-9), 1));
end

function draw_overlay(ax, R, lc, lw, isLog, isBiphasic)
% one sample: solid = stable, dashed = unstable, LP / BP / peak markers
for cb = R.curves(:)'
    A = cb.ATOT(:)'; X = cb.x2(:)';
    if isempty(A), continue; end
    st = logical(cb.stable(:)');
    ys = X; ys(~st) = NaN;
    yu = X; yu(st) = NaN;
    sw = find(diff(st) ~= 0);
    for j = sw, yu(j) = X(j); yu(j+1) = X(j+1); end
    if isLog, A(A <= 0) = NaN; ys(ys <= 0) = NaN; yu(yu <= 0) = NaN; end
    plot(ax, A, ys, '-',  'Color', lc, 'LineWidth', lw);
    plot(ax, A, yu, '--', 'Color', lc, 'LineWidth', lw);
end
if ~isempty(R.exact_LP_ATOT)
    plot(ax, R.exact_LP_ATOT, R.exact_LP_x2, 'o', 'MarkerFaceColor', [0.9 0 0], ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 5);
end
if ~isempty(R.exact_BP_ATOT)
    plot(ax, R.exact_BP_ATOT, R.exact_BP_x2, 's', 'MarkerFaceColor', [0 0.7 0], ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 5);
end
if isBiphasic && ~isempty(R.curves)
    types = arrayfun(@(cb) char(cb.type), R.curves, 'UniformOutput', false);
    b = find(strcmp(types, 'main'), 1); if isempty(b), b = 1; end
    [pk, i] = max(R.curves(b).x2);
    plot(ax, R.curves(b).ATOT(i), pk, '^', 'MarkerFaceColor', [0.85 0 0.85], ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 6);
end
end
