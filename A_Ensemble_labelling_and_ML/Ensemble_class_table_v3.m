% ======================= Ensemble_class_table_v3.m ==================================
% Class table figure: | Parameter sets | Labels | Bifurcations |
%   column 1: the first three parameter sets of each class, in the order of the
%             Excel sheets (iterations 1-4 one after the other), then "......"
%   column 2: the class label
%   column 3: one real bifurcation diagram per class (x2 = APP vs ATOT), drawn from
%             the first parameter set listed for that class (solid black = stable,
%             dashed red = unstable; red circle = LP, green square = BP); no text labels
% All texts are MATLAB text objects, so they can be edited in the saved .fig.
%
% Run from the Ensemble_v3 folder (paths set automatically). Reads
%   Ensemble_iter<i>_v3_1e6/MAT_Files/Complete_Bifurcation_Results_v3.mat, i = 1..4
% Saves Figures_combined_v3/Class_Table_v3.png/.pdf/.fig
% ===================================================================================
clear; close all; clc;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
addpath(fullfile(HERE, 'support')); cd(HERE);

%% ---------------- USER SETTINGS ----------------------------------------------------
resultDirs  = {'Ensemble_iter1_v3_1e6','Ensemble_iter2_v3_1e6', ...
               'Ensemble_iter3_v3_1e6','Ensemble_iter4_v3_1e6'};
classes     = 1:10;                % rows of the table (class numbers of class_info_v3)
nListed     = 3;                   % parameter sets listed per class before "......"
numbering   = 'global';            % 'global'    : 1..40000 (iteration 1 first, as in the
                                   %               pooled ML dataset)
                                   % 'iteration' : "it2: 5" (row number in that
                                   %               iteration's Excel sheet)
skipLPeqBP  = true;                % do not list / draw the flagged LP = BP sets
exampleRank = ones(1, 11);         % which listed set is drawn per class (1 = first)
titleText   = 'Bifurcation classes';   % e.g. 'MATCONT'; '' = no title
panelLabel  = '';                  % e.g. '[III]'
rowLabels = {'LP ONLY (S shaped)','LP (S shaped) + DETACHED','LP ONLY (Z shaped)', ...
             'LP (Z shaped) + DETACHED','LP DETACHED ONLY','BP ONLY','LP+BP (LP>BP)', ...
             'LP+BP (LP<BP)','BIPHASIC','MONOTONE','OTHER'};
% axes of the small diagrams (same choice as the class wheel)
iconScale = {'lin','log','lin','log','log','log','log','log','lin','lin','lin'};
logZoom   = 30;
showMarkers = true;
gridRC    = [5 2];                 % rows x columns of the diagram grid
LW.table  = 2.5;  LW.curve = 3.0;  LW.axes = 2.0;     % line widths
FS.header = 18;   FS.cell = 15;    FS.axes = 10;  FS.title = 18;   % font sizes
headerColor = [1.00 0.93 0.70];
outName   = 'Class_Table_v3';

%% ---------------- 1. labels and the first sets of every class ----------------------
cls = []; flg = []; src = [];
for d = 1:numel(resultDirs)
    f = fullfile(HERE, resultDirs{d}, 'MAT_Files', 'Complete_Bifurcation_Results_v3.mat');
    if ~exist(f, 'file'), error('Missing %s', f); end
    S = load(f, 'exact_class_v3', 'lp_eq_bp');
    n = numel(S.exact_class_v3);
    cls = [cls; S.exact_class_v3(:)]; flg = [flg; S.lp_eq_bp(:)]; %#ok<AGROW>
    src = [src; [d*ones(n,1), (1:n)']]; %#ok<AGROW>
end
nR = numel(classes);
listed = cell(1, nR); setText = cell(1, nR); example = zeros(1, nR);
for i = 1:nR
    c = classes(i);
    k = find(cls == c & ~(skipLPeqBP & flg));
    listed{i} = k(1:min(nListed, end));
    if isempty(k), setText{i} = '-'; continue; end
    switch lower(numbering)
        case 'global'
            s = arrayfun(@(g) sprintf('%d', g), listed{i}, 'UniformOutput', false);
        otherwise
            s = arrayfun(@(g) sprintf('it%d: %d', src(g,1), src(g,2)), listed{i}, 'UniformOutput', false);
    end
    setText{i} = [strjoin(s(:)', ', '), ', ......'];
    example(i) = listed{i}(min(exampleRank(c), numel(listed{i})));
end

% bifurcation data of the drawn examples only
Rex = cell(1, nR);
for d = 1:numel(resultDirs)
    need = find(example > 0 & arrayfun(@(g) g > 0 && src(max(g,1),1) == d, example));
    if isempty(need), continue; end
    f = fullfile(HERE, resultDirs{d}, 'MAT_Files', 'Complete_Bifurcation_Results_v3.mat');
    fprintf('Loading %s ...\n', resultDirs{d});
    S = load(f, 'res');
    for i = need, Rex{i} = S.res(src(example(i), 2)); end
    clear S
end

%% ---------------- 2. table -----------------------------------------------------------
fig = figure('Color','w', 'Units','pixels', 'Position',[40 40 1500 1050], ...
    'Name','Class table', 'NumberTitle','off');
ax0 = axes('Parent', fig, 'Position', [0 0 1 1], 'XLim', [0 1], 'YLim', [0 1], 'Visible', 'off');
hold(ax0, 'on');
x = [0.03 0.25 0.53 0.97];             % column edges
yTop = 0.92; yBot = 0.03; hHead = 0.065;
yRows = linspace(yTop - hHead, yBot, nR + 1);   % row edges below the header
% header background and outer frame
patch(ax0, [x(1) x(4) x(4) x(1)], [yTop yTop yTop-hHead yTop-hHead], headerColor, 'EdgeColor', 'none');
headers = {'Parameter Sets', 'Labels', 'Bifurcations'};
for j = 1:3
    text(ax0, mean(x(j:j+1)), yTop - hHead/2, headers{j}, 'FontSize', FS.header, ...
        'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle');
end
for i = 1:nR
    yc = mean(yRows(i:i+1));
    text(ax0, mean(x(1:2)), yc, setText{i}, 'FontSize', FS.cell, ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');
    text(ax0, x(2) + 0.015, yc, rowLabels{classes(i)}, 'FontSize', FS.cell, ...
        'HorizontalAlignment','left', 'VerticalAlignment','middle');
end
% grid lines (row lines only in the first two columns, as in the original layout)
k = [0 0 0];
plot(ax0, [x(1) x(4) x(4) x(1) x(1)], [yTop yTop yBot yBot yTop], '-', 'Color', k, 'LineWidth', LW.table);
plot(ax0, [x(1) x(4)], [1 1]*(yTop - hHead), '-', 'Color', k, 'LineWidth', LW.table);
for j = 2:3, plot(ax0, [x(j) x(j)], [yTop yBot], '-', 'Color', k, 'LineWidth', LW.table); end
for i = 2:nR, plot(ax0, [x(1) x(3)], [1 1]*yRows(i), '-', 'Color', k, 'LineWidth', LW.table*0.6); end
if ~isempty(titleText)
    text(ax0, 0.5, 0.965, titleText, 'FontSize', FS.title, 'FontWeight','bold', ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');
end
if ~isempty(panelLabel)
    text(ax0, 0.01, 0.985, panelLabel, 'FontSize', FS.title, 'FontWeight','bold', ...
        'VerticalAlignment','top');
end

%% ---------------- 3. small bifurcation diagrams ------------------------------------
gx = linspace(x(3), x(4), gridRC(2) + 1); gy = linspace(yTop - hHead, yBot, gridRC(1) + 1);
mL = 0.055; mR = 0.012; mB = 0.035; mT = 0.012;   % margins inside each grid cell
for i = 1:nR
    r = ceil(i/gridRC(2)); cc = i - (r-1)*gridRC(2);
    if r > gridRC(1), break; end
    pos = [gx(cc) + mL, gy(r+1) + mB, (gx(cc+1) - gx(cc)) - mL - mR, (gy(r) - gy(r+1)) - mB - mT];
    ax = axes('Parent', fig, 'Position', pos); hold(ax, 'on');
    R = Rex{i};
    if ~isempty(R), draw_diagram(ax, R, iconScale{classes(i)}, logZoom, LW, showMarkers); end
    set(ax, 'XTick', [], 'YTick', [], 'Box', 'off', 'LineWidth', LW.axes, 'Color', 'none', ...
        'XColor', [0 0 0], 'YColor', [0 0 0], 'FontSize', FS.axes);
    xlabel(ax, 'ATOT', 'FontWeight','bold', 'FontSize', FS.axes);
    ylabel(ax, {'Steady state','(APP)'}, 'FontWeight','bold', 'FontSize', FS.axes);
end

%% ---------------- 4. save ------------------------------------------------------------
outDir = fullfile(HERE, 'Figures_combined_v3');
if ~exist(outDir, 'dir'), mkdir(outDir); end
base = fullfile(outDir, outName);
savefig(fig, [base '.fig']);
try
    exportgraphics(fig, [base '.png'], 'Resolution', 400);
    exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
catch
    print(fig, base, '-dpng', '-r400');
end
fprintf('Saved %s (.png, .pdf, .fig)\n', base);
for i = 1:nR, fprintf('  %-26s %s\n', rowLabels{classes(i)}, setText{i}); end

%% ============== LOCAL FUNCTIONS =====================================================
function draw_diagram(ax, R, scaleMode, logZoom, LW, showMarkers)
% stable = solid black, unstable = dashed red; linear or log-log axes (diagram_window)
[xl, yl, isLog] = diagram_window(R, scaleMode, logZoom);
for cb = R.curves(:)'
    A = cb.ATOT(:)'; X = cb.x2(:)';
    if isempty(A), continue; end
    st = logical(cb.stable(:)');
    ys = X; ys(~st) = NaN;
    yu = X; yu(st) = NaN;
    sw = find(diff(st) ~= 0);
    for j = sw, yu(j) = X(j); yu(j+1) = X(j+1); end
    if isLog, A(A <= 0) = NaN; ys(ys <= 0) = NaN; yu(yu <= 0) = NaN; end
    plot(ax, A, ys, '-',  'Color', [0 0 0],        'LineWidth', LW.curve);
    plot(ax, A, yu, '--', 'Color', [0.9 0.05 0.05], 'LineWidth', LW.curve);
end
if showMarkers
    if ~isempty(R.exact_LP_ATOT)
        plot(ax, R.exact_LP_ATOT, R.exact_LP_x2, 'o', 'MarkerFaceColor', [0.9 0 0], ...
            'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1);
    end
    if ~isempty(R.exact_BP_ATOT)
        plot(ax, R.exact_BP_ATOT, R.exact_BP_x2, 's', 'MarkerFaceColor', [0 0.65 0], ...
            'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1);
    end
end
if isLog, set(ax, 'XScale', 'log', 'YScale', 'log'); end
xlim(ax, xl); ylim(ax, yl);
end
