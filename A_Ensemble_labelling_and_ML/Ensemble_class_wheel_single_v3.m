% =================== Ensemble_class_wheel_single_v3.m ===============================
% Class wheel with ONE representative bifurcation diagram per class.
%   * the ten classes of the v3 ensemble get EQUAL wedges (36 degrees each),
%     class 1 at the top, running clockwise (OTHER is empty and not drawn)
%   * every wedge holds exactly one REAL exact bifurcation diagram (x2 = APP vs ATOT)
%     of one parameter set of that class: solid black = stable, dashed red = unstable
%   * limit points are marked with a red dot and the text "LP", branch points with a
%     green square and the text "BP"
%   * coloured outer band with the class name as ONE editable text object per class,
%     written straight and tilted along the band (tangent to the circle at the middle
%     of the wedge); centre circle with text; numbered legend on the right, placed
%     close to the wheel
%   * diagrams slightly smaller than the largest frame that fits the wedge and moved
%     towards the centre (frameScale, inwardShift); larger fonts inside the diagrams
% (The previous version is kept in Ensemble_class_wheel_single_v3_old.m.)
%
% Which parameter set is shown (per class, see USER SETTINGS):
%   pickMode = 'clear'  : among the sets of the class, the ones whose singular points
%                         are well separated (x2 and ATOT spread of the LPs / BPs) are
%                         ranked, and one set between the pLo and pHi quantiles of
%                         this score is drawn at random (seeded): clear, but typical,
%                         not an extreme. Classes with one singular point or none:
%                         random set.
%   pickMode = 'random' : a random set of the class (seeded).
%   manualPick{c} = [iteration sample] overrides the choice for class c.
% Flagged sets (LP coincides with BP) are never drawn.
%
% Run from the Ensemble_v3 folder (paths set automatically). Reads
%   Ensemble_iter<i>_v3_1e6/MAT_Files/Complete_Bifurcation_Results_v3.mat
% Saves Figures_combined_v3/Class_Wheel_single_v3.png/.pdf/.fig
% Each results file is loaded once (about 3-5 GB in memory); only the ten chosen
% parameter sets are kept. The chosen sets are printed, so a choice can be fixed in
% manualPick.
% ===================================================================================
clear; close all; clc;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
addpath(fullfile(HERE, 'support')); cd(HERE);

%% ---------------- USER SETTINGS ----------------------------------------------------
resultDirs  = {'Ensemble_iter1_v3_1e6','Ensemble_iter2_v3_1e6', ...
               'Ensemble_iter3_v3_1e6','Ensemble_iter4_v3_1e6'};
classesShown = 1:10;       % classes on the wheel (order = clockwise from the top)
rngSeed     = 3;           % change to draw other representative sets
pickMode    = 'clear';     % 'clear' or 'random' (see header)
pLo = 0.60; pHi = 0.90;    % quantile window of the clarity score ('clear')
% iteration each class is taken from (0 = the iteration with most sets of the class)
fromIter    = zeros(1, 10);
% manual choice, e.g. manualPick{5} = [2 1234];  (iteration, sample number)
manualPick  = cell(1, 10);
startAngleDeg = 90;        % where class 1 starts (90 = top); clockwise
tiltFrac    = 0;           % 0 = every diagram upright (most readable); 1 = turned fully
                           % with its wedge; e.g. 0.3 = slight tilt towards the wedge
% axes of the diagrams (order of class_info_v3):
%   'lin' : linear axes, zoomed on the folds
%   'log' : log-log axes, ATOT from min(singular ATOT)/logZoom to max*logZoom
%           (branches on very different scales: detached branch, BP, LP + BP)
iconScale = {'lin','log','lin','log','log','log','log','log','lin','lin','lin'};
logZoom   = 30;
markLP    = true;  markBP = true;  labelMarkers = true;   % LP / BP markers and text
xLab = 'A_{TOT}';  yLab = 'Steady state (APP)';
logSuffix = '';            % e.g. ' (log)' appended to the axis labels of log diagrams
% line widths (points) and fonts
LW.curve  = 3.0;           % bifurcation curves
LW.axes   = 2.0;           % diagram axes
LW.sector = 2.5;           % grey lines between the wedges
LW.band   = 3.0;           % white lines between the band segments
LW.circle = 3.5;           % centre circle and outer circle
MS.LP = 9; MS.BP = 9;      % marker sizes
FS.band = 14; FS.axis = 13; FS.marker = 13; FS.centre = 20; FS.legend = 15;
frameScale  = 0.86;        % diagram size relative to the largest frame that fits a wedge
inwardShift = 0.6;         % 0 = keep the frame where it fits best; 1 = move it inwards by
                           % the full size reduction (more compact, away from the band)
legendGap   = 0.005;       % gap between the wheel and the legend (figure units)
centerText  = {'Codimension-1','Bifurcations'};   % {} = empty centre
showLegend  = true;        % numbered colour circles with the class names, right side
panelLabel  = '';          % e.g. '[IV]'
outName     = 'Class_Wheel_single_v3';

% names on the band and in the legend (order of class_info_v3)
bandLabels = {'LP ONLY (S shaped)','LP (S) + DETACHED','LP ONLY (Z shaped)', ...
              'LP (Z) + DETACHED','LP DETACHED ONLY','BP ONLY', ...
              'LP + BP (LP > BP)','LP + BP (LP < BP)','BIPHASIC','MONOTONE','OTHER'};
legendLabels = {'LP ONLY (S shaped)','LP (S shaped) + detached LP','LP ONLY (Z shaped)', ...
                'LP (Z shaped) + detached LP','LP on detached branch only','BP ONLY', ...
                'LP + BP (LP > BP)','LP + BP (LP < BP)','Biphasic','Monotone','Other'};

%% ---------------- GEOMETRY ---------------------------------------------------------
rC = 0.27;                 % centre circle
rB = 0.86;                 % inner edge of the coloured band
rO = 1.00;                 % outer edge
CI = class_info_v3();
nW = numel(classesShown);
ang = 360/nW;                                        % equal wedges
th1 = startAngleDeg - (0:nW-1)*ang;                  % upper edge of every wedge
th2 = th1 - ang;
% diagram frames: largest square (plot + axis labels) that fits every wedge at its
% rotation; the same size in all wedges
mids = (th1 + th2)/2;
rots = tiltFrac*(mod(mids - 90 + 90, 180) - 90);
hF = zeros(1, nW); rhoF = zeros(1, nW);
for w = 1:nW, [hF(w), rhoF(w)] = frame_size(rC, rB, th2(w), th1(w), rots(w), 0.02); end
hFrame0 = min(hF);
hFrame  = frameScale*hFrame0;
rhoF    = rhoF - inwardShift*(hFrame0 - hFrame);   % move the smaller frames inwards

%% ---------------- 1. class labels of all iterations ---------------------------------
nIt = numel(resultDirs);
cls = cell(1, nIt); flg = cell(1, nIt);
for d = 1:nIt
    f = fullfile(HERE, resultDirs{d}, 'MAT_Files', 'Complete_Bifurcation_Results_v3.mat');
    if ~exist(f, 'file'), error('Missing %s', f); end
    S = load(f, 'exact_class_v3', 'lp_eq_bp');
    cls{d} = S.exact_class_v3(:); flg{d} = logical(S.lp_eq_bp(:));
end
src = zeros(1, nW);                                  % iteration of every class
for w = 1:nW
    c = classesShown(w);
    if ~isempty(manualPick{c}), src(w) = manualPick{c}(1); continue; end
    if fromIter(c) > 0, src(w) = fromIter(c); continue; end
    n = cellfun(@(x, g) sum(x == c & ~g), cls, flg);
    if all(n == 0), error('Class %d (%s) has no parameter set.', c, CI.names{c}); end
    [~, src(w)] = max(n);
end

%% ---------------- 2. one representative set per class -------------------------------
rs = RandStream('mt19937ar', 'Seed', rngSeed);
Rsel = cell(1, nW); idxSel = zeros(1, nW);
for d = unique(src)
    f = fullfile(HERE, resultDirs{d}, 'MAT_Files', 'Complete_Bifurcation_Results_v3.mat');
    fprintf('Loading %s ...\n', resultDirs{d});
    S = load(f, 'res');
    for w = find(src == d)
        c = classesShown(w);
        if ~isempty(manualPick{c})
            k = manualPick{c}(2);
        else
            cand = find(cls{d} == c & ~flg{d} & [S.res.exact_ok]');
            if isempty(cand), error('Class %d has no usable set in %s.', c, resultDirs{d}); end
            k = choose_set(S.res(cand), cand, pickMode, pLo, pHi, rs);
        end
        Rsel{w} = S.res(k); idxSel(w) = k;
        fprintf('  class %2d %-20s : iteration %d, sample %d\n', c, CI.names{c}, d, k);
    end
    clear S
end

%% ---------------- 3. draw ------------------------------------------------------------
figH = 1100; if showLegend, figW = 1650; axPos = [0.01 0.02 0.64 0.96]; else, figW = 1100; axPos = [0.02 0.02 0.96 0.96]; end
fig = figure('Color','w', 'Units','pixels', 'Position',[40 40 figW figH], ...
    'Name','Class wheel (one diagram per class)', 'NumberTitle','off');
ax = axes('Parent', fig, 'Position', axPos);
hold(ax, 'on'); axis(ax, 'equal'); axis(ax, 'off');
lim = 1.02; xlim(ax, [-lim lim]); ylim(ax, [-lim lim]);
edgeGrey = [0.55 0.55 0.55];
opt = struct('markLP', markLP, 'markBP', markBP, 'labelMarkers', labelMarkers, ...
    'xLab', xLab, 'yLab', yLab, 'logSuffix', logSuffix, 'logZoom', logZoom, 'MS', MS, 'FS', FS, 'LW', LW);
for w = 1:nW
    c = classesShown(w); col = CI.colors(c,:);
    sector(ax, rC, rB, th2(w), th1(w), [1 1 1], edgeGrey, LW.sector);
    sector(ax, rB, rO, th2(w), th1(w), col, [1 1 1], LW.band);
    band_text(ax, bandLabels{c}, (rB + rO)/2, rO - rB, th2(w), th1(w), FS.band, [0.1 0.1 0.1]);
    ctr = rhoF(w)*[cosd(mids(w)), sind(mids(w))];
    diagram(ax, ctr, hFrame, rots(w), Rsel{w}, iconScale{c}, opt);
end
t = linspace(0, 2*pi, 400);
patch(ax, rC*cos(t), rC*sin(t), 'w', 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', LW.circle);
plot(ax, rO*cos(t), rO*sin(t), '-', 'Color', [0.2 0.2 0.2], 'LineWidth', LW.circle);
if ~isempty(centerText)
    text(ax, 0, 0, centerText, 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
        'FontSize', FS.centre, 'FontWeight','bold');
end
if ~isempty(panelLabel)
    text(ax, -lim, lim, panelLabel, 'FontSize', 18, 'FontWeight','bold', 'VerticalAlignment','top');
end
if showLegend
    % right edge of the wheel in figure units (square axes, axis equal)
    sq = min(axPos(3)*figW, axPos(4)*figH);                 % side of the square plot box (px)
    xRight = (axPos(1)*figW + axPos(3)*figW/2 + sq/2*rO/lim)/figW;
    hL = 0.80; wL = hL*figH*(5.2/(nW + 0.4))/figW;           % legend box matching its aspect
    axL = axes('Parent', fig, 'Position', [xRight + legendGap, 0.10, wL, hL]);
    hold(axL, 'on'); axis(axL, 'off'); axis(axL, 'equal');
    xlim(axL, [0 5.2]); ylim(axL, [0.3 nW + 0.7]);
    tt = linspace(0, 2*pi, 100); rr = 0.36;
    for w = 1:nW
        c = classesShown(w); y = nW + 1 - w;
        patch(axL, 0.45 + rr*cos(tt), y + rr*sin(tt), CI.colors(c,:), ...
            'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 2);
        text(axL, 0.45, y, sprintf('%d', c), 'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', 'FontSize', FS.legend + 2, 'FontWeight','bold');
        text(axL, 1.0, y, legendLabels{c}, 'VerticalAlignment','middle', ...
            'FontSize', FS.legend, 'FontWeight','bold');
    end
end

%% ---------------- 4. save -----------------------------------------------------------
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

%% ============== LOCAL FUNCTIONS =====================================================
function k = choose_set(RR, cand, pickMode, pLo, pHi, rs)
% one sample of the class: random, or 'clear' = random within a quantile window of a
% clarity score (spread of the singular points in x2 and ATOT, log10 units)
n = numel(cand);
if strcmpi(pickMode, 'clear')
    sc = NaN(n, 1);
    for i = 1:n
        A = [RR(i).exact_LP_ATOT(:); RR(i).exact_BP_ATOT(:)];
        X = [RR(i).exact_LP_x2(:);   RR(i).exact_BP_x2(:)];
        if numel(A) >= 2 && all(A > 0) && all(X > 0)
            sc(i) = log10(max(X)/min(X)) + log10(max(A)/min(A));
        end
    end
    ok = isfinite(sc);
    if nnz(ok) >= 5
        q = sort(sc(ok));
        lo = q(max(1, floor(pLo*numel(q)))); hi = q(max(1, ceil(pHi*numel(q))));
        win = find(ok & sc >= lo & sc <= hi);
        if ~isempty(win), k = cand(win(randperm(rs, numel(win), 1))); return; end
    end
end
k = cand(randperm(rs, n, 1));
end

function [h, rho] = frame_size(rIn, rOut, a1, a2, rot, m)
% largest square of half-size h, turned by rot degrees and centred on the bisector of
% the wedge rIn < r < rOut, a1 < angle < a2 (degrees), at radius rho; margin m.
% The whole outline (sampled) must lie inside the wedge.
mid = (a1 + a2)/2;
s = linspace(-1, 1, 25); o = ones(size(s));
U = [s, o, fliplr(s), -o]; V = [-o, s, o, fliplr(s)];          % square outline
Q = [cosd(rot) -sind(rot); sind(rot) cosd(rot)];
h = 0; rho = (rIn + rOut)/2;
for r0 = linspace(rIn, rOut, 120)
    c = r0*[cosd(mid); sind(mid)];
    lo = 0; hi = 0.6;
    for it = 1:30                                                % bisection on h
        hh = (lo + hi)/2;
        P = c + hh*Q*[U; V];
        rr = hypot(P(1,:), P(2,:)); aa = atan2d(P(2,:), P(1,:));
        da = mod(aa - mid + 180, 360) - 180;                     % angle from bisector
        marg = rad2deg(m./max(rr, eps));
        ok = all(rr >= rIn + m & rr <= rOut - m & abs(da) <= (a1 - a2)/2*sign(a1 - a2) - marg);
        if ok, lo = hh; else, hi = hh; end
    end
    if lo > h, h = lo; rho = r0; end
end
end

function sector(ax, r1, r2, a1, a2, faceCol, edgeCol, lw)
% annular sector between radii r1 < r2 and angles a1 < a2 (degrees)
t = linspace(deg2rad(a1), deg2rad(a2), max(8, ceil(abs(a2 - a1))));
x = [r2*cos(t), r1*cos(fliplr(t))];
y = [r2*sin(t), r1*sin(fliplr(t))];
patch(ax, x, y, faceCol, 'EdgeColor', edgeCol, 'LineWidth', lw);
end

function diagram(ax, ctr, h, rot, R, scaleMode, opt)
% one bifurcation diagram in a square frame of half-size h centred at ctr, turned by
% rot degrees. Local frame coordinates (u, v) in [-1, 1]; plot area [u0 u1] x [u0 u1].
u0 = -0.55; u1 = 0.95;
Q = [cosd(rot) -sind(rot); sind(rot) cosd(rot)];
W = @(u, v) deal(ctr(1) + h*(Q(1,1)*u + Q(1,2)*v), ctr(2) + h*(Q(2,1)*u + Q(2,2)*v));
% axes
[x, y] = W([u0 u0 u1], [u1 u0 u0]);
plot(ax, x, y, '-', 'Color', [0.1 0.1 0.1], 'LineWidth', opt.LW.axes);
isLog = strcmpi(scaleMode, 'log') && ~isempty([R.exact_LP_ATOT(:); R.exact_BP_ATOT(:)]);
sfx = ''; if isLog, sfx = opt.logSuffix; end
[x, y] = W((u0 + u1)/2, -0.80);
text(ax, x, y, [opt.xLab sfx], 'Rotation', rot, 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontSize', opt.FS.axis, 'FontWeight','bold');
[x, y] = W(-0.80, (u0 + u1)/2);
text(ax, x, y, [opt.yLab sfx], 'Rotation', rot + 90, 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontSize', opt.FS.axis, 'FontWeight','bold');
if isempty(R) || ~isfield(R, 'curves') || isempty(R.curves), return; end
[xl, yl] = icon_window(R, isLog, opt.logZoom);
if isLog, tr = @safelog10; else, tr = @(v) v; end
pu = @(a) u0 + (u1 - u0)*(tr(a) - xl(1))/(xl(2) - xl(1));
pv = @(v) u0 + (u1 - u0)*(tr(v) - yl(1))/(yl(2) - yl(1));
for cb = R.curves(:)'
    A = cb.ATOT(:)'; X = cb.x2(:)';
    if isempty(A), continue; end
    st = logical(cb.stable(:)');
    tA = tr(A); tX = tr(X);
    out = tA < xl(1) | tA > xl(2) | tX < yl(1) | tX > yl(2) | ~isfinite(tX) | ~isfinite(tA);
    ys = X; ys(~st | out) = NaN;                     % stable: solid black
    yu = X; yu(st) = NaN;                            % unstable: dashed red
    sw = find(diff(st) ~= 0);
    for j = sw, yu(j) = X(j); yu(j+1) = X(j+1); end
    yu(out) = NaN;
    [x, y] = W(pu(A), pv(ys)); plot(ax, x, y, '-',  'Color', [0 0 0],         'LineWidth', opt.LW.curve);
    [x, y] = W(pu(A), pv(yu)); plot(ax, x, y, '--', 'Color', [0.9 0.05 0.05], 'LineWidth', opt.LW.curve);
end
inBox = @(u, v) isfinite(u) & isfinite(v) & u >= u0 - 0.02 & u <= u1 + 0.02 & v >= u0 - 0.02 & v <= u1 + 0.02;
if opt.markBP && ~isempty(R.exact_BP_ATOT)
    uu = pu(R.exact_BP_ATOT(:)'); vv = pv(R.exact_BP_x2(:)');
    for j = find(inBox(uu, vv))
        [x, y] = W(uu(j), vv(j));
        plot(ax, x, y, 's', 'MarkerSize', opt.MS.BP, 'MarkerFaceColor', [0 0.62 0.2], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1);
        if opt.labelMarkers
            [x, y] = W(uu(j) + 0.10, vv(j) - 0.11);
            text(ax, x, y, 'BP', 'Rotation', rot, 'FontSize', opt.FS.marker, 'FontWeight','bold', ...
                'Color', [0 0.45 0.15], 'HorizontalAlignment','left', 'VerticalAlignment','middle');
        end
    end
end
if opt.markLP && ~isempty(R.exact_LP_ATOT)
    uu = pu(R.exact_LP_ATOT(:)'); vv = pv(R.exact_LP_x2(:)');
    done = zeros(0, 2);                               % labelled LPs (no crowding)
    for j = find(inBox(uu, vv))
        [x, y] = W(uu(j), vv(j));
        plot(ax, x, y, 'o', 'MarkerSize', opt.MS.LP, 'MarkerFaceColor', [0.9 0.05 0.05], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1);
        if opt.labelMarkers && (isempty(done) || min(hypot(done(:,1) - uu(j), done(:,2) - vv(j))) > 0.22)
            done(end+1, :) = [uu(j), vv(j)]; %#ok<AGROW>
            [x, y] = W(uu(j) - 0.08, vv(j) + 0.11);
            text(ax, x, y, 'LP', 'Rotation', rot, 'FontSize', opt.FS.marker, 'FontWeight','bold', ...
                'Color', [0.75 0 0], 'HorizontalAlignment','right', 'VerticalAlignment','middle');
        end
    end
end
end

function y = safelog10(v)
% log10 that maps v <= 0 to NaN (those points are not drawn)
v(v <= 0) = NaN; y = log10(v);
end

function [xl, yl] = icon_window(R, isLog, logZoom)
% plot window (log10 units when isLog), as in Ensemble_class_wheel_v3:
%  linear: zoomed on the folds / branch points; no singularity: automatic range
%  log   : ATOT from min(singular ATOT)/logZoom to max(singular ATOT)*logZoom
A = [R.exact_LP_ATOT(:); R.exact_BP_ATOT(:)]; X = [R.exact_LP_x2(:); R.exact_BP_x2(:)];
if isLog
    xl = log10([min(A)/logZoom, max(A)*logZoom]);
    v = [];
    for cb = R.curves(:)'
        in = cb.ATOT >= 10^xl(1) & cb.ATOT <= 10^xl(2) & cb.x2 > 0 & isfinite(cb.x2);
        v = [v, cb.x2(in)]; %#ok<AGROW>
    end
    if isempty(v), yl = [-3 3]; return; end
    yl = log10([min(v), max(v)]);
    pad = 0.06*max(diff(yl), 1); yl = yl + [-pad pad];
    return
end
xl = auto_xlim(R); yl = y_range(R.curves, xl);
if isempty(A), return; end
wA = max(A) - min(A); wX = max(X) - min(X);
if numel(A) >= 2 && wA > 0 && wA < 0.5*max(A)          % narrow loop: zoom on it
    xl = [max(0, min(A) - 1.5*wA), max(A) + 1.5*wA];
else
    xl = [0, 1.6*max(A)];
end
if numel(X) >= 2 && wX > 0
    yl = [max(0, min(X) - 1.2*wX), max(X) + 1.2*wX];
else
    yl = [0, 2.2*max(X)];
end
yAll = y_range(R.curves, xl);                            % never above the data
if yAll(2) < yl(2), yl(2) = yAll(2); end
if yl(2) <= yl(1), yl = y_range(R.curves, xl); end
end

function band_text(ax, str, rho, thick, a1, a2, fs, col)
% class name as ONE text object, straight, centred in the band at the middle of the
% wedge [a1, a2] (degrees) and tilted along the band (tangent to the circle there).
% Upright on the top half, turned by 180 degrees on the bottom half so that it is
% never upside down. The font shrinks (>= 6 pt) until the text fits the chord of the
% band and its thickness.
mid = (a1 + a2)/2;
rot = mid - 90; if sind(mid) < -0.05, rot = mid + 90; end
rot = mod(rot + 180, 360) - 180;
chord = 2*rho*sind(abs(a2 - a1)/2)*0.90;
h = text(ax, rho*cosd(mid), rho*sind(mid), str, 'Rotation', 0, 'FontSize', fs, ...
    'FontWeight','bold', 'Color', col, 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'Units','data');
for f = fs:-0.5:6
    set(h, 'FontSize', f);
    e = get(h, 'Extent');                 % unrotated extent in data units
    if e(3) <= chord && e(4) <= 0.80*thick, break; end
end
set(h, 'Rotation', rot);
end
