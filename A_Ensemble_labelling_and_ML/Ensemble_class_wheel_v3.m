% ======================= Ensemble_class_wheel_v3.m ==================================
% "Class wheel": one circle for the whole ensemble (all four iterations, 40,000
% parameter sets). Each class gets a wedge whose angle is its share of the ensemble
% (label: name and percentage on the coloured outer band). Inside each wedge, small
% discs show REAL bifurcation diagrams (x2 vs ATOT) of randomly chosen parameter sets
% of that class: solid black = stable, dashed red = unstable steady states.
%
% Run from the Ensemble_v3 folder (paths are set automatically). Reads
%   Ensemble_iter<i>_v3_1e6/MAT_Files/Complete_Bifurcation_Results_v3.mat, i = 1..4
% Saves Figures_combined_v3/Class_Wheel_v3.png/.pdf/.fig
% Memory: each results file is loaded once (about 3-5 GB in memory) and only the
% chosen parameter sets are kept.
% ===================================================================================
clear; close all; clc;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
addpath(fullfile(HERE, 'support')); cd(HERE);

%% ---------------- USER SETTINGS ----------------------------------------------------
resultDirs   = {'Ensemble_iter1_v3_1e6','Ensemble_iter2_v3_1e6', ...
                'Ensemble_iter3_v3_1e6','Ensemble_iter4_v3_1e6'};
rngSeed      = 7;          % which parameter sets are drawn (change for another choice)
angleMode    = 'min';      % 'proportional' : wedge angle = exact share of the ensemble
                           % 'min'          : same, but every class gets at least
                           %                  minAngleDeg (small classes stay readable);
                           %                  the printed percentages are always exact
minAngleDeg  = 16;
startAngleDeg= 90;         % where class 1 starts (90 = top); classes run clockwise
iconRadius   = 0.058;      % radius of the small discs (wheel radius = 1)
maxIcons     = 24;         % at most this many discs per class
excludeLPeqBP= true;       % do not draw the flagged LP = BP sets (they are still counted)
showMarkers  = false;      % red dot / green square at LP / BP inside the discs
showPercent  = false;      % true: add the percentage after the class name on the band
% line widths (points)
LW.curve   = 3.0;          % bifurcation curves inside the discs
LW.axes    = 2.0;          % small axes inside the discs
LW.disc    = 2.5;          % outline of the discs
LW.sector  = 2.5;          % grey lines between the wedges
LW.band    = 3.0;          % white lines between the coloured band segments
LW.circle  = 3.5;          % centre circle and outer circle
centerText   = {'Classification','and','Prediction'};   % {} for an empty centre
panelLabel   = '';         % e.g. '[IV]'
labelFont    = 11;         % font size of the class names on the band
letterSpace  = 1.08;       % spacing of the letters on the band (1 = tight)
% axes of the small diagrams, per class (order of class_info_v3):
%   'lin' : linear axes, zoomed on the folds (single-scale shapes: S, Z, biphasic, ...)
%   'log' : log-log axes, ATOT from min(singular ATOT)/logZoom to max*logZoom; needed
%           where branches live on very different scales (detached branch, pitchfork,
%           LP + BP)
iconScale = {'lin','log','lin','log','log','log','log','log','lin','lin','lin'};
logZoom   = 30;
layoutIter = 400;          % iterations of the even (irregular) spreading of the discs
fillFrac   = 0.62;         % fraction of a wedge's free area covered by discs
outName      = 'Class_Wheel_v3';

% short names on the band (same order as class_info_v3)
bandLabels = {'LP ONLY (S shaped)','LP (S shaped) + DETACHED','LP ONLY (Z shaped)', ...
              'LP (Z shaped) + DETACHED','LP DETACHED ONLY','BP ONLY', ...
              'LP + BP (LP > BP)','LP + BP (LP < BP)','BIPHASIC','MONOTONE','OTHER'};

%% ---------------- GEOMETRY ---------------------------------------------------------
rC = 0.30;                 % centre circle
rB = 0.86;                 % inner edge of the coloured band
rO = 1.00;                 % outer edge
CI = class_info_v3(); nC = CI.n;

%% ---------------- 1. class counts of all iterations ---------------------------------
cls = []; flg = []; src = [];
for d = 1:numel(resultDirs)
    f = fullfile(HERE, resultDirs{d}, 'MAT_Files', 'Complete_Bifurcation_Results_v3.mat');
    if ~exist(f, 'file'), error('Missing %s', f); end
    S = load(f, 'exact_class_v3', 'lp_eq_bp');
    cls = [cls; S.exact_class_v3(:)]; flg = [flg; S.lp_eq_bp(:)]; %#ok<AGROW>
    src = [src; [d*ones(numel(S.exact_class_v3),1), (1:numel(S.exact_class_v3))']]; %#ok<AGROW>
end
N = numel(cls);
counts = histc(cls(cls > 0), 1:nC); counts = counts(:)';
active = find(counts > 0);
pct = 100*counts/N;
fprintf('Ensemble: %d parameter sets\n', N);
for c = active, fprintf('  %2d %-20s %6d  %5.1f%%\n', c, CI.names{c}, counts(c), pct(c)); end

%% ---------------- 2. wedge angles ----------------------------------------------------
share = counts(active)/sum(counts(active));
ang = 360*share;
if strcmpi(angleMode, 'min')
    for it = 1:50                                   % lift small wedges, shrink the rest
        small = ang < minAngleDeg;
        if ~any(small), break; end
        ang(small) = minAngleDeg;
        free = 360 - sum(ang(small));
        ang(~small) = free*share(~small)/sum(share(~small));
    end
end
edges = startAngleDeg - [0 cumsum(ang)];            % clockwise
th1 = edges(1:end-1); th2 = edges(2:end);           % th1 > th2 (degrees)

%% ---------------- 3. disc positions in every wedge ----------------------------------
% evenly but irregularly spread: random start, then mutual repulsion inside the wedge
% (like particles that fill the space), so there are no rows or rings
rsL = RandStream('mt19937ar', 'Seed', rngSeed + 100);
pos = cell(1, numel(active));
for w = 1:numel(active)
    pos{w} = layout_discs(rC, rB, th2(w), th1(w), iconRadius, maxIcons, fillFrac, layoutIter, rsL);
end

%% ---------------- 4. random parameter sets for the discs ----------------------------
rs = RandStream('mt19937ar', 'Seed', rngSeed);
pick = cell(1, numel(active));
for w = 1:numel(active)
    c = active(w);
    cand = find(cls == c & ~(excludeLPeqBP & flg));
    k = min(size(pos{w},1), numel(cand));
    pick{w} = cand(randperm(rs, numel(cand), k));
end
allPick = vertcat(pick{:});
Rsel = cell(N, 1);                                   % only the chosen sets are kept
for d = 1:numel(resultDirs)
    need = allPick(src(allPick,1) == d);
    if isempty(need), continue; end
    f = fullfile(HERE, resultDirs{d}, 'MAT_Files', 'Complete_Bifurcation_Results_v3.mat');
    fprintf('Loading %s (%d sets) ...\n', resultDirs{d}, numel(need));
    S = load(f, 'res');
    for g = need(:)', Rsel{g} = S.res(src(g,2)); end
    clear S
end

%% ---------------- 5. draw ------------------------------------------------------------
fig = figure('Color','w', 'Units','pixels', 'Position',[60 40 1100 1100], ...
    'Name','Class wheel', 'NumberTitle','off');
ax = axes('Parent', fig, 'Position', [0.02 0.02 0.96 0.96]);
hold(ax, 'on'); axis(ax, 'equal'); axis(ax, 'off');
lim = 1.08; xlim(ax, [-lim lim]); ylim(ax, [-lim lim]);
edgeGrey = [0.55 0.55 0.55];

for w = 1:numel(active)
    c = active(w); col = CI.colors(c,:);
    % white inner sector with grey border, coloured outer band
    sector(ax, rC, rB, th2(w), th1(w), [1 1 1], edgeGrey, LW.sector);
    sector(ax, rB, rO, th2(w), th1(w), col, [1 1 1], LW.band);
    % class name and percentage along the band
    lab = bandLabels{c};
    if showPercent, lab = sprintf('%s  %.1f%%', lab, pct(c)); end
    arc_text(ax, lab, (rB + rO)/2, th2(w), th1(w), labelFont, [0.1 0.1 0.1], letterSpace);
    % discs with real diagrams
    tint = 1 - 0.45*(1 - col);                       % lighter version of the class colour
    for i = 1:numel(pick{w})
        R = Rsel{pick{w}(i)};
        icon(ax, pos{w}(i,:), iconRadius, R, tint, showMarkers, iconScale{c}, logZoom, LW);
    end
end
% centre and outline
t = linspace(0, 2*pi, 400);
patch(ax, rC*cos(t), rC*sin(t), 'w', 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', LW.circle);
plot(ax, rO*cos(t), rO*sin(t), '-', 'Color', [0.2 0.2 0.2], 'LineWidth', LW.circle);
if ~isempty(centerText)
    text(ax, 0, 0, centerText, 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
        'FontSize', 17, 'FontWeight','bold');
end
if ~isempty(panelLabel)
    text(ax, -lim, lim, panelLabel, 'FontSize', 18, 'FontWeight','bold', 'VerticalAlignment','top');
end
text(ax, lim, -lim, sprintf('n = %d parameter sets (4 iterations)', N), ...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize', 9, 'Color', [0.35 0.35 0.35]);

%% ---------------- 6. save -----------------------------------------------------------
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
function sector(ax, r1, r2, a1, a2, faceCol, edgeCol, lw)
% annular sector between radii r1 < r2 and angles a1 < a2 (degrees)
t = linspace(deg2rad(a1), deg2rad(a2), max(8, ceil(abs(a2 - a1))));
x = [r2*cos(t), r1*cos(fliplr(t))];
y = [r2*sin(t), r1*sin(fliplr(t))];
patch(ax, x, y, faceCol, 'EdgeColor', edgeCol, 'LineWidth', lw);
end

function icon(ax, c, r, R, tint, showMarkers, scaleMode, logZoom, LW)
% one disc with the bifurcation diagram of R (x2 vs ATOT); linear or log-log axes
t = linspace(0, 2*pi, 80);
patch(ax, c(1) + r*cos(t), c(2) + r*sin(t), tint, 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', LW.disc);
s = 0.60*r;                                        % half-width of the plot box
x0 = c(1) - s; y0 = c(2) - s;
plot(ax, [x0 x0 x0+2*s], [y0+2*s y0 y0], '-', 'Color', [0.15 0.15 0.15], 'LineWidth', LW.axes);   % axes
if isempty(R) || ~isfield(R, 'curves') || isempty(R.curves), return; end
isLog = strcmpi(scaleMode, 'log') && ~isempty([R.exact_LP_ATOT(:); R.exact_BP_ATOT(:)]);
[xl, yl] = icon_window(R, isLog, logZoom);
if isLog, tr = @safelog10; else, tr = @(v) v; end
mx = @(a) x0 + 2*s*(tr(a) - xl(1))/(xl(2) - xl(1));
my = @(v) y0 + 2*s*(tr(v) - yl(1))/(yl(2) - yl(1));
for cb = R.curves(:)'
    A = cb.ATOT(:)'; X = cb.x2(:)';
    if isempty(A), continue; end
    st = logical(cb.stable(:)');
    tA = tr(A); tX = tr(X);
    out = tA < xl(1) | tA > xl(2) | tX < yl(1) | tX > yl(2) | ~isfinite(tX) | (isLog & (A <= 0 | X <= 0));
    ys = X; ys(~st | out) = NaN;                    % stable: solid black
    yu = X; yu(st) = NaN;                           % unstable: dashed red
    sw = find(diff(st) ~= 0);
    for j = sw, yu(j) = X(j); yu(j+1) = X(j+1); end
    yu(out) = NaN;
    plot(ax, mx(A), my(ys), '-',  'Color', [0 0 0],        'LineWidth', LW.curve);
    plot(ax, mx(A), my(yu), '--', 'Color', [0.9 0.05 0.05], 'LineWidth', LW.curve);
end
if showMarkers
    if ~isempty(R.exact_LP_ATOT)
        plot(ax, mx(R.exact_LP_ATOT), my(R.exact_LP_x2), '.', 'Color', [0.9 0 0], 'MarkerSize', 7);
    end
    if ~isempty(R.exact_BP_ATOT)
        plot(ax, mx(R.exact_BP_ATOT), my(R.exact_BP_x2), 's', 'MarkerFaceColor', [0 0.6 0], ...
            'MarkerEdgeColor', 'none', 'MarkerSize', 3);
    end
end
end

function y = safelog10(v)
% log10 that keeps NaN as NaN and maps v <= 0 to NaN (so those points are not drawn)
v(v <= 0) = NaN; y = log10(v);
end

function [xl, yl] = icon_window(R, isLog, logZoom)
% plot window of one disc (in log10 units when isLog).
%  linear: zoomed on the folds / branch points (narrow hysteresis loops visible);
%          samples without singularity keep the automatic range (auto_xlim)
%  log   : ATOT from min(singular ATOT)/logZoom to max(singular ATOT)*logZoom, x2 over
%          the data inside that range
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
yAll = y_range(R.curves, xl);                           % never above the data
if yAll(2) < yl(2), yl(2) = yAll(2); end
if yl(2) <= yl(1), yl = y_range(R.curves, xl); end
end

function P = layout_discs(rIn, rOut, a1, a2, r, nMax, fillFrac, nIter, rs)
% Evenly but irregularly spread discs of radius r in the annular wedge
% rIn < rho < rOut, a1 < angle < a2 (degrees). Number of discs from the free area,
% random start (dart throwing), then mutual repulsion + projection into the wedge.
gap = 0.25*r; D = 2*r + gap;                       % minimum centre distance
t1 = deg2rad(a1); t2 = deg2rad(a2);
rLo = rIn + gap + r; rHi = rOut - gap - r;
if rHi < rLo, P = zeros(0,2); return; end
area = 0.5*(t2 - t1)*(rOut^2 - rIn^2);
n = max(1, min(nMax, floor(fillFrac*area/(pi*(r + gap/2)^2))));
inside = @(p) wedge_ok(p, rLo, rHi, t1, t2, r + gap/2);
for tryN = n:-1:1
    % dart throwing
    P = zeros(0,2); att = 0;
    while size(P,1) < tryN && att < 20000
        att = att + 1;
        rho = sqrt(rLo^2 + rand(rs)*(rHi^2 - rLo^2)); th = t1 + rand(rs)*(t2 - t1);
        q = [rho*cos(th), rho*sin(th)];
        if ~inside(q), continue; end
        if isempty(P) || min(sum((P - q).^2, 2)) >= (0.9*D)^2, P(end+1,:) = q; end %#ok<AGROW>
    end
    if size(P,1) < tryN, continue; end
    % repulsion towards an even spacing larger than D where room allows
    Dt = max(D, 0.98*sqrt(area/tryN));
    for it = 1:nIter
        F = zeros(size(P));
        for i = 1:size(P,1)
            d = P(i,:) - P; dist = sqrt(sum(d.^2, 2)); dist(i) = Inf;
            k = dist < Dt;
            if any(k), F(i,:) = sum(d(k,:)./dist(k).*(Dt - dist(k)), 1); end
        end
        step = 0.25*(1 - it/nIter) + 0.02;
        P = P + step*F + 0.004*r*(rand(rs, size(P)) - 0.5);   % small thermal jitter
        P = project_wedge(P, rLo, rHi, t1, t2, r + gap/2);
    end
    dmin = Inf;
    for i = 1:size(P,1)
        d = sqrt(sum((P - P(i,:)).^2, 2)); d(i) = Inf; dmin = min(dmin, min(d));
    end
    if dmin >= 2*r + 0.5*gap, return; end            % no overlaps: done
end
end

function ok = wedge_ok(p, rLo, rHi, t1, t2, m)
rho = hypot(p(1), p(2)); th = atan2(p(2), p(1));
th = t1 + mod(th - t1, 2*pi);
ok = rho >= rLo && rho <= rHi && th >= t1 + m/rho && th <= t2 - m/rho;
end

function P = project_wedge(P, rLo, rHi, t1, t2, m)
rho = hypot(P(:,1), P(:,2)); th = atan2(P(:,2), P(:,1));
th = t1 + mod(th - t1 + pi/2, 2*pi) - pi/2;          % angle near the wedge
rho = min(max(rho, rLo), rHi);
th = min(max(th, t1 + m./rho), t2 - m./rho);
P = [rho.*cos(th), rho.*sin(th)];
end

function arc_text(ax, str, rho, a1, a2, fs, col, letterSpace)
% text written along an arc, centred in the wedge [a1, a2] (degrees); readable from
% outside on the top half and from inside on the bottom half. The font shrinks to
% fit the wedge (down to 6 pt); if it still does not fit, the label goes outside the
% wheel, radially.
mid = (a1 + a2)/2;
bottom = sind(mid) < -0.05;
span = deg2rad(a2 - a1)*rho*0.92;                   % available arc length
for f = fs:-0.5:6
    [w, total] = char_widths(ax, str, f);
    w = w*letterSpace; total = total*letterSpace;
    if total <= span, break; end
end
if total > span                                      % does not fit: outside, radial
    r = 1.03; rot = mid; ha = 'left';
    if cosd(mid) < 0, rot = mid + 180; ha = 'right'; end
    text(ax, r*cosd(mid), r*sind(mid), str, 'Rotation', rot, 'HorizontalAlignment', ha, ...
        'VerticalAlignment','middle', 'FontSize', fs*0.85, 'FontWeight','bold', 'Color', col);
    return
end
cum = [0 cumsum(w)];
centres = (cum(1:end-1) + cum(2:end))/2 - total/2;  % arc position of each character
for k = 1:numel(str)
    if bottom
        a = mid + rad2deg(centres(k)/rho);          % left to right = increasing angle
        rot = a + 90;
    else
        a = mid - rad2deg(centres(k)/rho);          % left to right = decreasing angle
        rot = a - 90;
    end
    text(ax, rho*cosd(a), rho*sind(a), str(k), 'Rotation', rot, 'FontSize', f, ...
        'FontWeight','bold', 'Color', col, 'HorizontalAlignment','center', ...
        'VerticalAlignment','middle');
end
end

function [w, total] = char_widths(ax, str, fs)
% width (data units) of every character of str at font size fs (bold)
h = text(ax, 0, 0, '|', 'FontSize', fs, 'FontWeight','bold', 'Units','data', 'Visible','off');
e = get(h, 'Extent'); wBar = e(3);
w = zeros(1, numel(str)); prev = 0;
for k = 1:numel(str)
    set(h, 'String', [str(1:k) '|']);
    e = get(h, 'Extent'); cur = e(3) - wBar;
    w(k) = cur - prev; prev = cur;
end
total = prev;
delete(h);
end
