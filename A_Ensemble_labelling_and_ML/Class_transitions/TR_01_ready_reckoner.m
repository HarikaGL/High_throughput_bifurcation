% ======================= TR_01_ready_reckoner.m ========================================
% Fig 11 of the manuscript: "ready reckoner" of the ten bifurcation classes.
% For every class: its number and definition, the exact region it occupies, a one-line
% reading of what it means experimentally, the class (colour) against the whole ensemble
% (grey) in the plane of the two catalytic-dominance ratios
%     x = log10(k2 KTOT / k8 PTOT),   y = log10(k4 KTOT / k6 PTOT)
% (detached branch iff x < 0 and y > 0, pink), and the exact bifurcation diagram
% (x2 against ATOT) of one representative set.
% Open circles (Classes 1, 3, 9, 10): sets in the detached-branch quadrant whose detached
% fold lies above the labelling window ATOT <= 1e6, so the label counts only the main
% branch (tr_window_sets; full analysis in TR_05_window_open_circles.m). Star: example
% set (a) of TR_05, whose exact diagram extended beyond 1e6 is shown in panel (B).
%
% Representative set of each class (rule, reproducible): among the nCand sets of the class
% nearest to the class median in (x, y) [in (log10 k1/k5, log10 k11/k55) for Classes 7, 8],
% the set whose singularities all lie in 10^0.5 <= ATOT <= 10^5 and whose feature is the
% most visible (widest main-branch fold pair for Classes 1-4, with the detached fold of
% Classes 2 and 4 at least 0.7 decades from the main folds; for 7, 8 the largest
% |log10 ATOT(LP)/ATOT(BP)|); Classes 5, 6, 9, 10: the nearest such set to the median
% (Class 10 from Iterations 1-2, i.e. an asymmetric monotone set).
% Set repRows (iteration, row of ML_Ready_Dataset_v3_iter<i>.csv; row 1 = first data row)
% to fix the sets by hand.
%
% Output: ../Figures_combined_v3/Fig11_ready_reckoner/Fig11_ready_reckoner.png/.pdf/.fig and
%         Fig11_representative_sets.csv (the ten sets, their exact class and conditions)
% Needs ../support (exact labeller) and tr_lib/. Base MATLAB only; runs in GNU Octave too.
% ======================================================================================
close all;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
PARENT = fileparts(HERE);
addpath(fullfile(PARENT, 'support'), fullfile(HERE, 'tr_lib'));
outDir = fullfile(PARENT, 'Figures_combined_v3', 'Fig11_ready_reckoner');
if ~exist(outDir, 'dir'), mkdir(outDir); end
if ~exist('nCand', 'var'), nCand = 150; end
if ~exist('repRows', 'var'), repRows = []; end    % e.g. [1 9137; 1 4896; ...] (10 x 2)
nGrey = 9000; FS = 9;
resDir = fullfile(HERE, 'results'); if ~exist(resDir, 'dir'), mkdir(resDir); end

E = tr_load_ensemble(HERE);
use = E.use; Q = tr_conditions(E.p);
W = tr_window_sets(E, Q, resDir);          % open circles (cached in results/)
CI = class_info_v3();

%% ---------------- representative sets ---------------------------------------------------
rep = zeros(10, 1);
for c = 1:10
    if ~isempty(repRows)
        rep(c) = find(E.iteration == repRows(c,1) & E.row == repRows(c,2), 1); continue
    end
    k = find(use & E.class == c);
    if c == 10, k = k(E.iteration(k) <= 2); end   % a generic (asymmetric) monotone set
    if any(c == [7 8]), u = Q.a(k); v = Q.b(k); else, u = Q.x(k); v = Q.y(k); end
    d = (u - median(u)).^2 + (v - median(v)).^2;
    [~, o] = sort(d); k = k(o(1:min(nCand, numel(o))));
    best = -Inf; rep(c) = k(1);
    for j = k(:)'
        [cl, S] = tr_label(E.p(j, :));
        if cl ~= c, continue; end
        sing = [[S.folds.ATOT], [S.bp.ATOT]];
        if any(sing < 10^0.5 | sing > 1e5), continue; end
        isM = arrayfun(@(f) any(strcmp(S.branch(f.branch).type, {'main','sym'})), S.folds);
        switch c
            case {1, 3}
                fa = [S.folds(isM).ATOT]; score = log10(max(fa)/min(fa));
            case {2, 4}           % wide main fold pair, detached fold well separated from it
                fa = [S.folds(isM).ATOT]; fd = [S.folds(~isM).ATOT];
                if isempty(fd) || min(abs(log10(fd) - log10(fa(:)'))) < 0.7, continue; end
                score = log10(max(fa)/min(fa));
            case {7, 8}
                score = abs(log10(S.folds(1).ATOT/S.bp(1).ATOT));
            otherwise
                score = 0;
        end
        if score > best + 1e-12, best = score; rep(c) = j; end
        if any(c == [5 6 9 10]), break; end
    end
end
R = [(1:10)', E.iteration(rep), E.row(rep), E.class(rep), Q.x(rep), Q.y(rep), Q.side(rep), ...
     Q.bp(rep), Q.km(rep), E.p(rep, :)];
tr_write_csv(fullfile(outDir, 'Fig11_representative_sets.csv'), ['class,iteration,row,ensemble_class,' ...
    'x,y,side,log10_k4_over_k2,log10_KM1_over_KM5,KTOT,PTOT,k1,k11,k2,k3,k33,k4,k5,k55,k6,k7,k77,k8'], R);

%% ---------------- text of each class ----------------------------------------------------
TXT = {
 'LP only, S-shaped',          'region: x > 0, or detached fold > 10^6', {'Bistable between two folds on', 'the main branch. A dose ramp', 'shows a hysteresis loop.'}
 'S-shaped LP + detached LP',  'region: x < 0, y > 0; u > r',         {'S-hysteresis on the main branch', 'plus a detached branch. A ramp', 'sees only the main-branch loop.'}
 'LP only, Z-shaped',          'region: main branch on u < r',          {'Bistable between two folds;', 'output falls with A_{TOT}. A dose', 'ramp shows a hysteresis loop.'}
 'Z-shaped LP + detached LP',  'region: x < 0, y > 0; u < r',         {'Z-hysteresis on the main branch', 'plus a detached branch. A ramp', 'sees only the main-branch loop.'}
 'LP on detached branch only', 'region: x < 0, y > 0',                {'Hidden switch: graded main', 'branch; the second stable state', 'is on the detached branch only.'}
 'BP only (pitchfork)',        'region: u* = r, k_4 > k_2, K_{M1} = K_{M5}', {'Two stable arms beyond the', 'branch point. Needs exact symmetry;', 'any asymmetry removes it.'}
 'LP + BP, x_2(LP) > x_2(BP)', 'region: u* = r, k_4 > k_2, K_{M5} < K_{M1}', {'Crossing branch with a fold above', 'the branch point. Bistable;', 'fragile to asymmetry.'}
 'LP + BP, x_2(LP) < x_2(BP)', 'region: u* = r, k_4 > k_2, K_{M5} > K_{M1}', {'Crossing branch with a fold below', 'the branch point. Bistable;', 'fragile to asymmetry.'}
 'Biphasic, no singularity',   'region: x < 0, u < r; no fold',       {'Monostable; output rises,', 'peaks and falls. No hysteresis.', ''}
 'Monotone, no singularity',   'region: x > 0 or u > r; no fold',     {'Monostable, graded or', 'ultrasensitive. No hysteresis.', ''}};

%% ---------------- figure ----------------------------------------------------------------
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [1 1 42 44], 'PaperPositionMode', 'auto');
gi = find(use); gi = gi(tr_pick(numel(gi), min(nGrey, numel(gi)), 11));
rowH = 0.124; top = 0.905; star = zeros(1, 10);
for c = 1:10
    col = floor((c-1)/5); r = mod(c-1, 5);
    x0 = 0.015 + 0.5*col; y0 = top - (r+1)*rowH + 0.022;
    % badge and text
    axT = axes('Position', [x0, y0, 0.17, rowH - 0.03]); axis(axT, 'off'); hold(axT, 'on');
    set(axT, 'XLim', [0 1], 'YLim', [0 1]);
    plot(axT, 0.1, 0.515, 'o', 'MarkerSize', 34, 'MarkerFaceColor', tr_colors(c), 'MarkerEdgeColor', 'none');
    text(axT, 0.1, 0.515, sprintf('%d', c), 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'w');
    text(axT, 0.25, 0.92, TXT{c,1}, 'FontSize', FS+1, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    text(axT, 0.25, 0.72, TXT{c,2}, 'FontSize', FS-1, 'FontAngle', 'italic', 'VerticalAlignment', 'top');
    for q = 1:3
        text(axT, 0.25, 0.48 - 0.13*(q-1), TXT{c,3}{q}, 'FontSize', FS-1, 'VerticalAlignment', 'top');
    end
    % scatter
    axS = axes('Position', [x0 + 0.185, y0 + 0.012, 0.115, rowH - 0.04]); hold(axS, 'on'); box(axS, 'on');
    patch(axS, [-6 0 0 -6], [0 0 7 7], [0.98 0.87 0.85], 'EdgeColor', 'none');
    plot(axS, Q.x(gi), Q.y(gi), '.', 'Color', [0.80 0.80 0.80], 'MarkerSize', 3);
    k = find(use & E.class == c); k = k(tr_pick(numel(k), min(3000, numel(k)), 100 + c));
    plot(axS, Q.x(k), Q.y(k), '.', 'Color', tr_colors(c), 'MarkerSize', 4);
    plot(axS, [-6 6], [0 0], 'k:', 'LineWidth', 0.8); plot(axS, [0 0], [-6 7], 'k:', 'LineWidth', 0.8);
    set(axS, 'XLim', [-5.5 4.5], 'YLim', [-4.5 6], 'FontSize', FS-2, 'Layer', 'top');
    xlabel(axS, 'x = log_{10}(k_2K_{TOT}/k_8P_{TOT})', 'FontSize', FS-2);
    ylabel(axS, 'y = log_{10}(k_4K_{TOT}/k_6P_{TOT})', 'FontSize', FS-2);
    text(axS, 4.3, -4.2, sprintf('n = %s', tr_commas(sum(use & E.class == c))), 'FontSize', FS-2, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    if any(c == [1 3 9 10])                 % open circles: detached fold beyond the window
        kw = find(W.class == c);
        plot(axS, W.x(kw), W.y(kw), 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 0.6);
        ex = tr_window_examples(W, c); star(c) = ex(1);
        plot(axS, W.x(ex(1)), W.y(ex(1)), 'p', 'MarkerSize', 13, 'MarkerFaceColor', [0.98 0.80 0.20], 'MarkerEdgeColor', 'k');
        text(axS, -5.3, 5.7, sprintf('open: %.1f%% of class', 100*numel(kw)/W.nClass(c)), 'FontSize', FS-2, ...
             'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'Margin', 0.5);
    end
    % diagram
    axD = axes('Position', [x0 + 0.345, y0 + 0.012, 0.135, rowH - 0.04]);
    D = tr_curves(E.p(rep(c), :), 1e-2, 1e6, 600);
    opt = struct('lw', 1.5, 'ms', 5);
    if any(c == [7 8]) && ~isempty(D.lp) && abs(log10(D.lp(1,1)/D.bp(1,1))) < 0.3
        opt.zoom = tr_bp_zoom(E.p(rep(c), :), 0.3);   % inset only when the fold is close to the BP
    end
    tr_draw_diagram(axD, D, opt);
    set(axD, 'FontSize', FS-2);
    tr_tick_decades(axD);
    xlabel(axD, 'A_{TOT}', 'FontSize', FS-2); ylabel(axD, 'x_2 = [APP]', 'FontSize', FS-2);
end
% title and legend
axH = axes('Position', [0 0.925 1 0.075]); axis(axH, 'off'); set(axH, 'XLim', [0 1], 'YLim', [0 1]);
text(axH, 0.005, 0.8, '(A)', 'FontSize', 16, 'FontWeight', 'bold');
text(axH, 0.5, 0.78, 'The ten bifurcation classes of the dual phosphorylation cycle (exact diagrams of x_2 = [APP] against A_{TOT})', ...
     'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
text(axH, 0.5, 0.45, ['Middle: the class (colour) against the whole ensemble (grey; 39,968 sets) in the plane of the two catalytic-dominance ratios; ' ...
     'pink: detached-branch region x < 0, y > 0.  Right: exact diagram of one representative set.'], ...
     'HorizontalAlignment', 'center', 'FontSize', FS);
text(axH, 0.5, 0.18, ['Open circles (Classes 1, 3, 9, 10): sets inside the pink quadrant whose detached branch exists but whose fold lies above A_{TOT} = 10^6; ' ...
     'the label counts only the main branch. Star: the set shown in (B).'], 'HorizontalAlignment', 'center', 'FontSize', FS);

%% ---------------- (B) starred sets: diagrams extended beyond the window ---------------------
axB = axes('Position', [0 0.262 1 0.03]); axis(axB, 'off'); set(axB, 'XLim', [0 1], 'YLim', [0 1]);
text(axB, 0.005, 0.5, '(B)', 'FontSize', 16, 'FontWeight', 'bold');
text(axB, 0.5, 0.5, ['Starred sets of (A): the detached branch lies beyond the labelling window A_{TOT} \leq 10^6 (grey), ' ...
     'outside the experimentally accessible range'], 'HorizontalAlignment', 'center', 'FontSize', FS + 2, 'FontWeight', 'bold');
cs = [1 3 9 10];
for q = 1:4
    c = cs(q); j = W.idx(star(c));
    Amax = 10^(ceil(log10(W.detFold(star(c)))) + 1.5);
    D = tr_curves_ext(E.p(j, :), 1e-2, Amax, 700);
    a = axes('Position', [0.04 + 0.245*(q-1), 0.075, 0.195, 0.155]);
    tr_draw_diagram(a, D, struct('xlim', [1e-2 Amax], 'lw', 1.5, 'ms', 5));
    tr_shade_window(a, Amax); tr_tick_decades(a); set(a, 'FontSize', FS-2);
    xlabel(a, 'A_{TOT}', 'FontSize', FS-2); ylabel(a, 'x_2 = [APP]', 'FontSize', FS-2);
    title(a, {sprintf('Class %d (star): detached fold at A_{TOT} = 10^{%.1f}', c, log10(D.detFold)), ...
              sprintf('with the fold inside the window: Class %d', D.classWide)}, 'FontSize', FS, 'Color', tr_colors(c));
end
axL = axes('Position', [0 0.002 1 0.035]); axis(axL, 'off'); hold(axL, 'on'); set(axL, 'XLim', [0 1], 'YLim', [0 1]);
items = {'main branch', tr_colors('main'), '-'; 'detached branch', tr_colors('detached'), '-'; ...
         'branch-point line u = r', tr_colors('sym'), '-'; 'crossing branch', tr_colors('offsym'), '-'; ...
         'stable', [0 0 0], '-'; 'unstable', [0 0 0], '--'};
xpos = [0.03 0.12 0.22 0.35 0.45 0.51];
for q = 1:size(items, 1)
    plot(axL, xpos(q) + [0 0.025], [0.65 0.65], items{q,3}, 'Color', items{q,2}, 'LineWidth', 1.8);
    text(axL, xpos(q) + 0.03, 0.65, items{q,1}, 'FontSize', FS);
end
plot(axL, 0.585, 0.65, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.96 0.82 0.25], 'MarkerEdgeColor', 'k');
text(axL, 0.595, 0.65, 'fold (LP)', 'FontSize', FS);
plot(axL, 0.65, 0.65, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
text(axL, 0.66, 0.65, 'branch point (BP)', 'FontSize', FS);
plot(axL, 0.75, 0.65, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
text(axL, 0.76, 0.65, 'detached fold > 10^6', 'FontSize', FS);
patch(axL, [0.865 0.885 0.885 0.865], [0.45 0.45 0.85 0.85], [0.9 0.9 0.9], 'EdgeColor', [0.6 0.6 0.6]);
text(axL, 0.89, 0.65, 'A_{TOT} > 10^6 (outside window)', 'FontSize', FS);
text(axL, 0.5, 0.15, ['x = log_{10}(k_2K_{TOT}/k_8P_{TOT}), y = log_{10}(k_4K_{TOT}/k_6P_{TOT}); u = K/P, r = K_{TOT}/P_{TOT}, ' ...
     'u* = -\beta/\alpha (Eq. 7). Classes refer to A_{TOT} \leq 10^6.'], 'HorizontalAlignment', 'center', 'FontSize', FS-1);
tr_save(fig, fullfile(outDir, 'Fig11_ready_reckoner'));
fprintf('Fig 11 saved in %s\nRepresentative sets (iteration, row): %s\n', outDir, mat2str(R(:, 2:3)));
