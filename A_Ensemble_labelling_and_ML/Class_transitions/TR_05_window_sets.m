% ======================= TR_05_window_sets.m ===========================================
% Detached branches beyond the labelling window (Supplementary Note 18).
% PART 1 (Supplementary Fig 23; the open circles and stars of Fig 11): sets in the
%   detached-branch quadrant (x < 0, y > 0, criterion (10)) labelled 1, 3, 9 or 10. Their
%   detached branch exists, but its fold lies above ATOT = 1e6, so the label counts only
%   the main branch. For each class: the (x, y) panel with the open circles and two starred
%   examples (a: detached fold nearest 10^7, the set of Fig 11B; b: deepest inside the
%   quadrant), and the exact diagrams of a and b extended beyond 1e6.
% PART 2 (Supplementary Fig 24, Supplementary Table 19): misclassified held-out test sets
%   of the pooled LightGBM model (ML_v3_results_parallel/labels_Exact_Class_v3/
%   misclassified/Misclassified_test.csv) that are such sets: the classifier recognises
%   the detached branch from the signs of x and y, the label does not count it.
% Adapted from make_fig_window_open_circles.m and make_fig_window_misclassified.m
% (K. Sriram) to the ensemble and labeller of this study.
% Output: ../Figures_combined_v3/Class_transitions/S_window_open_circles.*,
%         S_window_misclassified.*; results/TR_window_sets.csv, TR_window_summary.txt,
%         TR_window_confusions.csv, TR_window_misclassified.txt
% ======================================================================================
close all;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
PARENT = fileparts(HERE);
addpath(fullfile(fileparts(HERE), 'support'), fullfile(HERE, 'tr_lib'));
resDir = fullfile(HERE, 'results'); if ~exist(resDir, 'dir'), mkdir(resDir); end
figDir = fullfile(fileparts(HERE), 'Figures_combined_v3', 'Class_transitions');
if ~exist(figDir, 'dir'), mkdir(figDir); end
FS = 9;
E = tr_load_ensemble(HERE); Q = tr_conditions(E.p);
W = tr_window_sets(E, Q, resDir);
CI = class_info_v3();
names = {1, 'LP only, S-shaped'; 3, 'LP only, Z-shaped'; 9, 'biphasic'; 10, 'monotone'};

% summary
fid = fopen(fullfile(resDir, 'TR_window_summary.txt'), 'w');
for f = [1 fid]
    fprintf(f, 'Sets in the detached-branch quadrant (x < 0, y > 0) labelled 1, 3, 9 or 10 (ML sets)\n');
    for q = 1:4
        c = names{q,1}; k = W.class == c;
        fprintf(f, '  Class %2d: %4d of %5d (%.1f%%); detached fold above 1e6 in %d; log10 ATOT of the fold %.2f-%.2f (median %.2f)\n', ...
            c, sum(k), W.nClass(c), 100*sum(k)/W.nClass(c), sum(W.detFold(k) > 1e6), min(log10(W.detFold(k))), max(log10(W.detFold(k))), median(log10(W.detFold(k))));
    end
    fprintf(f, '  All: %d of %d ML sets (%.2f%%); detached fold above 1e6 in %d\n', numel(W.idx), W.nUse, 100*numel(W.idx)/W.nUse, sum(W.detFold > 1e6));
    fprintf(f, '  Distance to the quadrant edges, min(-x, y): median %.2f decades (Classes 2, 4, 5: %.2f)\n', median(W.margin), median(W.marginDet));
end
fclose(fid);

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [1 1 40 56]);
gi = find(E.use); gi = gi(tr_pick(numel(gi), min(6000, numel(gi)), 11));
for q = 1:4
    c = names{q,1}; k = find(W.class == c); ex = tr_window_examples(W, c);
    y0 = 0.735 - 0.225*(q-1); hA = 0.14;
    ax = axes('Position', [0.05, y0, 0.25, hA]);
    kc = find(E.use & E.class == c); kc = kc(tr_pick(numel(kc), min(3000, numel(kc)), 100 + c));
    tr_plane_panel(ax, Q, gi, kc, tr_colors(c), FS-1, 4);
    plot(ax, W.x(k), W.y(k), 'o', 'MarkerSize', 4.5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 0.6);
    lab = 'ab';
    for e = 1:2
        plot(ax, W.x(ex(e)), W.y(ex(e)), 'p', 'MarkerSize', 15, 'MarkerFaceColor', [0.98 0.80 0.20], 'MarkerEdgeColor', 'k');
        text(ax, W.x(ex(e)) + 0.25, W.y(ex(e)) + 0.35, lab(e), 'FontSize', FS + 1, 'FontWeight', 'bold');
    end
    xlabel(ax, 'x = log_{10}(k_2K_{TOT}/k_8P_{TOT})', 'FontSize', FS-1); ylabel(ax, 'y = log_{10}(k_4K_{TOT}/k_6P_{TOT})', 'FontSize', FS-1);
    title(ax, {sprintf('Class %d: %s', c, names{q,2}), sprintf('open circles: %d sets, %.1f%% of the class', numel(k), 100*numel(k)/W.nClass(c))}, ...
          'FontSize', FS + 1, 'Color', tr_colors(c));
    for e = 1:2
        p = E.p(W.idx(ex(e)), :);
        Amax = 10^(ceil(log10(W.detFold(ex(e)))) + 1.5);
        D = tr_curves_ext(p, 1e-2, Amax, 700);
        a = axes('Position', [0.37 + 0.315*(e-1), y0, 0.27, hA]);
        tr_draw_diagram(a, D, struct('xlim', [1e-2 Amax], 'lw', 1.5, 'ms', 6));
        tr_shade_window(a, Amax); tr_tick_decades(a); set(a, 'FontSize', FS-1);
        xlabel(a, 'A_{TOT}', 'FontSize', FS-1); ylabel(a, 'x_2 = [APP]', 'FontSize', FS-1);
        title(a, {sprintf('set %s: detached fold at A_{TOT} = 10^{%.1f}  (x = %+.2f, y = %+.2f)', lab(e), log10(D.detFold), W.x(ex(e)), W.y(ex(e))), ...
                  sprintf('Class %d for A_{TOT} \\leq 10^6; with the fold inside the window: Class %d', D.class, D.classWide)}, 'FontSize', FS);
    end
end
axH = axes('Position', [0 0.955 1 0.045]); axis(axH, 'off'); set(axH, 'XLim', [0 1], 'YLim', [0 1]);
text(axH, 0.5, 0.65, 'Open circles: sets in the detached-branch quadrant labelled 1, 3, 9 or 10', 'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
text(axH, 0.5, 0.15, 'The detached branch exists (x < 0, y > 0), but its fold lies above A_{TOT} = 10^6, so within the labelling window only the main branch is seen', ...
     'HorizontalAlignment', 'center', 'FontSize', FS + 1);
axL = axes('Position', [0 0.003 1 0.022]); axis(axL, 'off'); hold(axL, 'on'); set(axL, 'XLim', [0 1], 'YLim', [0 1]);
plot(axL, 0.05, 0.6, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
text(axL, 0.06, 0.6, 'open circle: detached fold above A_{TOT} = 10^6', 'FontSize', FS);
plot(axL, 0.29, 0.6, 'p', 'MarkerSize', 11, 'MarkerFaceColor', [0.98 0.80 0.20], 'MarkerEdgeColor', 'k');
text(axL, 0.30, 0.6, 'example set (right)', 'FontSize', FS);
it = {'main branch', tr_colors('main'), '-'; 'detached branch', tr_colors('detached'), '-'; 'stable', [0 0 0], '-'; 'unstable', [0 0 0], '--'};
xp = [0.42 0.53 0.66 0.73];
for i = 1:4
    plot(axL, xp(i) + [0 0.025], [0.6 0.6], it{i,3}, 'Color', it{i,2}, 'LineWidth', 1.8); text(axL, xp(i) + 0.03, 0.6, it{i,1}, 'FontSize', FS);
end
plot(axL, 0.81, 0.6, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.96 0.82 0.25], 'MarkerEdgeColor', 'k'); text(axL, 0.82, 0.6, 'fold (LP)', 'FontSize', FS);
patch(axL, [0.88 0.9 0.9 0.88], [0.4 0.4 0.8 0.8], [0.9 0.9 0.9], 'EdgeColor', [0.6 0.6 0.6]); text(axL, 0.905, 0.6, 'A_{TOT} > 10^6', 'FontSize', FS);
tr_save(fig, fullfile(figDir, 'S_window_open_circles'));
fprintf('S_window_open_circles saved in %s\n', figDir);

%% ======================= PART 2: classifier errors from the window =====================
close all;
fMis = fullfile(PARENT, 'ML_v3_results_parallel', 'labels_Exact_Class_v3', 'misclassified', 'Misclassified_test.csv');

% --- read the misclassified test sets (Description column may contain commas) ---
fid = fopen(fMis, 'r'); fgetl(fid); M = zeros(0, 18);
while true
    s = fgetl(fid); if ~ischar(s), break; end
    f = strsplit(s, ',');
    if numel(f) < 20, continue; end
    M(end+1, :) = [str2double(f([1 2 3 5])), str2double(f(7:20))]; %#ok<SAGROW>
end
fclose(fid);
it = M(:,1); tru = M(:,3); prd = M(:,4); P = M(:, 5:18);
Qm = tr_conditions(P);
df = Inf(size(P, 1), 1);
cand = find(Qm.det & ~Qm.sym & ismember(tru, [1 3 9 10]));
for j = cand(:)'
    [~, S] = tr_label(P(j, :), 4000); if ~isempty(S), df(j) = S.detached_fold_ATOT; end
end
win = Qm.det & ~Qm.sym & ismember(tru, [1 3 9 10]) & isfinite(df) & df > 1e6;

% --- confusion table ---
pairs = unique([tru prd], 'rows'); T = zeros(size(pairs, 1), 4);
for q = 1:size(pairs, 1)
    k = tru == pairs(q,1) & prd == pairs(q,2);
    T(q, :) = [pairs(q, :), sum(k), sum(k & win)];
end
T = sortrows(T, [-3 1 2]);
tr_write_csv(fullfile(resDir, 'TR_window_confusions.csv'), 'true,predicted,errors,window_sets', T);
% (E and W from PART 1)
fo = fopen(fullfile(resDir, 'TR_window_misclassified.txt'), 'w');
for f = [1 fo]
    fprintf(f, 'Pooled LightGBM, held-out test sets: %d errors; window errors: %d (%.1f%%)\n', numel(tru), sum(win), 100*mean(win));
    fprintf(f, 'Window sets in the whole ensemble: %d of %d ML sets (%.2f%%)\n', numel(W.idx), W.nUse, 100*numel(W.idx)/W.nUse);
    fprintf(f, 'true -> predicted   errors   window sets\n');
    for q = 1:size(T, 1)
        fprintf(f, '%4d -> %-4d       %5d    %5d (%3.0f%%)\n', T(q,1), T(q,2), T(q,3), T(q,4), 100*T(q,4)/T(q,3));
    end
    fprintf(f, 'Window errors: log10 ATOT of the detached fold %.2f-%.2f (median %.2f)\n', ...
        min(log10(df(win))), max(log10(df(win))), median(log10(df(win))));
end
fclose(fo);

% --- figure: confusions with window errors, up to three sets each ---------------------
show = T(T(:,4) >= 2, :); show = show(1:min(4, end), :);
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [1 1 40 13*size(show,1) + 5]);
nR = size(show, 1); h = 0.87/nR;
for q = 1:nR
    k = find(tru == show(q,1) & prd == show(q,2) & win);
    [~, o] = sort(df(k)); k = k(o); k = k(unique(round(linspace(1, numel(k), min(3, numel(k))))));
    y0 = 0.91 - q*h;
    axT = axes('Position', [0.04, y0 + 0.86*h, 0.9, 0.04*h/0.2]); axis(axT, 'off'); set(axT, 'XLim', [0 1], 'YLim', [0 1]);
    text(axT, 0, 0.5, sprintf('Class %d predicted as Class %d: %d test sets, %d with the detached fold above A_{TOT} = 10^6', ...
        show(q,1), show(q,2), show(q,3), show(q,4)), 'FontSize', FS + 1, 'FontWeight', 'bold');
    for e = 1:numel(k)
        j = k(e); Amax = 10^(ceil(log10(df(j))) + 1.5);
        D = tr_curves_ext(P(j, :), 1e-2, Amax, 600);
        a = axes('Position', [0.06 + 0.315*(e-1), y0 + 0.22*h, 0.26, 0.50*h]);
        tr_draw_diagram(a, D, struct('xlim', [1e-2 Amax], 'lw', 1.5, 'ms', 6));
        tr_shade_window(a, Amax); tr_tick_decades(a); set(a, 'FontSize', FS-1);
        xlabel(a, 'A_{TOT}', 'FontSize', FS-1); ylabel(a, 'x_2 = [APP]', 'FontSize', FS-1);
        title(a, sprintf('true %d, predicted %d: detached fold at A_{TOT} = 10^{%.1f}  (x = %+.2f, y = %+.2f)', ...
            show(q,1), show(q,2), log10(df(j)), Qm.x(j), Qm.y(j)), 'FontSize', FS);
    end
end
axH = axes('Position', [0 0.93 1 0.07]); axis(axH, 'off'); set(axH, 'XLim', [0 1], 'YLim', [0 1]);
text(axH, 0.5, 0.7, 'Misclassified test sets whose detached branch lies beyond the labelling window', 'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
text(axH, 0.5, 0.25, 'The classifier detects the detached branch (x < 0, y > 0); the label ignores it because its fold lies above A_{TOT} = 10^6 (grey)', ...
     'HorizontalAlignment', 'center', 'FontSize', FS + 1);
tr_save(fig, fullfile(figDir, 'S_window_misclassified'));
fprintf('S_window_misclassified saved in %s\n', figDir);
