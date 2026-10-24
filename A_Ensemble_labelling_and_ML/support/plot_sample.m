function plot_sample(R, P, k, ax)
%PLOT_SAMPLE  Full bifurcation diagram of one sample (x2 = App vs ATOT), linear axes.
%   plot_sample(res(k), paramSamples(k,:), k)        opens a new figure
%   plot_sample(res(k), paramSamples(k,:), k, ax)    draws into existing axes
%   Solid = stable, dashed = unstable; blue = branch from ATOT = 0, red = detached
%   branch, purple = off-symmetric branch. Red circles = folds (LP), green squares =
%   branch points (BP). Thin grey = MATCONT branches, if they were computed.
%   The ATOT range is chosen from the singularities (auto_xlim.m).
if nargin < 4 || isempty(ax), figure('Color','w'); ax = gca; end
cla(ax); hold(ax, 'on'); box(ax, 'on');
xl = auto_xlim(R);
yl = y_range(R.curves, xl);
for j = 1:size(R.bistable_ATOT,1)          % shade bistable ATOT ranges (drawn first = bottom)
    patch(ax, R.bistable_ATOT(j,[1 2 2 1]), yl([1 1 2 2]), [1 0.9 0.6], ...
        'EdgeColor','none', 'FaceAlpha', 0.35, 'HandleVisibility','off');
end
draw_branches(ax, R.curves, 1.6);
if isfield(R, 'mc_curves')
    for c = R.mc_curves(:)'
        plot(ax, c.ATOT, c.x2, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.9);
    end
end
if ~isempty(R.exact_LP_ATOT)
    plot(ax, R.exact_LP_ATOT, R.exact_LP_x2, 'o', 'MarkerFaceColor','r', 'MarkerEdgeColor','k', 'MarkerSize', 6);
end
if ~isempty(R.exact_BP_ATOT)
    plot(ax, R.exact_BP_ATOT, R.exact_BP_x2, 's', 'MarkerFaceColor','g', 'MarkerEdgeColor','k', 'MarkerSize', 7);
end
xlim(ax, xl); ylim(ax, yl);
xlabel(ax, 'A_{TOT}'); ylabel(ax, 'x_2 (App)');
t = sprintf('#%d  %s  (max %d SS)', k, R.exact_label, R.max_ss);
if R.has_detached, t = sprintf('%s | detached fold at ATOT=%.3g', t, R.detached_fold_ATOT); end
if R.mushroom_candidate, t = [t ' | MUSHROOM?']; end
title(ax, t, 'Interpreter','none', 'FontSize', 9);
if nargin >= 2 && ~isempty(P)
    text(ax, 0.02, 0.97, sprintf('KTOT=%.3g  PTOT=%.3g', P(2), P(3)), 'Units','normalized', ...
        'VerticalAlignment','top', 'FontSize', 8);
end
end
