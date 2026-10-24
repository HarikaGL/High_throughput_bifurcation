function tr_draw_diagram(ax, D, opt)
%TR_DRAW_DIAGRAM  Draw a diagram from tr_curves on axes ax (log-log).
%   Colours: main branch blue, detached branch red, symmetric branch u = r dark grey,
%   crossing branch orange; solid stable, dashed unstable; yellow diamonds folds, black
%   dots branch points. opt.xlim (default [1e-2 1e6]), opt.lw (1.6), opt.zoom (struct from
%   tr_bp_zoom; draws an inset around the branch point), opt.ms (marker size, 6).
if nargin < 3, opt = struct(); end
if ~isfield(opt, 'xlim'), opt.xlim = [1e-2 1e6]; end
if ~isfield(opt, 'lw'), opt.lw = 1.6; end
if ~isfield(opt, 'ms'), opt.ms = 6; end
hold(ax, 'on'); box(ax, 'on');
yy = [];
for b = D.br(:)'
    plot_branch(ax, b, opt.lw);
    in = b.ATOT >= opt.xlim(1) & b.ATOT <= opt.xlim(2) & b.x2 > 0;
    yy = [yy, b.x2(in)]; %#ok<AGROW>
end
for i = 1:size(D.lp, 1)
    if D.lp(i,1) >= opt.xlim(1) && D.lp(i,1) <= opt.xlim(2)
        plot(ax, D.lp(i,1), D.lp(i,2), 'd', 'MarkerSize', opt.ms, 'MarkerFaceColor', [0.96 0.82 0.25], 'MarkerEdgeColor', 'k');
    end
end
for i = 1:size(D.bp, 1)
    plot(ax, D.bp(i,1), D.bp(i,2), 'o', 'MarkerSize', opt.ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
end
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', opt.xlim);
if ~isempty(yy)
    lo = max(min(yy), 1e-12); hi = max(yy);
    set(ax, 'YLim', 10.^[floor(log10(lo)) - 0.3, ceil(log10(hi)) + 0.3]);
end
if isfield(opt, 'zoom') && ~isempty(opt.zoom) && ~isempty(D.bp)
    draw_zoom(ax, D, opt);
end
end

function plot_branch(ax, b, lw)
col = tr_colors(b.type);
st = logical(b.stable);
ys = b.x2; ys(~(st | [st(2:end) false] | [false st(1:end-1)])) = NaN;
yu = b.x2; yu(st) = NaN;
uu = ~st | [~st(2:end) false] | [false ~st(1:end-1)]; yu(uu) = b.x2(uu);
plot(ax, b.ATOT, ys, '-', 'Color', col, 'LineWidth', lw);
plot(ax, b.ATOT, yu, '--', 'Color', col, 'LineWidth', lw*0.8);
end

function draw_zoom(ax, D, opt)
% inset around the branch point, wide enough to show the nearby fold (Classes 7, 8)
Z = opt.zoom; bp = D.bp(1, :);
lp = D.lp; dl = 3e-4;
if ~isempty(lp), dl = max(4*min(abs(lp(:,1) - bp(1)))/bp(1), dl); end
x0 = bp(1)*(1 - dl); x1 = bp(1)*(1 + dl);
pos = get(ax, 'Position');
ins = axes('Position', [pos(1) + 0.06*pos(3), pos(2) + 0.55*pos(4), 0.42*pos(3), 0.4*pos(4)]);
hold(ins, 'on'); box(ins, 'on');
for z = Z(:)'
    k = z.ATOT >= x0 & z.ATOT <= x1;
    if any(k), plot_branch(ins, struct('type', z.type, 'ATOT', z.ATOT(k), 'x2', z.x2(k), 'stable', z.stable(k)), 1.2); end
end
if ~isempty(lp)
    [~, j] = min(abs(lp(:,1) - bp(1)));
    plot(ins, lp(j,1), lp(j,2), 'd', 'MarkerSize', 5, 'MarkerFaceColor', [0.96 0.82 0.25], 'MarkerEdgeColor', 'k');
end
plot(ins, bp(1), bp(2), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
set(ins, 'XLim', [x0 x1], 'YScale', 'log', 'YLim', bp(2)*[1/4 4], 'XTick', [], 'YTick', []);
title(ins, 'zoom at BP', 'FontSize', 7, 'FontWeight', 'normal');
end
