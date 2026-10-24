function tr_plane_panel(ax, Q, gi, kc, colc, FS, ms)
%TR_PLANE_PANEL  One Fig-11-style panel in the (x, y) plane: pink detached-branch quadrant,
%   grey ensemble (rows gi), the class (rows kc, colour colc), dashed x = 0 and y = 0 and
%   the grey diagonal y = -x of the symmetric sets of Iterations 3-4.
if nargin < 7, ms = 4; end
hold(ax, 'on'); box(ax, 'on');
patch(ax, [-6 0 0 -6], [0 0 7 7], [0.98 0.87 0.85], 'EdgeColor', 'none');
plot(ax, Q.x(gi), Q.y(gi), '.', 'Color', [0.80 0.80 0.80], 'MarkerSize', 3);
if ~isempty(kc), plot(ax, Q.x(kc), Q.y(kc), '.', 'Color', colc, 'MarkerSize', ms); end
plot(ax, [-6 6], [0 0], 'k:', 'LineWidth', 0.8); plot(ax, [0 0], [-6 7], 'k:', 'LineWidth', 0.8);
set(ax, 'XLim', [-5.5 4.5], 'YLim', [-4.5 6], 'FontSize', FS, 'Layer', 'top');
end
