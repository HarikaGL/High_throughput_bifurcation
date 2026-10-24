function tr_shade_window(ax, Amax)
%TR_SHADE_WINDOW  Grey band ATOT > 1e6 (outside the labelling window) and a dotted line at 1e6.
yl = get(ax, 'YLim');
h = patch(ax, [1e6 Amax Amax 1e6], yl([1 1 2 2]), [0.90 0.90 0.90], 'EdgeColor', 'none');
try, uistack(h, 'bottom'); catch, ch = get(ax, 'Children'); set(ax, 'Children', [ch(ch ~= h); h]); end
plot(ax, [1e6 1e6], yl, 'k:', 'LineWidth', 0.9);
set(ax, 'YLim', yl, 'Layer', 'top');
end
