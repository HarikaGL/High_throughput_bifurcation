function tr_tick_decades(ax)
%TR_TICK_DECADES  Decade ticks 10^k on both log axes, at most 6 labels per axis.
% decade ticks 10^k, at most 6 labels per axis
for w = 'XY'
    L = get(ax, [w 'Lim']); k = ceil(log10(L(1))):floor(log10(L(2)));
    st = max(1, ceil(numel(k)/6)); k = k(mod(k - k(end), st) == 0);
    set(ax, [w 'Tick'], 10.^k, [w 'TickLabel'], arrayfun(@(v) sprintf('10^{%d}', v), k, 'UniformOutput', false));
end
end

