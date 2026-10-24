function h = draw_branches(ax, curves, lw)
%DRAW_BRANCHES  Plot steady-state branches x2 vs ATOT on axes ax.
%   Stable parts solid, unstable parts dashed. Colour by branch type:
%   main/sym = blue, detached = red, offsym (symmetric case) = purple.
if nargin < 3, lw = 1.4; end
col = struct('main',[0 0.35 0.75], 'sym',[0 0.35 0.75], 'detached',[0.85 0.2 0.1], 'offsym',[0.5 0.2 0.65]);
h = [];
hold(ax, 'on');
for c = curves(:)'
    if isempty(c.ATOT), continue; end
    if isfield(col, c.type), cc = col.(c.type); else, cc = [0.3 0.3 0.3]; end
    st = logical(c.stable);
    % extend each stable run by one point on each side so solid and dashed parts meet
    stS = st | [st(2:end) false] | [false st(1:end-1)];
    ys = c.x2; ys(~stS) = NaN;
    yu = c.x2; yu(st) = NaN;
    uu = ~st | [~st(2:end) false] | [false ~st(1:end-1)];
    yu(uu) = c.x2(uu);
    h(end+1) = plot(ax, c.ATOT, ys, '-',  'Color', cc, 'LineWidth', lw); %#ok<AGROW>
    h(end+1) = plot(ax, c.ATOT, yu, '--', 'Color', cc, 'LineWidth', max(0.6, lw*0.7)); %#ok<AGROW>
end
end