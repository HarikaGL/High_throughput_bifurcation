function [xl, yl, isLog] = diagram_window(R, scaleMode, logZoom)
%DIAGRAM_WINDOW  Axis limits that show the shape of one bifurcation diagram clearly.
%   [xl, yl, isLog] = diagram_window(R, scaleMode, logZoom)
%   scaleMode 'lin': linear axes zoomed on the folds / branch points (narrow
%                    hysteresis loops visible); no singularity -> auto_xlim range
%             'log': log-log axes, ATOT from min(singular ATOT)/logZoom to
%                    max(singular ATOT)*logZoom, x2 over the data in that range
%   Limits are returned in data units (not log10). isLog is false when 'log' was
%   asked for a sample without any singularity.
if nargin < 3 || isempty(logZoom), logZoom = 30; end
A = [R.exact_LP_ATOT(:); R.exact_BP_ATOT(:)]; X = [R.exact_LP_x2(:); R.exact_BP_x2(:)];
isLog = strcmpi(scaleMode, 'log') && ~isempty(A);
if isLog
    xl = [min(A)/logZoom, max(A)*logZoom];
    v = [];
    for cb = R.curves(:)'
        in = cb.ATOT >= xl(1) & cb.ATOT <= xl(2) & cb.x2 > 0 & isfinite(cb.x2);
        v = [v, cb.x2(in)]; %#ok<AGROW>
    end
    if isempty(v), yl = [1e-3 1e3]; return; end
    ly = log10([min(v), max(v)]); pad = 0.06*max(diff(ly), 1);
    yl = 10.^(ly + [-pad pad]);
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
