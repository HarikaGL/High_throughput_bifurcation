function xl = auto_xlim(rs, factor, defaultMax)
%AUTO_XLIM  Linear ATOT range for a bifurcation diagram: [0, factor * largest ATOT of a
%   fold, branch point or (biphasic) peak] over the samples rs; [0 defaultMax] if none.
if nargin < 2 || isempty(factor), factor = 2.5; end
if nargin < 3 || isempty(defaultMax), defaultMax = 200; end
m = 0;
for r = rs(:)'
    m = max([m; r.exact_LP_ATOT(:); r.exact_BP_ATOT(:)]);
    if strcmp(char(r.exact_label), 'NO_SING_biphasic') && ~isempty(r.curves)
        types = arrayfun(@(cb) char(cb.type), r.curves, 'UniformOutput', false);
        b = find(strcmp(types, 'main'), 1); if isempty(b), b = 1; end
        [~, i] = max(r.curves(b).x2); m = max(m, r.curves(b).ATOT(i));
    end
end
if m > 0, xl = [0 factor*m]; else, xl = [0 defaultMax]; end
end
