function D = tr_curves_ext(p, Amin, Amax, nPts)
%TR_CURVES_EXT  Like tr_curves, but the folds are located up to Amax (> 1e6), so that a
%   detached fold beyond the labelling window is drawn. D.class is the class in the
%   labelling window [0, 1e6]; D.classWide the class in [0, Amax]; D.detFold the ATOT of
%   the detached fold (Inf if none).
if nargin < 4, nPts = 600; end
c = tr_label(p);
[cw, S] = tr_label_win(p, Amax);
D = struct('class', c, 'classWide', cw, 'S', S, 'detFold', Inf, ...
           'br', struct('type', {}, 'ATOT', {}, 'x2', {}, 'stable', {}), 'lp', zeros(0, 2), 'bp', zeros(0, 2));
if isempty(S), return; end
D.detFold = S.detached_fold_ATOT;
for i = 1:numel(S.branch)
    b = S.branch(i);
    idx = find(b.ATOT >= Amin & b.ATOT <= Amax);
    if isempty(idx), continue; end
    idx = idx(unique(round(linspace(1, numel(idx), min(nPts, numel(idx))))));
    st = false(1, numel(idx));
    for q = 1:numel(idx)
        st(q) = max(real(eig(jac6(b.X(:, idx(q)), [b.ATOT(idx(q)) p(:)'])))) < 0;
    end
    D.br(end+1) = struct('type', b.type, 'ATOT', b.ATOT(idx), 'x2', b.X(2, idx), 'stable', st);
end
for f = S.folds, D.lp(end+1, :) = [f.ATOT, f.x(2)]; end
for b = S.bp,    D.bp(end+1, :) = [b.ATOT, b.x(2)]; end
end
