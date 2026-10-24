function D = tr_curves(p, Amin, Amax, nPts, N)
%TR_CURVES  Exact bifurcation diagram (x2 against ATOT) of p for plotting.
%   D.br(i): type ('main','detached','sym','offsym'), ATOT, x2, stable (eigenvalues of the
%            analytic Jacobian jac6), restricted to [Amin, Amax]
%   D.lp   : [ATOT x2] of the folds in the window [0, 1e6]; D.bp: [ATOT x2] of branch points
%   D.class: exact class; D.S: the ss_structure result
if nargin < 4, nPts = 500; end
if nargin < 5, N = 10000; end
[c, S] = tr_label(p, N);
D = struct('class', c, 'S', S, 'br', struct('type', {}, 'ATOT', {}, 'x2', {}, 'stable', {}), ...
           'lp', zeros(0, 2), 'bp', zeros(0, 2));
if isempty(S), return; end
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
