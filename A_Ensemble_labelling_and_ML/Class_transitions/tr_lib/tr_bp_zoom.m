function Z = tr_bp_zoom(p, w, n)
%TR_BP_ZOOM  Dense branches near the branch point of a symmetric set (u* = r):
%   the u = r branch for A in A_bp*10^[-w, w] and the crossing branch A = 1/(alpha u) for
%   u in r*10^[-w, w]; with stability. Used for the insets of Classes 7 and 8.
if nargin < 2, w = 0.25; end
if nargin < 3, n = 4000; end
P = [1 p(:)']; C = model_constants(P);
Abp = 1/(C.alpha*C.r);
s = logspace(-w, w, n);
Z = struct('type', {}, 'ATOT', {}, 'x2', {}, 'stable', {});
[X, G] = ss_state_from_uA(C.r*ones(1, n), Abp*s, C);   Z(1) = mk('sym', X, G);
u = C.r*s; [X, G] = ss_state_from_uA(u, 1./(C.alpha*u), C); Z(2) = mk('offsym', X, G);
    function z = mk(t, X, G)
        st = false(1, numel(G));
        for i = 1:numel(G), st(i) = max(real(eig(jac6(X(:, i), [G(i) p(:)'])))) < 0; end
        z = struct('type', t, 'ATOT', G, 'x2', X(2, :), 'stable', st);
    end
end
