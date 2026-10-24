function [L, LP, BP] = exact_labels(S, P0)
%EXACT_LABELS  Ten-class label and normal forms from the exact steady-state structure.
%   [L, LP, BP] = exact_labels(S, P0)   S = ss_structure(P0, ...)
%   LP, BP: struct arrays (ATOT, x2, type, x, nf) of the folds and branch points inside
%   the window; nf = normal_forms(...) of each point (empty if P0 is not given).
LP = struct('ATOT',{},'x2',{},'type',{},'x',{},'nf',{});
for f = S.folds
    nf = [];
    if nargin > 1
        Pf = P0; Pf(1) = f.ATOT;
        try, nf = normal_forms(f.x(:), Pf, 'LP'); catch, nf = []; end
    end
    LP(end+1) = struct('ATOT', f.ATOT, 'x2', f.x(2), 'type', S.branch(f.branch).type, ...
                       'x', f.x(:), 'nf', nf); %#ok<AGROW>
end
BP = struct('ATOT',{},'x2',{},'type',{},'x',{},'nf',{});
for b = S.bp
    nf = [];
    if nargin > 1
        Pb = P0; Pb(1) = b.ATOT;
        try, nf = normal_forms(b.x(:), Pb, 'BP'); catch, nf = []; end
    end
    BP(end+1) = struct('ATOT', b.ATOT, 'x2', b.x(2), 'type', 'sym', 'x', b.x(:), 'nf', nf); %#ok<AGROW>
end
L = classify_branches(LP, BP, S);
L.multistationary = S.multistationary;
L.bistable = S.bistable;
L.max_ss = S.max_ss;
L.has_detached = S.has_detached;
L.detached_fold_ATOT = S.detached_fold_ATOT;
L.degenerate = S.degenerate;
L.bistable_ATOT = S.bistable_ATOT;
end
