function L = classify_branches(LP, BP, S_exact)
%CLASSIFY_BRANCHES  Ten-class bifurcation label (v3) from fold (LP) and branch-point
%   (BP) lists. Used twice with the same rules, so the two labellings are directly
%   comparable:
%     - with the folds/BPs of the exact steady-state structure (ss_structure)
%     - with the LP/BP points that MATCONT finds on the seeded branches
%
%   LP, BP : struct arrays with fields ATOT, x2 (free doubly phosphorylated substrate,
%            state 2) and type (branch type: 'main','detached','sym','offsym').
%            Only points inside the ATOT window.
%   S_exact: output of ss_structure (used for the biphasic test on the main branch).
%
%   Classes (see class_info_v3.m):
%     1 LP_S_ONLY          two folds on the main branch, S (x2 of the lower-ATOT fold
%                          is the larger), no fold on a detached branch
%     2 LP_S_DETACHED      as 1, plus a fold on a detached branch
%     3 LP_Z_ONLY          two folds on the main branch, Z orientation, no detached fold
%     4 LP_Z_DETACHED      as 3, plus a fold on a detached branch
%     5 LP_DETACHED_ONLY   fold(s) only on a detached branch
%     6 BP_ONLY            branch point, no fold
%     7 LP_BP_LPgreaterBP  branch point and fold, mean x2(LP) > mean x2(BP)
%     8 LP_BP_LPlessBP     branch point and fold, mean x2(LP) < mean x2(BP)
%     9 NO_SING_biphasic   no singular point, x2 along the main branch has one maximum
%    10 NO_SING_monotone   no singular point, not biphasic
%    11 OTHER              anything else (e.g. a fold pair split by the window edge)
%
%   Also returned, for traceability with Ensemble_LHS_v2:
%     L.class_v2 / L.label_v2   the nine v2 classes, L.class_old7 the original seven.
%   Flags:
%     L.lp_coincides_bp   LP+BP set whose LP coincides with the BP (relative 1e-9):
%                         the turning point of the crossing branch at a pitchfork,
%                         not a separate fold. Reported, but excluded from the ML.
%     L.detached_in_window  a fold on a detached branch lies inside the window.

CI = class_info_v3();
isMain = @(t) any(strcmp(t, {'main','sym'}));
nLP = numel(LP); nBP = numel(BP);
if nLP > 0
    lpMain = LP(arrayfun(@(p) isMain(p.type), LP));
    lpDet  = LP(arrayfun(@(p) strcmp(p.type,'detached'), LP));
else
    lpMain = LP; lpDet = LP;
end
hasD = ~isempty(lpDet);

L = struct();
L.n_LP = nLP; L.n_BP = nBP; L.n_LP_main = numel(lpMain); L.n_LP_detached = numel(lpDet);

% ---- v2 rules (unchanged) ----------------------------------------------------------
if nBP > 0 && nLP > 0
    if mean([LP.x2]) > mean([BP.x2]), c2 = 4; else, c2 = 5; end
elseif nBP > 0
    c2 = 3;
elseif numel(lpMain) >= 2
    [~, o] = sort([lpMain.ATOT]);
    x2s = [lpMain(o).x2];
    if x2s(1) > x2s(2), c2 = 1; else, c2 = 2; end
elseif hasD && isempty(lpMain)
    c2 = 7;
elseif nLP > 0
    c2 = 9;
else
    if is_biphasic(S_exact), c2 = 6; else, c2 = 8; end
end

% ---- v3: regroup and split the main-branch hysteresis by the detached fold ---------
map = [1 3 6 7 8 9 5 10 11];          % v2 class -> v3 class
c = map(c2);
if c == 1 && hasD, c = 2; end
if c == 3 && hasD, c = 4; end

L.class      = c;
L.label      = CI.names{c};
L.class_v2   = c2;
L.label_v2   = CI.names_v2{c2};
L.class_old7 = min(c2, 7);
L.detached_in_window = hasD;
L.lp_coincides_bp = false;
if any(c == [7 8])
    a = [LP.ATOT]; b = BP(1).ATOT;
    L.lp_coincides_bp = min(abs(a - b)) <= 1e-9*abs(b);
end
end

function tf = is_biphasic(S)
% one interior maximum of x2 along the main branch, ordered by ATOT, inside the window
tf = false;
if isempty(S) || ~isfield(S,'branch'), return; end
for i = 1:numel(S.branch)
    b = S.branch(i);
    if ~any(strcmp(b.type, {'main','sym'})), continue; end
    in = b.ATOT >= S.opts.ATOT_min & b.ATOT <= S.opts.ATOT_max;
    a = b.ATOT(in); x2 = b.X(2,in);
    if numel(a) < 5, return; end
    [~, o] = sort(a); x2 = x2(o);
    d = diff(x2);
    tol = 1e-9*max(abs(x2));
    d(abs(d) <= tol) = 0;
    d = d(d ~= 0);
    nmax = sum(d(1:end-1) > 0 & d(2:end) < 0);
    tf = (nmax == 1);
    return
end
end
