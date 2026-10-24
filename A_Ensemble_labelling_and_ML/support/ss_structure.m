function S = ss_structure(P, opts)
%SS_STRUCTURE  Complete steady-state structure of the dual futile cycle (MAPKfull.m)
%   as a function of ATOT, for one parameter vector, with no initial conditions.
%
%   S = ss_structure(P)          P = [ATOT KTOT PTOT k1 k11 k2 k3 k33 k4 k5 k55 k6 k7 k77 k8]
%   S = ss_structure(P, opts)    opts.ATOT_max (window upper end, default 1e4)
%                                opts.ATOT_min (window lower end, default 0)
%                                opts.N        (grid points per branch, default 6000)
%
%   METHOD (exact, after Wang & Sontag 2008). At steady state the complexes are
%   x3 = A*K/Km1, x4 = x1*K/Km3, x5 = x2*P/Km5, x6 = x1*P/Km7, and flux balance gives
%   x1 = a1*u*A, x2 = a1*a2*u^2*A with u = K/P. The two enzyme conservation laws then
%   give A as an explicit function of u:
%         A(u) = (u - r) / ( u*(alpha*u + beta) ),   r = KTOT/PTOT.
%   So every positive steady state, for EVERY value of ATOT, is one point u > 0 with
%   A(u) > 0, and ATOT = G(u) is explicit. Counting steady states at a given ATOT means
%   counting solutions of G(u) = ATOT; folds (saddle-nodes, LP) are extrema of G.
%   The admissible u-set has one or two intervals:
%     'main'      interval with endpoint u = r: the branch that starts at ATOT = 0
%     'detached'  interval with both ends at A -> Inf: a U-shaped branch that begins at a
%                 fold at ATOT = min(G) and exists for all larger ATOT. It is NOT
%                 connected to the main branch, so continuation from low ATOT misses it.
%   Degenerate case alpha*r + beta = 0 (the kinase/phosphatase-symmetric ensembles):
%   u = r is itself a branch ('sym', any A) and A = 1/(alpha*u) is a second branch
%   ('offsym'); they cross at a branch point (BP).
%
%   Output fields (ATOT values are restricted to [ATOT_min, ATOT_max] unless noted):
%     S.branch(i)      .type, .ATOT, .X (6 x N states), .z (grid), .map (z -> u, A)
%     S.folds          table-like struct array: ATOT, x (6x1), branch, kind ('min'/'max')
%     S.bp             struct array: ATOT, x
%     S.intervals      [ATOT_lo ATOT_hi n_ss n_stable] for each interval between
%                      consecutive critical points in the window
%     S.max_ss         max number of steady states in the window
%     S.multistationary  any interval with >= 3 steady states
%     S.bistable       any interval with >= 2 stable steady states
%     S.has_detached   a detached branch exists whose fold lies in the window
%     S.detached_fold_ATOT  (anywhere, Inf if none)
%     S.degenerate     symmetric (alpha*r + beta = 0) case

if nargin < 2, opts = struct(); end
if ~isfield(opts,'ATOT_max'), opts.ATOT_max = 1e4; end
if ~isfield(opts,'ATOT_min'), opts.ATOT_min = 0;   end
if ~isfield(opts,'N'),        opts.N = 6000;       end

C = model_constants(P);
S = struct();
S.opts = opts;
S.constants = C;
S.degenerate = abs(C.alpha*C.r + C.beta) <= 1e-9*(abs(C.alpha)*C.r + abs(C.beta));

%% ---- build branches -------------------------------------------------------------
br = struct('type',{},'z',{},'map',{},'ATOT',{},'X',{});
if ~S.degenerate
    cands = [0, C.r];
    if C.alpha ~= 0
        us = -C.beta/C.alpha;
        if us > 0 && isfinite(us), cands(end+1) = us; end
    end
    cands = unique(cands);
    edges = [cands, Inf];
    for j = 1:numel(edges)-1
        lo = edges(j); hi = edges(j+1);
        if isfinite(hi), mid = (lo+hi)/2; else, mid = 2*lo + 1; end
        [~, Amid] = u_to_A(mid, mid - C.r, mid - C.ustar, C);
        if ~(Amid > 0 && isfinite(Amid)), continue; end
        if lo == C.r || hi == C.r, t = 'main'; else, t = 'detached'; end
        br(end+1) = new_branch(t, @(z) interval_map(z, lo, hi, C), linspace(-40, 40, opts.N)); %#ok<AGROW>
    end
else
    % symmetric branch: u = r, A free
    Asc = max(C.KTOT, C.PTOT);
    br(end+1) = new_branch('sym', @(z) deal(C.r + 0*z, Asc*exp(z)), linspace(-40, 25, opts.N));
    if C.alpha > 0
        br(end+1) = new_branch('offsym', @(z) deal(C.r*exp(z), 1./(C.alpha*C.r*exp(z))), ...
                               linspace(-40, 40, opts.N));
    end
end

for i = 1:numel(br)
    [u, A] = br(i).map(br(i).z);
    [X, G] = ss_state_from_uA(u, A, C);
    ok = isfinite(G) & G > 0 & all(isfinite(X),1) & all(X >= 0,1);
    br(i).z = br(i).z(ok); br(i).ATOT = G(ok); br(i).X = X(:,ok);
end
S.branch = br;

%% ---- folds (extrema of ATOT along each branch) and branch points ----------------
Gcap = 10*opts.ATOT_max;
folds = struct('ATOT',{},'x',{},'branch',{},'kind',{});
for i = 1:numel(br)
    G = br(i).ATOT; z = br(i).z;
    if numel(G) < 3, continue; end
    d = diff(G);
    sg = sign(d);
    for m = 2:numel(sg)
        if sg(m) == 0 || sg(m-1) == 0 || sg(m) == sg(m-1), continue; end
        if G(m) > Gcap, continue; end
        if sg(m-1) < 0, kind = 'min'; s = 1; else, kind = 'max'; s = -1; end
        fz = @(zz) s*G_of_z(zz, br(i).map, C);
        zf = fminbnd(fz, z(m-1), z(m+1), optimset('TolX',1e-12));
        [uf, Af] = br(i).map(zf);
        if strcmp(br(i).type,'offsym') && abs(uf - C.r) < 1e-6*C.r
            continue   % extremum at the symmetric point is the BP, not a fold
        end
        [xf, Gf] = ss_state_from_uA(uf, Af, C);
        folds(end+1) = struct('ATOT',Gf,'x',xf,'branch',i,'kind',kind); %#ok<AGROW>
    end
end
S.folds_all = folds;
inw = arrayfun(@(f) f.ATOT >= opts.ATOT_min && f.ATOT <= opts.ATOT_max, folds);
S.folds = folds(inw);

bp = struct('ATOT',{},'x',{});
if S.degenerate && C.alpha > 0
    [xb, Gb] = ss_state_from_uA(C.r, 1/(C.alpha*C.r), C);
    bp(1).ATOT = Gb; bp(1).x = xb;
end
S.bp_all = bp;
S.bp = bp(arrayfun(@(b) b.ATOT >= opts.ATOT_min && b.ATOT <= opts.ATOT_max, bp));

S.detached_fold_ATOT = Inf;
for i = 1:numel(br)
    if strcmp(br(i).type,'detached') && ~isempty(br(i).ATOT)
        S.detached_fold_ATOT = min([S.detached_fold_ATOT, folds([folds.branch]==i).ATOT, min(br(i).ATOT)]);
    end
end
S.has_detached = S.detached_fold_ATOT <= opts.ATOT_max && S.detached_fold_ATOT >= opts.ATOT_min;

%% ---- count steady states and stable ones between critical points ----------------
crit = [opts.ATOT_min, opts.ATOT_max, [folds.ATOT], [bp.ATOT]];
crit = unique(crit(crit >= opts.ATOT_min & crit <= opts.ATOT_max));
S.intervals = zeros(0,4);
for j = 1:numel(crit)-1
    a1 = crit(j); a2 = crit(j+1);
    if a1 > 0, a = sqrt(a1*a2); else, a = a2/2; end
    Xs = ss_states_at(S, a);
    ns = size(Xs,2);
    nst = 0;
    Pa = P; Pa(1) = a;
    for q = 1:ns
        ev = eig(jac6(Xs(:,q), Pa));
        if max(real(ev)) < 0, nst = nst + 1; end
    end
    S.intervals(end+1,:) = [a1 a2 ns nst];
end
if isempty(S.intervals)
    S.max_ss = 0; S.multistationary = false; S.bistable = false;
else
    S.max_ss = max(S.intervals(:,3));
    S.multistationary = any(S.intervals(:,3) >= 3);
    S.bistable = any(S.intervals(:,4) >= 2);
end
S.bistable_ATOT = S.intervals(S.intervals(:,4) >= 2, 1:2);
end

%% ===================================================================================
function b = new_branch(t, map, z)
b = struct('type',t,'z',z,'map',[],'ATOT',[],'X',[]);
b.map = map;
end

function G = G_of_z(z, map, C)
[u, A] = map(z);
[~, G] = ss_state_from_uA(u, A, C);
end

function [u, A] = interval_map(z, lo, hi, C)
% maps z in R to u in (lo,hi), keeping the distance to the nearer endpoint exact,
% so the singular factors (u - r) and (u - u*) are computed without cancellation.
u = zeros(size(z)); du_r = u; du_s = u;
if isfinite(hi)
    L = hi - lo;
    neg = z <= 0;
    d = L./(1 + exp(-z(neg)));            % distance from lo
    u(neg) = lo + d; du_r(neg) = (lo - C.r) + d; du_s(neg) = (lo - C.ustar) + d;
    d = L./(1 + exp(z(~neg)));            % distance from hi
    u(~neg) = hi - d; du_r(~neg) = (hi - C.r) - d; du_s(~neg) = (hi - C.ustar) - d;
else
    d = lo*exp(z);
    u = lo + d; du_r = (lo - C.r) + d; du_s = (lo - C.ustar) + d;
end
[~, A] = u_to_A(u, du_r, du_s, C);
end

function [u, A] = u_to_A(u, du_r, du_s, C)
if C.alpha ~= 0
    A = du_r ./ (u .* C.alpha .* du_s);
else
    A = du_r ./ (u .* C.beta);
end
end