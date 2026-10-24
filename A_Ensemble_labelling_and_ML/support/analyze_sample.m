function R = analyze_sample(P0, cfg, run_mc)
%ANALYZE_SAMPLE  All work for one parameter set (v3, ten classes); safe inside parfor.
%   1) exact steady-state structure over the ATOT window (ss_structure, inverse method)
%      -> ten-class label, all folds/branch points, their normal-form coefficients
%   2) optional MATCONT continuation, seeded on EVERY branch of steady states (including
%      detached ones), in both directions -> MATCONT label (same rules), MATCONT's own
%      fold coefficient a (s.data.a) for comparison with the computed one
%
%   Returns a flat struct (same fields for every sample, so parfor can collect it).
%   The exact analysis is identical to Ensemble_LHS_v2 (same ss_structure); only the
%   labels are refined to ten classes and the normal forms are added.

R = empty_result();
t0 = tic; c0 = cputime;                  % timing (core-hours report)
eopts = struct('ATOT_min', cfg.ATOT_window(1), 'ATOT_max', cfg.ATOT_window(2), 'N', cfg.exact_N);
if ~isfield(cfg, 'curve_points'), cfg.curve_points = 600; end

%% ---------------- 1. exact structure ----------------------------------------------
try
    S = ss_structure(P0, eopts);
    tnf = tic;
    [LE, LPe, BPe] = exact_labels(S, P0);
    R.t_nf = toc(tnf);                   % labels + normal forms
    R.exact_ok        = true;
    R.exact_class     = LE.class;          % ten classes (v3)
    R.exact_label     = LE.label;
    R.exact_class_v2  = LE.class_v2;       % nine classes (v2), traceability
    R.exact_label_v2  = LE.label_v2;
    R.exact_class_old7= LE.class_old7;
    R.lp_coincides_bp = LE.lp_coincides_bp;
    R.detached_in_window = LE.detached_in_window;
    R.multistationary = LE.multistationary;
    R.bistable        = LE.bistable;
    R.max_ss          = LE.max_ss;
    R.has_detached    = LE.has_detached;
    R.detached_fold_ATOT = LE.detached_fold_ATOT;
    R.degenerate      = LE.degenerate;
    R.exact_LP_ATOT   = [LPe.ATOT];
    R.exact_LP_x2     = [LPe.x2];
    R.exact_LP_type   = {LPe.type};
    R.exact_BP_ATOT   = [BPe.ATOT];
    R.exact_BP_x2     = [BPe.x2];
    R.bistable_ATOT   = LE.bistable_ATOT;
    R.intervals       = S.intervals;
    R.n_LP_main       = LE.n_LP_main;
    R.n_LP_detached   = LE.n_LP_detached;
    R.mushroom_candidate = LE.n_LP_main >= 4;   % two separate hysteresis loops on the main branch
    % normal-form coefficients of every fold and branch point (normal_forms.m)
    for i = 1:numel(LPe)
        nf = LPe(i).nf; if isempty(nf), nf = struct('a',NaN,'beta',NaN,'kappa',NaN,'kappa_rel',NaN); end
        R.exact_LP_a(i) = nf.a; R.exact_LP_beta(i) = nf.beta;
        R.exact_LP_kappa(i) = nf.kappa; R.exact_LP_kappa_rel(i) = nf.kappa_rel;
    end
    for i = 1:numel(BPe)
        nf = BPe(i).nf;
        if isempty(nf), nf = struct('a',NaN,'b',NaN,'c',NaN,'disc',NaN,'slopes',[NaN NaN],'pitchfork',false); end
        R.exact_BP_a(i) = nf.a; R.exact_BP_b(i) = nf.b; R.exact_BP_c(i) = nf.c;
        R.exact_BP_disc(i) = nf.disc; R.exact_BP_slopes(i,:) = nf.slopes;
        R.exact_BP_pitchfork(i) = nf.pitchfork;
    end
    % branch curves for plotting (ATOT, x2) inside the window, with the stability of
    % each point (eigenvalues of the analytic Jacobian): solid = stable, dashed = unstable
    for i = 1:numel(S.branch)
        b = S.branch(i);
        in = b.ATOT >= eopts.ATOT_min & b.ATOT <= eopts.ATOT_max;
        idx = find(in);
        if isempty(idx), continue; end
        idx = idx(unique(round(linspace(1, numel(idx), min(cfg.curve_points, numel(idx))))));
        st = false(1, numel(idx));
        for q = 1:numel(idx)
            Pq = P0; Pq(1) = b.ATOT(idx(q));
            st(q) = max(real(eig(jac6(b.X(:,idx(q)), Pq)))) < 0;
        end
        R.curves(end+1) = struct('type', b.type, 'ATOT', b.ATOT(idx), 'x2', b.X(2,idx), 'stable', st); %#ok<AGROW>
    end
catch ME
    R.exact_ok = false;
    R.message = ['exact: ' ME.message];
    R.t_exact = toc(t0); R.cpu_exact = cputime - c0;
    return
end
R.t_exact = toc(t0); R.cpu_exact = cputime - c0;

if nargin < 3, run_mc = cfg.run_matcont; end
if ~run_mc, return; end

%% ---------------- 2. MATCONT, one seed per branch ---------------------------------
t1 = tic; c1 = cputime;
R.mc_ran = true;
LPm = struct('ATOT',{},'x2',{},'type',{},'a',{});
BPm = struct('ATOT',{},'x2',{},'type',{},'a',{});
HBm = [];
nSeedOK = 0;
for i = 1:numel(S.branch)
    b = S.branch(i);
    a_seed = seed_ATOT(b, S, cfg);
    if isempty(a_seed), continue; end
    [Xs, bid] = ss_states_at(S, a_seed);
    Xs = Xs(:, bid == i);
    if isempty(Xs), continue; end
    xseed = Xs(:,1);
    Pseed = P0; Pseed(1) = a_seed;
    if norm(rhs6(xseed, Pseed)) > 1e-8*max(1, a_seed)
        R.message = sprintf('%s | seed residual too large on branch %d', R.message, i);
        continue
    end
    for dirn = [0 1]
        try
            x0 = init_EP_EP(@MAPKfull, xseed, Pseed, 1);
            opt = cfg.contset;
            opt = contset(opt, 'Backward', dirn);
            [x, ~, s] = quiet_cont(x0, opt);
        catch ME
            R.message = sprintf('%s | cont failed (branch %d, dir %d): %s', R.message, i, dirn, ME.message);
            R.mc_failed_runs = R.mc_failed_runs + 1;
            continue
        end
        nSeedOK = nSeedOK + 1;
        if isempty(x), continue; end
        ATOT = x(end,:);
        phys = physical_mask(x(1:6,:), ATOT, P0, cfg);
        R.mc_ATOT_reached = [min([R.mc_ATOT_reached, ATOT(phys)]), max([R.mc_ATOT_reached, ATOT(phys)])];
        R.mc_negative = R.mc_negative || any(~phys & ATOT > 0 & ATOT <= cfg.ATOT_window(2));
        for j = 1:numel(s)
            lab = strtrim(s(j).label);
            k = s(j).index;
            if k < 1 || k > size(x,2) || ~phys(k), continue; end
            a = ATOT(k);
            if a < cfg.ATOT_window(1) || a > cfg.ATOT_window(2), continue; end
            switch lab
                case 'LP'
                    aMC = NaN;      % MATCONT's fold normal-form coefficient (a_lp.m)
                    if isfield(s(j), 'data') && isstruct(s(j).data) && isfield(s(j).data, 'a')
                        aMC = s(j).data.a;
                    end
                    LPm = add_unique(LPm, a, x(2,k), b.type, aMC);
                case 'BP'
                    BPm = add_unique(BPm, a, x(2,k), 'sym', NaN);
                case {'H','HB'}
                    HBm(end+1) = a; %#ok<AGROW>
            end
        end
        if cfg.keep_mc_curves
            keep = find(phys);
            keep = keep(unique(round(linspace(1, numel(keep), min(400, numel(keep))))));
            R.mc_curves(end+1) = struct('type', b.type, 'ATOT', ATOT(keep), 'x2', x(2,keep)); %#ok<AGROW>
        end
    end
end
R.mc_ok = nSeedOK > 0;
R.t_mc = toc(t1); R.cpu_mc = cputime - c1;
% a turning point of the off-symmetric branch at its branch point is the BP itself
if ~isempty(BPm) && ~isempty(LPm)
    keep = arrayfun(@(p) all(abs(p.ATOT - [BPm.ATOT]) > 1e-4*max(1,p.ATOT)), LPm);
    LPm = LPm(keep);
end
if R.mc_ok
    LM = classify_branches(LPm, BPm, S);   % same rules as the exact label
    R.mc_class    = LM.class;
    R.mc_label    = LM.label;
    R.mc_class_v2 = LM.class_v2;
    R.mc_label_v2 = LM.label_v2;
    R.mc_LP_ATOT = [LPm.ATOT]; R.mc_LP_x2 = [LPm.x2]; R.mc_LP_type = {LPm.type};
    R.mc_LP_a    = [LPm.a];
    R.mc_BP_ATOT = [BPm.ATOT]; R.mc_BP_x2 = [BPm.x2];
    R.mc_H_ATOT  = HBm;   % Hopf or neutral saddle (MATCONT reports both as 'H')
    % |a| from MATCONT vs |a| computed at the exact fold (same normalisation)
    R.mc_LP_a_reldiff = NaN(1, numel(LPm));
    for i = 1:numel(LPm)
        if isempty(R.exact_LP_ATOT), break; end
        [d, j] = min(abs(R.exact_LP_ATOT - LPm(i).ATOT));
        if d <= 1e-3*max(1, LPm(i).ATOT)
            R.mc_LP_a_reldiff(i) = abs(abs(LPm(i).a) - abs(R.exact_LP_a(j)))/abs(R.exact_LP_a(j));
        end
    end
end
end

%% ===================================================================================
function R = empty_result()
R = struct('class_version',3,'t_exact',NaN,'cpu_exact',NaN,'t_nf',NaN,'t_mc',0,'cpu_mc',0,'exact_ok',false,'exact_class',NaN,'exact_label','', ...
    'exact_class_v2',NaN,'exact_label_v2','','exact_class_old7',NaN, ...
    'lp_coincides_bp',false,'detached_in_window',false, ...
    'multistationary',false,'bistable',false,'max_ss',NaN,'has_detached',false, ...
    'detached_fold_ATOT',Inf,'degenerate',false, ...
    'exact_LP_ATOT',[],'exact_LP_x2',[],'exact_LP_type',{{}}, ...
    'exact_BP_ATOT',[],'exact_BP_x2',[],'bistable_ATOT',zeros(0,2),'intervals',zeros(0,4), ...
    'n_LP_main',0,'n_LP_detached',0,'mushroom_candidate',false, ...
    'exact_LP_a',[],'exact_LP_beta',[],'exact_LP_kappa',[],'exact_LP_kappa_rel',[], ...
    'exact_BP_a',[],'exact_BP_b',[],'exact_BP_c',[],'exact_BP_disc',[], ...
    'exact_BP_slopes',zeros(0,2),'exact_BP_pitchfork',false(1,0), ...
    'curves',struct('type',{},'ATOT',{},'x2',{},'stable',{}), ...
    'mc_ran',false,'mc_ok',false,'mc_class',NaN,'mc_label','','mc_class_v2',NaN,'mc_label_v2','', ...
    'mc_LP_ATOT',[],'mc_LP_x2',[],'mc_LP_type',{{}},'mc_LP_a',[],'mc_LP_a_reldiff',[], ...
    'mc_BP_ATOT',[],'mc_BP_x2',[],'mc_H_ATOT',[], ...
    'mc_ATOT_reached',[],'mc_negative',false,'mc_failed_runs',0, ...
    'mc_curves',struct('type',{},'ATOT',{},'x2',{}),'message','');
end

function [x, v, s] = quiet_cont(x0, opt)
% cont() without the console chatter (MATCONT prints every singular point)
x = []; v = []; s = []; %#ok<NASGU>
[~] = evalc('[x, v, s] = cont(@equilibrium, x0, [], opt);');
end

function a = seed_ATOT(b, S, cfg)
% a point well inside each branch, where the branch has a regular (non-fold) state
a = [];
G = b.ATOT;
if isempty(G), return; end
switch b.type
    case {'main','sym'}
        a = cfg.ATOT_seed;                          % the original starting value (0.1)
        if a <= min(G) || a >= max(G), a = median(G); end
    case {'detached','offsym'}
        amin = min(G);                              % the fold (or the BP for 'offsym')
        if amin > cfg.ATOT_window(2), return; end   % branch lies beyond the window
        a = min(1.5*amin, 0.5*(amin + cfg.ATOT_window(2)));
end
% move off any fold or BP by a small relative margin
crit = [[S.folds_all.ATOT], [S.bp_all.ATOT]];
if ~isempty(crit) && min(abs(crit - a)) < 1e-4*a, a = a*1.01; end
end

function m = physical_mask(X, ATOT, P0, cfg)
tol = cfg.neg_tol*max([max(abs(ATOT)), P0(2), P0(3)]);
A = ATOT - sum(X,1); K = P0(2) - X(3,:) - X(4,:); Pf = P0(3) - X(5,:) - X(6,:);
m = ATOT >= 0 & all(X >= -tol,1) & A >= -tol & K >= -tol & Pf >= -tol;
end

function L = add_unique(L, a, x2, t, aMC)
if ~isempty(L) && any(abs([L.ATOT] - a) <= 1e-4*max(1,abs(a))), return; end
L(end+1) = struct('ATOT', a, 'x2', x2, 'type', t, 'a', aMC);
end
