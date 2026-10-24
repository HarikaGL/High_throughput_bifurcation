% ======================= TR_03_two_class_sweeps.m ======================================
% Two-class switches: one quantity is swept back and forth through the boundary between
% two neighbouring classes, in nPer sets of the source class, all other parameters fixed;
% the exact class is computed at nGrid points of the sweep.
%   clean   = only the two classes of the pair occur along the sweep (both of them; for the
%             7 <-> 8 sweep Class 6 is allowed exactly at the boundary KM1 = KM5)
%   switch  = the second class occurs somewhere along the sweep
%
% Swept quantities (value v on the grid linspace(-1, 1, nGrid)):
%   'x'   : x = log10(k2 KTOT/(k8 PTOT)) set to v by changing one parameter (k8, k2, KTOT or PTOT)
%   'y'   : y = log10(k4 KTOT/(k6 PTOT)) set to v by changing one parameter (k4, k6, KTOT or PTOT)
%   'bp'  : log10(k4/k2) = v with k4 = k8 (keeps the symmetry of Iterations 3-4)
%   'km'  : log10(KM1/KM5) = v by changing k5
%   'mul' : one parameter multiplied by 10^(2v) (fold of the main branch: k3 or k7;
%           k3 and k7 enter neither x nor y, so the detached-branch condition is fixed)
%
% Output (Class_transitions/results/): TR_sweeps.csv (spec, set, step, value, class),
%   TR_sweep_summary.csv / .txt, TR_two_class_examples.mat (for TR_04_figures.m)
% ======================================================================================
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
addpath(fullfile(fileparts(HERE), 'support'), fullfile(HERE, 'tr_lib'));
resDir = fullfile(HERE, 'results'); if ~exist(resDir, 'dir'), mkdir(resDir); end
if ~exist('nPer', 'var'), nPer = 80; end
if ~exist('nGrid', 'var'), nGrid = 41; end
if ~exist('specList', 'var'), specList = []; end   % run a subset of the specs, e.g. 1:9
g = linspace(-1, 1, nGrid);
% spec: id, iteration, source class, class B, type, parameter index, shown in figure (1/0)
%       parameter: 14 k8, 5 k2, 1 KTOT, 2 PTOT, 8 k4, 11 k6, 6 k3, 12 k7, 9 k5
SPEC = {
  1, 1,  1,  2, 'x',  14, 1
  2, 1, 10,  5, 'x',  14, 1
  3, 2,  3,  4, 'y',   8, 1
  4, 2,  9,  5, 'y',   8, 1
  5, 1,  1, 10, 'mul', 6, 1
  6, 1,  2,  5, 'mul', 6, 1
  7, 2,  3,  9, 'mul', 6, 1
  8, 2,  4,  5, 'mul', 6, 1
  9, 3, 10,  6, 'bp',  8, 1
 10, 4,  7,  8, 'km',  9, 1
 11, 1,  1,  2, 'x',   5, 0
 12, 1,  1,  2, 'x',   1, 0
 13, 1,  1,  2, 'x',   2, 0
 14, 2,  3,  4, 'y',  11, 0
 15, 2,  3,  4, 'y',   1, 0
 16, 2,  3,  4, 'y',   2, 0
 17, 1,  1, 10, 'mul',12, 0};
pn = {'KTOT','PTOT','k1','k11','k2','k3','k33','k4','k5','k55','k6','k7','k77','k8'};
if isempty(specList), specList = 1:size(SPEC, 1); end

E = tr_load_ensemble(HERE);
for sid = specList
    [id, it, A, B, typ, pj, shown] = SPEC{sid, :}; %#ok<ASGLU>
    f = fullfile(resDir, sprintf('TR_sweep_spec%02d.mat', id));
    if exist(f, 'file'), continue; end
    t0 = tic;
    pool = find(E.use & E.iteration == it & E.class == A);
    pick = pool(sort(tr_pick(numel(pool), 120, 1300 + 10*it + A)));   % the sets of TR_02_scan
    pick = pick(1:min(nPer, numel(pick)));
    CL = zeros(numel(pick), nGrid); VAL = zeros(numel(pick), nGrid);
    parfor s = 1:numel(pick)
        p0 = E.p(pick(s), :); cl = zeros(1, nGrid); val = zeros(1, nGrid);
        for k = 1:nGrid
            [p, val(k)] = tr_sweep_point(p0, typ, pj, g(k));
            cl(k) = tr_label(p);
        end
        CL(s, :) = cl; VAL(s, :) = val;
    end
    save(f, 'id', 'it', 'A', 'B', 'typ', 'pj', 'shown', 'pick', 'CL', 'VAL', 'g', '-v7');
    fprintf('spec %2d (%d <-> %d, %s via %s): %.0f s\n', id, A, B, typ, pn{pj}, toc(t0));
end

%% ---------------- summary ---------------------------------------------------------------
have = arrayfun(@(i) exist(fullfile(resDir, sprintf('TR_sweep_spec%02d.mat', i)), 'file') == 2, 1:size(SPEC,1));
if all(have)
    SUM = []; ROWS = []; EXM = struct();
    fid = fopen(fullfile(resDir, 'TR_sweep_summary.txt'), 'w');
    for f1 = [1 fid]
        fprintf(f1, 'Two-class sweeps: %d grid points, up to %d sets per pair\n', nGrid, nPer);
        fprintf(f1, 'spec pair     quantity  via    sets  clean   switch  other classes met\n');
    end
    for sid = 1:size(SPEC, 1)
        L = load(fullfile(resDir, sprintf('TR_sweep_spec%02d.mat', sid)));
        n = size(L.CL, 1); clean = false(n, 1); sw = false(n, 1); other = [];
        for s = 1:n
            u = unique(L.CL(s, :));
            sw(s) = any(u == L.B);
            ok = [L.A L.B];
            if strcmp(L.typ, 'km'), ok = [ok 6]; end   % Class 6 exactly at KM1 = KM5 (the boundary)
            clean(s) = sw(s) && any(u == L.A) && all(ismember(u, ok));
            other = [other, setdiff(u, ok)]; %#ok<AGROW>
        end
        uo = unique(other);
        os = strjoin(arrayfun(@(v) sprintf('%d', v), uo, 'UniformOutput', false), ',');
        for f1 = [1 fid]
            fprintf(f1, '%3d  %2d <-> %-2d %-8s  %-5s  %4d  %5.1f%%  %5.1f%%  %s\n', sid, L.A, L.B, L.typ, pn{L.pj}, n, ...
                100*mean(clean), 100*mean(sw), os);
        end
        SUM = [SUM; sid, L.it, L.A, L.B, L.pj, L.shown, n, mean(clean), mean(sw)]; %#ok<AGROW>
        ROWS = [ROWS; repmat(sid, numel(L.CL), 1), kron(L.pick(:), ones(nGrid, 1)), ...
                repmat((1:nGrid)', n, 1), reshape(L.VAL', [], 1), reshape(L.CL', [], 1)]; %#ok<AGROW>
        % example for the figure: the clean set whose switch lies nearest the middle of the sweep
        if L.shown
            best = Inf; bs = 0;
            for s = find(clean)'
                jb = find(L.CL(s, :) ~= L.CL(s, 1), 1); if isempty(jb), continue; end
                if abs(jb - (nGrid+1)/2) < best, best = abs(jb - (nGrid+1)/2); bs = s; end
            end
            if bs == 0, bs = find(sw, 1); end
            cl = L.CL(bs, :); val = L.VAL(bs, :);
            ka = find(cl == L.A); kb = find(cl == L.B);
            % one value inside each class, about a quarter of the sweep away from the boundary
            [~, ia] = min(abs(ka - mean(ka))); [~, ib] = min(abs(kb - mean(kb)));
            pA = tr_sweep_point(E.p(L.pick(bs), :), L.typ, L.pj, g(ka(ia)));
            pB = tr_sweep_point(E.p(L.pick(bs), :), L.typ, L.pj, g(kb(ib)));
            EXM(sid).A = L.A; EXM(sid).B = L.B; EXM(sid).typ = L.typ; EXM(sid).pj = L.pj;
            EXM(sid).pA = pA; EXM(sid).pB = pB; EXM(sid).vA = val(ka(ia)); EXM(sid).vB = val(kb(ib));
            EXM(sid).strip = [val(:), cl(:)]; EXM(sid).clean = mean(clean); EXM(sid).set = L.pick(bs);
        end
    end
    fclose(fid);
    tr_write_csv(fullfile(resDir, 'TR_sweep_summary.csv'), 'spec,iter,A,B,param,shown,sets,clean,switch', SUM);
    tr_write_csv(fullfile(resDir, 'TR_sweeps.csv'), 'spec,set,step,value,class', ROWS);
    save(fullfile(resDir, 'TR_two_class_examples.mat'), 'EXM', '-v7');
end
