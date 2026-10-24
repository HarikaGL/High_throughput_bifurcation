% =================== SingleStart_v3.m ================================================
% Single-start MATCONT labels ("MATCONT without the inverse method") for EXACTLY the
% 40,000 parameter sets of the ten-class ensemble, so that every comparison between
% single-start continuation and the exact labels refers to the same parameter sets.
%
% Procedure per parameter set: support/old_analyze_sample.m, the port of the original
% single-start pipeline (Ensemble_LHS.m) that is also used by the MATCONT check of
% Ensemble_LHS_v3.m:
%   ode15s from x = 0 at ATOT = 0.1 up to t = 5000 (RelTol 1e-8) -> one steady state
%   -> one forward MATCONT continuation in ATOT (1,500 points, settings of Table 11)
%   -> at every branch point: branch switching, continuation in both directions
%   -> 7 classes (1 LP S, 2 LP Z, 3 BP only, 4 LP+BP LP>BP, 5 LP+BP LP<BP,
%      6 biphasic, 7 other / no bifurcation); a state < 0 anywhere flags NEGATIVE.
%
% The parameter sets are READ from Ensemble_iter<i>_v3_1e6/MAT_Files/Parameter_Samples.mat
% (written by Ensemble_LHS_v3.m), not re-drawn, so they are identical by construction.
%
% Run from the Ensemble_v3 folder (inside MatCont7p6), like Ensemble_LHS_v3.m.
% Needs MATCONT 7p6, the Signal Processing Toolbox (findpeaks) and, for speed, the
% Parallel Computing Toolbox. Cost: about 5 s per parameter set, i.e. ~15 core-hours
% per iteration (~15 min per iteration on 64 workers).
%
% Output per iteration, in Ensemble_iter<i>_v3_1e6/SingleStart/:
%   ML_Ready_Dataset_singlestart.csv  14 parameters + single-start label and singular
%                                     points, in the column format of the original
%                                     ML_Ready_Dataset.csv (read by make_comparison_figures.m);
%                                     NEGATIVE / FAILED runs have Class_Number_7 = 9
%   SingleStart_vs_Exact.csv          single-start class against Exact_Class_v3
%   SingleStart_summary.txt           counts, agreement and timings: core-hours (sum of
%                                     per-set times), CPU time and wall clock, in the
%                                     same form as Core_Hours_iter<i>.txt of Ensemble_LHS_v3.m
%   (per-set times are also in the CSV: Run_seconds, CPU_seconds)
% ===================================================================================
clear; close all; clc;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
SUPPORT_DIR = fullfile(HERE, 'support'); addpath(SUPPORT_DIR); cd(HERE);
MATCONT_DIR = fileparts(HERE);
if ~exist(fullfile(MATCONT_DIR, 'init.m'), 'file')
    error('MATCONT not found in %s: put this folder inside MatCont7p6.', MATCONT_DIR);
end

%% ---------------- USER SETTINGS ----------------------------------------------------
ITERATIONS = 1:4;      % iterations to run
NumWorkers = 64;       % 0 = serial
nTest      = 0;        % > 0: run only the first nTest sets of each iteration (quick test)

%% ---------------- MATCONT and options of the original single-start pipeline -------
mc_init(MATCONT_DIR, SUPPORT_DIR);
oldOpt = contset;
oldOpt = contset(oldOpt,'InitStepSize',1);
oldOpt = contset(oldOpt,'MaxStepSize',1e-3);
oldOpt = contset(oldOpt,'MinStepSize',1e-5);
oldOpt = contset(oldOpt,'MaxNewtoniters',10);
oldOpt = contset(oldOpt,'MaxCorriters',15);
oldOpt = contset(oldOpt,'MaxTestiters',10);
oldOpt = contset(oldOpt,'VarTolerance',1e-7);
oldOpt = contset(oldOpt,'FunTolerance',1e-7);
oldOpt = contset(oldOpt,'TestTolerance',1e-6);
oldOpt = contset(oldOpt,'TestFunctions',1);
oldOpt = contset(oldOpt,'Locators',[1 1 1]);
oldOpt = contset(oldOpt,'Singularities',1);
oldOpt = contset(oldOpt,'Adapt',5);
oldOpt = contset(oldOpt,'MaxNumPoints',1500);
oldOpt = contset(oldOpt,'CheckClosed',50);
oldOpt = contset(oldOpt,'Backward',0);

useParallel = NumWorkers > 0 && license('test','Distrib_Computing_Toolbox');
if useParallel
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool('local', NumWorkers); end
    nW = pool.NumWorkers;
    wait(parfevalOnAll(pool, @cd, 0, HERE));
    wait(parfevalOnAll(pool, @mc_init, 0, MATCONT_DIR, SUPPORT_DIR));
    wait(parfevalOnAll(pool, @addpath, 0, SUPPORT_DIR));
else
    nW = 1;
end
[~, cpuName] = system('lscpu 2>/dev/null | grep "Model name" | cut -d: -f2'); cpuName = strtrim(cpuName);

names7 = {'LP_ONLY_S-shaped','LP_ONLY_Z-shaped','BP_ONLY','LP_BP_LPgreaterBP', ...
          'LP_BP_LPlessBP','NO_SING_biphasic','OTHER'};
COLS = {'KTOT','PTOT','k1','k11','k2','k3','k33','k4','k5','k55','k6','k7','k77','k8'};

for it = ITERATIONS
    iterDir = fullfile(HERE, sprintf('Ensemble_iter%d_v3_1e6', it));
    S = load(fullfile(iterDir, 'MAT_Files', 'Parameter_Samples.mat'), 'paramSamples', 'param_names');
    P = S.paramSamples; pn = S.param_names;
    n = size(P, 1); if nTest > 0, n = min(n, nTest); end
    fprintf('\n=== Iteration %d: single-start MATCONT on %d parameter sets ===\n', it, n);

    cls = zeros(n,1); neg = false(n,1); ok = false(n,1); reach = nan(n,1); tsec = nan(n,1); tcpu = nan(n,1);
    lpA = cell(n,1); lpX = cell(n,1); bpA = cell(n,1); bpX = cell(n,1); msg = cell(n,1);
    tWall = tic;
    if useParallel
        parfor k = 1:n
            [cls(k), neg(k), ok(k), reach(k), lpA{k}, lpX{k}, bpA{k}, bpX{k}, msg{k}, tsec(k), tcpu(k)] = ...
                run_one(P(k,:), oldOpt);
        end
    else
        for k = 1:n
            [cls(k), neg(k), ok(k), reach(k), lpA{k}, lpX{k}, bpA{k}, bpX{k}, msg{k}, tsec(k), tcpu(k)] = ...
                run_one(P(k,:), oldOpt);
            if mod(k,100) == 0, fprintf('  %d / %d\n', k, n); end
        end
    end
    wall = toc(tWall);

    %% ---- table in the format of the original ML_Ready_Dataset.csv ------------------
    T = array2table(P(1:n, cellfun(@(c) find(strcmp(pn, c)), COLS)), 'VariableNames', COLS);
    num2s = @(v) strjoin(arrayfun(@(a) sprintf('%.6f', a), v, 'UniformOutput', false), ' ');
    T.LP_Concentrations = strings(n,1); T.BP_Concentrations = strings(n,1);
    T.ATOT_LP = strings(n,1); T.ATOT_BP = strings(n,1);
    T.ATOT_LPBP_LP = repmat("N/A", n, 1); T.ATOT_LPBP_BP = repmat("N/A", n, 1);
    T.ATOT_NO_SING = repmat("N/A", n, 1);
    T.Multiclass_Label_7 = strings(n,1); T.Class_Number_7 = zeros(n,1);
    T.ATOT_reached = reach; T.Run_seconds = tsec; T.CPU_seconds = tcpu; T.Message = string(msg);
    for k = 1:n
        if isempty(lpA{k}), T.ATOT_LP(k) = "NO_SING"; T.LP_Concentrations(k) = "NO_SING";
        else, T.ATOT_LP(k) = num2s(lpA{k}); T.LP_Concentrations(k) = num2s(lpX{k}); end
        if isempty(bpA{k}), T.ATOT_BP(k) = "NONE"; T.BP_Concentrations(k) = "NONE";
        else, T.ATOT_BP(k) = num2s(bpA{k}); T.BP_Concentrations(k) = num2s(bpX{k}); end
        if ~ok(k) || neg(k)
            T.Class_Number_7(k) = 9;
            if neg(k), T.Multiclass_Label_7(k) = "NEGATIVE"; else, T.Multiclass_Label_7(k) = "FAILED"; end
            continue
        end
        T.Class_Number_7(k) = cls(k); T.Multiclass_Label_7(k) = names7{cls(k)};
        if any(cls(k) == [4 5])
            T.ATOT_LPBP_LP(k) = num2s(lpA{k}); T.ATOT_LPBP_BP(k) = num2s(bpA{k});
        end
        if cls(k) == 7 && isempty(lpA{k}) && isempty(bpA{k})
            T.ATOT_NO_SING(k) = sprintf('%.6f', reach(k));   % no singularity: range reached
        end
    end
    outDir = fullfile(iterDir, 'SingleStart'); if ~exist(outDir, 'dir'), mkdir(outDir); end
    writetable(T, fullfile(outDir, 'ML_Ready_Dataset_singlestart.csv'));

    %% ---- comparison with the exact labels ------------------------------------------
    E = readtable(fullfile(iterDir, 'Excel_Files', sprintf('ML_Ready_Dataset_v3_iter%d.csv', it)));
    ex = E.Exact_Class_v3(1:n);
    C = crosstab(categorical(ex), categorical(T.Class_Number_7));
    [ue, ~] = unique(ex); [us, ~] = unique(T.Class_Number_7);
    Ct = array2table(C, 'RowNames', compose('Exact_%d', ue), 'VariableNames', compose('Single_%d', us));
    writetable(Ct, fullfile(outDir, 'SingleStart_vs_Exact.csv'), 'WriteRowNames', true);

    fid = fopen(fullfile(outDir, 'SingleStart_summary.txt'), 'w');
    for f = [1 fid]
        fprintf(f, 'Single-start MATCONT on the v3 parameter sets, iteration %d (%s)\n', it, datestr(now));
        fprintf(f, 'parameter sets: %d | NEGATIVE: %d | FAILED: %d | used: %d\n', n, sum(neg), ...
            sum(~ok), sum(T.Class_Number_7 ~= 9));
        for c = 1:7, fprintf(f, '  class %d %-18s %6d\n', c, names7{c}, sum(T.Class_Number_7 == c)); end
        used = T.Class_Number_7 ~= 9;
        msS = any(T.Class_Number_7(used) == [1 2 3 4 5], 2);
        msE = ex(used) <= 8;
        fprintf(f, 'multistationary: single start %.1f%% | exact %.1f%% (same %d sets)\n', ...
            100*mean(msS), 100*mean(msE), sum(used));
        fprintf(f, 'ATOT reached (median, 5-95%%): %.0f (%.0f-%.0f)\n', median(reach(used)), ...
            prctile(reach(used), 5), prctile(reach(used), 95));
        fprintf(f, '\nCORE-HOURS (sum of per-set times = serial CPU time)\n');
        fprintf(f, '  single-start labelling, all %d sets : %.4f core-h (%.3f s/set mean, %.3f median, 5-95%% %.2f-%.2f)\n', ...
            n, sum(tsec,'omitnan')/3600, mean(tsec,'omitnan'), median(tsec,'omitnan'), prctile(tsec,5), prctile(tsec,95));
        fprintf(f, '  (CPU time reported by cputime: %.4f core-h)\n', sum(tcpu,'omitnan')/3600);
        fprintf(f, 'WALL CLOCK: %.1f s = %.3f h on %d workers -> allocated %.4f core-h\n', wall, wall/3600, nW, wall/3600*nW);
        fprintf(f, 'Machine: %s | CPU: %s | MATLAB %s\n', computer, cpuName, version);
    end
    fclose(fid);
    save(fullfile(outDir, 'SingleStart_results.mat'), 'T', 'cls', 'neg', 'ok', 'reach', ...
        'lpA', 'lpX', 'bpA', 'bpX', 'tsec', 'tcpu', 'wall', 'nW', 'oldOpt', '-v7.3');
    fprintf('saved %s\n', outDir);
end

%% =====================================================================================
function [c, neg, ok, reach, lpA, lpX, bpA, bpX, msg, t, tc] = run_one(P0, oldOpt)
t0 = tic; c0 = cputime;
c = 0; neg = false; ok = false; reach = NaN; lpA = []; lpX = []; bpA = []; bpX = []; msg = '';
try
    R = old_analyze_sample(P0, oldOpt, 50);
    ok = R.ok; neg = R.negative; c = R.class7; reach = R.ATOT_reached;
    lpA = R.LP_ATOT; lpX = R.LP_x2; bpA = R.BP_ATOT; bpX = R.BP_x2; msg = R.message;
catch ME
    msg = ME.message;
end
t = toc(t0); tc = cputime - c0;
end
