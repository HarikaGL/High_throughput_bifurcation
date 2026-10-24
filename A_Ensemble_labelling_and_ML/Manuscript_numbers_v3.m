% ======================= Manuscript_numbers_v3.m ===================================
% Reproduces the numbers of the manuscript that are computed from the output files
% of Ensemble_LHS_v3.m (no new analysis; reads CSV files only):
%   * class counts per iteration and multistationary fractions   (Tables 1 and 2)
%   * detached-branch criterion k8*PTOT > k2*KTOT & k4*KTOT > k6*PTOT against the
%     computed detached fold, kinase-/phosphatase-dominant quadrants (Results)
%   * dependence of the multistationary fraction on the ATOT window      (Table 3)
%   * branch point <=> k4 > k2 (Iterations 3, 4); Class 7 <=> K_M5 < K_M1 (Iter. 4)
%   * normal-form coefficients per singularity type, agreement with MATCONT
%                                                                      (Table 6)
% Run from the Ensemble_v3 folder after Ensemble_LHS_v3.m has been run for
% ITERATION = 1..4. Needs base MATLAB only. Output: Manuscript_numbers_v3.txt
% ===================================================================================
clear; clc;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
cd(HERE);
out = fullfile(HERE, 'Manuscript_numbers_v3.txt');
fid = fopen(out, 'w'); P = @(varargin) fprintfBoth(fid, varargin{:});

D = cell(1,4); NF = cell(1,4);
for it = 1:4
    dd = fullfile(HERE, sprintf('Ensemble_iter%d_v3_1e6', it), 'Excel_Files');
    D{it}  = readtable(fullfile(dd, sprintf('ML_Ready_Dataset_v3_iter%d.csv', it)), 'TextType','string');
    NF{it} = readtable(fullfile(dd, sprintf('Normal_Forms_v3_iter%d.csv', it)), 'TextType','string');
end

%% 1. class counts and multistationarity --------------------------------------------
P('1. CLASS COUNTS (Exact_Class_v3) and multistationary fraction\n');
tot = zeros(1,11);
for it = 1:4
    c = D{it}.Exact_Class_v3; n = histc(c, 1:11); tot = tot + n(:)';
    ms = mean(c >= 1 & c <= 8);
    P('  Iteration %d (n = %d): %s | multistationary (Classes 1-8) %.2f%% | flagged %d\n', ...
        it, numel(c), mat2str(n(1:10)'), 100*ms, sum(D{it}.LP_Coincides_BP));
end
P('  Total: %s | multistationary %.2f%%\n\n', mat2str(tot(1:10)), 100*sum(tot(1:8))/sum(tot));

%% 2. detached-branch criterion -------------------------------------------------------
P('2. DETACHED-BRANCH CRITERION (Iterations 1 and 2)\n');
for it = 1:2
    T = D{it};
    crit = T.k8.*T.PTOT > T.k2.*T.KTOT & T.k4.*T.KTOT > T.k6.*T.PTOT;
    det  = T.Has_Detached_Branch == 1;          % detached fold inside the window
    kin  = T.k8.*T.PTOT < T.k2.*T.KTOT & T.k4.*T.KTOT > T.k6.*T.PTOT;
    pho  = T.k8.*T.PTOT > T.k2.*T.KTOT & T.k4.*T.KTOT < T.k6.*T.PTOT;
    P('  Iteration %d: criterion true %d | detached fold in window %d | in window but criterion false %d | criterion true, fold above window %d\n', ...
        it, sum(crit), sum(det), sum(det & ~crit), sum(crit & ~det));
    P('     classes of the criterion-true sets: %s\n', classCounts(T.Exact_Class_v3(crit)));
    P('     kinase dominates both sites: %d sets, classes %s\n', sum(kin), classCounts(T.Exact_Class_v3(kin)));
    P('     phosphatase dominates both sites: %d sets, classes %s\n', sum(pho), classCounts(T.Exact_Class_v3(pho)));
end
P('\n');

%% 3. multistationary fraction against the ATOT window -----------------------------
P('3. MULTISTATIONARY FRACTION vs ATOT WINDOW (first singularity below the bound)\n');
for it = 1:4
    T = D{it}; n = height(T); m = inf(n,1);
    for k = 1:n
        v = [parseVec(T.LP_ATOT_values(k)), parseVec(T.BP_ATOT_values(k))];
        if ~isempty(v), m(k) = min(v); end
    end
    crit = T.k8.*T.PTOT > T.k2.*T.KTOT & T.k4.*T.KTOT > T.k6.*T.PTOT;
    bp = false(n,1); if it >= 3, bp = T.k4 > T.k2; end
    allA = isfinite(m) | crit | bp;
    P('  Iteration %d: <=1e2 %.1f%% | <=1e3 %.1f%% | <=1e4 %.1f%% | <=1e6 %.1f%% | all ATOT (lower bound) %.1f%%\n', ...
        it, 100*mean(m <= 1e2), 100*mean(m <= 1e3), 100*mean(m <= 1e4), 100*mean(m <= 1e6), 100*mean(allA));
end
P('\n');

%% 4. branch points and the orientation criterion ------------------------------------
P('4. BRANCH POINTS (Iterations 3, 4) and Class 7/8 criterion (Iteration 4)\n');
for it = 3:4
    T = D{it}; hasBP = ismember(T.Exact_Class_v3, [6 7 8]);
    P('  Iteration %d: k4 > k2 in %d sets | branch point in window %d | k4 > k2 without BP in window %d | BP with k4 <= k2 %d\n', ...
        it, sum(T.k4 > T.k2), sum(hasBP), sum(T.k4 > T.k2 & ~hasBP), sum(hasBP & T.k4 <= T.k2));
end
T = D{4}; KM1 = (T.k11 + T.k2)./T.k1; KM5 = (T.k55 + T.k6)./T.k5;
c7 = T.Exact_Class_v3 == 7; c8 = T.Exact_Class_v3 == 8;
P('  Iteration 4: Class 7 with KM5 < KM1: %d of %d | Class 8 with KM5 > KM1: %d of %d\n\n', ...
    sum(c7 & KM5 < KM1), sum(c7), sum(c8 & KM5 > KM1), sum(c8));

%% 5. normal-form coefficients ------------------------------------------------------
P('5. NORMAL-FORM COEFFICIENTS (|a|, singularities inside the window)\n');
A = []; for it = 1:4, N = NF{it}; N.Iter = repmat(it, height(N), 1); A = [A; N]; end %#ok<AGROW>
grp = { 'BP, Iteration 3',              A.Iter == 3 & A.Kind == "BP";
        'BP, Iteration 4',              A.Iter == 4 & A.Kind == "BP";
        'LP main branch, Iter. 1-2',    A.Iter <= 2 & A.Kind == "LP" & A.Branch == "main";
        'LP detached branch, Iter. 1-2',A.Iter <= 2 & A.Kind == "LP" & A.Branch == "detached";
        'LP crossing branch, Iter. 4',  A.Iter == 4 & A.Kind == "LP";
        'LP flagged, Iter. 3',          A.Iter == 3 & A.Kind == "LP"};
for g = 1:size(grp,1)
    a = abs(A.a(grp{g,2}));
    P('  %-30s n = %5d | |a| min %.2g  median %.2g  max %.2g\n', grp{g,1}, numel(a), min(a), median(a), max(a));
end
b4 = A.Iter == 4 & A.Kind == "BP"; s = A.slope1(b4).*A.slope2(b4); cl = A.Exact_Class_v3(b4);
P('  Iteration 4 BP slopes opposite sign: Class 7 %d of %d | same sign: Class 8 %d of %d\n', ...
    sum(s < 0 & cl == 7), sum(cl == 7), sum(s > 0 & cl == 8), sum(cl == 8));
m = ~isnan(A.a_MATCONT);
r = abs(abs(A.a(m)) - abs(A.a_MATCONT(m)))./abs(A.a_MATCONT(m));
P('  agreement with MATCONT (folds of the validation subsample): n = %d | median rel. diff %.2g | 95th percentile %.2g\n', ...
    sum(m), median(r), prctileSimple(r, 95));
P('  all singularities in the window: %d\n', height(A));
fclose(fid);
fprintf('\nWritten to %s\n', out);

%% ------------------------------------------------------------------------------------
function fprintfBoth(fid, varargin)
fprintf(varargin{:}); fprintf(fid, varargin{:});
end
function v = parseVec(s)
% '[40.97 40.56]' or 'zeros(1,0)' or a number -> numeric row vector
if ismissing(s), v = []; return; end
s = char(s); if contains(s, 'zeros'), v = []; return; end
v = sscanf(strrep(strrep(s, '[', ''), ']', ''), '%f')';
end
function s = classCounts(c)
u = unique(c(:))'; s = strjoin(arrayfun(@(x) sprintf('%d:%d', x, sum(c == x)), u, 'UniformOutput', false), ', ');
end
function p = prctileSimple(x, q)
x = sort(x(:)); p = x(max(1, ceil(q/100*numel(x))));
end
