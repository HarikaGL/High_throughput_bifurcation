%% ============================================================
%  ML ANALYSIS, TEN BIFURCATION CLASSES (v3)   -- PARALLEL VERSION --
%                                               ML_v3_10class_parallel.m
%  Same analyses, same splits, same settings and outputs as ML_v3_10class.m, but every
%  model fit is an independent job that runs on a pool of cfg.nWorkers MATLAB workers
%  (default 64). Each worker imports Python once and fits ONE model with ONE thread, so
%  64 fits run at the same time without oversubscribing the cores. Jobs are sent
%  longest-first (Gradient Boosting, CatBoost, ...) for good load balancing.
%  Needs: Parallel Computing Toolbox, Python with numpy / lightgbm / scikit-learn
%  (+ xgboost / catboost), and support/ml_v3_core.m (model core shared by the workers).
%  Memory: about 1.5-2 GB per worker (MATLAB + Python); 64 workers ~ 100-130 GB.
%
%  Same analyses as ML_RC_based2.m, now on the
%  ten-class labels of Ensemble_LHS_v3:
%     1 LP_S_ONLY          2 LP_S_DETACHED      3 LP_Z_ONLY
%     4 LP_Z_DETACHED      5 LP_DETACHED_ONLY   6 BP_ONLY
%     7 LP_BP_LPgreaterBP  8 LP_BP_LPlessBP     9 NO_SING_biphasic
%    10 NO_SING_monotone
%
%  HOW TO RUN: put this file in the Ensemble_v3 folder (next to
%  Ensemble_LHS_v3.m) and press Run. It reads
%     Ensemble_iter<i>_v3_1e6/Excel_Files/ML_Ready_Dataset_v3_iter<i>.csv
%  (falls back to the .xlsx) and writes everything to ML_v3_results/.
%  Python (numpy, lightgbm, scikit-learn) is used through pyenv, as before.
%
%  Models (as in the original manuscript, all with inverse-frequency sample weights):
%     Random Forest, XGBoost, LightGBM, CatBoost, Gradient Boosting
%  plus Gaussian Naive Bayes and the most-frequent-class baseline. The section
%  [MODELS] reproduces the manuscript's model table (held-out accuracy, k-fold CV
%  accuracy / macro-F1 / weighted-F1) pooled and within every iteration.
%  ROC / AUC (one-vs-rest, every class) and convergence curves (train / test log-loss and
%  error vs boosting rounds or number of trees, overall and per class) for every model:
%  section [AUC] / [CONV], figures ROC_*, convergence_*, convergence_per_class_*.
%  Run time: sklearn GradientBoosting is single-threaded and slow on 10 classes; set
%  cfg.gbFast = true (HistGradientBoosting) or cfg.quick = true for a fast test run.
%
%  FEATURES: only the 14 sampled parameters KTOT PTOT k1 k11 k2 k3 k33 k4 k5 k55 k6 k7 k77 k8
%  (log10), and quantities computed from them alone (pairwise log-ratios, CLR). No label,
%  descriptor (Bistable, Has_Detached_Branch, ...), singularity value or Iteration column is
%  ever an input. SPLITS: stratified 70% train / 10% validation / 20% test; models are fitted
%  on the train rows only (column filtering and JSD selection too); the validation rows are
%  used only for the convergence curves; every reported accuracy / ROC is on the test rows.
%  Model comparison CV folds are drawn inside the train rows only.
%
%  Changes relative to ML_RC_based2.m
%   * labels: Exact_Class_v3 (ten classes); Exact_Class_v2 / Class_Number_7
%     can be added in cfg.labelCols for comparison
%   * rows removed from the ML (reported in the summary):
%       - Use_For_ML == 0: LP+BP sets whose LP coincides with the BP
%         (31 in Iteration 3, 1 in Iteration 4 of the runs)
%       - any class with fewer than cfg.minClassPerIter rows inside an
%         iteration (too few for iteration-wise learning)
%   * decisive tests: several {iteration, class A, class B, axes}; axes may
%     be products, e.g. {'k2*KTOT','k8*PTOT'} = log10(k2 KTOT / (k8 PTOT))
%       Iter 1: 1 vs 2  (S only vs S + detached)  detached-branch planes
%       Iter 2: 3 vs 4  (Z only vs Z + detached)  detached-branch planes
%       Iter 4: 7 vs 8  (LP > BP vs LP < BP)      log10(k1/k5), log10(k11/k55)
%   * confusion audit for several class pairs (main-branch vs detached)
%   * small multiples in the class colours of class_info_v3 and on the
%     detached-branch coordinates
%   * thicker lines / larger fonts in all figures
%
%  Analyses: ML-a baselines with every accuracy; ML-b confound
%  audit; ML-c pooled feature-set ablation (A-F); ML-d iteration-wise
%  (stratified) classification, importances and confusion matrices; ML-e
%  decisive tests with linear boundaries; ML-f pooled confusion matrix;
%  FIG small multiples; EMB embedding audit; EX class x iteration and
%  bistable / multistationary fractions.
% =============================================================

clear; clc; close all;

%% ============================================================
%% 0. CONFIG
%% ============================================================

HERE = fileparts(mfilename('fullpath')); if isempty(HERE); HERE = pwd; end
cd(HERE);
if exist(fullfile(HERE,'support'),'dir'); addpath(fullfile(HERE,'support')); end   % support/ also holds feature_colors.m and cm_graded.m (figure colours)
if exist('ml_v3_core', 'file') ~= 2
    error(['support/ml_v3_core.m not found.\nThis script must sit in the Ensemble_v3 folder and ' ...
           'Ensemble_v3/support/ must contain ml_v3_core.m (copy it there from the package).\n' ...
           'Script folder: %s'], HERE);
end

cfg = struct();
cfg.nWorkers   = 64;                        % parallel workers (0 = run serially on the client)
cfg.nThreads   = -1;                        % threads per model on the client (workers use 1)
cfg.dataDir    = HERE;                      % the Ensemble_v3 folder
cfg.iterIds    = [1 2 3 4];
cfg.iterFiles  = arrayfun(@(k) fullfile(sprintf('Ensemble_iter%d_v3_1e6',k), 'Excel_Files', ...
                          sprintf('ML_Ready_Dataset_v3_iter%d.csv',k)), cfg.iterIds, 'UniformOutput', false);
cfg.pooledFile = '';                        % OR one pooled file with an 'Iteration' column

cfg.labelCols  = {'Exact_Class_v3'};        % add 'Exact_Class_v2','Class_Number_7' to compare
cfg.nameCols   = containers.Map({'Exact_Class_v3','Exact_Class_v2','Class_Number_7'}, ...
                                {'Exact_Label_v3','Exact_Label_v2','Multiclass_Label_7'});
cfg.flagCols   = {'Bistable','Multistationary','Has_Detached_Branch','Detached_Fold_In_Window'};

% rows used for the ML
cfg.useForMLCol     = 'Use_For_ML';         % 0 = flagged LP = BP set -> removed
cfg.minClassPerIter = 50;                   % classes smaller than this inside an iteration are removed

cfg.featNames  = {'KTOT','PTOT','k1','k11','k2','k3','k33','k4','k5',...
                  'k55','k6','k7','k77','k8'};

% Model comparison (as in the original manuscript, Table "ml_performance": five weighted
% tree ensembles) + Naive Bayes and the most-frequent-class baseline. Held-out accuracy and
% stratified k-fold CV accuracy / macro-F1 / weighted-F1, pooled AND within every iteration.
cfg.compareModels = {'RandomForest','XGBoost','LightGBM','CatBoost','GradBoost','NaiveBayes','MostFrequent'};
cfg.compareFeatureSets = [1 4];             % A (all log-ratios) and D (JSD + broken-stick, original pipeline)
cfg.cvFolds       = 5;
cfg.compareScopes = 'pooled+iterations';    % 'pooled' or 'pooled+iterations'
% models used in the (larger) ablation / stratum / decisive loops
cfg.models       = {'LightGBM','RandomForest','XGBoost','CatBoost','NaiveBayes'};
cfg.primaryModel = 'LightGBM';              % used for importances / confusion / figures
cfg.useClassWeights = true;
cfg.lgbRounds  = 300;   cfg.lgbLeaves = 63;      % original pipeline: 300 rounds, 63 leaves, lr 0.1
cfg.rfTrees    = 300;                           % original pipeline used 1000 (slower, same accuracy)
cfg.xgbRounds  = 300;   cfg.xgbDepth  = 6;       % original pipeline
cfg.catIter    = 300;   cfg.catDepth  = 6;       % original pipeline
cfg.gbRounds   = 200;   cfg.gbDepth   = 5;       % original pipeline (sklearn GradientBoosting, slow)
cfg.runCurves     = true;              % ROC / AUC and convergence (learning) curves
cfg.curveScopes   = 'pooled+iterations'; % 'pooled' or 'pooled+iterations'
cfg.curvePoints   = 30;                 % checkpoints along the boosting rounds / trees
cfg.curveTrainMax = 5000;               % training rows used for the training curve (speed)
cfg.gbFast     = false;  % true: sklearn HistGradientBoosting instead (much faster, multi-threaded)

cfg.nRepeats   = 3;                         % random stratified splits per experiment
cfg.testFrac   = 0.20;                      % held-out TEST set (reported accuracies)
cfg.valFrac    = 0.10;                      % VALIDATION set (convergence curves only)
                                            % -> 70% TRAIN: the only rows any model is fitted on
cfg.seed       = 42;
cfg.nBins      = 30;                        % JSD histogram bins
cfg.tolEq      = 1e-9;                      % |log10 ratio| < tol  <=>  ratio == 1
cfg.missingPolicy = 'drop';                 % rows with a missing feature are dropped

% Decisive tests: {iteration, classA, classB, {num1,den1; num2,den2}}  (one row per test)
% axis = log10(num/den); num/den may be products of parameters, e.g. 'k2*KTOT'
detAxes = {'k2*KTOT','k8*PTOT'; 'k4*KTOT','k6*PTOT'};     % detached-branch planes
kmAxes  = {'k1','k5'; 'k11','k55'};                       % first-binding asymmetry
cfg.decisiveTests = {1, '1', '2', detAxes; ...
                     2, '3', '4', detAxes; ...
                     4, '7', '8', kmAxes};

% Mechanistic coordinates for the small multiples: {numerator, denominator}
cfg.smAxes = detAxes;

% Pairs of classes whose mutual confusion is audited (ML-f), one pair per row
cfg.confusionCheck = {'1','2'; '3','4'; '2','5'; '4','5'};

% figure style
cfg.LW = 2.0;  cfg.FS = 13;

% Embedding audit
cfg.runEmbeddingAudit = true;
cfg.embedMethod   = 'tsne';                 % 'tsne' or 'umap' (needs umap-learn)
cfg.embedN        = 6000;                   % global embedding size
cfg.embedPerRegime= 3000;                   % per-regime embedding size
cfg.knnK          = 15;

cfg.quick = false;                          % true = 1 repeat, fewer trees, small embeddings
if cfg.quick
    cfg.nRepeats = 1; cfg.lgbRounds = 100; cfg.rfTrees = 100; cfg.xgbRounds = 100;
    cfg.catIter = 100; cfg.gbRounds = 50; cfg.cvFolds = 3;
    cfg.embedN = 2500; cfg.embedPerRegime = 1500;
end

outRoot = fullfile(HERE, 'ML_v3_results_parallel');
if ~exist(outRoot,'dir'); mkdir(outRoot); end

%% ============================================================
%% 1. PYTHON MODULES
%% ============================================================

fprintf('\n--- Installing / verifying Python packages ---\n');
try
    pyExec = char(pyenv().Executable);
    pkgs = 'lightgbm scikit-learn xgboost catboost';
    if strcmpi(cfg.embedMethod,'umap'); pkgs = [pkgs ' umap-learn']; end
    system(sprintf('"%s" -m pip install %s --quiet', pyExec, pkgs));
catch
    warning('Could not install Python packages automatically.');
end

mods = ml_v3_core('mods');                   % client-side modules (workers import their own)
if isempty(mods.xgb), cfg.models = setdiff(cfg.models, {'XGBoost'}, 'stable'); cfg.compareModels = setdiff(cfg.compareModels, {'XGBoost'}, 'stable'); warning('xgboost not available: XGBoost skipped.'); end
if isempty(mods.cb),  cfg.models = setdiff(cfg.models, {'CatBoost'}, 'stable'); cfg.compareModels = setdiff(cfg.compareModels, {'CatBoost'}, 'stable'); warning('catboost not available: CatBoost skipped.'); end

%% ============================================================
%% 1b. PARALLEL POOL (process workers: Python needs separate processes)
%% ============================================================
cfg.useParallel = cfg.nWorkers > 0 && license('test', 'Distrib_Computing_Toolbox');
if cfg.useParallel
    pool = gcp('nocreate');
    if ~isempty(pool) && (~isa(pool, 'parallel.ProcessPool') || pool.NumWorkers ~= cfg.nWorkers)
        delete(pool); pool = [];
    end
    if isempty(pool), pool = parpool('Processes', cfg.nWorkers); end
    % one thread per model on every worker (set before Python is loaded there)
    for v = {'OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'}
        wait(parfevalOnAll(pool, @setenv, 0, v{1}, '1'));
    end
    wait(parfevalOnAll(pool, @cd, 0, HERE));
    wait(parfevalOnAll(pool, @addpath, 0, fullfile(HERE, 'support')));
    fprintf('Parallel pool: %d workers\n', pool.NumWorkers);
    % warm-up: import the Python modules on every worker once
    wait(parfevalOnAll(pool, @ml_v3_core, 0, 'mods'));
end

%% ============================================================
%% 2. LOAD + POOL ALL ITERATIONS
%% ============================================================

D = loadIterations(cfg);
fprintf('\nPooled dataset: %d rows from iterations: %s\n', size(D.X,1), ...
    strjoin(string(unique(D.Iter)'), ', '));

%% ============================================================
%% 3. RUN EVERYTHING FOR EACH LABEL SCHEME
%% ============================================================

allSummary = {};
for li = 1:numel(cfg.labelCols)
    lc = cfg.labelCols{li};
    if ~isfield(D.lab, lc)
        fprintf('Label column %s not found in the data - skipped.\n', lc);
        continue;
    end
    [~, S] = analyseLabelScheme(D, lc, cfg, mods, outRoot);
    allSummary = [allSummary, S, {''}]; %#ok<AGROW>
end

%% ============================================================
%% 4. SUMMARY FILE (numbers to quote in the response letter)
%% ============================================================

summaryFile = fullfile(outRoot, 'ML_v3_summary.txt');
fid = fopen(summaryFile, 'w');
for i = 1:numel(allSummary)
    fprintf(fid, '%s\n', allSummary{i});
    fprintf('%s\n', allSummary{i});
end
fclose(fid);
fprintf('\nSummary written to %s\nAll outputs under %s\n', summaryFile, outRoot);

%% ============================================================
%% ============================================================
%% LOCAL FUNCTIONS
%% ============================================================
%% ============================================================

function D = loadIterations(cfg)
% Read every iteration file, keep ONLY the 14 features + label/descriptor
% columns, tag each row with its iteration, and pool.
    fn = cfg.featNames;
    labCols  = cfg.labelCols;
    nameCols = {};
    for i = 1:numel(labCols)
        if isKey(cfg.nameCols, labCols{i}); nameCols{end+1} = cfg.nameCols(labCols{i}); end %#ok<AGROW>
    end
    flagNames = cfg.flagCols;

    if ~isempty(cfg.pooledFile)
        srcFiles = {cfg.pooledFile}; srcIter = NaN;
    else
        srcFiles = cfg.iterFiles;    srcIter = cfg.iterIds;
    end

    Xc = {}; Ic = {}; Uc = {};
    labC = struct(); nameC = struct(); flagC = struct();
    for i = 1:numel(labCols);  labC.(labCols{i})   = {}; end
    for i = 1:numel(nameCols); nameC.(nameCols{i}) = {}; end
    for i = 1:numel(flagNames);flagC.(flagNames{i})= {}; end

    for k = 1:numel(srcFiles)
        f = fullfile(cfg.dataDir, srcFiles{k});
        if ~isfile(f)
            [pp, nn] = fileparts(f); fx = fullfile(pp, [nn '.xlsx']);
            if isfile(fx); f = fx; else
                warning('File not found, skipped: %s', f);
                continue;
            end
        end
        T = readtable(f);
        vars = T.Properties.VariableNames;
        miss = setdiff(fn, vars);
        if ~isempty(miss)
            error('%s lacks feature columns: %s', srcFiles{k}, strjoin(miss, ', '));
        end
        h = height(T);

        Tf = standardizeMissing(T(:, fn), {'N/A','NA','na','','NONE'});
        Xk = zeros(h, numel(fn));
        for j = 1:numel(fn)
            col = Tf.(fn{j});
            if ~isnumeric(col); col = str2double(strtrim(string(col))); end
            Xk(:,j) = col;
        end

        if isnan(srcIter(1))
            if ~ismember('Iteration', vars)
                error('pooledFile needs an "Iteration" column.');
            end
            itk = double(T.Iteration);
        else
            itk = repmat(srcIter(k), h, 1);
        end

        Xc{end+1} = Xk; Ic{end+1} = itk; %#ok<AGROW>
        if ismember(cfg.useForMLCol, vars)
            u = T.(cfg.useForMLCol);
            if ~islogical(u) && ~isnumeric(u); u = strcmpi(strtrim(string(u)),'true') | strtrim(string(u)) == "1"; end
            Uc{k} = logical(u(:));
        else
            Uc{k} = true(h,1);
        end

        for i = 1:numel(labCols)
            lc = labCols{i};
            if ismember(lc, vars); s = strtrim(string(T.(lc)));
            else;                  s = repmat(string(missing), h, 1); end
            labC.(lc){end+1} = s(:);
        end
        for i = 1:numel(nameCols)
            nc = nameCols{i};
            if ismember(nc, vars); s = strtrim(string(T.(nc)));
            else;                  s = repmat(string(missing), h, 1); end
            nameC.(nc){end+1} = s(:);
        end
        for i = 1:numel(flagNames)
            fl = flagNames{i};
            if ismember(fl, vars)
                v = T.(fl);
                if ~islogical(v) && ~isnumeric(v); v = strcmpi(strtrim(string(v)),'true'); end
                v = double(v);
            else
                v = NaN(h,1);
            end
            flagC.(fl){end+1} = v(:);
        end
        fprintf('  loaded %-40s : %6d rows (iteration %s)\n', srcFiles{k}, h, ...
            strjoin(string(unique(itk)'), ','));
    end

    if isempty(Xc); error('No data files could be loaded.'); end

    D.X    = vertcat(Xc{:});
    D.Iter = vertcat(Ic{:});
    D.UseForML = vertcat(Uc{~cellfun(@isempty, Uc)});
    D.lab = struct(); D.name = struct(); D.flags = struct();
    for i = 1:numel(labCols);  D.lab.(labCols{i})   = vertcat(labC.(labCols{i}){:});  end
    for i = 1:numel(nameCols); D.name.(nameCols{i}) = vertcat(nameC.(nameCols{i}){:}); end
    for i = 1:numel(flagNames);D.flags.(flagNames{i})= vertcat(flagC.(flagNames{i}){:});end

    % missing-feature policy
    bad = any(isnan(D.X), 2);
    if any(bad)
        fprintf('  %d rows have missing features -> dropped.\n', sum(bad));
        D.X = D.X(~bad,:); D.Iter = D.Iter(~bad); D.UseForML = D.UseForML(~bad);
        fl = fieldnames(D.lab);   for i=1:numel(fl); D.lab.(fl{i})  = D.lab.(fl{i})(~bad);  end
        fl = fieldnames(D.name);  for i=1:numel(fl); D.name.(fl{i}) = D.name.(fl{i})(~bad); end
        fl = fieldnames(D.flags); for i=1:numel(fl); D.flags.(fl{i})= D.flags.(fl{i})(~bad);end
    end
end

%% ------------------------------------------------------------
function [R, S] = analyseLabelScheme(D, labelCol, cfg, mods, outRoot)
% Run every analysis for one label column.
    S = {};
    R = struct();
    outDir = fullfile(outRoot, ['labels_' labelCol]);
    figDir = fullfile(outDir, 'figures');
    if ~exist(figDir,'dir'); mkdir(figDir); end
    S{end+1} = sprintf('==================== LABEL SCHEME: %s ====================', labelCol);

    %% ---- labels --------------------------------------------------------
    labRaw = D.lab.(labelCol);
    keep   = ~ismissing(labRaw) & labRaw ~= "" & labRaw ~= "0";
    if ~any(keep)
        S{end+1} = 'No labels found for this scheme.'; return;
    end
    % (1) flagged LP = BP sets
    nFlag = sum(keep & ~D.UseForML);
    S{end+1} = sprintf('[DATA] Removed %d flagged LP = BP sets (Use_For_ML = 0): %s', nFlag, ...
        strjoin(arrayfun(@(a) sprintf('iteration %d: %d', a, sum(keep & ~D.UseForML & D.Iter == a)), ...
        unique(D.Iter)', 'UniformOutput', false), ', '));
    keep = keep & D.UseForML;
    % (2) classes too small for iteration-wise learning
    for a = unique(D.Iter)'
        la = labRaw(keep & D.Iter == a);
        [u, ~, g] = unique(la);
        cnt = accumarray(g, 1);
        for q = find(cnt < cfg.minClassPerIter)'
            S{end+1} = sprintf('[DATA] Iteration %d: class %s has only %d sets (< %d) -> removed from the ML.', ...
                a, u(q), cnt(q), cfg.minClassPerIter); %#ok<AGROW>
            keep = keep & ~(D.Iter == a & labRaw == u(q));
        end
    end
    X  = D.X(keep,:);  It = D.Iter(keep);  labS = labRaw(keep);

    nameS = [];
    if isKey(cfg.nameCols, labelCol) && isfield(D.name, cfg.nameCols(labelCol))
        nameAll = D.name.(cfg.nameCols(labelCol));
        nameS = nameAll(keep);
    end
    flags = struct();
    ff = fieldnames(D.flags);
    for i = 1:numel(ff); flags.(ff{i}) = D.flags.(ff{i})(keep); end

    lblNum = str2double(labS);
    if ~any(isnan(lblNum))
        catStr = string(sort(unique(lblNum)));
        labS   = string(lblNum);
    else
        catStr = sort(unique(labS));
    end
    catStr = catStr(:)';
    cats   = cellstr(catStr);
    K      = numel(cats);
    [~, yCode] = ismember(labS, catStr);

    className = cats;
    if ~isempty(nameS)
        for c = 1:K
            nm = nameS(find(yCode==c, 1));
            if ~ismissing(nm); className{c} = sprintf('%s: %s', cats{c}, char(nm)); end
        end
    end

    %% ---- features -------------------------------------------------------
    fn = cfg.featNames;  p = numel(fn);
    L  = log10(X);
    if any(~isfinite(L(:)))
        error('Non-positive or non-finite parameter values: cannot take log10.');
    end
    idxK = find(strcmp(fn,'KTOT'));  idxP = find(strcmp(fn,'PTOT'));
    if isempty(idxK) || isempty(idxP); error('KTOT / PTOT columns required.'); end

    pairs = nchoosek(1:p, 2);  nFP = size(pairs,1);
    Rlog  = L(:,pairs(:,1)) - L(:,pairs(:,2));
    ratioNames = cell(1,nFP);
    for fp = 1:nFP
        ratioNames{fp} = sprintf('%s/%s', fn{pairs(fp,1)}, fn{pairs(fp,2)});
    end

    iters = unique(It)';
    [~, itIdx] = ismember(It, iters);

    % design-constrained ratio = equal to 1 in >=99% of rows of some iteration
    fracOne = zeros(numel(iters), nFP);
    for a = 1:numel(iters)
        fracOne(a,:) = mean(abs(Rlog(itIdx==a,:)) < cfg.tolEq, 1);
    end
    isCons = any(fracOne >= 0.99, 1);

    % CLR representations (row-centred log-parameters)
    Lclr = L - mean(L, 2);
    rateIdx = setdiff(1:p, [idxK idxP]);
    Lr = L(:, rateIdx);  Lr = Lr - mean(Lr, 2);
    LclrSplit = [Lr, L(:,idxK) - L(:,idxP)];
    clrSplitNames = [strcat('CLR_', fn(rateIdx)), {'log(KTOT/PTOT)'}];

    FS = {};
    FS{end+1} = makeFS('A: 14 params + 91 log-ratios',            [L, Rlog],            [fn, ratioNames], false);
    FS{end+1} = makeFS('B: A minus design-constrained ratios',    [L, Rlog(:,~isCons)], [fn, ratioNames(~isCons)], false);
    FS{end+1} = makeFS('C: raw 14 parameters only',               L,                    fn, false);
    FS{end+1} = makeFS('D: JSD+broken-stick (original pipeline)', [L, Rlog],            [fn, ratioNames], true);
    FS{end+1} = makeFS('E: CLR of 14 log-params',                 Lclr,                 strcat('CLR_', fn), false);
    FS{end+1} = makeFS('F: CLR rates + log(KTOT/PTOT)',           LclrSplit,            clrSplitNames, false);

    key = It*1000 + yCode;
    splits = cell(1, cfg.nRepeats);
    for r = 1:cfg.nRepeats
        splits{r} = stratSplit(key, cfg.testFrac, cfg.seed + r, cfg.valFrac);
    end

    C = struct('L',L, 'Rlog',Rlog, 'pairs',pairs, 'ratioNames',{ratioNames}, 'names14',{fn}, ...
        'isCons',isCons, 'yCode',yCode, 'cats',{cats}, 'className',{className}, 'It',It, ...
        'iters',iters, 'itIdx',itIdx, 'K',K, 'FS',{FS}, 'key',key, 'splits',{splits}, ...
        'cfg',cfg, 'mods',mods, 'outDir',outDir, 'figDir',figDir, 'idxK',idxK, 'idxP',idxP, ...
        'fracOne',fracOne, 'flags',flags);

    % data shared with the workers once (feature sets + labels)
    Dshared = struct('FS', {FS}, 'yCode', yCode, 'itIdx', itIdx);
    if cfg.useParallel
        C.Dc = parallel.pool.Constant(Dshared);
    else
        C.Dc = Dshared;
    end
    C.useParallel = cfg.useParallel;

    S{end+1} = sprintf('Rows: %d | classes: %d (%s) | iterations: %s', numel(yCode), K, ...
        strjoin(cats, ','), strjoin(string(iters), ','));

    %% ---- sections --------------------------------------------------------
    s = secSplitFigure(C);                     S = [S, s];
    [s, R.confound] = secConfound(C);          S = [S, s];
    [s, R.ablation] = secAblation(C);          S = [S, s];
    [s, R.strata]   = secStrata(C);            S = [S, s];
    [s, R.decisive] = secDecisive(C);          S = [S, s];
    [s, R.confusion]= secConfusion(C);         S = [S, s];
    [s, R.models]   = secModelComparison(C);   S = [S, s];
    if cfg.runCurves
        [s, R.curves] = secCurves(C);           S = [S, s];
    end
    s = secSmallMultiples(C);                  S = [S, s];
    if cfg.runEmbeddingAudit
        [s, R.embedding] = secEmbedding(C);    S = [S, s];
    end
end

%% ------------------------------------------------------------
function [S, out] = secConfound(C)
% ML-b / EX: design confound audit and class x iteration structure.
    S = {};
    cfg = C.cfg; n = numel(C.yCode); iters = C.iters; nI = numel(iters); K = C.K;
    tol = cfg.tolEq;
    consIdx = find(C.isCons);
    consNames = C.ratioNames(consIdx);
    S{end+1} = sprintf('Design-constrained ratios (== 1 in >=99%% of some iteration): %d -> %s', ...
        numel(consIdx), strjoin(consNames, ', '));

    % fraction of rows equal to 1, per iteration
    vn = matlab.lang.makeValidName(strrep(consNames, '/', '_over_'));
    T1 = [table(iters(:), 'VariableNames', {'Iteration'}), ...
          array2table(C.fracOne(:,consIdx)*100, 'VariableNames', vn)];
    writetable(T1, fullfile(C.outDir,'confound_tables.xlsx'), 'Sheet', 'PctRowsEqualOne_ByIter');
    for a = 1:nI
        S{end+1} = sprintf('   Iteration %d: %s', iters(a), ...
            strjoin(arrayfun(@(j) sprintf('%s=1 in %.1f%%', consNames{j}, C.fracOne(a,consIdx(j))*100), ...
            1:numel(consIdx), 'UniformOutput', false), '; ')); %#ok<AGROW>
    end

    nEq = sum(abs(C.Rlog(:,consIdx)) < tol, 2);
    S{end+1} = sprintf('   Rows with >=4 constrained ratios exactly 1: %.1f%% | with all %d: %.1f%%', ...
        mean(nEq>=4)*100, numel(consIdx), mean(nEq==numel(consIdx))*100);

    % constraint signature -> sampling iteration
    sig = [abs(C.Rlog(:,consIdx)) < tol, sign(C.L(:,C.idxK) - C.L(:,C.idxP))];
    [~,~,g] = unique(sig, 'rows');
    M = accumarray([g, C.itIdx], 1);
    sigAcc = sum(max(M,[],2)) / n * 100;
    S{end+1} = sprintf('Constraint signature (%d distinct patterns) recovers the sampling iteration for %.2f%% of rows.', ...
        size(M,1), sigAcc);

    % class x iteration
    CT = accumarray([C.yCode, C.itIdx], 1, [K nI]);
    [pur, modal] = max(CT, [], 2);
    purity = pur ./ sum(CT,2) * 100;
    vnI = cellstr(matlab.lang.makeValidName(strcat('Iter', string(iters))));
    T2 = [table(C.className(:), 'VariableNames', {'Class'}), array2table(CT, 'VariableNames', vnI), ...
          table(sum(CT,2), purity, iters(modal)', 'VariableNames', {'Total','PctInModalIteration','ModalIteration'})];
    writetable(T2, fullfile(C.outDir,'confound_tables.xlsx'), 'Sheet', 'Class_by_Iteration');
    for c = 1:K
        S{end+1} = sprintf('   Class %-28s n=%6d | %.1f%% from a single iteration (iteration %d)', ...
            C.className{c}, sum(CT(c,:)), purity(c), iters(modal(c))); %#ok<AGROW>
    end
    S{end+1} = sprintf('Class-by-iteration table written. Majority-class baseline for the pooled data = %.2f%% (class %s).', ...
        max(sum(CT,2))/n*100, C.cats{find(sum(CT,2)==max(sum(CT,2)),1)});

    % Newton-verified descriptor flags (EX)
    fl = C.flags; names = fieldnames(fl);
    hasFlag = false;
    for i = 1:numel(names); if any(~isnan(fl.(names{i}))); hasFlag = true; end; end
    if hasFlag
        TF = table(iters(:), accumarray(C.itIdx,1), 'VariableNames', {'Iteration','N'});
        for i = 1:numel(names)
            v = fl.(names{i});
            pct = accumarray(C.itIdx, v, [], @(x) mean(x,'omitnan')) * 100;
            TF.(['Pct_' names{i}]) = pct;
        end
        writetable(TF, fullfile(C.outDir,'confound_tables.xlsx'), 'Sheet', 'Flags_by_Iteration');
        for a = 1:nI
            parts = arrayfun(@(i) sprintf('%s %.1f%%', names{i}, TF.(['Pct_' names{i}])(a)), ...
                1:numel(names), 'UniformOutput', false);
            S{end+1} = sprintf('Iteration %d (n=%d): %s', iters(a), TF.N(a), strjoin(parts, ' | ')); %#ok<AGROW>
        end
    end

    % figures
    fig = figure('Position',[50 50 1200 420],'Color','w');
    subplot(1,2,1);
    imagesc(C.fracOne(:,consIdx)*100); colorbar; caxis([0 100]);
    xticks(1:numel(consIdx)); xticklabels(consNames); xtickangle(45);
    yticks(1:nI); yticklabels(strcat('Iter ', string(iters)));
    title('% rows where ratio == 1 (design-constrained ratios)','FontWeight','bold');
    for a = 1:nI
        for j = 1:numel(consIdx)
            text(j, a, sprintf('%.0f', C.fracOne(a,consIdx(j))*100), 'HorizontalAlignment','center', 'Color','w','FontSize',8);
        end
    end
    subplot(1,2,2);
    imagesc(CT ./ sum(CT,2) * 100); colorbar; caxis([0 100]);
    xticks(1:nI); xticklabels(strcat('Iter ', string(iters)));
    yticks(1:K); yticklabels(C.cats);
    ylabel('Class'); title('% of each class by sampling iteration','FontWeight','bold');
    for c = 1:K
        for a = 1:nI
            text(a, c, sprintf('%.0f', CT(c,a)/sum(CT(c,:))*100), 'HorizontalAlignment','center','Color','w','FontSize',8);
        end
    end
    sgtitle('Design confound audit','FontWeight','bold');
    saveBothFormats(fig, C.figDir, 'confound_audit');

    out = struct('consNames',{consNames}, 'sigAcc',sigAcc, 'classByIter',CT, 'purity',purity);
end

%% ------------------------------------------------------------
function [S, T] = secAblation(C)
% ML-a / ML-c / CLR: pooled feature-set ablation with baselines.
    S = {};
    cfg = C.cfg; nFS = numel(C.FS); models = cfg.models; nM = numel(models); nR = cfg.nRepeats;
    results = cell(nFS, nM, nR);
    fprintf('\n[%s] Pooled ablation: %d feature sets x %d models x %d repeats\n', C.outDir, nFS, nM, nR);
    jobs = mkJob(); map = zeros(0,3);
    for r = 1:nR
        sp = C.splits{r};
        for f = 1:nFS
            for m = 1:nM
                jobs(end+1) = mkJob('eval', f, [], {}, false, [], 'yCode', sp.tr, sp.te, [], models{m}, cfg.seed + r); %#ok<AGROW>
                map(end+1,:) = [f m r]; %#ok<AGROW>
            end
        end
    end
    % leak check jobs (sampling iteration as the target) in the same batch
    for r = 1:nR
        sp = C.splits{r};
        jobs(end+1) = mkJob('eval', 1, [], {}, false, [], 'itIdx', sp.tr, sp.te, [], cfg.primaryModel, cfg.seed + r); %#ok<AGROW>
        map(end+1,:) = [0 0 r]; %#ok<AGROW>
    end
    out = runBatch(C, jobs, 'pooled ablation');
    itRes = cell(1,nR);
    for q = 1:numel(jobs)
        if map(q,1) > 0, results{map(q,1),map(q,2),map(q,3)} = out{q}; else, itRes{map(q,3)} = out{q}; end
    end

    agg = @(f,m,fld) arrayfun(@(r) results{f,m,r}.(fld), 1:nR);
    rows = {};
    accM = zeros(nFS,nM); accS = zeros(nFS,nM); baseM = zeros(nFS,nM);
    for f = 1:nFS
        for m = 1:nM
            a  = agg(f,m,'acc');   bT = agg(f,m,'baseTest'); bR = agg(f,m,'baseTrain');
            ba = agg(f,m,'balacc');mf = agg(f,m,'macroF1');  kp = agg(f,m,'kappa');
            nf = agg(f,m,'nFeat');
            accM(f,m) = mean(a); accS(f,m) = std(a); baseM(f,m) = mean(bT);
            rows(end+1,:) = {C.FS{f}.name, models{m}, mean(nf), mean(a), std(a), mean(bT), mean(bR), ...
                mean(a)-mean(bT), mean(ba), mean(mf), mean(kp)}; %#ok<AGROW>
            S{end+1} = sprintf('%-44s | %-12s | acc %.2f +/- %.2f | majority baseline %.2f | lift %+.2f | bal.acc %.2f | macroF1 %.2f | %g feats', ...
                C.FS{f}.name, models{m}, mean(a), std(a), mean(bT), mean(a)-mean(bT), mean(ba), mean(mf), mean(nf)); %#ok<AGROW>
        end
    end
    T = cell2table(rows, 'VariableNames', {'FeatureSet','Model','nFeatures','Acc_mean','Acc_sd', ...
        'Baseline_TestMajority','Baseline_TrainMajority','Lift_pts','BalancedAcc','MacroF1','Kappa'});
    writetable(T, fullfile(C.outDir,'ablation_pooled.xlsx'));

    % sampling-iteration recoverability from the parameters (leak check)
    itAcc = mean(cellfun(@(x) x.acc, itRes)); itBase = mean(cellfun(@(x) x.baseTest, itRes));
    S{end+1} = sprintf('%s predicting the SAMPLING ITERATION from feature set A: %.2f%% (majority baseline %.2f%%).', ...
        cfg.primaryModel, itAcc, itBase);

    fig = figure('Position',[50 50 1100 520],'Color','w');
    b = bar(accM); hold on;
    for m = 1:nM
        errorbar(b(m).XEndPoints, accM(:,m), accS(:,m), 'k.', 'LineWidth',1);
    end
    yline(mean(baseM(:)), 'r--', 'Majority baseline', 'LineWidth',1.5);
    xticks(1:nFS); xticklabels(cellfun(@(x) x.name, C.FS, 'UniformOutput', false)); xtickangle(25);
    ylabel('Accuracy (%)'); ylim([0 105]); legend(models, 'Location','southeast');
    title('Pooled accuracy by feature set (held-out, mean \pm sd)','FontWeight','bold');
    grid on; box on;
    saveBothFormats(fig, C.figDir, 'ablation_pooled_accuracy');
end

%% ------------------------------------------------------------
function [S, T] = secStrata(C)
% ML-d: classification within each sampling iteration.
    S = {};
    cfg = C.cfg; models = cfg.models; nM = numel(models); nR = cfg.nRepeats;
    nFS = numel(C.FS); iters = C.iters; nI = numel(iters);
    rows = {};
    impStore = cell(1,nI); confStore = cell(1,nI);

    % all fits of all iterations in one parallel batch
    jobs = mkJob(); map = zeros(0,4);
    for a = 1:nI
        idx = find(C.It == iters(a)); ys = C.yCode(idx);
        if numel(unique(ys)) < 2, continue; end
        for r = 1:nR
            sp = stratSplit(ys, cfg.testFrac, cfg.seed + 100*a + r, cfg.valFrac);
            for f = 1:nFS
                for m = 1:nM
                    jobs(end+1) = mkJob('eval', f, [], {}, false, [], 'yCode', idx(sp.tr), idx(sp.te), [], models{m}, cfg.seed + r); %#ok<AGROW>
                    map(end+1,:) = [a f m r]; %#ok<AGROW>
                end
            end
        end
    end
    outAll = runBatch(C, jobs, 'iteration-wise ablation');

    for a = 1:nI
        idx = find(C.It == iters(a));
        ys  = C.yCode(idx);
        nCls = numel(unique(ys));
        if nCls < 2
            S{end+1} = sprintf('Iteration %d: only one class present - baseline is 100%%, nothing to learn.', iters(a)); %#ok<AGROW>
            continue;
        end
        fprintf('  stratum iteration %d (%d rows, %d classes)\n', iters(a), numel(idx), nCls);
        res = cell(nFS, nM, nR);
        for q = find(map(:,1) == a)'
            res{map(q,2), map(q,3), map(q,4)} = outAll{q};
        end
        agg = @(f,m,fld) arrayfun(@(r) res{f,m,r}.(fld), 1:nR);
        for f = 1:nFS
            for m = 1:nM
                acc = agg(f,m,'acc'); bT = agg(f,m,'baseTest');
                rows(end+1,:) = {iters(a), numel(idx), nCls, C.FS{f}.name, models{m}, mean(agg(f,m,'nFeat')), ...
                    mean(acc), std(acc), mean(bT), mean(acc)-mean(bT), mean(agg(f,m,'balacc')), ...
                    mean(agg(f,m,'macroF1'))}; %#ok<AGROW>
                if strcmp(models{m}, cfg.primaryModel)
                    S{end+1} = sprintf('Iteration %d | %-44s | %s acc %.2f +/- %.2f vs baseline %.2f (lift %+.2f)', ...
                        iters(a), C.FS{f}.name, models{m}, mean(acc), std(acc), mean(bT), mean(acc)-mean(bT)); %#ok<AGROW>
                end
            end
        end
        mP = find(strcmp(models, cfg.primaryModel), 1);
        impStore{a}  = res{1,mP,1};
        confStore{a} = res{1,mP,1};
    end
    if isempty(rows); T = table(); return; end
    T = cell2table(rows, 'VariableNames', {'Iteration','N','nClasses','FeatureSet','Model','nFeatures', ...
        'Acc_mean','Acc_sd','Baseline_TestMajority','Lift_pts','BalancedAcc','MacroF1'});
    writetable(T, fullfile(C.outDir,'ablation_by_iteration.xlsx'));

    % importance within each stratum (feature set A: constants are dropped automatically)
    fig = figure('Position',[20 20 420*nI 520],'Color','w');
    tl = tiledlayout(1, nI, 'TileSpacing','compact','Padding','compact'); %#ok<NASGU>
    for a = 1:nI
        nexttile;
        rs = impStore{a};
        if isempty(rs); axis off; title(sprintf('Iteration %d (single class)', iters(a))); continue; end
        [v, o] = sort(rs.imp, 'descend'); nT = min(12, numel(v));
        hb = barh(v(1:nT), 'FaceColor', 'flat', 'EdgeColor', [0.15 0.15 0.15]); hb.CData = feature_colors(rs.impNames(o(1:nT)));   % one colour per feature
        set(gca,'YDir','reverse','FontSize',8);
        yticks(1:nT); yticklabels(rs.impNames(o(1:nT)));
        set(gca,'TickLabelInterpreter','none');
        xlabel('Importance (gain, normalised)');
        title(sprintf('Iteration %d', iters(a)), 'FontWeight','bold');
        S{end+1} = sprintf('Iteration %d top features (%s, set A): %s', iters(a), cfg.primaryModel, ...
            strjoin(rs.impNames(o(1:min(6,numel(o)))), ', ')); %#ok<AGROW>
    end
    sgtitle(sprintf('Within-iteration feature importance (%s)', cfg.primaryModel), 'FontWeight','bold');
    saveBothFormats(fig, C.figDir, 'importance_by_iteration');

    % confusion matrices per stratum (plain axes: exportable with exportgraphics / print)
    fig2 = figure('Position',[20 20 480*nI 520],'Color','w');
    tl2 = tiledlayout(fig2, 1, nI, 'TileSpacing','compact','Padding','compact');
    for a = 1:nI
        ax = nexttile(tl2);
        rs = confStore{a};
        if isempty(rs); axis(ax, 'off'); title(ax, sprintf('Iteration %d (single class)', iters(a))); continue; end
        present = unique([rs.yte; rs.pred]);
        cm = confusionmat(rs.yte, rs.pred, 'Order', present);
        drawCM(ax, cm, C.cats(present), sprintf('Iteration %d: acc %.1f%% (baseline %.1f%%)', iters(a), rs.acc, rs.baseTest));
    end
    title(tl2, 'Confusion matrices within each sampling iteration (row %, counts)', 'FontWeight','bold');
    saveBothFormats(fig2, C.figDir, 'confusion_by_iteration');
end

%% ------------------------------------------------------------
function [S, T] = secDecisive(C)
% ML-e: pairwise decisive test inside one iteration + linear boundary figure.
    S = {}; T = table();
    cfg = C.cfg; models = cfg.models; nM = numel(models); nR = cfg.nRepeats;
    tests = cfg.decisiveTests;
    rowsT = {};
    for t = 1:size(tests,1)
        itr = tests{t,1}; cA = string(tests{t,2}); cB = string(tests{t,3});
        if size(tests,2) >= 4 && ~isempty(tests{t,4}); ax = tests{t,4}; else; ax = cfg.smAxes; end
        [x1, xl1] = mechCoordSpec(C, ax{1,1}, ax{1,2});  [x2, xl2] = mechCoordSpec(C, ax{2,1}, ax{2,2});
        ca = find(strcmp(C.cats, cA)); cb = find(strcmp(C.cats, cB));
        if isempty(ca) || isempty(cb)
            S{end+1} = sprintf('Iteration %d, class %s vs %s: class not present in this label scheme - skipped.', itr, cA, cB); %#ok<AGROW>
            continue;
        end
        rows = find(C.It == itr & (C.yCode == ca | C.yCode == cb));
        yb = 1 + (C.yCode(rows) == cb);
        nA = sum(yb==1); nB = sum(yb==2);
        if min(nA,nB) < 10
            S{end+1} = sprintf('Iteration %d, class %s vs %s: too few samples (%d / %d) - skipped.', itr, cA, cB, nA, nB); %#ok<AGROW>
            continue;
        end
        tag = sprintf('Iter%d_Class%s_vs_%s', itr, cA, cB);
        fsAll  = struct('name','all non-constant features', 'X', C.FS{1}.X(rows,:), 'names', {C.FS{1}.names}, 'jsd', false);
        fsMech = struct('name','2 mechanistic log-ratios',  'X', [x1(rows), x2(rows)], 'names', {{xl1, xl2}}, 'jsd', false);

        accAll = zeros(nM,nR); accMech = zeros(nM,nR); accLin = zeros(1,nR); base = zeros(1,nR);
        res1 = []; svm1 = []; sp1 = [];
        sps = cell(1,nR); jobs = mkJob();
        for r = 1:nR
            sps{r} = stratSplit(yb, cfg.testFrac, cfg.seed + 200 + r, cfg.valFrac);
            for m = 1:nM
                jobs(end+1) = mkJob('eval', 0, fsAll.X,  fsAll.names,  false, yb, '', sps{r}.tr, sps{r}.te, [], models{m}, cfg.seed + r); %#ok<AGROW>
                jobs(end+1) = mkJob('eval', 0, fsMech.X, fsMech.names, false, yb, '', sps{r}.tr, sps{r}.te, [], models{m}, cfg.seed + r); %#ok<AGROW>
            end
        end
        outD = runBatch(C, jobs, ['decisive ' tag]);
        q = 0;
        for r = 1:nR
            sp = sps{r};
            for m = 1:nM
                q = q + 2; ra = outD{q-1}; rm = outD{q};
                accAll(m,r) = ra.acc; accMech(m,r) = rm.acc;
                if r == 1 && strcmp(models{m}, cfg.primaryModel); res1 = ra; end
                if m == 1; base(r) = ra.baseTest; end
            end
            svm = fitcsvm(fsMech.X(sp.tr,:), yb(sp.tr), 'KernelFunction','linear', ...
                'Standardize', false, 'BoxConstraint', 1);
            pl = predict(svm, fsMech.X(sp.te,:));
            accLin(r) = mean(pl == yb(sp.te)) * 100;
            if r == 1; svm1 = svm; sp1 = sp; end
        end

        for m = 1:nM
            rowsT(end+1,:) = {tag, models{m}, 'all non-constant features', mean(accAll(m,:)), std(accAll(m,:)), mean(base)}; %#ok<AGROW>
            rowsT(end+1,:) = {tag, models{m}, '2 mechanistic log-ratios',  mean(accMech(m,:)), std(accMech(m,:)), mean(base)}; %#ok<AGROW>
            S{end+1} = sprintf('%s | %-12s | all features %.2f +/- %.2f | 2 mech. log-ratios %.2f +/- %.2f | baseline %.2f', ...
                tag, models{m}, mean(accAll(m,:)), std(accAll(m,:)), mean(accMech(m,:)), std(accMech(m,:)), mean(base)); %#ok<AGROW>
        end
        rowsT(end+1,:) = {tag, 'Linear SVM', '2 mechanistic log-ratios', mean(accLin), std(accLin), mean(base)}; %#ok<AGROW>
        S{end+1} = sprintf('%s | LINEAR boundary on %s vs %s: %.2f +/- %.2f (baseline %.2f)', ...
            tag, xl1, xl2, mean(accLin), std(accLin), mean(base)); %#ok<AGROW>

        [v, o] = sort(res1.imp, 'descend');
        S{end+1} = sprintf('%s top features (%s): %s', tag, cfg.primaryModel, ...
            strjoin(cellfun(@(nm,iv) sprintf('%s (%.3f)', nm, iv), res1.impNames(o(1:min(6,numel(o)))), ...
            num2cell(v(1:min(6,numel(o)))), 'UniformOutput', false), ', ')); %#ok<AGROW>

        % ---- figure: mechanistic coordinates, colour-blind safe ----
        cOrange = classColor(ca, C, [0.90 0.60 0.00]); cBlue = classColor(cb, C, [0.00 0.45 0.70]);
        fig = figure('Position',[60 60 760 640],'Color','w');
        hold on;
        scatter(x1(rows(yb==1)), x2(rows(yb==1)), 16, cOrange, 'filled', 'MarkerFaceAlpha', 0.55);
        scatter(x1(rows(yb==2)), x2(rows(yb==2)), 16, cBlue,   'filled', 'MarkerFaceAlpha', 0.55);
        w = svm1.Beta; b0 = svm1.Bias;
        xr = xlim;
        if abs(w(2)) > 1e-12
            plot(xr, -(w(1)*xr + b0)/w(2), 'k-', 'LineWidth', 3);
        else
            xline(-b0/w(1), 'k-', 'LineWidth', 3);
        end
        xlim(xr);
        xlabel(xl1, 'Interpreter','none'); ylabel(xl2, 'Interpreter','none');
        legend({sprintf('Class %s (n=%d)', cA, nA), sprintf('Class %s (n=%d)', cB, nB), 'Linear boundary'}, ...
            'Location','best');
        title(sprintf('Iteration %d, class %s vs %s\nlinear boundary %.1f%% | boosted trees (2 feats) %.1f%% | baseline %.1f%%', ...
            itr, cA, cB, mean(accLin), mean(accMech(strcmp(models,cfg.primaryModel),:)), mean(base)), ...
            'FontWeight','bold');
        box on; grid on; set(gca, 'LineWidth', cfg.LW, 'FontSize', cfg.FS);
        saveBothFormats(fig, C.figDir, ['mechanistic_scatter_' tag]);
    end
    if ~isempty(rowsT)
        T = cell2table(rowsT, 'VariableNames', {'Test','Model','Features','Acc_mean','Acc_sd','Baseline'});
        writetable(T, fullfile(C.outDir,'decisive_tests.xlsx'));
    end
end

%% ------------------------------------------------------------
function [S, T] = secModelComparison(C)
% Original-manuscript model comparison: every model in cfg.compareModels, feature sets
% cfg.compareFeatureSets, held-out accuracy (split 1) + stratified k-fold CV on the
% training part (accuracy, macro-F1, weighted-F1; mean +/- sd), pooled and per iteration.
% Also one figure with the row-normalised confusion matrices of all models and one with
% the top-20 importances of all models (pooled, feature set A).
    S = {}; cfg = C.cfg; models = cfg.compareModels; nM = numel(models);
    scopes = {'pooled'};
    if contains(cfg.compareScopes, 'iteration')
        scopes = [scopes, arrayfun(@(a) sprintf('Iteration %d', a), C.iters, 'UniformOutput', false)];
    end
    rows = {}; keepPooled = cell(1, nM);
    % build every fit (held-out + k CV folds inside the training rows) as a parallel job
    jobs = mkJob(); map = zeros(0,4);              % scope, feature set, model, fold (0 = held-out)
    nIdx = zeros(1, numel(scopes)); k = cfg.cvFolds;
    for sc = 1:numel(scopes)
        if sc == 1
            idx = (1:numel(C.yCode))'; sp = C.splits{1}; tr = sp.tr; te = sp.te;
        else
            idx = find(C.It == C.iters(sc-1));
            if numel(unique(C.yCode(idx))) < 2; continue; end
            spl = stratSplit(C.yCode(idx), cfg.testFrac, cfg.seed + 500 + sc, cfg.valFrac);
            tr = idx(spl.tr); te = idx(spl.te);
        end
        nIdx(sc) = numel(idx);
        ytr = C.yCode(tr); fold = zeros(numel(tr), 1);
        rsF = RandStream('mt19937ar', 'Seed', cfg.seed + 7);
        for c = unique(ytr)'
            ii = find(ytr == c); ii = ii(randperm(rsF, numel(ii)));
            fold(ii) = mod(0:numel(ii)-1, k) + 1;
        end
        for f = cfg.compareFeatureSets
            for m = 1:nM
                jobs(end+1) = mkJob('eval', f, [], {}, false, [], 'yCode', tr, te, [], models{m}, cfg.seed); %#ok<AGROW>
                map(end+1,:) = [sc f m 0]; %#ok<AGROW>
                for q = 1:k
                    jobs(end+1) = mkJob('eval', f, [], {}, false, [], 'yCode', tr(fold ~= q), tr(fold == q), [], models{m}, cfg.seed + q); %#ok<AGROW>
                    map(end+1,:) = [sc f m q]; %#ok<AGROW>
                end
            end
        end
    end
    out = runBatch(C, jobs, 'model comparison');
    for sc = 1:numel(scopes)
        for f = cfg.compareFeatureSets
            fs = C.FS{f};
            for m = 1:nM
                h = find(map(:,1) == sc & map(:,2) == f & map(:,3) == m & map(:,4) == 0);
                if isempty(h), continue; end
                rh = out{h};
                cvq = find(map(:,1) == sc & map(:,2) == f & map(:,3) == m & map(:,4) > 0);
                cvA = cellfun(@(x) x.acc, out(cvq))'; cvM = cellfun(@(x) x.macroF1, out(cvq))';
                cvW = cellfun(@(x) x.weightedF1, out(cvq))';
                rows(end+1,:) = {scopes{sc}, fs.name, models{m}, nIdx(sc), rh.acc, rh.baseTest, ...
                    mean(cvA), std(cvA), mean(cvM), std(cvM), mean(cvW), std(cvW), rh.balacc, rh.kappa}; %#ok<AGROW>
                S{end+1} = sprintf('[MODELS] %-12s | %-40s | %-12s | held-out %.2f (baseline %.2f) | CV acc %.2f +/- %.2f | CV macro-F1 %.2f +/- %.2f | CV wt-F1 %.2f +/- %.2f', ...
                    scopes{sc}, fs.name, models{m}, rh.acc, rh.baseTest, mean(cvA), std(cvA), ...
                    mean(cvM), std(cvM), mean(cvW), std(cvW)); %#ok<AGROW>
                if sc == 1 && f == cfg.compareFeatureSets(1); keepPooled{m} = rh; end
            end
        end
    end
    T = cell2table(rows, 'VariableNames', {'Scope','FeatureSet','Model','N','HeldOut_Acc', ...
        'Baseline_TestMajority','CV_Acc_mean','CV_Acc_sd','CV_MacroF1_mean','CV_MacroF1_sd', ...
        'CV_WeightedF1_mean','CV_WeightedF1_sd','HeldOut_BalancedAcc','HeldOut_Kappa'});
    writetable(T, fullfile(C.outDir, 'model_comparison.xlsx'));

    % ---- figure 1: held-out accuracy per model and scope (feature set A) ----
    fA = C.FS{cfg.compareFeatureSets(1)}.name;
    uS = unique(T.Scope, 'stable');
    M = NaN(numel(uS), nM); B = NaN(numel(uS), 1);
    for i = 1:numel(uS)
        for m = 1:nM
            r = strcmp(T.Scope, uS{i}) & strcmp(T.Model, models{m}) & strcmp(T.FeatureSet, fA);
            if any(r); M(i,m) = T.HeldOut_Acc(find(r,1)); B(i) = T.Baseline_TestMajority(find(r,1)); end
        end
    end
    fig = figure('Position',[50 50 1250 520],'Color','w');
    bar(M); hold on;
    plot((1:numel(uS))', B, 'k_', 'MarkerSize', 40, 'LineWidth', 3);
    xticks(1:numel(uS)); xticklabels(uS); ylabel('Held-out accuracy (%)'); ylim([0 105]);
    legend([models, {'majority baseline'}], 'Location','southoutside','Orientation','horizontal','Interpreter','none');
    title(sprintf('Model comparison, feature set %s', fA), 'FontWeight','bold'); grid on; box on;
    saveBothFormats(fig, C.figDir, 'model_comparison_accuracy');

    % ---- figure 2: confusion matrices of all models (pooled, row-normalised) ----
    nc = min(4, nM); nr = ceil(nM/nc);
    fig2 = figure('Position',[20 20 420*nc 400*nr],'Color','w');
    tiledlayout(nr, nc, 'TileSpacing','compact','Padding','compact');
    for m = 1:nM
        rh = keepPooled{m}; if isempty(rh); continue; end
        cm = confusionmat(rh.yte, rh.pred, 'Order', 1:C.K);
        cmn = cm ./ max(sum(cm,2), 1) * 100;
        nexttile; imagesc(cmn, [0 100]); colormap(gca, cm_graded(256)); axis square;
        xticks(1:C.K); yticks(1:C.K); xticklabels(C.cats); yticklabels(C.cats);
        xlabel('Predicted'); ylabel('True');
        for i = 1:C.K, for j = 1:C.K
            if cm(i,j) > 0
                text(j, i, sprintf('%.0f', cmn(i,j)), 'HorizontalAlignment','center', 'FontSize', 8, ...
                    'Color', [1 1 1]*(cmn(i,j) > 50), 'FontWeight','bold');
            end
        end, end
        title(sprintf('%s: %.1f%%', models{m}, rh.acc), 'FontWeight','bold', 'Interpreter','none');
    end
    sgtitle('Pooled held-out confusion matrices (row %), feature set A', 'FontWeight','bold');
    saveBothFormats(fig2, C.figDir, 'confusion_all_models');

    % ---- figure 3: top-20 importances of all tree models (pooled) ----
    tm = find(~ismember(models, {'NaiveBayes','MostFrequent'}));
    nc = min(3, numel(tm)); nr = ceil(numel(tm)/nc);
    fig3 = figure('Position',[20 20 460*nc 560*nr],'Color','w');
    tiledlayout(nr, nc, 'TileSpacing','compact','Padding','compact');
    for q = tm
        rh = keepPooled{q}; if isempty(rh) || isempty(rh.imp) || all(rh.imp == 0); continue; end
        [v, o] = sort(rh.imp, 'descend'); nT = min(20, numel(v));
        nexttile; hb = barh(v(1:nT), 'FaceColor', 'flat', 'EdgeColor', [0.15 0.15 0.15]); hb.CData = feature_colors(rh.impNames(o(1:nT)));   % one colour per feature
        set(gca, 'YDir','reverse', 'TickLabelInterpreter','none');
        yticks(1:nT); yticklabels(rh.impNames(o(1:nT)));
        xlabel('Importance (normalised)'); title(models{q}, 'FontWeight','bold');
        S{end+1} = sprintf('[MODELS] %s pooled top features: %s', models{q}, ...
            strjoin(rh.impNames(o(1:min(6,nT))), ', ')); %#ok<AGROW>
    end
    sgtitle('Top-20 feature importances, pooled, feature set A', 'FontWeight','bold');
    saveBothFormats(fig3, C.figDir, 'importance_all_models');
end

%% ------------------------------------------------------------
function S = secSplitFigure(C)
% Train / test / validation split (pooled, repeat 1 = the split used for the confusion
% matrices, model comparison and curves): pie of the split, class distribution, and the
% number of sets of every class in each part (as in the original pipeline figure).
    S = {}; sp = C.splits{1}; K = C.K; n = numel(C.yCode);
    cTr = accumarray(C.yCode(sp.tr), 1, [K 1])';
    cTe = accumarray(C.yCode(sp.te), 1, [K 1])';
    cVa = accumarray(C.yCode(sp.va), 1, [K 1])';
    cAll = accumarray(C.yCode, 1, [K 1])';
    T = array2table([cTr; cTe; cVa; cAll], 'VariableNames', matlab.lang.makeValidName(strcat('Class_', C.cats)), ...
        'RowNames', {'Train','Test','Validation','Total'});
    writetable(T, fullfile(C.outDir, 'split_counts.xlsx'), 'WriteRowNames', true);
    S{end+1} = sprintf('[SPLIT] Train %d (%.1f%%) | Test %d (%.1f%%) | Validation %d (%.1f%%) | Total %d', ...
        numel(sp.tr), 100*numel(sp.tr)/n, numel(sp.te), 100*numel(sp.te)/n, numel(sp.va), 100*numel(sp.va)/n, n);
    for c = 1:K
        S{end+1} = sprintf('   class %-28s train %6d | test %5d | validation %5d', C.className{c}, cTr(c), cTe(c), cVa(c)); %#ok<AGROW>
    end

    colTr = [0.62 0.45 0.75]; colTe = [0.98 0.80 0.38]; colVa = [0.58 0.80 0.38];
    fig = figure('Position', [40 40 1500 980], 'Color', 'w');
    ax0 = axes('Parent', fig, 'Position', [0 0 1 1], 'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1]);
    hold(ax0, 'on');
    text(ax0, 0.5, 0.965, 'SPLITTING DATASET INTO TRAIN-TEST-VALIDATION', 'FontSize', 22, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    % pie
    axP = axes('Parent', fig, 'Position', [0.02 0.50 0.27 0.40]);
    hp = pie(axP, [numel(sp.tr) numel(sp.te) numel(sp.va)]);
    set(hp(1), 'FaceColor', colTr, 'LineWidth', 2); set(hp(3), 'FaceColor', colTe, 'LineWidth', 2);
    set(hp(5), 'FaceColor', colVa, 'LineWidth', 2);
    for j = [2 4 6], set(hp(j), 'FontSize', 15, 'FontWeight', 'bold'); end
    % legend block and total
    text(ax0, 0.30, 0.86, 'Per class split', 'FontSize', 20, 'FontWeight', 'bold');
    labs = {sprintf('Train (%.0f%%)', 100*numel(sp.tr)/n), sprintf('Test (%.0f%%)', 100*numel(sp.te)/n), ...
            sprintf('Validation (%.0f%%)', 100*numel(sp.va)/n)};
    cols = {colTr, colTe, colVa};
    for j = 1:3
        y = 0.80 - 0.045*(j-1);
        rectangle(ax0, 'Position', [0.30 y-0.013 0.022 0.026], 'FaceColor', cols{j}, 'EdgeColor', 'k', 'LineWidth', 1.5);
        text(ax0, 0.33, y, labs{j}, 'FontSize', 15, 'FontWeight', 'bold');
    end
    text(ax0, 0.30, 0.60, sprintf('Total: %d', n), 'FontSize', 20, 'FontWeight', 'bold');
    % class distribution
    axB = axes('Parent', fig, 'Position', [0.56 0.53 0.41 0.34]); hold(axB, 'on'); box(axB, 'on'); grid(axB, 'on');
    for c = 1:K
        bar(axB, c, cAll(c), 0.7, 'FaceColor', classColor(c, C, [0.1 0.4 0.75]), 'EdgeColor', 'k', 'LineWidth', 1.5);
        text(axB, c, cAll(c), {sprintf('%d', cAll(c)), sprintf('(%.1f%%)', 100*cAll(c)/n)}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9, 'FontWeight', 'bold');
    end
    set(axB, 'XTick', 1:K, 'XTickLabel', C.cats, 'XLim', [0.4 K+0.6], 'YLim', [0 1.25*max(cAll)], ...
        'LineWidth', 2, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel(axB, 'Class'); ylabel(axB, 'Number of Samples');
    title(axB, 'Original Class Distribution', 'FontSize', 14);
    % table: class circles + three rows
    x0 = 0.02; x1 = 0.98; wHead = 0.16; wc = (x1 - x0 - wHead)/K;
    yC = 0.43; rC = min(0.032, 0.42*wc);
    for c = 1:K
        xc = x0 + wHead + (c - 0.5)*wc;
        figPos = get(fig, 'Position'); ar = figPos(3)/figPos(4);
        rectangle(ax0, 'Position', [xc - rC, yC - rC*ar, 2*rC, 2*rC*ar], 'Curvature', [1 1], ...
            'FaceColor', classColor(c, C, [0.8 0.8 0.8]), 'EdgeColor', 'k', 'LineWidth', 2);
        text(ax0, xc, yC, C.cats{c}, 'FontSize', 20, 'HorizontalAlignment', 'center', 'FontWeight', 'normal');
    end
    yT = 0.36; hR = 0.105;
    rowsN = {sprintf('Train\n(%d)', numel(sp.tr)), sprintf('Test\n(%d)', numel(sp.te)), sprintf('Validation\n(%d)', numel(sp.va))};
    vals = {cTr, cTe, cVa}; rcol = {colTr, colTe, colVa};
    for r = 1:3
        yb = yT - r*hR;
        patch(ax0, [x0 x1 x1 x0], [yb yb yb+hR yb+hR], 1 - 0.55*(1 - rcol{r}), 'EdgeColor', 'none');
        patch(ax0, [x0 x0+wHead x0+wHead x0], [yb yb yb+hR yb+hR], 1 - 0.85*(1 - rcol{r}), 'EdgeColor', 'none');
        text(ax0, x0 + wHead/2, yb + hR/2, strsplit(rowsN{r}, newline), 'FontSize', 16, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center');
        for c = 1:K
            text(ax0, x0 + wHead + (c - 0.5)*wc, yb + hR/2, sprintf('%d', vals{r}(c)), 'FontSize', 16, ...
                'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end
    end
    yB = yT - 3*hR;
    for r = 0:3, plot(ax0, [x0 x1], [1 1]*(yT - r*hR), 'k-', 'LineWidth', 2); end
    plot(ax0, [x0 x0], [yB yT], 'k-', 'LineWidth', 2); plot(ax0, [x1 x1], [yB yT], 'k-', 'LineWidth', 2);
    for c = 0:K, plot(ax0, [1 1]*(x0 + wHead + c*wc), [yB yT], 'k-', 'LineWidth', 1.5); end
    try
        exportgraphics(fig, fullfile(C.figDir, 'train_test_validation_split.png'), 'Resolution', 300);
    catch
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, fullfile(C.figDir, 'train_test_validation_split'), '-dpng', '-r300');
    end
    savefig(fig, fullfile(C.figDir, 'train_test_validation_split.fig'));
end

%% ------------------------------------------------------------
function [S, T] = secCurves(C)
% ROC / AUC (one-vs-rest, test rows) and convergence curves (log-loss / error vs rounds on
% the validation rows and a training subsample, overall and per class) for every model,
% feature set A, pooled and (cfg.curveScopes) within every iteration. All fits run in one
% parallel batch; plotting is done afterwards on the client.
    S = {}; cfg = C.cfg; T = table();
    models = setdiff(cfg.compareModels, {'MostFrequent'}, 'stable');
    scopes = {'pooled'};
    if contains(cfg.curveScopes, 'iteration')
        scopes = [scopes, arrayfun(@(a) sprintf('Iteration %d', a), C.iters, 'UniformOutput', false)];
    end
    CI = []; try, CI = class_info_v3(); catch, end
    nM = numel(models);
    jobs = mkJob(); map = zeros(0,2);
    for sc = 1:numel(scopes)
        if sc == 1
            sp = C.splits{1}; tr = sp.tr; te = sp.te; va = sp.va;
        else
            idx = find(C.It == C.iters(sc-1));
            if numel(unique(C.yCode(idx))) < 2; continue; end
            spl = stratSplit(C.yCode(idx), cfg.testFrac, cfg.seed + 500 + sc, cfg.valFrac);
            tr = idx(spl.tr); te = idx(spl.te); va = idx(spl.va);
        end
        for m = 1:nM
            jobs(end+1) = mkJob('curve', 1, [], {}, false, [], 'yCode', tr, te, va, models{m}, cfg.seed); %#ok<AGROW>
            map(end+1,:) = [sc m]; %#ok<AGROW>
        end
    end
    out = runBatch(C, jobs, 'ROC + convergence curves');

    rowsA = {}; rowsC = {};
    colOf = @(c, cls, K) pickColor(c, cls, K, C, CI);
    for sc = 1:numel(scopes)
        qs = find(map(:,1) == sc)'; if isempty(qs), continue; end
        tag = matlab.lang.makeValidName(scopes{sc});
        nr = ceil(nM/4); nc = min(4, nM);
        figR = figure('Position',[20 20 430*nc 400*nr],'Color','w');
        tlR = tiledlayout(figR, nr, nc, 'TileSpacing','compact','Padding','compact');
        figL = figure('Position',[20 20 430*nc 380*nr],'Color','w');
        tlL = tiledlayout(figL, nr, nc, 'TileSpacing','compact','Padding','compact');
        figP = figure('Position',[20 20 430*nc 380*nr],'Color','w');
        tlP = tiledlayout(figP, nr, nc, 'TileSpacing','compact','Padding','compact');
        for q = qs
            o = out{q}; mdl = models{map(q,2)};
            if isempty(o.auc), continue; end
            cls = o.cls; K = numel(cls);
            % ROC
            ax = nexttile(tlR); hold(ax, 'on'); box(ax, 'on');
            for c = 1:K
                if isempty(o.fpr{c}), continue; end
                plot(ax, o.fpr{c}, o.tpr{c}, '-', 'Color', colOf(c, cls, K), 'LineWidth', 2, ...
                    'DisplayName', sprintf('%s (%.3f)', C.cats{cls(c)}, o.auc(c)));
                rowsA(end+1,:) = {scopes{sc}, mdl, C.className{cls(c)}, o.auc(c)}; %#ok<AGROW>
            end
            plot(ax, [0 1], [0 1], 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');
            axis(ax, [0 1 0 1]); axis(ax, 'square');
            xlabel(ax, 'False positive rate'); ylabel(ax, 'True positive rate');
            title(ax, sprintf('%s: macro AUC %.3f', mdl, mean(o.auc,'omitnan')), 'Interpreter','none');
            legend(ax, 'Location','southeast', 'FontSize', 7, 'Interpreter','none');
            S{end+1} = sprintf('[AUC] %-12s | %-12s | macro AUC %.4f | per class: %s', scopes{sc}, mdl, ...
                mean(o.auc,'omitnan'), strjoin(arrayfun(@(c) sprintf('%s %.3f', C.cats{cls(c)}, o.auc(c)), 1:K, ...
                'UniformOutput', false), ', ')); %#ok<AGROW>
            if isempty(o.steps), continue; end
            steps = o.steps;
            % overall convergence
            ax = nexttile(tlL); hold(ax, 'on'); box(ax, 'on');
            yyaxis(ax, 'left');
            plot(ax, steps, o.llTr, 'b-', 'LineWidth', 2, 'DisplayName', 'train log-loss');
            plot(ax, steps, o.llVa, 'r-', 'LineWidth', 2, 'DisplayName', 'validation log-loss');
            ylabel(ax, 'Multiclass log-loss');
            yyaxis(ax, 'right');
            plot(ax, steps, o.erTr, 'b--', 'LineWidth', 1.5, 'DisplayName', 'train error %');
            plot(ax, steps, o.erVa, 'r--', 'LineWidth', 1.5, 'DisplayName', 'validation error %');
            ylabel(ax, 'Error (%)'); xlabel(ax, roundLabel(mdl));
            last = steps >= steps(end) - 0.2*(steps(end) - steps(1)); i0 = find(last, 1);
            rel = (o.llVa(i0) - o.llVa(end)) / max(o.llVa(i0), eps) * 100;
            [~, ib] = min(o.llVa);
            title(ax, sprintf('%s: \\Delta last 20%% = %.2f%%, best at %d/%d', mdl, rel, steps(ib), steps(end)), ...
                'Interpreter','tex', 'FontSize', 9);
            legend(ax, 'Location','northeast', 'FontSize', 7);
            % per class
            ax = nexttile(tlP); hold(ax, 'on'); box(ax, 'on');
            for c = 1:K
                plot(ax, steps, o.errC(c,:), '-', 'Color', colOf(c, cls, K), 'LineWidth', 2, 'DisplayName', C.cats{cls(c)});
                rowsC(end+1,:) = {scopes{sc}, mdl, C.className{cls(c)}, o.errC(c,1), o.errC(c,end), ...
                    o.errC(c,i0) - o.errC(c,end), steps(end)}; %#ok<AGROW>
            end
            xlabel(ax, roundLabel(mdl)); ylabel(ax, 'Validation error of the class (%)');
            title(ax, mdl, 'Interpreter','none');
            legend(ax, 'Location','northeast', 'FontSize', 7, 'NumColumns', 2);
            [we, wc] = max(o.errC(:,end));
            S{end+1} = sprintf('[CONV] %-12s | %-12s | validation log-loss %.4f -> %.4f | change over last 20%% of rounds %.2f%% | min at %d of %d | worst class error at end: %s %.1f%%', ...
                scopes{sc}, mdl, o.llVa(1), o.llVa(end), rel, steps(ib), steps(end), C.cats{cls(wc)}, we); %#ok<AGROW>
        end
        title(tlR, sprintf('ROC curves (one-vs-rest), %s, feature set A', scopes{sc}), 'FontWeight','bold');
        title(tlL, sprintf('Convergence: log-loss and error vs rounds, %s', scopes{sc}), 'FontWeight','bold');
        title(tlP, sprintf('Convergence per class: validation error vs rounds, %s', scopes{sc}), 'FontWeight','bold');
        saveBothFormats(figR, C.figDir, ['ROC_' tag]);
        saveBothFormats(figL, C.figDir, ['convergence_' tag]);
        saveBothFormats(figP, C.figDir, ['convergence_per_class_' tag]);
        close([figR figL figP]);
    end
    if ~isempty(rowsA)
        T = cell2table(rowsA, 'VariableNames', {'Scope','Model','Class','AUC'});
        writetable(T, fullfile(C.outDir, 'roc_auc.xlsx'));
    end
    if ~isempty(rowsC)
        writetable(cell2table(rowsC, 'VariableNames', {'Scope','Model','Class','ValError_first', ...
            'ValError_last','ErrorDrop_last20pct','Rounds'}), fullfile(C.outDir, 'convergence_per_class.xlsx'));
    end
end

function col = pickColor(c, cls, K, C, CI)
    col = lines(K); col = col(c,:);
    if ~isempty(CI)
        k = str2double(C.cats{cls(c)});
        if ~isnan(k) && k >= 1 && k <= size(CI.colors,1), col = CI.colors(k,:); end
    end
end

function s = roundLabel(model)
    switch model
        case 'RandomForest', s = 'Number of trees';
        case 'NaiveBayes',   s = '-';
        otherwise,           s = 'Boosting rounds';
    end
end





%% ------------------------------------------------------------
function [S, out] = secConfusion(C)
% ML-f: pooled confusion matrix with baseline; audit of a chosen class pair.
    S = {}; cfg = C.cfg; K = C.K;
    sp = C.splits{1};
    res = evalFS(C.FS{1}, C.yCode, sp.tr, sp.te, cfg, C.mods, cfg.primaryModel, cfg.seed + 1);
    cm = confusionmat(res.yte, res.pred, 'Order', 1:K);
    rowSum = sum(cm, 2);
    S{end+1} = sprintf('Pooled %s (feature set A): held-out accuracy %.2f%% vs majority baseline %.2f%% (train-majority %.2f%%); balanced acc %.2f%%, macro-F1 %.2f%%.', ...
        cfg.primaryModel, res.acc, res.baseTest, res.baseTrain, res.balacc, res.macroF1);

    % per-class precision / recall / F1 (pooled, primary model)
    rowsP = {};
    for i = 1:K
        tp = cm(i,i); fp = sum(cm(:,i)) - tp; fn = rowSum(i) - tp;
        pr = tp/max(tp+fp,1)*100; re = tp/max(tp+fn,1)*100; f1 = 2*pr*re/max(pr+re,eps);
        rowsP(end+1,:) = {C.className{i}, rowSum(i), pr, re, f1}; %#ok<AGROW>
        S{end+1} = sprintf('   class %-28s n_test=%5d | precision %.1f | recall %.1f | F1 %.1f', C.className{i}, rowSum(i), pr, re, f1); %#ok<AGROW>
    end
    writetable(cell2table(rowsP, 'VariableNames', {'Class','N_test','Precision','Recall','F1'}), ...
        fullfile(C.outDir,'per_class_metrics_pooled.xlsx'));

    % top confusions
    rowsC = {};
    for i = 1:K
        for j = 1:K
            if i ~= j && cm(i,j) > 0
                rowsC(end+1,:) = {C.className{i}, C.className{j}, cm(i,j), rowSum(i), cm(i,j)/rowSum(i)*100}; %#ok<AGROW>
            end
        end
    end
    if ~isempty(rowsC)
        TC = cell2table(rowsC, 'VariableNames', {'TrueClass','PredictedAs','Count','TrueClassN','PctOfTrueClass'});
        TC = sortrows(TC, 'PctOfTrueClass', 'descend');
        writetable(TC, fullfile(C.outDir,'top_confusions.xlsx'));
        for i = 1:min(8, height(TC))
            S{end+1} = sprintf('   confusion: %s -> %s : %d of %d (%.1f%%)', TC.TrueClass{i}, TC.PredictedAs{i}, ...
                TC.Count(i), TC.TrueClassN(i), TC.PctOfTrueClass(i)); %#ok<AGROW>
        end
    end

    % explicit audit of the pairs named in cfg.confusionCheck
    for pc = 1:size(cfg.confusionCheck, 1)
    ia = find(strcmp(C.cats, cfg.confusionCheck{pc,1})); ib = find(strcmp(C.cats, cfg.confusionCheck{pc,2}));
    if ~isempty(ia) && ~isempty(ib)
        a2b = cm(ia,ib); b2a = cm(ib,ia);
        S{end+1} = sprintf('Class %s -> %s: %d of %d (%.2f%%) | Class %s -> %s: %d of %d (%.2f%%).', ...
            C.cats{ia}, C.cats{ib}, a2b, rowSum(ia), a2b/max(rowSum(ia),1)*100, ...
            C.cats{ib}, C.cats{ia}, b2a, rowSum(ib), b2a/max(rowSum(ib),1)*100);
        for q = [ia ib]
            err = cm(q,:); err(q) = 0;
            if sum(err) > 0
                [mx, dest] = max(err);
                S{end+1} = sprintf('   Errors of class %s go mostly to class %s (%d of %d errors, %.1f%%).', ...
                    C.cats{q}, C.cats{dest}, mx, sum(err), mx/sum(err)*100); %#ok<AGROW>
            else
                S{end+1} = sprintf('   Class %s has no misclassified test samples.', C.cats{q}); %#ok<AGROW>
            end
        end
    else
        S{end+1} = sprintf('Pair %s/%s not both present in this label scheme.', cfg.confusionCheck{pc,1}, cfg.confusionCheck{pc,2}); %#ok<AGROW>
    end
    end

    fig = figure('Position',[60 60 820 760],'Color','w');
    drawCM(axes(fig), cm, C.cats, sprintf('Pooled %s: acc %.1f%% | majority baseline %.1f%%', cfg.primaryModel, res.acc, res.baseTest));
    saveBothFormats(fig, C.figDir, 'confusion_pooled');

    out = struct('cm', cm, 'acc', res.acc, 'baseline', res.baseTest);
end

%% ------------------------------------------------------------
function S = secSmallMultiples(C)
% Figure suggestion: one panel per class (accent colour) against the grey ensemble.
    S = {};
    [x, xl] = mechCoordSpec(C, C.cfg.smAxes{1,1}, C.cfg.smAxes{1,2});
    [y, yl] = mechCoordSpec(C, C.cfg.smAxes{2,1}, C.cfg.smAxes{2,2});
    K = C.K; n = numel(x);
    nCols = min(4, K); nRows = ceil(K / nCols);
    accent = [0.84 0.37 0.00];                    % vermillion (Okabe-Ito)
    rng(C.cfg.seed);
    gi = randperm(n, min(n, 8000));
    fig = figure('Position',[40 40 nCols*330 nRows*300],'Color','w');
    tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');
    for c = 1:K
        nexttile; hold on;
        scatter(x(gi), y(gi), 4, [0.78 0.78 0.78], 'filled', 'MarkerFaceAlpha', 0.35);
        ic = C.yCode == c;
        scatter(x(ic), y(ic), 8, classColor(c, C, accent), 'filled', 'MarkerFaceAlpha', 0.75);
        title(sprintf('%s (n=%d)', C.className{c}, nnz(ic)), 'Interpreter','none', 'FontSize', 10);
        set(gca,'FontSize',10,'LineWidth',C.cfg.LW); box on;
        xline(0, 'k--', 'LineWidth', 1); yline(0, 'k--', 'LineWidth', 1);
        if mod(c-1, nCols) == 0; ylabel(yl, 'Interpreter','none'); end
        if c > K - nCols;       xlabel(xl, 'Interpreter','none'); end
    end
    sgtitle('Mechanistic coordinates: each class (colour) against the full ensemble (grey)', 'FontWeight','bold');
    saveBothFormats(fig, C.figDir, 'small_multiples_mechanistic');
    S{end+1} = sprintf('[FIG] Small multiples saved (axes: %s vs %s).', xl, yl);
end

%% ------------------------------------------------------------
function [S, T] = secEmbedding(C)
% Embedding audit: what does a 2-D embedding of the parameters encode?
    S = {}; cfg = C.cfg; n = numel(C.yCode);
    rng(cfg.seed);
    dKP = C.L(:,C.idxK) - C.L(:,C.idxP);
    regime = ones(n,1);  regime(abs(dKP) < cfg.tolEq) = 3;  regime(dKP < -cfg.tolEq) = 2;
    % regime 1: KTOT > PTOT | 2: KTOT < PTOT | 3: KTOT == PTOT
    regNames = {'KTOT>PTOT','KTOT<PTOT','KTOT=PTOT'};

    repNames = {'raw log-params', 'CLR (14)', 'CLR rates + log(KTOT/PTOT)'};
    repX     = {C.L, C.FS{5}.X, C.FS{6}.X};

    % stratified-by-regime subsample for the global embedding
    sub = [];
    for g = 1:3
        ig = find(regime == g);
        if isempty(ig); continue; end
        take = min(numel(ig), ceil(cfg.embedN * numel(ig) / n));
        sub = [sub; ig(randperm(numel(ig), take))]; %#ok<AGROW>
    end

    rows = {};
    % ---- global embedding: iteration vs class ----
    for r = 1:numel(repNames)
        fprintf('  embedding (%s) of %d points...\n', repNames{r}, numel(sub));
        Y = embed2d(repX{r}(sub,:), cfg.embedMethod, cfg.seed);
        accIt  = knnCV(Y, C.It(sub),    cfg.knnK, cfg.seed);
        accCl  = knnCV(Y, C.yCode(sub), cfg.knnK, cfg.seed);
        bIt = majorityShare(C.It(sub)); bCl = majorityShare(C.yCode(sub));
        S{end+1} = sprintf('[EMB] Global embedding of %s: sampling iteration recoverable at %.1f%% (baseline %.1f%%); class at %.1f%% (baseline %.1f%%).', ...
            repNames{r}, accIt, bIt, accCl, bCl); %#ok<AGROW>
        rows(end+1,:) = {'global', repNames{r}, numel(sub), numel(unique(C.yCode(sub))), bCl, accCl, NaN, NaN, accIt, bIt}; %#ok<AGROW>
        if r == 1
            fig = figure('Position',[60 60 1250 520],'Color','w');
            subplot(1,2,1);
            scatter(Y(:,1), Y(:,2), 6, C.itIdx(sub), 'filled'); colormap(gca, lines(numel(C.iters))); colorbar;
            title(sprintf('Embedding of raw log-params coloured by SAMPLING ITERATION\n(kNN recovers it at %.0f%%, baseline %.0f%%)', accIt, bIt));
            subplot(1,2,2);
            scatter(Y(:,1), Y(:,2), 6, C.yCode(sub), 'filled'); colormap(gca, lines(C.K)); colorbar;
            title(sprintf('Same embedding coloured by BIFURCATION CLASS\n(kNN recovers it at %.0f%%, baseline %.0f%%)', accCl, bCl));
            saveBothFormats(fig, C.figDir, 'embedding_iteration_vs_class');
        end
    end

    % ---- per-regime re-embedding vs supervised ----
    for g = 1:3
        ig = find(regime == g);
        if numel(ig) < 50 || numel(unique(C.yCode(ig))) < 2; continue; end
        if numel(ig) > cfg.embedPerRegime; ig = ig(randperm(numel(ig), cfg.embedPerRegime)); end
        yg = C.yCode(ig);
        bCl = majorityShare(yg);
        sp = stratSplit(yg, cfg.testFrac, cfg.seed + 300 + g, cfg.valFrac);
        for r = 1:numel(repNames)
            Y = embed2d(repX{r}(ig,:), cfg.embedMethod, cfg.seed);
            accKnn = knnCV(Y, yg, cfg.knnK, cfg.seed);
            fsT = struct('name', repNames{r}, 'X', repX{r}(ig,:), 'names', {C.FS{1}.names(1:size(repX{r},2))}, 'jsd', false);
            rs = evalFS(fsT, yg, sp.tr, sp.te, cfg, C.mods, cfg.primaryModel, cfg.seed);
            lift = NaN;
            if rs.acc > bCl; lift = (accKnn - bCl) / (rs.acc - bCl) * 100; end
            rows(end+1,:) = {regNames{g}, repNames{r}, numel(ig), numel(unique(yg)), bCl, accKnn, rs.acc, lift, NaN, NaN}; %#ok<AGROW>
            S{end+1} = sprintf('[EMB] Regime %-10s | %-28s | embedding kNN %.1f%% | supervised %.1f%% | baseline %.1f%% | lift recovered %.0f%%', ...
                regNames{g}, repNames{r}, accKnn, rs.acc, bCl, lift); %#ok<AGROW>
        end
    end

    T = cell2table(rows, 'VariableNames', {'Scope','Representation','N','nClasses','Baseline_Class', ...
        'Embedding_kNN_Class','Supervised_Class','PctLiftRecovered','Embedding_kNN_Iteration','Baseline_Iteration'});
    writetable(T, fullfile(C.outDir,'embedding_audit.xlsx'));
    S{end+1} = '[EMB] Reading: an embedding that encodes iteration but not class cannot be used as evidence about class boundaries (revise the t-SNE passage accordingly).';
end

%% ------------------------------------------------------------
function [x, lab] = mechCoord(C, r)
% log10 ratio for the r-th mechanistic axis in cfg.smAxes.
    a = C.cfg.smAxes{r,1}; b = C.cfg.smAxes{r,2};
    ia = find(strcmp(C.names14, a)); ib = find(strcmp(C.names14, b));
    if isempty(ia) || isempty(ib)
        error('cfg.smAxes entry %s/%s not found among the feature names.', a, b);
    end
    x = C.L(:,ia) - C.L(:,ib);
    lab = sprintf('log10(%s/%s)', a, b);
end

function [x, lab] = mechCoordSpec(C, numSpec, denSpec)
% log10( prod(numerator parameters) / prod(denominator parameters) ); specs like 'k2*KTOT'
    x = termLog(C, numSpec) - termLog(C, denSpec);
    lab = sprintf('log10(%s / %s)', strrep(numSpec,'*',' '), strrep(denSpec,'*',' '));
end

function v = termLog(C, spec)
    parts = strtrim(strsplit(spec, '*'));
    v = zeros(size(C.L,1), 1);
    for i = 1:numel(parts)
        j = find(strcmp(C.names14, parts{i}));
        if isempty(j); error('Parameter %s (in %s) not among the features.', parts{i}, spec); end
        v = v + C.L(:,j);
    end
end

function col = classColor(c, C, fallback)
% colour of class code c from class_info_v3 (support folder), if available
    col = fallback;
    try
        CI = class_info_v3();
        k = str2double(C.cats{c});
        if ~isnan(k) && k >= 1 && k <= size(CI.colors,1); col = CI.colors(k,:); end
    catch
    end
end

function res = evalFS(fs, y, tr, te, cfg, ~, modelName, seed)
% single fit on the client (all threads); same core as the workers
    res = ml_v3_core('evalFS', fs, y, tr, te, cfg, modelName, seed);
end

function b = majorityShare(y)
    b = ml_v3_core('majorityShare', y);
end

function j = mkJob(kind, fsId, X, names, jsd, y, yName, tr, te, va, model, seed)
% one model fit for the parallel batch (mkJob() = empty job array)
    j = struct('kind', {}, 'fsId', {}, 'X', {}, 'names', {}, 'jsd', {}, 'y', {}, 'yName', {}, ...
               'tr', {}, 'te', {}, 'va', {}, 'model', {}, 'seed', {});
    if nargin == 0, return; end
    j(1).kind = kind; j.fsId = fsId; j.X = X; j.names = names; j.jsd = jsd; j.y = y; j.yName = yName;
    j.tr = uint32(tr(:)); j.te = uint32(te(:)); j.va = uint32(va(:)); j.model = model; j.seed = seed;
end

function c = jobCost(j)
% rough relative run time, used to send the longest jobs first
    w = struct('GradBoost', 40, 'CatBoost', 8, 'RandomForest', 5, 'XGBoost', 4, 'LightGBM', 2, ...
               'NaiveBayes', 0.05, 'MostFrequent', 0.01);
    if isfield(w, j.model), c = w.(j.model); else, c = 1; end
    c = c * numel(j.tr);
    if strcmp(j.kind, 'curve'), c = c * 3; end
end

function out = runBatch(C, jobs, label)
% run all jobs: in parallel on the pool (one thread per model) or serially on the client
    n = numel(jobs); out = cell(n, 1);
    if n == 0, return; end
    t0 = tic;
    fprintf('  [%s] %d model fits ', label, n);
    if C.useParallel
        cost = arrayfun(@jobCost, jobs); [~, ord] = sort(cost, 'descend');
        jobsO = jobs(ord); cfgW = C.cfg; cfgW.nThreads = 1; Dc = C.Dc; res = cell(n, 1);
        pool = gcp(); fprintf('on %d workers ... ', pool.NumWorkers);
        opts = parforOptions(pool, 'RangePartitionMethod', 'fixed', 'SubrangeSize', 1);
        parfor (i = 1:n, opts)
            res{i} = ml_v3_core('job', jobsO(i), Dc.Value, cfgW);
        end
        out(ord) = res;
    else
        fprintf('serially ... ');
        for i = 1:n, out{i} = ml_v3_core('job', jobs(i), C.Dc, C.cfg); end
    end
    fprintf('done in %.1f min\n', toc(t0)/60);
end

function fs = makeFS(name, X, names, useJSD)
    fs = struct('name', name, 'X', X, 'names', {names}, 'jsd', useJSD);
end

function sp = stratSplit(key, testFrac, seed, valFrac)
% Stratified TRAIN / VALIDATION / TEST split on an integer key (class, or iteration x class).
% sp.tr = train (models are fitted ONLY on these rows), sp.va = validation (convergence
% curves), sp.te = held-out test (all reported accuracies). Every stratum keeps >= 1 train row.
    if nargin < 4; valFrac = 0; end
    rng(seed);
    n = numel(key); part = zeros(n,1);                % 0 train, 1 validation, 2 test
    u = unique(key);
    for i = 1:numel(u)
        idx = find(key == u(i)); m = numel(idx);
        if m < 2; continue; end                       % singletons stay in train
        perm = idx(randperm(m));
        nTe = min(max(1, round(testFrac*m)), m-1);
        nVa = min(round(valFrac*m), m - 1 - nTe);
        part(perm(1:nTe)) = 2;
        part(perm(nTe+1:nTe+nVa)) = 1;
    end
    sp.tr = find(part == 0); sp.va = find(part == 1); sp.te = find(part == 2);
end







function acc = knnCV(Y, lab, k, seed)
% 5-fold cross-validated kNN accuracy of `lab` inside the embedding Y.
    rng(seed);
    lab = lab(:);
    try
        cvp = cvpartition(lab, 'KFold', 5);
    catch
        cvp = cvpartition(numel(lab), 'KFold', 5);
    end
    mdl = fitcknn(Y, lab, 'NumNeighbors', k);
    cvm = crossval(mdl, 'CVPartition', cvp);
    acc = (1 - kfoldLoss(cvm)) * 100;
end

function Y = embed2d(X, method, seed)
% 2-D embedding: UMAP if requested and installed, otherwise t-SNE.
    rng(seed);
    if strcmpi(method, 'umap')
        try
            um  = py.importlib.import_module('umap');
            red = um.UMAP(pyargs('n_neighbors', int32(15), 'min_dist', 0.1, 'random_state', int32(seed)));
            Y = double(red.fit_transform(py.numpy.array(double(X))));
            return;
        catch ME
            warning('UMAP unavailable (%s); using t-SNE.', ME.message);
        end
    end
    perp = max(5, min(30, floor((size(X,1) - 1) / 3)));
    Y = tsne(X, 'Perplexity', perp, 'Standardize', false);
end




function drawCM(ax, cm, cats, ttl)
% confusion matrix as a heatmap: colour = row-normalised %, text = count (row %)
    K = size(cm,1); cmn = cm ./ max(sum(cm,2), 1) * 100;
    imagesc(ax, cmn, [0 100]); colormap(ax, cm_graded(256)); axis(ax, 'square'); hold(ax, 'on');
    for i = 1:K
        for j = 1:K
            if cm(i,j) > 0
                text(ax, j, i, {sprintf('%d', cm(i,j)), sprintf('%.0f%%', cmn(i,j))}, 'HorizontalAlignment','center', ...
                    'FontSize', 7, 'FontWeight','bold', 'Color', [1 1 1]*(cmn(i,j) > 50));
            end
        end
    end
    set(ax, 'XTick', 1:K, 'YTick', 1:K, 'XTickLabel', cats, 'YTickLabel', cats, 'TickLabelInterpreter','none', ...
        'LineWidth', 1.8, 'FontWeight','bold');
    xlabel(ax, 'Predicted class'); ylabel(ax, 'True class');
    title(ax, ttl, 'Interpreter','none');
    cb = colorbar(ax); ylabel(cb, 'Row %');
end

function saveBothFormats(fig, figDir, name)
    % thicker axes and lines, larger fonts in every figure
    ax = findall(fig, 'Type', 'axes');
    for i = 1:numel(ax)
        set(ax(i), 'LineWidth', 1.8, 'FontWeight', 'bold');
        if get(ax(i), 'FontSize') < 11; set(ax(i), 'FontSize', 11); end
        box(ax(i), 'on');
    end
    ln = findall(fig, 'Type', 'line');
    for i = 1:numel(ln)
        if get(ln(i), 'LineWidth') < 2; set(ln(i), 'LineWidth', 2); end
    end
    % exportgraphics cannot export figures with several containers (e.g. confusion charts in
    % uipanels): fall back to print, which saves the whole figure
    try
        exportgraphics(fig, fullfile(figDir,[name,'.png']), 'Resolution', 300);
    catch
        try
            set(fig, 'PaperPositionMode', 'auto');
            print(fig, fullfile(figDir, name), '-dpng', '-r300');
        catch ME
            warning('PNG of %s not written (%s); the .fig is saved.', name, ME.message);
        end
    end
    savefig(fig, fullfile(figDir,[name,'.fig']));
end