%% ============================================================
%  ML ANALYSIS OF THE SINGLE-START AND EXACT LABELS
%
%  Reads the per-iteration datasets (Newton-verified labels), pools
%  them, and produces the numbers and figures of the analyses below:
%
%   ML-a  Majority-class baseline reported with EVERY accuracy
%         (test-majority share and train-majority accuracy, plus
%          balanced accuracy, macro-F1, Cohen's kappa, lift).
%   ML-b  Confound audit: which ratios are fixed by the sampling design,
%         how often they equal exactly 1, whether the constraint
%         signature recovers the sampling iteration, and how much of each
%         class is determined by iteration alone.
%   ML-c  Feature-set ablation on the pooled data:
%           A  14 params + all 91 pairwise ratios
%           B  A with the design-constrained ratios deleted
%           C  raw 14 parameters only
%           D  JSD + broken-stick selected features (original pipeline)
%           E  CLR (centred log-ratio) of the 14 log-parameters      [CLR]
%           F  CLR of the 12 rates + log(KTOT/PTOT)                  [CLR]
%   ML-d  Stratified-by-iteration classification (constrained ratios have
%         zero variance inside a stratum) + per-stratum importances and
%         confusion matrices.
%   ML-e  "Decisive test": Class 4 vs Class 5 inside Iteration 4
%         (configurable), with baseline, ranked features, and a LINEAR
%         boundary on log10(k1/k5) vs log10(k11/k55).
%   ML-f  Pooled confusion matrix re-examined (row-normalised, top
%         confusions, explicit check of the Class 1 <-> Class 2 claim).
%   FIG   Mechanistic-coordinate scatter (Okabe-Ito colour-blind safe) and
%         small multiples: one panel per class against the grey ensemble.
%   EMB   Embedding audit (t-SNE, or UMAP if installed): what does the
%         embedding encode — iteration or class?  Per-regime re-embedding,
%         raw log-params vs CLR, and % of the supervised lift recovered.
%   EX    Class x iteration table + bistable / multistationary fractions
%         from the Newton-verified flag columns.
%
%  Every analysis is run once per label scheme in cfg.labelCols:
%     'Exact_Class_v2'  — new (Newton-verified) 9-class labels
%     'Class_Number_7'  — single-start 7-class labels
%
%  Notes on choices
%   * Features are log10(parameters) and log10 pairwise ratios
%     (= differences of log-params). Tree models are invariant to this
%     monotone transform; it makes JSD binning and CLR well behaved and
%     lets "ratio == 1" be tested exactly (log-ratio == 0).
%   * Only the 14 kinetic parameters are ever used as inputs. The
%     Bistable / Multistationary / Has_Detached_Branch columns are label
%     descriptors and are used ONLY in the descriptive tables.
%   * Inside a stratum, exactly duplicated feature columns are collapsed (keep first),
%     so importances are not diluted across identical ratios (Iteration 3: 98 non-constant
%     columns but only 43 distinct; Iteration 4: 100 -> 59). A listed feature stands for
%     all of its aliases (e.g. k1/k2 also means k5/k6 there).
%   * Every split is stratified (pooled: by iteration x class) and
%     repeated cfg.nRepeats times; results are mean +/- sd.
% =============================================================

clear; clc; close all;

%% ============================================================
%% 0. CONFIG
%% ============================================================

cfg = struct();
cfg.dataDir    = '.';                       % folder holding the iteration files
cfg.iterIds    = [1 2 3 4];
cfg.iterFiles  = arrayfun(@(k) sprintf('ML_Ready_Dataset_v2_iter%d.xlsx',k), ...
                          cfg.iterIds, 'UniformOutput', false);
cfg.pooledFile = '';                        % OR one pooled file with an 'Iteration' column

cfg.labelCols  = {'Exact_Class_v2','Class_Number_7'};
cfg.nameCols   = containers.Map({'Exact_Class_v2','Class_Number_7'}, ...
                                {'Exact_Label_v2','Multiclass_Label_7'});
cfg.flagCols   = {'Bistable','Multistationary','Has_Detached_Branch'};

cfg.featNames  = {'KTOT','PTOT','k1','k11','k2','k3','k33','k4','k5',...
                  'k55','k6','k7','k77','k8'};

cfg.models       = {'LightGBM','RandomForest','NaiveBayes'};
cfg.primaryModel = 'LightGBM';              % used for importances / confusion / figures
cfg.useClassWeights = true;
cfg.lgbRounds  = 300;
cfg.rfTrees    = 300;

cfg.nRepeats   = 3;                         % random stratified splits per experiment
cfg.testFrac   = 0.30;
cfg.seed       = 42;
cfg.nBins      = 30;                        % JSD histogram bins
cfg.tolEq      = 1e-9;                      % |log10 ratio| < tol  <=>  ratio == 1
cfg.missingPolicy = 'drop';                 % rows with a missing feature are dropped

% Decisive tests: {iteration, classA, classB}  (one row per test)
cfg.decisiveTests = {4, '4', '5'};

% Mechanistic coordinates for the figures: {numerator, denominator}
cfg.smAxes = {'k1','k5'; 'k11','k55'};

% Pair of classes whose mutual confusion is audited (ML-f)
cfg.confusionCheck = {'1','2'};

% Embedding audit
cfg.runEmbeddingAudit = true;
cfg.embedMethod   = 'tsne';                 % 'tsne' or 'umap' (needs umap-learn)
cfg.embedN        = 6000;                   % global embedding size
cfg.embedPerRegime= 3000;                   % per-regime embedding size
cfg.knnK          = 15;

cfg.quick = false;                          % true = 1 repeat, fewer trees, small embeddings
if cfg.quick
    cfg.nRepeats = 1; cfg.lgbRounds = 100; cfg.rfTrees = 100;
    cfg.embedN = 2500; cfg.embedPerRegime = 1500;
end

outRoot = 'ML_label_comparison';
if ~exist(outRoot,'dir'); mkdir(outRoot); end

%% ============================================================
%% 1. PYTHON MODULES
%% ============================================================

fprintf('\n--- Installing / verifying Python packages ---\n');
try
    pyExec = char(pyenv().Executable);
    pkgs = 'lightgbm scikit-learn';
    if strcmpi(cfg.embedMethod,'umap'); pkgs = [pkgs ' umap-learn']; end
    system(sprintf('"%s" -m pip install %s --quiet', pyExec, pkgs));
catch
    warning('Could not install Python packages automatically.');
end

mods = struct();
mods.np     = py.importlib.import_module('numpy');
mods.lgb    = py.importlib.import_module('lightgbm');
mods.sk_ens = py.importlib.import_module('sklearn.ensemble');
mods.sk_nb  = py.importlib.import_module('sklearn.naive_bayes');

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

summaryFile = fullfile(outRoot, 'ML_label_comparison_summary.txt');
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

    Xc = {}; Ic = {};
    labC = struct(); nameC = struct(); flagC = struct();
    for i = 1:numel(labCols);  labC.(labCols{i})   = {}; end
    for i = 1:numel(nameCols); nameC.(nameCols{i}) = {}; end
    for i = 1:numel(flagNames);flagC.(flagNames{i})= {}; end

    for k = 1:numel(srcFiles)
        f = fullfile(cfg.dataDir, srcFiles{k});
        if ~isfile(f)
            warning('File not found, skipped: %s', f);
            continue;
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
    D.lab = struct(); D.name = struct(); D.flags = struct();
    for i = 1:numel(labCols);  D.lab.(labCols{i})   = vertcat(labC.(labCols{i}){:});  end
    for i = 1:numel(nameCols); D.name.(nameCols{i}) = vertcat(nameC.(nameCols{i}){:}); end
    for i = 1:numel(flagNames);D.flags.(flagNames{i})= vertcat(flagC.(flagNames{i}){:});end

    % missing-feature policy
    bad = any(isnan(D.X), 2);
    if any(bad)
        fprintf('  %d rows have missing features -> dropped.\n', sum(bad));
        D.X = D.X(~bad,:); D.Iter = D.Iter(~bad);
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
    keep   = ~ismissing(labRaw) & labRaw ~= "";
    if ~any(keep)
        S{end+1} = 'No labels found for this scheme.'; return;
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
        splits{r} = stratSplit(key, cfg.testFrac, cfg.seed + r);
    end

    C = struct('L',L, 'Rlog',Rlog, 'pairs',pairs, 'ratioNames',{ratioNames}, 'names14',{fn}, ...
        'isCons',isCons, 'yCode',yCode, 'cats',{cats}, 'className',{className}, 'It',It, ...
        'iters',iters, 'itIdx',itIdx, 'K',K, 'FS',{FS}, 'key',key, 'splits',{splits}, ...
        'cfg',cfg, 'mods',mods, 'outDir',outDir, 'figDir',figDir, 'idxK',idxK, 'idxP',idxP, ...
        'fracOne',fracOne, 'flags',flags);

    S{end+1} = sprintf('Rows: %d | classes: %d (%s) | iterations: %s', numel(yCode), K, ...
        strjoin(cats, ','), strjoin(string(iters), ','));

    %% ---- sections --------------------------------------------------------
    [s, R.confound] = secConfound(C);          S = [S, s];
    [s, R.ablation] = secAblation(C);          S = [S, s];
    [s, R.strata]   = secStrata(C);            S = [S, s];
    [s, R.decisive] = secDecisive(C);          S = [S, s];
    [s, R.confusion]= secConfusion(C);         S = [S, s];
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
    for r = 1:nR
        sp = C.splits{r};
        for f = 1:nFS
            for m = 1:nM
                results{f,m,r} = evalFS(C.FS{f}, C.yCode, sp.tr, sp.te, cfg, C.mods, models{m}, cfg.seed + r);
            end
        end
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
    itRes = cell(1,nR);
    for r = 1:nR
        sp = C.splits{r};
        itRes{r} = evalFS(C.FS{1}, C.itIdx, sp.tr, sp.te, cfg, C.mods, cfg.primaryModel, cfg.seed + r);
    end
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
        for r = 1:nR
            sp = stratSplit(ys, cfg.testFrac, cfg.seed + 100*a + r);
            trG = idx(sp.tr); teG = idx(sp.te);
            for f = 1:nFS
                for m = 1:nM
                    res{f,m,r} = evalFS(C.FS{f}, C.yCode, trG, teG, cfg, C.mods, models{m}, cfg.seed + r);
                end
            end
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
        barh(v(1:nT)); set(gca,'YDir','reverse','FontSize',8);
        yticks(1:nT); yticklabels(rs.impNames(o(1:nT)));
        set(gca,'TickLabelInterpreter','none');
        xlabel('Importance (gain, normalised)');
        title(sprintf('Iteration %d', iters(a)), 'FontWeight','bold');
        S{end+1} = sprintf('Iteration %d top features (%s, set A): %s', iters(a), cfg.primaryModel, ...
            strjoin(rs.impNames(o(1:min(6,numel(o)))), ', ')); %#ok<AGROW>
    end
    sgtitle(sprintf('Within-iteration feature importance (%s)', cfg.primaryModel), 'FontWeight','bold');
    saveBothFormats(fig, C.figDir, 'importance_by_iteration');

    % confusion matrices per stratum
    fig2 = figure('Position',[20 20 480*nI 520],'Color','w');
    for a = 1:nI
        rs = confStore{a};
        hp = uipanel(fig2,'Position',[(a-1)/nI 0 1/nI 0.93],'BorderType','none','BackgroundColor','w');
        if isempty(rs); continue; end
        present = unique([rs.yte; rs.pred]);
        yt = categorical(C.cats(rs.yte), C.cats(present));
        yp = categorical(C.cats(rs.pred), C.cats(present));
        cm = confusionchart(hp, yt, yp);
        cm.Title = sprintf('Iteration %d: acc %.1f%% (baseline %.1f%%)', iters(a), rs.acc, rs.baseTest);
        cm.RowSummary = 'row-normalized'; cm.FontSize = 8;
    end
    annotation(fig2,'textbox',[0 0.94 1 0.06],'String','Confusion matrices within each sampling iteration', ...
        'FontSize',12,'FontWeight','bold','EdgeColor','none','HorizontalAlignment','center');
    saveBothFormats(fig2, C.figDir, 'confusion_by_iteration');
end

%% ------------------------------------------------------------
function [S, T] = secDecisive(C)
% ML-e: pairwise decisive test inside one iteration + linear boundary figure.
    S = {}; T = table();
    cfg = C.cfg; models = cfg.models; nM = numel(models); nR = cfg.nRepeats;
    tests = cfg.decisiveTests;
    rowsT = {};
    [x1, xl1] = mechCoord(C, 1);  [x2, xl2] = mechCoord(C, 2);

    for t = 1:size(tests,1)
        itr = tests{t,1}; cA = string(tests{t,2}); cB = string(tests{t,3});
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
        for r = 1:nR
            sp = stratSplit(yb, cfg.testFrac, cfg.seed + 200 + r);
            for m = 1:nM
                ra = evalFS(fsAll,  yb, sp.tr, sp.te, cfg, C.mods, models{m}, cfg.seed + r);
                rm = evalFS(fsMech, yb, sp.tr, sp.te, cfg, C.mods, models{m}, cfg.seed + r);
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
        cOrange = [0.90 0.60 0.00]; cBlue = [0.00 0.45 0.70];
        fig = figure('Position',[60 60 760 640],'Color','w');
        hold on;
        scatter(x1(rows(yb==1)), x2(rows(yb==1)), 14, cOrange, 'filled', 'MarkerFaceAlpha', 0.45);
        scatter(x1(rows(yb==2)), x2(rows(yb==2)), 14, cBlue,   'filled', 'MarkerFaceAlpha', 0.45);
        w = svm1.Beta; b0 = svm1.Bias;
        xr = xlim;
        if abs(w(2)) > 1e-12
            plot(xr, -(w(1)*xr + b0)/w(2), 'k-', 'LineWidth', 2);
        else
            xline(-b0/w(1), 'k-', 'LineWidth', 2);
        end
        xlim(xr);
        xlabel(xl1, 'Interpreter','none'); ylabel(xl2, 'Interpreter','none');
        legend({sprintf('Class %s (n=%d)', cA, nA), sprintf('Class %s (n=%d)', cB, nB), 'Linear boundary'}, ...
            'Location','best');
        title(sprintf('Iteration %d, class %s vs %s\nlinear boundary %.1f%% | boosted trees (2 feats) %.1f%% | baseline %.1f%%', ...
            itr, cA, cB, mean(accLin), mean(accMech(strcmp(models,cfg.primaryModel),:)), mean(base)), ...
            'FontWeight','bold');
        box on; grid on;
        saveBothFormats(fig, C.figDir, ['mechanistic_scatter_' tag]);
    end
    if ~isempty(rowsT)
        T = cell2table(rowsT, 'VariableNames', {'Test','Model','Features','Acc_mean','Acc_sd','Baseline'});
        writetable(T, fullfile(C.outDir,'decisive_tests.xlsx'));
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

    % explicit audit of the pair named in cfg.confusionCheck
    ia = find(strcmp(C.cats, cfg.confusionCheck{1})); ib = find(strcmp(C.cats, cfg.confusionCheck{2}));
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
        S{end+1} = sprintf('Pair %s/%s not both present in this label scheme.', cfg.confusionCheck{1}, cfg.confusionCheck{2});
    end

    fig = figure('Position',[60 60 760 700],'Color','w');
    yt = categorical(C.cats(res.yte), C.cats);
    yp = categorical(C.cats(res.pred), C.cats);
    ch = confusionchart(fig, yt, yp);
    ch.Title = sprintf('Pooled %s: acc %.1f%% | majority baseline %.1f%%', cfg.primaryModel, res.acc, res.baseTest);
    ch.RowSummary = 'row-normalized'; ch.ColumnSummary = 'column-normalized'; ch.FontSize = 9;
    saveBothFormats(fig, C.figDir, 'confusion_pooled');

    out = struct('cm', cm, 'acc', res.acc, 'baseline', res.baseTest);
end

%% ------------------------------------------------------------
function S = secSmallMultiples(C)
% Figure suggestion: one panel per class (accent colour) against the grey ensemble.
    S = {};
    [x, xl] = mechCoord(C, 1);  [y, yl] = mechCoord(C, 2);
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
        scatter(x(ic), y(ic), 8, accent, 'filled', 'MarkerFaceAlpha', 0.65);
        title(sprintf('%s (n=%d)', C.className{c}, nnz(ic)), 'Interpreter','none', 'FontSize', 9);
        set(gca,'FontSize',8); box on;
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
        sp = stratSplit(yg, cfg.testFrac, cfg.seed + 300 + g);
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

function fs = makeFS(name, X, names, useJSD)
    fs = struct('name', name, 'X', X, 'names', {names}, 'jsd', useJSD);
end

function sp = stratSplit(key, testFrac, seed)
% Stratified train/test split on an arbitrary integer key.
    rng(seed);
    n = numel(key); isTe = false(n,1);
    u = unique(key);
    for i = 1:numel(u)
        idx = find(key == u(i)); m = numel(idx);
        if m < 2; continue; end                       % singletons stay in train
        perm = idx(randperm(m));
        nTe  = min(max(1, round(testFrac*m)), m-1);
        isTe(perm(1:nTe)) = true;
    end
    sp.te = find(isTe); sp.tr = find(~isTe);
end

function res = evalFS(fs, y, tr, te, cfg, mods, modelName, seed)
% Train on rows tr, score on rows te for one feature set / model.
    Xtr = fs.X(tr,:); Xte = fs.X(te,:); ytr = y(tr); yte = y(te);
    keepC = find(std(Xtr, 0, 1) > 0);                 % drop zero-variance columns
    % drop exact duplicate columns (keep the first): inside Iterations 3/4 many ratios are
    % identical (e.g. k1/k2 == k5/k6 because k1 = k5, k2 = k6), which would split importance
    if numel(keepC) > 1
        [~, ia] = unique(round(Xtr(:,keepC), 9)', 'rows', 'stable');
        keepC = keepC(sort(ia));
    end
    if fs.jsd && numel(unique(ytr)) > 1 && ~isempty(keepC)
        selLoc = selectJSDcols(Xtr(:,keepC), ytr, cfg.nBins);
        cols = keepC(selLoc);
    else
        cols = keepC;
    end
    [pred, imp] = fitPredictPy(modelName, Xtr(:,cols), ytr, Xte(:,cols), cfg, mods, seed);
    res = scoreRun(yte, pred, ytr);
    res.nFeat = numel(cols);
    res.imp = imp;
    res.impNames = fs.names(cols);
    res.pred = pred;
    res.yte = yte;
end

function [pred, imp] = fitPredictPy(modelName, Xtr, ytr, Xte, cfg, mods, seed)
% Weighted sklearn / LightGBM fit via the MATLAB-Python interface.
    [cls, ~, yi] = unique(ytr(:));
    nTe = size(Xte,1);
    if numel(cls) < 2 || size(Xtr,2) == 0
        pred = repmat(mode(ytr(:)), nTe, 1); imp = zeros(1, size(Xtr,2)); return;
    end
    cnt = accumarray(yi, 1);
    if cfg.useClassWeights
        sw = (numel(yi) ./ (numel(cls) .* cnt));  sw = sw(yi);
    else
        sw = ones(numel(yi), 1);
    end
    np = mods.np;
    % MATLAB hands an n-by-1 matrix (or a 1-by-p row) to Python as a 1-D array; sklearn needs
    % 2-D. Force the shape explicitly (happens when a feature set keeps only ONE column).
    Xtr_np = np.array(double(Xtr));
    Xtr_np = Xtr_np.reshape(int32(size(Xtr,1)), int32(size(Xtr,2)));
    Xte_np = np.array(double(Xte));
    Xte_np = Xte_np.reshape(int32(size(Xte,1)), int32(size(Xte,2)));
    y_np   = np.array(double(yi - 1));   y_np  = y_np.ravel();
    sw_np  = np.array(double(sw));       sw_np = sw_np.ravel();
    switch modelName
        case 'LightGBM'
            clf = mods.lgb.LGBMClassifier(pyargs('n_estimators', int32(cfg.lgbRounds), ...
                'num_leaves', int32(31), 'learning_rate', 0.1, 'min_child_samples', int32(3), ...
                'importance_type', 'gain', 'random_state', int32(seed), 'verbose', int32(-1)));
        case 'RandomForest'
            clf = mods.sk_ens.RandomForestClassifier(pyargs('n_estimators', int32(cfg.rfTrees), ...
                'max_features', 0.3, 'min_samples_leaf', int32(2), ...
                'random_state', int32(seed), 'n_jobs', int32(-1)));
        case 'NaiveBayes'
            clf = mods.sk_nb.GaussianNB();
        otherwise
            error('Unknown model %s', modelName);
    end
    clf.fit(Xtr_np, y_np, pyargs('sample_weight', sw_np));
    p = double(clf.predict(Xte_np));
    pred = cls(p(:) + 1);
    if strcmp(modelName, 'NaiveBayes')
        imp = zeros(1, size(Xtr,2));
    else
        imp = double(py.array.array('d', clf.feature_importances_));
        imp = imp(:)';  s = sum(imp);  if s > 0; imp = imp / s; end
    end
end

function m = scoreRun(yte, pred, ytr)
% Accuracy with majority baselines, balanced accuracy, macro-F1 and kappa.
    yte = yte(:); pred = pred(:);
    m.acc = mean(yte == pred) * 100;
    cls = unique(yte);
    rec = zeros(numel(cls),1); f1 = zeros(numel(cls),1);
    for i = 1:numel(cls)
        c = cls(i);
        tp = sum(pred==c & yte==c); fp = sum(pred==c & yte~=c); fn = sum(pred~=c & yte==c);
        pr = tp / max(tp+fp, 1); re = tp / max(tp+fn, 1);
        rec(i) = re; f1(i) = 2*pr*re / max(pr+re, eps);
    end
    m.balacc  = mean(rec) * 100;
    m.macroF1 = mean(f1) * 100;
    m.baseTest  = majorityShare(yte);                       % share of the commonest class in the test set
    m.baseTrain = mean(yte == mode(ytr(:))) * 100;          % always predict the TRAIN majority class
    allc = unique([yte; pred]);
    pe = 0;
    for i = 1:numel(allc); pe = pe + mean(yte==allc(i)) * mean(pred==allc(i)); end
    po = mean(yte == pred);
    m.kappa = (po - pe) / max(1 - pe, eps);
end

function b = majorityShare(y)
    [~,~,g] = unique(y(:));
    b = max(accumarray(g, 1)) / numel(y) * 100;
end

function sel = selectJSDcols(X, y, nBins)
% Original pipeline selection: per class pair JSD scores + strict
% broken-stick elbow; union over pairs. Returns column indices into X.
    cls = unique(y(:)); nC = numel(cls); nF = size(X,2);
    if nC < 2; sel = 1:nF; return; end
    keep = false(1, nF);
    pairs = nchoosek(1:nC, 2);
    for cp = 1:size(pairs,1)
        XA = X(y == cls(pairs(cp,1)), :);
        XB = X(y == cls(pairs(cp,2)), :);
        sc = zeros(1, nF);
        for j = 1:nF
            sc(j) = jsdScore(XA(:,j), XB(:,j), nBins);
        end
        [~, ord] = sort(sc, 'descend');
        ne = brokenStickElbow(sc);
        keep(ord(1:ne)) = true;
    end
    sel = find(keep);
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

function jsd = jsdScore(a, b, nBins)
% Jensen-Shannon Divergence (base-e) between two sample vectors.
    allV  = [a(:); b(:)];
    edges = linspace(min(allV)-eps, max(allV)+eps, nBins+1);
    pA    = histcounts(a, edges, 'Normalization','probability') + eps;
    pB    = histcounts(b, edges, 'Normalization','probability') + eps;
    pA    = pA / sum(pA);
    pB    = pB / sum(pB);
    M     = 0.5*(pA + pB);
    kl_pm = sum(pA .* log(pA ./ M));
    kl_qm = sum(pB .* log(pB ./ M));
    jsd   = max(0, 0.5*kl_pm + 0.5*kl_qm);
end

function nElbow = brokenStickElbow(values)
% STRICT broken-stick elbow (first crossing), minimum 1 feature.
    values = sort(values(:)', 'descend');
    n      = numel(values);
    if n == 0;  nElbow = 1; return; end
    total  = sum(values);
    if total <= 0; nElbow = 1; return; end
    proportions = values / total;
    expected    = brokenStickExpected(n);
    first_below = find(proportions < expected, 1, 'first');
    if isempty(first_below)
        nElbow = n;
    else
        nElbow = max(first_below - 1, 1);
    end
end

function expProp = brokenStickExpected(n)
% Broken-stick expected proportions for n components (row vector).
    if n <= 0; expProp = zeros(0,1); return; end
    expProp = zeros(1, n);
    for k = 1:n
        expProp(k) = (1/n) * sum(1./(k:n));
    end
end

function saveBothFormats(fig, figDir, name)
    exportgraphics(fig, fullfile(figDir,[name,'.png']), 'Resolution', 300);
    savefig(fig, fullfile(figDir,[name,'.fig']));
end