function varargout = ml_v3_core(action, varargin)
%ML_V3_CORE  Model fitting / scoring core of ML_v3_10class_parallel.m (runs on the client
%   and on every parallel worker; each process imports its own Python modules once).
%
%   M   = ml_v3_core('mods')                         Python modules (persistent per process)
%   res = ml_v3_core('evalFS', fs, y, tr, te, cfg, model, seed)
%   out = ml_v3_core('job', job, D, cfg)             one job of the parallel batch
%   cols= ml_v3_core('selectCols', fs, Xtr, ytr, cfg)
%   b   = ml_v3_core('majorityShare', y)
%
%   Job kinds: 'eval'  fit on job.tr, score on job.te
%              'cv'    'eval' + stratified k-fold CV inside job.tr (acc, macro-F1, weighted-F1)
%              'curve' fit on job.tr; ROC on job.te; convergence curves on job.va and a
%                      subsample of job.tr
%   Models are only ever fitted on the training rows. cfg.nThreads = threads per model
%   (1 on parallel workers, -1 = all cores on the client).
switch action
    case 'mods',          varargout{1} = getMods();
    case 'evalFS',        varargout{1} = evalFS(varargin{:});
    case 'job',           varargout{1} = runJob(varargin{:});
    case 'selectCols',    varargout{1} = selectCols(varargin{:});
    case 'majorityShare', varargout{1} = majorityShare(varargin{:});
    otherwise, error('ml_v3_core: unknown action %s', action);
end
end

%% =====================================================================================
function M = getMods()
persistent P
if isempty(P)
    P = struct();
    P.np     = py.importlib.import_module('numpy');
    P.lgb    = py.importlib.import_module('lightgbm');
    P.sk_ens = py.importlib.import_module('sklearn.ensemble');
    P.sk_nb  = py.importlib.import_module('sklearn.naive_bayes');
    P.xgb = []; P.cb = [];
    try, P.xgb = py.importlib.import_module('xgboost');  catch, end
    try, P.cb  = py.importlib.import_module('catboost'); catch, end
end
M = P;
end

%% =====================================================================================
function out = runJob(job, D, cfg)
% D: struct with FS (cell of feature sets), yCode, itIdx
if job.fsId > 0
    fs = D.FS{job.fsId};
else
    fs = struct('X', job.X, 'names', {job.names}, 'jsd', job.jsd);
end
if ~isempty(job.y), y = job.y; else, y = D.(job.yName); end
t0 = tic;
switch job.kind
    case 'eval'
        out = evalFS(fs, y, job.tr, job.te, cfg, job.model, job.seed);
    case 'cv'
        out = evalFS(fs, y, job.tr, job.te, cfg, job.model, job.seed);
        k = cfg.cvFolds; tr = job.tr(:); ytr = y(tr);
        fold = zeros(numel(tr), 1); s = RandStream('mt19937ar', 'Seed', job.seed + 7);
        for c = unique(ytr)'
            ii = find(ytr == c); ii = ii(randperm(s, numel(ii)));
            fold(ii) = mod(0:numel(ii)-1, k) + 1;
        end
        out.cvA = zeros(1,k); out.cvM = zeros(1,k); out.cvW = zeros(1,k);
        for q = 1:k
            rq = evalFS(fs, y, tr(fold ~= q), tr(fold == q), cfg, job.model, job.seed + q);
            out.cvA(q) = rq.acc; out.cvM(q) = rq.macroF1; out.cvW(q) = rq.weightedF1;
        end
    case 'curve'
        out = curveJob(fs, y, job, cfg);
    otherwise
        error('unknown job kind %s', job.kind);
end
out.seconds = toc(t0);
end

%% =====================================================================================
function res = evalFS(fs, y, tr, te, cfg, modelName, seed)
% Train on rows tr, score on rows te for one feature set / model.
Xtr = fs.X(tr,:); Xte = fs.X(te,:); ytr = y(tr); yte = y(te);
cols = selectCols(fs, Xtr, ytr, cfg);
[pred, imp] = fitPredictPy(modelName, Xtr(:,cols), ytr, Xte(:,cols), cfg, seed);
res = scoreRun(yte, pred, ytr);
res.nFeat = numel(cols);
res.imp = imp;
res.cols = cols;
res.impNames = fs.names(cols);
res.pred = pred;
res.yte = yte;
end

function cols = selectCols(fs, Xtr, ytr, cfg)
% columns used for training (computed from the TRAINING rows only): non-constant,
% de-duplicated, JSD + broken-stick selected for feature set D
keepC = find(std(Xtr, 0, 1) > 0);
if numel(keepC) > 1
    [~, ia] = unique(round(Xtr(:,keepC), 9)', 'rows', 'stable');
    keepC = keepC(sort(ia));
end
if fs.jsd && numel(unique(ytr)) > 1 && ~isempty(keepC)
    cols = keepC(selectJSDcols(Xtr(:,keepC), ytr, cfg.nBins));
else
    cols = keepC;
end
end

function [pred, imp, clf, cls] = fitPredictPy(modelName, Xtr, ytr, Xte, cfg, seed)
% Weighted fit (inverse-frequency sample weights) through the MATLAB-Python interface.
M = getMods(); np = M.np;
[cls, ~, yi] = unique(ytr(:));
nTe = size(Xte,1); clf = [];
if numel(cls) < 2 || size(Xtr,2) == 0
    pred = repmat(mode(ytr(:)), nTe, 1); imp = zeros(1, size(Xtr,2)); return;
end
cnt = accumarray(yi, 1);
if cfg.useClassWeights, sw = numel(yi) ./ (numel(cls) .* cnt); sw = sw(yi);
else, sw = ones(numel(yi), 1); end
if strcmp(modelName, 'MostFrequent')
    [~, im] = max(cnt); pred = repmat(cls(im), nTe, 1); imp = zeros(1, size(Xtr,2)); return;
end
nt = int32(cfg.nThreads);
Xtr_np = np.array(double(Xtr)); Xtr_np = Xtr_np.reshape(int32(size(Xtr,1)), int32(size(Xtr,2)));
Xte_np = np.array(double(Xte)); Xte_np = Xte_np.reshape(int32(size(Xte,1)), int32(size(Xte,2)));
y_np  = np.array(double(yi - 1)); y_np  = y_np.ravel();
sw_np = np.array(double(sw));     sw_np = sw_np.ravel();
switch modelName
    case 'LightGBM'
        clf = M.lgb.LGBMClassifier(pyargs('n_estimators', int32(cfg.lgbRounds), ...
            'num_leaves', int32(cfg.lgbLeaves), 'learning_rate', 0.1, 'min_child_samples', int32(3), ...
            'importance_type', 'gain', 'random_state', int32(seed), 'verbose', int32(-1), 'n_jobs', nt));
    case 'RandomForest'
        clf = M.sk_ens.RandomForestClassifier(pyargs('n_estimators', int32(cfg.rfTrees), ...
            'max_features', 0.3, 'min_samples_leaf', int32(2), 'random_state', int32(seed), 'n_jobs', nt));
    case 'NaiveBayes'
        clf = M.sk_nb.GaussianNB();
    case 'XGBoost'
        clf = M.xgb.XGBClassifier(pyargs('n_estimators', int32(cfg.xgbRounds), ...
            'max_depth', int32(cfg.xgbDepth), 'learning_rate', 0.1, 'eval_metric', 'mlogloss', ...
            'tree_method', 'hist', 'random_state', int32(seed), 'n_jobs', nt));
    case 'CatBoost'
        clf = M.cb.CatBoostClassifier(pyargs('iterations', int32(cfg.catIter), ...
            'depth', int32(cfg.catDepth), 'learning_rate', 0.1, 'random_seed', int32(seed), ...
            'verbose', int32(0), 'thread_count', nt, 'allow_writing_files', false));
    case 'GradBoost'
        if cfg.gbFast
            clf = M.sk_ens.HistGradientBoostingClassifier(pyargs('max_iter', int32(cfg.gbRounds), ...
                'max_depth', int32(cfg.gbDepth), 'learning_rate', 0.1, 'random_state', int32(seed)));
        else
            clf = M.sk_ens.GradientBoostingClassifier(pyargs('n_estimators', int32(cfg.gbRounds), ...
                'max_depth', int32(cfg.gbDepth), 'learning_rate', 0.1, 'random_state', int32(seed)));
        end
    otherwise
        error('Unknown model %s', modelName);
end
clf.fit(Xtr_np, y_np, pyargs('sample_weight', sw_np));
p = double(py.array.array('d', np.ravel(clf.predict(Xte_np))));
pred = cls(p(:) + 1);
if strcmp(modelName, 'NaiveBayes') || (strcmp(modelName, 'GradBoost') && cfg.gbFast)
    imp = zeros(1, size(Xtr,2));
else
    imp = double(py.array.array('d', np.ravel(py.getattr(clf, 'feature_importances_'))));
    imp = imp(:)'; s = sum(imp); if s > 0, imp = imp / s; end
end
end

%% =====================================================================================
function out = curveJob(fs, y, job, cfg)
% ROC (test rows) and convergence curves (validation rows + training subsample)
tr = job.tr(:); te = job.te(:); va = job.va(:); if isempty(va), va = te; end
Xtr = fs.X(tr,:); ytr = y(tr);
cols = selectCols(fs, Xtr, ytr, cfg);
[~, ~, clf, cls] = fitPredictPy(job.model, Xtr(:,cols), ytr, fs.X(te,cols), cfg, job.seed);
out = struct('model', job.model, 'cls', cls, 'auc', [], 'fpr', {{}}, 'tpr', {{}}, 'steps', [], ...
    'llTr', [], 'llVa', [], 'erTr', [], 'erVa', [], 'errC', []);
if isempty(clf), return; end
M = getMods(); np = M.np; K = numel(cls);
yte = y(te);
Pte = probaPy(clf, fs.X(te,cols), np, K);
out.auc = NaN(1,K); out.fpr = cell(1,K); out.tpr = cell(1,K);
for c = 1:K
    pos = (yte == cls(c));
    if ~any(pos) || all(pos), continue; end
    [f, t, ~, out.auc(c)] = perfcurve(pos, Pte(:,c), true);
    keep = unique(round(linspace(1, numel(f), min(400, numel(f)))));   % thin for storage
    out.fpr{c} = f(keep); out.tpr{c} = t(keep);
end
s = RandStream('mt19937ar', 'Seed', job.seed);
trSub = randperm(s, numel(tr), min(numel(tr), cfg.curveTrainMax));
[steps, PvaS, PtrS] = stagedProbaPy(clf, job.model, fs.X(va,cols), Xtr(trSub,cols), np, K, cfg);
if isempty(steps), return; end
[~, yvI] = ismember(y(va), cls); [~, ytI] = ismember(ytr(trSub), cls);
nS = numel(steps);
out.steps = steps; out.llTr = zeros(1,nS); out.llVa = zeros(1,nS);
out.erTr = zeros(1,nS); out.erVa = zeros(1,nS); out.errC = NaN(K, nS);
for q = 1:nS
    [out.llVa(q), out.erVa(q)] = lossErr(PvaS{q}, yvI);
    [out.llTr(q), out.erTr(q)] = lossErr(PtrS{q}, ytI);
    [~, pr] = max(PvaS{q}, [], 2);
    for c = 1:K
        if any(yvI == c), out.errC(c,q) = mean(pr(yvI == c) ~= c) * 100; end
    end
end
end

function [ll, er] = lossErr(P, yI)
P = max(min(P, 1 - 1e-15), 1e-15); P = P ./ sum(P, 2);
ll = -mean(log(P(sub2ind(size(P), (1:numel(yI))', yI(:)))));
[~, pr] = max(P, [], 2); er = mean(pr(:) ~= yI(:)) * 100;
end

function P = probaPy(clf, X, np, K)
Xn = np.array(double(X)); Xn = Xn.reshape(int32(size(X,1)), int32(size(X,2)));
P = pyToMat(clf.predict_proba(Xn), np, size(X,1), K);
end

function P = pyToMat(pyArr, np, n, K)
v = double(py.array.array('d', np.ravel(pyArr)));
P = reshape(v, K, n)';
end

function [steps, PvaS, PtrS] = stagedProbaPy(clf, model, Xva, Xtr, np, K, cfg)
steps = []; PvaS = {}; PtrS = {};
Xe = np.array(double(Xva)); Xe = Xe.reshape(int32(size(Xva,1)), int32(size(Xva,2)));
Xr = np.array(double(Xtr)); Xr = Xr.reshape(int32(size(Xtr,1)), int32(size(Xtr,2)));
switch model
    case 'LightGBM',     nMax = cfg.lgbRounds;
    case 'XGBoost',      nMax = cfg.xgbRounds;
    case 'CatBoost',     nMax = double(py.getattr(clf, 'tree_count_'));
    case 'GradBoost',    nMax = cfg.gbRounds;
    case 'RandomForest', nMax = cfg.rfTrees;
    otherwise,           return;
end
steps = unique(round(linspace(1, nMax, cfg.curvePoints)));
switch model
    case {'LightGBM','XGBoost','CatBoost'}
        for q = 1:numel(steps)
            t = int32(steps(q));
            switch model
                case 'LightGBM', a = pyargs('num_iteration', t);
                case 'XGBoost',  a = pyargs('iteration_range', py.tuple({int32(0), t}));
                case 'CatBoost', a = pyargs('ntree_end', t);
            end
            PvaS{q} = pyToMat(clf.predict_proba(Xe, a), np, size(Xva,1), K); %#ok<AGROW>
            PtrS{q} = pyToMat(clf.predict_proba(Xr, a), np, size(Xtr,1), K); %#ok<AGROW>
        end
    case 'GradBoost'
        ge = cell(py.list(clf.staged_predict_proba(Xe))); gr = cell(py.list(clf.staged_predict_proba(Xr)));
        steps = steps(steps <= numel(ge));
        for q = 1:numel(steps)
            PvaS{q} = pyToMat(ge{steps(q)}, np, size(Xva,1), K); %#ok<AGROW>
            PtrS{q} = pyToMat(gr{steps(q)}, np, size(Xtr,1), K); %#ok<AGROW>
        end
    case 'RandomForest'
        est = cell(py.list(py.getattr(clf, 'estimators_'))); nT = numel(est);
        steps = steps(steps <= nT);
        Ae = zeros(size(Xva,1), K); Ar = zeros(size(Xtr,1), K); q = 1;
        for t = 1:nT
            Ae = Ae + pyToMat(est{t}.predict_proba(Xe), np, size(Xva,1), K);
            Ar = Ar + pyToMat(est{t}.predict_proba(Xr), np, size(Xtr,1), K);
            if q <= numel(steps) && t == steps(q)
                PvaS{q} = Ae / t; PtrS{q} = Ar / t; q = q + 1; %#ok<AGROW>
            end
        end
end
end

%% =====================================================================================
function m = scoreRun(yte, pred, ytr)
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
m.baseTest  = majorityShare(yte);
m.baseTrain = mean(yte == mode(ytr(:))) * 100;
allc = unique([yte; pred]); pe = 0;
for i = 1:numel(allc), pe = pe + mean(yte==allc(i)) * mean(pred==allc(i)); end
m.kappa = (mean(yte == pred) - pe) / max(1 - pe, eps);
w = arrayfun(@(c) mean(yte == c), cls);
m.weightedF1 = sum(w(:) .* f1(:)) * 100;
end

function b = majorityShare(y)
[~,~,g] = unique(y(:));
b = max(accumarray(g, 1)) / numel(y) * 100;
end

function sel = selectJSDcols(X, y, nBins)
cls = unique(y(:)); nC = numel(cls); nF = size(X,2);
if nC < 2, sel = 1:nF; return; end
keep = false(1, nF);
pairs = nchoosek(1:nC, 2);
for cp = 1:size(pairs,1)
    XA = X(y == cls(pairs(cp,1)), :); XB = X(y == cls(pairs(cp,2)), :);
    sc = zeros(1, nF);
    for j = 1:nF, sc(j) = jsdScore(XA(:,j), XB(:,j), nBins); end
    [~, ord] = sort(sc, 'descend');
    keep(ord(1:brokenStickElbow(sc))) = true;
end
sel = find(keep);
end

function jsd = jsdScore(a, b, nBins)
allV  = [a(:); b(:)];
edges = linspace(min(allV)-eps, max(allV)+eps, nBins+1);
pA = histcounts(a, edges, 'Normalization','probability') + eps;
pB = histcounts(b, edges, 'Normalization','probability') + eps;
pA = pA / sum(pA); pB = pB / sum(pB);
Mx = 0.5*(pA + pB);
jsd = max(0, 0.5*sum(pA .* log(pA ./ Mx)) + 0.5*sum(pB .* log(pB ./ Mx)));
end

function nElbow = brokenStickElbow(values)
values = sort(values(:)', 'descend'); n = numel(values);
if n == 0 || sum(values) <= 0, nElbow = 1; return; end
proportions = values / sum(values);
expected = arrayfun(@(k) (1/n) * sum(1./(k:n)), 1:n);
first_below = find(proportions < expected, 1, 'first');
if isempty(first_below), nElbow = n; else, nElbow = max(first_below - 1, 1); end
end
