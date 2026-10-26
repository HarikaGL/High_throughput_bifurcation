%% ML_v3_inference_time.m
% Time per parameter set for the trained LightGBM classifier to assign a class
% (main text: "The class of a construct can be assigned in milliseconds", Discussion
% and Methods, Software, hardware and computational cost).
%
%   1. Builds feature set A (14 log-parameters + 91 pairwise log-ratios) for the 39,968
%      sets used for machine learning.
%   2. Trains LightGBM with the settings of ML_v3_10class_parallel.m (300 rounds, 63 leaves,
%      learning rate 0.1, inverse-frequency class weights) on a stratified 80% split.
%   3. Times the prediction on the held-out 20%, including the computation of the 105
%      features from the 14 raw parameters:
%        - in one batch (all held-out sets at once), time per set;
%        - one set at a time (500 sets), median time per set.
%   One thread is used throughout, so the times are per core.
%
% Input : Ensemble_iter<i>_v3_1e6/Excel_Files/ML_Ready_Dataset_v3_iter<i>.csv
% Output: ML_v3_results_parallel/inference_time.txt
% Needs : Python with numpy and lightgbm, reachable through pyenv (as for
%         ML_v3_10class_parallel.m).

clear;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
addpath(fullfile(HERE, 'support'));
OUT = fullfile(HERE, 'ML_v3_results_parallel'); if ~exist(OUT, 'dir'), mkdir(OUT); end

E = load_ensemble_v3(HERE);
keep = E.use;
P = E.p(keep,:); y = E.class(keep); it = E.iteration(keep);
pairs = nchoosek(1:14, 2);
feats = @(Q) [log10(Q), log10(Q(:,pairs(:,1))) - log10(Q(:,pairs(:,2)))];

% stratified 80/20 split by iteration and class
rng(42);
key = it*100 + y; te = false(size(y));
for k = unique(key)'
    idx = find(key == k); idx = idx(randperm(numel(idx)));
    te(idx(1:round(0.2*numel(idx)))) = true;
end
tr = ~te;

% LightGBM through the MATLAB-Python interface
np  = py.importlib.import_module('numpy');
lgb = py.importlib.import_module('lightgbm');
toNp = @(M) np.array(M(:)').reshape(int32(size(M,2)), int32(size(M,1))).T;
[cls, ~, yi] = unique(y(tr));
cnt = accumarray(yi, 1); w = numel(yi) ./ (numel(cls) .* cnt);
clf = lgb.LGBMClassifier(pyargs('n_estimators', int32(300), 'num_leaves', int32(63), ...
      'learning_rate', 0.1, 'min_child_samples', int32(3), 'random_state', int32(1), ...
      'verbose', int32(-1), 'n_jobs', int32(1)));
tFit = tic;
clf.fit(toNp(feats(P(tr,:))), np.array(yi' - 1), pyargs('sample_weight', np.array(w(yi)')));
tFit = toc(tFit);

Pte = P(te,:); nTe = size(Pte, 1);
pred = double(py.array.array('d', np.ravel(clf.predict(toNp(feats(Pte))))));
acc = mean(cls(pred(:) + 1) == y(te)) * 100;

% batch: all held-out sets at once (second run reported, after warm-up)
for rep = 1:2
    tB = tic; clf.predict_proba(toNp(feats(Pte))); tB = toc(tB);
end
% one set at a time
nS = 500; tS = zeros(nS, 1);
for i = 1:nS
    t1 = tic; clf.predict_proba(toNp(feats(Pte(i,:)))); tS(i) = toc(t1);
end

fid = fopen(fullfile(OUT, 'inference_time.txt'), 'w');
for f = [1 fid]
    fprintf(f, 'LightGBM, feature set A (105 features), one thread\n');
    fprintf(f, 'training: %d sets, %.1f s | held-out accuracy %.2f%% (%d sets)\n', sum(tr), tFit, acc, nTe);
    fprintf(f, 'batch prediction of %d sets: %.3f s in total, %.3f ms per set\n', nTe, tB, tB/nTe*1e3);
    fprintf(f, 'single-set prediction: median %.3f ms per set (90th percentile %.3f ms, %d sets)\n', ...
            median(tS)*1e3, prctile(tS, 90)*1e3, nS);
    fprintf(f, 'times include the computation of the 105 features from the 14 parameters\n');
end
fclose(fid);
