%% ============================================================
%  MULTI-CLASS CLASSIFICATION (1–7) — FULL PIPELINE
%
%  Feature filtering (JSD score + Broken-Stick Elbow method):
%    - For every pair of classes (21 pairs), score each
%      individual feature and each ratio feature using JSD.
%    - Apply the broken-stick elbow method per class pair to
%      identify the naturally separating features.
%    - STRICTER ELBOW: use the FIRST crossing point where the
%      observed proportion drops BELOW the broken-stick expected,
%      rather than keeping all features above the line.
%      This meaningfully reduces the selected feature count.
%    - Keep a feature / ratio only if it survives the elbow
%      in AT LEAST ONE class-pair comparison.
%    - Final feature set = union of elbow-selected features
%      across all class-pair comparisons.
%    - Broken-stick elbow plots shown for BOTH individual
%      features AND ratio features.
%
%  All confusion matrices   → one figure (subplots)
%  All top-20 importance    → one figure (subplots)
%  All class-wise perf      → one figure (subplots)
%
%  *** PYTHON-ONLY MODELS (via MATLAB–Python interface) ***
%  Models:
%    1. Random Forest       (sklearn, weighted)
%    2. XGBoost             (xgboost, weighted)
%    3. LightGBM            (lightgbm, weighted)
%    4. CatBoost            (catboost, weighted)
%    5. Gradient Boosting   (sklearn, weighted)
%
%  Splits: 70% train | 20% test | 10% validation (per class)
% =============================================================

clear; clc; close all;

%% ============================================================
%% 1. SETUP
%% ============================================================

outDir = 'ML_JSD_BrokenStick_PythonOnly_final';
figDir = fullfile(outDir,'figures');
if ~exist(outDir,'dir'); mkdir(outDir); end
if ~exist(figDir,'dir'); mkdir(figDir); end

saveFig = @(name) saveBothFormats(gcf, figDir, name);

nBins = 30;   % histogram bins for JSD estimation

%% ============================================================
%% 2. LOAD DATA
%% ============================================================

T = readtable('ML_Ready_Dataset_ALL_1to7_10k_full.xlsx');
T = standardizeMissing(T, {'N/A','NA','na','','NONE'});

keepCols  = {'KTOT','PTOT','k1','k11','k2','k3','k33','k4','k5',...
             'k55','k6','k7','k77','k8','Class_Number_7'};
availCols = T.Properties.VariableNames;
keepCols  = keepCols(ismember(keepCols, availCols));
T         = T(:, keepCols);
T         = rmmissing(T);
targetVar = 'Class_Number_7';

writetable(T, fullfile(outDir,'Dataset_KTOT_to_k8.xlsx'));
fprintf('Dataset: %d rows, %d cols\n', height(T), width(T));

%% ============================================================
%% 3. CLASS DISTRIBUTION
%% ============================================================

Y_all        = categorical(T.(targetVar));
classLabels  = categories(Y_all);
classCounts  = countcats(Y_all);
classPercent = 100 * classCounts / sum(classCounts);

figure;
bar(classCounts);
xticklabels(classLabels);
ylabel('Number of Samples');
title('Original Class Distribution');
for i = 1:length(classCounts)
    text(i, classCounts(i), sprintf('%d (%.1f%%)', classCounts(i), classPercent(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold');
end
saveFig('class_distribution_original');

%% ============================================================
%% 4. STRATIFIED SPLIT: 70% TRAIN | 20% TEST | 10% VALIDATION
%% ============================================================

rng(42);
trainIdx = []; testIdx = []; valIdx = [];
for c = 1:numel(classLabels)
    cls    = classLabels{c};
    clsIdx = find(Y_all == cls);
    n      = numel(clsIdx);
    perm   = clsIdx(randperm(n));
    nTrain = round(0.70 * n);
    nTest  = round(0.20 * n);
    trainIdx = [trainIdx; perm(1:nTrain)];
    testIdx  = [testIdx;  perm(nTrain+1 : nTrain+nTest)];
    valIdx   = [valIdx;   perm(nTrain+nTest+1 : end)];
end

T_train = T(trainIdx, :);
T_test  = T(testIdx,  :);
T_val   = T(valIdx,   :);

fprintf('Train: %d | Test: %d | Val: %d\n', height(T_train), height(T_test), height(T_val));
writetable(T_val, fullfile(outDir,'Validation_Dataset_10pct_perClass.xlsx'));

%% Pie chart of split sizes
figure;
splits = [height(T_train), height(T_test), height(T_val)];
labels = {'Train (70%)','Test (20%)','Validation (10%)'};
pie(splits, labels);
title(sprintf('Data Split: %d total samples', sum(splits)), 'FontWeight','bold');
saveFig('data_split_pie');

%% ============================================================
%% 5. BASE FEATURE MATRICES
%% ============================================================

Y_train = categorical(T_train.(targetVar));
Y_test  = categorical(T_test.(targetVar));
Y_val   = categorical(T_val.(targetVar));

T_train(:, targetVar) = [];
T_test(:,  targetVar) = [];
T_val(:,   targetVar) = [];

X_train      = table2array(T_train);
X_test       = table2array(T_test);
X_val        = table2array(T_val);
featureNames = T_train.Properties.VariableNames;
nFeat        = numel(featureNames);

%% ============================================================
%% 6. CORRELATION MATRIX
%% ============================================================

R = corr(X_train,'Rows','pairwise');
figure; imagesc(R); colorbar; axis square;
xticks(1:nFeat); yticks(1:nFeat);
xticklabels(featureNames); yticklabels(featureNames);
xtickangle(45);
title('Feature Correlation Matrix (Raw Features)');
saveFig('correlation_matrix');

%% ============================================================
%% 7A. JSD-BASED SEPARABILITY — ALL CLASS PAIRS
%%     STRICT Broken-Stick Elbow selection (first-crossing, min 1)
%%     Uses FIRST CROSSING POINT: elbow = index just before the
%%     observed proportion first falls below the broken-stick expected.
%%     At least 1 feature is always kept per class pair (even if the
%%     top feature itself falls below the broken-stick line).
%% ============================================================

fprintf('\n--- JSD Separability: Stricter Broken-Stick Elbow filtering ---\n');

nC         = numel(classLabels);
classPairs = nchoosek(1:nC, 2);
nCP        = size(classPairs,1);

featPairs = nchoosek(1:nFeat, 2);
nFP       = size(featPairs,1);

% Track which features ever survive the elbow across ANY class pair
indiv_everSelected = false(nFeat, 1);
ratio_everSelected = false(nFP,   1);

% Per-class-pair bookkeeping for logging
cpIndivSelected = cell(nCP, 1);
cpRatioSelected = cell(nCP, 1);
cpIndivCount    = zeros(nCP, 1);
cpRatioCount    = zeros(nCP, 1);

% Store all JSD score matrices for heat-maps and elbow plots
jsdIndivMatrix = zeros(nCP, nFeat);
jsdRatioMatrix = zeros(nCP, nFP);

fprintf('Scoring %d class pairs x (%d indiv + %d ratio) features  [JSD, nBins=%d]\n', ...
    nCP, nFeat, nFP, nBins);

for cp = 1:nCP
    cA = classPairs(cp,1);
    cB = classPairs(cp,2);
    XA = X_train(Y_train == classLabels{cA}, :);
    XB = X_train(Y_train == classLabels{cB}, :);

    %% --- Score individual features ---
    jsdIndiv = zeros(1, nFeat);
    for fi = 1:nFeat
        jsdIndiv(fi) = jsdScore(XA(:,fi), XB(:,fi), nBins);
    end
    jsdIndivMatrix(cp,:) = jsdIndiv;

    % Sort descending → stricter broken-stick elbow
    [sortedIndiv, sortIdxIndiv] = sort(jsdIndiv, 'descend');
    nElbowIndiv  = brokenStickElbow(sortedIndiv);   % first-crossing elbow
    if nElbowIndiv > 0
        selIndivLocal = sortIdxIndiv(1:nElbowIndiv);
    else
        selIndivLocal = [];
    end

    for fi = selIndivLocal
        indiv_everSelected(fi) = true;
    end
    cpIndivSelected{cp} = featureNames(selIndivLocal);
    cpIndivCount(cp)    = nElbowIndiv;

    %% --- Score ratio features ---
    jsdRatio = zeros(1, nFP);
    for fp = 1:nFP
        fi = featPairs(fp,1);
        fj = featPairs(fp,2);
        rA = XA(:,fi) ./ (XA(:,fj) + eps);
        rB = XB(:,fi) ./ (XB(:,fj) + eps);
        jsdRatio(fp) = jsdScore(rA, rB, nBins);
    end
    jsdRatioMatrix(cp,:) = jsdRatio;

    [sortedRatio, sortIdxRatio] = sort(jsdRatio, 'descend');
    nElbowRatio  = brokenStickElbow(sortedRatio);   % first-crossing elbow
    if nElbowRatio > 0
        selRatioLocal = sortIdxRatio(1:nElbowRatio);
    else
        selRatioLocal = [];
    end

    for fp = selRatioLocal
        ratio_everSelected(fp) = true;
    end
    cpRatioSelected{cp} = cellfun(@(a,b) sprintf('%s/%s',a,b), ...
        featureNames(featPairs(selRatioLocal,1)), ...
        featureNames(featPairs(selRatioLocal,2)), ...
        'UniformOutput', false);
    cpRatioCount(cp) = nElbowRatio;

    fprintf('  Cls %d vs Cls %d : %d indiv selected, %d ratios selected  (elbow)\n', ...
        cA, cB, nElbowIndiv, nElbowRatio);
end

%% --- Summary table ---
cpLabels = cell(nCP,1);
for cp = 1:nCP
    cpLabels{cp} = sprintf('Cls%d_vs_Cls%d', classPairs(cp,1), classPairs(cp,2));
end
PairSummary = table(cpLabels, cpIndivCount, cpRatioCount, ...
    'VariableNames',{'ClassPair','Elbow_Indiv','Elbow_Ratio'});
disp(PairSummary);
writetable(PairSummary, fullfile(outDir,'per_classpair_elbow_counts.xlsx'));

%% --- JSD heatmap: individual features ---
figure('Position',[100 100 max(600, nFeat*40+100) max(400, nCP*22+100)]);
imagesc(jsdIndivMatrix);
colorbar; colormap('hot');
xticks(1:nFeat); xticklabels(featureNames); xtickangle(45);
yticks(1:nCP);   yticklabels(cpLabels);
title('JSD Scores — Individual Features × Class Pairs','FontWeight','bold');
saveFig('jsd_heatmap_individual');

%% --- Final selected sets ---
selectedIndivIdx   = find(indiv_everSelected);
selectedIndivNames = featureNames(selectedIndivIdx);
nSelIndiv          = numel(selectedIndivIdx);

selectedRatioIdx   = find(ratio_everSelected);
nSelRatio          = numel(selectedRatioIdx);

ratioLabelsFull = cell(nFP,1);
for fp = 1:nFP
    ratioLabelsFull{fp} = sprintf('%s/%s', ...
        featureNames{featPairs(fp,1)}, featureNames{featPairs(fp,2)});
end

fprintf('\n=== FINAL SELECTION (JSD + Stricter Broken-Stick Elbow, union across all class pairs) ===\n');
fprintf('Individual features selected : %d / %d\n', nSelIndiv, nFeat);
fprintf('Ratio features selected      : %d / %d\n', nSelRatio, nFP);
fprintf('\nIndividual features:\n'); disp(selectedIndivNames);
fprintf('Ratio features:\n');       disp(ratioLabelsFull(selectedRatioIdx));

IndivSelTable = table(selectedIndivNames(:), ...
    'VariableNames',{'SelectedIndivFeature'});
RatioSelTable = table(ratioLabelsFull(selectedRatioIdx), ...
    'VariableNames',{'SelectedRatioFeature'});
writetable(IndivSelTable, fullfile(outDir,'selected_features_JSD.xlsx'),'Sheet','Individual');
writetable(RatioSelTable, fullfile(outDir,'selected_features_JSD.xlsx'),'Sheet','Ratios');

%% ============================================================
%% 7A-PLOT. BROKEN-STICK ELBOW — INDIVIDUAL FEATURES (per class pair)
%% ============================================================

nCPcols = 4;
nCProws = ceil(nCP / nCPcols);

figure('Position',[50 50 nCPcols*380 nCProws*300],'Color','w');
for cp = 1:nCP
    subplot(nCProws, nCPcols, cp);
    scores = sort(jsdIndivMatrix(cp,:), 'descend');
    n      = numel(scores);
    total  = sum(scores);
    if total > 0
        propScores = scores / total;
    else
        propScores = zeros(1,n);
    end
    expProp = brokenStickExpected(n);

    bar(propScores, 'FaceColor',[0.4 0.6 0.9],'EdgeColor','none'); hold on;
    plot(expProp, 'r-', 'LineWidth',1.5);

    nE = brokenStickElbow(scores);
    if nE > 0
        xline(nE + 0.5, 'k--', 'LineWidth',1.2);
        title(sprintf('%s  [k=%d]', cpLabels{cp}, nE), ...
            'FontSize',7,'FontWeight','bold','Interpreter','none');
    else
        title(cpLabels{cp},'FontSize',7,'FontWeight','bold','Interpreter','none');
    end
    xlabel('Feature rank','FontSize',6);
    ylabel('Prop. JSD','FontSize',6);
    set(gca,'FontSize',6);
    if cp == 1
        legend({'Actual','Broken stick','Elbow'},'FontSize',5,'Location','northeast');
    end
end
sgtitle('Broken-Stick Elbow — JSD INDIVIDUAL Features per Class Pair', ...
    'FontSize',12,'FontWeight','bold');
saveFig('broken_stick_elbow_indiv_all_pairs');

%% ============================================================
%% 7A-PLOT2. BROKEN-STICK ELBOW — RATIO FEATURES (per class pair)
%% ============================================================

figure('Position',[50 50 nCPcols*380 nCProws*300],'Color','w');
for cp = 1:nCP
    subplot(nCProws, nCPcols, cp);
    scores = sort(jsdRatioMatrix(cp,:), 'descend');
    n      = numel(scores);
    total  = sum(scores);
    if total > 0
        propScores = scores / total;
    else
        propScores = zeros(1,n);
    end
    expProp = brokenStickExpected(n);

    bar(propScores, 'FaceColor',[0.6 0.4 0.9],'EdgeColor','none'); hold on;
    plot(expProp, 'r-', 'LineWidth',1.5);

    nE = brokenStickElbow(scores);
    if nE > 0
        xline(nE + 0.5, 'k--', 'LineWidth',1.2);
        title(sprintf('%s  [k=%d]', cpLabels{cp}, nE), ...
            'FontSize',7,'FontWeight','bold','Interpreter','none');
    else
        title(cpLabels{cp},'FontSize',7,'FontWeight','bold','Interpreter','none');
    end
    xlabel('Ratio rank','FontSize',6);
    ylabel('Prop. JSD','FontSize',6);
    set(gca,'FontSize',6);
    if cp == 1
        legend({'Actual','Broken stick','Elbow'},'FontSize',5,'Location','northeast');
    end
end
sgtitle('Broken-Stick Elbow — JSD RATIO Features per Class Pair', ...
    'FontSize',12,'FontWeight','bold');
saveFig('broken_stick_elbow_ratio_all_pairs');

%% ============================================================
%% 7B. BUILD FINAL FEATURE MATRIX
%% ============================================================

fprintf('\n--- Building feature matrix ---\n');

X_indiv_train = X_train(:, selectedIndivIdx);
X_indiv_test  = X_test(:,  selectedIndivIdx);
X_indiv_val   = X_val(:,   selectedIndivIdx);

X_ratio_train = zeros(size(X_train,1), nSelRatio);
X_ratio_test  = zeros(size(X_test,1),  nSelRatio);
X_ratio_val   = zeros(size(X_val,1),   nSelRatio);
ratioNamesSelected = cell(nSelRatio,1);

for k = 1:nSelRatio
    fp = selectedRatioIdx(k);
    fi = featPairs(fp,1);
    fj = featPairs(fp,2);
    X_ratio_train(:,k) = X_train(:,fi) ./ (X_train(:,fj) + eps);
    X_ratio_test(:,k)  = X_test(:,fi)  ./ (X_test(:,fj)  + eps);
    X_ratio_val(:,k)   = X_val(:,fi)   ./ (X_val(:,fj)   + eps);
    ratioNamesSelected{k} = sprintf('%s_div_%s', featureNames{fi}, featureNames{fj});
end

X_train_exp     = [X_indiv_train, X_ratio_train];
X_test_exp      = [X_indiv_test,  X_ratio_test];
X_val_exp       = [X_indiv_val,   X_ratio_val];
featureNamesExp = [selectedIndivNames, ratioNamesSelected(:)'];

fprintf('Total features: %d  (individual: %d | ratios: %d)\n', ...
    size(X_train_exp,2), nSelIndiv, nSelRatio);

if size(X_train_exp,2) == 0
    error('No features survived the elbow filter for ANY class pair. Cannot proceed.');
end

R_exp = corr(X_train_exp, 'Rows', 'pairwise');
figure('Position', [100 100 800 700]);
imagesc(R_exp); colorbar; axis square;
xticks(1:numel(featureNamesExp));
yticks(1:numel(featureNamesExp));
xticklabels(featureNamesExp); yticklabels(featureNamesExp);
xtickangle(90);
title('Correlation Matrix — JSD-Elbow Selected Features');
saveFig('correlation_matrix_expanded');

%% ============================================================
%% 8. MIN-MAX SCALING (fit on train, transform all)
%% ============================================================

min_tr   = min(X_train_exp,[],1);
range_tr = max(X_train_exp,[],1) - min_tr;
range_tr(range_tr==0) = 1;
X_train_sc = (X_train_exp - min_tr) ./ range_tr;
X_test_sc  = (X_test_exp  - min_tr) ./ range_tr;
X_val_sc   = (X_val_exp   - min_tr) ./ range_tr;
fprintf('Scaling done. Features: %d\n', size(X_train_sc,2));

%% ============================================================
%% 9. PCA + t-SNE (using TRAINING data)
%% ============================================================

[~,score,~,~,explained] = pca(normalize(X_train_sc));
figure;
gscatter(score(:,1),score(:,2),Y_train);
xlabel(sprintf('PC1 (%.1f%%)',explained(1)));
ylabel(sprintf('PC2 (%.1f%%)',explained(2)));
title('PCA Projection (JSD-Elbow Features)');
saveFig('pca_projection');

rng(1);
Ytsne = tsne(zscore(X_train_sc),'Perplexity',30,'Standardize',true);
figure;
gscatter(Ytsne(:,1),Ytsne(:,2),Y_train);
title('t-SNE Projection (JSD-Elbow Features)');
saveFig('tsne_all');

classes  = categories(Y_train);
nClasses = numel(classes);
figure; set(gcf,'Position',[100 100 1600 800],'Color','w');
plotIdx = 1;
for i = 1:nClasses-1
    for j = i+1:nClasses
        subplot(3,7,plotIdx); hold on
        idxPair = (Y_train==classes{i}) | (Y_train==classes{j});
        rng(1);
        Ytp = tsne(zscore(X_train_sc(idxPair,:)),'Perplexity',30);
        gscatter(Ytp(:,1),Ytp(:,2),Y_train(idxPair));
        title(sprintf('%s vs %s',classes{i},classes{j}),'FontSize',9,'FontWeight','bold');
        set(gca,'XTick',[],'YTick',[],'LineWidth',2); box on;
        plotIdx = plotIdx+1;
    end
end
sgtitle('t-SNE — All Two-Class Combinations','FontSize',14,'FontWeight','bold');
saveFig('tsne_all_pairs_3x7');

%% ============================================================
%% 9B. SEPARABILITY METRIC COMPARISON — JSD vs BD
%%     Independent block: scores every individual feature and
%%     every ratio feature using two metrics:
%%       1. JSD — Jensen-Shannon Divergence
%%       2. BD  — Bhattacharyya Distance
%%     For each metric: broken-stick elbow per class pair,
%%     union of selected features across all pairs.
%%     Final output: comparison table of top features per metric.
%% ============================================================

fprintf('\n=== SECTION 9B: Separability Metric Comparison (JSD | BD) ===\n');

%% --- Storage matrices: rows = class pairs, cols = features/ratios ---
jsd_indiv_9b = zeros(nCP, nFeat);
bd_indiv_9b  = zeros(nCP, nFeat);

jsd_ratio_9b = zeros(nCP, nFP);
bd_ratio_9b  = zeros(nCP, nFP);

fprintf('Scoring %d class pairs x %d indiv + %d ratio features...\n', nCP, nFeat, nFP);

for cp = 1:nCP
    cA   = classPairs(cp,1);
    cB   = classPairs(cp,2);
    idxA = find(Y_train == classLabels{cA});
    idxB = find(Y_train == classLabels{cB});
    XA   = X_train(idxA, :);
    XB   = X_train(idxB, :);

    %% Individual features
    for fi = 1:nFeat
        jsd_indiv_9b(cp,fi) = jsdScore(XA(:,fi),   XB(:,fi),   nBins);
        bd_indiv_9b(cp,fi)  = bdScore_9b(XA(:,fi), XB(:,fi),   nBins);
    end

    %% Ratio features
    for fp = 1:nFP
        fi2 = featPairs(fp,1);
        fj2 = featPairs(fp,2);
        rA  = XA(:,fi2) ./ (XA(:,fj2) + eps);
        rB  = XB(:,fi2) ./ (XB(:,fj2) + eps);
        jsd_ratio_9b(cp,fp) = jsdScore(rA,    rB,    nBins);
        bd_ratio_9b(cp,fp)  = bdScore_9b(rA,  rB,    nBins);
    end

    if mod(cp,5)==0 || cp==nCP
        fprintf('  Completed %d / %d class pairs\n', cp, nCP);
    end
end

%% ---------------------------------------------------------------
%% Elbow selection per metric
%% ---------------------------------------------------------------
metrics9b      = {'JSD','BD'};
nMetrics9b     = numel(metrics9b);
indivMats9b    = {jsd_indiv_9b,  bd_indiv_9b};
ratioMats9b    = {jsd_ratio_9b,  bd_ratio_9b};
metricColors9b = {[0.25 0.55 0.85], [0.85 0.35 0.25]};

selected_indiv_9b = cell(nMetrics9b,1);
selected_ratio_9b = cell(nMetrics9b,1);
nSel_indiv_9b     = zeros(nMetrics9b,1);
nSel_ratio_9b     = zeros(nMetrics9b,1);

for m9 = 1:nMetrics9b
    indivMat = indivMats9b{m9};
    ratioMat = ratioMats9b{m9};

    sel_indiv = false(nFeat,1);
    for cp = 1:nCP
        [~, sIdx] = sort(indivMat(cp,:), 'descend');
        nE = brokenStickElbow(indivMat(cp,:));
        sel_indiv(sIdx(1:nE)) = true;
    end
    selected_indiv_9b{m9} = featureNames(sel_indiv);
    nSel_indiv_9b(m9)     = sum(sel_indiv);

    sel_ratio = false(nFP,1);
    for cp = 1:nCP
        [~, sIdx] = sort(ratioMat(cp,:), 'descend');
        nE = brokenStickElbow(ratioMat(cp,:));
        sel_ratio(sIdx(1:nE)) = true;
    end
    selected_ratio_9b{m9} = ratioLabelsFull(sel_ratio);
    nSel_ratio_9b(m9)     = sum(sel_ratio);

    fprintf('  %-5s  Indiv selected: %d / %d   Ratio selected: %d / %d\n', ...
        metrics9b{m9}, nSel_indiv_9b(m9), nFeat, nSel_ratio_9b(m9), nFP);
end

%% ---------------------------------------------------------------
%% ELBOW PLOTS — Individual features, one figure per metric (2 figs)
%% ---------------------------------------------------------------
for m9 = 1:nMetrics9b
    indivMat = indivMats9b{m9};
    fc = metricColors9b{m9};
    figure('Position',[50 50 nCPcols*380 nCProws*300],'Color','w');
    for cp = 1:nCP
        subplot(nCProws, nCPcols, cp);
        scores = sort(indivMat(cp,:), 'descend');
        total  = sum(scores);
        propS  = scores / max(total, eps);
        expP   = brokenStickExpected(numel(scores));
        bar(propS, 'FaceColor', fc, 'EdgeColor','none'); hold on;
        plot(expP, 'r-', 'LineWidth',1.5);
        nE = brokenStickElbow(scores);
        if nE > 0; xline(nE+0.5,'k--','LineWidth',1.2); end
        title(sprintf('%s  [k=%d]', cpLabels{cp}, nE), ...
            'FontSize',7,'FontWeight','bold','Interpreter','none');
        xlabel('Feature rank','FontSize',6);
        ylabel(sprintf('Prop. %s',metrics9b{m9}),'FontSize',6);
        set(gca,'FontSize',6);
        if cp==1
            legend({'Actual','Broken stick','Elbow'},'FontSize',5,'Location','northeast');
        end
    end
    sgtitle(sprintf('Broken-Stick Elbow — %s  |  INDIVIDUAL Features per Class Pair', ...
        metrics9b{m9}), 'FontSize',12,'FontWeight','bold');
    saveFig(sprintf('9b_elbow_indiv_%s', lower(metrics9b{m9})));
end

%% ---------------------------------------------------------------
%% ELBOW PLOTS — Ratio features, one figure per metric (2 figs)
%% ---------------------------------------------------------------
for m9 = 1:nMetrics9b
    ratioMat = ratioMats9b{m9};
    fc = metricColors9b{m9} * 0.75 + 0.25;
    figure('Position',[50 50 nCPcols*380 nCProws*300],'Color','w');
    for cp = 1:nCP
        subplot(nCProws, nCPcols, cp);
        scores = sort(ratioMat(cp,:), 'descend');
        total  = sum(scores);
        propS  = scores / max(total, eps);
        expP   = brokenStickExpected(numel(scores));
        bar(propS, 'FaceColor', fc, 'EdgeColor','none'); hold on;
        plot(expP, 'r-', 'LineWidth',1.5);
        nE = brokenStickElbow(scores);
        if nE > 0; xline(nE+0.5,'k--','LineWidth',1.2); end
        title(sprintf('%s  [k=%d]', cpLabels{cp}, nE), ...
            'FontSize',7,'FontWeight','bold','Interpreter','none');
        xlabel('Ratio rank','FontSize',6);
        ylabel(sprintf('Prop. %s',metrics9b{m9}),'FontSize',6);
        set(gca,'FontSize',6);
        if cp==1
            legend({'Actual','Broken stick','Elbow'},'FontSize',5,'Location','northeast');
        end
    end
    sgtitle(sprintf('Broken-Stick Elbow — %s  |  RATIO Features per Class Pair', ...
        metrics9b{m9}), 'FontSize',12,'FontWeight','bold');
    saveFig(sprintf('9b_elbow_ratio_%s', lower(metrics9b{m9})));
end

%% ---------------------------------------------------------------
%% MEAN SCORE BARS — one figure (2x2 grid: metric x indiv/ratio)
%% ---------------------------------------------------------------
figure('Position',[50 50 1000 700],'Color','w');
for m9 = 1:nMetrics9b
    subplot(2, nMetrics9b, m9);
    meanScoreIndiv = mean(indivMats9b{m9}, 1);
    bar(meanScoreIndiv, 'FaceColor', metricColors9b{m9}, 'EdgeColor','none');
    xticks(1:nFeat); xticklabels(featureNames); xtickangle(45);
    ylabel(sprintf('Mean %s Score', metrics9b{m9}));
    title(sprintf('%s — Individual Features (mean over pairs)', metrics9b{m9}), ...
        'FontWeight','bold');
    set(gca,'FontSize',7);

    subplot(2, nMetrics9b, nMetrics9b + m9);
    meanScoreRatio = mean(ratioMats9b{m9}, 1);
    [~, sortR] = sort(meanScoreRatio, 'descend');
    topR = min(20, nFP);
    bar(meanScoreRatio(sortR(1:topR)), 'FaceColor', metricColors9b{m9}*0.75+0.25, 'EdgeColor','none');
    xticks(1:topR);
    xticklabels(ratioLabelsFull(sortR(1:topR)));
    xtickangle(60);
    ylabel(sprintf('Mean %s Score', metrics9b{m9}));
    title(sprintf('%s — Top 20 Ratio Features (mean over pairs)', metrics9b{m9}), ...
        'FontWeight','bold');
    set(gca,'FontSize',6);
end
sgtitle('Mean Separability Score per Feature — JSD | BD', ...
    'FontSize',13,'FontWeight','bold');
saveFig('9b_mean_scores_all_metrics');

%% ---------------------------------------------------------------
%% COMPARISON TABLE — Individual features
%% ---------------------------------------------------------------
allSelectedIndiv = union(selected_indiv_9b{1}, selected_indiv_9b{2});
allSelectedIndiv = allSelectedIndiv(:);
nRows  = numel(allSelectedIndiv);
sel_by = false(nRows, nMetrics9b);
meanS  = zeros(nRows, nMetrics9b);
rankS  = zeros(nRows, nMetrics9b);

for m9 = 1:nMetrics9b
    meanScoreAll = mean(indivMats9b{m9}, 1);
    [~, globalRank] = sort(meanScoreAll, 'descend');
    rankLookup = zeros(1,nFeat);
    rankLookup(globalRank) = 1:nFeat;
    for r = 1:nRows
        fname = allSelectedIndiv{r};
        fi    = find(strcmp(featureNames, fname), 1);
        if ~isempty(fi)
            sel_by(r,m9) = ismember(fname, selected_indiv_9b{m9});
            meanS(r,m9)  = meanScoreAll(fi);
            rankS(r,m9)  = rankLookup(fi);
        end
    end
end

sel_str = repmat({''}, nRows, nMetrics9b);
for m9 = 1:nMetrics9b
    for r = 1:nRows
        if sel_by(r,m9); sel_str{r,m9} = 'YES'; else; sel_str{r,m9} = '-'; end
    end
end

nMetricsSel = sum(sel_by, 2);
[~, sortOrder] = sortrows([-nMetricsSel, -meanS(:,1)]);

CompTable = table( ...
    allSelectedIndiv(sortOrder), ...
    sel_str(sortOrder,1), sel_str(sortOrder,2), ...
    round(meanS(sortOrder,1),4), round(meanS(sortOrder,2),4), ...
    rankS(sortOrder,1), rankS(sortOrder,2), ...
    'VariableNames', { ...
        'Feature', ...
        'Sel_JSD','Sel_BD', ...
        'MeanScore_JSD','MeanScore_BD', ...
        'Rank_JSD','Rank_BD'});

fprintf('\n=== INDIVIDUAL FEATURE COMPARISON TABLE (JSD | BD) ===\n');
disp(CompTable);
writetable(CompTable, fullfile(outDir,'9b_feature_comparison_table.xlsx'), 'Sheet','Individual');

%% ---------------------------------------------------------------
%% COMPARISON TABLE — Ratio features
%% ---------------------------------------------------------------
allSelectedRatio = union(selected_ratio_9b{1}, selected_ratio_9b{2});
allSelectedRatio = allSelectedRatio(:);
nRowsR  = numel(allSelectedRatio);
sel_byR = false(nRowsR, nMetrics9b);
meanSR  = zeros(nRowsR, nMetrics9b);
rankSR  = zeros(nRowsR, nMetrics9b);

for m9 = 1:nMetrics9b
    meanScoreAll = mean(ratioMats9b{m9}, 1);
    [~, globalRank] = sort(meanScoreAll, 'descend');
    rankLookup = zeros(1,nFP);
    rankLookup(globalRank) = 1:nFP;
    for r = 1:nRowsR
        rname = allSelectedRatio{r};
        fp2   = find(strcmp(ratioLabelsFull, rname), 1);
        if ~isempty(fp2)
            sel_byR(r,m9) = ismember(rname, selected_ratio_9b{m9});
            meanSR(r,m9)  = meanScoreAll(fp2);
            rankSR(r,m9)  = rankLookup(fp2);
        end
    end
end

sel_strR = repmat({''}, nRowsR, nMetrics9b);
for m9 = 1:nMetrics9b
    for r = 1:nRowsR
        if sel_byR(r,m9); sel_strR{r,m9} = 'YES'; else; sel_strR{r,m9} = '-'; end
    end
end

nMetricsSelR = sum(sel_byR, 2);
[~, sortOrderR] = sortrows([-nMetricsSelR, -meanSR(:,1)]);

CompTableRatio = table( ...
    allSelectedRatio(sortOrderR), ...
    sel_strR(sortOrderR,1), sel_strR(sortOrderR,2), ...
    round(meanSR(sortOrderR,1),4), round(meanSR(sortOrderR,2),4), ...
    rankSR(sortOrderR,1), rankSR(sortOrderR,2), ...
    'VariableNames', { ...
        'RatioFeature', ...
        'Sel_JSD','Sel_BD', ...
        'MeanScore_JSD','MeanScore_BD', ...
        'Rank_JSD','Rank_BD'});

writetable(CompTableRatio, fullfile(outDir,'9b_feature_comparison_table.xlsx'), 'Sheet','Ratios');
fprintf('Comparison tables saved to 9b_feature_comparison_table.xlsx\n');

%% ---------------------------------------------------------------
%% SUMMARY BAR — number of features selected by each metric
%% ---------------------------------------------------------------
figure('Position',[100 100 600 400],'Color','w');
subplot(1,2,1);
b1 = bar(nSel_indiv_9b);
b1.FaceColor = 'flat';
for m9=1:nMetrics9b; b1.CData(m9,:) = metricColors9b{m9}; end
xticks(1:nMetrics9b); xticklabels(metrics9b);
ylabel('# Individual Features Selected');
title('Individual Features Selected','FontWeight','bold');
ylim([0 nFeat+1]);
for m9=1:nMetrics9b
    text(m9, nSel_indiv_9b(m9)+0.1, num2str(nSel_indiv_9b(m9)), ...
        'HorizontalAlignment','center','FontWeight','bold');
end

subplot(1,2,2);
b2 = bar(nSel_ratio_9b);
b2.FaceColor = 'flat';
for m9=1:nMetrics9b; b2.CData(m9,:) = metricColors9b{m9}; end
xticks(1:nMetrics9b); xticklabels(metrics9b);
ylabel('# Ratio Features Selected');
title('Ratio Features Selected','FontWeight','bold');
for m9=1:nMetrics9b
    text(m9, nSel_ratio_9b(m9)+0.5, num2str(nSel_ratio_9b(m9)), ...
        'HorizontalAlignment','center','FontWeight','bold');
end
sgtitle('Broken-Stick Elbow Feature Count: JSD vs BD', ...
    'FontSize',12,'FontWeight','bold');
saveFig('9b_feature_count_comparison');

fprintf('\n=== Section 9B complete ===\n');

%% ============================================================
%% 10. CLASS WEIGHTS (on training set)
%% ============================================================

trainCats    = categories(Y_train);
trainCounts  = countcats(Y_train);
nTrain       = numel(Y_train);
nCls         = numel(trainCats);
classWeights = nTrain ./ (nCls .* trainCounts);
sampleWeights= zeros(nTrain,1);
for c = 1:nCls
    sampleWeights(Y_train == trainCats{c}) = classWeights(c);
end

%% ============================================================
%% 11. INSTALL PYTHON PACKAGES
%% ============================================================

fprintf('\n--- Installing / verifying Python packages ---\n');
try
    pyExec = char(pyenv().Executable);
    system(sprintf('"%s" -m pip install xgboost lightgbm catboost scikit-learn joblib --quiet', pyExec));
catch
    warning('Could not install Python packages automatically.');
end

%% ============================================================
%% 12. IMPORT PYTHON MODULES
%% ============================================================

np     = py.importlib.import_module('numpy');
xgb    = py.importlib.import_module('xgboost');
lgb    = py.importlib.import_module('lightgbm');
cb     = py.importlib.import_module('catboost');
sk_ens = py.importlib.import_module('sklearn.ensemble');

Y_tr_int = double(grp2idx(Y_train));
catOrder = categories(Y_train);

Xtr_np   = np.array(double(X_train_sc));
ytr_np   = np.array(Y_tr_int - 1);
Xtest_np = np.array(double(X_test_sc));
ytest_np = np.array(double(grp2idx(Y_test)) - 1);
sw_np    = np.array(sampleWeights);

p2c = @(p) categorical(catOrder(double(p)+1));

%% ============================================================
%% 13. TRAIN ALL PYTHON MODELS (WEIGHTED)
%% ============================================================

fprintf('\n--- Training all Python models (weighted) ---\n');

%% --- Random Forest (sklearn, W) ---
% NOTE: class_weight is NOT set — sample_weight already encodes
% inverse-frequency weighting.  Passing both simultaneously
% double-counts the correction and degrades minority-class boundaries
% (classes 1 & 2 drop to ~35-45%).  max_features=0.3 (30% of
% features per split) gives each tree more signal, matching MATLAB
% TreeBagger behaviour more closely.  min_samples_leaf=2 stabilises
% leaf probability estimates for minority classes.
fprintf('  Random Forest (sklearn, W)...\n');
clf_rf_w = sk_ens.RandomForestClassifier(pyargs(...
    'n_estimators',     int32(1000), ...
    'max_features',     0.3, ...
    'min_samples_leaf', int32(2), ...
    'random_state',     int32(42), ...
    'n_jobs',           int32(-1)));
clf_rf_w.fit(Xtr_np, ytr_np, pyargs('sample_weight', sw_np));

%% --- XGBoost (W) ---
fprintf('  XGBoost (W)...\n');
xgbArgs = pyargs('n_estimators', int32(300), 'max_depth', int32(6), ...
    'learning_rate', 0.1, 'eval_metric', 'mlogloss', 'random_state', int32(42));
clf_xgb_w = xgb.XGBClassifier(xgbArgs);
clf_xgb_w.fit(Xtr_np, ytr_np, pyargs(...
    'sample_weight', sw_np, ...
    'eval_set', py.list({py.tuple({Xtr_np, ytr_np}), ...
                         py.tuple({Xtest_np, ytest_np})}), ...
    'verbose', false));

%% --- LightGBM (W) ---
fprintf('  LightGBM (W)...\n');
lgbArgs = pyargs('n_estimators', int32(300), 'num_leaves', int32(63), ...
    'learning_rate', 0.1, 'random_state', int32(42), 'verbose', int32(-1));
clf_lgb_w = lgb.LGBMClassifier(lgbArgs);
clf_lgb_w.fit(Xtr_np, ytr_np, pyargs(...
    'sample_weight', sw_np, ...
    'eval_set', py.list({py.tuple({Xtr_np, ytr_np}), ...
                         py.tuple({Xtest_np, ytest_np})}), ...
    'eval_metric', 'multi_logloss', ...
    'callbacks', py.list({lgb.log_evaluation(py.int(0))})));

%% --- CatBoost (W) ---
fprintf('  CatBoost (W)...\n');
catArgs = pyargs('iterations', int32(300), 'depth', int32(6), ...
    'learning_rate', 0.1, 'random_seed', int32(42), 'verbose', int32(0));
clf_cat_w = cb.CatBoostClassifier(catArgs);
clf_cat_w.fit(Xtr_np, ytr_np, pyargs(...
    'sample_weight', sw_np, ...
    'eval_set', py.tuple({Xtest_np, ytest_np}), ...
    'verbose', int32(0)));

%% --- Gradient Boosting (sklearn, W) ---
fprintf('  Gradient Boosting (W)...\n');
clf_gb_w = sk_ens.GradientBoostingClassifier(pyargs(...
    'n_estimators', int32(200), 'max_depth', int32(5), ...
    'learning_rate', 0.1, 'random_state', int32(42)));
clf_gb_w.fit(Xtr_np, ytr_np, pyargs('sample_weight', sw_np));

fprintf('  All models trained.\n');

%% ============================================================
%% 14. LOSS CURVES — ALL MODELS (train vs test)
%% ============================================================

fprintf('\n--- Plotting loss curves ---\n');

fig_loss = figure('Position',[10 10 1600 700],'Color','w');
nLoss_cols = 3; nLoss_rows = 2;

%% ---------- Random Forest: warm_start Train vs Test log-loss ----------
subplot(nLoss_rows, nLoss_cols, 1);
% Random Forest has no sequential loss like boosting models.
% We use warm_start to incrementally add trees, recording train
% log-loss and test log-loss at each step — giving a proper
% Train vs Test curve comparable to the other models.
rfScriptPath  = fullfile(outDir, 'rf_warmstart_curve.py');
rfTrLossPath  = fullfile(outDir, 'rf_train_loss.npy');
rfTeLossPath  = fullfile(outDir, 'rf_test_loss.npy');
rfStepsPath   = fullfile(outDir, 'rf_steps.npy');

jl = py.importlib.import_module('joblib');
np.save(fullfile(outDir,'Xtr_rf.npy'),  Xtr_np);
np.save(fullfile(outDir,'ytr_rf.npy'),  ytr_np);
np.save(fullfile(outDir,'Xte_rf.npy'),  Xtest_np);
np.save(fullfile(outDir,'yte_rf.npy'),  ytest_np);
np.save(fullfile(outDir,'sw_rf.npy'),   sw_np);

rfLines = { ...
    'import sys, numpy as np', ...
    'from sklearn.ensemble import RandomForestClassifier', ...
    'from sklearn.metrics import log_loss', ...
    'Xtr = np.load(sys.argv[1])', ...
    'ytr = np.load(sys.argv[2]).astype(int)', ...
    'Xte = np.load(sys.argv[3])', ...
    'yte = np.load(sys.argv[4]).astype(int)', ...
    'sw  = np.load(sys.argv[5])', ...
    'n_total = 500', ...
    'step_size = 10', ...
    'steps = list(range(step_size, n_total + 1, step_size))', ...
    'tr_loss = []', ...
    'te_loss = []', ...
    'classes = np.unique(ytr)', ...
    'rf = RandomForestClassifier(n_estimators=step_size, max_features=0.3,', ...
    '    min_samples_leaf=2, warm_start=True,', ...
    '    random_state=42, n_jobs=-1)', ...
    'for n in steps:', ...
    '    rf.n_estimators = n', ...
    '    rf.fit(Xtr, ytr, sample_weight=sw)', ...
    '    tr_prob = rf.predict_proba(Xtr)', ...
    '    te_prob = rf.predict_proba(Xte)', ...
    '    tr_loss.append(log_loss(ytr, tr_prob, labels=classes))', ...
    '    te_loss.append(log_loss(yte, te_prob, labels=classes))', ...
    'np.save(sys.argv[6], np.array(tr_loss))', ...
    'np.save(sys.argv[7], np.array(te_loss))', ...
    'np.save(sys.argv[8], np.array(steps))', ...
};
fid = fopen(rfScriptPath,'w');
for ln = 1:numel(rfLines); fprintf(fid,'%s\n',rfLines{ln}); end
fclose(fid);

pyExec = char(pyenv().Executable);
cmd = sprintf('"%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s"', ...
    pyExec, rfScriptPath, ...
    fullfile(outDir,'Xtr_rf.npy'), fullfile(outDir,'ytr_rf.npy'), ...
    fullfile(outDir,'Xte_rf.npy'), fullfile(outDir,'yte_rf.npy'), ...
    fullfile(outDir,'sw_rf.npy'),  ...
    rfTrLossPath, rfTeLossPath, rfStepsPath);
[status, cmdout] = system(cmd);
if status == 0
    rfTrLoss = double(np.load(rfTrLossPath));
    rfTeLoss = double(np.load(rfTeLossPath));
    rfSteps  = double(np.load(rfStepsPath));
    plot(rfSteps, rfTrLoss, 'b-',  'LineWidth',1.5); hold on;
    plot(rfSteps, rfTeLoss, 'r--', 'LineWidth',1.5);
    legend('Train','Test','Location','best');
    xlabel('Number of Trees'); ylabel('Log Loss');
else
    warning('RF warm_start curve script failed:\n%s', cmdout);
    text(0.5,0.5,'Loss curve unavailable','Units','normalized', ...
        'HorizontalAlignment','center','FontSize',9);
end
title('Random Forest (W)'); grid on;

%% ---------- XGBoost ----------
subplot(nLoss_rows, nLoss_cols, 2);
try
    evalRes  = clf_xgb_w.evals_result_;
    trainXGB = double(py.array.array('d', evalRes{'validation_0'}{'mlogloss'}));
    testXGB  = double(py.array.array('d', evalRes{'validation_1'}{'mlogloss'}));
    plot(trainXGB, 'b-',  'LineWidth',1.5); hold on;
    plot(testXGB,  'r--', 'LineWidth',1.5);
    legend('Train','Test','Location','best');
catch
    try
        evalRes  = clf_xgb_w.evals_result();
        trainXGB = double(py.array.array('d', evalRes{'validation_0'}{'mlogloss'}));
        testXGB  = double(py.array.array('d', evalRes{'validation_1'}{'mlogloss'}));
        plot(trainXGB, 'b-',  'LineWidth',1.5); hold on;
        plot(testXGB,  'r--', 'LineWidth',1.5);
        legend('Train','Test','Location','best');
    catch ME2
        fprintf('  XGBoost: using subprocess fallback (%s)\n', ME2.message);
        [trainXGB, testXGB] = extractXGBLoss_subprocess(clf_xgb_w, outDir, pyenv().Executable);
        if ~isempty(trainXGB)
            plot(trainXGB, 'b-',  'LineWidth',1.5); hold on;
            plot(testXGB,  'r--', 'LineWidth',1.5);
            legend('Train','Test','Location','best');
        else
            text(0.5,0.5,'Loss unavailable','Units','normalized', ...
                'HorizontalAlignment','center','FontSize',9);
        end
    end
end
xlabel('Boosting Round'); ylabel('Multi-Log Loss');
title('XGBoost (W)'); grid on;

%% ---------- LightGBM ----------
subplot(nLoss_rows, nLoss_cols, 3);
try
    evalRes_lgb = clf_lgb_w.evals_result_;
    trainLGB = double(py.array.array('d', evalRes_lgb{'training'}{'multi_logloss'}));
    testLGB  = double(py.array.array('d', evalRes_lgb{'valid_1'}{'multi_logloss'}));
    plot(trainLGB, 'b-',  'LineWidth',1.5); hold on;
    plot(testLGB,  'r--', 'LineWidth',1.5);
    legend('Train','Test','Location','best');
catch ME
    warning('LightGBM loss curve: %s', ME.message);
    text(0.5,0.5,'Loss unavailable','Units','normalized', ...
        'HorizontalAlignment','center','FontSize',9);
end
xlabel('Boosting Round'); ylabel('Multi-Log Loss');
title('LightGBM (W)'); grid on;

%% ---------- CatBoost ----------
subplot(nLoss_rows, nLoss_cols, 4);
try
    evalsResult_cb = clf_cat_w.get_evals_result();
    trainCB = double(py.array.array('d', evalsResult_cb{'learn'}{'MultiClass'}));
    testCB  = double(py.array.array('d', evalsResult_cb{'validation'}{'MultiClass'}));
    plot(trainCB, 'b-',  'LineWidth',1.5); hold on;
    plot(testCB,  'r--', 'LineWidth',1.5);
    legend('Train','Test','Location','best');
catch ME
    warning('CatBoost loss curve: %s', ME.message);
    text(0.5,0.5,'Loss unavailable','Units','normalized', ...
        'HorizontalAlignment','center','FontSize',9);
end
xlabel('Iteration'); ylabel('Multi-Class Loss');
title('CatBoost (W)'); grid on;

%% ---------- Gradient Boosting: subprocess ----------
subplot(nLoss_rows, nLoss_cols, 5);

gbScriptPath  = fullfile(outDir, 'gb_staged_loss.py');
gbModelPath   = fullfile(outDir, 'gb_model.pkl');
XtrNpyPath    = fullfile(outDir, 'Xtr_gb.npy');
XteNpyPath    = fullfile(outDir, 'Xte_gb.npy');
ytrNpyPath    = fullfile(outDir, 'ytr_gb.npy');
yteNpyPath    = fullfile(outDir, 'yte_gb.npy');
trLossNpyPath = fullfile(outDir, 'gb_train_loss.npy');
teLossNpyPath = fullfile(outDir, 'gb_test_loss.npy');

jl = py.importlib.import_module('joblib');
jl.dump(clf_gb_w, gbModelPath);
np.save(XtrNpyPath,  Xtr_np);
np.save(XteNpyPath,  Xtest_np);
np.save(ytrNpyPath,  ytr_np);
np.save(yteNpyPath,  ytest_np);

pyLines = { ...
    'import numpy as np, sys', ...
    'from joblib import load', ...
    'from sklearn.metrics import log_loss', ...
    'clf  = load(sys.argv[1])', ...
    'Xtr  = np.load(sys.argv[2])', ...
    'Xte  = np.load(sys.argv[3])', ...
    'ytr  = np.load(sys.argv[4]).astype(int)', ...
    'yte  = np.load(sys.argv[5]).astype(int)', ...
    'tr_ll = [log_loss(ytr, p) for p in clf.staged_predict_proba(Xtr)]', ...
    'te_ll = [log_loss(yte, p) for p in clf.staged_predict_proba(Xte)]', ...
    'np.save(sys.argv[6], np.array(tr_ll))', ...
    'np.save(sys.argv[7], np.array(te_ll))', ...
};
fid = fopen(gbScriptPath, 'w');
for ln = 1:numel(pyLines); fprintf(fid, '%s\n', pyLines{ln}); end
fclose(fid);

pyExec = char(pyenv().Executable);
cmd = sprintf('"%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s"', ...
    pyExec, gbScriptPath, ...
    gbModelPath, XtrNpyPath, XteNpyPath, ...
    ytrNpyPath,  yteNpyPath, ...
    trLossNpyPath, teLossNpyPath);
[status, cmdout] = system(cmd);
if status ~= 0
    warning('GB staged loss script failed:\n%s', cmdout);
    text(0.5, 0.5, 'Loss unavailable', 'Units','normalized', ...
        'HorizontalAlignment','center','FontSize',9);
else
    trainLL = double(np.load(trLossNpyPath));
    testLL  = double(np.load(teLossNpyPath));
    plot(trainLL, 'b-',  'LineWidth',1.5); hold on;
    plot(testLL,  'r--', 'LineWidth',1.5);
    legend('Train','Test','Location','best');
end
xlabel('Boosting Stage'); ylabel('Log Loss');
title('Grad. Boosting (W)'); grid on;

%% Hide unused subplot
subplot(nLoss_rows, nLoss_cols, 6); axis off;

annotation(fig_loss, 'textbox', [0 0.95 1 0.05], ...
    'String','Loss Curves — Training vs Test  |  All Python Weighted Models', ...
    'FontSize',13,'FontWeight','bold','EdgeColor','none', ...
    'HorizontalAlignment','center');
saveFig('loss_curves_all_models');

%% ============================================================
%% 15. FEATURE IMPORTANCE — TOP 20 (ALL MODELS)
%% ============================================================

fprintf('\n--- Feature Importance ---\n');
topN = 20;

pyMdls   = {clf_rf_w,   clf_xgb_w,  clf_lgb_w,  clf_cat_w,  clf_gb_w};
pyNames  = {'RF_W','XGBoost_W','LightGBM_W','CatBoost_W','GradBoost_W'};
pyLabels = {'Random Forest (W)','XGBoost (W)','LightGBM (W)','CatBoost (W)','Grad. Boosting (W)'};
nModels  = numel(pyMdls);

impVecs  = cell(nModels,1);
impValid = false(nModels,1);

for p = 1:nModels
    try
        if contains(pyNames{p},'CatBoost')
            impVecs{p} = double(py.array.array('d', pyMdls{p}.get_feature_importance()));
        else
            impVecs{p} = double(py.array.array('d', pyMdls{p}.feature_importances_));
        end
        impValid(p) = true;
    catch ME
        fprintf('  Importance skip for %s: %s\n', pyLabels{p}, ME.message);
    end
end

nImp_cols = 3; nImp_rows = 2;
fig_imp   = figure('Position',[10 10 nImp_cols*420 nImp_rows*400],'Color','w');
spIdx = 0;
for p = 1:nModels
    if ~impValid(p); continue; end
    spIdx = spIdx + 1;
    ax_i  = subplot(nImp_rows, nImp_cols, spIdx, 'Parent', fig_imp);
    imp   = impVecs{p};
    [imp_s, idx_s] = sort(imp(:)', 'descend');
    n = min(topN, numel(imp_s));
    barh(ax_i, imp_s(1:n));
    yticks(ax_i,1:n);
    yticklabels(ax_i, featureNamesExp(idx_s(1:n)));
    set(ax_i,'YDir','reverse','FontSize',7,'FontWeight','bold');
    xlabel(ax_i,'Importance','FontSize',7);
    title(ax_i, sprintf('Top %d — %s', topN, pyLabels{p}), ...
        'FontSize',8,'FontWeight','bold','Interpreter','none');
    box(ax_i,'on');
end
for sp = spIdx+1 : nImp_rows*nImp_cols
    subplot(nImp_rows,nImp_cols,sp); axis off;
end
sgtitle(fig_imp, sprintf('Top %d Feature Importance — All Python Weighted Models', topN), ...
    'FontSize',13,'FontWeight','bold');
saveFig('feat_imp_all_models');

%% ============================================================
%% 16. 5-FOLD CROSS-VALIDATION (on TRAINING set only)
%% ============================================================

fprintf('\n--- 5-Fold CV (on training set, Python models only) ---\n');
kFold = 5;
cvp   = cvpartition(Y_train,'KFold',kFold,'Stratify',true);
cv_acc     = zeros(kFold, nModels);
cv_macroF1 = zeros(kFold, nModels);
cv_wtF1    = zeros(kFold, nModels);

% Model labels aligned with pyMdls order
allImpLabels = pyLabels;

for fold = 1:kFold
    fprintf('  Fold %d/%d\n', fold, kFold);
    idxTr = training(cvp,fold); idxTe = test(cvp,fold);
    Xf_tr = X_train_sc(idxTr,:); Yf_tr = Y_train(idxTr);
    Xf_te = X_train_sc(idxTe,:); Yf_te = Y_train(idxTe);

    fCats  = categories(Yf_tr); fCounts = countcats(Yf_tr);
    nFTr   = numel(Yf_tr);      nFCls   = numel(fCats);
    fCW    = nFTr ./ (nFCls .* fCounts);
    sw_f   = zeros(nFTr,1);
    for c = 1:nFCls; sw_f(Yf_tr==fCats{c}) = fCW(c); end

    Xftr_np = np.array(double(Xf_tr));
    Xfte_np = np.array(double(Xf_te));
    yftr_np = np.array(double(grp2idx(Yf_tr))-1);
    sw_fnp  = np.array(sw_f);

    fx_rf_w = sk_ens.RandomForestClassifier(pyargs(...
        'n_estimators',     int32(300), ...
        'max_features',     0.3, ...
        'min_samples_leaf', int32(2), ...
        'random_state',     int32(42), ...
        'n_jobs',           int32(-1)));
    fx_rf_w.fit(Xftr_np, yftr_np, pyargs('sample_weight', sw_fnp));

    fx_xgb_w = xgb.XGBClassifier(pyargs(...
        'n_estimators',int32(150),'max_depth',int32(6),...
        'learning_rate',0.1,'eval_metric','mlogloss','random_state',int32(42)));
    fx_xgb_w.fit(Xftr_np, yftr_np, pyargs('sample_weight',sw_fnp));

    fx_lgb_w = lgb.LGBMClassifier(pyargs(...
        'n_estimators',int32(150),'num_leaves',int32(63),...
        'learning_rate',0.1,'random_state',int32(42),'verbose',int32(-1)));
    fx_lgb_w.fit(Xftr_np, yftr_np, pyargs('sample_weight',sw_fnp));

    fx_cat_w = cb.CatBoostClassifier(pyargs(...
        'iterations',int32(150),'depth',int32(6),...
        'learning_rate',0.1,'random_seed',int32(42),'verbose',int32(0)));
    fx_cat_w.fit(Xftr_np, yftr_np, pyargs('sample_weight',sw_fnp));

    fx_gb_w = sk_ens.GradientBoostingClassifier(pyargs(...
        'n_estimators',int32(100),'max_depth',int32(5),...
        'learning_rate',0.1,'random_state',int32(42)));
    fx_gb_w.fit(Xftr_np, yftr_np, pyargs('sample_weight',sw_fnp));

    fold_preds = { ...
        p2c(fx_rf_w.predict(Xfte_np)), ...
        p2c(fx_xgb_w.predict(Xfte_np)), ...
        p2c(fx_lgb_w.predict(Xfte_np)), ...
        p2c(fx_cat_w.predict(Xfte_np)), ...
        p2c(fx_gb_w.predict(Xfte_np))  };

    for m = 1:nModels
        Yp = fold_preds{m}; Yt = Yf_te;
        cv_acc(fold,m) = mean(Yp==Yt)*100;
        fCatsV = categories(Yt); nFC = numel(fCatsV);
        f1_c = zeros(nFC,1); sup_c = zeros(nFC,1);
        for c = 1:nFC
            tp   = sum(Yp==fCatsV{c} & Yt==fCatsV{c});
            fp_n = sum(Yp==fCatsV{c} & Yt~=fCatsV{c});
            fn   = sum(Yp~=fCatsV{c} & Yt==fCatsV{c});
            pr   = tp/max(tp+fp_n,1); re = tp/max(tp+fn,1);
            f1_c(c)  = 2*pr*re/max(pr+re,eps);
            sup_c(c) = sum(Yt==fCatsV{c});
        end
        cv_macroF1(fold,m) = mean(f1_c)*100;
        cv_wtF1(fold,m)    = sum(f1_c.*sup_c)/sum(sup_c)*100;
    end
end

cv_acc_mean     = mean(cv_acc,1)';     cv_acc_sd     = std(cv_acc,0,1)';
cv_macroF1_mean = mean(cv_macroF1,1)'; cv_macroF1_sd = std(cv_macroF1,0,1)';
cv_wtF1_mean    = mean(cv_wtF1,1)';    cv_wtF1_sd    = std(cv_wtF1,0,1)';

fprintf('\n=== 5-FOLD CV RESULTS ===\n');
for m = 1:nModels
    fprintf('%-25s  Acc=%.1f+/-%.1f  MacroF1=%.1f+/-%.1f  WtF1=%.1f+/-%.1f\n', ...
        allImpLabels{m}, cv_acc_mean(m), cv_acc_sd(m), ...
        cv_macroF1_mean(m), cv_macroF1_sd(m), cv_wtF1_mean(m), cv_wtF1_sd(m));
end

%% ============================================================
%% 17. VALIDATION PREDICTIONS (using the 10% validation set)
%% ============================================================

fprintf('\n--- Validation Predictions ---\n');
Xval_np = np.array(double(X_val_sc));

predictions = { ...
    p2c(clf_rf_w.predict(Xval_np)), ...
    p2c(clf_xgb_w.predict(Xval_np)), ...
    p2c(clf_lgb_w.predict(Xval_np)), ...
    p2c(clf_cat_w.predict(Xval_np)), ...
    p2c(clf_gb_w.predict(Xval_np))  };

modelLabelsShort = {'RandomForest_W','XGBoost_W','LightGBM_W','CatBoost_W','GradBoost_W'};
modelLabelsFull  = {'Random Forest (W)','XGBoost (W)','LightGBM (W)','CatBoost (W)','Grad. Boosting (W)'};

%% ============================================================
%% 18. MODEL ACCURACY (VALIDATION SET)
%% ============================================================

accAll = zeros(nModels,1);
for m = 1:nModels
    accAll(m) = mean(predictions{m} == Y_val)*100;
end

fprintf('\n=== VALIDATION ACCURACY ===\n');
for m = 1:nModels
    fprintf('%-25s : %.2f%%\n', modelLabelsFull{m}, accAll(m));
end

%% ============================================================
%% 19. CONFUSION MATRICES — ONE FIGURE (3 cols x 2 rows)
%% ============================================================

allCats  = categories(Y_val);
nClasses = numel(allCats);
nCM_cols = 3; nCM_rows = 2;

fig_cm = figure('Position',[10 10 1800 nCM_rows*450],'Color','w');
padL = 0.01; padR = 0.01; padT = 0.07; padB = 0.02;
pw   = (1 - padL - padR) / nCM_cols;
ph   = (1 - padT - padB) / nCM_rows;

for m = 1:nModels
    col  = mod(m-1, nCM_cols);
    row  = floor((m-1) / nCM_cols);
    xpos = padL + col * pw;
    ypos = 1 - padT - (row+1)*ph;
    hp   = uipanel(fig_cm,'Position',[xpos ypos pw ph], ...
               'BorderType','none','BackgroundColor','w');
    cm   = confusionchart(hp, Y_val, predictions{m});
    cm.Title         = sprintf('%s  (%.1f%%)', modelLabelsFull{m}, accAll(m));
    cm.RowSummary    = 'row-normalized';
    cm.ColumnSummary = 'column-normalized';
    cm.FontSize      = 7;
end
% Hide unused panel slot (5 models in 3x2 = 6 slots, one empty)
col  = mod(nModels, nCM_cols);
row  = floor(nModels / nCM_cols);
xpos = padL + col * pw;
ypos = 1 - padT - (row+1)*ph;
uipanel(fig_cm,'Position',[xpos ypos pw ph],'BorderType','none','BackgroundColor','w');

annotation(fig_cm,'textbox',[0 1-padT 1 padT-0.005], ...
    'String','Confusion Matrices — All Python Weighted Models (Validation, Row-Norm)', ...
    'FontSize',13,'FontWeight','bold','EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle');
saveFig('confusion_matrices_all_models');

%% ============================================================
%% 20. CLASS-WISE STACKED BARS — ONE FIGURE (3x2)
%% ============================================================

correct_all     = zeros(nClasses, nModels);
incorrect_all   = zeros(nClasses, nModels);
pct_correct_all = zeros(nClasses, nModels);

for m = 1:nModels
    for i = 1:nClasses
        cls = allCats{i};
        idx = (Y_val == cls);
        correct_all(i,m)   = sum(predictions{m}(idx) == cls);
        incorrect_all(i,m) = sum(predictions{m}(idx) ~= cls);
    end
    pct_correct_all(:,m) = 100*correct_all(:,m) ./ ...
        (correct_all(:,m) + incorrect_all(:,m));
end

fig_bar = figure('Position',[10 10 1200 nCM_rows*350],'Color','w');
for m = 1:nModels
    ax_b = subplot(nCM_rows, nCM_cols, m, 'Parent', fig_bar);
    bar(ax_b, 1:nClasses, [pct_correct_all(:,m), 100-pct_correct_all(:,m)], 'stacked');
    xticks(ax_b,1:nClasses); xticklabels(ax_b,allCats);
    ylim(ax_b,[0 100]);
    title(ax_b, sprintf('%s\n%.1f%%', modelLabelsFull{m}, accAll(m)), ...
        'FontSize',7,'FontWeight','bold');
    set(ax_b,'FontSize',7,'LineWidth',1); box(ax_b,'on');
end
subplot(nCM_rows,nCM_cols,nModels+1); axis off;   % hide 6th panel
axLeg = axes('Parent',fig_bar,'Position',[0.01 0.965 0.98 0.025],'Visible','off');
hold(axLeg,'on');
fill(axLeg,NaN,NaN,[0.2 0.5 0.8],'DisplayName','Correct %');
fill(axLeg,NaN,NaN,[0.9 0.3 0.3],'DisplayName','Incorrect %');
legend(axLeg,'Location','north','Orientation','horizontal','FontSize',10,'Box','off');
sgtitle(fig_bar,'Class-wise Performance — All Python Weighted Models (Validation)', ...
    'FontSize',13,'FontWeight','bold');
saveFig('classwise_performance_all_models');

%% ============================================================
%% 21. OVERALL ACCURACY BAR
%% ============================================================

figure('Position',[100 100 900 500]);
bar(accAll);
xticks(1:nModels); xticklabels(modelLabelsShort); xtickangle(30);
ylabel('Accuracy (%)','FontWeight','bold');
title('Overall Validation Accuracy — All Python Weighted Models (JSD-Elbow Features)',...
    'FontWeight','bold');
ylim([0 115]);
for m = 1:nModels
    text(m, accAll(m)+1, sprintf('%.1f%%',accAll(m)), ...
        'HorizontalAlignment','center','FontWeight','bold','FontSize',9);
end
set(gca,'FontWeight','bold','LineWidth',1.5); box on;
saveFig('overall_accuracy_all_models');

%% ============================================================
%% 22. EXPORT TABLES
%% ============================================================

ClassTable = table(allCats,'VariableNames',{'Class'});
for m = 1:nModels
    tag = modelLabelsShort{m};
    ClassTable.([tag,'_Correct_N'])   = correct_all(:,m);
    ClassTable.([tag,'_Incorrect_N']) = incorrect_all(:,m);
    ClassTable.([tag,'_Correct_pct']) = pct_correct_all(:,m);
end

SummaryTable_ML = table(modelLabelsFull(:), accAll, ...
    cv_acc_mean, cv_acc_sd, cv_macroF1_mean, cv_macroF1_sd, ...
    cv_wtF1_mean, cv_wtF1_sd, ...
    'VariableNames',{'Model','Val_Accuracy_pct', ...
    'CV_Acc_Mean','CV_Acc_SD','CV_MacroF1_Mean','CV_MacroF1_SD', ...
    'CV_WtF1_Mean','CV_WtF1_SD'});

writetable(ClassTable,      fullfile(outDir,'classwise_performance_validation.xlsx'));
writetable(SummaryTable_ML, fullfile(outDir,'overall_performance_CV.xlsx'));

disp(SummaryTable_ML);
fprintf('\nAll outputs saved to: %s\n', outDir);

%% ============================================================
%% LOCAL FUNCTIONS
%% ============================================================

function jsd = jsdScore(a, b, nBins)
% Jensen-Shannon Divergence (base-e) between two sample vectors.
% Returns a value in [0, log(2)] (~0.693 max).
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
% STRICT broken-stick elbow (first-crossing), minimum 1 feature.
%
% Finds the FIRST rank where the observed proportion drops BELOW the
% broken-stick expected proportion and keeps everything before it.
%
% At least 1 feature is always returned, even when the top feature
% itself falls below the broken-stick line.
%
% Input:  values — any numeric vector (sorted descending internally).
% Output: nElbow — number of features to keep (integer >= 1).

    values = sort(values(:)', 'descend');
    n      = numel(values);
    if n == 0;  nElbow = 1; return; end
    total  = sum(values);
    if total <= 0; nElbow = 1; return; end

    proportions = values / total;
    expected    = brokenStickExpected(n);

    % Find the first index where observed drops BELOW expected
    below       = proportions < expected;
    first_below = find(below, 1, 'first');

    if isempty(first_below)
        % Never crosses below → keep all
        nElbow = n;
    else
        % Keep features BEFORE the first crossing, minimum 1
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

function [trainLoss, testLoss] = extractXGBLoss_subprocess(clf_xgb, outDir, pyEnvExec)
% Fallback: retrieve XGBoost evals_result via subprocess.
    trainLoss = []; testLoss = [];
    try
        pyExec     = char(pyEnvExec);
        jl         = py.importlib.import_module('joblib');
        np_mod     = py.importlib.import_module('numpy');

        mdlPath    = fullfile(outDir, 'xgb_tmp.pkl');
        trPath     = fullfile(outDir, 'xgb_tr_loss.npy');
        tePath     = fullfile(outDir, 'xgb_te_loss.npy');
        scriptPath = fullfile(outDir, 'xgb_loss.py');

        jl.dump(clf_xgb, mdlPath);

        pyLines = { ...
            'import sys, joblib, numpy as np', ...
            'clf = joblib.load(sys.argv[1])', ...
            'try:', ...
            '    er = clf.evals_result_', ...
            'except AttributeError:', ...
            '    er = clf.evals_result()', ...
            'tr = list(er[''validation_0''][''mlogloss''])', ...
            'te = list(er[''validation_1''][''mlogloss''])', ...
            'np.save(sys.argv[2], np.array(tr))', ...
            'np.save(sys.argv[3], np.array(te))', ...
        };
        fid = fopen(scriptPath,'w');
        for ln = 1:numel(pyLines); fprintf(fid,'%s\n', pyLines{ln}); end
        fclose(fid);

        cmd = sprintf('"%s" "%s" "%s" "%s" "%s"', ...
            pyExec, scriptPath, mdlPath, trPath, tePath);
        [status, ~] = system(cmd);
        if status == 0
            trainLoss = double(np_mod.load(trPath));
            testLoss  = double(np_mod.load(tePath));
        end
    catch
    end
end

function bd = bdScore_9b(a, b, nBins)
% Bhattacharyya Distance between two sample vectors (for Section 9B).
% BD = -ln( sum(sqrt(pA .* pB)) )  in [0, +inf).
    allV  = [a(:); b(:)];
    edges = linspace(min(allV)-eps, max(allV)+eps, nBins+1);
    pA    = histcounts(a, edges, 'Normalization','probability') + eps;
    pB    = histcounts(b, edges, 'Normalization','probability') + eps;
    pA    = pA / sum(pA);
    pB    = pB / sum(pB);
    BC    = sum(sqrt(pA .* pB));
    BC    = min(max(BC, eps), 1);
    bd    = -log(BC);
end

function saveBothFormats(fig, figDir, name)
    exportgraphics(fig, fullfile(figDir,[name,'.png']),'Resolution',300);
    savefig(fig, fullfile(figDir,[name,'.fig']));
end