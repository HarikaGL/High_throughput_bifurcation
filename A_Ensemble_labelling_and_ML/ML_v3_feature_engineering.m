%% ML_v3_feature_engineering.m
% Feature engineering and feature selection (Supplementary Note 9, Supplementary
% Figs 9-12 and Supplementary Tables 6-7; main-text Fig 12, block 4).
%
%   Feature engineering. Each parameter set is described by its 14 sampled parameters
%   (log10) and all 91 pairwise log-ratios log10(p_i/p_j) = log10 p_i - log10 p_j
%   (feature set A, 105 features). Feature sets B-F of ML_v3_10class_parallel.m are
%   derived from A (B: without the 7 design-constrained ratios; C: the 14 parameters
%   only; D: JSD + broken-stick selection; E, F: centred log-ratios).
%
%   Feature selection (feature set D). For every pair of classes, every feature is
%   scored by the Jensen-Shannon divergence (JSD, natural log, 30-bin histograms over the
%   common range) between its two class-conditional distributions; the features are
%   ranked, and those before the broken-stick elbow (first rank at which the normalised
%   JSD falls below the broken-stick expectation) are kept. The union over all class
%   pairs is the selected set. Constant and duplicated columns are removed first. This is
%   the procedure of selectJSDcols in support/ml_v3_core.m, applied here to all sets of a
%   scope (the classifiers apply it to the training rows of each split).
%
%   Scopes: the pooled data (45 class pairs) and each sampling iteration (Iteration 1:
%   Classes 1, 2, 5, 10; 2: 3, 4, 5, 9; 3: 6, 10; 4: 7, 8, 10).
%
% Outputs (ML_v3_results_parallel/labels_Exact_Class_v3/feature_engineering/):
%   SFig_JSD_pooled.png          broken-stick profiles, all 45 class pairs (pooled)
%   SFig_JSD_within.png          broken-stick profiles, the 16 class pairs within iterations
%   SFig_JSD_summary.png         features kept per class pair; selection frequency
%   SFig_feature_pairs.png       two-parameter feature matrices for two class pairs
%   JSD_selected_features.txt    selected features per scope and per class pair
%   JSD_scores_<scope>.csv       JSD of every feature for every class pair
% Runs in MATLAB (R2020b or newer) and GNU Octave; base functions only.

clear; close all;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
addpath(fullfile(HERE, 'support'));
OUT = fullfile(HERE, 'ML_v3_results_parallel', 'labels_Exact_Class_v3', 'feature_engineering');
if ~exist(OUT, 'dir'), mkdir(OUT); end
nBins = 30;

%% data and features (feature set A)
E = load_ensemble_v3(HERE);
keep = E.use; y = E.class(keep); it = E.iteration(keep);
fn = E.names; L = log10(E.p(keep,:));
pairs = nchoosek(1:14, 2);
X = [L, L(:,pairs(:,1)) - L(:,pairs(:,2))];
names = [fn, arrayfun(@(k) [fn{pairs(k,1)} '/' fn{pairs(k,2)}], 1:size(pairs,1), 'UniformOutput', false)];
consNames = {'KTOT/PTOT','k1/k5','k11/k55','k2/k6','k3/k7','k33/k77','k4/k8'};
isCons = ismember(names, consNames);
isPar  = [true(1,14), false(1,91)];
FC = feature_colors(names);
fprintf('%d parameter sets, %d features (14 parameters + %d log-ratios)\n', numel(y), numel(names), size(pairs,1));

%% JSD + broken stick in every scope
scopes = {'pooled', 'Iteration 1', 'Iteration 2', 'Iteration 3', 'Iteration 4'};
R = struct('scope', {}, 'cols', {}, 'pairs', {}, 'score', {}, 'nElbow', {}, 'sel', {});
for s = 1:numel(scopes)
    if s == 1, rows = true(size(y)); else, rows = it == s-1; end
    Xs = X(rows,:); ys = y(rows);
    cols = find(std(Xs, 0, 1) > 0);                         % non-constant
    [~, ia] = unique(round(Xs(:,cols)*1e9)', 'rows', 'stable');  % de-duplicated
    cols = cols(sort(ia));
    cls = unique(ys)'; cp = nchoosek(cls, 2);
    sc = zeros(size(cp,1), numel(cols)); ne = zeros(size(cp,1), 1); sel = false(1, numel(cols));
    for q = 1:size(cp,1)
        A = Xs(ys == cp(q,1), cols); B = Xs(ys == cp(q,2), cols);
        for j = 1:numel(cols), sc(q,j) = jsd_score(A(:,j), B(:,j), nBins); end
        ne(q) = broken_stick_elbow(sc(q,:));
        [~, o] = sort(sc(q,:), 'descend'); sel(o(1:ne(q))) = true;
    end
    R(s).scope = scopes{s}; R(s).cols = cols; R(s).pairs = cp; R(s).score = sc;
    R(s).nElbow = ne; R(s).sel = cols(sel);
    fprintf('%-12s %3d features after removal of constant/duplicated columns, %2d class pairs, %3d selected (%d parameters, %d ratios, %d design-constrained)\n', ...
        scopes{s}, numel(cols), size(cp,1), sum(sel), sum(isPar(cols(sel))), sum(~isPar(cols(sel))), sum(isCons(cols(sel))));
end

%% text output
fid = fopen(fullfile(OUT, 'JSD_selected_features.txt'), 'w');
for s = 1:numel(R)
    fprintf(fid, '== %s: %d of %d features selected (%d parameters, %d log-ratios, %d design-constrained ratios)\n', ...
        R(s).scope, numel(R(s).sel), numel(R(s).cols), sum(isPar(R(s).sel)), sum(~isPar(R(s).sel)), sum(isCons(R(s).sel)));
    fprintf(fid, '   %s\n', strjoin(names(R(s).sel), ', '));
    for q = 1:size(R(s).pairs,1)
        [v, o] = sort(R(s).score(q,:), 'descend'); k = 1:R(s).nElbow(q);
        lst = arrayfun(@(i) sprintf('%s (%.3f)', names{R(s).cols(o(i))}, v(i)), k, 'UniformOutput', false);
        fprintf(fid, '   Class %d vs %d: %d kept: %s\n', R(s).pairs(q,1), R(s).pairs(q,2), numel(k), strjoin(lst, ', '));
    end
    f2 = fopen(fullfile(OUT, sprintf('JSD_scores_%s.csv', strrep(R(s).scope, ' ', ''))), 'w');
    fprintf(f2, 'pair,%s\n', strjoin(names(R(s).cols), ','));
    for q = 1:size(R(s).pairs,1)
        fprintf(f2, '%d vs %d%s\n', R(s).pairs(q,1), R(s).pairs(q,2), sprintf(',%.5f', R(s).score(q,:)));
    end
    fclose(f2);
end
fclose(fid);
type(fullfile(OUT, 'JSD_selected_features.txt'));

%% Supplementary Fig: broken-stick profiles (pooled, 45 pairs; within iterations, 16 pairs)
plot_profiles(R(1), names, FC, isCons, 5, 9, 'Pooled data: all 45 class pairs', ...
              fullfile(OUT, 'SFig_JSD_pooled.png'), [36 25]);
G = struct('R', {}, 'q', {});
for s = 2:5
    for q = 1:size(R(s).pairs,1), G(end+1).R = s; G(end).q = q; end %#ok<SAGROW>
end
plot_profiles_within(R, G, names, FC, isCons, fullfile(OUT, 'SFig_JSD_within.png'));

%% Supplementary Fig: summary of the selection
fig = new_figure([32 13]);
% (A) features kept per class pair, pooled
axes('Position',[0.05 0.14 0.27 0.72]);
M = nan(10); for q = 1:size(R(1).pairs,1), M(R(1).pairs(q,2), R(1).pairs(q,1)) = R(1).nElbow(q); end
hold on; mx = max(M(:));
for i = 1:10, for j = 1:10, if ~isnan(M(i,j))   % grey level: number of features kept
    g = 1 - 0.7*sqrt(M(i,j)/mx);
    patch(j + [-0.5 0.5 0.5 -0.5], i + [-0.5 -0.5 0.5 0.5], [g g g], 'EdgeColor', [0.8 0.8 0.8]);
    text(j, i, sprintf('%d', M(i,j)), 'HorizontalAlignment','center','FontSize',8, 'Color', [1 1 1]*(g < 0.5));
end, end, end
axis([0.5 10.5 0.5 10.5]); axis square; set(gca, 'YDir', 'reverse');
same = false(10); for s = 2:5, c = unique(y(it == s-1)); same(c, c) = true; end
for i = 1:10, for j = 1:10, if ~isnan(M(i,j)) && same(i,j), rectangle('Position',[j-0.5 i-0.5 1 1],'EdgeColor',[0.85 0.2 0.2],'LineWidth',1.5); end, end, end
set(gca, 'XTick', 1:10, 'YTick', 1:10, 'FontSize', 9, 'TickLength', [0 0]); xlabel('class'); ylabel('class');
title({'(A) Features kept per class pair', '(pooled; red frame: both classes in one iteration)'}, 'FontSize', 10);
% (B) selection frequency, pooled
axes('Position',[0.40 0.14 0.25 0.72]);
cnt = zeros(1, numel(names));
for q = 1:size(R(1).pairs,1)
    [~, o] = sort(R(1).score(q,:), 'descend'); c = R(1).cols(o(1:R(1).nElbow(q))); cnt(c) = cnt(c) + 1;
end
[v, o] = sort(cnt, 'descend'); o = o(v > 0); o = o(1:min(20, numel(o)));
hold on;
for i = 1:numel(o)
    barh(i, cnt(o(i)), 0.75, 'FaceColor', FC(o(i),:), 'EdgeColor', 'k', 'LineWidth', 0.3 + 1.2*isCons(o(i)));
end
set(gca, 'YDir', 'reverse', 'YTick', 1:numel(o), 'YTickLabel', strcat(names(o), repmat({''}, 1, numel(o))), 'FontSize', 8);
lab = get(gca, 'YTickLabel'); for i = 1:numel(o), if isCons(o(i)), lab{i} = [lab{i} ' *']; end, end; set(gca, 'YTickLabel', lab);
xlabel('class pairs in which the feature is kept (of 45)'); box on; ylim([0.3 numel(o)+0.7]);
title({'(B) Most often kept features (pooled)', '* design-constrained ratio'}, 'FontSize', 10);
% (C) within iterations: selected features and the JSD of the best one
axes('Position',[0.73 0.14 0.25 0.72]); hold on;
yy = 0; lbl = {};
for s = 2:5
    sel = R(s).sel; [~, ord] = sort(max(R(s).score(:, ismember(R(s).cols, sel)), [], 1), 'descend');
    sel = sel(ord); sel = sel(1:min(6, numel(sel)));
    for k = 1:numel(sel)
        yy = yy + 1; j = find(R(s).cols == sel(k));
        barh(yy, max(R(s).score(:, j)), 0.75, 'FaceColor', FC(sel(k),:), 'EdgeColor', 'k', 'LineWidth', 0.3);
        lbl{yy} = names{sel(k)}; %#ok<SAGROW>
    end
    text(0.93, yy - numel(sel)/2 + 0.5, sprintf('Iteration %d\n(%d selected)', s-1, numel(R(s).sel)), 'FontSize', 8, 'HorizontalAlignment', 'right');
    yy = yy + 1; lbl{yy} = '';
end
set(gca, 'YDir', 'reverse', 'YTick', find(~cellfun(@isempty, lbl)), 'YTickLabel', lbl(~cellfun(@isempty, lbl)), 'FontSize', 8);
xlim([0 0.95]); ylim([0.3 yy+0.2]); box on; xlabel('largest JSD over the class pairs of the iteration');
title({'(C) Selected within each iteration', '(top six by JSD)'}, 'FontSize', 10);
print(fig, fullfile(OUT, 'SFig_JSD_summary.png'), '-dpng', '-r300');

%% Supplementary Fig: two-parameter feature matrices
% the four features with the largest JSD for the pair, plus the ratio that completes the
% analytical boundary (criterion (10): k8*PTOT/k2*KTOT; Classes 7/8: Km5 against Km1)
ex = {2, 1, 5, 'KTOT/PTOT', {'(A) Iteration 1: Class 1, LP S-shaped only (green)', 'against Class 5, LP on detached branch only (brown)'}; ...
      5, 7, 8, 'k11/k55',   {'(B) Iteration 4: Class 7, LP + BP with LP > BP (orange)', 'against Class 8, LP + BP with LP < BP (violet)'}};
fig = new_figure([36 19]);
CI = class_info_v3();
for e = 1:size(ex,1)
    s = ex{e,1}; c1 = ex{e,2}; c2 = ex{e,3};
    rows = it == s-1 & ismember(y, [c1 c2]);
    q = find(R(s).pairs(:,1) == c1 & R(s).pairs(:,2) == c2);
    [~, o] = sort(R(s).score(q,:), 'descend'); f = R(s).cols(o(1:4));
    f = [f, find(strcmp(names, ex{e,4}))]; nf = numel(f);
    jj = arrayfun(@(c) find(R(s).cols == c), f);   % columns of f in the JSD table
    rng(1); r1 = find(rows & y == c1); r2 = find(rows & y == c2);
    r1 = r1(randperm(numel(r1), min(1500, numel(r1)))); r2 = r2(randperm(numel(r2), min(1500, numel(r2))));
    x0 = 0.055 + 0.5*(e-1); y0 = 0.085; w = 0.43/nf; h = 0.79/nf;
    annotation('textbox', [x0-0.04 0.885 0.47 0.11], 'String', ex{e,5}, ...
               'EdgeColor','none','FontWeight','bold','FontSize',10,'VerticalAlignment','top');
    for a = 1:nf
        for b = 1:nf
            axes('Position', [x0+(b-1)*w+0.005, y0+(nf-a)*h+0.005, w-0.012, h-0.012]); hold on; box on;
            if a == b
                ed = linspace(min(X(rows,f(a))), max(X(rows,f(a))), 31);
                h1 = histc(X(y==c1 & rows, f(a)), ed); h2 = histc(X(y==c2 & rows, f(a)), ed);
                stairs(ed, h1/sum(h1), 'Color', CI.colors(c1,:)*0.8, 'LineWidth', 1.4);
                stairs(ed, h2/sum(h2), 'Color', CI.colors(c2,:)*0.8, 'LineWidth', 1.4);
                text(0.03, 0.92, sprintf('JSD %.2f', R(s).score(q, jj(a))), 'Units','normalized','FontSize',7);
                set(gca, 'YTick', []);
            else
                plot(X(r1, f(b)), X(r1, f(a)), '.', 'Color', CI.colors(c1,:)*0.85, 'MarkerSize', 3);
                plot(X(r2, f(b)), X(r2, f(a)), '.', 'Color', CI.colors(c2,:)*0.85, 'MarkerSize', 3);
            end
            set(gca, 'FontSize', 6);
            if a == nf, xlabel(['log_{10} ' names{f(b)}], 'FontSize', 8); else, set(gca, 'XTickLabel', []); end
            if b == 1 && a ~= 1, ylabel(['log_{10} ' names{f(a)}], 'FontSize', 8); elseif b > 1, set(gca, 'YTickLabel', []); end
        end
    end
end
print(fig, fullfile(OUT, 'SFig_feature_pairs.png'), '-dpng', '-r250');
fprintf('Feature-engineering figures written to %s\n', OUT);

%% ---------------------------------------------------------------------------------
function plot_profiles(Rs, names, FC, isCons, nr, nc, ttl, fout, sizeCm)
fig = new_figure(sizeCm);
nTop = 12; dx = 0.94/nc; dy = 0.95/nr;
for q = 1:size(Rs.pairs,1)
    r = ceil(q/nc); c = q - (r-1)*nc;
    axes('Position', [0.035+(c-1)*dx, 0.955-r*dy+0.40*dy, 0.82*dx, 0.53*dy]); hold on; box on;
    draw_profile(Rs.score(q,:), Rs.cols, Rs.nElbow(q), names, FC, isCons, nTop);
    title(sprintf('Class %d vs %d  (%d kept)', Rs.pairs(q,1), Rs.pairs(q,2), Rs.nElbow(q)), 'FontSize', 7);
end
annotation('textbox', [0.02 0.975 0.96 0.02], 'String', ttl, 'EdgeColor','none', 'FontWeight','bold', 'FontSize', 10);
print(fig, fout, '-dpng', '-r250');
end

function plot_profiles_within(R, G, names, FC, isCons, fout)
fig = new_figure([36 25]);
nTop = 15; nc = 4;
for g = 1:numel(G)
    s = G(g).R; q = G(g).q; Rs = R(s);
    r = ceil(g/nc); c = g - (r-1)*nc;
    axes('Position', [0.05+(c-1)*0.24, 0.995-r*0.245+0.09, 0.20, 0.125]); hold on; box on;
    draw_profile(Rs.score(q,:), Rs.cols, Rs.nElbow(q), names, FC, isCons, nTop);
    title(sprintf('%s: Class %d vs %d (%d kept)', Rs.scope, Rs.pairs(q,1), Rs.pairs(q,2), Rs.nElbow(q)), 'FontSize', 8);
end
print(fig, fout, '-dpng', '-r250');
end

function draw_profile(sc, cols, ne, names, FC, isCons, nTop)
n = numel(sc); [v, o] = sort(sc, 'descend');
p = v / max(sum(v), eps);
ex = arrayfun(@(k) (1/n) * sum(1./(k:n)), 1:n);
k = min(nTop, n);
for i = 1:k
    j = cols(o(i));
    bar(i, p(i), 0.8, 'FaceColor', FC(j,:), 'EdgeColor', 'k', 'LineWidth', 0.2 + 1.0*isCons(j));
end
plot(1:k, ex(1:k), 'r-', 'LineWidth', 1.2);
plot([ne+0.5 ne+0.5], [0 max([p(1:k), ex(1)])*1.05], 'k--', 'LineWidth', 0.8);
lab = names(cols(o(1:k)));
for i = 1:k, if isCons(cols(o(i))), lab{i} = [lab{i} '*']; end, end
yl = max([p(1:k), ex(1)])*1.08;
set(gca, 'XTick', 1:k, 'XTickLabel', {}, 'FontSize', 5.5, 'XLim', [0.3 k+0.7], 'YLim', [0 yl]);
for i = 1:k   % rotated feature names under the bars
    text(i, -0.02*yl, lab{i}, 'Rotation', 90, 'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'middle', 'FontSize', 5.5);
end
end

function fig = new_figure(sizeCm)
% figure of a fixed printed size (independent of the screen size)
fig = figure('Color','w','Units','centimeters','Position',[1 1 sizeCm], ...
             'PaperUnits','centimeters','PaperPositionMode','manual', ...
             'PaperPosition',[0 0 sizeCm],'PaperSize',sizeCm);
end

function jsd = jsd_score(a, b, nBins)
% Jensen-Shannon divergence of two samples (histograms on the common range), as in
% support/ml_v3_core.m
allV = [a(:); b(:)];
edges = linspace(min(allV)-eps, max(allV)+eps, nBins+1);
pA = histc(a(:), edges)'; pA = pA(1:nBins); pA(end) = pA(end) + sum(a(:) == edges(end));
pB = histc(b(:), edges)'; pB = pB(1:nBins); pB(end) = pB(end) + sum(b(:) == edges(end));
pA = pA/numel(a) + eps; pB = pB/numel(b) + eps;
pA = pA/sum(pA); pB = pB/sum(pB);
Mx = 0.5*(pA + pB);
jsd = max(0, 0.5*sum(pA .* log(pA ./ Mx)) + 0.5*sum(pB .* log(pB ./ Mx)));
end

function nElbow = broken_stick_elbow(values)
% number of features before the first rank at which the normalised score falls below the
% broken-stick expectation (as in support/ml_v3_core.m)
values = sort(values(:)', 'descend'); n = numel(values);
if n == 0 || sum(values) <= 0, nElbow = 1; return; end
proportions = values / sum(values);
expected = arrayfun(@(k) (1/n) * sum(1./(k:n)), 1:n);
first_below = find(proportions < expected, 1, 'first');
if isempty(first_below), nElbow = n; else, nElbow = max(first_below - 1, 1); end
end
