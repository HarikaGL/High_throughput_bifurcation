function W = tr_window_sets(E, Q, resDir)
%TR_WINDOW_SETS  Sets in the detached-branch quadrant (x < 0, y > 0, criterion (10)) that
%   are labelled 1, 3, 9 or 10 (symmetric sets of Iterations 3-4, u* = r, excluded), i.e. whose detached branch is not counted by the label.
%   For each: the ATOT of its detached fold (exact; expected above the window 1e6) and the
%   distance to the quadrant edges, margin = min(-x, y) (decades). Cached in
%   resDir/TR_window_sets.mat; also written as TR_window_sets.csv.
%   W.idx (rows of E), W.class, W.x, W.y, W.detFold, W.margin; W.nClass(c) = number of
%   ML sets of class c; W.marginDet = margins of the Class 2, 4, 5 sets (for comparison).
f = fullfile(resDir, 'TR_window_sets.mat');
if exist(f, 'file') == 2, L = load(f); W = L.W; return; end
k = find(E.use & ismember(E.class, [1 3 9 10]) & Q.det & ~Q.sym);   % symmetric sets have no detached branch
df = nan(numel(k), 1); P = E.p(k, :);
parfor j = 1:numel(k)
    [~, S] = tr_label(P(j, :), 4000);
    if ~isempty(S), df(j) = S.detached_fold_ATOT; end
end
W = struct('idx', k, 'class', E.class(k), 'iteration', E.iteration(k), 'row', E.row(k), ...
           'x', Q.x(k), 'y', Q.y(k), 'detFold', df, 'margin', min(-Q.x(k), Q.y(k)));
W.nClass = arrayfun(@(c) sum(E.use & E.class == c), 1:10);
kd = E.use & ismember(E.class, [2 4 5]);
W.marginDet = min(-Q.x(kd), Q.y(kd));
W.nUse = sum(E.use);
save(f, 'W', '-v7');
tr_write_csv(fullfile(resDir, 'TR_window_sets.csv'), 'iteration,row,class,x,y,margin,log10_detached_fold_ATOT', ...
    [W.iteration, W.row, W.class, W.x, W.y, W.margin, log10(W.detFold)]);
end
