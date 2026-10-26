function E = load_ensemble_v3(folderA)
%LOAD_ENSEMBLE_V3  The 40,000 parameter sets of the v3 ensemble with their exact labels.
%   E = load_ensemble_v3(folderA) reads
%       <folderA>/Ensemble_iter<i>_v3_1e6/Excel_Files/ML_Ready_Dataset_v3_iter<i>.csv, i = 1..4
%   and returns
%       E.p          n x 14 parameters [KTOT PTOT k1 k11 k2 k3 k33 k4 k5 k55 k6 k7 k77 k8]
%       E.names      the 14 parameter names
%       E.class      exact class (Exact_Class_v3, 1..10)
%       E.use        Use_For_ML flag (false for the 32 flagged sets)
%       E.iteration  sampling iteration (1..4)
%       E.row        row in the iteration file (1 = first data row)
%   Works in MATLAB and in GNU Octave (list-valued columns in quotes are skipped).
names = {'KTOT','PTOT','k1','k11','k2','k3','k33','k4','k5','k55','k6','k7','k77','k8'};
E = struct('p', [], 'names', {names}, 'class', [], 'use', false(0,1), 'iteration', [], 'row', []);
for it = 1:4
    f = fullfile(folderA, sprintf('Ensemble_iter%d_v3_1e6', it), 'Excel_Files', ...
                 sprintf('ML_Ready_Dataset_v3_iter%d.csv', it));
    [H, M] = read_numeric_csv(f);
    col = @(nm) M(:, find(strcmp(H, nm), 1));
    X = zeros(size(M,1), 14);
    for j = 1:14, X(:, j) = col(names{j}); end
    n = size(X, 1);
    E.p = [E.p; X];
    E.class = [E.class; col('Exact_Class_v3')];
    E.use = [E.use; col('Use_For_ML') == 1];
    E.use = logical(E.use);
    E.iteration = [E.iteration; it*ones(n,1)];
    E.row = [E.row; (1:n)'];
end
end

function [H, M] = read_numeric_csv(f)
% header names and a numeric matrix (text and quoted list fields become NaN)
fid = fopen(f, 'r'); txt = fread(fid, Inf, '*char')'; fclose(fid);
txt = regexprep(txt, '"[^"]*"', 'NaN');
L = regexp(txt, '\r?\n', 'split'); L = L(~cellfun(@isempty, L));
H = strtrim(strsplit(L{1}, ','));
nC = numel(H); M = nan(numel(L)-1, nC);
for i = 2:numel(L)
    v = str2double(strsplit(L{i}, ','));
    M(i-1, 1:min(nC, numel(v))) = v(1:min(nC, numel(v)));
end
end
