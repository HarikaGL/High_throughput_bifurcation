function E = tr_load_ensemble(HERE)
%TR_LOAD_ENSEMBLE  The 40,000 parameter sets of the v3 ensemble with their exact labels.
%   Reads ../Ensemble_iter<i>_v3_1e6/Excel_Files/ML_Ready_Dataset_v3_iter<i>.csv (folder A);
%   falls back to data/Ensemble_v3_parameters_and_classes.csv in this folder.
%   E.p (n x 14), E.class (Exact_Class_v3), E.use (Use_For_ML), E.iteration, E.row (row in
%   the iteration file, 1 = first data row).
names = {'KTOT','PTOT','k1','k11','k2','k3','k33','k4','k5','k55','k6','k7','k77','k8'};
PARENT = fileparts(HERE);
E = struct('p', [], 'class', [], 'use', [], 'iteration', [], 'row', [], 'names', {names});
f1 = fullfile(PARENT, 'Ensemble_iter1_v3_1e6', 'Excel_Files', 'ML_Ready_Dataset_v3_iter1.csv');
if exist(f1, 'file') && exist('readtable', 'file')
    for it = 1:4
        T = readtable(fullfile(PARENT, sprintf('Ensemble_iter%d_v3_1e6', it), 'Excel_Files', ...
                      sprintf('ML_Ready_Dataset_v3_iter%d.csv', it)));
        n = height(T); X = zeros(n, 14);
        for j = 1:14, X(:, j) = double(T.(names{j})); end
        E.p = [E.p; X]; E.class = [E.class; double(T.Exact_Class_v3)];
        E.use = [E.use; double(T.Use_For_ML) == 1]; E.iteration = [E.iteration; it*ones(n,1)];
        E.row = [E.row; (1:n)'];
    end
else
    M = dlmread(fullfile(HERE, 'data', 'Ensemble_v3_parameters_and_classes.csv'), ',', 1, 0);
    E.iteration = M(:,1); E.row = M(:,2); E.class = M(:,3); E.use = M(:,4) == 1; E.p = M(:,5:18);
end
end
