% ======================= Fig_normal_forms_v3.m =====================================
% Normal-form coefficients of every singularity of the ten-class ensemble
% (manuscript figure "Normal-form coefficients at every branch point and fold").
%   (a) quadratic branching coefficient |a| at all branch points, Iteration 3
%       (pitchforks, a = 0 to round-off) and Iteration 4 (transcritical, a ~= 0)
%   (b) Iteration 4: slope s of the crossing branch at the branch point, Eq. (13),
%       against log10(K_M5/K_M1), coloured by exact class (7: x2_LP > x2_BP,
%       8: x2_LP < x2_BP); sign(s) = sign(K_M5 - K_M1) decides the class, Eq. (14)
%   (c) fold coefficient |a_LP| at every fold: main branch (Iter. 1-2), detached
%       branch (Iter. 1-2), crossing branch (Iter. 4), flagged folds (Iter. 3)
% Reads only the CSV output of Ensemble_LHS_v3.m:
%   Ensemble_iter<i>_v3_1e6/Excel_Files/Normal_Forms_v3_iter<i>.csv
%   Ensemble_iter4_v3_1e6/Excel_Files/ML_Ready_Dataset_v3_iter4.csv
% Run from the Ensemble_v3 folder. Base MATLAB (also runs in GNU Octave).
% Output: Figures_combined_v3/Normal_Forms_v3_summary.png / .fig
% ===================================================================================
clear; close all; clc;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
cd(HERE);

LW = 3; FS = 15; floorA = 1e-20;              % |a| below floorA drawn at floorA
c3 = [0.000 0.447 0.741]; c4 = [0.850 0.325 0.098];          % Iteration 3 / 4
c7 = [0.965 0.741 0.431]; c8 = [0.843 0.557 0.878];          % Classes 7 / 8 (class_info_v3)
cM = [0.000 0.447 0.741]; cD = [0.494 0.184 0.556]; cX = [0.850 0.325 0.098]; cF = [0.4 0.4 0.4];

%% data -------------------------------------------------------------------------------
NF = cell(1,4);
for it = 1:4
    NF{it} = readCSV(fullfile(HERE, sprintf('Ensemble_iter%d_v3_1e6', it), 'Excel_Files', ...
        sprintf('Normal_Forms_v3_iter%d.csv', it)));
end
M4 = readCSV(fullfile(HERE, 'Ensemble_iter4_v3_1e6', 'Excel_Files', 'ML_Ready_Dataset_v3_iter4.csv'));

aBP3 = abs(num(NF{3}, 'a', strcmp(NF{3}.Kind, 'BP')));
aBP4 = abs(num(NF{4}, 'a', strcmp(NF{4}.Kind, 'BP')));
aLPm = []; aLPd = [];
for it = 1:2
    isLP = strcmp(NF{it}.Kind, 'LP');
    aLPm = [aLPm; abs(num(NF{it}, 'a', isLP & strcmp(NF{it}.Branch, 'main')))];     %#ok<AGROW>
    aLPd = [aLPd; abs(num(NF{it}, 'a', isLP & strcmp(NF{it}.Branch, 'detached')))]; %#ok<AGROW>
end
aLPx = abs(num(NF{4}, 'a', strcmp(NF{4}.Kind, 'LP')));
aLPf = abs(num(NF{3}, 'a', strcmp(NF{3}.Kind, 'LP')));

% slope s, Eq. (13), for the Iteration-4 sets of Classes 7 and 8
p = @(n) num(M4, n, true(size(M4.KTOT)));
k1 = p('k1'); k11 = p('k11'); k2 = p('k2'); k3 = p('k3'); k33 = p('k33'); k4 = p('k4');
k5 = p('k5'); k55 = p('k55'); k6 = p('k6'); KT = p('KTOT'); cls = p('Exact_Class_v3');
KM1 = (k11 + k2)./k1; KM3 = (k33 + k4)./k3; KM5 = (k55 + k6)./k5;      % k6 = k2 in Iteration 4
s = k4.^2.*(KM5 - KM1) ./ (KT.*(k4.^2 - k2.^2) + k4.^2.*(KM1 + KM5) + k2.*k4.*KM3);
x = log10(KM5./KM1);
i7 = cls == 7; i8 = cls == 8;
fprintf('BP Iter 3: n = %d, max |a| = %.2g | BP Iter 4: n = %d, min |a| = %.2g\n', ...
    numel(aBP3), max(aBP3), numel(aBP4), min(aBP4));
fprintf('Class 7: %d sets, s < 0 in %d | Class 8: %d sets, s > 0 in %d\n', ...
    sum(i7), sum(i7 & s < 0), sum(i8), sum(i8 & s > 0));
fprintf('Folds: main %d, detached %d, crossing %d, flagged %d (total %d)\n', numel(aLPm), ...
    numel(aLPd), numel(aLPx), numel(aLPf), numel(aLPm)+numel(aLPd)+numel(aLPx)+numel(aLPf));

%% figure -----------------------------------------------------------------------------
fig = figure('Color','w', 'Position', [50 50 1800 720]);
% (a)
ax = subplot(1,3,1); hold(ax, 'on'); box(ax, 'on');
e = -20:0.5:0;
h3 = histc(log10(max(aBP3, floorA)), e); h4 = histc(log10(max(aBP4, floorA)), e);
stairs(ax, e, h3, '-', 'Color', c3, 'LineWidth', LW);
stairs(ax, e, h4, '-', 'Color', c4, 'LineWidth', LW);
set(ax, 'YScale', 'log', 'FontSize', FS, 'LineWidth', 2);
xlabel(ax, 'log_{10} |a|  (branch points)'); ylabel(ax, 'Number of branch points');
legend(ax, {sprintf('Iteration 3 (n = %d): pitchfork', numel(aBP3)), ...
    sprintf('Iteration 4 (n = %d): transcritical', numel(aBP4))}, 'Location', 'southoutside');
title(ax, '(a)'); xlim(ax, [-20.5 0.5]);
% (b)
ax = subplot(1,3,2); hold(ax, 'on'); box(ax, 'on');
plot(ax, x(i7), s(i7), 'o', 'MarkerSize', 4, 'MarkerFaceColor', c7, 'MarkerEdgeColor', c7);
plot(ax, x(i8), s(i8), 'o', 'MarkerSize', 4, 'MarkerFaceColor', c8, 'MarkerEdgeColor', c8);
yl = get(ax, 'YLim'); xl = get(ax, 'XLim');
plot(ax, [0 0], yl, 'k--', 'LineWidth', 2); plot(ax, xl, [0 0], 'k--', 'LineWidth', 2);
set(ax, 'FontSize', FS, 'LineWidth', 2);
xlabel(ax, 'log_{10}(K_{M5}/K_{M1})'); ylabel(ax, 'slope s of the crossing branch, Eq. (13)');
legend(ax, {sprintf('Class 7 (n = %d): s < 0', sum(i7)), sprintf('Class 8 (n = %d): s > 0', sum(i8))}, ...
    'Location', 'southoutside');
title(ax, '(b)  Iteration 4');
% (c)
ax = subplot(1,3,3); hold(ax, 'on'); box(ax, 'on');
e = -20:0.5:0;
H = {aLPm, aLPd, aLPx, aLPf}; C = {cM, cD, cX, cF};
for k = 1:4
    stairs(ax, e, histc(log10(max(H{k}, floorA)), e), '-', 'Color', C{k}, 'LineWidth', LW);
end
set(ax, 'YScale', 'log', 'FontSize', FS, 'LineWidth', 2);
xlabel(ax, 'log_{10} |a_{LP}|  (folds)'); ylabel(ax, 'Number of folds');
legend(ax, {sprintf('main branch, Iter. 1-2 (n = %d)', numel(aLPm)), ...
    sprintf('detached branch, Iter. 1-2 (n = %d)', numel(aLPd)), ...
    sprintf('crossing branch, Iter. 4 (n = %d)', numel(aLPx)), ...
    sprintf('flagged, Iter. 3 (n = %d)', numel(aLPf))}, 'Location', 'southoutside');
title(ax, '(c)'); xlim(ax, [-20.5 0.5]);

outDir = fullfile(HERE, 'Figures_combined_v3'); if ~exist(outDir, 'dir'), mkdir(outDir); end
base = fullfile(outDir, 'Normal_Forms_v3_summary');
print(fig, [base '.png'], '-dpng', '-r300');
try, savefig(fig, [base '.fig']); catch, end
fprintf('Saved %s.png\n', base);

%% ------------------------------------------------------------------------------------
function T = readCSV(f)
% minimal CSV reader (MATLAB and Octave): struct of columns, text as cellstr
txt = fileread(f);
L = regexp(txt, '\r?\n', 'split'); L = L(~cellfun(@isempty, L));
hdr = strsplit(L{1}, ',', 'CollapseDelimiters', false);
n = numel(L) - 1; C = cell(n, numel(hdr));
for i = 1:n
    r = L{i+1}; q = mod(cumsum(r == '"'), 2) == 1;  % commas inside quotes, e.g. "zeros(1,0)"
    r(r == ',' & q) = ';';
    v = strsplit(strrep(r, '"', ''), ',', 'CollapseDelimiters', false);
    C(i, 1:min(numel(v), numel(hdr))) = v(1:min(numel(v), numel(hdr)));
end
T = struct();
for j = 1:numel(hdr), T.(strtrim(hdr{j})) = C(:, j); end
end

function v = num(T, name, mask)
v = str2double(T.(name)); v = v(mask);
end
