%% Fig_criterion_v3.m
% Main-text Fig 7: criterion (10) against the full steady-state computation.
%
%   (A) the 20,000 parameter sets of Iterations 1 and 2 on the two catalytic-dominance
%       ratios x = log10(k8*PTOT/k2*KTOT) and y = log10(k4*KTOT/k6*PTOT), coloured by
%       whether a detached branch exists at any ATOT (computed from the interval structure
%       of A(u), independently of the criterion);
%   (B) the same plane coloured by the exact class; open circles: detached branch whose
%       fold lies above the labelling window (ATOT > 1e6);
%   (C) criterion (10) as a decision tree with the number of sets at every node.
%
% Input : Ensemble_iter<i>_v3_1e6/Excel_Files/ML_Ready_Dataset_v3_iter<i>.csv (i = 1, 2)
% Output: Figures_combined_v3/Fig7_criterion.png/.pdf and Fig7_criterion_counts.txt
% Runs in MATLAB (R2020b or newer) and GNU Octave; base functions only.

clear; close all;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
addpath(fullfile(HERE, 'support'));
OUT = fullfile(HERE, 'Figures_combined_v3'); if ~exist(OUT, 'dir'), mkdir(OUT); end

E  = load_ensemble_v3(HERE);
in = E.iteration <= 2;
P  = E.p(in,:);  cls = E.class(in);  n = sum(in);
KT = P(:,1); PT = P(:,2); k2 = P(:,5); k4 = P(:,8); k6 = P(:,11); k8 = P(:,14);
x  = log10(k8.*PT ./ (k2.*KT));          % > 0: phosphatase dominates site 1
y  = log10(k4.*KT ./ (k6.*PT));          % > 0: kinase dominates site 2

% detached branch at any ATOT, from the sign pattern of A(u) (ss_structure.m logic)
det = false(n,1);
for i = 1:n
    C = model_constants([1, P(i,:)]);
    c = [0, C.r];
    if C.alpha ~= 0 && -C.beta/C.alpha > 0, c(end+1) = -C.beta/C.alpha; end %#ok<AGROW>
    c = [unique(c), Inf];
    for j = 1:numel(c)-1
        lo = c(j); hi = c(j+1);
        if isinf(hi), u = 10*max(lo,1); elseif lo == 0, u = hi/10; else, u = sqrt(lo*hi); end
        A = (u - C.r) / (u*(C.alpha*u + C.beta));
        if A > 0 && lo ~= C.r && hi ~= C.r, det(i) = true; end
    end
end
crit  = x > 0 & y > 0;                   % criterion (10)
inWin = ismember(cls, [2 4 5]);          % detached fold inside ATOT <= 1e6
above = det & ~inWin;                    % detached fold above the window

%% counts (panel C and caption)
fid = fopen(fullfile(OUT, 'Fig7_criterion_counts.txt'), 'w');
pr = @(varargin) fprintf(fid, varargin{:});
pr('Iterations 1-2: %d sets\n', n);
pr('criterion (10) agrees with the computation: %d of %d\n', sum(crit == det), n);
pr('detached branch: %d | fold in window: %d | fold above window: %d\n', sum(det), sum(det & inWin), sum(above));
pr('  in window, by class: 5 %d, 2 %d, 4 %d\n', sum(det & cls==5), sum(det & cls==2), sum(det & cls==4));
pr('  above window, by class: 1 %d, 3 %d, 9 %d, 10 %d\n', sum(above & cls==1), sum(above & cls==3), sum(above & cls==9), sum(above & cls==10));
kin = x < 0 & y > 0;  pho = x > 0 & y < 0;
pr('kinase dominates both: %d (class 1 %d, class 10 %d)\n', sum(kin), sum(kin & cls==1), sum(kin & cls==10));
pr('phosphatase dominates both: %d (class 3 %d, class 9 %d)\n', sum(pho), sum(pho & cls==3), sum(pho & cls==9));
pr('lower-left quadrant: %d\n', sum(x < 0 & y < 0));
fclose(fid);
type(fullfile(OUT, 'Fig7_criterion_counts.txt'));

%% figure
CI = class_info_v3(); col = CI.colors;
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 34]);
lim = [-4 6 -2 6];
axA = axes('Position',[0.07 0.60 0.40 0.34]);  axB = axes('Position',[0.56 0.60 0.40 0.34]);
for ax = [axA axB]
    axes(ax); hold on; box on;
    patch([-4 0 0 -4], [-2 -2 0 0], [0.93 0.93 0.93], 'EdgeColor','none');
    text(-2, -1, 'not sampled', 'HorizontalAlignment','center', 'Color',[0.45 0.45 0.45], 'FontSize',10);
end
axes(axA);
scatter(x(~det), y(~det), 4, [0.72 0.72 0.70], 'filled');
scatter(x(det),  y(det),  4, [0.84 0.24 0.24], 'filled');
text(5.8, 5.6, {'k_8P_{TOT} > k_2K_{TOT}', 'k_4K_{TOT} > k_6P_{TOT}', sprintf('detached branch: %s', commas(sum(det)))}, ...
     'HorizontalAlignment','right','VerticalAlignment','top','FontSize',10,'EdgeColor',[0.84 0.24 0.24],'BackgroundColor','w');
text(-3.8, 5.6, sprintf('no detached branch: %s', commas(sum(kin))), 'VerticalAlignment','top','FontSize',10,'EdgeColor',[0.6 0.6 0.6],'BackgroundColor','w');
text(5.8, -1.6, sprintf('no detached branch: %s', commas(sum(pho))), 'HorizontalAlignment','right','FontSize',10,'EdgeColor',[0.6 0.6 0.6],'BackgroundColor','w');
title(sprintf('(A) Criterion (10) vs full computation (n = %s)', commas(n)), 'FontSize',12);
xlabel(sprintf('agreement: %s of %s sets', commas(sum(crit == det)), commas(n)), 'FontSize',10);
axes(axB);
order = [10 1 2 5 4 3 9];
for c = order
    k = cls == c; scatter(x(k), y(k), 4, col(c,:), 'filled');
end
plot(x(above), y(above), 'ko', 'MarkerSize', 3.5, 'LineWidth', 0.5);
text(-1.7, 5.4, {'kinase dominates', 'both sites'}, 'HorizontalAlignment','center','FontWeight','bold','FontSize',10);
text(3.0, 5.4, 'crossed dominance', 'HorizontalAlignment','center','FontWeight','bold','FontSize',10);
text(3.0, -1.6, 'phosphatase dominates both', 'HorizontalAlignment','center','FontWeight','bold','FontSize',10);
title('(B) Exact class in the same plane', 'FontSize',12);
for ax = [axA axB]
    axes(ax);
    plot([0 0], lim(3:4), 'k-', 'LineWidth', 1.2); plot(lim(1:2), [0 0], 'k-', 'LineWidth', 1.2);
    axis(lim); set(ax, 'FontSize', 10, 'Layer', 'top');
    ylabel('log_{10}(k_4K_{TOT}/k_6P_{TOT})  (site 2: > 0, kinase dominates)');
end
annotation('textbox', [0.07 0.545 0.40 0.03], 'String', ...
    'log_{10}(k_8P_{TOT}/k_2K_{TOT})  (site 1: > 0, phosphatase dominates)', ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',10);
annotation('textbox', [0.56 0.545 0.40 0.03], 'String', ...
    'log_{10}(k_8P_{TOT}/k_2K_{TOT})  (site 1: > 0, phosphatase dominates)', ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',10);

% (C) decision tree
axC = axes('Position',[0.03 0.02 0.94 0.50]); axis(axC, [0 100 0 100]); axis(axC, 'off'); hold(axC, 'on');
text(0, 98, '(C) Criterion (10) as a decision tree: Iterations 1-2', 'FontWeight','bold','FontSize',12);
node = @(xc, yc, w, s, ec) rectangle('Position',[xc-w/2 yc-3 w 6], 'Curvature',0.15, 'EdgeColor',ec, 'LineWidth',1.2, 'FaceColor','w');
lab  = @(xc, yc, s) text(xc, yc, s, 'HorizontalAlignment','center', 'FontSize',9);
red = [0.84 0.24 0.24]; grey = [0.55 0.55 0.55];
node(8, 48, 14, '', 'k');      lab(8, 48, {'Iterations 1-2', sprintf('%s sets', commas(n))});
node(30, 70, 18, '', red);     lab(30, 70, {'crossed dominance', sprintf('detached branch: %s', commas(sum(det)))});
node(30, 30, 18, '', grey);    lab(30, 30, {'kinase dominates both', sprintf('%s', commas(sum(kin)))});
node(30, 12, 18, '', grey);    lab(30, 12, {'phosphatase dominates both', sprintf('%s', commas(sum(pho)))});
node(54, 82, 16, '', red);     lab(54, 82, {'fold at A_{TOT} \leq 10^6', commas(sum(det & inWin))});
node(54, 58, 16, '', red);     lab(54, 58, {'fold at A_{TOT} > 10^6', commas(sum(above))});
leaves = {5, det & cls==5, 90; 2, det & cls==2, 84; 4, det & cls==4, 78; ...
          1, above & cls==1, 66; 3, above & cls==3, 60; 9, above & cls==9, 54; 10, above & cls==10, 48; ...
          1, kin & cls==1, 34; 10, kin & cls==10, 28; 3, pho & cls==3, 16; 9, pho & cls==9, 10};
parent = [82 82 82 58 58 58 58 30 30 12 12];
for i = 1:size(leaves,1)
    c = leaves{i,1}; yy = leaves{i,3};
    rectangle('Position',[70 yy-2.4 22 4.8], 'Curvature',0.2, 'EdgeColor',col(c,:)*0.8, 'FaceColor',1-0.25*(1-col(c,:)), 'LineWidth',1.2);
    text(71, yy, sprintf('Class %d', c), 'FontSize',9, 'FontWeight','bold');
    text(91, yy, commas(sum(leaves{i,2})), 'FontSize',9, 'HorizontalAlignment','right', 'FontWeight','bold');
    xp = 62 - 24*(parent(i) <= 30);  % parent box right edge
    plot([xp 66 66 70], [parent(i) parent(i) yy yy], '-', 'Color', grey);
end
plot([15 18 18 21], [48 48 70 70], 'k-'); plot([18 18 21], [48 30 30], 'k-'); plot([18 18 21], [30 12 12], 'k-');
plot([39 42 42 46], [70 70 82 82], '-', 'Color', red); plot([42 42 46], [70 58 58], '-', 'Color', red);
text(99, 2, 'With K_{TOT}, P_{TOT} free, criterion (10) can be met iff k_4k_8 > k_2k_6.', ...
     'HorizontalAlignment','right','FontSize',9,'Color',[0.35 0.35 0.35]);

print(fig, fullfile(OUT, 'Fig7_criterion.png'), '-dpng', '-r300');
print(fig, fullfile(OUT, 'Fig7_criterion.pdf'), '-dpdf');
fprintf('Fig 7 written to %s\n', OUT);

function s = commas(v)
s = sprintf('%d', round(v));
k = numel(s) - 3;
while k > 0, s = [s(1:k) ',' s(k+1:end)]; k = k - 3; end
end
