%% Fig_inverse_method_v3.m
% Main-text Fig 3: the inverse (u-parametrisation) method for one parameter set
% (Iteration 1, row 96: S-dominant, Class 5, detached fold at ATOT = 19.49).
%
% Computes every number shown in the figure from the exact steady-state structure
% (support/ss_structure.m and the same equations, Eqs (7)-(9)):
%   r = KTOT/PTOT, alpha, beta, u* = -beta/alpha;
%   the grid in u (6 points per decade), A(u) and the excluded band r < u < u*;
%   the state at u = 1e4 (A, K, P, all substrate forms and their sum ATOT);
%   the fold of the detached branch (minimum of ATOT(u));
%   the three steady states at ATOT = ATOT(u = 1e4), with their stability.
% and draws the panels of Fig 3: (A) the branch followed by single-start continuation
% from ATOT = 0.1 (the main branch) and the detached branch it misses; (B1) the two
% factors of A(u) on the grid in u; (B2) A(u); (B3) the state at u = 1e4; (B4) ATOT(u)
% with the fold and the three roots; (B5) the complete diagram APP against ATOT.
%
% Input : Ensemble_iter1_v3_1e6/Excel_Files/ML_Ready_Dataset_v3_iter1.csv
% Output: Figures_combined_v3/Fig3_inverse_method.png/.pdf, Fig3_inverse_method_numbers.txt
% Runs in MATLAB (R2020b or newer) and GNU Octave; base functions only.

clear; close all;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
addpath(fullfile(HERE, 'support'));
OUT = fullfile(HERE, 'Figures_combined_v3'); if ~exist(OUT, 'dir'), mkdir(OUT); end
IT = 1; ROW = 96; UEX = 1e4;

E = load_ensemble_v3(HERE);
p = E.p(E.iteration == IT & E.row == ROW, :);
P = [1, p];                                    % [ATOT KTOT PTOT k1 ... k8]
C = model_constants(P);
f = @(u) (u - C.r) ./ (u .* (C.alpha*u + C.beta));   % A(u), Eq. (8)

% grid in u, 6 points per decade
u = 10.^(-1:1/6:6);
A = f(u); ok = A > 0 & isfinite(A);
[X, G] = ss_state_from_uA(u(ok), A(ok), C);

% example point
Aex = f(UEX); [Xex, Gex] = ss_state_from_uA(UEX, Aex, C);
Kex = C.KTOT/(1 + Aex/C.Km1 + C.a1*UEX*Aex/C.Km3); Pex = Kex/UEX;

% exact structure, detached fold and the steady states at ATOT = Gex
S = ss_structure(P, struct('ATOT_max', 1e3, 'N', 6000));
isdet = arrayfun(@(k) strcmp(S.branch(S.folds(k).branch).type, 'detached'), 1:numel(S.folds));
fold = S.folds(find(isdet, 1));
[Xs, bid] = ss_states_at(S, Gex);
stab = false(1, size(Xs,2));
for k = 1:size(Xs,2), stab(k) = all(real(eig(jac6(Xs(:,k), [Gex p]))) < 0); end
% u of each root (u = K/P)
uroot = zeros(1, size(Xs,2));
for k = 1:size(Xs,2)
    Kk = C.KTOT - Xs(3,k) - Xs(4,k); Pk = C.PTOT - Xs(5,k) - Xs(6,k); uroot(k) = Kk/Pk;
end
[uroot, o] = sort(uroot); Xs = Xs(:,o); stab = stab(o);
ufold = (C.KTOT - fold.x(3) - fold.x(4)) / (C.PTOT - fold.x(5) - fold.x(6));

fid = fopen(fullfile(OUT, 'Fig3_inverse_method_numbers.txt'), 'w');
for fo = [1 fid]
    fprintf(fo, 'Iteration %d, row %d: KTOT = %.2f, PTOT = %.2f, class %d\n', IT, ROW, p(1), p(2), E.class(E.iteration == IT & E.row == ROW));
    fprintf(fo, 'r = KTOT/PTOT = %.4g | u* = -beta/alpha = %.4g | alpha = %.4g | beta = %.4g\n', C.r, C.ustar, C.alpha, C.beta);
    fprintf(fo, 'grid: %d points (6 per decade), %d excluded (A < 0)\n', numel(u), sum(~ok));
    fprintf(fo, 'u = %g: A = %.4g, K = %.5g, P = %.4g\n', UEX, Aex, Kex, Pex);
    nm = {'AP','APP','AK','APK','APPC','APC'};
    for k = 1:6, fprintf(fo, '   %-5s = %.4g\n', nm{k}, Xex(k)); end
    fprintf(fo, '   sum = ATOT = %.4g\n', Gex);
    fprintf(fo, 'detached fold: u = %.4g, ATOT = %.4g, APP = %.4g\n', ufold, fold.ATOT, fold.x(2));
    for k = 1:numel(uroot)
        fprintf(fo, 'steady state %d at ATOT = %.4g: u = %.4g, APP = %.4g, %s\n', k, Gex, uroot(k), Xs(2,k), stab_word(stab(k)));
    end
end
fclose(fid);

% ---- figure (layout of main Fig 3)
W = 19; H = 23;
fig = figure('Color','w','Units','centimeters','Position',[1 1 W H], ...
             'PaperUnits','centimeters','PaperPositionMode','manual','PaperPosition',[0 0 W H],'PaperSize',[W H]);
lu = log10(u(ok)); cmap = jet(256);
co = cmap(1 + round(255*(lu - (-1))/7), :);          % marker colour = log10(u), 0.1 ... 1e6
uo = u(ok); Ao = A(ok); mo = uo < C.r; md = uo > C.ustar;
blue = [0.15 0.40 0.70]; grey = [0.55 0.55 0.55]; red = [0.84 0.24 0.24];
ib = find(strcmp({S.branch.type}, 'main'), 1); jb = find(strcmp({S.branch.type}, 'detached'), 1);

% overlay axes for headings, flow chart and footnote (normalised figure units)
ov = axes('Position',[0 0 1 1],'Visible','off'); hold(ov,'on'); axis(ov,[0 1 0 1]);
txt = @(x, y, s, varargin) text(ov, x, y, s, varargin{:});
txt(0.01, 0.98, '(A) Single-start continuation (MATCONT without the inverse method)', 'FontWeight','bold','FontSize',10);
txt(0.01, 0.74, '(B) Inverse method (exact labels): every steady state from a grid in u = K/P', 'FontWeight','bold','FontSize',10);
steps = {{'Fix','A_{TOT} = 0.1'}, {'Integrate ODEs','from X_0 = 0'}, {'Continue in A_{TOT}','(one run)'}, {'Label from the','branch followed'}};
for k = 1:4
    x0 = 0.02 + 0.145*(k-1);
    rectangle(ov,'Position',[x0 0.875 0.12 0.06],'Curvature',0.2,'EdgeColor',grey,'FaceColor',[0.96 0.96 0.96]);
    txt(x0 + 0.06, 0.905, steps{k}, 'HorizontalAlignment','center','FontSize',7);
    if k < 4, plot(ov, [x0+0.122 x0+0.143], [0.905 0.905], 'k-'); plot(ov, x0+0.139, 0.905, 'k>', 'MarkerFaceColor','k', 'MarkerSize',4); end
end
txt(0.02, 0.845, 'Outcome for this parameter set:', 'FontSize',7,'FontWeight','bold');
txt(0.03, 0.82, '- follows the main branch only (1 steady state for every A_{TOT})', 'FontSize',7);
txt(0.03, 0.797, sprintf('- misses the detached fold at A_{TOT} = %.2f and the bistable range A_{TOT} > %.2f', fold.ATOT, fold.ATOT), 'FontSize',7);
txt(0.5, 0.012, sprintf(['Iteration %d, row %d (K_{TOT} = %.2f, P_{TOT} = %.2f). Marker colour: log_{10} u (0.1 ... 10^6). ' ...
     'Filled, stable; open, unstable (analytic Jacobian). x, grid value with A < 0. a, b, c: steady states at A_{TOT} = %.2f.'], ...
     IT, ROW, p(1), p(2), Gex), 'HorizontalAlignment','center','FontSize',6);
axes('Position',[0.66 0.80 0.31 0.15]); hold on; box on;
Gm = S.branch(ib).ATOT; km = Gm >= 0.1;
plot(Gm(km), S.branch(ib).X(2,km), '-', 'Color', blue, 'LineWidth', 2);
plot(S.branch(jb).ATOT, S.branch(jb).X(2,:), ':', 'Color', grey, 'LineWidth', 1.2);
plot(fold.ATOT, fold.x(2), 'p', 'MarkerSize', 8, 'MarkerEdgeColor', grey);
text(fold.ATOT + 5, fold.x(2)/30, {'fold missed', sprintf('A_{TOT} = %.2f', fold.ATOT)}, 'FontSize',6,'Color',grey);
text(5, 3e-5, 'start, A_{TOT} = 0.1', 'FontSize',6,'Color',blue);
text(150, 4e-3, 'branch followed', 'FontSize',6,'Color',blue,'HorizontalAlignment','right');
text(150, 3e2, 'detached branch (not reached)', 'FontSize',6,'Color',grey,'HorizontalAlignment','right');
set(gca,'YScale','log','FontSize',7); xlim([0 150]); ylim([1e-6 1e3]); xlabel('A_{TOT}'); ylabel('APP');

% (1) the grid in u: the two factors of A(u)
T = @(y) sign(y).*log10(1 + abs(y));                  % signed log scale
axes('Position',[0.07 0.45 0.25 0.22]); hold on; box on;
ug = logspace(-1, 6, 400);
patch([C.r C.ustar C.ustar C.r], [-1 -1 7 7], [0.9 0.9 0.9], 'EdgeColor','none');
h1 = plot(log10(ug), T((ug - C.r)/C.r), 'k-', 'LineWidth', 1.2);
h2 = plot(log10(ug), T((C.alpha*ug + C.beta)/abs(C.beta)), '--', 'Color', red, 'LineWidth', 1.2);
for k = 1:numel(uo), plot(lu(k), 0, 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', co(k,:), 'MarkerEdgeColor', 'none'); end
plot(log10(UEX), 0, 'ko', 'MarkerSize', 8);
yt = [-1 0 1e2 1e4 1e6]; set(gca,'YTick',T(yt),'YTickLabel',{'-1','0','10^2','10^4','10^6'}, ...
    'XTick',-1:6,'XTickLabel',{'0.1','1','10','10^2','10^3','10^4','10^5','10^6'},'FontSize',7,'Layer','top');
xlim([-1 6]); ylim([T(-1.2) T(1e6)]);
xlabel('u = K/P (grid: 6 points per decade)'); ylabel('scaled factor of A(u)');
legend([h1 h2], {'(u - r)/r','(\alphau + \beta)/|\beta|'},'Location','northwest','FontSize',6,'Box','off');
text(5.9, T(5), {sprintf('r = K_{TOT}/P_{TOT} = %.2f', C.r), sprintf('u^* = -\\beta/\\alpha = %.2f', C.ustar), ...
     sprintf('\\alpha = %.3g', C.alpha), sprintf('\\beta = %.3g', C.beta)}, 'HorizontalAlignment','right','FontSize',6,'BackgroundColor','w','EdgeColor',grey);
title('(1) Choose the grid in u','FontSize',8);
% (2) A in closed form
axes('Position',[0.40 0.45 0.25 0.22]); hold on; box on;
patch([C.r C.ustar C.ustar C.r], [1e-4 1e-4 1e5 1e5], [0.9 0.9 0.9], 'EdgeColor','none');
plot(uo(mo), Ao(mo), '-', 'Color', grey); plot(uo(md), Ao(md), '-', 'Color', grey);
for k = 1:numel(uo), plot(uo(k), Ao(k), 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', co(k,:), 'MarkerEdgeColor', 'none'); end
plot(u(~ok), 2e-4*ones(1, sum(~ok)), 'kx', 'MarkerSize', 7);
plot(UEX, Aex, 'ko', 'MarkerSize', 8);
text(UEX*1.5, Aex*20, {'u = 10^4', sprintf('A = %.4f', Aex)}, 'FontSize',6);
text(0.15, 1e-2, {'main', 'branch (u < r)'}, 'FontSize',6); text(300, 3e-3, {'detached', 'branch (u > u^*)'}, 'FontSize',6,'HorizontalAlignment','center');
text(20, 2e-4, 'A < 0: no steady state', 'FontSize',6);
set(gca,'XScale','log','YScale','log','FontSize',7,'Layer','top'); xlim([0.1 1e6]); ylim([1e-4 1e5]);
xlabel('u'); ylabel('free substrate A'); title('(2) Solve for A in closed form','FontSize',8);
% (3) the state at u = 1e4
axes('Position',[0.77 0.45 0.20 0.22]); hold on; box on;
vals = [Aex; Xex(:); Gex]; lbl = {'A','AP','APP','AK','APK','APPC','APC','sum = A_{TOT}'};
for k = 1:8
    yk = 9 - k - (k == 8)*0.5; c = [0.35 0.75 0.45]; if k == 8, c = [0.15 0.15 0.15]; end
    patch(log10([1e-4 vals(k) vals(k) 1e-4]), yk + [-0.35 -0.35 0.35 0.35], c, 'EdgeColor','none');
    text(log10(vals(k)) + 0.15, yk, sprintf('%.3g', vals(k)), 'FontSize',6);
end
set(gca,'YTick',[0.5 2:8],'YTickLabel',fliplr(lbl),'XTick',-4:2:2, ...
    'XTickLabel',{'10^{-4}','10^{-2}','1','10^2'},'FontSize',7);
xlim([-4 3]); ylim([0 8.5]); xlabel('steady-state amount');
title({'(3) Rebuild the state (u = 10^4)', sprintf('free K = %.2f, free P = %.4f', Kex, Pex)},'FontSize',8);
% (4) ATOT(u) and its roots
axes('Position',[0.07 0.07 0.40 0.30]); hold on; box on;
patch([C.r C.ustar C.ustar C.r], [1 1 1e4 1e4], [0.9 0.9 0.9], 'EdgeColor','none');
plot(uo(mo), G(mo), '-', 'Color', grey); plot(uo(md), G(md), '-', 'Color', grey);
for k = 1:numel(uo), plot(uo(k), G(k), 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', co(k,:), 'MarkerEdgeColor', 'none'); end
plot([0.1 1e6], [Gex Gex], 'k--');
mk = {'s','s','s'}; ab = 'abc';
for k = 1:numel(uroot)
    fc = 'k'; if ~stab(k), fc = 'w'; end
    plot(uroot(k), Gex, mk{k}, 'MarkerSize', 7, 'MarkerFaceColor', fc, 'MarkerEdgeColor', 'k');
    text(uroot(k), Gex/1.6, ab(k), 'FontSize',7,'FontWeight','bold','HorizontalAlignment','center');
end
plot(ufold, fold.ATOT, 'p', 'MarkerSize', 9, 'MarkerFaceColor', red, 'MarkerEdgeColor', red);
text(ufold, fold.ATOT/2.2, {'fold (min of A_{TOT}(u))', sprintf('u = %.0f, A_{TOT} = %.2f', ufold, fold.ATOT)}, ...
     'Color', red, 'FontSize',6,'HorizontalAlignment','center');
tt = {sprintf('steady states at A_{TOT} = %.2f', Gex)};
for k = 1:numel(uroot), tt{end+1} = sprintf('%s  u = %.4g   APP = %.3g   %s', ab(k), uroot(k), Xs(2,k), stab_word(stab(k))); end %#ok<AGROW>
text(1.5e6/1.2, 9e3, tt, 'FontSize',6,'VerticalAlignment','top','HorizontalAlignment','right','BackgroundColor','w','EdgeColor',grey);
set(gca,'XScale','log','YScale','log','FontSize',7,'Layer','top'); xlim([0.1 1e6]); ylim([1 1e4]);
xlabel('u'); ylabel('A_{TOT}(u) = sum of all forms');
title('(4) Sum gives A_{TOT}(u); its roots are the steady states','FontSize',8);
% (5) the complete diagram
axes('Position',[0.57 0.07 0.40 0.30]); hold on; box on;
patch([fold.ATOT 150 150 fold.ATOT], [1e-7 1e-7 1e3 1e3], [0.96 0.93 0.88], 'EdgeColor','none');
for b = 1:numel(S.branch)
    Gb = S.branch(b).ATOT; Xb = S.branch(b).X;
    step = max(1, floor(numel(Gb)/600)); idx = 1:step:numel(Gb); st = false(size(idx));
    for k = 1:numel(idx), st(k) = all(real(eig(jac6(Xb(:,idx(k)), [Gb(idx(k)) p]))) < 0); end
    Gs = Gb(idx); Xp = Xb(2, idx); Xu = Xp; Xp(~st) = NaN; Xu(st) = NaN;
    plot(Gs, Xp, '-', 'Color', blue, 'LineWidth', 1.5); plot(Gs, Xu, '--', 'Color', blue, 'LineWidth', 1.2);
end
for k = 1:numel(uo)
    if G(k) <= 150, plot(G(k), X(2,k), 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', co(k,:), 'MarkerEdgeColor', 'none'); end
end
for k = 1:numel(uroot)
    fc = 'k'; if ~stab(k), fc = 'w'; end
    plot(Gex, Xs(2,k), 's', 'MarkerSize', 7, 'MarkerFaceColor', fc, 'MarkerEdgeColor', 'k');
    text(Gex - 4, Xs(2,k), ab(k), 'FontSize',7,'FontWeight','bold');
end
plot(fold.ATOT, fold.x(2), 'p', 'MarkerSize', 9, 'MarkerFaceColor', red, 'MarkerEdgeColor', red);
text(fold.ATOT + 8, fold.x(2)/2.5, sprintf('LP: A_{TOT} = %.2f, APP = %.2f', fold.ATOT, fold.x(2)), 'Color', red, 'FontSize',6);
plot([Gex Gex], [1e-7 1e3], 'k:');
text(3, 3e-7, '1 stable', 'FontSize',6); text(85, 3e-7, 'bistable: 3 steady states (2 stable)', 'FontSize',6,'HorizontalAlignment','center');
set(gca,'YScale','log','FontSize',7,'Layer','top'); xlim([0 150]); ylim([1e-7 1e3]);
xlabel('A_{TOT}'); ylabel('APP'); title('(5) Place every grid point: the complete diagram','FontSize',8);
print(fig, fullfile(OUT, 'Fig3_inverse_method.png'), '-dpng', '-r300');
print(fig, fullfile(OUT, 'Fig3_inverse_method.pdf'), '-dpdf');
fprintf('Fig 3 written to %s\n', OUT);

function s = stab_word(b)
if b, s = 'stable'; else, s = 'unstable'; end
end
