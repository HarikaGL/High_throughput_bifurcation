% ======================= TR_04_figures.m ===============================================
% Supplementary figures and table of the class transitions (needs TR_02 and TR_03 results).
%   S_transition_map      transition graph: classes as (x, y) scatter panels (as in
%                         Fig 11), edges = exact conditions
%                         (solid), main-branch fold (dashed), symmetry breaking (dotted),
%                         labelled with the controlling quantity, the two-class sweep
%                         result and the best single-parameter move of the scan
%   S_two_class_switches  for every edge: exact diagrams in the two classes and the class
%                         along the whole sweep of the controlling quantity (one set)
%   S_symmetry_breaking   Class-6 sets before and after one move that breaks u* = r
%   TR_edges_table.tex    LaTeX rows of the Supplementary table of transitions
% Output: ../Figures_combined_v3/Class_transitions/*.png/.pdf/.fig and results/
% ======================================================================================
close all;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
addpath(fullfile(fileparts(HERE), 'support'), fullfile(HERE, 'tr_lib'));
resDir = fullfile(HERE, 'results');
figDir = fullfile(fileparts(HERE), 'Figures_combined_v3', 'Class_transitions');
if ~exist(figDir, 'dir'), mkdir(figDir); end
FS = 10;
pn = {'K_{TOT}','P_{TOT}','k_1','k_{11}','k_2','k_3','k_{33}','k_4','k_5','k_{55}','k_6','k_7','k_{77}','k_8','k_4&k_8','k_1&k_{11}'};
SC   = tr_read_csv(fullfile(resDir, 'TR_scan.csv'));          % iter,src,set,move,d,dst,x,y,crit10,detfold
EX   = tr_read_csv(fullfile(resDir, 'TR_exact_moves.csv'));   % iter,src,set,move_id,value,dst,detfold
BEST = tr_read_csv(fullfile(resDir, 'TR_best_moves.csv'));    % iter,src,dst,move,d,fraction
SW   = tr_read_csv(fullfile(resDir, 'TR_sweep_summary.csv')); % spec,iter,A,B,param,shown,sets,clean,switch
load(fullfile(resDir, 'TR_two_class_examples.mat'), 'EXM');
E = tr_load_ensemble(HERE);

%% ======================= edges ==========================================================
% kind: 1 exact, 2 fold, 3 symmetry breaking. spec: two-class sweep (TR_03) or 0
%       from to iter kind spec  quantity
EDG = {  1,  2, 1, 1, 1, 'x < 0';
        10,  5, 1, 1, 2, 'x < 0';
         3,  4, 2, 1, 3, 'y > 0';
         9,  5, 2, 1, 4, 'y > 0';
         1, 10, 1, 2, 5, 'fold (S)';
         2,  5, 1, 2, 6, 'fold (S)';
         3,  9, 2, 2, 7, 'fold (Z)';
         4,  5, 2, 2, 8, 'fold (Z)';
        10,  6, 3, 1, 9, 'k_4 > k_2';
         6,  7, 3, 4, 10, 'K_{M5} < K_{M1}';
         6,  8, 3, 4, 10, 'K_{M5} > K_{M1}';
         6,  5, 3, 3, 0, 'u^* \neq r';
         6,  9, 3, 3, 0, 'u^* \neq r';
         6, 10, 3, 3, 0, 'u^* \neq r'};
nE = size(EDG, 1); best = cell(nE, 2); bestR = cell(nE, 2); clean = nan(nE, 1);
for e = 1:nE
    [a, b, it, kind, sp] = EDG{e, 1:5};
    % best single move in both directions (same iteration)
    for dirn = 1:2
        if dirn == 1, s = a; t = b; else, s = b; t = a; end
        k = BEST(:,1) == it & BEST(:,2) == s & BEST(:,3) == t;
        if ~any(k) && it == 3, k = BEST(:,1) == 4 & BEST(:,2) == s & BEST(:,3) == t; end
        if any(k)
            r = BEST(find(k, 1), :); bestR{e, dirn} = r;
            m = sprintf('%g', 10^abs(r(5))); if abs(abs(r(5)) - 0.5) < 1e-9, m = '10^{0.5}'; end
            if r(5) > 0, op = '\times'; else, op = '\div'; end
            best{e, dirn} = sprintf('%s %s%s: %.0f%%', pn{r(4)}, op, m, 100*r(6));
        else
            best{e, dirn} = '--';
        end
    end
    if sp > 0, clean(e) = 100*SW(SW(:,1) == sp, 8); end
end

%% ======================= figure 1: transition map =======================================
% Every class is drawn as its own (x, y) scatter panel (as in Fig 11: the class in colour
% against the ensemble in grey, pink = detached-branch quadrant); edges between the panels.
pos = containers.Map('KeyType', 'double', 'ValueType', 'any');
pos(1) = [1.5 3.0];  pos(2) = [5.2 3.0];  pos(4) = [8.9 3.0];  pos(3) = [12.6 3.0];
pos(10) = [1.5 0.0]; pos(5) = [7.05 0.0]; pos(9) = [12.6 0.0];
pos(106) = [16.2 3.0]; pos(6) = [16.2 0.0]; pos(7) = [15.0 -3.0]; pos(8) = [17.4 -3.0];   % 106: Class 10 of Iterations 3-4
XL = [0 18.6]; YL = [-4.4 4.55]; nwD = 1.9; nhD = 1.55;        % node size in data units
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [1 1 46 32]);
axP = [0.01 0.065 0.98 0.89];
ax = axes('Position', axP); hold(ax, 'on'); axis(ax, 'off');
set(ax, 'XLim', XL, 'YLim', YL);
toFig = @(P) [axP(1) + (P(1) - XL(1))/diff(XL)*axP(3), axP(2) + (P(2) - YL(1))/diff(YL)*axP(4)];
patch(ax, [0.2 14.2 14.2 0.2], [-1.25 -1.25 4.45 4.45], [0.965 0.965 0.98], 'EdgeColor', [0.8 0.8 0.85]);
patch(ax, [14.4 18.5 18.5 14.4], [-4.35 -4.35 4.45 4.45], [0.98 0.965 0.965], 'EdgeColor', [0.85 0.8 0.8]);
text(ax, 7.2, 4.62, 'Iterations 1-2: catalytic dominance at the two sites', 'HorizontalAlignment', 'center', 'FontSize', FS+2, 'FontWeight', 'bold');
text(ax, 16.45, 4.62, 'Iterations 3-4: symmetry', 'HorizontalAlignment', 'center', 'FontSize', FS+2, 'FontWeight', 'bold');
sty = {'-', '--', ':', '-'}; ecol = {[0.15 0.15 0.15], [0.35 0.35 0.35], [0.55 0.2 0.65], [0.15 0.15 0.15]};   % kind 4: exact, one way
% label offset (data units) from the middle of each edge
LO = [0 0.55; 0 -0.62; 0 0.55; 0 -0.62; 0.95 0; -0.75 0.15; 0.95 0; 0.75 0.15; 0.95 0; -0.95 0.2; 0.95 0.2; 0 0; 0 0; 0 0];
for e = 1:nE
    [a, b, it, kind] = EDG{e, 1:4};
    ka = a; kb = b; if it >= 3 && a == 10, ka = 106; end; if it >= 3 && b == 10, kb = 106; end
    if kind == 3, kb = b; end
    P = pos(ka); Qp = pos(kb); v = Qp - P; u = v/norm(v);
    cut = @(w) min((nwD/2 + 0.08)/max(abs(w(1)), eps), (nhD/2 + 0.2)/max(abs(w(2)), eps));
    if kind == 3 && any(b == [5 10])          % long symmetry-breaking edges: curved below
        C = [(P(1) + Qp(1))/2, -2.9 - 0.7*(b == 10)];
        t = linspace(0, 1, 80)';
        B = (1-t).^2*P + 2*(1-t).*t*C + t.^2*Qp;
        inside = abs(B(:,1) - P(1)) < nwD/2 + 0.08 & abs(B(:,2) - P(2)) < nhD/2 + 0.2 | ...
                 abs(B(:,1) - Qp(1)) < nwD/2 + 0.08 & abs(B(:,2) - Qp(2)) < nhD/2 + 0.2;
        B = B(~inside, :);
        plot(ax, B(:,1), B(:,2), sty{kind}, 'Color', ecol{kind}, 'LineWidth', 2.2);
        ut = B(end, :) - B(end-1, :); ut = ut/norm(ut); tr_arrowhead(ax, B(end, :), ut, ecol{kind});
        if b == 5, M = B(round(0.40*size(B,1)), :) + [0 -0.32]; else, M = B(round(0.5*size(B,1)), :) + [0 -0.32]; end
    else
        P1 = P + cut(u)*u; Q1 = Qp - cut(u)*u;
        plot(ax, [P1(1) Q1(1)], [P1(2) Q1(2)], sty{kind}, 'Color', ecol{kind}, 'LineWidth', 2.2);
        tr_arrowhead(ax, Q1, u, ecol{kind});
        if kind <= 2, tr_arrowhead(ax, P1, -u, ecol{kind}); end
        M = (P1 + Q1)/2;
    end
    lab = {['\bf' EDG{e, 6}]};
    if kind ~= 3 && isfinite(clean(e)), lab{end+1} = sprintf('two-class sweep: %.0f%%', clean(e)); end
    if kind <= 2, lab{end+1} = sprintf('%d\\rightarrow%d %s', a, b, best{e,1}); lab{end+1} = sprintf('%d\\rightarrow%d %s', b, a, best{e,2});
    else, lab{end+1} = sprintf('%d\\rightarrow%d %s', a, b, best{e,1}); end
    M = M + LO(e, :);
    text(ax, M(1), M(2), lab, 'HorizontalAlignment', 'center', 'FontSize', FS-2, ...
         'BackgroundColor', 'w', 'Margin', 1, 'Color', ecol{kind});
end
% nodes: one scatter panel per class
CI = class_info_v3(); Q = tr_conditions(E.p);
gi = find(E.use); gi = gi(tr_pick(numel(gi), min(4000, numel(gi)), 11));
for kk = cell2mat(keys(pos))
    c = kk; if kk == 106, c = 10; end
    P = pos(kk); F = toFig(P);
    w = nwD/diff(XL)*axP(3); h = nhD/diff(YL)*axP(4);
    axN = axes('Position', [F(1) - w/2, F(2) - h/2, w, h]);
    kc = find(E.use & E.class == c);
    if kk == 10, kc = kc(E.iteration(kc) <= 2); elseif kk == 106, kc = kc(E.iteration(kc) >= 3); end
    kc = kc(tr_pick(numel(kc), min(2500, numel(kc)), 100 + kk));
    msz = 3; if any(kk == [106 6 7 8]), msz = 6; end   % symmetric sets lie on the line y = -x
    tr_plane_panel(axN, Q, gi, kc, tr_colors(c), FS-4, msz);
    set(axN, 'XTick', [-4 -2 0 2 4], 'YTick', [-4 -2 0 2 4 6], 'LineWidth', 1.4, 'XColor', tr_colors(c), 'YColor', tr_colors(c));
    nm = strrep(CI.labels{c}, '_2^{LP}', '_2(LP)'); nm = strrep(nm, '_2^{BP}', '_2(BP)');
    sub = ''; if kk == 10, sub = ' (Iter. 1-2)'; elseif kk == 106, sub = ' (Iter. 3-4)'; end
    title(axN, sprintf('%d: %s%s', c, nm, sub), 'FontSize', FS-1, 'FontWeight', 'bold', 'Color', tr_colors(c));
    plot(axN, -4.6, 5.0, 'o', 'MarkerSize', 13, 'MarkerFaceColor', tr_colors(c), 'MarkerEdgeColor', 'w');
    text(axN, -4.6, 5.0, sprintf('%d', c), 'HorizontalAlignment', 'center', 'FontSize', FS-1, 'FontWeight', 'bold', 'Color', 'w');
end
axL = axes('Position', [0 0 1 0.05]); axis(axL, 'off'); hold(axL, 'on'); set(axL, 'XLim', [0 1], 'YLim', [0 1]);
lt = {'exact condition (closed form)', 'fold on the main branch (no closed form)', 'symmetry breaking (any single change with u^* \neq r)'};
for q = 1:3
    plot(axL, 0.05 + 0.31*(q-1) + [0 0.04], [0.65 0.65], sty{q}, 'Color', ecol{q}, 'LineWidth', 2.2);
    text(axL, 0.1 + 0.31*(q-1), 0.65, lt{q}, 'FontSize', FS);
end
text(axL, 0.5, 0.18, ['Panels: each class (colour) against the ensemble (grey) in the plane x = log_{10}(k_2K_{TOT}/k_8P_{TOT}), y = log_{10}(k_4K_{TOT}/k_6P_{TOT}) (pink: x < 0, y > 0). ', ...
     'Labels: controlling quantity; two-class sweep (TR\_03); best single-parameter change in each direction (TR\_02).'], 'HorizontalAlignment', 'center', 'FontSize', FS-2);
tr_save(fig, fullfile(figDir, 'S_transition_map'));

%% ======================= figure 2: two-class switches ====================================
show = find(SW(:,6) == 1)';
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [1 1 42 52]);
qName = struct('x', 'x = log_{10}(k_2K_{TOT}/k_8P_{TOT})', 'y', 'y = log_{10}(k_4K_{TOT}/k_6P_{TOT})', ...
               'mul', 'log_{10}(k_3/k_7)', 'bp', 'log_{10}(k_4/k_2), k_4 = k_8', 'km', 'log_{10}(K_{M1}/K_{M5})');
for q = 1:numel(show)
    s = show(q); X = EXM(s);
    col = mod(q-1, 2); row = floor((q-1)/2);
    x0 = 0.04 + 0.5*col; y0 = 0.81 - 0.192*row; w = 0.13;
    via = pn{X.pj}; if strcmp(X.typ, 'mul'), via = pn{X.pj}; end
    ttl = sprintf('%d \\leftrightarrow %d   %s, via %s   (two-class in %.0f%% of sets)', X.A, X.B, qName.(X.typ), via, 100*X.clean);
    axT = axes('Position', [x0, y0 + 0.155, 0.45, 0.02]); axis(axT, 'off');
    text(axT, 0, 0.5, ttl, 'FontSize', FS, 'FontWeight', 'bold');
    for side = 1:2
        if side == 1, p = X.pA; c = X.A; v = X.vA; else, p = X.pB; c = X.B; v = X.vB; end
        a = axes('Position', [x0 + (side-1)*0.155, y0 + 0.012, w, 0.125]);
        D = tr_curves(p, 1e-2, 1e6, 500);
        opt = struct('lw', 1.4, 'ms', 5);
        if any(c == [7 8]) && ~isempty(D.lp) && ~isempty(D.bp) && abs(log10(D.lp(1,1)/D.bp(1,1))) < 0.3
            opt.zoom = tr_bp_zoom(p, 0.3);
        end
        tr_draw_diagram(a, D, opt); tr_tick_decades(a);
        set(a, 'FontSize', FS-3);
        title(a, sprintf('Class %d (value %+.2f)', D.class, v), 'FontSize', FS-1, 'Color', tr_colors(c));
        xlabel(a, 'A_{TOT}', 'FontSize', FS-2); if side == 1, ylabel(a, 'x_2', 'FontSize', FS-2); end
    end
    a = axes('Position', [x0 + 0.315, y0 + 0.012, w, 0.125]); hold(a, 'on'); box(a, 'on');
    st = sortrows(X.strip, 1); vv = st(:,1); cc = st(:,2);
    mids = [vv(1); (vv(1:end-1) + vv(2:end))/2; vv(end)];
    for i = 1:numel(vv)
        patch(a, [mids(i) mids(i+1) mids(i+1) mids(i)], [0 0 1 1], tr_colors(min(max(cc(i),1),10)), 'EdgeColor', 'none');
    end
    run = [1; find(diff(cc) ~= 0) + 1; numel(cc) + 1];
    for i = 1:numel(run)-1
        xm = (mids(run(i)) + mids(run(i+1)))/2;
        text(a, xm, 0.55, sprintf('%d', cc(run(i))), 'HorizontalAlignment', 'center', 'FontSize', FS+1, 'FontWeight', 'bold', 'Color', 'w');
    end
    plot(a, [X.vA X.vB], [0.1 0.1], 'k^', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    set(a, 'XLim', [vv(1) vv(end)], 'YLim', [0 1], 'YTick', [], 'FontSize', FS-3, 'Layer', 'top');
    xlabel(a, qName.(X.typ), 'FontSize', FS-2); title(a, 'class along the sweep', 'FontSize', FS-1);
end
axL = axes('Position', [0 0.005 1 0.02]); axis(axL, 'off'); set(axL, 'XLim', [0 1], 'YLim', [0 1]);
text(axL, 0.5, 0.5, ['Each row: one parameter set; only the named quantity is changed. Diagrams at the two marked values (\^); ', ...
     'the bar shows the exact class over the whole sweep. Solid stable, dashed unstable; diamonds folds, dots branch points.'], ...
     'HorizontalAlignment', 'center', 'FontSize', FS-1);
tr_save(fig, fullfile(figDir, 'S_two_class_switches'));

%% ======================= figure 3: symmetry breaking =====================================
MV = {6, -0.5, 5, 'k_3 \times 10^{-0.5}'; 1, -1, 9, 'K_{TOT} \div 10'; 1, 1, 10, 'K_{TOT} \times 10'};
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [1 1 42 26]);
for m = 1:size(MV, 1)
    yA = 0.69 - 0.32*(m-1);
    k = SC(:,1) == 3 & SC(:,2) == 6 & SC(:,4) == MV{m,1} & SC(:,5) == MV{m,2};
    dist = SC(k, 6); u = unique(dist)';
    ds = strjoin(arrayfun(@(v) sprintf('%d: %.0f%%', v, 100*mean(dist == v)), u, 'UniformOutput', false), ', ');
    axT = axes('Position', [0.03, yA + 0.245, 0.9, 0.03]); axis(axT, 'off');
    text(axT, 0, 0.5, sprintf('6 \\rightarrow %d:  %s (symmetry broken)    120 Class-6 sets after the move: %s', MV{m,3}, MV{m,4}, ds), ...
         'FontSize', FS+1, 'FontWeight', 'bold');
    sets = SC(k & SC(:,6) == MV{m,3}, 3); sets = sets(1:min(3, end));
    for s = 1:numel(sets)
        p0 = E.p(sets(s), :);
        for side = 1:2
            if side == 1, p = p0; else, p = tr_move(p0, 'mul', MV{m,1}, MV{m,2}); end
            a = axes('Position', [0.04 + 0.16*(2*(s-1) + side - 1), yA, 0.13, 0.19]);
            D = tr_curves(p, 1e-2, 1e6, 500);
            tr_draw_diagram(a, D, struct('lw', 1.4, 'ms', 5)); tr_tick_decades(a); set(a, 'FontSize', FS-3);
            if side == 1, t = 'before'; else, t = 'after'; end
            title(a, sprintf('set %d %s: Class %d', s, t, D.class), 'FontSize', FS-1, 'Color', tr_colors(D.class));
            if m == size(MV,1), xlabel(a, 'A_{TOT}', 'FontSize', FS-2); end
        end
    end
end
tr_save(fig, fullfile(figDir, 'S_symmetry_breaking'));

%% ======================= table rows =======================================================
fid = fopen(fullfile(resDir, 'TR_edges_table.tex'), 'w');
kn = {'exact', 'fold', 'symmetry', 'exact'};
pt = {'\KTOT','\PTOT','k_1','k_{11}','k_2','k_3','k_{33}','k_4','k_5','k_{55}','k_6','k_7','k_{77}','k_8','k_4, k_8','k_1, k_{11}'};
for e = 1:nE
    q = EDG{e, 6}; q = strrep(q, 'K_{M5}', '\KM{5}'); q = strrep(q, 'K_{M1}', '\KM{1}');
    if strncmp(q, 'fold', 4), q = ['main-branch ' q]; else, q = ['$' q '$']; end
    cs = '--'; if isfinite(clean(e)), cs = sprintf('%.0f', clean(e)); end
    if EDG{e, 4} <= 2, ar = '\leftrightarrow'; else, ar = '\to'; end
    mv = cell(1, 2);
    for dirn = 1:2
        r = bestR{e, dirn};
        if isempty(r) || EDG{e, 4} > 2 && dirn == 2, mv{dirn} = '--'; continue; end
        if abs(abs(r(5)) - 0.5) < 1e-9, m = '10^{0.5}'; else, m = sprintf('%g', 10^abs(r(5))); end
        if r(5) > 0, op = '\times'; else, op = '\div'; end
        mv{dirn} = sprintf('$%s\\,%s\\,%s$ (%.0f\\%%)', pt{r(4)}, op, m, 100*r(6));
    end
    fprintf(fid, '$%d %s %d$ & %s & %s & %s & %s & %s \\\\\n', EDG{e,1}, ar, EDG{e,2}, kn{EDG{e,4}}, q, cs, mv{1}, mv{2});
end
fclose(fid);
fprintf('Figures in %s\n', figDir);
