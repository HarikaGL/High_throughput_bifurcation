function C = matcont_check(paramSamples, res, cfg, oldOpt, nCheck, seed, figDir, useParallel)
%MATCONT_CHECK  Visual performance check of MATCONT on randomly chosen samples.
%
%   C = matcont_check(paramSamples, res, cfg, oldOpt, nCheck, seed, figDir, useParallel)
%
%   For nCheck randomly chosen samples (exact analysis successful), MATCONT is run
%   twice and drawn in grey OVER the exact branches (colour; solid = stable,
%   dashed = unstable):
%
%   Figure A  "MATCONT_check_WITH_inverse_method"
%       MATCONT seeded on EVERY branch of steady states found by the inverse method
%       (analyze_sample, both directions). Shows that MATCONT reproduces all branches,
%       folds and branch points when it is told where the branches are.
%   Figure B  "MATCONT_check_WITHOUT_inverse_method"
%       The original single-start procedure (old_analyze_sample: ode15s from zero
%       initial conditions at ATOT = 0.1 -> one steady state -> one forward
%       continuation, plus continuation from any BP). Shows what it misses:
%       detached branches (never reached) and anything beyond the ATOT it reaches.
%
%   The samples are drawn with their own random stream (RandStream with 'seed'), so
%   the main LHS design and the choice of MATCONT-validated samples are unchanged.
%   Figures are saved as .png and .fig, in pages of up to 25 panels. A summary table
%   (one row per sample) is written as MATCONT_check_summary.csv and returned in C.

if nCheck <= 0, C = []; return; end
if ~exist(figDir, 'dir'), mkdir(figDir); end
CI = class_info_v3();

ok = find(arrayfun(@(r) r.exact_ok, res));
rs = RandStream('mt19937ar', 'Seed', seed);
pick = ok(randperm(rs, numel(ok), min(nCheck, numel(ok))));
pick = sort(pick(:));
n = numel(pick);
fprintf('\n=== MATCONT check on %d random samples (seed %d) ===\n', n, seed);

cfgC = cfg; cfgC.keep_mc_curves = true;
Rnew = cell(n,1); Rold = cell(n,1); P = paramSamples(pick,:);
npts = cfg.curve_points;
if useParallel
    parfor i = 1:n
        Rnew{i} = analyze_sample(P(i,:), cfgC, true);
        Rold{i} = old_analyze_sample(P(i,:), oldOpt, npts);
    end
else
    for i = 1:n
        Rnew{i} = analyze_sample(P(i,:), cfgC, true);
        Rold{i} = old_analyze_sample(P(i,:), oldOpt, npts);
        fprintf('  %d / %d\n', i, n);
    end
end
Rnew = [Rnew{:}]'; Rold = [Rold{:}]';

%% ---- per-sample summary -------------------------------------------------------------
W = cfg.ATOT_window;
match = @(ex, found) sum(arrayfun(@(a) any(abs(found - a) <= 1e-3*max(1, a)), ex));
rows = cell(n, 16);
for i = 1:n
    r = Rnew(i); o = Rold(i);
    exLP = r.exact_LP_ATOT; exBP = r.exact_BP_ATOT;
    oLP = o.LP_ATOT(o.LP_ATOT >= W(1) & o.LP_ATOT <= W(2));
    oBP = o.BP_ATOT(o.BP_ATOT >= W(1) & o.BP_ATOT <= W(2));
    nEx = numel(exLP) + numel(exBP);
    nNew = match(exLP, r.mc_LP_ATOT) + match(exBP, r.mc_BP_ATOT);
    nOld = match(exLP, oLP) + match(exBP, oBP);
    if r.mc_ok, mcLab = r.mc_label; else, mcLab = 'N/A'; end
    rows(i,:) = {pick(i), r.exact_class, r.exact_label, r.lp_coincides_bp, numel(exLP), numel(exBP), ...
        r.mc_class, mcLab, nNew, isequal(r.mc_class, r.exact_class), ...
        o.class7, CI.names_7{min(o.class7, 7)}, nOld, nEx, ...
        numel(r.exact_LP_ATOT(strcmp(r.exact_LP_type, 'detached'))), o.ATOT_reached};
    if o.class7 == 8, rows{i,12} = 'FAILED_or_NEGATIVE'; end
end
C = cell2table(rows, 'VariableNames', {'Sample','Exact_Class_v3','Exact_Label_v3','LP_Coincides_BP', ...
    'Exact_n_LP','Exact_n_BP','MC_Seeded_Class_v3','MC_Seeded_Label_v3','MC_Seeded_singularities_found', ...
    'MC_Seeded_agrees','SingleStart_Class7','SingleStart_Label7','SingleStart_singularities_found', ...
    'Exact_singularities','Exact_n_LP_detached','SingleStart_ATOT_reached'});
writetable(C, fullfile(figDir, 'MATCONT_check_summary.csv'));
save(fullfile(figDir, 'MATCONT_check_results.mat'), 'C', 'Rnew', 'Rold', 'pick', 'seed', 'cfgC', 'oldOpt');
nTot = sum(C.Exact_singularities);
fprintf('Singularities found: seeded on every branch %d / %d, single start %d / %d\n', ...
    sum(C.MC_Seeded_singularities_found), nTot, sum(C.SingleStart_singularities_found), nTot);
fprintf('Class agreement (seeded MATCONT vs exact): %d / %d\n', sum(C.MC_Seeded_agrees), n);

%% ---- figures ----------------------------------------------------------------------------
perPage = 25; nPages = ceil(n/perPage);
for pg = 1:nPages
    ii = (pg-1)*perPage + 1 : min(n, pg*perPage);
    make_page(Rnew(ii), [], pick(ii), CI, 'A', pg, nPages, ...
        sprintf(['MATCONT WITH the inverse method: seeded on every exact branch (%d of %d singularities found, ' ...
                 'class agreement %d/%d)'], sum(C.MC_Seeded_singularities_found), nTot, sum(C.MC_Seeded_agrees), n), ...
        fullfile(figDir, sprintf('MATCONT_check_WITH_inverse_method_p%02d', pg)));
    make_page(Rnew(ii), Rold(ii), pick(ii), CI, 'B', pg, nPages, ...
        sprintf(['MATCONT WITHOUT the inverse method: single start from ode15s at A_{TOT} = %g ' ...
                 '(%d of %d singularities found)'], cfg.ATOT_seed, sum(C.SingleStart_singularities_found), nTot), ...
        fullfile(figDir, sprintf('MATCONT_check_WITHOUT_inverse_method_p%02d', pg)));
end
end

%% =========================================================================================
function make_page(Rn, Ro, ks, CI, mode, pg, nPages, ttl, base)
m = numel(Rn);
nc = ceil(sqrt(m)); nr = ceil(m/nc);
f = figure('Color','w', 'Position', [30 30 min(1900, 380*nc) min(1050, 300*nr + 90)], ...
    'Name', sprintf('MATCONT check %s, page %d/%d', mode, pg, nPages), 'NumberTitle','off');
grey = [0.45 0.45 0.45];
for i = 1:m
    r = Rn(i);
    ax = subplot(nr, nc, i); hold(ax, 'on'); box(ax, 'on');
    xl = auto_xlim(r); yl = y_range(r.curves, xl);
    draw_branches(ax, r.curves, 2.6);                       % exact (colour)
    nFound = 0; nEx = numel(r.exact_LP_ATOT) + numel(r.exact_BP_ATOT);
    if mode == 'A'                                          % seeded MATCONT (grey, on top)
        for c = r.mc_curves(:)'
            plot(ax, c.ATOT, c.x2, '-', 'Color', grey, 'LineWidth', 1.0);
        end
        lpA = r.mc_LP_ATOT; lpX = r.mc_LP_x2; bpA = r.mc_BP_ATOT; bpX = r.mc_BP_x2;
    else                                                    % single-start MATCONT (grey)
        o = Ro(i);
        if ~isempty(o.main), plot(ax, o.main(1,:), o.main(2,:), '-', 'Color', grey, 'LineWidth', 1.0); end
        for e = 1:numel(o.ext)
            plot(ax, o.ext{e}(1,:), o.ext{e}(2,:), '-', 'Color', grey, 'LineWidth', 1.0);
        end
        lpA = o.LP_ATOT; lpX = o.LP_x2; bpA = o.BP_ATOT; bpX = o.BP_x2;
    end
    % exact singularities (filled) and the ones MATCONT reported (grey crosses)
    if ~isempty(r.exact_LP_ATOT)
        plot(ax, r.exact_LP_ATOT, r.exact_LP_x2, 'o', 'MarkerFaceColor','r', 'MarkerEdgeColor','k', 'MarkerSize', 7);
    end
    if ~isempty(r.exact_BP_ATOT)
        plot(ax, r.exact_BP_ATOT, r.exact_BP_x2, 's', 'MarkerFaceColor','g', 'MarkerEdgeColor','k', 'MarkerSize', 8);
    end
    if ~isempty(lpA), plot(ax, lpA, lpX, 'x', 'Color', [0.15 0.15 0.15], 'MarkerSize', 10, 'LineWidth', 1.6); end
    if ~isempty(bpA), plot(ax, bpA, bpX, '+', 'Color', [0.15 0.15 0.15], 'MarkerSize', 11, 'LineWidth', 1.6); end
    for a = r.exact_LP_ATOT, nFound = nFound + any(abs(lpA - a) <= 1e-3*max(1, a)); end
    for a = r.exact_BP_ATOT, nFound = nFound + any(abs(bpA - a) <= 1e-3*max(1, a)); end
    xlim(ax, xl); ylim(ax, yl);
    t1 = sprintf('#%d  %d %s', ks(i), r.exact_class, r.exact_label);
    if r.lp_coincides_bp, t1 = [t1 ' (LP=BP)']; end
    t2 = sprintf('MATCONT found %d/%d singularities', nFound, nEx);
    if mode == 'B' && Ro(i).class7 == 8, t2 = [t2 ' | run failed/negative']; end
    title(ax, {t1, t2}, 'Interpreter','none', 'FontSize', 8);
    if i > m - nc, xlabel(ax, 'A_{TOT}', 'FontSize', 8); end
    if mod(i-1, nc) == 0, ylabel(ax, 'x_2', 'FontSize', 8); end
    set(ax, 'FontSize', 7);
end
% legend in the first panel
ax1 = subplot(nr, nc, 1);
h = [plot(ax1, nan, nan, '-', 'Color', [0 0.35 0.75], 'LineWidth', 2.6), ...
     plot(ax1, nan, nan, '-', 'Color', [0.85 0.2 0.1], 'LineWidth', 2.6), ...
     plot(ax1, nan, nan, '-', 'Color', grey, 'LineWidth', 1.0), ...
     plot(ax1, nan, nan, 'o', 'MarkerFaceColor','r', 'MarkerEdgeColor','k'), ...
     plot(ax1, nan, nan, 'x', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.6)];
legend(ax1, h, {'exact, main branch','exact, detached branch','MATCONT','exact LP / BP','MATCONT LP / BP'}, ...
    'Location','best', 'FontSize', 6, 'Box','off');
try
    sgtitle(f, {ttl, 'colour = exact branches (inverse method; solid stable, dashed unstable), grey = MATCONT'}, ...
        'FontSize', 11, 'FontWeight', 'bold');
catch
end
drawnow;
saveas(f, [base '.png']);
savefig(f, [base '.fig']);
fprintf('Saved %s (.png, .fig)\n', base);
end
