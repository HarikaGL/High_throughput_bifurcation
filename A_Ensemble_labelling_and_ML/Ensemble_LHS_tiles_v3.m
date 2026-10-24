% ======================= Ensemble_LHS_tiles_v3.m ====================================
% Representative figure of the Latin hypercube sampling ("gradient tile view"):
% rows = parameter sets, columns = the 14 sampled kinetic parameters, colour = where
% the value lies inside its sampling range on a log scale
%     scaled = (log10 p - log10 lower) / (log10 upper - log10 lower)   in [0, 1]
% which is exactly the LHS coordinate used by Ensemble_LHS_v3 (lhsdesign, maximin).
% The parameter sets are the REAL ones of the ensemble (Parameter_Samples.mat).
%
% Run from the Ensemble_v3 folder (paths set automatically). Reads
%   Ensemble_iter<i>_v3_1e6/MAT_Files/Parameter_Samples.mat
% Saves Figures_combined_v3/LHS_Tiles_v3.png/.pdf/.fig
% ===================================================================================
clear; close all; clc;
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
cd(HERE);

%% ---------------- USER SETTINGS ----------------------------------------------------
iterations = 1;            % iteration(s) the rows are taken from, e.g. 1 or [1 2 3 4]
nRows      = 15;           % parameter sets shown
pickMode   = 'first';      % 'first'  : the first sets of each iteration (Excel order)
                           % 'random' : random sets (seed below)
rngSeed    = 1;
showIterationTag = false;  % row labels "it2: 5" instead of 1..nRows
titleText  = 'Latin Hypercube Sampling (Gradient Tile View)';
panelLabel = '';           % e.g. '[II]'
headerText = 'Parameter Set Generation';   % '' = none
% tick labels of the 14 sampled parameters (TeX), same order as param_names(2:end)
paramLabels = {'K_{TOT}','P_{TOT}','k_{1}','k_{11}','k_{2}','k_{3}','k_{33}','k_{4}', ...
               'k_{5}','k_{55}','k_{6}','k_{7}','k_{77}','k_{8}'};
cmap       = parula(256);
FS.ticks = 16; FS.label = 20; FS.title = 20; FS.header = 20;
LW.axes  = 2.5;            % axes and colour-bar outline
gridLines = false;         % thin white lines between the tiles
outName  = 'LHS_Tiles_v3';

%% ---------------- data --------------------------------------------------------------
Z = []; tags = {};
nPer = ceil(nRows/numel(iterations));
for it = iterations(:)'
    f = fullfile(HERE, sprintf('Ensemble_iter%d_v3_1e6', it), 'MAT_Files', 'Parameter_Samples.mat');
    if ~exist(f, 'file'), error('Missing %s', f); end
    S = load(f, 'paramSamples', 'param_names', 'bounds');
    names = S.param_names(2:end);                       % ATOT is not sampled
    lo = cellfun(@(n) S.bounds.(n)(1), names); hi = cellfun(@(n) S.bounds.(n)(2), names);
    P = S.paramSamples(:, 2:end);
    sc = (log10(P) - log10(lo)) ./ (log10(hi) - log10(lo));
    switch lower(pickMode)
        case 'random'
            rs = RandStream('mt19937ar', 'Seed', rngSeed + it);
            k = sort(randperm(rs, size(P,1), nPer));
        otherwise
            k = 1:nPer;
    end
    Z = [Z; sc(k,:)]; %#ok<AGROW>
    tags = [tags, arrayfun(@(j) sprintf('it%d: %d', it, j), k, 'UniformOutput', false)]; %#ok<AGROW>
end
Z = Z(1:nRows, :); tags = tags(1:nRows);
Z = min(max(Z, 0), 1);

%% ---------------- figure ------------------------------------------------------------
fig = figure('Color','w', 'Units','pixels', 'Position',[60 60 1000 900], ...
    'Name','LHS tiles', 'NumberTitle','off');
ax = axes('Parent', fig, 'Position', [0.14 0.17 0.66 0.68]);
imagesc(ax, Z, [0 1]); colormap(ax, cmap);
axis(ax, 'tight'); set(ax, 'YDir', 'reverse');
if gridLines
    hold(ax, 'on');
    for j = 0.5:1:size(Z,2)+0.5, plot(ax, [j j], [0.5 size(Z,1)+0.5], 'w-', 'LineWidth', 0.8); end
    for i = 0.5:1:size(Z,1)+0.5, plot(ax, [0.5 size(Z,2)+0.5], [i i], 'w-', 'LineWidth', 0.8); end
end
if showIterationTag, yl = tags; else, yl = arrayfun(@num2str, 1:nRows, 'UniformOutput', false); end
set(ax, 'XTick', 1:size(Z,2), 'XTickLabel', paramLabels, 'XTickLabelRotation', 45, ...
    'YTick', 1:nRows, 'YTickLabel', yl, 'TickLabelInterpreter', 'tex', ...
    'FontSize', FS.ticks, 'FontWeight', 'bold', 'LineWidth', LW.axes, 'TickLength', [0 0], 'Box', 'on');
xlabel(ax, 'Parameters', 'FontSize', FS.label, 'FontWeight', 'bold');
ylabel(ax, 'Samples',    'FontSize', FS.label, 'FontWeight', 'bold');
title(ax, titleText, 'FontSize', FS.title, 'FontWeight', 'bold');
cb = colorbar(ax); set(cb, 'Position', [0.84 0.17 0.035 0.68]);
set(cb, 'FontSize', FS.ticks, 'FontWeight', 'bold', 'LineWidth', LW.axes, 'Ticks', 0:0.2:1);
ylabel(cb, 'Scaled value in sampling range (log)', 'FontSize', FS.ticks, 'FontWeight', 'bold');
ax0 = axes('Parent', fig, 'Position', [0 0 1 1], 'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1]);
if ~isempty(headerText)
    text(ax0, 0.5, 0.955, headerText, 'FontSize', FS.header, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
end
if ~isempty(panelLabel)
    text(ax0, 0.02, 0.975, panelLabel, 'FontSize', FS.header, 'FontWeight', 'bold');
end

%% ---------------- save --------------------------------------------------------------
outDir = fullfile(HERE, 'Figures_combined_v3');
if ~exist(outDir, 'dir'), mkdir(outDir); end
base = fullfile(outDir, outName);
savefig(fig, [base '.fig']);
try
    exportgraphics(fig, [base '.png'], 'Resolution', 400);
    exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
catch
    print(fig, base, '-dpng', '-r400');
end
fprintf('Saved %s (.png, .pdf, .fig)\n', base);
