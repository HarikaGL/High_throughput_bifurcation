% =========================================================================
%  PUBLICATION-READY ENSEMBLE BIFURCATION DIAGRAMS – ALL 7 CLASSES
%  (Extended branches included)
% =========================================================================
clear; close all; clc;

%% ---------- USER SETTINGS ----------
outputDir = 'ABBA_iter4';
matDir    = fullfile(outputDir, 'MAT_Files');
resultFile = fullfile(matDir, 'Complete_Bifurcation_Results.mat');

maxCurves    = 1;                  % max representative curves per class
outFigName   = 'Ensemble_Bifurcation_All_7_Classes_1_iter4';
outFigFormat = '-dpng';
outFigResolution = '-r300';

%% ---------- CLASS COLOUR MAP (1–7) ----------
classColors = [
    0.00  0.60  0.00;   % 1 - LP S-shaped
    0.00  0.30  0.80;   % 2 - LP Z-shaped
    0.80  0.00  0.00;   % 3 - BP Only
    1.00  0.50  0.00;   % 4 - LP+BP (LP>BP)
    0.70  0.00  0.70;   % 5 - LP+BP (LP<BP)
    0.00  0.70  0.70;   % 6 - Biphasic
    0.50  0.50  0.50;   % 7 - OTHER
];

% Class labels (subplot titles)
classLabels = {
    'LP S-shaped';
    'LP Z-shaped';
    'BP Only';
    'LP+BP (LP>BP)';
    'LP+BP (LP<BP)';
    'Biphasic';
    'OTHER'
};

% EXACT bifurcation_type strings as stored in the .mat file
classTypeStrings = {
    'LP_ONLY_S-shaped';
    'LP_ONLY_Z-shaped';
    'BP_ONLY';
    'LP_BP_LPgreaterBP';
    'LP_BP_LPlessBP';
    'NO_SING_biphasic';
    'NO_SING_other'                     % explicitly used
};

%% ---------- LOAD DATA ----------
if ~exist(resultFile, 'file')
    error('File not found: %s\nPlease check the output directory.', resultFile);
end
fprintf('Loading results from: %s\n', resultFile);
load(resultFile, 'all_results');

sampleDetails = all_results.sample_details;
numSamples    = all_results.num_samples;

% Flag invalid samples (failed or negative concentration)
valid_mask = true(numSamples, 1);
for k = 1:numSamples
    if sampleDetails(k).failed || sampleDetails(k).negative_in_continuation
        valid_mask(k) = false;
    end
end

% Build a cell array of bifurcation type strings
bifType = cell(numSamples, 1);
for k = 1:numSamples
    if ~valid_mask(k)
        bifType{k} = 'INVALID';
    else
        if isfield(sampleDetails(k), 'bifurcation_type') && ...
                ~isempty(sampleDetails(k).bifurcation_type)
            bifType{k} = sampleDetails(k).bifurcation_type;
        else
            bifType{k} = 'UNKNOWN';
        end
    end
end

%% ---------- DEFINE CLASS INDICES (all 7 using exact labels) ----------
classIndices = cell(7,1);
for c = 1:7
    classIndices{c} = find(strcmp(bifType, classTypeStrings{c}));
end

% Determine which classes actually have samples
hasSamples = false(7,1);
for c = 1:7
    hasSamples(c) = ~isempty(classIndices{c});
end
activeClasses = find(hasSamples);

fprintf('Classes with samples:\n');
for c = activeClasses'
    fprintf('  Class %d (%s): %d samples\n', c, classLabels{c}, length(classIndices{c}));
end

%% ---------- SELECT REPRESENTATIVE CURVES ----------
rng(42);   % reproducible selection
sel = cell(7,1);
for c = activeClasses'
    sel{c} = pickRepresentatives(classIndices{c}, maxCurves);
    fprintf('  Class %d selected %d curves\n', c, length(sel{c}));
end

%% ---------- CREATE FIGURE ----------
nPlots = length(activeClasses);
nCols = min(4, nPlots);
nRows = ceil(nPlots / nCols);

fig = figure('Color','w', ...
    'Position',[100 100 400*nCols 350*nRows], ...
    'Name','Ensemble Bifurcation – All 7 Classes', ...
    'NumberTitle','off');

xLabelStr = 'Total ATOT';
yLabelStr = 'Concentration (State 2)';

% Plot each active class
for idx = 1:nPlots
    c = activeClasses(idx);
    subplot(nRows, nCols, idx);
    hold on; box on; grid on;

    % **** MODIFIED: Use new function that also plots extended branches ****
    plotClassCurvesExtended(sampleDetails, sel{c}, classColors(c,:), 2);

    % Add LP, BP, peak markers
    addClassMarkers(sampleDetails, sel{c}, bifType);

    title(sprintf('%s (n=%d)', classLabels{c}, length(sel{c})), ...
        'FontSize',12, 'FontWeight','bold');
    xlabel(xLabelStr, 'FontSize',10);
    ylabel(yLabelStr, 'FontSize',10);
end

sgtitle('Ensemble Bifurcation Diagrams – All 7 Classes', ...
    'FontSize',16, 'FontWeight','bold');

drawnow;

%% ---------- SAVE FIGURE ----------
saveDir = fullfile(outputDir, 'Figures');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end
savePath = fullfile(saveDir, outFigName);
fprintf('Saving figure to: %s\n', savePath);
print(fig, savePath, outFigFormat, outFigResolution);
savefig(fig, [savePath '.fig']);
fprintf('Done.\n');

%% ============== LOCAL FUNCTIONS ==============
function idxOut = pickRepresentatives(idxIn, maxN)
    if isempty(idxIn)
        idxOut = [];
        return;
    end
    n = length(idxIn);
    if n <= maxN
        idxOut = idxIn;
    else
        p = randperm(n);
        idxOut = idxIn(p(1:maxN));
    end
end

% -------------------------------------------------------------------------
% NEW: plot both main and extended branches for each sample
% -------------------------------------------------------------------------
function plotClassCurvesExtended(sampleDetails, indices, color, lw)
    % For each sample index, draw the main continuation curve (solid) and
    % any extended branches (dashed with the same colour).
    for i = 1:length(indices)
        k = indices(i);
        sd = sampleDetails(k);

        % --- main branch ---
        if isfield(sd.continuation_data, 'x') && ~isempty(sd.continuation_data.x)
            xMain = sd.continuation_data.x;
            paramRow = size(xMain, 1);
            stateRow = 2;
            plot(xMain(paramRow,:), xMain(stateRow,:), ...
                'Color', color, 'LineWidth', lw, 'LineStyle','-');
        end

        % --- extended branches (from BP continuation) ---
        if isfield(sd, 'extended_branches') && ~isempty(sd.extended_branches)
            x_ext = sd.extended_branches.x;   % cell array of branches
            if ~isempty(x_ext)
                for b = 1:length(x_ext)
                    xB = x_ext{b};
                    if ~isempty(xB)
                        paramRow = size(xB, 1);
                        % Use dashed line to distinguish from the main branch
                        plot(xB(paramRow,:), xB(2,:), ...
                            'Color', color, 'LineWidth', lw*0.8, 'LineStyle','--');
                    end
                end
            end
        end
    end
end

% -------------------------------------------------------------------------
% (Unchanged) add markers for LP, BP, biphasic peak
% -------------------------------------------------------------------------
function addClassMarkers(sampleDetails, indices, bifType)
    for i = 1:length(indices)
        k = indices(i);
        sd = sampleDetails(k);
        if ~isfield(sd.continuation_data, 'x')
            continue;
        end

        % 1) Limit points (LP)
        if isfield(sd, 'ATOT_LP_values') && ~isempty(sd.ATOT_LP_values) ...
                && isfield(sd, 'LP_concentrations') && ~isempty(sd.LP_concentrations)
            plot(sd.ATOT_LP_values, sd.LP_concentrations, ...
                'ro', 'MarkerSize', 6, 'MarkerFaceColor','r', ...
                'MarkerEdgeColor','k');
        end

        % 2) Branch points (BP)
        if isfield(sd, 'ATOT_BP_values') && ~isempty(sd.ATOT_BP_values) ...
                && isfield(sd, 'BP_concentrations') && ~isempty(sd.BP_concentrations)
            plot(sd.ATOT_BP_values, sd.BP_concentrations, ...
                'gs', 'MarkerSize', 6, 'MarkerFaceColor','g', ...
                'MarkerEdgeColor','k');
        end

        % 3) Biphasic peak (only for class "NO_SING_biphasic")
        if strcmp(bifType{k}, 'NO_SING_biphasic')
            xData = sd.continuation_data.x;
            stateRow = 2;
            paramRow = size(xData, 1);
            yData = xData(stateRow, :);
            try
                [pks, locs] = findpeaks(yData);
                if ~isempty(pks)
                    [~, maxIdx] = max(pks);
                    peakATOT = xData(paramRow, locs(maxIdx));
                    peakConc = pks(maxIdx);
                    plot(peakATOT, peakConc, 'm^', ...
                        'MarkerSize', 8, 'MarkerFaceColor','m', ...
                        'MarkerEdgeColor','k');
                end
            catch
                % skip if findpeaks fails
            end
        end
    end
end