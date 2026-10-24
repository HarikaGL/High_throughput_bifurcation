% ================= COMPLETE BIFURCATION ANALYSIS CODE =================
% For the MATCONT code below, store all the data in the folder. 
% Store only in the Excel and .mat formats. 
clear 
close all
clc
%%
init
%%
% to freely vary, select type 1, to add constraints select type 2
type =1; 
tot_list = {}; % to break symmetry (use 'k5','k55' for SNTC, Leave it blank for PF)

%% ================= BLOCK 1: CREATE OUTPUT DIRECTORY =================
% outputDir = sprintf('Bifurcation_Analysis_TEST');
outputDir = sprintf('Ensemble_iter3');
mkdir(outputDir);
fprintf('Output directory created: %s\n', outputDir);

% Create subdirectories
figuresDir = fullfile(outputDir, 'Figures'); 
excelDir = fullfile(outputDir, 'Excel_Files');
matDir = fullfile(outputDir, 'MAT_Files');
mkdir(figuresDir);
mkdir(excelDir);
mkdir(matDir);
fprintf('Subdirectories created:\n');
fprintf('  - %s\n', figuresDir);
fprintf('  - %s\n', excelDir);
fprintf('  - %s\n', matDir);

%% ================= BLOCK 2: START TIMING =================
fprintf('\n=== Starting bifurcation analysis ===\n');
total_start_time = tic;

%% ================= BLOCK 3: BASE PARAMETERS SETUP =================
param_names = {'ATOT','KTOT','PTOT','k1','k11','k2','k3','k33','k4','k5','k55','k6','k7','k77','k8'};
Param_Base_Vec =[0.1;100;100;0.02;1;0.01;0.0024;1;0.21;0.198;1;0.01;0.0024;1;0.2]';
numParams = numel(Param_Base_Vec);

%% ================= BLOCK 4: LHS SETUP AND SAMPLING =================
%% ================= LHS SETUP WITH LOG-UNIFORM SAMPLING =================
lhs_start_time = tic;

numSamples = 10000;  % Ensemble size



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 vary_list = {'KTOT','PTOT','k1','k11','k2','k3','k33','k4','k5','k55',...
             'k6','k7','k77','k8'};

%tot_list = {}; % k6, k7, k77, k8 created only nonsingularities
%  bounds.KTOT = [10 500];
%  bounds.PTOT = [10 500];
% bounds.k1 = [1e-3 1e1];
% bounds.k11 = [1e-3 1e1];
% bounds.k2 = [1e-3 1e1];
% bounds.k3 = [1e-3 1e1];
% bounds.k33 = [1e-3 1e1];
% bounds.k4 = [1e-3 1e1];
% 
% bounds.k5 = [1e-3 1e1];
% bounds.k55 = [1e-3 1e1];
% bounds.k6 = [1e-3 1e1];
% bounds.k7 = [1e-3 1e1];
% bounds.k77 =[1e-3 1e1];
% bounds.k8 = [1e-3 1e1];

% % Bounds for Z dominant 
%  bounds.KTOT = [10 100];
%  bounds.PTOT = [100 500];
% bounds.k1 = [1e-1 1e1];
% bounds.k11 = [1e-3 1e-1];
% bounds.k2 = [1e-3 1e-1];
% bounds.k3 = [1e-1 1e1];
% bounds.k33 = [1e-3 1e-1];
% bounds.k4 = [1e-1 1e1];
% 
% bounds.k5 = [1e-3 1e-2];
% bounds.k55 = [1e-3 1e1];
% bounds.k6 = [1e-3 1e-1];
% bounds.k7 = [1e-3 1e1];
% bounds.k77 =[1e-3 1e1];
% bounds.k8 = [1e-1 1e1];

% Bounds for S Dominant
 bounds.KTOT = [100 500];
 bounds.PTOT = [10 100];
bounds.k1 = [1e-3 1e-2];
bounds.k11 = [1e-3 1e1];
bounds.k2 = [1e-3 1e-1];
bounds.k3 = [1e-3 1];
bounds.k33 = [1e-3 1e1];
bounds.k4 = [1e-1 1e1];

bounds.k5 = [1e-1 1e1];
bounds.k55 = [1e-3 1e-1];
bounds.k6 = [1e-3 1e-1];
bounds.k7 = [1e-1 1e1];
bounds.k77 =[1e-3 1e-1];
bounds.k8 = [1e-3 1e1];
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert bounds to log10 space
bounds_log = struct();
for p = 1:numel(vary_list)
    param_name = vary_list{p};
    bounds_log.(param_name) = log10(bounds.(param_name));
end

% Initialize parameter matrices
paramSamples = zeros(numSamples, numParams);
paramSamples_log = zeros(numSamples, numParams);

% Generate LHS design in log space
lhsNorm = lhsdesign(numSamples, numParams, ...
                   'criterion','maximin','iterations',50);
if type == 1 % No constraints 
% Sample each parameter
for p = 1:numParams
    param_name = param_names{p};

    if strcmp(param_name, 'ATOT')
        % ATOT is fixed at baseline (not sampled)
        paramSamples(:, p) = Param_Base_Vec(p);
        paramSamples_log(:, p) = log10(Param_Base_Vec(p));

    elseif ismember(param_name, vary_list)
        % Sample in log space
        log_lower = bounds_log.(param_name)(1);
        log_upper = bounds_log.(param_name)(2);
        log_samples = log_lower + lhsNorm(:, p) * (log_upper - log_lower);

        % Store log values
        paramSamples_log(:, p) = log_samples;

        % Convert to original space (guarantees positivity)
        paramSamples(:, p) = 10.^log_samples;

    else
        % Non-varying parameters use baseline values
        paramSamples(:, p) = Param_Base_Vec(p);
        paramSamples_log(:, p) = log10(Param_Base_Vec(p));
    end
end
elseif type ==2 % To apply symmetry
% Sample each parameter
for p = 1:numParams
    param_name = param_names{p};
    
    if strcmp(param_name, 'ATOT')
        % ATOT is fixed at baseline (not sampled)
        paramSamples(:, p) = Param_Base_Vec(p);
        paramSamples_log(:, p) = log10(Param_Base_Vec(p));
    elseif strcmp(param_name, 'KTOT')
        % Sample in log space
        log_lower = bounds_log.(param_name)(1);
        log_upper = bounds_log.(param_name)(2);
        log_samples = log_lower + lhsNorm(:, p) * (log_upper - log_lower);
        
        % Store log values
        paramSamples_log(:, p) = log_samples;
        
        % Convert to original space (guarantees positivity)
        paramSamples(:, p) = 10.^log_samples;    
        % PTOT
        % Store log values
        paramSamples_log(:, 3) = paramSamples_log(:, p);
        
        % Convert to original space (guarantees positivity)
        paramSamples(:, 3) = paramSamples(:, p);    

    elseif strcmp(param_name, 'k1')
        % Sample in log space
        log_lower = bounds_log.(param_name)(1);
        log_upper = bounds_log.(param_name)(2);
        log_samples = log_lower + lhsNorm(:, p) * (log_upper - log_lower);
        
        % Store log values
        paramSamples_log(:, p) = log_samples;
        
        % Convert to original space (guarantees positivity)
        paramSamples(:, p) = 10.^log_samples;    
        % k5
        % Store log values
        paramSamples_log(:, 10) = paramSamples_log(:, p);
        
        % Convert to original space (guarantees positivity)
        paramSamples(:, 10) = paramSamples(:, p);    
    
    elseif strcmp(param_name, 'k11')
        % Sample in log space
        log_lower = bounds_log.(param_name)(1);
        log_upper = bounds_log.(param_name)(2);
        log_samples = log_lower + lhsNorm(:, p) * (log_upper - log_lower);
        
        % Store log values
        paramSamples_log(:, p) = log_samples;
        
        % Convert to original space (guarantees positivity)
        paramSamples(:, p) = 10.^log_samples;    
        % k55
        % Store log values
        paramSamples_log(:, 11) = paramSamples_log(:, p);
        
        % Convert to original space (guarantees positivity)
        paramSamples(:, 11) = paramSamples(:, p);   
    
    elseif strcmp(param_name, 'k2')
        % Sample in log space
        log_lower = bounds_log.(param_name)(1);
        log_upper = bounds_log.(param_name)(2);
        log_samples = log_lower + lhsNorm(:, p) * (log_upper - log_lower);
        
        % Store log values
        paramSamples_log(:, p) = log_samples;
        
        % Convert to original space (guarantees positivity)
        paramSamples(:, p) = 10.^log_samples;    
        % k6
        % Store log values
        paramSamples_log(:, 12) = paramSamples_log(:, p);
        
        % Convert to original space (guarantees positivity)
        paramSamples(:, 12) = paramSamples(:, p);    
    
    elseif strcmp(param_name, 'k3')
        % Sample in log space
        log_lower = bounds_log.(param_name)(1);
        log_upper = bounds_log.(param_name)(2);
        log_samples = log_lower + lhsNorm(:, p) * (log_upper - log_lower);
        
        % Store log values
        paramSamples_log(:, p) = log_samples;
        
        % Convert to original space (guarantees positivity)
        paramSamples(:, p) = 10.^log_samples;    
        % k7
        % Store log values
        paramSamples_log(:, 13) = paramSamples_log(:, p);
        
        % Convert to original space (guarantees positivity)
        paramSamples(:, 13) = paramSamples(:, p);    

    elseif strcmp(param_name, 'k33')
        % Sample in log space
        log_lower = bounds_log.(param_name)(1);
        log_upper = bounds_log.(param_name)(2);
        log_samples = log_lower + lhsNorm(:, p) * (log_upper - log_lower);
        
        % Store log values
        paramSamples_log(:, p) = log_samples;
        
        % Convert to original space (guarantees positivity)
        paramSamples(:, p) = 10.^log_samples;    
        % k5
        % Store log values
        paramSamples_log(:, 14) = paramSamples_log(:, p);
        
        % Convert to original space (guarantees positivity)
        paramSamples(:, 14) = paramSamples(:, p);    

    elseif strcmp(param_name, 'k4')
        % Sample in log space
        log_lower = bounds_log.(param_name)(1);
        log_upper = bounds_log.(param_name)(2);
        log_samples = log_lower + lhsNorm(:, p) * (log_upper - log_lower);
        
        % Store log values
        paramSamples_log(:, p) = log_samples;
        
        % Convert to original space (guarantees positivity)
        paramSamples(:, p) = 10.^log_samples;    
        % k8
        % Store log values
        paramSamples_log(:, 15) = paramSamples_log(:, p);
        
        % Convert to original space (guarantees positivity)
        paramSamples(:, 15) = paramSamples(:, p);    

    elseif ismember(param_name, tot_list)
        % Sample in log space
        log_lower = bounds_log.(param_name)(1);
        log_upper = bounds_log.(param_name)(2);
        log_samples = log_lower + lhsNorm(:, p) * (log_upper - log_lower);

        % Store log values
        paramSamples_log(:, p) = log_samples;

        % Convert to original space (guarantees positivity)
        paramSamples(:, p) = 10.^log_samples;
        
    else
        % % Non-varying parameters use baseline values
        % paramSamples(:, p) = Param_Base_Vec(p);
        % paramSamples_log(:, p) = log10(Param_Base_Vec(p));
        fprintf('done\n');

    end
end
end
lhs_time = toc(lhs_start_time);

% Display information about the sampling
fprintf('Generated %d parameter sets using log-uniform LHS in %.2f seconds\n', numSamples, lhs_time);
fprintf('Parameter ranges:\n');
for p = 1:numParams
    param_name = param_names{p};
    if ismember(param_name, vary_list)
        fprintf('  %s: [%.4g, %.4g] (log10: [%.2f, %.2f])\n', ...
                param_name, bounds.(param_name)(1), bounds.(param_name)(2), ...
                bounds_log.(param_name)(1), bounds_log.(param_name)(2));
    elseif strcmp(param_name, 'ATOT')
        fprintf('  %s: fixed at %.4g\n', param_name, Param_Base_Vec(p));
    else
        fprintf('  %s: fixed at baseline value %.4g\n', param_name, Param_Base_Vec(p));
    end
end

% Create lb and ub vectors for reference
lb = zeros(1, numParams);
ub = zeros(1, numParams);
for p = 1:numParams
    param_name = param_names{p};
    if ismember(param_name, vary_list)
        lb(p) = bounds.(param_name)(1);
        ub(p) = bounds.(param_name)(2);
    else
        lb(p) = Param_Base_Vec(p);
        ub(p) = Param_Base_Vec(p);
    end
end
%% ================= BLOCK 5: INITIALIZE BASIC FLAGS =================
LP_flag   = false(numSamples,1);
BP_flag   = false(numSamples,1);
HB_flag   = false(numSamples,1);
NO_flag   = false(numSamples,1);
FAIL_flag = false(numSamples,1);

%% ================= BLOCK 5.5: INITIALIZE ALL CLASSIFICATION VARIABLES =================
% Initialize ALL classification variables early
LPBP_LPgreaterBP = false(numSamples, 1);
LPBP_LPlessBP = false(numSamples, 1);
NO_SING_biphasic = false(numSamples, 1);
NO_SING_other = false(numSamples, 1);
S_shaped_LP = false(numSamples, 1);
Z_shaped_LP = false(numSamples, 1);
LP_shape = cell(numSamples, 1);
LP_BP_comparison = cell(numSamples, 1);
Bifurcation_Type = cell(numSamples, 1);
Maxima_Count = zeros(numSamples, 1);
LP_count = zeros(numSamples, 1);
BP_count = zeros(numSamples, 1);

% ================= ADDED: BP COEFFICIENT STORAGE =================
BP_coefficient_a = cell(numSamples, 1);  % Store the 'a' coefficient for BP
BP_type = cell(numSamples, 1);           % Store BP type: 'Transcritical' or 'Pitchfork'
% ================================================================

% Initialize all cells with default values
for k = 1:numSamples
    LP_BP_comparison{k} = 'N/A';
    Bifurcation_Type{k} = 'Not Classified';
    LP_shape{k} = 'N/A';
    BP_coefficient_a{k} = [];  % Initialize as empty
    BP_type{k} = 'N/A';        % Initialize as N/A
end

%% ================= BLOCK 5.6: INITIALIZE NEGATIVE CONCENTRATION DETECTION =================
% Add flags for negative concentrations during continuation
NEGATIVE_in_continuation_flag = false(numSamples, 1);
NEGATIVE_concentration_points = cell(numSamples, 1); % Store points where negative concentrations occur
NEGATIVE_state_indices = cell(numSamples, 1); % Which states go negative
NEGATIVE_ATOT_values = cell(numSamples, 1); % ATOT values where negative occurs

% Initialize with default values
for k = 1:numSamples
    NEGATIVE_concentration_points{k} = [];
    NEGATIVE_state_indices{k} = [];
    NEGATIVE_ATOT_values{k} = [];
end

%% ================= BLOCK 6: ADD CATEGORY STORAGE =================
category = cell(numSamples,1);

%% ================= BLOCK 7: ADD STORAGE FOR BIFURCATION VALUES =================
ATOT_LP_values = cell(numSamples, 1);
ATOT_BP_values = cell(numSamples, 1);
LP_concentrations = cell(numSamples, 1);
BP_concentrations = cell(numSamples, 1);

%% ================= BLOCK 8: CREATE DATA STRUCTURES FOR SAVING =================
all_results = struct();
all_results.parameter_names = param_names;
all_results.parameter_bounds = bounds;
all_results.parameter_samples = paramSamples;
all_results.num_samples = numSamples;
all_results.sample_details = struct();

for k = 1:numSamples
    all_results.sample_details(k).sample_index = k;
    all_results.sample_details(k).parameters = paramSamples(k,:);
    all_results.sample_details(k).category = '';
    all_results.sample_details(k).hasLP = false;
    all_results.sample_details(k).hasBP = false;
    all_results.sample_details(k).hasHB = false;
    all_results.sample_details(k).failed = false;
    all_results.sample_details(k).continuation_data = [];
    all_results.sample_details(k).singularity_points = [];
    all_results.sample_details(k).extended_branches = [];
    all_results.sample_details(k).ATOT_LP_values = [];
    all_results.sample_details(k).ATOT_BP_values = [];
    all_results.sample_details(k).LP_concentrations = [];
    all_results.sample_details(k).BP_concentrations = [];
    all_results.sample_details(k).extended_LP_values = [];
    all_results.sample_details(k).extended_LP_concentrations = [];
    all_results.sample_details(k).negative_in_continuation = false;
    all_results.sample_details(k).negative_points = [];
    all_results.sample_details(k).negative_states = [];
    all_results.sample_details(k).negative_ATOT_values = [];
    % ================= ADDED: BP COEFFICIENTS =================
    all_results.sample_details(k).BP_coefficient_a = [];
    all_results.sample_details(k).BP_type = 'N/A';
    % =========================================================
end

%% ================= BLOCK 9: CREATE FIGURE 1 =================
fig1 = figure('Color','w','Position',[100 100 1200 800], 'Name', 'Bifurcation Categories', 'NumberTitle', 'off');

subplot(2,2,1); hold on; box on; title('LP only', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Parameter', 'FontSize', 10); ylabel('State 2', 'FontSize', 10);
subplot(2,2,2); hold on; box on; title('LP + BP', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Parameter', 'FontSize', 10); ylabel('State 2', 'FontSize', 10);
subplot(2,2,3); hold on; box on; title('BP only', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Parameter', 'FontSize', 10); ylabel('State 2', 'FontSize', 10);
subplot(2,2,4); hold on; box on; title('No singularities', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Parameter', 'FontSize', 10); ylabel('State 2', 'FontSize', 10);

current_category_index = 0;
fprintf('\n=== Starting main bifurcation analysis loop ===\n');
analysis_start_time = tic;

%% ================= BLOCK 10: MAIN LOOP - BIFURCATION ANALYSIS =================
for k = 1:numSamples
    fprintf('\n=== Sample %d / %d ===\n',k,numSamples)
    
    X0 = zeros(6,1);
    P0 = paramSamples(k,:);
    ap1 = 1;
    
    %% Integrate ODE
    hls = feval(@MAPKfull);
    opts = odeset('RelTol',1e-8);
    Pcell = num2cell(P0);
    
    try
        [t,y] = ode15s(hls{2},[0 5000],X0,opts,Pcell{:});
        IC = y(end,:)';
    catch
        FAIL_flag(k) = true;
        all_results.sample_details(k).failed = true;
        all_results.sample_details(k).category = 'FAILED';
        current_category_index = current_category_index + 1;
        category{current_category_index} = 'FAILED';
        ATOT_LP_values{k} = 'FAILED';
        ATOT_BP_values{k} = 'FAILED';
        LP_concentrations{k} = 'FAILED';
        BP_concentrations{k} = 'FAILED';
        continue
    end
    
    %% Initialize continuation
    try
        [x0,v0] = init_EP_EP(@MAPKfull,IC,P0,ap1);
    catch
        FAIL_flag(k) = true;
        all_results.sample_details(k).failed = true;
        all_results.sample_details(k).category = 'FAILED';
        current_category_index = current_category_index + 1;
        category{current_category_index} = 'FAILED';
        ATOT_LP_values{k} = 'FAILED';
        ATOT_BP_values{k} = 'FAILED';
        LP_concentrations{k} = 'FAILED';
        BP_concentrations{k} = 'FAILED';
        continue
    end
    
    %% Setup options and run continuation
    opt=contset;
    opt=contset(opt,'InitStepSize',1);
    opt=contset(opt,'MaxStepSize',1e-3);
    opt=contset(opt,'MinStepSize',1e-5);
    opt=contset(opt,'MaxNewtoniters',10);
    opt=contset(opt,'MaxCorriters',15);
    opt=contset(opt,'MaxTestiters',10);
    opt=contset(opt,'VarTolerance',1e-7);
    opt=contset(opt,'FunTolerance',1e-7);
    opt=contset(opt,'TestTolerance',1e-6);
    opt = contset(opt,'TestFunctions',1);
    opt = contset(opt,'Locators',[1 1 1]);
    opt=contset(opt,'Singularities',1);
    opt=contset(opt,'Adapt',5);
    opt=contset(opt,'MaxNumPoints',1500);
    opt=contset(opt,'CheckClosed',50);
    opt=contset(opt, 'Backward', 0);
    
    cont_failed = false;
    try
        [x,v,s,~,~] = cont(@equilibrium,x0,[],opt);
        
        % ===== CHECK FOR NEGATIVE CONCENTRATIONS DURING CONTINUATION =====
        nStates = length(IC);
        NV = size(x,1);
        param_row = NV; % ATOT parameter row
        
        % Check each state variable for negative values
        negative_detected = false;
        negative_points = [];
        negative_states = [];
        negative_ATOTs = [];
        
        for state_idx = 1:nStates
            state_values = x(state_idx, :);
            neg_indices = find(state_values < 0);
            
            if ~isempty(neg_indices)
                negative_detected = true;
                NEGATIVE_in_continuation_flag(k) = true;
                
                % Store negative points
                for idx = 1:length(neg_indices)
                    point_idx = neg_indices(idx);
                    negative_points = [negative_points; point_idx];
                    negative_states = [negative_states; state_idx];
                    negative_ATOTs = [negative_ATOTs; x(param_row, point_idx)];
                end
                
                fprintf('WARNING: Sample %d - State %d goes negative during continuation\n', k, state_idx);
                fprintf('  Number of negative points: %d\n', length(neg_indices));
                fprintf('  Range of negative values: [%.2e, %.2e]\n', ...
                    min(state_values(neg_indices)), max(state_values(neg_indices)));
                fprintf('  ATOT range where negative occurs: [%.4f, %.4f]\n', ...
                    min(x(param_row, neg_indices)), max(x(param_row, neg_indices)));
            end
        end
        
        if negative_detected
            % Store negative concentration information
            NEGATIVE_concentration_points{k} = negative_points;
            NEGATIVE_state_indices{k} = negative_states;
            NEGATIVE_ATOT_values{k} = negative_ATOTs;
            
            % Store in all_results
            all_results.sample_details(k).negative_in_continuation = true;
            all_results.sample_details(k).negative_points = negative_points;
            all_results.sample_details(k).negative_states = negative_states;
            all_results.sample_details(k).negative_ATOT_values = negative_ATOTs;
            
            % Create detailed report
            fprintf('=== Negative Concentration Report for Sample %d ===\n', k);
            fprintf('Total negative points: %d\n', length(negative_points));
            
            % Group by state
            unique_states = unique(negative_states);
            for state = unique_states'
                state_mask = (negative_states == state);
                state_neg_points = negative_points(state_mask);
                state_ATOTs = negative_ATOTs(state_mask);
                
                fprintf('  State %d: %d negative points\n', state, sum(state_mask));
                if ~isempty(state_neg_points)
                    state_vals = x(state, state_neg_points);
                    fprintf('    Min negative value: %.2e\n', min(state_vals));
                    fprintf('    Max negative value: %.2e\n', max(state_vals));
                    fprintf('    ATOT range: [%.4f, %.4f]\n', ...
                        min(state_ATOTs), max(state_ATOTs));
                end
            end
            fprintf('============================================\n');
        end
        
    catch
        cont_failed = true;
    end
    
    if cont_failed
        FAIL_flag(k) = true;
        all_results.sample_details(k).failed = true;
        all_results.sample_details(k).category = 'FAILED';
        category{k} = 'FAILED';
        ATOT_LP_values{k} = 'FAILED';
        ATOT_BP_values{k} = 'FAILED';
        LP_concentrations{k} = 'FAILED';
        BP_concentrations{k} = 'FAILED';
        continue
    end
    
    %% Store continuation data
    all_results.sample_details(k).continuation_data.x = x;
    all_results.sample_details(k).continuation_data.v = v;
    all_results.sample_details(k).continuation_data.s = s;
    
%% ================= ADDED: EXTRACT BP COEFFICIENTS =================
% Extract normal form coefficients for BP points
labels = {s.label};
idxBP = find(strcmp(labels, 'BP'));

if ~isempty(idxBP)
    fprintf('BP points detected: %d\n', length(idxBP));
    
    bp_coeffs = [];
    bp_types = {};
    extracted_count = 0;
    
    for j = 1:length(idxBP)
        bp_idx = idxBP(j);
        s_point = s(bp_idx);
        
        fprintf('  Examining BP point %d...\n', j);
        
        % Debug: Check the structure of s_point
        fprintf('    s_point fields: %s\n', strjoin(fieldnames(s_point), ', '));
        
        % In MATCONT, BP type information is usually in s.data
        % Let's check what's available
        if isfield(s_point, 'data')
            data_info = s_point.data;
            fprintf('    s_point.data type: %s\n', class(data_info));
            
            % Check if it's a struct with coefficient information
            if isstruct(data_info)
                fprintf('    s_point.data fields: %s\n', strjoin(fieldnames(data_info), ', '));
                
                % Try to find coefficient information
                found_coeff = false;
                
                % Check common field names in MATCONT
                if isfield(data_info, 'a')
                    a = data_info.a;
                    if isnumeric(a)
                        bp_coeffs = [bp_coeffs, a];
                        if abs(a) > 1e-6
                            bp_type = 'Transcritical';
                        else
                            bp_type = 'Pitchfork';
                        end
                        bp_types{end+1} = bp_type;
                        found_coeff = true;
                        fprintf('    Found coefficient a = %.4e -> %s\n', a, bp_type);
                    end
                elseif isfield(data_info, 'normalform')
                    if isfield(data_info.normalform, 'a')
                        a = data_info.normalform.a;
                        if isnumeric(a)
                            bp_coeffs = [bp_coeffs, a];
                            if abs(a) > 1e-6
                                bp_type = 'Transcritical';
                            else
                                bp_type = 'Pitchfork';
                            end
                            bp_types{end+1} = bp_type;
                            found_coeff = true;
                            fprintf('    Found coefficient a = %.4e in normalform -> %s\n', a, bp_type);
                        end
                    end
                end
                
                if ~found_coeff
                    % Try to print the entire data structure for debugging
                    fprintf('    Data structure contents:\n');
                    disp(data_info);
                    
                    bp_coeffs = [bp_coeffs, NaN];
                    bp_types{end+1} = 'Unknown (struct)';
                end
            elseif isnumeric(data_info) && ~isempty(data_info)
                % s.data might be a numeric array
                a = data_info(1);
                bp_coeffs = [bp_coeffs, a];
                if abs(a) > 1e-6
                    bp_type = 'Transcritical';
                else
                    bp_type = 'Pitchfork';
                end
                bp_types{end+1} = bp_type;
                fprintf('    Found numeric coefficient a = %.4e -> %s\n', a, bp_type);
            else
                fprintf('    s.data is not struct or numeric\n');
                bp_coeffs = [bp_coeffs, NaN];
                bp_types{end+1} = 'No coefficient data';
            end
        else
            fprintf('    No data field in s_point\n');
            bp_coeffs = [bp_coeffs, NaN];
            bp_types{end+1} = 'No data field';
        end
    end
    
    % Store coefficients and types
    if ~isempty(bp_coeffs)
        BP_coefficient_a{k} = bp_coeffs;
        BP_type{k} = bp_types;
        all_results.sample_details(k).BP_coefficient_a = bp_coeffs;
        all_results.sample_details(k).BP_type = bp_types;
        
        % Print summary
        fprintf('BP Classification Summary for sample %d:\n', k);
        for j = 1:length(bp_coeffs)
            if ~isnan(bp_coeffs(j))
                fprintf('  BP %d: Coefficient a = %.4e -> %s bifurcation\n', j, bp_coeffs(j), bp_types{j});
            else
                fprintf('  BP %d: No coefficient extracted, Type: %s\n', j, bp_types{j});
            end
        end
    end
else
    fprintf('  No BP points detected\n');
    BP_coefficient_a{k} = [];
    BP_type{k} = {'No BP'};
    all_results.sample_details(k).BP_coefficient_a = [];
    all_results.sample_details(k).BP_type = {'No BP'};
end
%% ================================================================
    %% Initial classification
    labels = {s.label};
    hasLP = any(strcmp(labels,'LP'));
    hasBP = any(strcmp(labels,'BP'));
    hasHB = any(strcmp(labels,'HB'));
    
    LP_flag(k) = hasLP;
    BP_flag(k) = hasBP;
    HB_flag(k) = hasHB;
    
    all_results.sample_details(k).hasLP = hasLP;
    all_results.sample_details(k).hasBP = hasBP;
    all_results.sample_details(k).hasHB = hasHB;
    
    if ~(hasLP || hasBP || hasHB)
        NO_flag(k) = true;
    end
    
    nStates = length(IC);
    NV = size(x,1);
    param_row = NV;
    state_row = 2;
    
    %% Extract bifurcation values from ORIGINAL branch
    lp_atot_vals = [];
    bp_atot_vals = [];
    lp_conc_vals = [];
    bp_conc_vals = [];
    
    if hasLP
        idxLP = find(strcmp(labels, 'LP'));
        for j = 1:length(idxLP)
            idx = s(idxLP(j)).index;
            atot_val = x(param_row, idx);
            conc_val = x(state_row, idx);
            lp_atot_vals = [lp_atot_vals, atot_val];
            lp_conc_vals = [lp_conc_vals, conc_val];
        end
    end
    
    if hasBP
        idxBP = find(strcmp(labels, 'BP'));
        for j = 1:length(idxBP)
            idx = s(idxBP(j)).index;
            atot_val = x(param_row, idx);
            conc_val = x(state_row, idx);
            bp_atot_vals = [bp_atot_vals, atot_val];
            bp_conc_vals = [bp_conc_vals, conc_val];
        end
    end
    
    %% ========== CRITICAL: BP EXTENSION TO FIND LP POINTS ==========
    x_ext = {};
    s_ext = {};
    branch_colors = {};
    extended_lp_atot_vals = [];
    extended_lp_conc_vals = [];
    
    if hasBP
        fprintf('BP detected — extending branches to find LP points.\n');
        bp_idx = find(strcmp(labels,'BP'));
        
        for ii = 1:length(bp_idx)
            kk = bp_idx(ii);
            fprintf('Extending from BP at index %d\n',kk);
            
            try
                xBP = x(1:nStates,s(kk).index);
                param_temp = P0;
                param_temp(ap1) = x(param_row,s(kk).index);
                
                [xx,vx] = init_BP_EP(@MAPKfull,xBP,param_temp,s(kk),1e-8);
                
                opt_ext = opt;
                opt_ext = contset(opt_ext, 'Backward', 1);
                opt_ext = contset(opt_ext, 'MaxNumPoints', 1500);
                [x2, v2, s2, h2, f2] = cont(@equilibrium, xx, vx, opt_ext);
                
                if ~isempty(x2)
                    x_ext{end+1} = x2;
                    s_ext{end+1} = s2;
                    branch_colors{end+1} = 'b';
                    
                    % Check for LP points in extended branch
                    if ~isempty(s2)
                        labels_ext = {s2.label};
                        idxLP_ext = find(strcmp(labels_ext, 'LP'));
                        for j = 1:length(idxLP_ext)
                            idx = s2(idxLP_ext(j)).index;
                            atot_val = x2(end, idx);
                            conc_val = x2(state_row, idx);
                            extended_lp_atot_vals = [extended_lp_atot_vals, atot_val];
                            extended_lp_conc_vals = [extended_lp_conc_vals, conc_val];
                            fprintf('Found LP on extended branch: ATOT=%.4f, Conc=%.4f\n', atot_val, conc_val);
                        end
                    end
                end
                
                opt_ext = contset(opt_ext, 'Backward', 0);
                 opt_ext = contset(opt_ext, 'MaxNumPoints', 1500);
                [x3, v3, s3, h3, f3] = cont(@equilibrium, xx, vx, opt_ext);
                
                if ~isempty(x3)
                    x_ext{end+1} = x3;
                    s_ext{end+1} = s3;
                    branch_colors{end+1} = 'r';
                    
                    % Check for LP points in extended branch
                    if ~isempty(s3)
                        labels_ext = {s3.label};
                        idxLP_ext = find(strcmp(labels_ext, 'LP'));
                        for j = 1:length(idxLP_ext)
                            idx = s3(idxLP_ext(j)).index;
                            atot_val = x3(end, idx);
                            conc_val = x3(state_row, idx);
                            extended_lp_atot_vals = [extended_lp_atot_vals, atot_val];
                            extended_lp_conc_vals = [extended_lp_conc_vals, conc_val];
                            fprintf('Found LP on extended branch: ATOT=%.4f, Conc=%.4f\n', atot_val, conc_val);
                        end
                    end
                end
                
            catch ME
                fprintf('Warning: BP extension failed: %s\n', ME.message);
                continue;
            end
        end
    end
    
    % Store extended branches data
    all_results.sample_details(k).extended_branches.x = x_ext;
    all_results.sample_details(k).extended_branches.s = s_ext;
    all_results.sample_details(k).extended_branches.colors = branch_colors;
    all_results.sample_details(k).extended_LP_values = extended_lp_atot_vals;
    all_results.sample_details(k).extended_LP_concentrations = extended_lp_conc_vals;
    
    % Combine LP values from original branch AND extended branches
    all_lp_atot_vals = [lp_atot_vals, extended_lp_atot_vals];
    all_lp_conc_vals = [lp_conc_vals, extended_lp_conc_vals];
    
    %% Store values (using combined LP values from original + extended)
    if NO_flag(k)
        if NEGATIVE_in_continuation_flag(k)
            ATOT_LP_values{k} = 'NO_SING_NEGATIVE';
            ATOT_BP_values{k} = 'NO_SING_NEGATIVE';
            LP_concentrations{k} = 'NO_SING_NEGATIVE';
            BP_concentrations{k} = 'NO_SING_NEGATIVE';
        else
            ATOT_LP_values{k} = 'NO_SING';
            ATOT_BP_values{k} = 'NO_SING';
            LP_concentrations{k} = 'NO_SING';
            BP_concentrations{k} = 'NO_SING';
        end
    elseif NEGATIVE_in_continuation_flag(k)
        ATOT_LP_values{k} = 'NEGATIVE_IN_CONT';
        ATOT_BP_values{k} = 'NEGATIVE_IN_CONT';
        LP_concentrations{k} = 'NEGATIVE_IN_CONT';
        BP_concentrations{k} = 'NEGATIVE_IN_CONT';
    elseif ~isempty(all_lp_atot_vals) || ~isempty(bp_atot_vals)
        if ~isempty(all_lp_atot_vals)
            ATOT_LP_values{k} = mat2str(all_lp_atot_vals, 6);
            LP_concentrations{k} = mat2str(all_lp_conc_vals, 6);
        else
            ATOT_LP_values{k} = 'NONE';
            LP_concentrations{k} = 'NONE';
        end
        
        if ~isempty(bp_atot_vals)
            ATOT_BP_values{k} = mat2str(bp_atot_vals, 6);
            BP_concentrations{k} = mat2str(bp_conc_vals, 6);
        else
            ATOT_BP_values{k} = 'NONE';
            BP_concentrations{k} = 'NONE';
        end
    else
        ATOT_LP_values{k} = 'NONE';
        ATOT_BP_values{k} = 'NONE';
        LP_concentrations{k} = 'NONE';
        BP_concentrations{k} = 'NONE';
    end
    
    all_results.sample_details(k).ATOT_LP_values = all_lp_atot_vals;
    all_results.sample_details(k).ATOT_BP_values = bp_atot_vals;
    all_results.sample_details(k).LP_concentrations = all_lp_conc_vals;
    all_results.sample_details(k).BP_concentrations = bp_conc_vals;
    
    %% Compare LP and BP values
    lp_bp_comparison = 'N/A';
    if ~isempty(all_lp_conc_vals) && ~isempty(bp_conc_vals)
        avg_lp = mean(all_lp_conc_vals);
        avg_bp = mean(bp_conc_vals);
        
        if avg_lp > avg_bp
            lp_bp_comparison = 'LP>BP';
        else
            lp_bp_comparison = 'LP<BP';
        end
    end
    all_results.sample_details(k).lp_bp_comparison = lp_bp_comparison;
    
    %% ========== CRITICAL FIX: CORRECT CLASSIFICATION LOGIC ==========
    % Store original LP detection (before considering extended branches)
    original_hasLP = hasLP;
    original_hasBP = hasBP;
    
    % Determine final category - FIXED LOGIC
    if NEGATIVE_in_continuation_flag(k)
        final_category = 'NEGATIVE_IN_CONTINUATION';
    elseif NO_flag(k)
        final_category = 'NO_SING';
    else
        % Check for extended LP points
        has_any_LP = original_hasLP || ~isempty(extended_lp_atot_vals);
        
        if original_hasBP && has_any_LP
            final_category = 'LP_BP';
            % Update hasLP flag to include extended LPs
            hasLP = true;
            LP_flag(k) = true; % Update LP_flag to include extended LP cases
        elseif original_hasBP && ~has_any_LP
            final_category = 'BP_ONLY';
        elseif has_any_LP && ~original_hasBP
            final_category = 'LP_ONLY';
            hasLP = true;
            LP_flag(k) = true;
        else
            final_category = 'OTHER';
        end
    end
    
    fprintf('Classification: %s (origLP=%d, origBP=%d, extLP=%d, anyLP=%d, Negative=%d)\n', ...
            final_category, original_hasLP, original_hasBP, ~isempty(extended_lp_atot_vals), ...
            (original_hasLP || ~isempty(extended_lp_atot_vals)), NEGATIVE_in_continuation_flag(k));
    
    % Store in category array
    category{k} = final_category;
    all_results.sample_details(k).category = final_category;
    % ================================================================
    
    %% ========== SKIP PLOTTING FOR NEGATIVE CONCENTRATION CASES ==========
    if NEGATIVE_in_continuation_flag(k)
        fprintf('Sample %d: Skipping plot due to negative concentrations in continuation\n', k);
        continue; % Skip the plotting section
    end
    
    %% Plotting (only for samples without negative concentrations)
    if NO_flag(k)
        subplot(2,2,4); hold on;
        plot(x(param_row,:), x(state_row,:), 'k-', 'LineWidth', 1)
        continue
    end
    
    if strcmp(final_category, 'LP_ONLY')
        plot_subplot = 1;
    elseif strcmp(final_category, 'BP_ONLY')
        plot_subplot = 3;
    elseif strcmp(final_category, 'LP_BP')
        plot_subplot = 2;
    else
        plot_subplot = 4;
    end
    
    subplot(2,2,plot_subplot); hold on;
    plot(x(param_row,:), x(state_row,:), 'k-', 'LineWidth', 1.5);
    
    % Mark LP points (original + extended)
    if ~isempty(all_lp_atot_vals)
        for j = 1:length(all_lp_atot_vals)
            plot(all_lp_atot_vals(j), all_lp_conc_vals(j), 'ro', ...
                 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        end
    end
    
    % Mark BP points
    if ~isempty(bp_atot_vals)
        for j = 1:length(bp_atot_vals)
            plot(bp_atot_vals(j), bp_conc_vals(j), 'gs', ...
                 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        end
    end
    
    % Plot extended branches and their LP points
    if hasBP && ~isempty(x_ext)
        for branch_idx = 1:length(x_ext)
            x_branch = x_ext{branch_idx};
            s_branch = s_ext{branch_idx};
            branch_color = branch_colors{branch_idx};
            
            if ~isempty(x_branch)
                param_row_branch = size(x_branch,1);
                
                if strcmp(branch_color, 'b')
                    plot(x_branch(param_row_branch,:), x_branch(state_row,:), 'b--', 'LineWidth', 1);
                else
                    plot(x_branch(param_row_branch,:), x_branch(state_row,:), 'r--', 'LineWidth', 1);
                end
                
                % Plot LP points on extended branches
                if ~isempty(s_branch)
                    labels_branch = {s_branch.label};
                    idxLP_branch = find(strcmp(labels_branch, 'LP'));
                    
                    for j = 1:length(idxLP_branch)
                        idx = s_branch(idxLP_branch(j)).index;
                        plot(x_branch(param_row_branch, idx), x_branch(state_row, idx), 'm*', ...
                             'MarkerSize', 12, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
                    end
                end
            end
        end
    end
end

analysis_time = toc(analysis_start_time);
fprintf('\n=== Main analysis completed in %.2f seconds ===\n', analysis_time);

%% ================= BLOCK 11: SAVE FIGURE 1 =================
fig1Path = fullfile(figuresDir, 'Bifurcation_Categories.png');
saveas(fig1, fig1Path);
fprintf('Saved bifurcation categories figure: %s\n', fig1Path);

%% ================= BLOCK 12: CALCULATE SUMMARY FLAGS =================
% Update summary flags to exclude negative continuation cases
% FIXED: Use the correct logic for classification
LP_only_flag    = false(numSamples, 1);
BP_only_flag    = false(numSamples, 1);
LPBP_flag       = false(numSamples, 1);
NO_only_flag    = false(numSamples, 1);

for k = 1:numSamples
    if FAIL_flag(k) || NEGATIVE_in_continuation_flag(k)
        continue;
    end
    
    cat = category{k};
    if strcmp(cat, 'LP_ONLY')
        LP_only_flag(k) = true;
    elseif strcmp(cat, 'BP_ONLY')
        BP_only_flag(k) = true;
    elseif strcmp(cat, 'LP_BP')
        LPBP_flag(k) = true;
    elseif strcmp(cat, 'NO_SING')
        NO_only_flag(k) = true;
    end
end

FAILED_flag = FAIL_flag;
NEGATIVE_CONT_flag = NEGATIVE_in_continuation_flag;

%% ================= BLOCK 13: PRINT SUMMARY =================
fprintf('\n========== SUMMARY ==========\n')
fprintf('Total samples                    : %d\n', numSamples)
fprintf('LP only (no negative)           : %d\n', sum(LP_only_flag))
fprintf('BP only (no negative)           : %d\n', sum(BP_only_flag))
fprintf('LP + BP (no negative)           : %d\n', sum(LPBP_flag))
fprintf('No singularity (no negative)    : %d\n', sum(NO_only_flag))
fprintf('Negative in continuation        : %d\n', sum(NEGATIVE_CONT_flag))
fprintf('Continuation failed             : %d\n', sum(FAILED_flag))

%% ================= BLOCK 14: PERFORM ADDITIONAL CLASSIFICATION =================
fprintf('\n=== Performing additional classification ===\n');

for k = 1:numSamples
    if FAIL_flag(k)
        Bifurcation_Type{k} = 'FAILED';
        LP_BP_comparison{k} = 'FAILED';
        LP_shape{k} = 'FAILED';
        continue;
    end
    
    % Check for negative concentrations in continuation
    if NEGATIVE_in_continuation_flag(k)
        Bifurcation_Type{k} = 'NEGATIVE_IN_CONTINUATION';
        LP_BP_comparison{k} = 'NEGATIVE_CONT';
        LP_shape{k} = 'NEGATIVE_CONT';
        
        % Store in all_results
        all_results.sample_details(k).lp_shape = 'NEGATIVE_CONT';
        all_results.sample_details(k).lp_bp_comparison = 'NEGATIVE_CONT';
        all_results.sample_details(k).bifurcation_type = 'NEGATIVE_IN_CONTINUATION';
        all_results.sample_details(k).maxima_count = 0;
        all_results.sample_details(k).lp_count = 0;
        all_results.sample_details(k).bp_count = 0;
        all_results.sample_details(k).is_biphasic = false;
        
        continue;
    end
    
    % Get the final category
    final_category = category{k};
    
    if strcmp(final_category, 'NO_SING')
        Bifurcation_Type{k} = 'NO_SINGULARITY';
        LP_BP_comparison{k} = 'NO_SINGULARITY';
        LP_shape{k} = 'NO_SINGULARITY';
        
        % Check for biphasic
        if isfield(all_results.sample_details(k).continuation_data, 'x')
            x_data = all_results.sample_details(k).continuation_data.x;
            if ~isempty(x_data)
                state_row = 2;
                state_vals = x_data(state_row, :);
                
                if length(state_vals) > 10
                    try
                        [max_vals, ~] = findpeaks(state_vals);
                        Maxima_Count(k) = length(max_vals);
                        
                        if length(max_vals) == 1
                            NO_SING_biphasic(k) = true;
                            Bifurcation_Type{k} = 'NO_SING_biphasic';
                            LP_BP_comparison{k} = 'Biphasic (1 max)';
                        else
                            NO_SING_other(k) = true;
                            Bifurcation_Type{k} = 'NO_SING_other';
                            LP_BP_comparison{k} = sprintf('Other (%d max)', length(max_vals));
                        end
                    catch
                        NO_SING_other(k) = true;
                        Bifurcation_Type{k} = 'NO_SING_other';
                        LP_BP_comparison{k} = 'Other (findpeaks failed)';
                    end
                else
                    NO_SING_other(k) = true;
                    Bifurcation_Type{k} = 'NO_SING_other';
                    LP_BP_comparison{k} = 'Other (insufficient points)';
                end
            end
        end
        
        % Store results
        all_results.sample_details(k).lp_shape = LP_shape{k};
        all_results.sample_details(k).lp_bp_comparison = LP_BP_comparison{k};
        all_results.sample_details(k).bifurcation_type = Bifurcation_Type{k};
        all_results.sample_details(k).maxima_count = Maxima_Count(k);
        all_results.sample_details(k).lp_count = LP_count(k);
        all_results.sample_details(k).bp_count = BP_count(k);
        
        if NO_SING_biphasic(k)
            all_results.sample_details(k).is_biphasic = true;
        elseif NO_SING_other(k)
            all_results.sample_details(k).is_biphasic = false;
        else
            all_results.sample_details(k).is_biphasic = false;
        end
        
        continue;
    end
    
    % For samples with singularities (NOT NO_SING)
    % Get ALL LP values (original + extended)
    all_lp_atot_vals = all_results.sample_details(k).ATOT_LP_values;
    all_lp_conc_vals = all_results.sample_details(k).LP_concentrations;
    bp_atot_vals = all_results.sample_details(k).ATOT_BP_values;
    bp_conc_vals = all_results.sample_details(k).BP_concentrations;
    
    % Count LP and BP points (including extended)
    if ~isempty(all_lp_conc_vals)
        LP_count(k) = length(all_lp_conc_vals);
    end
    if ~isempty(bp_conc_vals)
        BP_count(k) = length(bp_conc_vals);
    end
    
    fprintf('Sample %d: LP_count=%d, BP_count=%d, Category=%s\n', ...
            k, LP_count(k), BP_count(k), final_category);
    
    % ========== CLASSIFY BASED ON FINAL CATEGORY ==========
    % This is the key fix: Use the already determined final_category
    if strcmp(final_category, 'LP_ONLY')
        Bifurcation_Type{k} = 'LP_ONLY';
        
        if LP_count(k) >= 2
            if length(all_lp_conc_vals) >= 2 && length(all_lp_atot_vals) >= 2
                [atot_sorted, atot_indices] = sort(all_lp_atot_vals);
                conc_sorted = all_lp_conc_vals(atot_indices);
                
                if conc_sorted(1) > conc_sorted(2)
                    S_shaped_LP(k) = true;
                    LP_shape{k} = 'S-shaped';
                    LP_BP_comparison{k} = 'S-shaped';
                    Bifurcation_Type{k} = 'LP_ONLY_S-shaped';
                else
                    Z_shaped_LP(k) = true;
                    LP_shape{k} = 'Z-shaped';
                    LP_BP_comparison{k} = 'Z-shaped';
                    Bifurcation_Type{k} = 'LP_ONLY_Z-shaped';
                end
            else
                LP_shape{k} = 'Multiple LP (unclassified)';
                LP_BP_comparison{k} = 'Multiple LP';
            end
        else
            LP_shape{k} = 'Single LP';
            LP_BP_comparison{k} = 'Single LP';
        end
    
    elseif strcmp(final_category, 'BP_ONLY')
        Bifurcation_Type{k} = 'BP_ONLY';
        LP_BP_comparison{k} = 'BP_ONLY';
        LP_shape{k} = 'N/A';
    
    elseif strcmp(final_category, 'LP_BP')
        Bifurcation_Type{k} = 'LP_BP';
        
        if ~isempty(all_lp_conc_vals) && ~isempty(bp_conc_vals)
            avg_lp = mean(all_lp_conc_vals);
            avg_bp = mean(bp_conc_vals);
            
            if avg_lp > avg_bp
                LPBP_LPgreaterBP(k) = true;
                Bifurcation_Type{k} = 'LP_BP_LPgreaterBP';
                LP_BP_comparison{k} = 'LP>BP';
            else
                LPBP_LPlessBP(k) = true;
                Bifurcation_Type{k} = 'LP_BP_LPlessBP';
                LP_BP_comparison{k} = 'LP<BP';
            end
        else
            LP_BP_comparison{k} = 'LP+BP (no comparison)';
        end
        LP_shape{k} = 'N/A';
    else
        Bifurcation_Type{k} = final_category;
        LP_BP_comparison{k} = final_category;
        LP_shape{k} = 'N/A';
    end
    
    % Store results
    all_results.sample_details(k).lp_shape = LP_shape{k};
    all_results.sample_details(k).lp_bp_comparison = LP_BP_comparison{k};
    all_results.sample_details(k).bifurcation_type = Bifurcation_Type{k};
    all_results.sample_details(k).maxima_count = Maxima_Count(k);
    all_results.sample_details(k).lp_count = LP_count(k);
    all_results.sample_details(k).bp_count = BP_count(k);
    all_results.sample_details(k).is_biphasic = false;
end

fprintf('\n=== Additional Classification Summary ===\n');
fprintf('LP_ONLY: %d\n', sum(LP_only_flag));
fprintf('  -> S-shaped: %d\n', sum(S_shaped_LP));
fprintf('  -> Z-shaped: %d\n', sum(Z_shaped_LP));
fprintf('BP_ONLY: %d\n', sum(BP_only_flag));
fprintf('LP+BP: %d\n', sum(LPBP_flag));
fprintf('  -> LP > BP: %d\n', sum(LPBP_LPgreaterBP));
fprintf('  -> LP < BP: %d\n', sum(LPBP_LPlessBP));
fprintf('NO_SINGULARITY: %d\n', sum(NO_only_flag));
fprintf('  -> Biphasic: %d\n', sum(NO_SING_biphasic));
fprintf('  -> Other: %d\n', sum(NO_SING_other));
fprintf('NEGATIVE_IN_CONTINUATION: %d\n', sum(NEGATIVE_CONT_flag));

%% ================= BLOCK 15: 7-CLASS CLASSIFICATION =================
fprintf('\n=== Performing 7-class classification ===\n');

% ================= DEFINE 7 CLASSES =================
% 1. LP_ONLY_S-shaped (S-shaped LP)
% 2. LP_ONLY_Z-shaped (Z-shaped LP)
% 3. BP_ONLY (Only BP points)
% 4. LP_BP_LPgreaterBP (LP + BP with LP > BP)
% 5. LP_BP_LPlessBP (LP + BP with LP < BP)
% 6. NO_SING_biphasic (Biphasic without singularities)
% 7. OTHER (US - Unstable and other cases)

% Initialize class flags
Class1_S_shaped_LP = false(numSamples, 1);       % LP_ONLY_S-shaped
Class2_Z_shaped_LP = false(numSamples, 1);       % LP_ONLY_Z-shaped
Class3_BP_only = false(numSamples, 1);           % BP_ONLY
Class4_LPBP_LPgreaterBP = false(numSamples, 1);  % LP_BP_LPgreaterBP
Class5_LPBP_LPlessBP = false(numSamples, 1);     % LP_BP_LPlessBP
Class6_NO_SING_biphasic = false(numSamples, 1);  % NO_SING_biphasic
Class7_OTHER = false(numSamples, 1);             % OTHER

% Classify samples
for k = 1:numSamples
    if FAIL_flag(k)
        % Failed samples go to OTHER
        Class7_OTHER(k) = true;
        continue;
    end
    
    if NEGATIVE_in_continuation_flag(k)
        % Negative cases go to OTHER
        Class7_OTHER(k) = true;
        continue;
    end
    
    % Get the final category
    final_category = category{k};
    
    if strcmp(final_category, 'LP_ONLY')
        % Check if S-shaped or Z-shaped
        if S_shaped_LP(k)
            Class1_S_shaped_LP(k) = true;
        elseif Z_shaped_LP(k)
            Class2_Z_shaped_LP(k) = true;
        else
            % Other LP_ONLY cases go to OTHER
            Class7_OTHER(k) = true;
        end
        
    elseif strcmp(final_category, 'BP_ONLY')
        Class3_BP_only(k) = true;
        
    elseif strcmp(final_category, 'LP_BP')
        % Check LP vs BP comparison
        if LPBP_LPgreaterBP(k)
            Class4_LPBP_LPgreaterBP(k) = true;
        elseif LPBP_LPlessBP(k)
            Class5_LPBP_LPlessBP(k) = true;
        else
            % Other LP_BP cases go to OTHER
            Class7_OTHER(k) = true;
        end
        
    elseif strcmp(final_category, 'NO_SING')
        if NO_SING_biphasic(k)
            Class6_NO_SING_biphasic(k) = true;
        else
            % Other NO_SING cases (US) go to OTHER
            Class7_OTHER(k) = true;
        end
    else
        % All other categories go to OTHER
        Class7_OTHER(k) = true;
    end
end

% Calculate counts for 7 classes
class_counts_7 = [
    sum(Class1_S_shaped_LP);
    sum(Class2_Z_shaped_LP);
    sum(Class3_BP_only);
    sum(Class4_LPBP_LPgreaterBP);
    sum(Class5_LPBP_LPlessBP);
    sum(Class6_NO_SING_biphasic);
    sum(Class7_OTHER)
];

class_labels_7 = {
    'LP_ONLY_S-shaped';
    'LP_ONLY_Z-shaped';
    'BP_ONLY';
    'LP_BP_LPgreaterBP';
    'LP_BP_LPlessBP';
    'NO_SING_biphasic';
    'OTHER'
};

fprintf('\n=== 7-Class Classification Summary ===\n');
for i = 1:length(class_labels_7)
    fprintf('  %s: %d samples (%.1f%%)\n', ...
        class_labels_7{i}, class_counts_7(i), ...
        class_counts_7(i)/numSamples*100);
end

% Store for later use
all_results.classification_7.class_labels = class_labels_7;
all_results.classification_7.class_counts = class_counts_7;
all_results.classification_7.percentages = (class_counts_7/numSamples * 100)';
all_results.classification_7.proportions = (class_counts_7/numSamples)';

%% ================= BLOCK 16: CLASSIFICATION ANALYSIS =================
fprintf('\n=== Classification Analysis ===\n');
fprintf('Total samples for classification: %d\n', numSamples);
fprintf('Distribution of classes:\n');
fprintf('  LP_ONLY:        %d samples (%.1f%%)\n', sum(LP_only_flag), sum(LP_only_flag)/numSamples*100);
fprintf('  BP_ONLY:        %d samples (%.1f%%)\n', sum(BP_only_flag), sum(BP_only_flag)/numSamples*100);
fprintf('  LP_BP:          %d samples (%.1f%%)\n', sum(LPBP_flag), sum(LPBP_flag)/numSamples*100);
fprintf('  NO_SINGULARITY: %d samples (%.1f%%)\n', sum(NO_only_flag), sum(NO_only_flag)/numSamples*100);
fprintf('  NEGATIVE_CONT:  %d samples (%.1f%%)\n', sum(NEGATIVE_CONT_flag), sum(NEGATIVE_CONT_flag)/numSamples*100);
fprintf('  FAILED:         %d samples (%.1f%%)\n', sum(FAILED_flag), sum(FAILED_flag)/numSamples*100);

% Breakdown of LP+BP
lp_bp_indices = find(LPBP_flag);
if ~isempty(lp_bp_indices)
    lp_gt_bp = 0;
    lp_lt_bp = 0;
    for idx = lp_bp_indices'
        comparison = all_results.sample_details(idx).lp_bp_comparison;
        if strcmp(comparison, 'LP>BP')
            lp_gt_bp = lp_gt_bp + 1;
        elseif strcmp(comparison, 'LP<BP')
            lp_lt_bp = lp_lt_bp + 1;
        end
    end
    
    fprintf('    -> LP > BP:   %d samples (%.1f%%)\n', lp_gt_bp, lp_gt_bp/numSamples*100);
    fprintf('    -> LP < BP:   %d samples (%.1f%%)\n', lp_lt_bp, lp_lt_bp/numSamples*100);
end

% Create classification summary table
class_labels = {'LP_ONLY'; 'BP_ONLY'; 'LP_BP'; 'NO_SINGULARITY'; 'NEGATIVE_IN_CONTINUATION'; 'FAILED'};
class_counts = [sum(LP_only_flag); sum(BP_only_flag); sum(LPBP_flag); ...
                sum(NO_only_flag); sum(NEGATIVE_CONT_flag); sum(FAILED_flag)];
percentages = (class_counts/numSamples * 100);
proportions = class_counts/numSamples;

classification_summary = table(class_labels, class_counts, percentages, proportions, ...
                              'VariableNames', {'Class', 'Count', 'Percentage', 'Proportion'});

all_results.classification.class_labels = class_labels';
all_results.classification.class_counts = class_counts';
all_results.classification.percentages = percentages';
all_results.classification.proportions = proportions';

%% ================= BLOCK 17: BIFURCATION VALUE STATISTICS =================
fprintf('\n=== Bifurcation Value Statistics ===\n');

% Calculate number of LP and BP points
num_LP_points = zeros(numSamples, 1);
num_BP_points = zeros(numSamples, 1);

for k = 1:numSamples
    if NEGATIVE_in_continuation_flag(k)
        num_LP_points(k) = -1;  % Mark as negative case
        num_BP_points(k) = -1;  % Mark as negative case
    elseif ~strcmp(ATOT_LP_values{k}, 'FAILED') && ~strcmp(ATOT_LP_values{k}, 'NONE') && ~strcmp(ATOT_LP_values{k}, 'NO_SING')
        try
            lp_vals = str2num(ATOT_LP_values{k});
            num_LP_points(k) = length(lp_vals);
        catch
            num_LP_points(k) = 0;
        end
    end
    
    if NEGATIVE_in_continuation_flag(k)
        num_BP_points(k) = -1;  % Already marked
    elseif ~strcmp(ATOT_BP_values{k}, 'FAILED') && ~strcmp(ATOT_BP_values{k}, 'NONE') && ~strcmp(ATOT_BP_values{k}, 'NO_SING')
        try
            bp_vals = str2num(ATOT_BP_values{k});
            num_BP_points(k) = length(bp_vals);
        catch
            num_BP_points(k) = 0;
        end
    end
end

% Count samples with LP/BP values
has_lp_values = ~strcmp(ATOT_LP_values, 'NONE') & ~strcmp(ATOT_LP_values, 'FAILED') & ~strcmp(ATOT_LP_values, 'NO_SING') & ~NEGATIVE_in_continuation_flag;
has_bp_values = ~strcmp(ATOT_BP_values, 'NONE') & ~strcmp(ATOT_BP_values, 'FAILED') & ~strcmp(ATOT_BP_values, 'NO_SING') & ~NEGATIVE_in_continuation_flag;

% Count samples with LP from extended branches
has_extended_LP_values = zeros(numSamples, 1);
for k = 1:numSamples
    if ~FAIL_flag(k) && ~NEGATIVE_in_continuation_flag(k) && ~isempty(all_results.sample_details(k).extended_LP_values)
        has_extended_LP_values(k) = 1;
    end
end

fprintf('ATOT values at singular points:\n');
fprintf('  Samples with LP values: %d (%.1f%%)\n', sum(has_lp_values), sum(has_lp_values)/numSamples*100);
fprintf('  Samples with BP values: %d (%.1f%%)\n', sum(has_bp_values), sum(has_bp_values)/numSamples*100);
fprintf('  Samples with LP from extended branches: %d (%.1f%%)\n', sum(has_extended_LP_values), sum(has_extended_LP_values)/numSamples*100);
fprintf('  Samples with negative concentrations: %d (%.1f%%)\n', sum(NEGATIVE_CONT_flag), sum(NEGATIVE_CONT_flag)/numSamples*100);

if sum(num_LP_points > 0) > 0
    valid_lp_counts = num_LP_points(num_LP_points > 0);
    fprintf('  Average number of LP points per sample: %.2f\n', mean(valid_lp_counts));
end
if sum(num_BP_points > 0) > 0
    valid_bp_counts = num_BP_points(num_BP_points > 0);
    fprintf('  Average number of BP points per sample: %.2f\n', mean(valid_bp_counts));
end

% Store extended LP statistics
all_results.summary.samples_with_extended_LP = sum(has_extended_LP_values);
all_results.summary.percentage_extended_LP = sum(has_extended_LP_values)/numSamples*100;
all_results.summary.samples_with_negative_concentrations = sum(NEGATIVE_CONT_flag);
all_results.summary.percentage_negative_concentrations = sum(NEGATIVE_CONT_flag)/numSamples*100;

%% ================= BLOCK 18: CREATE FIGURE 3 - BAR PLOT =================
fig3_bar = figure('Color','w','Position',[300 300 800 500], 'Name', 'Bifurcation Statistics', 'NumberTitle', 'off');
counts_original = [sum(LP_only_flag), sum(BP_only_flag), sum(LPBP_flag), sum(NO_only_flag), sum(NEGATIVE_CONT_flag), sum(FAILED_flag)];
bar(counts_original);
set(gca,'XTickLabel',{'LP','BP','LP+BP','No sing.','Negative','Failed'}, 'FontSize',11, 'FontWeight', 'bold');
ylabel('Count', 'FontSize', 12, 'FontWeight', 'bold');
title('Bifurcation Statistics', 'FontSize', 14, 'FontWeight', 'bold');
box on; grid on;

for i = 1:length(counts_original)
    percentage = (counts_original(i)/sum(counts_original))*100;
    text(i, counts_original(i) + max(counts_original)*0.02, ...
        sprintf('%.1f%%', percentage), 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

%% ================= BLOCK 19: SAVE FIGURE 3 =================
fig3Path = fullfile(figuresDir, 'Bifurcation_Statistics.png');
saveas(fig3_bar, fig3Path);
fprintf('Saved bifurcation statistics figure: %s\n', fig3Path);

%%%% ================= BLOCK 20: SAVE DATA IN EXCEL FORMAT =================
fprintf('\n=== Saving data in Excel format ===\n');

% 1. Create the main data table
excelData = table();
for i = 1:numParams
    excelData.(param_names{i}) = paramSamples(:,i);
end

% Add basic classification columns
excelData.LP_detected = LP_flag;
excelData.BP_detected = BP_flag;
excelData.NO_singularity = NO_flag;
excelData.FAILED = FAIL_flag;
excelData.Has_Negative_In_Continuation = NEGATIVE_in_continuation_flag;

% Add final category - IMPORTANT: Use the category determined in the main loop
excelData.Category = category(1:numSamples);

% Create classification label - Use the same logic as in main loop
classification_label = cell(numSamples, 1);
for k = 1:numSamples
    if FAIL_flag(k)
        classification_label{k} = 'FAILED';
    elseif NEGATIVE_in_continuation_flag(k)
        classification_label{k} = 'NEGATIVE_IN_CONTINUATION';
    else
        classification_label{k} = category{k};
    end
end
excelData.Classification_Label = classification_label;

% ================= ADD BIFURCATION PARAMETER VALUES =================
% Add ATOT values at LP and BP points
excelData.ATOT_LP_Values = ATOT_LP_values;
excelData.ATOT_BP_Values = ATOT_BP_values;
excelData.LP_Concentrations = LP_concentrations;
excelData.BP_Concentrations = BP_concentrations;

% ================= ADD LP vs BP COMPARISON COLUMN =================
% Create LP vs BP comparison column from all_results
lp_bp_comparison_cell = cell(numSamples, 1);
for k = 1:numSamples
    if FAIL_flag(k)
        lp_bp_comparison_cell{k} = 'FAILED';
    elseif NEGATIVE_in_continuation_flag(k)
        lp_bp_comparison_cell{k} = 'NEGATIVE_CONT';
    elseif strcmp(category{k}, 'NO_SING')
        lp_bp_comparison_cell{k} = 'NO_SING';
    elseif strcmp(category{k}, 'LP_BP')
        % Get the comparison from all_results
        if isfield(all_results.sample_details(k), 'lp_bp_comparison') && ...
           ~isempty(all_results.sample_details(k).lp_bp_comparison)
            comparison = all_results.sample_details(k).lp_bp_comparison;
            if strcmp(comparison, 'LP>BP')
                lp_bp_comparison_cell{k} = 'LP > BP';
            elseif strcmp(comparison, 'LP<BP')
                lp_bp_comparison_cell{k} = 'LP < BP';
            else
                lp_bp_comparison_cell{k} = comparison;
            end
        else
            lp_bp_comparison_cell{k} = 'LP+BP';
        end
    else
        lp_bp_comparison_cell{k} = category{k};
    end
end
excelData.LP_vs_BP_Comparison = lp_bp_comparison_cell;

% ================= ADD BP COEFFICIENTS AND TYPE TO EXCEL =================
% Add BP coefficient 'a' and BP type to Excel table
BP_coefficient_str = cell(numSamples, 1);
BP_type_str = cell(numSamples, 1);

for k = 1:numSamples
    if FAIL_flag(k)
        BP_coefficient_str{k} = 'FAILED';
        BP_type_str{k} = 'FAILED';
    elseif NEGATIVE_in_continuation_flag(k)
        BP_coefficient_str{k} = 'NEGATIVE_CONT';
        BP_type_str{k} = 'NEGATIVE_CONT';
    elseif ~isempty(BP_coefficient_a{k}) && ~all(isnan(BP_coefficient_a{k}))
        % Format coefficient values
        coeff_vals = BP_coefficient_a{k};
        valid_coeffs = coeff_vals(~isnan(coeff_vals));
        
        if ~isempty(valid_coeffs)
            if length(valid_coeffs) == 1
                BP_coefficient_str{k} = sprintf('%.4e', valid_coeffs(1));
            else
                coeff_str = '[';
                for j = 1:length(valid_coeffs)
                    coeff_str = [coeff_str sprintf('%.4e', valid_coeffs(j))];
                    if j < length(valid_coeffs)
                        coeff_str = [coeff_str ', '];
                    end
                end
                coeff_str = [coeff_str ']'];
                BP_coefficient_str{k} = coeff_str;
            end
        else
            BP_coefficient_str{k} = 'No valid coefficients';
        end
        
        % Format BP types
        if ~isempty(BP_type{k}) && iscell(BP_type{k})
            % Filter out placeholder types
            valid_types = {};
            for t = 1:length(BP_type{k})
                if ~strcmp(BP_type{k}{t}, 'No data') && ~strcmp(BP_type{k}{t}, 'Extraction failed')
                    valid_types{end+1} = BP_type{k}{t};
                end
            end
            
            if ~isempty(valid_types)
                if length(valid_types) == 1
                    BP_type_str{k} = valid_types{1};
                else
                    type_str = strjoin(valid_types, '; ');
                    BP_type_str{k} = type_str;
                end
            else
                BP_type_str{k} = 'N/A';
            end
        else
            BP_type_str{k} = 'N/A';
        end
    else
        BP_coefficient_str{k} = 'N/A';
        BP_type_str{k} = 'N/A';
    end
end

excelData.BP_Coefficient_a = BP_coefficient_str;
excelData.BP_Type = BP_type_str;
% ======================================================================
% ======================================================================

% ================= CREATE DETAILED CATEGORY WITH ALL CLASSIFICATIONS =================
detailed_category = cell(numSamples, 1);
biphasic_classification = cell(numSamples, 1);

for k = 1:numSamples
    if FAIL_flag(k)
        detailed_category{k} = 'FAILED';
        biphasic_classification{k} = 'FAILED';
    elseif NEGATIVE_in_continuation_flag(k)
        detailed_category{k} = 'NEGATIVE_IN_CONTINUATION';
        biphasic_classification{k} = 'NEGATIVE_CONT';
    elseif strcmp(category{k}, 'NO_SING')
        % Check if biphasic
        if NO_SING_biphasic(k)
            detailed_category{k} = 'NO_SING_biphasic';
            biphasic_classification{k} = 'BIPHASIC';
        elseif NO_SING_other(k)
            detailed_category{k} = 'NO_SING_other';
            biphasic_classification{k} = 'NOT_BIPHASIC';
        else
            detailed_category{k} = 'NO_SINGULARITY';
            biphasic_classification{k} = 'UNKNOWN';
        end
    elseif strcmp(category{k}, 'LP_BP')
        % Check LP vs BP comparison
        if isfield(all_results.sample_details(k), 'lp_bp_comparison') && ...
           ~isempty(all_results.sample_details(k).lp_bp_comparison)
            comparison = all_results.sample_details(k).lp_bp_comparison;
            if strcmp(comparison, 'LP>BP')
                detailed_category{k} = 'LP_BP_LPgreaterBP';
            elseif strcmp(comparison, 'LP<BP')
                detailed_category{k} = 'LP_BP_LPlessBP';
            else
                detailed_category{k} = 'LP_BP';
            end
        else
            detailed_category{k} = 'LP_BP';
        end
        biphasic_classification{k} = 'N/A (has singularities)';
    elseif strcmp(category{k}, 'LP_ONLY')
        % Check if it's S-shaped or Z-shaped
        if isfield(all_results.sample_details(k), 'lp_shape') && ...
           ~isempty(all_results.sample_details(k).lp_shape)
            lp_shape = all_results.sample_details(k).lp_shape;
            if strcmp(lp_shape, 'S-shaped')
                detailed_category{k} = 'LP_ONLY_S-shaped';
            elseif strcmp(lp_shape, 'Z-shaped')
                detailed_category{k} = 'LP_ONLY_Z-shaped';
            else
                detailed_category{k} = 'LP_ONLY';
            end
        else
            detailed_category{k} = 'LP_ONLY';
        end
        biphasic_classification{k} = 'N/A (has singularities)';
    elseif strcmp(category{k}, 'BP_ONLY')
        detailed_category{k} = 'BP_ONLY';
        biphasic_classification{k} = 'N/A (has singularities)';
    else
        detailed_category{k} = category{k};
        biphasic_classification{k} = 'UNKNOWN';
    end
end

excelData.Detailed_Category = detailed_category;
excelData.Biphasic_Classification = biphasic_classification;

% ================= ADD 7-CLASS CLASSIFICATION =================
% Create 7-class multiclass label and class number
multiclass_label_7 = cell(numSamples, 1);
class_number_7 = zeros(numSamples, 1);

for k = 1:numSamples
    if FAIL_flag(k)
        multiclass_label_7{k} = 'FAILED';
        class_number_7(k) = 8;  % Use 8 for failed
    elseif NEGATIVE_in_continuation_flag(k)
        multiclass_label_7{k} = 'NEGATIVE';
        class_number_7(k) = 9;  % Use 9 for negative
    elseif Class1_S_shaped_LP(k)
        multiclass_label_7{k} = 'LP_ONLY_S-shaped';
        class_number_7(k) = 1;
    elseif Class2_Z_shaped_LP(k)
        multiclass_label_7{k} = 'LP_ONLY_Z-shaped';
        class_number_7(k) = 2;
    elseif Class3_BP_only(k)
        multiclass_label_7{k} = 'BP_ONLY';
        class_number_7(k) = 3;
    elseif Class4_LPBP_LPgreaterBP(k)
        multiclass_label_7{k} = 'LP_BP_LPgreaterBP';
        class_number_7(k) = 4;
    elseif Class5_LPBP_LPlessBP(k)
        multiclass_label_7{k} = 'LP_BP_LPlessBP';
        class_number_7(k) = 5;
    elseif Class6_NO_SING_biphasic(k)
        multiclass_label_7{k} = 'NO_SING_biphasic';
        class_number_7(k) = 6;
    elseif Class7_OTHER(k)
        multiclass_label_7{k} = 'OTHER';
        class_number_7(k) = 7;
    else
        multiclass_label_7{k} = 'UNKNOWN';
        class_number_7(k) = 0;
    end
end

% Add 7-class classification to Excel
excelData.Multiclass_Label_7 = multiclass_label_7;
excelData.Class_Number_7 = class_number_7;

% ================= ADD NEGATIVE CONTINUATION INFORMATION =================
% Add detailed negative continuation information
negative_continuation_info = cell(numSamples, 1);
negative_states_info = cell(numSamples, 1);
negative_ATOT_info = cell(numSamples, 1);

for k = 1:numSamples
    if NEGATIVE_in_continuation_flag(k)
        % Create detailed info string
        if ~isempty(NEGATIVE_concentration_points{k})
            num_points = length(NEGATIVE_concentration_points{k});
            unique_states = unique(NEGATIVE_state_indices{k});
            
            info_str = sprintf('Negative in %d points', num_points);
            states_str = sprintf('States: %s', mat2str(unique_states'));
            
            if ~isempty(NEGATIVE_ATOT_values{k})
                atot_range = sprintf('ATOT range: [%.4f, %.4f]', ...
                    min(NEGATIVE_ATOT_values{k}), max(NEGATIVE_ATOT_values{k}));
                info_str = sprintf('%s | %s | %s', info_str, states_str, atot_range);
            else
                info_str = sprintf('%s | %s', info_str, states_str);
            end
            
            negative_continuation_info{k} = info_str;
            negative_states_info{k} = mat2str(unique(unique_states));
            
            if ~isempty(NEGATIVE_ATOT_values{k})
                negative_ATOT_info{k} = sprintf('[%.4f, %.4f]', ...
                    min(NEGATIVE_ATOT_values{k}), max(NEGATIVE_ATOT_values{k}));
            else
                negative_ATOT_info{k} = 'N/A';
            end
        else
            negative_continuation_info{k} = 'Negative detected (no details)';
            negative_states_info{k} = 'N/A';
            negative_ATOT_info{k} = 'N/A';
        end
    else
        negative_continuation_info{k} = 'NONE';
        negative_states_info{k} = 'NONE';
        negative_ATOT_info{k} = 'NONE';
    end
end

excelData.Negative_Continuation_Info = negative_continuation_info;
excelData.Negative_States = negative_states_info;
excelData.Negative_ATOT_Range = negative_ATOT_info;

% Add extended branch information
has_extended_LP_cell = cell(numSamples, 1);
for k = 1:numSamples
    if FAIL_flag(k)
        has_extended_LP_cell{k} = 'FAILED';
    elseif NEGATIVE_in_continuation_flag(k)
        has_extended_LP_cell{k} = 'NEGATIVE_CONT';
    elseif strcmp(category{k}, 'NO_SING')
        has_extended_LP_cell{k} = 'NO_SING';
    elseif ~isempty(all_results.sample_details(k).extended_LP_values)
        has_extended_LP_cell{k} = 'TRUE';
    else
        has_extended_LP_cell{k} = 'FALSE';
    end
end
excelData.Has_Extended_LP = has_extended_LP_cell;

% Add the base ATOT value used in the simulation
atot_idx = find(strcmp(param_names, 'ATOT'));
excelData.Base_ATOT = paramSamples(:, atot_idx);

% Add calculated number of points
excelData.Num_LP_Points = num_LP_points;
excelData.Num_BP_Points = num_BP_points;

% Add presence/absence flags for classification
excelData.Has_LP_Points = num_LP_points > 0;
excelData.Has_BP_Points = num_BP_points > 0;

% ================= ADD LP_SHAPE INFORMATION =================
% Add LP shape information
lp_shape_cell = cell(numSamples, 1);
for k = 1:numSamples
    if FAIL_flag(k)
        lp_shape_cell{k} = 'FAILED';
    elseif NEGATIVE_in_continuation_flag(k)
        lp_shape_cell{k} = 'NEGATIVE_CONT';
    elseif strcmp(category{k}, 'NO_SING')
        lp_shape_cell{k} = 'NO_SING';
    elseif isfield(all_results.sample_details(k), 'lp_shape') && ...
           ~isempty(all_results.sample_details(k).lp_shape)
        lp_shape_cell{k} = all_results.sample_details(k).lp_shape;
    else
        lp_shape_cell{k} = 'N/A';
    end
end
excelData.LP_Shape = lp_shape_cell;

% ================= ADD MAXIMA COUNT FOR NO_SING CASES =================
% Add maxima count for biphasic analysis
maxima_count_cell = cell(numSamples, 1);
for k = 1:numSamples
    if FAIL_flag(k)
        maxima_count_cell{k} = 'FAILED';
    elseif NEGATIVE_in_continuation_flag(k)
        maxima_count_cell{k} = 'NEGATIVE_CONT';
    elseif strcmp(category{k}, 'NO_SING')
        maxima_count_cell{k} = num2str(Maxima_Count(k));
    else
        maxima_count_cell{k} = 'N/A';
    end
end
excelData.Maxima_Count = maxima_count_cell;

% Save to main Excel file
excelFilePath = fullfile(excelDir, 'Bifurcation_Analysis_Results.xlsx');
writetable(excelData, excelFilePath);
fprintf('1. Saved main file: %s\n', excelFilePath);

% 2. NOW CREATE Classification_Datasets.xlsx SEPARATELY
classificationFilePath = fullfile(excelDir, 'Classification_Datasets.xlsx');
fprintf('2. Creating classification file: %s\n', classificationFilePath);

% Create a NEW workbook
writetable(excelData, classificationFilePath, 'Sheet', 'Complete_Dataset');

% Create Features_Only sheet
featureData = excelData(:, 1:numParams);
writetable(featureData, classificationFilePath, 'Sheet', 'Features_Only');

% Create Labels_Only sheet - IMPROVED VERSION
labelData = table((1:numSamples)', classification_label, detailed_category, ...
                  lp_bp_comparison_cell, biphasic_classification, lp_shape_cell, maxima_count_cell, ...
                  BP_coefficient_str, BP_type_str, multiclass_label_7, class_number_7, ...
                  'VariableNames', {'Sample_Index', 'Main_Category', 'Detailed_Category', ...
                                    'LP_vs_BP_Comparison', 'Biphasic_Classification', ...
                                    'LP_Shape', 'Maxima_Count', 'BP_Coefficient_a', 'BP_Type', ...
                                    'Multiclass_Label_7', 'Class_Number_7'});
writetable(labelData, classificationFilePath, 'Sheet', 'Labels_Only');

% Create Bifurcation_Values sheet with detailed info
bifurcationData = table((1:numSamples)', ATOT_LP_values, ATOT_BP_values, ...
                       LP_concentrations, BP_concentrations, ...
                       num_LP_points, num_BP_points, ...
                       has_extended_LP_cell, lp_bp_comparison_cell, ...
                       biphasic_classification, maxima_count_cell, ...
                       BP_coefficient_str, BP_type_str, multiclass_label_7, class_number_7, ...
                       negative_continuation_info, ...
                       'VariableNames', {'Sample_Index', 'ATOT_LP_Values', 'ATOT_BP_Values', ...
                                         'LP_Concentrations', 'BP_Concentrations', ...
                                         'Num_LP_Points', 'Num_BP_Points', ...
                                         'Has_Extended_LP', 'LP_vs_BP_Comparison', ...
                                         'Biphasic_Classification', 'Maxima_Count', ...
                                         'BP_Coefficient_a', 'BP_Type', ...
                                         'Multiclass_Label_7', 'Class_Number_7', ...
                                         'Negative_Continuation_Info'});
writetable(bifurcationData, classificationFilePath, 'Sheet', 'Bifurcation_Values');

% Verify the file was created
if exist(classificationFilePath, 'file')
    fprintf('✓ SUCCESS: Classification_Datasets.xlsx created at: %s\n', classificationFilePath);
    
    % List what's in the Excel file
    [~, sheetNames] = xlsfinfo(classificationFilePath);
    fprintf('Sheets in the file:\n');
    for i = 1:length(sheetNames)
        fprintf('  - %s\n', sheetNames{i});
    end
else
    fprintf('✗ ERROR: File not created. Trying alternative...\n');
    
    % Try saving directly to current directory
    altPath = 'Classification_Datasets_Direct.xlsx';
    writetable(excelData, altPath, 'Sheet', 'Complete_Dataset');
    if exist(altPath, 'file')
        fprintf('✓ Saved to current directory: %s\n', altPath);
        fprintf('Full path: %s\\%s\n', pwd, altPath);
    end
end

% 3. Save summary with LP>BP/LP<BP breakdown
fprintf('\n3. Saving summary statistics...\n');

% Count LP>BP and LP<BP
lp_gt_bp_count = sum(strcmp(lp_bp_comparison_cell, 'LP > BP'));
lp_lt_bp_count = sum(strcmp(lp_bp_comparison_cell, 'LP < BP'));

% Count biphasic cases
biphasic_count = sum(strcmp(biphasic_classification, 'BIPHASIC'));
not_biphasic_count = sum(strcmp(biphasic_classification, 'NOT_BIPHASIC'));

% Create enhanced summary table
summary_categories = {
    'LP_ONLY', 
    'BP_ONLY', 
    'LP+BP_LP>BP', 
    'LP+BP_LP<BP', 
    'NO_SING_BIPHASIC', 
    'NO_SING_OTHER', 
    'NEGATIVE_IN_CONTINUATION',
    'FAILED'};
summary_counts = [
    sum(LP_only_flag);
    sum(BP_only_flag);
    lp_gt_bp_count;
    lp_lt_bp_count;
    biphasic_count;
    not_biphasic_count;
    sum(NEGATIVE_CONT_flag);
    sum(FAILED_flag)];

summaryTable = table(summary_categories, summary_counts, ...
                     (summary_counts/numSamples * 100), ...
                     'VariableNames', {'Category', 'Count', 'Percentage'});
writetable(summaryTable, excelFilePath, 'Sheet', 'Summary');

% 4. Save extended branch statistics
extended_stats = table(sum(has_extended_LP_values), ...
                       sum(has_extended_LP_values)/numSamples*100, ...
                       'VariableNames', {'Samples_with_Extended_LP', 'Percentage'});
writetable(extended_stats, excelFilePath, 'Sheet', 'Extended_Branch_Stats');

% 5. Save biphasic statistics
biphasic_stats = table(biphasic_count, not_biphasic_count, ...
                       biphasic_count/numSamples*100, not_biphasic_count/numSamples*100, ...
                       'VariableNames', {'Biphasic_Count', 'Not_Biphasic_Count', ...
                                         'Biphasic_Percentage', 'Not_Biphasic_Percentage'});
writetable(biphasic_stats, excelFilePath, 'Sheet', 'Biphasic_Stats');

% 6. Save negative continuation statistics
negative_cont_stats = table(sum(NEGATIVE_CONT_flag), ...
                           sum(NEGATIVE_CONT_flag)/numSamples*100, ...
                           'VariableNames', {'Negative_Continuation_Cases', 'Percentage'});
writetable(negative_cont_stats, excelFilePath, 'Sheet', 'Negative_Continuation_Stats');

% 7. Save 7-class summary
summary_7class = table();
summary_7class.Class_Number = (1:7)';
summary_7class.Class_Label = class_labels_7;
summary_7class.Count = class_counts_7;
summary_7class.Percentage = (class_counts_7/numSamples * 100);

writetable(summary_7class, excelFilePath, 'Sheet', '7_Class_Summary');

fprintf('\n=== Excel Files Created ===\n');
fprintf('1. Main file: %s\\Bifurcation_Analysis_Results.xlsx\n', excelDir);
fprintf('2. Classification file: %s\\Classification_Datasets.xlsx\n', excelDir);
fprintf('3. Check folder: %s\n', excelDir);

% Print biphasic summary
fprintf('\n=== Biphasic Classification Summary ===\n');
fprintf('Biphasic curves: %d samples (%.1f%%)\n', biphasic_count, biphasic_count/numSamples*100);
fprintf('Not biphasic (other NO_SING): %d samples (%.1f%%)\n', not_biphasic_count, not_biphasic_count/numSamples*100);
fprintf('Total NO_SINGULARITY: %d samples\n', sum(NO_only_flag));

% Print negative continuation summary
fprintf('\n=== Negative Continuation Summary ===\n');
fprintf('Negative concentration cases: %d samples (%.1f%%)\n', sum(NEGATIVE_CONT_flag), sum(NEGATIVE_CONT_flag)/numSamples*100);

% Print 7-class summary
fprintf('\n=== 7-Class Classification Summary ===\n');
for i = 1:length(class_labels_7)
    fprintf('  %s: %d samples (%.1f%%)\n', class_labels_7{i}, class_counts_7(i), class_counts_7(i)/numSamples*100);
end

%% ================= BLOCK 21: CREATE SEPARATE FILE FOR NEGATIVE CONTINUATION CASES =================
% Create separate Excel file for negative continuation cases
negative_cont_indices = find(NEGATIVE_in_continuation_flag);
if ~isempty(negative_cont_indices)
    negativeContData = excelData(negative_cont_indices, :);
    
    % Save negative continuation cases separately
    negativeContExcelPath = fullfile(excelDir, 'Negative_Continuation_Cases.xlsx');
    writetable(negativeContData, negativeContExcelPath);
    fprintf('\n=== Saved negative continuation cases separately ===\n');
    fprintf('Negative continuation cases file: %s\n', negativeContExcelPath);
    fprintf('Number of negative continuation cases: %d\n', length(negative_cont_indices));
    
    % Create a summary sheet for negative continuation cases
    neg_cont_summary = table();
    neg_cont_summary.Metric = {
        'Total negative continuation cases';
        'Percentage of total samples';
        'Average negative points per case';
        };
    
    % Calculate statistics
    total_neg_cont = length(negative_cont_indices);
    percentage_neg_cont = (total_neg_cont / numSamples) * 100;
    
    % Calculate average negative points
    neg_points_counts = zeros(total_neg_cont, 1);
    for i = 1:total_neg_cont
        idx = negative_cont_indices(i);
        if ~isempty(NEGATIVE_concentration_points{idx})
            neg_points_counts(i) = length(NEGATIVE_concentration_points{idx});
        end
    end
    avg_neg_points = mean(neg_points_counts(neg_points_counts > 0));
    
    neg_cont_summary.Value = {
        total_neg_cont;
        percentage_neg_cont;
        avg_neg_points;
        };
    
    writetable(neg_cont_summary, negativeContExcelPath, 'Sheet', 'Negative_Cont_Summary');
    
    % Save detailed negative point information
    if ~isempty(negative_cont_indices)
        detailed_neg_info = table();
        detailed_neg_info.Sample_Index = negative_cont_indices;
        
        num_neg_points = zeros(total_neg_cont, 1);
        affected_states = cell(total_neg_cont, 1);
        atot_ranges = cell(total_neg_cont, 1);
        
        for i = 1:total_neg_cont
            idx = negative_cont_indices(i);
            if ~isempty(NEGATIVE_concentration_points{idx})
                num_neg_points(i) = length(NEGATIVE_concentration_points{idx});
                affected_states{i} = mat2str(unique(NEGATIVE_state_indices{idx}));
                
                if ~isempty(NEGATIVE_ATOT_values{idx})
                    atot_ranges{i} = sprintf('[%.4f, %.4f]', ...
                        min(NEGATIVE_ATOT_values{idx}), max(NEGATIVE_ATOT_values{idx}));
                else
                    atot_ranges{i} = 'N/A';
                end
            end
        end
        
        detailed_neg_info.Num_Negative_Points = num_neg_points;
        detailed_neg_info.Affected_States = affected_states;
        detailed_neg_info.ATOT_Range = atot_ranges;
        
        writetable(detailed_neg_info, negativeContExcelPath, 'Sheet', 'Detailed_Negative_Info');
    end
end

%% ================= BLOCK 22: ADD ENGINEERED FEATURES FOR ML =================
% Add normalized features (important for some ML algorithms)
for i = 1:numParams
    param_data = paramSamples(:, i);
    if std(param_data) > 0
        excelData.(sprintf('%s_norm', param_names{i})) = ...
            (param_data - mean(param_data)) / std(param_data);
    else
        excelData.(sprintf('%s_norm', param_names{i})) = zeros(size(param_data));
    end
end

% === REPLACE: Hardcoded ratios with comprehensive biologically meaningful ratios ===
% Remove hardcoded ratios:
% excelData.k4_times_k5 = paramSamples(:,9) .* paramSamples(:,10);
% excelData.k4_div_k5 = paramSamples(:,9) ./ (paramSamples(:,10) + eps);
% excelData.ratio_k4_k8 = paramSamples(:,9) ./ (paramSamples(:,15) + eps);
% excelData.ratio_k5_k6 = paramSamples(:,10) ./ (paramSamples(:,12) + eps);

% === ADD: Comprehensive biologically meaningful ratios ===
% 1. Total amounts ratios (VERY important for bistability)
excelData.KTOT_div_PTOT = paramSamples(:,2) ./ (paramSamples(:,3) + eps);
excelData.ATOT_div_KTOT = paramSamples(:,1) ./ (paramSamples(:,2) + eps);
excelData.ATOT_div_PTOT = paramSamples(:,1) ./ (paramSamples(:,3) + eps);

% Find indices for parameter names
k1_idx = find(strcmp(param_names, 'k1'));
k11_idx = find(strcmp(param_names, 'k11'));
k3_idx = find(strcmp(param_names, 'k3'));
k33_idx = find(strcmp(param_names, 'k33'));
k4_idx = find(strcmp(param_names, 'k4'));
k5_idx = find(strcmp(param_names, 'k5'));
k6_idx = find(strcmp(param_names, 'k6'));
k7_idx = find(strcmp(param_names, 'k7'));
k77_idx = find(strcmp(param_names, 'k77'));
k8_idx = find(strcmp(param_names, 'k8'));

% 2. Phosphorylation/dephosphorylation pairs
excelData.k1_div_k11 = paramSamples(:,k1_idx) ./ (paramSamples(:,k11_idx) + eps);
excelData.k3_div_k33 = paramSamples(:,k3_idx) ./ (paramSamples(:,k33_idx) + eps);
excelData.k7_div_k77 = paramSamples(:,k7_idx) ./ (paramSamples(:,k77_idx) + eps);

% 3. Interaction terms between key phosphorylation rates
excelData.k4_times_k5 = paramSamples(:,k4_idx) .* paramSamples(:,k5_idx);
excelData.k4_div_k5 = paramSamples(:,k4_idx) ./ (paramSamples(:,k5_idx) + eps);
excelData.k4_div_k8 = paramSamples(:,k4_idx) ./ (paramSamples(:,k8_idx) + eps);
excelData.k5_div_k6 = paramSamples(:,k5_idx) ./ (paramSamples(:,k6_idx) + eps);

% 4. Cascade progression ratios (signal flow through the cascade)
excelData.k1_div_k3 = paramSamples(:,k1_idx) ./ (paramSamples(:,k3_idx) + eps);
excelData.k3_div_k4 = paramSamples(:,k3_idx) ./ (paramSamples(:,k4_idx) + eps);
excelData.k4_div_k7 = paramSamples(:,k4_idx) ./ (paramSamples(:,k7_idx) + eps);

% 5. Sum of phosphorylation rates (total forward flux potential)
excelData.sum_phosphorylation = paramSamples(:,k1_idx) + paramSamples(:,k3_idx) + ...
                                paramSamples(:,k4_idx) + paramSamples(:,k7_idx);
excelData.sum_dephosphorylation = paramSamples(:,k11_idx) + paramSamples(:,k33_idx) + ...
                                  paramSamples(:,k5_idx) + paramSamples(:,k77_idx) + ...
                                  paramSamples(:,k8_idx);

% 6. Phosphorylation to dephosphorylation ratio (overall balance)
excelData.total_phosph_to_dephosph = excelData.sum_phosphorylation ./ ...
                                     (excelData.sum_dephosphorylation + eps);

fprintf('Added comprehensive engineered features:\n');
fprintf('  - 3 total amount ratios\n');
fprintf('  - 3 phosphorylation/dephosphorylation pairs\n');
fprintf('  - 4 interaction terms between key rates\n');
fprintf('  - 3 cascade progression ratios\n');
fprintf('  - 3 flux summation terms\n');

% Save feature importance analysis placeholder
feature_importance = table(param_names', zeros(numParams,1), ...
                          'VariableNames', {'Parameter', 'Importance_Score'});
writetable(feature_importance, excelFilePath, 'Sheet', 'Feature_Importance_Template');

% Update the Excel file with new features
writetable(excelData, excelFilePath);
fprintf('Updated Excel file with comprehensive engineered features.\n');

%% ================= BLOCK 23: EXTRACT ATOT VALUES FOR DIFFERENT CATEGORIES =================
fprintf('\n=== Extracting ATOT values for different categories ===\n');

% ========== CRITICAL: INITIALIZE ALL VARIABLES FIRST ==========
% Initialize with "N/A" strings for all samples
ML_ATOT_LP = cell(numSamples, 1);
ML_ATOT_BP = cell(numSamples, 1);
ML_ATOT_LPBP_LP = cell(numSamples, 1);
ML_ATOT_LPBP_BP = cell(numSamples, 1);
ML_ATOT_MAX_BIPHASIC = cell(numSamples, 1);
ML_ATOT_NO_SING = cell(numSamples, 1);
ML_ATOT_NO_BIPHASIC = cell(numSamples, 1);
ML_CONC_MAX_BIPHASIC = cell(numSamples, 1);

% Initialize all with "N/A"
for k = 1:numSamples
    ML_ATOT_LP{k} = 'N/A';
    ML_ATOT_BP{k} = 'N/A';
    ML_ATOT_LPBP_LP{k} = 'N/A';
    ML_ATOT_LPBP_BP{k} = 'N/A';
    ML_ATOT_MAX_BIPHASIC{k} = 'N/A';
    ML_ATOT_NO_SING{k} = 'N/A';
    ML_ATOT_NO_BIPHASIC{k} = 'N/A';
    ML_CONC_MAX_BIPHASIC{k} = 'N/A';
end

% Now extract ATOT values for each sample
for k = 1:numSamples
    if FAIL_flag(k) || NEGATIVE_in_continuation_flag(k)
        continue;
    end
    
    % Get category
    cat = category{k};
    
    % Get ATOT values from bifurcation analysis
    if isfield(all_results.sample_details(k), 'ATOT_LP_values')
        lp_atot_vals = all_results.sample_details(k).ATOT_LP_values;
    else
        lp_atot_vals = [];
    end
    
    if isfield(all_results.sample_details(k), 'ATOT_BP_values')
        bp_atot_vals = all_results.sample_details(k).ATOT_BP_values;
    else
        bp_atot_vals = [];
    end
    
    % Store values based on category
    if strcmp(cat, 'LP_ONLY')
        if ~isempty(lp_atot_vals) && isnumeric(lp_atot_vals)
            % For LP only, take the average of all LP ATOT values
            avg_value = mean(lp_atot_vals);
            ML_ATOT_LP{k} = sprintf('%.6f', avg_value);
        end
        
    elseif strcmp(cat, 'BP_ONLY')
        if ~isempty(bp_atot_vals) && isnumeric(bp_atot_vals)
            % For BP only, take the average of all BP ATOT values
            avg_value = mean(bp_atot_vals);
            ML_ATOT_BP{k} = sprintf('%.6f', avg_value);
        end
        
    elseif strcmp(cat, 'LP_BP')
        if ~isempty(lp_atot_vals) && isnumeric(lp_atot_vals)
            avg_value = mean(lp_atot_vals);
            ML_ATOT_LPBP_LP{k} = sprintf('%.6f', avg_value);
        end
        if ~isempty(bp_atot_vals) && isnumeric(bp_atot_vals)
            avg_value = mean(bp_atot_vals);
            ML_ATOT_LPBP_BP{k} = sprintf('%.6f', avg_value);
        end
        
    elseif strcmp(cat, 'NO_SING')
        % Check if biphasic
        if isfield(all_results.sample_details(k).continuation_data, 'x')
            x_data = all_results.sample_details(k).continuation_data.x;
            if ~isempty(x_data)
                state_row = 2;
                param_row = size(x_data, 1);
                state_vals = x_data(state_row, :);
                atot_vals = x_data(param_row, :);
                
                % Find maxima for biphasic analysis
                try
                    [max_vals, max_indices] = findpeaks(state_vals);
                    if length(max_vals) == 1
                        % This is biphasic - store the ATOT at the maximum
                        atot_at_max = atot_vals(max_indices(1));
                        conc_at_max = max_vals(1);
                        
                        ML_ATOT_MAX_BIPHASIC{k} = sprintf('%.6f', atot_at_max);
                        ML_CONC_MAX_BIPHASIC{k} = sprintf('%.6f', conc_at_max);
                        ML_ATOT_NO_BIPHASIC{k} = sprintf('%.6f', atot_at_max);
                    else
                        % Not biphasic - take the mean ATOT value
                        mean_atot = mean(atot_vals);
                        ML_ATOT_NO_SING{k} = sprintf('%.6f', mean_atot);
                        ML_ATOT_NO_BIPHASIC{k} = sprintf('%.6f', mean_atot);
                    end
                catch
                    % If findpeaks fails, take mean
                    mean_atot = mean(atot_vals);
                    ML_ATOT_NO_SING{k} = sprintf('%.6f', mean_atot);
                    ML_ATOT_NO_BIPHASIC{k} = sprintf('%.6f', mean_atot);
                end
            end
        end
    end
end

%% ================= BLOCK 24: UPDATE EXCEL FILES WITH ATOT VALUES =================

% Add ATOT values to the main excelData table
excelData.ATOT_LP = ML_ATOT_LP;
excelData.ATOT_BP = ML_ATOT_BP;
excelData.ATOT_LPBP_LP = ML_ATOT_LPBP_LP;
excelData.ATOT_LPBP_BP = ML_ATOT_LPBP_BP;
excelData.ATOT_MAX_BIPHASIC = ML_ATOT_MAX_BIPHASIC;
excelData.ATOT_NO_SING = ML_ATOT_NO_SING;
excelData.ATOT_NO_BIPHASIC = ML_ATOT_NO_BIPHASIC;
excelData.CONC_MAX_BIPHASIC = ML_CONC_MAX_BIPHASIC;

% Save updated main Excel file
writetable(excelData, excelFilePath);
fprintf('Updated main Excel file with ATOT values (as strings with N/A for missing)\n');

%% ================= BLOCK 25: CREATE ML-READY DATASET WITH 7 CLASSES =================
fprintf('\n=== Creating ML-ready dataset with 7 classes ===\n');

% Create ML-ready dataset
ML_Ready_FilePath = fullfile(excelDir, 'ML_Ready_Dataset.xlsx');

% ================= CREATE THE SIMPLIFIED ML_READY_DATA TABLE =================
% Create a table with ONLY the requested columns
ML_Ready_Data_Simple = table();

% Add parameter columns in the EXACT order specified
ML_Ready_Data_Simple.KTOT = paramSamples(:, strcmp(param_names, 'KTOT'));
ML_Ready_Data_Simple.PTOT = paramSamples(:, strcmp(param_names, 'PTOT'));
ML_Ready_Data_Simple.k1 = paramSamples(:, strcmp(param_names, 'k1'));
ML_Ready_Data_Simple.k11 = paramSamples(:, strcmp(param_names, 'k11'));
ML_Ready_Data_Simple.k2 = paramSamples(:, strcmp(param_names, 'k2'));
ML_Ready_Data_Simple.k3 = paramSamples(:, strcmp(param_names, 'k3'));
ML_Ready_Data_Simple.k33 = paramSamples(:, strcmp(param_names, 'k33'));
ML_Ready_Data_Simple.k4 = paramSamples(:, strcmp(param_names, 'k4'));
ML_Ready_Data_Simple.k5 = paramSamples(:, strcmp(param_names, 'k5'));
ML_Ready_Data_Simple.k55 = paramSamples(:, strcmp(param_names, 'k55'));
ML_Ready_Data_Simple.k6 = paramSamples(:, strcmp(param_names, 'k6'));
ML_Ready_Data_Simple.k7 = paramSamples(:, strcmp(param_names, 'k7'));
ML_Ready_Data_Simple.k77 = paramSamples(:, strcmp(param_names, 'k77'));
ML_Ready_Data_Simple.k8 = paramSamples(:, strcmp(param_names, 'k8'));

% Add concentration data
ML_Ready_Data_Simple.LP_Concentrations = LP_concentrations;
ML_Ready_Data_Simple.BP_Concentrations = BP_concentrations;

% ================= CRITICAL FIX: Handle SN/BP bifurcations =================
% For SN/BP bifurcations, ATOT_LP should contain both LP and BP concentrations
ATOT_LP_combined = cell(numSamples, 1);
for k = 1:numSamples
    if FAIL_flag(k) || NEGATIVE_in_continuation_flag(k)
        ATOT_LP_combined{k} = 'N/A';
    elseif strcmp(category{k}, 'LP_BP')
        % For LP_BP cases, combine LP and BP concentrations
        lp_conc_str = LP_concentrations{k};
        bp_conc_str = BP_concentrations{k};
        
        % Extract numeric values
        if ~strcmp(lp_conc_str, 'NONE') && ~strcmp(lp_conc_str, 'FAILED') && ~strcmp(lp_conc_str, 'NO_SING')
            lp_vals = str2num(lp_conc_str); %#ok<ST2NM>
        else
            lp_vals = [];
        end
        
        if ~strcmp(bp_conc_str, 'NONE') && ~strcmp(bp_conc_str, 'FAILED') && ~strcmp(bp_conc_str, 'NO_SING')
            bp_vals = str2num(bp_conc_str); %#ok<ST2NM>
        else
            bp_vals = [];
        end
        
        % Combine and format
        combined_vals = [lp_vals, bp_vals];
        if ~isempty(combined_vals)
            ATOT_LP_combined{k} = mat2str(combined_vals, 6);
        else
            ATOT_LP_combined{k} = 'NONE';
        end
    else
        % For other cases, use LP concentrations as is
        ATOT_LP_combined{k} = LP_concentrations{k};
    end
end

% Add ATOT columns with the combined values for LP_BP cases
ML_Ready_Data_Simple.ATOT_LP = ATOT_LP_combined;
ML_Ready_Data_Simple.ATOT_BP = ATOT_BP_values;
ML_Ready_Data_Simple.ATOT_LPBP_LP = ML_ATOT_LPBP_LP;
ML_Ready_Data_Simple.ATOT_LPBP_BP = ML_ATOT_LPBP_BP;
ML_Ready_Data_Simple.ATOT_MAX_BIPHASIC = ML_ATOT_MAX_BIPHASIC;
ML_Ready_Data_Simple.ATOT_NO_SING = ML_ATOT_NO_SING;
ML_Ready_Data_Simple.ATOT_NO_BIPHASIC = ML_ATOT_NO_BIPHASIC;
ML_Ready_Data_Simple.CONC_MAX_BIPHASIC = ML_CONC_MAX_BIPHASIC;

% Add 7-class classification labels
ML_Ready_Data_Simple.Multiclass_Label_7 = multiclass_label_7;
ML_Ready_Data_Simple.Class_Number_7 = class_number_7;

% ================= SAVE THE SIMPLIFIED VERSION =================
% Save to Excel
writetable(ML_Ready_Data_Simple, ML_Ready_FilePath, 'Sheet', 'ML_Ready_Data');

% Also create CSV for easy ML import
csvFilePath = fullfile(excelDir, 'ML_Ready_Dataset.csv');
writetable(ML_Ready_Data_Simple, csvFilePath);

fprintf('Created ML-ready files with 7 classes:\n');
fprintf('  - %s\n', ML_Ready_FilePath);
fprintf('  - %s\n', csvFilePath);

% Display summary of columns
fprintf('\n=== Columns in ML_Ready_Dataset ===\n');
column_names = ML_Ready_Data_Simple.Properties.VariableNames;
for i = 1:length(column_names)
    fprintf('%d. %s\n', i, column_names{i});
end
fprintf('\nTotal columns: %d\n', length(column_names));

% Display 7-class distribution
fprintf('\n=== 7-Class Distribution in ML_Ready_Dataset ===\n');
for i = 1:length(class_labels_7)
    count = sum(class_number_7 == i);
    fprintf('  %s: %d samples (%.1f%%)\n', ...
        class_labels_7{i}, count, count/numSamples*100);
end

% Display special handling for SN/BP bifurcations
fprintf('\n=== Special handling for SN/BP bifurcations ===\n');
lp_bp_indices = find(strcmp(category, 'LP_BP'));
if ~isempty(lp_bp_indices)
    fprintf('Found %d LP_BP cases. ATOT_LP contains:\n', length(lp_bp_indices));
    for i = 1:min(3, length(lp_bp_indices))
        idx = lp_bp_indices(i);
        fprintf('  Sample %d: %s\n', idx, ML_Ready_Data_Simple.ATOT_LP{idx});
    end
    if length(lp_bp_indices) > 3
        fprintf('  ... and %d more\n', length(lp_bp_indices)-3);
    end
end

%% ================= BLOCK 26: SAVE DATA IN .MAT FORMAT =================
fprintf('\n=== Saving data in .mat format ===\n');

% Save complete results with ATOT values
all_results.ATOT_values.LP = ATOT_LP_combined;  % Use the combined version
all_results.ATOT_values.BP = ATOT_BP_values;
all_results.ATOT_values.LPBP_LP = ML_ATOT_LPBP_LP;
all_results.ATOT_values.LPBP_BP = ML_ATOT_LPBP_BP;
all_results.ATOT_values.MAX_BIPHASIC = ML_ATOT_MAX_BIPHASIC;
all_results.ATOT_values.NO_SING = ML_ATOT_NO_SING;
all_results.ATOT_values.NO_BIPHASIC = ML_ATOT_NO_BIPHASIC;
all_results.ATOT_values.CONC_MAX_BIPHASIC = ML_CONC_MAX_BIPHASIC;
all_results.ML_dataset = ML_Ready_Data_Simple;  % Use the simplified version
all_results.ML_ready_data = ML_Ready_Data_Simple;  % Also save as ML_ready_data for compatibility

% Store 7-class classification results
all_results.classification_7.class_labels = class_labels_7;
all_results.classification_7.class_counts = class_counts_7;
all_results.classification_7.percentages = (class_counts_7/numSamples * 100)';
all_results.classification_7.proportions = (class_counts_7/numSamples)';

matFilePath = fullfile(matDir, 'Complete_Bifurcation_Results.mat');
save(matFilePath, 'all_results', '-v7.3');
fprintf('Saved complete results structure with ATOT values: %s\n', matFilePath);

% Save bifurcation values separately
save(fullfile(matDir, 'Bifurcation_Values.mat'), 'ATOT_LP_values', 'ATOT_BP_values', ...
     'LP_concentrations', 'BP_concentrations', 'num_LP_points', 'num_BP_points');

% Save individual components
save(fullfile(matDir, 'Parameter_Samples.mat'), 'paramSamples', 'param_names', 'bounds', 'lb', 'ub', 'numSamples');
save(fullfile(matDir, 'Bifurcation_Flags.mat'), 'LP_flag', 'BP_flag', 'HB_flag', 'NO_flag', 'FAIL_flag', ...
     'LP_only_flag', 'BP_only_flag', 'LPBP_flag', 'NO_only_flag', 'FAILED_flag', 'NEGATIVE_CONT_flag', 'category', 'numSamples');

% Now save it with the correct variable name
save(fullfile(matDir, 'Summary_Statistics.mat'), 'counts_original', 'classification_summary', 'analysis_time');

% Save negative continuation detection data
save(fullfile(matDir, 'Negative_Continuation_Detection.mat'), ...
    'NEGATIVE_in_continuation_flag', ...
    'NEGATIVE_concentration_points', ...
    'NEGATIVE_state_indices', ...
    'NEGATIVE_ATOT_values');

% Save compact version
compactResults = struct();
compactResults.parameter_samples = paramSamples;
compactResults.parameter_names = param_names;
compactResults.categories = category(1:numSamples);
compactResults.ATOT_LP_values = ATOT_LP_values;
compactResults.ATOT_BP_values = ATOT_BP_values;
compactResults.LP_concentrations = LP_concentrations;
compactResults.BP_concentrations = BP_concentrations;
compactResults.num_LP_points = num_LP_points;
compactResults.num_BP_points = num_BP_points;
compactResults.has_extended_LP = has_extended_LP_values;
compactResults.negative_continuation.has_negative = NEGATIVE_in_continuation_flag;
compactResults.negative_continuation.negative_points = NEGATIVE_concentration_points;
compactResults.negative_continuation.negative_states = NEGATIVE_state_indices;
compactResults.negative_continuation.negative_ATOTs = NEGATIVE_ATOT_values;
compactResults.summary = all_results.summary;
compactResults.classification = all_results.classification;
compactResults.flags.LP = LP_flag;
compactResults.flags.BP = BP_flag;
compactResults.flags.HB = HB_flag;
compactResults.flags.NO = NO_flag;
compactResults.flags.FAIL = FAIL_flag;
compactResults.flags.NEGATIVE_CONT = NEGATIVE_CONT_flag;
compactResults.ATOT_values = all_results.ATOT_values;
compactResults.ML_dataset = ML_Ready_Data_Simple;

% Save 7-class classification to compact results
compactResults.classification_7 = all_results.classification_7;

save(fullfile(matDir, 'Compact_Results.mat'), 'compactResults');

% Save separate MAT file for negative continuation cases only
if ~isempty(negative_cont_indices)
    negative_cont_cases = struct();
    negative_cont_cases.indices = negative_cont_indices;
    negative_cont_cases.parameters = paramSamples(negative_cont_indices, :);
    negative_cont_cases.parameter_names = param_names;
    negative_cont_cases.continuation_data = cell(length(negative_cont_indices), 1);
    
    % Store continuation data for negative cases
    for i = 1:length(negative_cont_indices)
        idx = negative_cont_indices(i);
        negative_cont_cases.continuation_data{i} = all_results.sample_details(idx).continuation_data;
    end
    
    negative_cont_cases.negative_info = struct();
    negative_cont_cases.negative_info.points = NEGATIVE_concentration_points(negative_cont_indices);
    negative_cont_cases.negative_info.states = NEGATIVE_state_indices(negative_cont_indices);
    negative_cont_cases.negative_info.ATOTs = NEGATIVE_ATOT_values(negative_cont_indices);
    
    save(fullfile(matDir, 'Negative_Continuation_Cases.mat'), 'negative_cont_cases', '-v7.3');
    fprintf('Saved negative continuation cases to: %s\\Negative_Continuation_Cases.mat\n', matDir);
end

fprintf('Saved all .mat files in directory: %s\n', matDir);

%% ================= BLOCK 27: CREATE FIGURE 4 - 7-CLASS DISTRIBUTION =================
%% ================= BLOCK 27: CREATE FIGURE 4 - 7-CLASS DISTRIBUTION =================
fig7class = figure('Color','w','Position',[300 300 1200 600], 'Name', '7-Class Distribution', 'NumberTitle', 'off');

% Calculate percentages
percentages = (class_counts_7/sum(class_counts_7))*100;

% Subplot 1: Regular bar chart
subplot(1,2,1);
bar(class_counts_7, 'FaceColor', [0.2, 0.6, 0.4]);
set(gca,'XTickLabel',class_labels_7, 'FontSize',9, 'FontWeight','bold', 'XTickLabelRotation',45);
ylabel('Count', 'FontSize',11, 'FontWeight','bold');
title('7-Class Distribution (Counts)', 'FontSize',13, 'FontWeight','bold');
grid on;

% Add count labels
for i = 1:length(class_counts_7)
    if class_counts_7(i) > 0
        text(i, class_counts_7(i) + max(class_counts_7)*0.02, ...
            sprintf('%d', class_counts_7(i)), ...
            'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize',9);
    end
end

% Subplot 2: PERCENTAGE bar chart (replacing pie chart)
subplot(1,2,2);
bar(percentages, 'FaceColor', [0.4, 0.2, 0.6]);
set(gca,'XTickLabel',class_labels_7, 'FontSize',9, 'FontWeight','bold', 'XTickLabelRotation',45);
ylabel('Percentage (%)', 'FontSize',11, 'FontWeight','bold');
title('7-Class Distribution (Percentages)', 'FontSize',13, 'FontWeight','bold');
grid on;

% Add percentage labels
for i = 1:length(percentages)
    if class_counts_7(i) > 0
        text(i, percentages(i) + max(percentages)*0.02, ...
            sprintf('%.1f%%', percentages(i)), ...
            'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize',9);
    end
end

% Save figure
fig7classPath = fullfile(figuresDir, '7_Class_Distribution.png');
saveas(fig7class, fig7classPath);
fprintf('Saved 7-class distribution figure: %s\n', fig7classPath);

%% ================= BLOCK 28: CREATE FIGURE 5 - ADDITIONAL CLASSIFICATION =================
fig5 = figure('Color','w','Position',[200 200 1400 800], 'Name', 'Figure 5: Additional Classification', 'NumberTitle', 'off');

% Create 6 subplots
for i = 1:6
    subplot(2,3,i); 
    hold on; box on; grid on;
    xlabel('ATOT', 'FontSize', 10);
    ylabel('State 2', 'FontSize', 10);
end

% Set initial titles
subplot(2,3,1); title('LP+BP: LP > BP', 'FontSize', 12, 'FontWeight', 'bold');
subplot(2,3,2); title('LP+BP: LP < BP', 'FontSize', 12, 'FontWeight', 'bold');
subplot(2,3,3); title('LP-only: S-shaped', 'FontSize', 12, 'FontWeight', 'bold');
subplot(2,3,4); title('LP-only: Z-shaped', 'FontSize', 12, 'FontWeight', 'bold');
subplot(2,3,5); title('No Singularity: Biphasic', 'FontSize', 12, 'FontWeight', 'bold');
subplot(2,3,6); title('Other Categories', 'FontSize', 12, 'FontWeight', 'bold');

% Track counts
counts = zeros(1,6);

% Plot all curves (skip negative cases)
for k = 1:numSamples
    if FAIL_flag(k) || NEGATIVE_in_continuation_flag(k)
        continue;
    end
    
    % Get data
    if ~isfield(all_results.sample_details(k).continuation_data, 'x')
        continue;
    end
    
    x_data = all_results.sample_details(k).continuation_data.x;
    if isempty(x_data)
        continue;
    end
    
    param_row = size(x_data, 1);
    state_row = 2;
    
    % Get the bifurcation type
    bif_type = Bifurcation_Type{k};
    
    % Determine which subplot
    subplot_idx = 0;
    
    if strcmp(bif_type, 'LP_BP_LPgreaterBP')
        subplot_idx = 1;  % LP > BP
    elseif strcmp(bif_type, 'LP_BP_LPlessBP')
        subplot_idx = 2;  % LP < BP
    elseif strcmp(bif_type, 'LP_ONLY_S-shaped')
        subplot_idx = 3;  % S-shaped
    elseif strcmp(bif_type, 'LP_ONLY_Z-shaped')
        subplot_idx = 4;  % Z-shaped
    elseif strcmp(bif_type, 'NO_SING_biphasic')
        subplot_idx = 5;  % Biphasic
    elseif strcmp(bif_type, 'LP_BP') || strcmp(bif_type, 'LP_ONLY') || strcmp(bif_type, 'BP_ONLY') || ...
           strcmp(bif_type, 'NO_SING_other') || strcmp(bif_type, 'NO_SINGULARITY')
        subplot_idx = 6;  % Other categories
    else
        subplot_idx = 6;  % Default to other
    end
    
    % Plot
    if subplot_idx > 0 && subplot_idx <= 6
        subplot(2,3,subplot_idx);
        
        % Choose color based on category
        if subplot_idx == 1
            color = [0, 0.4470, 0.7410]; % Blue
            ext_color1 = [0.3, 0.6, 0.9]; % Lighter blue for extended branches
            ext_color2 = [0.1, 0.3, 0.6]; % Darker blue for extended branches
        elseif subplot_idx == 2
            color = [0.8500, 0.3250, 0.0980]; % Orange
            ext_color1 = [0.95, 0.5, 0.2]; % Lighter orange
            ext_color2 = [0.7, 0.2, 0.05]; % Darker orange
        elseif subplot_idx == 3
            color = [0.4660, 0.6740, 0.1880]; % Green
            ext_color1 = [0.6, 0.8, 0.3]; % Lighter green
            ext_color2 = [0.3, 0.5, 0.1]; % Darker green
        elseif subplot_idx == 4
            color = [0.4940, 0.1840, 0.5560]; % Purple
            ext_color1 = [0.6, 0.3, 0.7]; % Lighter purple
            ext_color2 = [0.3, 0.1, 0.4]; % Darker purple
        elseif subplot_idx == 5
            color = [0.9290, 0.6940, 0.1250]; % Yellow
            ext_color1 = [0.95, 0.8, 0.3]; % Lighter yellow
            ext_color2 = [0.8, 0.5, 0.1]; % Darker yellow
        else
            color = [0.5, 0.5, 0.5]; % Gray
            ext_color1 = [0.7, 0.7, 0.7]; % Lighter gray
            ext_color2 = [0.3, 0.3, 0.3]; % Darker gray
        end
        
        % ========== PLOT MAIN BRANCH ==========
        plot(x_data(param_row,:), x_data(state_row,:), 'Color', color, 'LineWidth', 1.5);
        counts(subplot_idx) = counts(subplot_idx) + 1;
        
        % ========== PLOT EXTENDED BRANCHES ==========
        % Check if we have extended branches data
        if isfield(all_results.sample_details(k), 'extended_branches') && ...
           ~isempty(all_results.sample_details(k).extended_branches)
            
            x_ext = all_results.sample_details(k).extended_branches.x;
            branch_colors = all_results.sample_details(k).extended_branches.colors;
            
            if ~isempty(x_ext)
                for branch_idx = 1:length(x_ext)
                    x_branch = x_ext{branch_idx};
                    
                    if ~isempty(x_branch)
                        param_row_branch = size(x_branch, 1);
                        
                        % Plot extended branch with appropriate style
                        if branch_idx == 1
                            plot(x_branch(param_row_branch,:), x_branch(state_row,:), '--', ...
                                 'Color', ext_color1, 'LineWidth', 1);
                        else
                            plot(x_branch(param_row_branch,:), x_branch(state_row,:), ':', ...
                                 'Color', ext_color2, 'LineWidth', 1);
                        end
                    end
                end
            end
        end
        
        % ========== MARK LP POINTS ==========
        if LP_flag(k) && ~isempty(all_results.sample_details(k).LP_concentrations)
            lp_vals = all_results.sample_details(k).LP_concentrations;
            lp_atot = all_results.sample_details(k).ATOT_LP_values;
            if ~isempty(lp_vals)
                % Plot all LP points
                plot(lp_atot, lp_vals, 'ro', 'MarkerSize', 10, ...
                     'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
            end
        end
        
        % ========== MARK BP POINTS ==========
        if BP_flag(k) && ~isempty(all_results.sample_details(k).BP_concentrations)
            bp_vals = all_results.sample_details(k).BP_concentrations;
            bp_atot = all_results.sample_details(k).ATOT_BP_values;
            if ~isempty(bp_vals)
                % Plot all BP points
                plot(bp_atot, bp_vals, 'gs', 'MarkerSize', 10, ...
                     'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
            end
        end
        
        % For biphasic curves, mark the peak
        if strcmp(bif_type, 'NO_SING_biphasic') && subplot_idx == 5
            state_vals = x_data(state_row, :);
            try
                [max_vals, max_locs] = findpeaks(state_vals);
                if ~isempty(max_locs)
                    max_x = x_data(param_row, max_locs);
                    plot(max_x, max_vals, 'k*', 'MarkerSize', 12, 'LineWidth', 1.5);
                end
            catch
                % Skip if findpeaks fails
            end
        end
    end
end

% Update titles with counts
subplot(2,3,1); title(sprintf('LP+BP: LP > BP (n=%d)', counts(1)), 'FontSize', 12, 'FontWeight', 'bold');
subplot(2,3,2); title(sprintf('LP+BP: LP < BP (n=%d)', counts(2)), 'FontSize', 12, 'FontWeight', 'bold');
subplot(2,3,3); title(sprintf('LP-only: S-shaped (n=%d)', counts(3)), 'FontSize', 12, 'FontWeight', 'bold');
subplot(2,3,4); title(sprintf('LP-only: Z-shaped (n=%d)', counts(4)), 'FontSize', 12, 'FontWeight', 'bold');
subplot(2,3,5); title(sprintf('No Singularity: Biphasic (n=%d)', counts(5)), 'FontSize', 12, 'FontWeight', 'bold');
subplot(2,3,6); title(sprintf('Other Categories (n=%d)', counts(6)), 'FontSize', 12, 'FontWeight', 'bold');

% If no curves were plotted in a subplot, add a text message
for i = 1:6
    subplot(2,3,i);
    if counts(i) == 0
        text(0.5, 0.5, 'No curves in this category', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', 12, 'Color', 'red');
    end
end

fprintf('\n=== Figure 5 Summary ===\n');
fprintf('1. LP+BP LP>BP: %d\n', counts(1));
fprintf('2. LP+BP LP<BP: %d\n', counts(2));
fprintf('3. LP-only S-shaped: %d\n', counts(3));
fprintf('4. LP-only Z-shaped: %d\n', counts(4));
fprintf('5. Biphasic: %d\n', counts(5));
fprintf('6. Other: %d\n', counts(6));

%% ================= BLOCK 29: SAVE FIGURE 5 =================
fig5Path = fullfile(figuresDir, 'Figure5_Additional_Classification.png');
saveas(fig5, fig5Path);
fprintf('Saved Figure 5 (Additional Classification): %s\n', fig5Path);

%% ================= BLOCK 30: CREATE FIGURE 6 - SUMMARY STATISTICS =================
fig6 = figure('Color','w','Position',[300 300 1000 600], 'Name', 'Figure 6: Classification Summary', 'NumberTitle', 'off');

% Define main categories for Figure 6
main_categories = {
    'LP_ONLY', 
    'BP_ONLY', 
    'LP+BP', 
    'NO_SINGULARITY', 
    'NEGATIVE_CONT', 
    'FAILED'
};

% Calculate counts for main categories
main_counts = [
    sum(LP_only_flag);
    sum(BP_only_flag);
    sum(LPBP_flag);
    sum(NO_only_flag);
    sum(NEGATIVE_CONT_flag);
    sum(FAILED_flag)
];

% Create bar plot for main categories
subplot(1,2,1);
bar(main_counts, 'FaceColor', [0.2, 0.4, 0.6]);
set(gca, 'XTickLabel', main_categories, 'FontSize', 10, 'FontWeight', 'bold', 'XTickLabelRotation', 45);
ylabel('Count', 'FontSize', 12, 'FontWeight', 'bold');
title('Main Bifurcation Categories', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Add count numbers
for i = 1:length(main_counts)
    percentage = (main_counts(i)/numSamples)*100;
    text(i, main_counts(i), sprintf('%d\n(%.1f%%)', main_counts(i), percentage), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
         'FontSize', 9, 'FontWeight', 'bold');
end

% Create detailed breakdown for additional classification
subplot(1,2,2);
detailed_categories = {
    'LP+BP LP>BP', 
    'LP+BP LP<BP', 
    'LP-only S-shaped', 
    'LP-only Z-shaped', 
    'Biphasic', 
    'Other NO_SING'
};

detailed_counts = [
    sum(LPBP_LPgreaterBP);
    sum(LPBP_LPlessBP);
    sum(S_shaped_LP);
    sum(Z_shaped_LP);
    sum(NO_SING_biphasic);
    sum(NO_SING_other)
];

bar(detailed_counts, 'FaceColor', [0.6, 0.2, 0.4]);
set(gca, 'XTickLabel', detailed_categories, 'FontSize', 9, 'FontWeight', 'bold', 'XTickLabelRotation', 45);
ylabel('Count', 'FontSize', 12, 'FontWeight', 'bold');
title('Detailed Classification', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Add count numbers
for i = 1:length(detailed_counts)
    if detailed_counts(i) > 0
        percentage = (detailed_counts(i)/numSamples)*100;
        text(i, detailed_counts(i), sprintf('%d\n(%.1f%%)', detailed_counts(i), percentage), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontSize', 8, 'FontWeight', 'bold');
    end
end

% Add overall statistics text
annotation('textbox', [0.02, 0.02, 0.96, 0.05], ...
    'String', sprintf('Total Samples: %d | Negative in Continuation: %d (%.1f%%) | Analysis Time: %.1f seconds', ...
    numSamples, sum(NEGATIVE_CONT_flag), sum(NEGATIVE_CONT_flag)/numSamples*100, analysis_time), ...
    'FontSize', 10, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'BackgroundColor', [0.9, 0.9, 0.9]);

fprintf('\n=== Figure 6 Classification Summary ===\n');
fprintf('Main Categories:\n');
for i = 1:length(main_categories)
    fprintf('  %s: %d (%.1f%%)\n', main_categories{i}, main_counts(i), (main_counts(i)/numSamples)*100);
end

fprintf('\nDetailed Breakdown:\n');
for i = 1:length(detailed_categories)
    if detailed_counts(i) > 0
        fprintf('  %s: %d (%.1f%%)\n', detailed_categories{i}, detailed_counts(i), (detailed_counts(i)/numSamples)*100);
    end
end

%% ================= BLOCK 31: SAVE FIGURE 6 =================
fig6Path = fullfile(figuresDir, 'Figure6_Summary_Statistics.png');
saveas(fig6, fig6Path);
fprintf('Saved Figure 6 (Summary Statistics): %s\n', fig6Path);

%% ================= BLOCK 32: END TIMING AND FINAL SUMMARY =================
total_time = toc(total_start_time);
fprintf('\n=== Bifurcation Analysis Complete ===\n');
fprintf('Total analysis time: %.2f seconds (%.2f minutes)\n', total_time, total_time/60);
fprintf('All results saved in: %s\n', outputDir);
fprintf('Figures saved in: %s\n', figuresDir);
fprintf('Excel files saved in: %s\n', excelDir);
fprintf('MAT files saved in: %s\n', matDir);
fprintf('\n=== FINAL 7-CLASS CLASSIFICATION SUMMARY ===\n');
fprintf('1. LP_ONLY_S-shaped:     %d (%.1f%%)\n', class_counts_7(1), class_counts_7(1)/numSamples*100);
fprintf('2. LP_ONLY_Z-shaped:     %d (%.1f%%)\n', class_counts_7(2), class_counts_7(2)/numSamples*100);
fprintf('3. BP_ONLY:              %d (%.1f%%)\n', class_counts_7(3), class_counts_7(3)/numSamples*100);
fprintf('4. LP_BP_LPgreaterBP:    %d (%.1f%%)\n', class_counts_7(4), class_counts_7(4)/numSamples*100);
fprintf('5. LP_BP_LPlessBP:       %d (%.1f%%)\n', class_counts_7(5), class_counts_7(5)/numSamples*100);
fprintf('6. NO_SING_biphasic:     %d (%.1f%%)\n', class_counts_7(6), class_counts_7(6)/numSamples*100);
fprintf('7. OTHER (US):           %d (%.1f%%)\n', class_counts_7(7), class_counts_7(7)/numSamples*100);
fprintf('Total samples:           %d\n', numSamples);
fprintf('\n=== Files Created ===\n');
fprintf('1. Bifurcation_Analysis_Results.xlsx - Main results with all classifications\n');
fprintf('2. Classification_Datasets.xlsx - Detailed classification data\n');
fprintf('3. ML_Ready_Dataset.xlsx - ML-ready dataset with 7 classes\n');
fprintf('4. Negative_Continuation_Cases.xlsx - All negative continuation cases\n');
fprintf('5. Complete_Bifurcation_Results.mat - Complete MAT file with all data\n');

%% Classification of BP etc.

function print_BP_classification(all_results, sample_indices)
    % Print BP classification for specified samples
    fprintf('\n=== BP Classification Report ===\n');
    
    for k = sample_indices
        if k <= length(all_results.sample_details)
            fprintf('\nSample %d:\n', k);
            
            if isfield(all_results.sample_details(k), 'BP_coefficient_a')
                bp_coeffs = all_results.sample_details(k).BP_coefficient_a;
                bp_types = all_results.sample_details(k).BP_type;
                
                if ~isempty(bp_coeffs)
                    for j = 1:length(bp_coeffs)
                        if ~isnan(bp_coeffs(j))
                            fprintf('  BP %d: Coefficient a = %.4e -> %s\n', j, bp_coeffs(j), bp_types{j});
                        else
                            fprintf('  BP %d: %s\n', j, bp_types{j});
                        end
                    end
                else
                    fprintf('  No BP coefficients stored\n');
                end
            else
                fprintf('  No BP information available\n');
            end
        else
            fprintf('Sample %d does not exist\n', k);
        end
    end
end
