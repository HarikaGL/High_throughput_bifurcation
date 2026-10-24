function R = old_analyze_sample(P0, opt, npts)
% OLD_ANALYZE_SAMPLE  Legacy method, one parameter set.
% Faithful port of the per-sample body of the simple ensemble script:
%   ode15s to t = 5000 -> ONE steady state -> init_EP_EP -> single forward MATCONT
%   continuation in ATOT from that state -> for every BP, init_BP_EP and continue the
%   emanating branch in BOTH directions -> classify (7 classes).
% Detached LP branches are never reached, because only one steady state is continued.
%
% Kept exactly as in the legacy code (so the comparison is against what was published):
%   * any state value < 0 (even round-off) flags the sample NEGATIVE
%   * Hopf is tested as 'HB' (MATCONT writes 'H ')
%   * biphasic = exactly one findpeaks maximum of x2 along the main branch
% Not ported (not needed for the classes): BP normal-form coefficient extraction.
%
% Class numbers: 1 LP_ONLY_S, 2 LP_ONLY_Z, 3 BP_ONLY, 4 LP_BP LP>BP, 5 LP_BP LP<BP,
%                6 NO_SING_biphasic, 7 OTHER, 8 FAILED or NEGATIVE
if nargin < 3, npts = 600; end
nStates = 6; ap1 = 1; state_row = 2;

R = struct('ok',false,'status','FAILED','category','FAILED','class7',8, ...
    'LP_ATOT',[],'LP_x2',[],'BP_ATOT',[],'BP_x2',[],'ext_LP_n',0, ...
    'negative',false,'min_state',NaN,'ATOT_reached',NaN, ...
    'main',[],'ext',{{}},'message','');
P0 = P0(:)';

%% initial steady state and main continuation
try
    hls   = feval(@MAPKfull);
    Pcell = num2cell(P0);
    [~,y] = ode15s(hls{2}, [0 5000], zeros(nStates,1), odeset('RelTol',1e-8), Pcell{:});
    IC    = y(end,:)';
    x0    = init_EP_EP(@MAPKfull, IC, P0, ap1);
catch ME
    R.message = ['init: ' ME.message]; return
end
try
    [x,~,s] = quiet_cont(x0, [], opt);
catch ME
    R.message = ['cont: ' ME.message]; return
end
param_row = size(x,1);
st = x(1:nStates,:);
R.min_state = min(st(:));
R.negative  = any(st(:) < 0);

labels = {s.label};
hasLP = any(strcmp(labels,'LP'));
hasBP = any(strcmp(labels,'BP'));
hasHB = any(strcmp(labels,'HB'));
noSing = ~(hasLP || hasBP || hasHB);

lpA = []; lpX = []; bpA = []; bpX = [];
for j = find(strcmp(labels,'LP'))
    ii = s(j).index; lpA(end+1) = x(param_row,ii); lpX(end+1) = x(state_row,ii); %#ok<AGROW>
end
for j = find(strcmp(labels,'BP'))
    ii = s(j).index; bpA(end+1) = x(param_row,ii); bpX(end+1) = x(state_row,ii); %#ok<AGROW>
end

%% extension from every BP (both directions) to find LPs on the emanating branches
extA = []; extX = []; ext = {};
reached = max(x(param_row,:));
for kk = find(strcmp(labels,'BP'))
    try
        xBP = x(1:nStates, s(kk).index);
        pt  = P0; pt(ap1) = x(param_row, s(kk).index);
        [xx,vx] = init_BP_EP(@MAPKfull, xBP, pt, s(kk), 1e-8);
        for bw = [1 0]
            o2 = contset(opt,'Backward',bw);
            o2 = contset(o2,'MaxNumPoints',1500);
            [x2,~,s2] = quiet_cont(xx, vx, o2);
            if isempty(x2), continue; end
            ext{end+1} = thin(x2([end, state_row],:), npts); %#ok<AGROW>
            reached = max(reached, max(x2(end,:)));
            if ~isempty(s2)
                l2 = {s2.label};
                for j = find(strcmp(l2,'LP'))
                    ii = s2(j).index;
                    extA(end+1) = x2(end,ii); extX(end+1) = x2(state_row,ii); %#ok<AGROW>
                end
            end
        end
    catch
        % extension failed: skip, as in the legacy code
    end
end
allLPa = [lpA extA]; allLPx = [lpX extX];

%% classification (legacy logic)
cls = 7; cat = 'OTHER';
if R.negative
    cat = 'NEGATIVE'; cls = 8;
elseif noSing
    cat = 'NO_SING';
    sv = x(state_row,:);
    if numel(sv) > 10
        try
            if numel(findpeaks(sv)) == 1, cls = 6; end
        catch
        end
    end
else
    hasAnyLP = hasLP || ~isempty(extA);
    if hasBP && hasAnyLP
        cat = 'LP_BP';
        if mean(allLPx) > mean(bpX), cls = 4; else, cls = 5; end
    elseif hasBP
        cat = 'BP_ONLY'; cls = 3;
    elseif hasAnyLP
        cat = 'LP_ONLY';
        if numel(allLPx) >= 2
            [~,o] = sort(allLPa); cs = allLPx(o);
            if cs(1) > cs(2), cls = 1; else, cls = 2; end
        else
            cls = 7;
        end
    end
end

R.ok = true;  R.status = 'OK'; if R.negative, R.status = 'NEGATIVE'; end
R.category = cat; R.class7 = cls;
R.LP_ATOT = allLPa; R.LP_x2 = allLPx; R.BP_ATOT = bpA; R.BP_x2 = bpX;
R.ext_LP_n = numel(extA);
R.ATOT_reached = reached;
R.main = thin(x([param_row, state_row],:), npts);
R.ext  = ext;
end

function [x,v,s,h,f] = quiet_cont(x0, v0, opt)
% cont() without the console chatter
[~] = evalc('[x,v,s,h,f] = cont(@equilibrium, x0, v0, opt);');
end

function c = thin(c, n)
m = size(c,2);
if m > n, c = c(:, unique(round(linspace(1,m,n)))); end
end