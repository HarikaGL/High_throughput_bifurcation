function [c, S, L] = tr_label(p, N)
%TR_LABEL  Exact ten-class label (v3) of p = [KTOT PTOT k1 k11 k2 k3 k33 k4 k5 k55 k6 k7 k77 k8].
%   The labeller and settings of Ensemble_LHS_v3: ss_structure + exact_labels, ATOT window
%   [0, 1e6], N grid points per branch (default 10000). Normal forms are not computed.
%   c = class 1..11 (11 = OTHER), NaN if the labeller fails.
if nargin < 2, N = 10000; end
try
    S = ss_structure([1 p(:)'], struct('ATOT_min', 0, 'ATOT_max', 1e6, 'N', N));
    L = exact_labels(S); c = L.class;
catch
    S = []; L = []; c = NaN;
end
end
