function [c, S] = tr_label_win(p, Amax, N)
%TR_LABEL_WIN  Exact class of p with the labelling window widened to [0, Amax]
%   (same labeller as tr_label; tr_label uses Amax = 1e6). Used to show which class a set
%   would have if the window included its detached fold.
if nargin < 3, N = 10000; end
try
    S = ss_structure([1 p(:)'], struct('ATOT_min', 0, 'ATOT_max', Amax, 'N', N));
    L = exact_labels(S); c = L.class;
catch
    S = []; c = NaN;
end
end
