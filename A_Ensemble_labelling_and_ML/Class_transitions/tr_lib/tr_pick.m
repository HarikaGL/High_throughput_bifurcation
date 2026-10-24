function idx = tr_pick(n, k, state)
%TR_PICK  k distinct random indices from 1..n (portable Park-Miller generator: MATLAB and
%   GNU Octave pick the same sets).
a = 16807; M = 2147483647; u = zeros(n, 1);
for i = 1:n, state = mod(a*state, M); u(i) = state/M; end
[~, o] = sort(u); idx = o(1:min(k, n));
end
