function [Xs, br_id] = ss_states_at(S, a)
%SS_STATES_AT  All positive steady states at ATOT = a, from an ss_structure result.
%   Xs is 6 x n (n <= 3 for this network), br_id the branch index of each state.
C = S.constants; Xs = zeros(6,0); br_id = zeros(1,0);
for i = 1:numel(S.branch)
    G = S.branch(i).ATOT; z = S.branch(i).z;
    h = G - a;
    idx = find(h(1:end-1).*h(2:end) < 0);
    for m = idx(:)'
        map = S.branch(i).map;
        g = @(zz) Gz(zz, map, C) - a;
        zr = fzero(g, [z(m) z(m+1)], optimset('TolX',1e-14));
        [u, A] = map(zr);
        Xs(:,end+1) = ss_state_from_uA(u, A, C); %#ok<AGROW>
        br_id(end+1) = i; %#ok<AGROW>
    end
end
end

function G = Gz(z, map, C)
[u, A] = map(z);
[~, G] = ss_state_from_uA(u, A, C);
end