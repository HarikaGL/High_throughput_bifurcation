function [X, G] = ss_state_from_uA(u, A, C)
%SS_STATE_FROM_UA  Steady state (6 x n, MAPKfull.m ordering) and its ATOT from (u = K/P, A).
K  = C.KTOT ./ (1 + A./C.Km1 + C.a1.*u.*A./C.Km3);
Pf = K ./ u;
x1 = C.a1.*u.*A;
x2 = C.a1.*C.a2.*u.^2.*A;
x3 = A.*K./C.Km1;
x4 = x1.*K./C.Km3;
x5 = x2.*Pf./C.Km5;
x6 = x1.*Pf./C.Km7;
X = [x1; x2; x3; x4; x5; x6];
G = A + sum(X,1);
end