function C = model_constants(P)
%MODEL_CONSTANTS  Michaelis constants and reduction coefficients of the dual futile cycle.
%   P = [ATOT KTOT PTOT k1 k11 k2 k3 k33 k4 k5 k55 k6 k7 k77 k8] (same order as MAPKfull.m)
C.KTOT = P(2); C.PTOT = P(3);
k1=P(4); k11=P(5); k2=P(6); k3=P(7); k33=P(8); k4=P(9);
k5=P(10); k55=P(11); k6=P(12); k7=P(13); k77=P(14); k8=P(15);
C.Km1 = (k11+k2)/k1;   % A  + K <-> x3 -> x1 + K
C.Km3 = (k33+k4)/k3;   % x1 + K <-> x4 -> x2 + K
C.Km5 = (k55+k6)/k5;   % x2 + P <-> x5 -> x1 + P
C.Km7 = (k77+k8)/k7;   % x1 + P <-> x6 -> A  + P
C.a1 = k2*C.Km7/(k8*C.Km1);
C.a2 = k4*C.Km5/(k6*C.Km3);
C.r  = C.KTOT/C.PTOT;
C.alpha = C.r*C.a1*C.a2/C.Km5 - C.a1/C.Km3;
C.beta  = C.r*C.a1/C.Km7 - 1/C.Km1;
if C.alpha ~= 0, C.ustar = -C.beta/C.alpha; else, C.ustar = Inf; end
end