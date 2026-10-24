function Q = tr_conditions(P)
%TR_CONDITIONS  The exact quantities that separate neighbouring classes, for each row of P
%   (n x 14, [KTOT PTOT k1 k11 k2 k3 k33 k4 k5 k55 k6 k7 k77 k8]).
%   Q.x    = log10(k2 KTOT / (k8 PTOT))   site 1: detached branch needs x < 0
%   Q.y    = log10(k4 KTOT / (k6 PTOT))   site 2: detached branch needs y > 0
%   Q.det  = x < 0 & y > 0                criterion (10): a detached branch exists
%   Q.bp   = log10(k4 / k2)               symmetric sets: branch point iff > 0
%   Q.km   = log10(KM1 / KM5)             branch-point sets: Class 7 iff > 0 (KM5 < KM1)
%   Q.side = sign(alpha r + beta)         +1: main branch on u > r (S side), -1: u < r
%   Q.sym  = |u* - r| <= 1e-9 r           symmetric (degenerate) set, as in ss_structure
%   Q.a, Q.b = log10(k1/k5), log10(k11/k55)
KT=P(:,1); PT=P(:,2); k1=P(:,3); k11=P(:,4); k2=P(:,5); k3=P(:,6); k33=P(:,7); k4=P(:,8);
k5=P(:,9); k55=P(:,10); k6=P(:,11); k7=P(:,12); k77=P(:,13); k8=P(:,14);
KM1=(k11+k2)./k1; KM3=(k33+k4)./k3; KM5=(k55+k6)./k5; KM7=(k77+k8)./k7;
r=KT./PT; a1=k2.*KM7./(k8.*KM1); a2=k4.*KM5./(k6.*KM3);
al=r.*a1.*a2./KM5 - a1./KM3; be=r.*a1./KM7 - 1./KM1;        % as model_constants
Q.x = log10(k2.*KT./(k8.*PT)); Q.y = log10(k4.*KT./(k6.*PT));
Q.det = Q.x < 0 & Q.y > 0;
Q.bp = log10(k4./k2); Q.km = log10(KM1./KM5);
Q.side = sign(al.*r + be);
Q.sym = abs(al.*r + be) <= 1e-9*(abs(al).*r + abs(be));
Q.a = log10(k1./k5); Q.b = log10(k11./k55);
end
