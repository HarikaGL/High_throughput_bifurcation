function J = jac6(x, P)
%JAC6  Analytic Jacobian of the reduced 6-variable system (identical to MAPKfull.m).
ATOT=P(1); KTOT=P(2); PTOT=P(3); k1=P(4); k11=P(5); k2=P(6); k3=P(7); k33=P(8); k4=P(9);
k5=P(10); k55=P(11); k6=P(12); k7=P(13); k77=P(14); k8=P(15);
s  = sum(x) - ATOT;          % = -A
cK = x(3) + x(4) - KTOT;     % = -K
cP = x(5) + x(6) - PTOT;     % = -P
J = [ k3*cK + k7*cP, 0, k2 + x(1)*k3, k33 + x(1)*k3, k6 + x(1)*k7, k77 + x(1)*k7;
      0, k5*cP, 0, k4, k55 + x(2)*k5, x(2)*k5;
      k1*cK, k1*cK, k1*s - k11 - k2 + k1*cK, k1*s + k1*cK, k1*cK, k1*cK;
     -k3*cK, 0, -x(1)*k3, -k4 - k33 - x(1)*k3, 0, 0;
      0, -k5*cP, 0, 0, -k6 - k55 - x(2)*k5, -x(2)*k5;
     -k7*cP, 0, 0, 0, -x(1)*k7, -k8 - k77 - x(1)*k7 ];
end