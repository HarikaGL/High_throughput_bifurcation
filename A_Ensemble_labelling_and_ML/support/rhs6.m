function [dx, J] = rhs6(x, P)
%RHS6  Right-hand side of MAPKfull.m (for residual checks outside MATCONT).
ATOT=P(1); KTOT=P(2); PTOT=P(3); k1=P(4); k11=P(5); k2=P(6); k3=P(7); k33=P(8); k4=P(9);
k5=P(10); k55=P(11); k6=P(12); k7=P(13); k77=P(14); k8=P(15);
A = ATOT - sum(x); K = KTOT - x(3) - x(4); Pf = PTOT - x(5) - x(6);
v1=k1*A*K-k11*x(3); v2=k2*x(3); v3=k3*x(1)*K-k33*x(4); v4=k4*x(4);
v5=k5*x(2)*Pf-k55*x(5); v6=k6*x(5); v7=k7*x(1)*Pf-k77*x(6); v8=k8*x(6);
dx = [v2-v3+v6-v7; v4-v5; v1-v2; v3-v4; v5-v6; v7-v8];
if nargout > 1, J = jac6(x, P); end
end