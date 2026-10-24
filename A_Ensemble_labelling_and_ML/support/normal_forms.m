function NF = normal_forms(x, P, kind)
%NORMAL_FORMS  Normal-form coefficients of a fold (LP) or branch point (BP) of the dual
%   futile cycle (MAPKfull.m, 6 state variables, parameter ATOT = P(1)).
%
%   NF = normal_forms(x, P, 'LP')   x = 6x1 state at the fold, P(1) = ATOT at the fold
%   NF = normal_forms(x, P, 'BP')   x = 6x1 state at the branch point
%
%   The exact location comes from the inverse method (ss_structure), so the
%   coefficients are computed for EVERY singularity of every sample. The right-hand
%   side is quadratic in (x, ATOT), so all second derivatives below are exact
%   (polarisation identities, no step-size error).
%
%   LP (fold, saddle-node). J q = 0, p' J = 0, normalised as in MATCONT
%   (||q|| = 1, p'q = 1; sign fixed by q(2) > 0, MATCONT's sign is arbitrary):
%       a     = 1/2 p' B(q,q)        quadratic coefficient. SAME definition and
%                                    normalisation as MATCONT's s.data.a (a_lp.m,
%                                    nf_LP.m), so |a| can be compared directly
%       beta  = p' dF/dATOT          transversality coefficient
%       kappa = -a/beta              ATOT - ATOT_LP ~ kappa*xi^2 along the branch
%                                    (sign-independent; > 0 fold opens to larger ATOT)
%       kappa_rel = kappa*ATOT_LP    dimensionless curvature
%   Non-degenerate fold: a ~= 0 and beta ~= 0.
%
%   BP (branch point). J q = 0, p' J = 0 (||q|| = 1, p'q = 1, q(2) > 0); p'dF/dATOT = 0.
%   v0 solves J v0 = -dF/dATOT with q'v0 = 0, Q = (q,0), V = (v0,1) in (x, ATOT).
%   Algebraic branching equation  a xi^2 + 2 b xi eta + c eta^2 = 0  with
%       a = p' B(Q,Q),  b = p' B(Q,V),  c = p' B(V,V)   (B: second derivative in (x,ATOT))
%       disc = b^2 - a c  (> 0: two branches cross transversally)
%       slopes = dx2/dATOT of the two crossing branches at the BP (Inf = vertical,
%                i.e. a branch that turns back at the BP)
%       pitchfork = |a| <= 1e-6 max(|b|,|c|): one branch turns at the BP (symmetric
%                   case); otherwise a transcritical crossing
%   MATCONT reports no normal-form coefficient for BP, hence this computation.

F = @(xx, AT) rhs6(xx, [AT, P(2:end)]);
A0 = P(1);
J = jac6(x, P);
[U, Sg, V] = svd(J);
q = V(:,end); p = U(:,end);
q = q/norm(q); if q(2) < 0, q = -q; end
p = p/(p'*q);
sv = diag(Sg);
dA = max(1, abs(A0));
FA = (F(x, A0 + dA) - F(x, A0))/dA;          % exact: F is linear in ATOT
NF = struct();
NF.sigma_rel = sv(end)/sv(1);                 % ~0 at a genuine singularity
switch upper(kind)
    case 'LP'
        Bqq = F(x+q, A0) + F(x-q, A0) - 2*F(x, A0);   % exact: F quadratic in x
        NF.a     = 0.5*p'*Bqq;
        NF.beta  = p'*FA;
        NF.kappa = -NF.a/NF.beta;
        NF.kappa_rel = NF.kappa*A0;
        NF.q = q; NF.p = p;
    case 'BP'
        NF.pFA = p'*FA/max(norm(FA), eps);     % ~0 at a branch point
        M = [J, p; q', 0];
        sol = M \ [-FA; 0];
        v0 = sol(1:6);
        Qy = [q; 0]; Vy = [v0; 1]; y0 = [x; A0];
        Fy = @(y) rhs6(y(1:6), [y(7), P(2:end)]);
        B = @(u, w) (Fy(y0+u+w) - Fy(y0+u-w) - Fy(y0-u+w) + Fy(y0-u-w))/4;
        NF.a = p'*B(Qy, Qy);
        NF.b = p'*B(Qy, Vy);
        NF.c = p'*B(Vy, Vy);
        NF.disc = NF.b^2 - NF.a*NF.c;
        NF.pitchfork = abs(NF.a) <= 1e-6*max(abs(NF.b), abs(NF.c));
        % tangent directions xi*Q + eta*V -> slopes in the (ATOT, x2) plane
        if NF.pitchfork || NF.a == 0
            r = -NF.c/(2*NF.b);                % second branch (eta ~= 0)
            NF.slopes = [Inf, r*q(2) + v0(2)];
        elseif NF.disc >= 0
            r = (-NF.b + [1 -1]*sqrt(NF.disc))/NF.a;   % xi/eta
            NF.slopes = r*q(2) + v0(2);
        else
            NF.slopes = [NaN NaN];
        end
        NF.q = q; NF.p = p; NF.v0 = v0;
    otherwise
        error('normal_forms: kind must be LP or BP');
end
end
