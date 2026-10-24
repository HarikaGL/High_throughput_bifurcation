function tr_arrowhead(ax, P, u, col)
%TR_ARROWHEAD  Filled arrowhead at point P pointing along the unit vector u.
s = 0.16; n = [-u(2) u(1)];
X = [P(1), P(1) - s*u(1) + 0.5*s*n(1), P(1) - s*u(1) - 0.5*s*n(1)];
Y = [P(2), P(2) - s*u(2) + 0.5*s*n(2), P(2) - s*u(2) - 0.5*s*n(2)];
patch(ax, X, Y, col, 'EdgeColor', col);
end
