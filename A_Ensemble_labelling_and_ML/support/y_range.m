function yl = y_range(curves, xl)
%Y_RANGE  y-limits [0, 1.08*max x2] of the branch curves inside the ATOT range xl.
m = 0;
for c = curves(:)'
    in = c.ATOT >= xl(1) & c.ATOT <= xl(2) & isfinite(c.x2);
    if any(in), m = max(m, max(c.x2(in))); end
end
if m <= 0, m = 1; end
yl = [0 1.08*m];
end
