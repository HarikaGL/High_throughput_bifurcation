function C = cm_graded(n)
%CM_GRADED  Graded colour scale for the confusion matrices (white -> pale teal -> teal ->
%   deep blue), replacing the black-and-white scale. C: n x 3 (default n = 256).
if nargin < 1, n = 256; end
K = [1.00 1.00 1.00; 0.87 0.95 0.95; 0.62 0.85 0.85; 0.33 0.68 0.76; 0.16 0.47 0.66; 0.08 0.27 0.50; 0.04 0.13 0.33];
t = linspace(0, 1, size(K, 1)); ti = linspace(0, 1, n);
C = [interp1(t, K(:,1), ti)', interp1(t, K(:,2), ti)', interp1(t, K(:,3), ti)'];
end
