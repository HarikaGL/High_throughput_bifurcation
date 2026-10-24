function ex = tr_window_examples(W, c)
%TR_WINDOW_EXAMPLES  Two example window sets of class c (indices into W):
%   (a) the set whose detached fold is nearest to 10^7 (just beyond the window);
%   (b) the set deepest inside the detached-branch quadrant (largest min(-x, y)).
%   Fixed rule, so that MATLAB and Octave pick the same sets.
k = find(W.class == c & isfinite(W.detFold) & W.detFold > 1e6);
[~, j] = min(abs(log10(W.detFold(k)) - 7)); a = k(j);
k2 = k(k ~= a); [~, j] = max(W.margin(k2)); b = k2(j);
ex = [a b];
end
