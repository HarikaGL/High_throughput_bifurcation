function s = tr_commas(n)
%TR_COMMAS  Integer with thousands separators.
s = sprintf('%d', n); k = numel(s) - 3;
while k > 0, s = [s(1:k) ',' s(k+1:end)]; k = k - 3; end
end
