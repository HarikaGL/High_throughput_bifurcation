function s = tr_texify(s)
%TR_TEXIFY  Label of TR_04 (TeX of MATLAB text) -> LaTeX for the table.
s = strrep(s, '%', '\%');
s = strrep(s, char(215), '$\times$'); s = strrep(s, '\div', '$\div$');
s = regexprep(s, '(K_\{TOT\}|P_\{TOT\}|k_\{?\d+\}?)', '$$1$');
s = strrep(s, '$K_{TOT}$', '$\KTOT$'); s = strrep(s, '$P_{TOT}$', '$\PTOT$');
s = strrep(s, '$$', '');
end
