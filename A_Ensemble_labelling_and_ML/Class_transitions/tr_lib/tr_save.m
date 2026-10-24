function tr_save(fig, base)
%TR_SAVE  PNG (300 dpi) and, in MATLAB, vector PDF and .fig.
if exist('OCTAVE_VERSION', 'builtin')
    print(fig, [base '.png'], '-dpng', '-r300');
else
    try
        exportgraphics(fig, [base '.png'], 'Resolution', 300);
        exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
    catch
        print(fig, [base '.png'], '-dpng', '-r300'); print(fig, [base '.pdf'], '-dpdf', '-painters');
    end
    savefig(fig, [base '.fig']);
end
end
