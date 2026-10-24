function tr_write_csv(f, header, M)
%TR_WRITE_CSV  Numeric matrix with one header line (MATLAB and Octave).
fid = fopen(f, 'w'); fprintf(fid, '%s\n', header);
if ~isempty(M), fprintf(fid, [repmat('%.10g,', 1, size(M,2)-1) '%.10g\n'], M'); end
fclose(fid);
end
