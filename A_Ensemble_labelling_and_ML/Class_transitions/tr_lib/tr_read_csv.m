function M = tr_read_csv(f)
%TR_READ_CSV  Numeric CSV with one header line (MATLAB: readmatrix; Octave: dlmread).
if exist('OCTAVE_VERSION', 'builtin'), M = dlmread(f, ',', 1, 0);
else, M = readmatrix(f, 'NumHeaderLines', 1, 'Delimiter', ','); end
end
