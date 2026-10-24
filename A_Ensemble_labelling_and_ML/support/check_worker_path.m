function ok = check_worker_path(supportDir)
%CHECK_WORKER_PATH  true if this MATLAB process (client or worker) uses the v3 helpers
%   in supportDir (and not an older copy of analyze_sample.m from another folder).
ok = strcmp(fileparts(which('analyze_sample')), supportDir) && ...
     strcmp(fileparts(which('classify_branches')), supportDir);
end
