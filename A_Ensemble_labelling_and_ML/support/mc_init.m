function mc_init(matcont_dir, project_dir)
%MC_INIT  Initialise MATCONT in the calling MATLAB process (client or parallel worker).
%   matcont_dir : folder with MATCONT's init.m ('' = init.m is already on the path,
%                 as in the original script, which simply called `init`)
%   project_dir : folder with MAPKfull.m and the helper functions (added to the path)
if nargin >= 2 && ~isempty(project_dir), addpath(project_dir); end
here = pwd;
cleanup = onCleanup(@() cd(here));
if nargin >= 1 && ~isempty(matcont_dir)
    cd(matcont_dir);
end
init;              % MATCONT's init.m
end