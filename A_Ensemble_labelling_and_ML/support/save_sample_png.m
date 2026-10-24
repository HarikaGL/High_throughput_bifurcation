function save_sample_png(R, P, k, folder)
%SAVE_SAMPLE_PNG  Write the bifurcation diagram of one sample to a PNG and a .fig
%   (works on parallel workers: the figure is invisible).
if ~exist(folder, 'dir'), mkdir(folder); end
f = figure('Visible','off', 'Color','w', 'Position',[100 100 640 460]);
cleanup = onCleanup(@() close(f));
plot_sample(R, P, k, axes(f));
base = fullfile(folder, sprintf('sample_%05d_%s', k, R.exact_label));
saveas(f, [base '.png']);
set(f, 'Visible', 'on');   % so the .fig opens visible
savefig(f, [base '.fig']);
end
