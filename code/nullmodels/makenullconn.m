cd(homedir); addpath(genpath(pwd));
load(fullfile(savedir,'W.mat'));

%% degree sequence preserving rewiring
Wrw_in = dir_generate_srand(W',1e8)';
Wrw_out = dir_generate_srand(W,1e8);

save(fullfile(savedir,'Wrewire.mat'),'Wrw_in','Wrw_out');