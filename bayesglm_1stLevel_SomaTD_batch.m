%% Specify and estimate Bayesian 1st level GLMs 
clc
close all
clear all
                  
%% 

ana_dir     = 'my_ana_dir';
src_dir     = 'my_src_dir';  

SJs         = { 'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S08' 'S09' 'S11' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S23' 'S24' 'S25' 'S26' 'S28' 'S29' 'S30' 'S32'};
runs        = 1:4;

ffx_names   = {'BayesFFXmpctd_IntRT'
               'BayesFFXmpctd_PFRT'
               'BayesFFXmpctd_DetRT'
               'BayesFFXmpctd_UncRT'
               'BayesFFXmpctd_RepRT'
               'BayesFFXmpctd_10RT' };

% cycle over subjects
for m = 1:numel(ffx_names)
    for s = 1:numel(SJs)
        cd(ana_dir);
        data_dir = fullfile(src_dir, SJs{s});
        bayesglm_1stLevel_SomaTD(data_dir,ffx_names{m},runs,SJs{s});
    end
end

cd(ana_dir);
