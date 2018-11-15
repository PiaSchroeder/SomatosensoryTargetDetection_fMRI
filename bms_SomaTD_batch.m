%% Bayesian Model Selection
clc
close all
clear all

%%

data_dir = 'my_data_dir';    

SJs      = { 'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S08' 'S09' 'S11' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S23' 'S24' 'S25' 'S26' 'S28' 'S29' 'S30' 'S32'};

name     = 'IntDetPFUncRep_RT_2x2x2';

models   = {'mpctd_IntRT'
            'mpctd_DetRT'
            'mpctd_PFRT'
            'mpctd_UncRT'
            'mpctd_RepRT'};
       
level = 'RFX';
disp(name)
bms_SomaTD(data_dir,name,SJs,models,level);
