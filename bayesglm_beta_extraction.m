%% Extract Bayesian beta estimates
close all
clear all
clc

%%

SJs         = { 'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S08' 'S09' 'S11' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S23' 'S24' 'S25' 'S26' 'S28' 'S29' 'S30' 'S32'};

data_dir    = 'my_data_dir';
model_dir   = 'BayesGLM';

voi_dirs    = { 'VOIs\FFXmpctd_IntRT'
                'VOIs\FFXmpctd_IntRT'
                'VOIs\FFXmpctd_IntRT'
                'VOIs\FFXmpctd_IntRT'
                'VOIs\FFXmpctd_PFRT'
                'VOIs\FFXmpctd_PFRT'
                'VOIs\FFXmpctd_PFRT'
                'VOIs\FFXmpctd_DetRT'
                'VOIs\FFXmpctd_DetRT'
                'VOIs\FFXmpctd_DetRT'
                'VOIs\FFXmpctd_DetRT'
                'VOIs\FFXmpctd_DetRT'
                'VOIs\FFXmpctd_DetRT'
                'VOIs\FFXmpctd_UncRT'
                'VOIs\FFXmpctd_UncRT'
                'VOIs\FFXmpctd_UncRT'
                'VOIs\FFXmpctd_RepRT'
                'VOIs\FFXmpctd_RepRT'
                'VOIs\FFXmpctd_RepRT'};

trg_dir = fullfile(data_dir, '2nd level', model_dir,'betas');
if ~exist(trg_dir,'dir')
    mkdir(trg_dir)
end
name = 'Betas_4mmBMSPeakVOIs_mean';
info = '4mm VOIs centered on BMS peak within group level model clusters, mean estimates';

% Bayesian 1st level model to extract parameters from
ffx = 'BayesFFXmpctd_10RT';

% VOIs must be created beforehand!
VOIs = {'BMS_Int_rSI_4mm'
        'BMS_Int_rSIIa_4mm'
        'BMS_Int_rSIIp_4mm'
        'BMS_Int_lSII_4mm'
        'BMS_PF_rSI_4mm'
        'BMS_PF_rSII_4mm'
        'BMS_PF_lSII_4mm'
        'BMS_Det_rSIIs_4mm'
        'BMS_Det_rSIIi_4mm'
        'BMS_Det_lSII_4mm'
        'BMS_Det_lIPS_4mm'
        'BMS_Det_lMFG_4mm'
        'BMS_Det_lLG_4mm'
        'BMS_UncMSF_4mm'
        'BMS_UncrAIC_4mm'
        'BMS_UnclAIC_4mm'
        'BMS_ReplSMA_4mm'
        'BMS_ReplThal_4mm'
        'BMS_ReprSMG_4mm'
    };

nVOI = size(VOIs,1);
nSubs = numel(SJs);
nInt = 10;
nRuns = 4;

prefix_Cbeta = 'rCbeta_'; 
b = 1:4:37; % index beta images of regressors of interest
r = 1:nRuns;

%% Get estimates
P = cell2struct(cell(size(VOIs,1),1),VOIs(:,1),1);
for v = 1:nVOI
    
    voi = VOIs{v,1};
    disp(voi)
    
    if exist(fullfile(trg_dir,[name '_' voi '.mat']),'file')
        load(fullfile(trg_dir,[name '_' voi '.mat']))
        fns = fieldnames(betas);
        val = struct2cell(betas);
        P.(VOIs{v,1}) = cell2struct(val,fns);
        disp(['Loaded ' name 'betas into P struct!'])

    else
        
        % Prepare data structure
        P.(VOIs{v,1}).pM = nan(nInt, nRuns, nSubs);
        P.(VOIs{v,1}).info = info;
        P.(VOIs{v,1}).voi = VOIs{v};
        
        % Assemble mean estimates for all voxels in VOIs, subjects, runs, and intensities
        for s = 1:nSubs
            disp(SJs{s})
            % get betas
            sub_dir = fullfile(data_dir, SJs{s}, model_dir, ffx);
            voi_file = fullfile(data_dir, SJs{s}, voi_dirs{v}, ['VOI_' SJs{s} '_' voi '_mask.nii']);
            for r = 1:nRuns
                disp(r)
                for i = 1:nInt
                    beta_idx = b(i)+(r-1)*56; 
                    Cbeta_img = spm_select('FPList', sub_dir, ['^' prefix_Cbeta '0{1,3}' num2str(beta_idx) '.nii$']);
                    P.(VOIs{v,1}).pM(i,r,s) = spm_summarise(Cbeta_img, voi_file, @mean);
                end
            end
        end
        
        % Aggregate runs
        P.(voi).pM_r = squeeze(mean(P.(voi).pM,2))';
              
        % Aggregate subjects
        P.(voi).pM_rs = mean(P.(voi).pM_r,1);
        P.(voi).SD_rs = std(P.(voi).pM_r,1);
        P.(voi).SE_rs = P.(voi).SD_rs/sqrt(nSubs);
    end
end



%% Plot

% Regressors
% intensity
INT = (1:10)';
% detection
det = fullfile('mydir','normDet.mat');
load(det)
DET = mean(det)';
% pfs
pfs = fullfile('mydir','normPFs.mat');
load(pfs)
PF = mean(pfs)';
% uncertainty
unc = fullfile('mydir','normUnc.mat');
load(unc)
UNC = mean(unc)';
% report 
REP = ones(10,1); 

reg = {INT,INT,INT,INT,PF,PF,PF,DET,DET,DET,DET,DET,DET,UNC,UNC,UNC,REP,REP,REP};

cols = {[0 1 0]
        [0 1 0]
        [0 1 0]
        [0 1 0]
        [0 0 1]
        [0 0 1]
        [0 0 1]
        [1 0 0]
        [1 0 0]
        [1 0 0]
        [1 0 0]
        [1 0 0]
        [1 0 0]
        [0 1 1]
        [0 1 1]
        [0 1 1]
        [1 0 1]
        [1 0 1]
        [1 0 1]};
    
for v = 1:nVOI
    
    voi = VOIs{v,1};
        
    beta = P.(voi).pM_rs;
    SEbeta = P.(voi).SE_rs;
    
    [b,dev,stats] = glmfit(reg{v},beta','normal');
    
    mod = b(2)*reg{v}+b(1);
    
    f = figure;
    hold on
        
    plot(mod,'LineWidth',6,'Color', cols{v})
    errorbar(beta,SEbeta,'ko','LineWidth',4)
        
    xlim([0.5 10.5])
    set(gca,'FontSize',18)
    title(voi)
    xlabel('Intensity level')
    ylabel('Beta estimates')
    box off
    
    % Axes
    switch voi
        case 'BMS_Int_rSI_4mm'
            ylim([-.9 .15]) % SI
            set(gca,'YTick',[-.9 -.4 .1])
        case 'BMS_Int_rSIIa_4mm'
            ylim([-.8 .6])
            set(gca,'YTick',[-.8 -.1 .6])
        case 'BMS_Int_rSIIp_4mm'
            ylim([-.6 1])
            set(gca,'YTick',[-.6 .2 1])
        case 'BMS_Int_lSII_4mm'
            ylim([-.7 .7])
            set(gca,'YTick',[-.7 0 .7])
        case 'BMS_PF_rSI_4mm'
            ylim([-.5 .55])
            set(gca,'YTick',[-.5 0 .5])
        case 'BMS_PF_rSII_4mm'
            ylim([-1 .8]) % rSII
            set(gca,'YTick',[-1 -.1 .8])
        case 'BMS_PF_lSII_4mm'
            ylim([-.65 .6])
            set(gca,'YTick',[-.6 0 .6])
        case 'BMS_Det_rSIIs_4mm'
            ylim([-.2 1])
            set(gca,'YTick',[-.2 .4 1])
        case 'BMS_Det_rSIIi_4mm'
            ylim([-1.05 .2])
            set(gca,'YTick',[-1 -.4 .2])
        case 'BMS_Det_lSII_4mm'
            ylim([-.05 .8])
            set(gca,'YTick',[0 .4 .8])
        case 'BMS_Det_lIPS_4mm'
            ylim([-.3 .3]) % lAG
            set(gca,'YTick',[-.3 0 .3])
        case 'BMS_Det_lMFG_4mm'
            ylim([-1.4 .6]) % lMFG
            set(gca,'YTick',[-1.4 -.4 .6])
        case 'BMS_Det_lLG_4mm'
            ylim([.8 1.8]) 
            set(gca,'YTick',[.8 1.3 1.8])
        case 'BMS_UnclAIC_4mm'
            ylim([-.2 1.2]) % ACC
            set(gca,'YTick',[-.2 .5 1.2])
        case 'BMS_UncrAIC_4mm'
            ylim([-.4 1.25])
            set(gca,'YTick',[-.4 .4 1.2])
        case 'BMS_UncMSF_4mm'
            ylim([-.6 1.6]) % ACC
            set(gca,'YTick',[-.4 .5 1.4])
        case 'BMS_ReplSMA_4mm'
            ylim([.8 2.2]) % SMA
            set(gca,'YTick',[.8 1.5 2.2])
        case 'BMS_ReplThal_4mm'
            ylim([.2 1]) % Thal
            set(gca,'YTick',[.2 .6 1])
        case  'BMS_ReprSMG_4mm'
            ylim([-.4 .8]) % SMG
            set(gca,'YTick',[-.4 .2 .8])
    end
 
    x0 = 0;
    y0 = 0;
    width = 15;
    height = 15;
    set(gcf,'units','centimeters','position',[x0,y0,width,height])
    
    set(gca,'XTick',1:10)
    saveas(f,fullfile(trg_dir,[name,'_',VOIs{v,1},'.fig']))
    betas = P.(voi);
    save(fullfile(trg_dir,[name,'_',VOIs{v,1},'.mat']),'betas')
end



