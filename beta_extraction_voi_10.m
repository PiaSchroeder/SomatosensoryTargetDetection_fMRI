%% Extract beta estimates for 10 Intensity levels
close all
clear all
clc

%% Data info

SJs         = { 'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S08' 'S09' 'S11' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S23' 'S24' 'S25' 'S26' 'S28' 'S29' 'S30' 'S32'};

data_dir    = 'my_data_dir';
model_dir   = 'BayesGLM';

trg_dir = fullfile(data_dir, '10Int betas');
if ~exist(trg_dir,'dir')
    mkdir(trg_dir)
end

name = 'Betas_4mmBMSPeakVOIs_10Int';
info = 'Beta estimates for 10 intensity levels from 4mm VOIs centered on individual BMS peaks within group level BMS ROIs, mean estimates';

% Bayesian 1st level models that define VOIs
models   = {'Int'
            'PF'
            'Det'
            'Unc'
            'Rep'};
        
% Bayesian 1st level model to extract parameters from
ffx = 'BayesFFXmpctd_10RT';
  
%                VOI                     name       model
VOIs    =   {   'BMS_Int_rSI_4mm'       'R SIa'     1   
                'BMS_Int_rSIIa_4mm'     'R SIIa'    1   
                'BMS_Int_rSIIp_4mm'     'R SIIp'    1   
                'BMS_Int_lSII_4mm'      'L SIIm'    1   
                'BMS_PF_rSI_4mm'        'R SIp'     2   
                'BMS_PF_rSII_4mm'       'R SII'     2   
                'BMS_PF_lSII_4mm'       'L SII'     2   
                'BMS_Det_rSIIs_4mm'     'R SIIs'    3   
                'BMS_Det_rSIIi_4mm'     'R SIIi'    3   
                'BMS_Det_lSII_4mm'      'L SIIl'    3   
                'BMS_Det_lIPS_4mm'      'L IPL'     3   
                'BMS_Det_lMFG_4mm'      'L SFG'     3   
                'BMS_Det_lLG_4mm'       'L V3'      3   
                'BMS_UncMSF_4mm'        'SMG/ACC'   4   
                'BMS_UncrAIC_4mm'       'R AIC'     4   
                'BMS_UnclAIC_4mm'       'L AIC'     4   
                'BMS_ReplSMA_4mm'       'L SMA'     5   
                'BMS_ReplThal_4mm'      'L Thal'    5   
                'BMS_ReprSMG_4mm'       'R SMaG'    5  };   

nVOI = size(VOIs,1);
nSubs = numel(SJs);
nInt = 10;
nRuns = 4;

voi_dirs = cellfun(@strcat,repmat({'VOIs\FFXmpctd_'},nVOI,1),models([VOIs{:,3}]),repmat({'RT'},nVOI,1),'UniformOutput',0);

% beta images
prefix_Cbeta = 'rCbeta_'; 
b = 1:4:37; % index beta images of regressors of interest

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
        P.(VOIs{v,1}).voi = VOIs{v,1};
        
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

reg = { INT
        PF
        DET
        UNC
        REP };

cols = {[0 1 0]     % Green: Intensity
        [0 0 1]     % Blue: P(Detection)
        [1 0 0]     % Red: Detection
        [0 1 1]     % Cyan: Uncertainty
        [1 0 1]};   % Magenta: Report
    
for v = 1:nVOI
    
    voi = VOIs{v,1};
    
    % Get data
    beta = P.(voi).pM_rs;
    SEbeta = P.(voi).SE_rs;
    
    % Get model
    [b,dev,stats] = glmfit(reg{VOIs{v,3}},beta','normal');
    mod = b(2)*reg{VOIs{v,3}}+b(1);
    
    f = figure;
    hold on
    set(gca,'FontSize',10,'FontName','Calibri')
    
    % Plot data and model
    plot(mod,'LineWidth',2,'Color', cols{VOIs{v,3}})
    ha = errorbar(beta,SEbeta,'ko','MarkerSize',3,'MarkerFaceColor',[0 0 0]);
    
    % Remove horizonzal lines from error bars
    hb = get(ha,'children');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    xleft = temp; xright = temp+1;
    Xdata(xleft) = 0;
    Xdata(xright) = 0;
    set(hb(2),'Xdata',Xdata)
    
    % Axes and Ticks
    top = beta+SEbeta;
    bottom = beta-SEbeta;
    min_b = min(bottom);
    max_b = max(top);
    range_b = max_b - min_b;
    margin = range_b/20;
    ylim([(min_b - margin) (max_b + margin)])
    ylim_min = (ceil(min_b*10)/10);
    ylim_max = (floor(max_b*10)/10);
    clear mod
    if mod(round((ylim_max-ylim_min)*10),2) ~= 0
        ylim_min = ylim_min+.1;
    end
    ylabs = [ylim_min  (ylim_max - (ylim_max-ylim_min)/2) ylim_max];
    set(gca,'YTick',ylabs)
    ylabs_str = cellstr(num2str(ylabs'));
    ylabs_str = cellfun(@strrep,ylabs_str, repmat({'0.'},numel(ylabs_str),1), repmat({'.'},numel(ylabs_str),1),'UniformOutput',0);
    set(gca,'YTickLabel',{})
    xlims = xlim;
    text(repmat(xlims(1)-0.01,numel(ylabs_str),1),ylabs,ylabs_str,'HorizontalAlignment','right','FontName','Calibri','FontSize',8)
    set(gca,'XTick',[])
    set(gca,'xcolor',[1 1 1])
    set(gca,'TickLength',[0.02,0.5])
    
    title(VOIs{v,2})
    box off
    
    % Figure size
    x0 = 1;
    y0 = 1;
    width = 1.5;
    height = 2;
    set(gcf,'Units','centimeters','position',[0,0,14,21])
    set(gca,'units','centimeters','position',[x0,y0,width,height])
    
    % Save
    saveas(f,fullfile(trg_dir,[name,'_',VOIs{v,1},'.fig']))
    betas = P.(voi);
    save(fullfile(trg_dir,[name,'_',VOIs{v,1},'.mat']),'betas')
end



