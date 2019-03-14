%% Extract beta estimates for experimental regressors
clear all
close all
clc

%% Data info

SJs         = { 'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S08' 'S09' 'S11' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S23' 'S24' 'S25' 'S26' 'S28' 'S29' 'S30' 'S32'};

bms_name    = 'BMS_FFX_IntDetPFUncRep_RT_2x2x2';

data_dir    = 'my_data_dir';                            % Data directory
roi_dir     = 'my_roi_dir';                             % Directory containnig masks of ROIs
model_dir   = 'BayesGLM';                               % Directory containing first level beta images
map_dir     = fullfile(model_dir, 'BMS', bms_name);     % Directory containing single subject BMS results (PPMs)

trg_dir     = fullfile(data_dir, 'ExpReg betas');
if ~exist(trg_dir,'dir')
    mkdir(trg_dir)
end

name        = 'BayesGLM_Betas_BMSPeaks_ExpRegs';
info        = 'Beta estimates of experimental regressors from individual BMS peaks within group level BMS ROIs';

% Bayesian 1st level models to extract parameters from
models   = {'Int'
            'PF'
            'Det'
            'Unc'
            'Rep'};

mod = cellfun(@strcat,repmat({'mpctd_'},numel(models),1),models,repmat({'RT'},numel(models),1),'UniformOutput',0); % Complete file names

%                ROI             name       model
ROIs    =   {   'Int_rSI'       'R SIa'     1   
                'Int_rSIIp'     'R SIIp'    1   
                'Int_rSIIa'     'R SIIa'    1   
                'Int_lSII'      'L SIIm'    1   
                'PF_rSI'        'R SIp'     2   
                'PF_rSII'       'R SII'     2   
                'PF_lSII'       'L SII'     2   
                'Det_rSIIs'     'R SIIs'    3   
                'Det_rSIIi'     'R SIIi'    3   
                'Det_lSII'      'L SIIl'    3   
                'Det_lMFG'      'L SFG'     3   
                'Det_lIPS'      'L IPL'     3   
                'Det_lLG'       'L V3'      3   
                'Unc_MSF'       'SMG/ACC'   4   
                'Unc_lAIC'      'L AIC'     4   
                'Unc_rAIC'      'R AIC'     4   
                'Rep_lSMA'      'L SMA'     5   
                'Rep_lThal'     'L Thal'    5   
                'Rep_rSMG'      'R SMaG'    5  };                                        

nROI = size(ROIs,1);
nSub = numel(SJs);
nRuns = 4;

% beta images
prefix_Cbeta = 'srCbeta_';
prefix_SDbeta = 'srSDbeta_';
b = 3; % beta idx expreg

%% Get estimates
if exist(fullfile(trg_dir,[name '.mat']),'file')
    load(fullfile(trg_dir,[name '.mat']))
    disp(['Loaded ' name 'betas!'])
else
    P = cell2struct(cell(size(ROIs,1),1),ROIs(:,1),1);
    P.info = info;
    for r = 1:nROI
        
        roi = ROIs{r,1};
        m = ROIs{r,3};
        model = ['BayesFFX' mod{m}];
        map = mod{m};
        disp(roi)
        
        % Get peaks
        ROImax = find_bms_peaks(map, map_dir, roi, roi_dir, SJs);
        coords = ROImax.coords';
        peaks = cell2mat(coords(:,2))';
        
        % Prepare data structure
        P.(ROIs{r,1}).beta = nan(nSub, nRuns);
        P.(ROIs{r,1}).sd = nan(nSub, nRuns);
        P.(ROIs{r,1}).roi = ROIs{r,1};
        
        % Assemble beta estimates for BMS peaks in ROIs per subject and run
        for s = 1:nSub
            disp(SJs{s})
            % get betas
            sub_dir = fullfile(data_dir, SJs{s}, model_dir, model);
            for r = 1:nRuns
                beta_idx = b+(r-1)*22;
                Cbeta_img = spm_select('FPList', sub_dir, ['^' prefix_Cbeta '0{1,3}' num2str(beta_idx) '.nii$']);
                SDbeta_img = spm_select('FPList', sub_dir, ['^' prefix_SDbeta '0{1,3}' num2str(beta_idx) '.nii$']);
                P.(ROIs{r,1}).beta(s,r) = spm_summarise(Cbeta_img, peaks(:,s));
                P.(ROIs{r,1}).sd(s,r) = spm_summarise(SDbeta_img, peaks(:,s));
            end
        end
        
        % Aggregate runs
        P.(roi).beta_r = squeeze(mean(P.(roi).beta,2))';
        P.(roi).sd_r = squeeze(mean(P.(roi).sd,2))';
        
        % Aggregate subjects
        P.(roi).beta_rs = mean(P.(roi).beta_r,1);
        P.(roi).sd_rs = mean(P.(roi).sd_r,1);
    end
end

save(fullfile(trg_dir,[name '.mat']),'P');


%% Perform Bayesian t-tests on betas

bf10 = nan(nROI,1);
pValue = nan(nROI,1);
fprintf('\n\n Bayes factors and p-values \n\n')

for r = 1:nROI
    [bf10(r),pValue(r)] = bf_ttest(P.(ROIs{r,1}).beta_r);
    bf01 = 1./bf10;
    fprintf('%s mean beta = %.3f BF10 = %.3f BF01 = %.3f p = %.3f\n',ROI_names{r},mean(P.(ROIs{r,1}).beta_r),bf10(r),bf01(r),pValue(r));
end

%% Scatter plot per ROI

betas = nan(nSub,nROI);
for r = 1:nROI
    betas(:,r) = P.(ROIs{r,1}).beta_r;
end

cols = {[0 1 0]     % Green: Intensity
        [0 0 1]     % Blue: P(Detection)
        [1 0 0]     % Red: Detection
        [0 1 1]     % Cyan: Uncertainty
        [1 0 1]};   % Magenta: Report

for r = 1:nROI
    
    f = figure;
    hold on
    set(gca,'FontSize',8,'FontName','Calibri')
    
    % Scatter plot
    scatter(ones(nSub,1),betas(:,r),5, 'MarkerFaceColor','none','MarkerEdgeColor',cols{ROIs{r,3},:},'LineWidth',1.5) 
    scatter(1,mean(betas(:,r)),5, 'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'LineWidth',1.5) 
    
    % Axes and Ticks
    xlim([.999999 1.000001])
    line(xlim,[0 0],'LineWidth',2/3,'Color','k')
    min_b = min(betas(:,r));
    max_b = max(betas(:,r));
    range_b = max_b - min_b;
    margin = range_b/20;
    bf_topup = 3*margin;
    ylim([(min_b - margin) (max_b + margin + bf_topup )])
    ylim_min = (ceil(min_b*10)/10);
    ylim_max = (floor(max_b*10)/10);
    if ylim_min < 0 && ylim_max > 0
        ylabs = [ylim_min 0 ylim_max];
    else
        ylabs = [ylim_min ylim_max];
    end
    set(gca,'YTick',ylabs)
    set(gca,'XTick',[])
    ylabs_str = cellstr(num2str(ylabs'));
    ylabs_str = cellfun(@strrep,ylabs_str, repmat({'0.'},numel(ylabs_str),1), repmat({'.'},numel(ylabs_str),1),'UniformOutput',0);
    set(gca,'YTickLabel',{})
    xlims = xlim;
    text(repmat(xlims(1)-0.0000002,numel(ylabs_str),1),ylabs,ylabs_str,'HorizontalAlignment','right','FontName','Calibri','FontSize',8)
    set(gca,'xcolor',[1 1 1])
    
    % Figure size
    set(gcf,'Units','centimeters');
    x0 = 1;
    y0 = 1;
    width = .4;
    height = 2;
    set(gca,'Units','centimeters','position', [x0,y0,width,height])
    set(gcf,'Units','centimeters','position',[0,0,10.5,8])
    
    title(ROI_names{r})
    
    saveas(f,fullfile(trg_dir, 'model betas', [name,'_',ROIs{r,1},'.emf']))
end


%% Test unimodality and normality

nboot = 1000;
dip = nan(nROI,1);
p_dip = nan(nROI,1);
lillie = nan(nROI,1);
p_lillie = nan(nROI,1);

for r = 1:nROI
    
    % Unimodality
    [dip(r), p_dip(r)] = HartigansDipSignifTest(betas(:,r), nboot);
    
    subplot(4,5,r)
    hist(betas(:,r))
    title([ROIs{r,2} ': dip=',num2str(dip(r),3), ', p=',num2str(p_dip(r),3)])
    
    % Normality
    [lillie(r), p_lillie(r)] = lillietest(betas(:,r));
end




