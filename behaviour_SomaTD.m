%% SomaTD behavioural analysis 
% PFs --> data exclusion
% detection rates
% RTs
% Associations between detection and report
%

clear all
close all
clc

%%
% =========================================================================
% 0. Load data
% =========================================================================

SJs         = { 'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S08' 'S09' 'S10' 'S11' 'S12' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S23' 'S24' 'S25' 'S26' 'S27' 'S28' 'S29' 'S30' 'S31' 'S32'};

data_dir = 'my_data_dir';   % Data directory
ana_dir = 'my_ana_dir';     % Analysis directory
trg_dir = 'my_trg_dir';     % Target directory
log_dir = 'my_log_dir';     % Log file directory

nSubs = length(SJs);
nRuns = 4;
nInts = 10;
nTrials = 110;

runs = cell(1,nSubs);
runs(:) = {1:nRuns};

cd(ana_dir)
[Data, PFs] = load_logs(SJs,runs,data_dir,log_dir);

%%
% =========================================================================
% 1. Psychometric functions
% =========================================================================

logistic = @(c,x) (1./(1+exp(-c(2)*(x-c(1)))));

% 1.1 Normalise PFs to intensity range 1-10
iniT = 5.5;
ints = cell(nSubs,nRuns);
resp = cell(nSubs,nRuns);
norm_PFs = cell(nSubs,nRuns);
norm_slopes = nan(nSubs,nRuns);
norm_threshs = nan(nSubs,nRuns);

for s = 1:nSubs
    for r = 1:nRuns
        [~,ints{s,r}] = ismember(Data{s,r}.behaviour.PF.Intensities,Data{s,r}.Exp.Intensities); % Get intensity levels
        resp{s,r} = Data{s,r}.behaviour.PF.Responses;                                           % Get responses 
        norm_PFs{s,r} = fit_logistic(iniT,ints{s,r},resp{s,r},SJs{s},0);                        % Fit logistic function --> normalised PF
        norm_slopes(s,r) = norm_PFs{s,r}.fit_logistic(1,2);                                     % Get normalised slope
        norm_threshs(s,r) = norm_PFs{s,r}.fit_logistic(1,1);                                    % Get normalised T50
    end
end

mean_norm_slopes = mean(norm_slopes,2);
mean_norm_threshs = mean(norm_threshs,2);

normPFs.norm_PFs = norm_PFs;
normPFs.norm_slopes = norm_slopes;
normPFs.norm_threshs = norm_threshs;
normPFs.mean_norm_slopes = mean_norm_slopes;
normPFs.mean_norm_threshs = mean_norm_threshs;

% 1.2 Fit mean PF for every participant
int_range = 1:nInts;
mean_fit = nan(nSubs,nInts);
for s = 1:nSubs
    mean_fit(s,:) = logistic([mean_norm_threshs(s), mean_norm_slopes(s)],int_range);
end

normPFs.mean_fit = mean_fit;

% 1.3 Apply exclusion criterion: mean fitted detection probability min>10% or max<90%
excl = mean_fit(:,1) > .1 | mean_fit(:,end) < .9;

% 1.4 Plot PFs
plot_range = 0:0.01:nInts;

figure 
hold on
set(gca,'FontSize',22)

% Exclusion areas
area([0 10],[1 1],'FaceColor',[1 1 1], 'LineStyle', 'none')
area([0 10],[.9 .9],'FaceColor',[.7 .7 .7], 'LineStyle', 'none')
area([0 10],[.1 .1],'FaceColor',[1 1 1], 'LineStyle', 'none')

for s = 1:nSubs
    sub_pf = logistic([mean_norm_threshs(s), mean_norm_slopes(s)],plot_range);
    if excl(s)
        line(plot_range,sub_pf,'Color',[1 0 0],'LineWidth',2,'LineStyle','--');
    else
        line(plot_range,sub_pf,'Color',[0 0 0],'LineWidth',2);
    end
end
set(gca, 'XTick', 1:10)
set(gca, 'YTick', [0 .5 1])
axis([1 10 0 1]);
ylabel('Detection probability');
xlabel('Intensity level');

% 1.5 Exclude outliers from further analyses
SJs(excl) = [];
nSub = length(SJs);
runs = cell(1,nSub);
runs(:) = {1:nRuns};
Data(excl,:) = [];
PFs(excl,:) = []; 

%%
% =========================================================================
% 2. Detection rates and threshold intensities
% =========================================================================

% 2.1 Detection rate
det_rates = nan(nSub,1);
for s = 1:nSub    
    run_dr = nan(nRuns,1);
    for r = runs{s}
        run_dr(r) = Data{s,r}.behaviour.det_rate;
    end
    det_rates(s) = mean(run_dr);
end

mean_det_rates = mean(det_rates);

DetRates.det_rates = det_rates;
DetRates.mean_det_rates = mean_det_rates;

% 2.2 Initially estimated threshold intensities (T01, T50, T99)
T50s = nan(nSub,1);
T01s = nan(nSub,1);
T99s = nan(nSub,1);

for s = 1:nSub
    % load log file for threshold run
    sub_dir = fullfile(data_dir,SJs{s},'Logs');
    log = dir(fullfile(sub_dir, ['psychFunc_' SJs{s} '.mat']));
    log = log.name;
    load(fullfile(sub_dir, log));
    
    T50s(s) = PF.T50;
    T01s(s) = PF.T01;
    T99s(s) = PF.T99;
end

mean_T50 = mean(T50s);
mean_T01 = mean(T01s);
mean_T99 = mean(T99s);

Ts.T50s = T50s;
Ts.mean_T50 = mean_T50;
Ts.T01s = T01s;
Ts.mean_T01 = mean_T01;
Ts.T99s = T99s;
Ts.mean_T99 = mean_T99;

%%
% =========================================================================
% 3. Reaction times
% =========================================================================

RT.yes = nan(nSub,1);
RT.no = nan(nSub,1);
for s = 1:nSub
    run_yes = nan(nRuns,1);
    run_no = nan(nRuns,1);
    for r = 1:nRuns
        run_yes(r) = mean(Data{s,r}.behaviour.resp_times(Data{s,r}.behaviour.detection == 1));
        run_no(r) = mean(Data{s,r}.behaviour.resp_times(Data{s,r}.behaviour.detection == 0));
    end
    RT.yes(s) = mean(run_yes);
    RT.no(s) = mean(run_no);
end

RT.mean_yes = mean(RT.yes);
RT.mean_no = mean(RT.no);
RT.std_yes = std(RT.yes);
RT.std_no = std(RT.no);
RT.se_yes = std(RT.yes)/sqrt(nSub);
RT.se_no = std(RT.no)/sqrt(nSub);
RT.diff = RT.yes - RT.no;
RT.mean_diff = mean(RT.diff);
RT.se_diff = std(RT.diff)/sqrt(nSub);

[RT.h_normal,RT.p_normal] = lillietest(RT.diff);
[RT.h_diff,RT.p_diff,RT.ci_diff, RT.stats_diff] = ttest(RT.yes,RT.no);

% Plot group result det vs nodet
% figure
% set(gca,'fontsize',18)
% f = barwitherr([RT.se_yes RT.se_no],[RT.mean_yes RT.mean_no]);
% xlim([0.5 2.5])
% ylim([.33 .38])
% set(f,'FaceColor',[0.7 0.7 0.7])
% set(gca,'XTickLabel',{'Detected', 'Not detected'})
% ylabel('Reaction time (sec)');
% box off

% Plot individual RTs det vs no det
% X = [ones(nSub,1) 2*ones(nSub,1)];
% Y = [RT.yes RT.no];
% 
% figure;
% hold on
% plot(X',Y','-*','MarkerSize',10,'LineWidth',2)
% plot([1;2],mean(Y)','k-o','MarkerSize',10,'LineWidth',4)
% set(gca, 'fontsize',18)
% xlim([0.5 2.5])
% set(gca, 'XTick',[1 2])
% set(gca, 'XTickLabel', {'Detected','Not detected'})
% ylabel('Reaction time [s]')
% title('Individual reaction times for detected and undetected trials')
% box off


%%
% =========================================================================
% 4. Response associations
% =========================================================================

Resp.Det = nan(nSub,nTrials*nRuns);     % Detected/Not detected
Resp.Match = nan(nSub,nTrials*nRuns);   % Match/Mismatch

Det_Match_p = nan(nSub,1);

for s = 1:nSub
    for r = 1:nRuns
        Resp.Det(s,(1+(r-1)*nTrials):(nTrials+(r-1)*nTrials)) = Data{s,r}.behaviour.detection;
        Resp.Match(s,(1+(r-1)*nTrials):(nTrials+(r-1)*nTrials)) = Data{s,r}.behaviour.match; 
    end
    [~,~,Det_Match_p(s)] = crosstab(Resp.Det(s,:),Resp.Match(s,:));
end

Resp.Det_Match_p = Det_Match_p;

%%
% =========================================================================
% Assemble results
% =========================================================================

Behaviour.Data = Data;
Behaviour.PFs = PFs;
Behaviour.normPFs = normPFs;
Behaviour.DetRates = DetRates;
Behaviour.Ts = Ts;
Behaviour.RT = RT;
Behaviour.Resp = Resp;

save(fullfile(trg_dir,'Behaviour.mat'),'Behaviour');
