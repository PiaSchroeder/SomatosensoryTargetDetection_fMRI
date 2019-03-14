function [PF, fig] = fit_logistic(iniT, int, resp, name, plotit)
% Fits a logistic function to the responses (1|0) stored in resp given the
% intensities stored in int and the initial estimate of threshold iniT
% Parameters:
% a = threshold
% b = slope

if nargin < 5
    plotit = 1;
end

PF.name = name;
PF.Intensities = int;
PF.Responses = resp;
PF.iniT = iniT;

useInts = unique(int);
nInts = length(useInts);

%% Compute detection probabilities
PF.prob = zeros(2,nInts);

for i = 1:nInts
    ints = [int(int==useInts(i)); ...
        resp(int==useInts(i))];
    PF.prob(1,i) = useInts(i);
    if isempty(ints(2,~isnan(ints(2,:))))
        PF.prob(2,i) = nan;
    else
        ints(:,isnan(ints(2,:))) = [];
        PF.prob(2,i) = sum(ints(2,:))/length(ints(2,:));
    end
end



%% Fit sigmoidal

% DEFINITON OF LOGISTIC FUNCTION
logistic = @(c,x) (1./(1+exp(-c(2)*(x-c(1)))));

% FIT PSYCHOMETRIC FUNCTION
[PF.fit_logistic, ~, ~, ~, PF.mse_fit_logistic] = nlinfit(PF.prob(1,:), PF.prob(2,:), logistic, [iniT, 1]);

%% Determine Thresholds

% 1%
logistic_01 =  @(x) (1./(1+exp(-PF.fit_logistic(2)*(x-PF.fit_logistic(1))))) - 0.01;
PF.T01 =  fzero(logistic_01,2);

% 99%
logistic_99 =  @(x) (1./(1+exp(-PF.fit_logistic(2)*(x-PF.fit_logistic(1))))) - 0.99;
PF.T99 =  fzero(logistic_99,2);

% 50%
logistic_50 =  @(x) (1./(1+exp(-PF.fit_logistic(2)*(x-PF.fit_logistic(1))))) - 0.5;
PF.T50 =  fzero(logistic_50,2);

%% Plot
if plotit
    range = useInts(1):0.01:useInts(end);
    fit_log(1,:) = range;
    fit_log(2,:) = logistic(PF.fit_logistic,range);
    
    fig = figure;
    set(gca,'fontsize',18);
    hold on
    
    scatter(PF.prob(1,:),PF.prob(2,:),[],[0,0,0], ...
        'MarkerEdgeColor',[0 .5 .5],...
        'MarkerFaceColor',[0 .7 .7],...
        'LineWidth',1.5);
    line(fit_log(1,:),fit_log(2,:),'Color',[.8,0,.8],'LineWidth',2);
    axis([useInts(1) useInts(end) 0 1]);
    
    ylabel('P(Det)');
    xlabel('Intensity (mA)');
    % title(['Psychometric function: ', name]);
    
    fprintf('T50 = %1.2f\nT01 = %1.2f\nT99 = %1.2f\n', PF.T50, PF.T01, PF.T99);
else
    fig = [];
end
