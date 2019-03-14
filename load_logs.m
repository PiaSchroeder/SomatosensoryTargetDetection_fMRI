function [Data, PF] = load_logs(SJs,runs,data_dir,log_dir)
%% load all runs of all subjects into one cell array *Data*

Data = cell(length(SJs), 4);
PF = cell(length(SJs),1);

for s = 1:length(SJs)
    sub_log_dir = fullfile(data_dir, SJs{s}, log_dir);
    for r = 1:length(runs{s})
        % Retrieve file names
        log = dir(fullfile(sub_log_dir, ['log_periT_' SJs{s} '_' num2str(runs{s}(r)) '*.mat']));
        log = log.name;
        % load log file
        Log = load(fullfile(sub_log_dir, log));
        Data{s,r} = Log.log_periT;
    end
    pf_dir = fullfile(data_dir, SJs{s}, 'Logs');
    pf = dir(fullfile(pf_dir,['psychFunc_' SJs{s} '.mat']));
    pf = pf.name;
    pf = load(fullfile(pf_dir, pf));
    PF{s} = pf.PF;
end