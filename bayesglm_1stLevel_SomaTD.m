function bayesglm_1stLevel_SomaTD(data_dir,ffx_name, runs, sub)
%% Specifies and estimates 1st level GLM with Bayesian parameter estimation to produce PPMs

rng('shuffle');

log_dir = fullfile(data_dir,'Logs');
tgt_dir = fullfile(data_dir, 'BayesGLM', ffx_name);
if ~exist(tgt_dir, 'dir')
    mkdir(tgt_dir)
    disp(['Created directory: ' tgt_dir])
else
    disp(['Directory already exists: ' tgt_dir])
    disp(['Skipping ' sub ' ' ffx_name])
    return
end

filt = '^warf.*\.nii$';

%% Model Specification
% -------------------------------------------------------------------------
% output directory
matlabbatch{1}.stats{1}.fmri_spec.dir = cellstr(tgt_dir);
% timing parameters
matlabbatch{1}.stats{1}.fmri_spec.timing.units     = 'secs';
matlabbatch{1}.stats{1}.fmri_spec.timing.RT        = 2;
matlabbatch{1}.stats{1}.fmri_spec.timing.fmri_t    = 16;
matlabbatch{1}.stats{1}.fmri_spec.timing.fmri_t0   = 8;

%% Start looping through sessions
%% Onsets + Intensity levels
if regexp(ffx_name, '_IntRT$')
    
    for r = 1:length(runs)
        
        source_dir  = fullfile(data_dir, ['run0' num2str(runs(r))]);
        f           = spm_select('FPList',source_dir, filt);     

        % load the log_file
        log         = dir([log_dir '\log_periT_' sub '_' num2str(runs(r)) '*.mat']);
        log         = log.name;
        load(fullfile(log_dir, log));
        
        % Condition specifications        
        filter              = ~isnan(log_periT.behaviour.detection);                            
        onsets              = log_periT.Design(3,filter);                                       
        [~,intensities]     = ismember(log_periT.Design(1,filter),log_periT.Exp.Intensities);   
        intensities         = zscore(intensities);
        RTs                 = log_periT.behaviour.resp_times(filter);
        RTs                 = zscore(RTs);
        
        % Fill model specification
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).scans = cellstr([f, repmat(',1', size(f,1),1)]);
    
        % All onsets
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).name        = 'Onsets';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).onset       = onsets;
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).duration    = 0;
        
        % Parametric modulation
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).tmod = 0;
        
        % 1: Intensity
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(1).name    = 'Intensity';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(1).param   = intensities';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(1).poly    = 1;
        
        % 1: RTs
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(2).name    = 'RTs';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(2).param   = RTs';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(2).poly    = 1;
        
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).orth = 0;
        
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).hpf = 128;
        
        % Covariates
        % CSF and White matter parameters
        PC      = dir(fullfile(source_dir, ['nPCAs_' sub '_run0' num2str(runs(r)) '.mat']));
        PC      = fullfile(source_dir,PC.name);
        load(PC);
        PC      = R;
                
        % Motion Parameters
        del     = ' ';
        mfile   = dir(fullfile(source_dir, 'rp_f*.txt'));
        mfile   = fullfile(source_dir, mfile.name);
        M       = importdata(mfile, del);
        
        % Fill covariates
        cov = [M PC];
        cov_file = fullfile(source_dir, ['Cov_' sub '.txt']);
        if ~exist(cov_file,'file')
            dlmwrite(cov_file, cov, 'delimiter', '\t', 'newline','pc');
        end
        mf = spm_select('FPList', source_dir, '^Cov_.*\.txt$');
               
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).multi_reg = {mf};
    end

    
%% Onsets + Detectability (PF)
elseif regexp(ffx_name, '_PFRT$')
       
    for r = 1:length(runs)
        
        source_dir  = fullfile(data_dir, ['run0' num2str(runs(r))]);
        f           = spm_select('FPList',source_dir, filt);     

        % load the log_file
        log         = dir([log_dir '\log_periT_' sub '_' num2str(runs(r)) '*.mat']);
        log         = log.name;
        load(fullfile(log_dir, log));
        
        % Retrieve psychometric function
        logistic    = @(c,x) (1./(1+exp(-c(2)*(x-c(1)))));
        fitted      = logistic(log_periT.behaviour.PF.fit_logistic, log_periT.Exp.Intensities);
        
        % Condition specifications
        filter              = ~isnan(log_periT.behaviour.detection);  
        onsets              = log_periT.Design(3,filter);
        [~,intensities]     = ismember(log_periT.Design(1,filter),log_periT.Exp.Intensities);   
        PF                  = fitted(intensities);
        PF                  = zscore(PF);
        RTs                 = log_periT.behaviour.resp_times(filter);
        RTs                 = zscore(RTs);
        
        % Fill model specification
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).scans   = cellstr([f, repmat(',1', size(f,1),1)]);
    
        % All onsets
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).name        = 'Onsets';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).onset       = onsets;
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).duration    = 0;
        
        % Parametric modulation
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).tmod = 0;
        
        % 1: Intensity
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(1).name    = 'PF';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(1).param   = PF';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(1).poly    = 1;
        
        % 1: RTs
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(2).name    = 'RTs';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(2).param   = RTs';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(2).poly    = 1;
        
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).orth = 0;
        
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).hpf = 128;
        
        % Covariates
        % CSF and White matter parameters
        PC      = dir(fullfile(source_dir, ['nPCAs_' sub '_run0' num2str(runs(r)) '.mat']));
        PC      = fullfile(source_dir,PC.name);
        load(PC);
        PC      = R;
                
        % Motion Parameters
        del     = ' ';
        mfile   = dir(fullfile(source_dir, 'rp_f*.txt'));
        mfile   = fullfile(source_dir, mfile.name);
        M       = importdata(mfile, del);
        
        % Fill covariates
        cov = [M PC];
        cov_file = fullfile(source_dir, ['Cov_' sub '.txt']);
        if ~exist(cov_file,'file')
            dlmwrite(cov_file, cov, 'delimiter', '\t', 'newline','pc');
        end
        mf = spm_select('FPList', source_dir, '^Cov_.*\.txt$');
               
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).multi_reg = {mf};
    end    
    
%% Onsets + Detection 
elseif regexp(ffx_name, '_DetRT$')
       
    for r = 1:length(runs)
        
        source_dir  = fullfile(data_dir, ['run0' num2str(runs(r))]);
        f           = spm_select('FPList',source_dir, filt);      

        % load the log_file
        log         = dir([log_dir '\log_periT_' sub '_' num2str(runs(r)) '*.mat']);
        log         = log.name;
        load(fullfile(log_dir, log));
        
        % Condition specification
        filter      = ~isnan(log_periT.behaviour.detection); 
        onsets      = log_periT.Design(3,filter);
        detection   = log_periT.behaviour.detection(filter);
        detection   = zscore(detection);
        RTs         = log_periT.behaviour.resp_times(filter);
        RTs         = zscore(RTs);
        
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).scans   = cellstr([f, repmat(',1', size(f,1),1)]);
        
        % All onsets
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).name        = 'Onsets';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).onset       = onsets;
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).duration    = 0;
        
        % Parametric modulation
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).tmod = 0;
        
        % 1: Detection
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(1).name   = 'Detection';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(1).param  = detection';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(1).poly   = 1;
        
        % 1: RTs
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(2).name    = 'RTs';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(2).param   = RTs';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(2).poly    = 1;
        
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).orth = 0;
        
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).hpf = 128;
        
        % Covariates
        % CSF and White matter parameters
        PC      = dir(fullfile(source_dir, ['nPCAs_' sub '_run0' num2str(runs(r)) '.mat']));
        PC      = fullfile(source_dir,PC.name);
        load(PC);
        PC      = R;
                
        % Motion Parameters
        del     = ' ';
        mfile   = dir(fullfile(source_dir, 'rp_f*.txt'));
        mfile   = fullfile(source_dir, mfile.name);
        M       = importdata(mfile, del);
        
        % Fill covariates
        cov = [M PC];
        cov_file = fullfile(source_dir, ['Cov_' sub '.txt']);
        if ~exist(cov_file,'file')
            dlmwrite(cov_file, cov, 'delimiter', '\t', 'newline','pc');
        end
        mf = spm_select('FPList', source_dir, '^Cov_.*\.txt$');
               
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).multi_reg = {mf};
    end
        
%% Onsets + Uncertainty
% Uncertainty defined as the slope of individual PFs
elseif regexp(ffx_name, 'UncRT$')
    
    for r = 1:length(runs)
        
        source_dir  = fullfile(data_dir, ['run0' num2str(runs(r))]);
        f           = spm_select('FPList',source_dir, filt);     

        % load the log_file
        log         = dir([log_dir '\log_periT_' sub '_' num2str(runs(r)) '*.mat']);
        log         = log.name;
        load(fullfile(log_dir, log));
        
        % Get uncertainty from PF: derivative
        diff_log = @(c,x) (c(2)*exp(c(2)*(x - c(1))))./(exp(c(2)*(x - c(1))) + 1).^2;               % Derivative of logistic function
        diff_fitted = diff_log(log_periT.behaviour.PF.fit_logistic, log_periT.Exp.Intensities);     % Use parameters from fitted logistic
        
        % Condition specification
        filter              = ~isnan(log_periT.behaviour.detection); 
        onsets              = log_periT.Design(3,filter);
        [~,intensities]     = ismember(log_periT.Design(1,filter),log_periT.Exp.Intensities);  
        uncertainty         = diff_fitted(intensities);
        uncertainty         = zscore(uncertainty);
        RTs                 = log_periT.behaviour.resp_times(filter);
        RTs                 = zscore(RTs);
        
        % Fill model specification
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).scans   = cellstr([f, repmat(',1', size(f,1),1)]);
    
        % All onsets
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).name        = 'Onsets';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).onset       = onsets;
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).duration    = 0;
        
        % Parametric modulation
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).tmod = 0;
        
        % 1: Uncertainty
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(1).name    = 'Uncertainty';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(1).param   = uncertainty';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(1).poly    = 1;
        
        % 1: RTs
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(2).name    = 'RTs';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(2).param   = RTs';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(2).poly    = 1;
        
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).orth = 0;
        
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).hpf = 128;
        
        % Covariates
        % CSF and White matter parameters
        PC      = dir(fullfile(source_dir, ['nPCAs_' sub '_run0' num2str(runs(r)) '.mat']));
        PC      = fullfile(source_dir,PC.name);
        load(PC);
        PC      = R;
                
        % Motion Parameters
        del     = ' ';
        mfile   = dir(fullfile(source_dir, 'rp_f*.txt'));
        mfile   = fullfile(source_dir, mfile.name);
        M       = importdata(mfile, del);
        
        % Fill covariates
        cov = [M PC];
        cov_file = fullfile(source_dir, ['Cov_' sub '.txt']);
        if ~exist(cov_file,'file')
            dlmwrite(cov_file, cov, 'delimiter', '\t', 'newline','pc');
        end
        mf = spm_select('FPList', source_dir, '^Cov_.*\.txt$');
               
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).multi_reg = {mf};
    end

%% Report: Match vs Mismatch
elseif regexp(ffx_name, '_RepRT$')
       
    for r = 1:length(runs)
        
        source_dir  = fullfile(data_dir, ['run0' num2str(runs(r))]);
        f           = spm_select('FPList',source_dir, filt);      

        % load the log file
        log         = dir([log_dir '\log_periT_' sub '_' num2str(runs(r)) '*.mat']);
        log         = log.name;
        load(fullfile(log_dir, log));
        
        % Condition specification
        filter      = ~isnan(log_periT.behaviour.detection);  
        onsets      = log_periT.Design(3,filter);
        report      = log_periT.behaviour.match(filter);
        report      = zscore(report);
        RTs         = log_periT.behaviour.resp_times(filter);
        RTs         = zscore(RTs);
        
        % Fill model specification
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).scans   = cellstr([f, repmat(',1', size(f,1),1)]);
        
        % All onsets
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).name        = 'Onsets';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).onset       = onsets;
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).duration    = 0;
        
        % Parametric modulation
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).tmod = 0;
        
        % 1: Report
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(1).name   = 'Report';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(1).param  = report';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(1).poly   = 1;
                
        % 1: RTs
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(2).name   = 'RTs';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(2).param  = RTs';
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).pmod(2).poly   = 1;
        
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(1).orth = 0;
        
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).hpf = 128;
        
        % Covariates
        % CSF and White matter parameters
        PC      = dir(fullfile(source_dir, ['nPCAs_' sub '_run0' num2str(runs(r)) '.mat']));
        PC      = fullfile(source_dir,PC.name);
        load(PC);
        PC      = R;
                
        % Motion Parameters
        del     = ' ';
        mfile   = dir(fullfile(source_dir, 'rp_f*.txt'));
        mfile   = fullfile(source_dir, mfile.name);
        M       = importdata(mfile, del);
        
        % Fill covariates
        cov = [M PC];
        cov_file = fullfile(source_dir, ['Cov_' sub '.txt']);
        if ~exist(cov_file,'file')
            dlmwrite(cov_file, cov, 'delimiter', '\t', 'newline','pc');
        end
        mf = spm_select('FPList', source_dir, '^Cov_.*\.txt$');
               
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).multi_reg = {mf};
        
    end
    
%% 10 conditions: 1 per intensity, modulated by RT    
elseif regexp(ffx_name, '_10RT$')
    
    for r = 1:length(runs)
        
        source_dir  = fullfile(data_dir, ['run0' num2str(runs(r))]);
        f           = spm_select('FPList',source_dir, filt);     

        % load the log_file
        log         = dir([log_dir '\log_periT_' sub '_' num2str(runs(r)) '*.mat']);
        log         = log.name;
        load(fullfile(log_dir, log));
        
        % Condition specification: one condition per intensity 
        for i = 1:10
            filter = ~isnan(log_periT.behaviour.detection) & ...
                log_periT.Design(1,:) == log_periT.Exp.Intensities(i);
            
            onsets = log_periT.Design(3,filter);
            RTs    = log_periT.behaviour.resp_times(filter);
            RTs    = zscore(RTs);
            
            matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(i).name        = ['Int' num2str(i)];
            matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(i).onset       = onsets';
            matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(i).duration    = 0;
            
            % Parametric modulation
            matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(i).tmod = 0;
            matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(i).pmod(1).name    = 'RTs';
            matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(i).pmod(1).param   = RTs';
            matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(i).pmod(1).poly    = 1;
            matlabbatch{1}.stats{1}.fmri_spec.sess(r).cond(i).orth = 0;
            
            clear onsets RTs
        end
        
        % Fill model specification
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).scans = cellstr([f, repmat(',1', size(f,1),1)]);
        
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).hpf = 128;
        
        % Covariates
        % CSF and White matter parameters
        PC      = dir(fullfile(source_dir, ['nPCAs_' sub '_run0' num2str(runs(r)) '.mat']));
        PC      = fullfile(source_dir,PC.name);
        load(PC);
        PC      = R;
                
        % Motion Parameters
        del     = ' ';
        mfile   = dir(fullfile(source_dir, 'rp_f*.txt'));
        mfile   = fullfile(source_dir, mfile.name);
        M       = importdata(mfile, del);
        
        % Fill covariates
        cov = [M PC];
        cov_file = fullfile(source_dir, ['Cov_' sub '.txt']);
        if ~exist(cov_file,'file')
            dlmwrite(cov_file, cov, 'delimiter', '\t', 'newline','pc');
        end
        mf = spm_select('FPList', source_dir, '^Cov_.*\.txt$');
               
        matlabbatch{1}.stats{1}.fmri_spec.sess(r).multi_reg = {mf};
    end
end

%% additional model parameters
matlabbatch{1}.stats{1}.fmri_spec.fact = struct('name', {}, 'levels', {}); 
matlabbatch{1}.stats{1}.fmri_spec.bases.hrf.derivs     = [1 0];            % model hrf and first temporal derivative
matlabbatch{1}.stats{1}.fmri_spec.volt                 = 1;                
matlabbatch{1}.stats{1}.fmri_spec.global               = 'None';           
matlabbatch{1}.stats{1}.fmri_spec.mask                 = {''};             
matlabbatch{1}.stats{1}.fmri_spec.cvi                  = 'AR(1)';          

% create the model
fprintf('Creating BayesGLM\n')
spm_jobman('run', matlabbatch);

clear matlabbatch


%%  Bayesian Model Estimation
% -------------------------------------------------------------------------
spm_file = fullfile(tgt_dir, 'SPM.mat');

matlabbatch{1}.spm.stats.fmri_est.spmmat = {spm_file};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.space.volume.block_type = 'Slices';
matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.signal = 'UGL';
matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.ARP = 3;
matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.noise.UGL = 1;
matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.LogEv = 'Yes';
matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.anova.first = 'No';
matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.anova.second = 'No';

fprintf('Estimating BayesGLM: %s \n', sub);
spm_jobman('run', matlabbatch);

clear SPM;


