function bms_SomaTD(data_dir, name, SJs, models, level)

name = ['BMS_' level '_' name];
filt = 'srLogEv.nii'; % smoothed and resampled log evidence maps

trg_dir = fullfile(data_dir, '2nd level', 'BayesGLM', name);
if ~exist(trg_dir, 'dir')
    mkdir(trg_dir)
end

% set SPM defaults
spm('defaults','fmri')
spm_jobman('initcfg');
spm_get_defaults('cmdline',true)

matlabbatch{1}.spm.stats.bms_map.inference.dir = {trg_dir};

for s = 1:numel(SJs)
    f = cell(numel(models),1);
    % Assemble logEv maps
    for m = 1:numel(models)
        logEv_dir   = fullfile(data_dir, SJs{s}, 'BayesGLM', ['BayesFFX' models{m}]);
        f{m,1}      = [spm_select('FPList', logEv_dir, filt), ',1'];
        if length(f{m,1}) < 3
            error(['LogEv map not found: ' SJs{s} ' ' models{m}]);
        end
    end
    
    matlabbatch{1}.spm.stats.bms_map.inference.sess_map{s}.mod_map = f;
end

mask_file = 'mydir\mask_ICV.nii';
matlabbatch{1}.spm.stats.bms_map.inference.mod_name = models';
matlabbatch{1}.spm.stats.bms_map.inference.method_maps = level;
matlabbatch{1}.spm.stats.bms_map.inference.out_file = 1;
matlabbatch{1}.spm.stats.bms_map.inference.mask = {mask_file};
matlabbatch{1}.spm.stats.bms_map.inference.nsamp = '1e6';

% Run model selection
fprintf('Running BMS\n')
spm_jobman('run', matlabbatch);
clear matlabbatch


