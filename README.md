# SomatosensoryTargetDetection_fMRI
Matlab scripts to analyse data from an fMRI study on somatosensory target detection

Requires SPM12 software package: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/

Model log evidence maps and exceedance probabilities can be downloaded at https://figshare.com/articles/Neural_basis_of_somatosensory_target_detection_-_Data/7347167 and directly inspected at https://neurovault.org/collections/4496/

## Scripts
### Behavioural analysis:
* behaviour_SomaTD.m
* load_logs.m
* fit_logistic.m

### Bayesian 1st level GLMs:
* bayesglm_1stlevel_SomaTD_batch.m
* bayesglm_1stlevel_SomaTD.m

### Bayesian Model Selection:
* bms_SomaTD_batch.m
* bms_SomaTD.m

### Beta extraction:
* find_bms_peaks.m
* beta_extraction_peak_ExpReg.m
* beta_extraction_voi_10.m
