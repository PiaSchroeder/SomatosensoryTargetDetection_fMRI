function VOI_max = find_bms_peaks(map, maps_dir, mask, mask_dir, SJs)
%% Find  BMS peak voxels in ROIs. 
%
% INPUTS:
%
% map:      BMS mask from single subject BMS winning model
% map_dir:  BMS model mask directory
% mask:     Name of group mask that should be used to search peaks (in LogEv
%           space!!)
% mask_dir: Group mask directory
% SJs:      Participants to be included 
%

%%
map = [map '_mask.nii'];
data_dir = 'my_data_dir';
disp(mask)

%% Read in VOI coordinates
mask_hdr        = spm_vol(fullfile(mask_dir, [mask '.nii']));
[mask_Y,~]      = spm_read_vols(mask_hdr);
mask_Y          = reshape(mask_Y,[1 prod(mask_hdr.dim)]);
mask_idx        = find(mask_Y);

VOI_max.mask    = mask;
VOI_max.map     = map;
VOI_max.coords  = cell(3,numel(SJs));

%% Run through participants
for s = 1:length(SJs)
    
    disp(['Extracting peaks: ' SJs{s}])
    
    bmsmap          = [map(1:end-8) 'model_ppm.nii'];
    bmsmap_dir      = fullfile(data_dir, SJs{s}, maps_dir, bmsmap);
    hdr             = spm_vol(bmsmap_dir); 
    [Y,XYZ]         = spm_read_vols(hdr);  
    Y               = reshape(Y,[1 prod(hdr.dim)]);
    Y_masked        = Y(mask_idx);
    XYZ_masked      = XYZ(:,mask_idx);
    [max_Y,max_ind] = max(Y_masked); 
    xyz             = XYZ_masked(:,max_ind);
    VOI_max.map     = bmsmap;
    
    VOI_max.coords(1,s) = SJs(s);   % subject ID
    VOI_max.coords(2,s) = {xyz'};   % peak coordinates
    VOI_max.coords(3,s) = {max_Y};  % posterior model probability
    
end
