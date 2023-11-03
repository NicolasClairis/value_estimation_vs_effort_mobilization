function [] = binarize_grey_mask( mask_path, proba_threshold )
%binarize_grey_mask( mask_path, proba_threshold ) binarizes grey masks in
%the mask_path folder based on the proba_threshold.
%
% INPUTS
% mask_path: path where mask(s) to binarize is/are located
% 
% proba_threshold: numeric value in percentage (0 to 1) containing the
% threshold on which to binarize the mask (recommended to set to 0.25 at
% least=25% sure it is grey matter)
%

%% check that value is ok
if ~exist('proba_threshold','var') || isempty(proba_threshold)
    proba_threshold = 0.25;
    warning(['proba_threshold set at ',num2str(proba_threshold)]);
else
    if proba_threshold <= 0 || proba_threshold > 1
        error('proba threshold should be included in the ]0, 1] range.');
    end
end

%% file names
if exist('mask_path','var') && ~isempty(mask_path)
    grey_mask_file_nm           = ls([mask_path, filesep, 'mwc1*.nii']);
else
    error('please specify mask_path in input');
end
% grey_mask_file_nm           = ls([mask_path, filesep, 'grey.nii']); % to binarize SPM mask
% check that path is ok
if isempty(grey_mask_file_nm)
    error(['no grey mask found in ',mask_path,', please check your path is correct.']);
end
nMasks_to_binarize          = size(grey_mask_file_nm,1);

%% prepare SPM batch for binarization
matlabbatch = cell(nMasks_to_binarize,1);
for iMask = 1:nMasks_to_binarize
    
    grey_mask_nm = strrep(grey_mask_file_nm(iMask,:),' ',''); % kill potential parasite spaces in the file name
    
    binarized_grey_mask_file_nm = ['b',grey_mask_nm];
    
    matlabbatch{iMask}.spm.util.imcalc.input            = {[mask_path, filesep, grey_mask_nm]};
    matlabbatch{iMask}.spm.util.imcalc.output           = binarized_grey_mask_file_nm;
    matlabbatch{iMask}.spm.util.imcalc.outdir           = {mask_path};
    matlabbatch{iMask}.spm.util.imcalc.expression       = ['i1>',num2str(proba_threshold)]; % keep only voxels with probability of being in the grey matter higher than 25%
    matlabbatch{iMask}.spm.util.imcalc.var              = struct('name', {}, 'value', {});
    matlabbatch{iMask}.spm.util.imcalc.options.dmtx     = 0;
    matlabbatch{iMask}.spm.util.imcalc.options.mask     = 0;
    matlabbatch{iMask}.spm.util.imcalc.options.interp   = 1;
    matlabbatch{iMask}.spm.util.imcalc.options.dtype    = 4;
end

%% perform the batch
spm_jobman('run',matlabbatch)


end

