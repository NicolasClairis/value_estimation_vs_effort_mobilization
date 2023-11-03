function[] = First_level_MS2_megaconcatenation_NicoC_batch(GLM, checking,...
    bayesian_mdl)
%[] = First_level_MS2_megaconcatenation_NicoC_batch(GLM, checking, subject_id, pc_cluster, bayesian_mdl)
% Batch script for first level fMRI for multiband (1.1s) fMRI sequence of
% MotiScan-2
% Use of SPM-12
% define subjects manually, ask the GLM you want to use and then displays
% it into SPM GUI ("interactive") or runs it directly ("run") depending on
% the parameter entered in the last line
%
% INPUTS
% GLM: GLM identification number (see which_GLM_MS2.m for details)
%
% checking
% (0) directly runs the first level GLM
% (1) only displays interactive window on one subject to see if everything
% is fine
%
% bayesian_mdl
% (0) use classical First level model (by default if not entered in the inputs)
% (1) use bayesian model
%
% Makes a global GLM where all three tasks of MS2 (learning, grip and stroop) are pooled together
%
% See also: preprocessing_NicoC_batch, onsets_for_fMRI_MS2,
% group_onsets_for_fMRI_MS2, which_GLM_MS2

%% working directories
root = 'define path here';
scripts_folder = 'path here';

%% add script folder to the path
addpath(scripts_folder);

%% by default checking = 0 if not selected
if ~exist('checking','var') || isempty(checking)
    checking = 0;
end

%% classical or bayesian model (only if you need to compare models with BMS)
if ~exist('bayesian_mdl','var') || isempty(bayesian_mdl)
    bayesian_mdl = 0; % by default classical model
end

%% GLM
if ~exist('GLM','var') || isempty(GLM)
    GLM = input(sprintf('GLM number? (check which_GLM_MS2.m if doubt on which one to use) \n '));
end
% load specific GLM parameters
[GLMprm] = which_GLM_MS2(GLM);
% extract relevant parameters for this script
grey_mask   = GLMprm.gal.grey_mask;
add_drv     = GLMprm.gal.add_drv;
FIR         = GLMprm.gal.FIR;
if FIR == 1
    FIR_dur     = GLMprm.gal.FIR_dur;
    FIR_nBins   = GLMprm.gal.FIR_nBins;
end

%% load subjects
subject_id = {'enter list of subjects here'};
NS = length(subject_id);

%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% main parameters
learningRuns    = 3;
gripRuns        = 2;
stroopRuns      = 2;
nbRuns = learningRuns + gripRuns + stroopRuns;

% particular sequence parameters
TR = 1.1; % TR for multiband

nb_batch_per_subj = 2; % (model + estimate)

matlabbatch = cell(nb_batch_per_subj*NS,1);

for iSub = 1:NS
    
    sub_nm = subject_id{iSub};
    % extract subject number
    if strcmp(sub_nm(3),'_')
        subid   = sub_nm(2);
    elseif strcmp(sub_nm(3),'_') == 0 && strcmp(sub_nm(4),'_')
        subid   = sub_nm(2:3);
    end
    % define working folders
    subj_folder             = [root, filesep, sub_nm];
    subj_analysis_folder    = [subj_folder, filesep, 'fMRI_analysis' filesep];
    subj_anat_folder        = [subj_analysis_folder, filesep, 'anatomical' filesep];
    subj_scans_folder       = [subj_folder, filesep, 'fMRI_scans' filesep];
    subj_behavior_folder    = [subj_folder, filesep, 'behavior' filesep];
    % folder names
    subj_scan_folders_names = ls([subj_scans_folder, filesep, '*TR1100_2iso_PA_RUN*']); % takes all functional runs folders (if TR = 1.10s, for multiband seq in particular)
    
    % create folder for storing data for this subject
    switch bayesian_mdl
        case 0
            filename = [subj_analysis_folder 'functional', filesep,...
                'GLM',num2str(GLM),'_megaconcatenation'];
        case 1
            filename = [subj_analysis_folder 'functional', filesep,...
                'GLM',num2str(GLM),'_bayesian'];
    end
    mkdir(filename);
    
    
    %% starting 1st level GLM batch
    sub_idx = nb_batch_per_subj*(iSub-1) + 1 ;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.dir            = {filename};
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.units   = 'secs';
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.RT      = TR;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.fmri_t  = 16;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
    
    % loop through runs and tasks
    for iRun = 1:nbRuns
        runname = num2str(iRun);
        % erase useless spaces from folder with run name
        n_char = size(subj_scan_folders_names(iRun,:),2);
        for iLetter = 1:n_char
            if strcmp(subj_scan_folders_names(iRun,iLetter),' ') == 0 % erase space
                subj_runFoldername(iLetter) = subj_scan_folders_names(iRun,iLetter);
            end
        end
        
        % load scans in the GLM
        cd([subj_scans_folder filesep subj_runFoldername, filesep]); % go to run folder
        preprocessed_filenames = cellstr(spm_select('ExtFPList',pwd,'^swrf.*\.img$')); % extracts all the preprocessed swrf files (smoothed, normalized, realigned)
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).scans = preprocessed_filenames;
        
        % identify task corresponding to this run
        switch pc_cluster
            case 'pc'
                behaviorFile = ls([subj_behavior_folder, filesep,...
                    'global_sub_',subid,'_session_',runname,'_*']);
            case 'cluster'
                behaviorFile = cell2mat(get_folder_list(['global_sub_',subid,'_session_',runname,'_*'],...
                    subj_behavior_folder));
        end
        taskName = getfield(load([subj_behavior_folder, filesep, behaviorFile], 'taskName'),'taskName'); % load taskName to see what task is corresponding to this run
        
        switch taskName
            %% learning task
            case {'learning'}
                matlabbatch = First_level_MS2_learning_prm(matlabbatch, GLMprm,...
                    subid, sub_idx, iRun, subj_analysis_folder);
                %% grip task
            case {'gripRP'}
                matlabbatch = First_level_MS2_gripRP_prm(matlabbatch, GLMprm,...
                    subid, sub_idx, iRun, subj_analysis_folder);
                %% stroop task
            case {'mentalRP'}
                matlabbatch = First_level_MS2_mentalRP_prm(matlabbatch, GLMprm,...
                    subid, sub_idx, iRun, subj_analysis_folder);
        end % task
        
        %% global run parameters (rp movement file, etc.)
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).multi = {''};
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).regress = struct('name', {}, 'val', {});
        mvmtFolder = [subj_scans_folder, filesep, subj_runFoldername, filesep];
        switch pc_cluster
            case 'pc'
                movement_file = ls([mvmtFolder, 'rp*']);
            case 'cluster'
                movement_file = cell2mat(get_folder_list('rp*',mvmtFolder));
        end
        movement_filePath = [mvmtFolder, movement_file];
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).multi_reg = {movement_filePath};
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).hpf = 128;
        
        % clear folder variable at the end of each run loop to avoid bugs
        clear('subj_runFoldername');
    end % run loop
    
    
    %% global parameters for subject batch
    matlabbatch{sub_idx}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
   
    if FIR == 0
        % add temporal derivative or not
        if add_drv == 0
            matlabbatch{sub_idx}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        elseif add_drv == 1
            matlabbatch{sub_idx}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
        elseif add_drv == 2
            matlabbatch{sub_idx}.spm.stats.fmri_spec.bases.hrf.derivs = [1 1];
        end
    elseif FIR == 1
        % SPM will extract 'FIR_bins' datapoints from the onset until
        % onset+FIR_dur duration
        matlabbatch{sub_idx}.spm.stats.fmri_spec.bases.fir.length   = FIR_dur; % total duration checked (from the onset until onset + FIR_dur)
        matlabbatch{sub_idx}.spm.stats.fmri_spec.bases.fir.order    = FIR_nBins; % number of bins
    end
    
    matlabbatch{sub_idx}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.global = 'None';
    switch grey_mask
        case 0 % filter with threshold here if no mask entered
            matlabbatch{sub_idx}.spm.stats.fmri_spec.mthresh = 0.8; % default value
        case {1,2,3} % no filter here since the mask will do the job
            matlabbatch{sub_idx}.spm.stats.fmri_spec.mthresh = -Inf; % implicitly masks the first level depending on the probability that the voxel is "relevant"
    end
    
    % add grey mask or not
    switch grey_mask
        case 0
            matlabbatch{sub_idx}.spm.stats.fmri_spec.mask = {''};
        case {1,2,3,4}
            
            % find grey matter mask
            switch grey_mask
                case 1 % grey mask per subject
                    switch pc_cluster
                        case 'pc'
                            mask_file = ls([subj_anat_folder filesep 'bmwc1s*']); % modulated grey matter mask
                        case 'cluster' % path stored when using ls with the path
                            mask_file = cell2mat(get_folder_list('bmwc1s*',subj_anat_folder));
                    end
                    mask_file_path = [subj_anat_folder, mask_file];
                case 2 % grey matter filter across subs
                    mask_file_path = [root, filesep, 'Second_level', filesep,...
                        'mean_anatomy', filesep,'bsmean_greyMatter_22_subjects.nii'];
                case 3 % SPM template
                    mask_file_path = [root, filesep, 'Second_level', filesep,...
                        'mean_anatomy', filesep,'10percent', filesep,...
                        'bgrey_10perc_SPM_template.nii'];
                case 4 % SPM template bis (305 subjects, 10%)
                    mask_file_path = [root, filesep, 'Second_level', filesep,...
                        'mean_anatomy',filesep,'10percent',...
                        filesep,'greyMatter_SPM305template_10perc.nii'];
            end
            % load grey mask (or check if file missing)
            if exist(mask_file_path,'file')
                matlabbatch{sub_idx}.spm.stats.fmri_spec.mask = {mask_file_path};
            else
                error('problem: mean anatomy file not found for this number of subjects');
            end
        case 5 % SPM template, exclusive mask
            mask_file_path = [root, filesep, 'Second_level', filesep,...
                'mean_anatomy',filesep,'10percent', filesep,...
                'LCR-exclusive_SPM305template_10perc.nii'];
            % load grey mask (or check if file missing)
            if exist(mask_file_path,'file')
                matlabbatch{sub_idx}.spm.stats.fmri_spec.mask = {mask_file_path};
            else
                error('problem: mean anatomy file not found for this number of subjects');
            end
    end
    
    matlabbatch{sub_idx}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    
    %% estimate model
    estimate_mdl_rtg_idx = nb_batch_per_subj*iSub;
    matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{sub_idx}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.write_residuals = 0;
    switch bayesian_mdl
        case 0
            % classical model
            matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.method.Classical = 1;
        case 1
            % bayesian model (to compare with BMS which one is the best model afterwards)
            matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.method.Bayesian.space.volume.block_type = 'Slices';
            matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.method.Bayesian.signal = 'UGL';
            matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.method.Bayesian.ARP = 3;
            matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.method.Bayesian.noise.UGL = 1;
            matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.method.Bayesian.LogEv = 'Yes'; % crucial part = write log model evidence (per voxel) to compute model comparison afterwards
            matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.method.Bayesian.anova.first = 'No';
            matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.method.Bayesian.anova.second = 'No';
            matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.method.Bayesian.gcon = struct('name', {}, 'convec', {});
    end
end % subject loop

cd(scripts_folder)

%% display spm batch before running it or run it directly
switch checking
    case 0
        spm_jobman('run',matlabbatch);
    case 1
        spm_jobman('interactive',matlabbatch);
end

end % function