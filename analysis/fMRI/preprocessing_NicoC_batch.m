% preprocessing NicoC batch
%% preprocessing for subjects with runs all with with multiband
% enter subject identification in 'subject_id' (sXX_ddMMyy dd: day,
% MM: month, yy: year)
%
% all data is stored in the fMRI_scans folders + the script also creates
% fMRI_analysis folder where the anatomical files are copied for the
% preprocessing + functional folder created for 1st level analysis
%
% See also First_level_NicoC_batch, contrasts_NicoC_batch and
% Second_level_NicoC_batch

clear;

% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

% define subjects and folder where to work
subject_id = {'DEFINE LIST OF SUBJECTS'}; % actual subjects
root = 'path where data is stored to enter here';
NS = length(subject_id); % nber of subjects
scripts_folder = 'enter here path where analysis scripts are stored';

for iSubject = 1:NS % loop through subjects
    
    % create working directories and copy anat. file inside
    % \fMRI_analysis\anatomical folder
    cd([root,subject_id{iSubject}]);
    if exist('fMRI_analysis','dir') ~= 7
        mkdir fMRI_analysis;
    end
    subj_analysis_folder = [root,subject_id{iSubject},'\fMRI_analysis'];
    cd(subj_analysis_folder);
    if exist('anatomical','dir') ~= 7
        mkdir anatomical;
    end
    if exist('functional','dir') ~= 7
        mkdir functional;
    end
    subj_scans_folder = [root,subject_id{iSubject},'\fMRI_scans'];
    cd(subj_scans_folder);
    anat_folder = ls('*t1*');
    copyfile(anat_folder,[subj_analysis_folder,'\anatomical\']); % copies the contempt of the anatomical folder
    
    %% extract files
    cd(subj_scans_folder);
    subj_scan_folders_names = ls('*TR1100_2iso_PA_RUN*'); % takes all functional runs folders (if TR = 1.10s, for multiband seq in particular)
    %%
    cd(subj_analysis_folder)
    nb_runs = 7; % 7 runs for second pilot subject and for main ones (2 grip, 2 mental, 3 learning)
    nb_preprocessing_steps = 6;
 % nb of preprocessing steps per subject: 1)realign-reslice, 2)coregistration functional - anatomical, 3) segmentation, 4) normalization functional, 5) normalization anatomical, 6) smoothing functional
    cd(subj_scans_folder);
    for run = 1:nb_runs % loop through runs for 3 ratings, 3 choices 1D, 3 choices 2D runs
        cd(subj_scan_folders_names(run,:)); % go to run folder
        filenames = cellstr(spm_select('ExtFPList',pwd,'^f.*\.img$')); % extracts all the f-files
        if run == 1
            run1 = filenames;
        elseif run == 2
            run2 = filenames;
        elseif run == 3
            run3 = filenames;
        elseif run == 4
            run4 = filenames;
        elseif run == 5
            run5 = filenames;
        elseif run == 6
            run6 = filenames;
        elseif run == 7
            run7 = filenames;
        elseif run == 8
            run8 = filenames;
        elseif run == 9
            run9 = filenames;
        end
        cd(subj_scans_folder);
    end
    % create file for storing all preprocessed files
    
    %% realignement
    cd(subj_scans_folder);
    preproc_step = 1; % first step of preprocessing
    realign_step = nb_preprocessing_steps*(iSubject-1) + preproc_step; % number depends also on the number of subjects
    if nb_runs == 7
        matlabbatch{realign_step}.spm.spatial.realign.estwrite.data = {run1; run2; run3; run4; run5; run6; run7}';
    else
        warning('problem with nb_runs');
        return;
    end
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    %% coregistration
    preproc_step = 2;
    coreg_step = nb_preprocessing_steps*(iSubject-1) + preproc_step;
    matlabbatch{coreg_step}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{realign_step}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
    cd([root,subject_id{iSubject},'\fMRI_analysis\anatomical\']);
    anat_file = ls(['s*.img']);
    matlabbatch{coreg_step}.spm.spatial.coreg.estimate.source = {[root,subject_id{iSubject},'\fMRI_analysis\anatomical\',anat_file]};
    cd(subj_scans_folder);
    matlabbatch{coreg_step}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{coreg_step}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{coreg_step}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{coreg_step}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{coreg_step}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    %% segmentation
    cd(subj_scans_folder);
    preproc_step = 3;
    segm_step = nb_preprocessing_steps*(iSubject-1) + preproc_step;
    matlabbatch{segm_step}.spm.spatial.preproc.channel.vols = {[root,subject_id{iSubject},'\fMRI_analysis\anatomical\',anat_file]};
    matlabbatch{segm_step}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{segm_step}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{segm_step}.spm.spatial.preproc.channel.write = [1 1];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(1).tpm = {'C:\Users\nicolas.clairis\Desktop\spm12\tpm\TPM.nii,1'};
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(1).warped = [1 1]; % normalise grey matter and record modulated and unmodulated warped tissue
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(2).tpm = {'C:\Users\nicolas.clairis\Desktop\spm12\tpm\TPM.nii,2'};
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(3).tpm = {'C:\Users\nicolas.clairis\Desktop\spm12\tpm\TPM.nii,3'};
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(4).tpm = {'C:\Users\nicolas.clairis\Desktop\spm12\tpm\TPM.nii,4'};
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(5).tpm = {'C:\Users\nicolas.clairis\Desktop\spm12\tpm\TPM.nii,5'};
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(6).tpm = {'C:\Users\nicolas.clairis\Desktop\spm12\tpm\TPM.nii,6'};
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{segm_step}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{segm_step}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{segm_step}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{segm_step}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{segm_step}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{segm_step}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{segm_step}.spm.spatial.preproc.warp.write = [1 1];
    %% normalization
    preproc_step = 4;
    normf_step = nb_preprocessing_steps*(iSubject-1) + preproc_step;
    matlabbatch{normf_step}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{segm_step}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    for run = 1:nb_runs
        matlabbatch{normf_step}.spm.spatial.normalise.write.subj.resample(run) = cfg_dep(['Realign: Estimate & Reslice: Resliced Images (Sess ',num2str(run),')'], substruct('.','val', '{}',{realign_step}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{run}, '.','rfiles'));
    end
    matlabbatch{normf_step}.spm.spatial.normalise.write.subj.resample(nb_runs+1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{realign_step}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
    matlabbatch{normf_step}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{normf_step}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{normf_step}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{normf_step}.spm.spatial.normalise.write.woptions.prefix = 'w';
    %% normalization anatomical scan
    preproc_step = 5;
    norma_step = nb_preprocessing_steps*(iSubject-1) + preproc_step;
    matlabbatch{norma_step}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{segm_step}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    matlabbatch{norma_step}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{segm_step}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
    matlabbatch{norma_step}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{norma_step}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{norma_step}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{norma_step}.spm.spatial.normalise.write.woptions.prefix = 'w';
    %% smoothing
    preproc_step = 6;
    smooth_step = nb_preprocessing_steps*(iSubject-1) + preproc_step;
    matlabbatch{smooth_step}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{normf_step}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{smooth_step}.spm.spatial.smooth.fwhm = [8 8 8];
    matlabbatch{smooth_step}.spm.spatial.smooth.dtype = 0;
    matlabbatch{smooth_step}.spm.spatial.smooth.im = 0;
    matlabbatch{smooth_step}.spm.spatial.smooth.prefix = 's';

    cd(root);
    
end

cd(scripts_folder);
% display spm batch before running it
spm_jobman('interactive',matlabbatch);
% run the batch without checking first
% spm_jobman('run',mtatlabbatch);