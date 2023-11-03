function [] = contrasts_megaconcatenation_MS2_NicoC(GLM, checking)
% contrasts_megaconcatenation_NicoC_batch(GLM, checking) 
% concatenates together all tasks in the same GLM
% For Motiscan1 study for april-may 2016 (multisequence) or june 2016 (EPI
% only)
%
% INPUTS
% GLM: GLM number to be used (check which_GLM_MS2.m for more details about each
% GLM specifications)
%
% checking: if 0 or empty, run the GLM with all the subjects, if equal to
% 1, displays the interactive window on only one subject
%
% See also which_GLM_MS2.m, First_level_MS2_megaconcatenation_NicoC_batch

close all; clc;

%% working directories
root = 'enter path here';
scripts_folder = 'enter path here';

%% by default checking = 0 if not selected
if ~exist('checking','var') || isempty(checking)
    checking = 0;
end

%% GLM parameters
if ~exist('GLM','var') || isempty(GLM)
    GLM = input(sprintf('GLM number? (check which_GLM_MS2.m if doubt on which one to use) \n '));
end
[GLMprm] = which_GLM_MS2(GLM);

%% subjects identification
subject_id = {'enter list of subjects here'};
NS = length(subject_id);

%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% loop through subjects to extract all the regressors
matlabbatch = cell(NS,1);
batch_idx = 0;
for iSub = 1:NS
    
    %% extract current subject name
    sub_nm = subject_id{iSub};
    
    batch_idx = batch_idx + 1;
    
    run_foldername = ['GLM',num2str(GLM) filesep];
    
    matlabbatch{batch_idx}.spm.stats.con.spmmat = {fullfile(root,sub_nm,'fMRI_analysis','functional',run_foldername,'SPM.mat')};
    
    %% extract contrasts list (vectors + corresponding names
    [con_names, con_vec] = MS2_load_con(GLMprm, sub_nm);
    n_con = length(con_names);
    
    %% add each contrast to the list
    for iCon = 1:n_con
        matlabbatch{batch_idx}.spm.stats.con.consess{iCon}.tcon.name     = con_names{iCon};
        matlabbatch{batch_idx}.spm.stats.con.consess{iCon}.tcon.weights  = con_vec(iCon,:);
        matlabbatch{batch_idx}.spm.stats.con.consess{iCon}.tcon.sessrep  = 'none';
    end
    
    matlabbatch{batch_idx}.spm.stats.con.delete = 1; % deletes previous contrasts
    
end % subject loop

cd(scripts_folder);

%% display spm batch before running it or run it directly
switch checking
    case 0
        spm_jobman('run',matlabbatch);
    case 1
        spm_jobman('interactive',matlabbatch);
end

end % function end