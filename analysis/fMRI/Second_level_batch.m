function[] = Second_level_batch(GLM)
%Second_level_batch(GLM)
% Second_level_batch: concatenate MBB EPI subjects different
% tasks together
% based on first levels with all tasks (ratings/choice1D/choiceRE) pooled together
%
% INPUTS
% GLM: GLM number
% See also: First_level_batch, contrasts_MS2_batch

close all; clc;

% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% define study and subjects
% which GLM?
if ~exist('GLM','var') || isempty(GLM) || GLM <= 0
    GLM = spm_input('GLM number?',1,'e',[]);
end
GLM_str = num2str(GLM);
% separate per sex category?
% subSplit = spm_input('separate per sex ?',1,'b','All | Males | Females',[0 1 2]);

% load GLM parameters
[GLMprm] = which_GLM_MS2(GLM);
% define main folder
root = 'enter path here';
% define subjects
subject_id = {'enter list of subjects here'};
NS = length(subject_id);

% load contrasts (names + order)
[list_con, ~] = MS2_load_con(GLMprm, subject_id{1});
n_con = length(list_con);
NS_str = num2str(NS);

% create results folder
results_folder = [root,filesep,'Second_level',filesep,...
    'GLM',GLM_str,'_',NS_str,'subs',filesep];
if exist(results_folder,'dir') ~= 7
    mkdir(results_folder);
end

batch_idx = 0;

%% 1) take mean anatomy
% extract parameters for mean calculation
mean_filename = ['mean_anat_',NS_str,'_subjects'];
wms_anat = cell(NS,1);
% add all the anat files for all the subjects
% extract anat EPI
for iSubject = 1:NS
    sub_anat_folder = [root,subject_id{iSubject},filesep,...
        'fMRI_analysis',filesep,'anatomical',filesep];
    wms_anat_name = ls([sub_anat_folder,'wms*']);
    wms_anat(iSubject) = {[sub_anat_folder, wms_anat_name]};
end

% spm calculation of the mean
batch_idx = batch_idx + 1;
matlabbatch{batch_idx}.spm.util.imcalc.input = wms_anat;
matlabbatch{batch_idx}.spm.util.imcalc.output = mean_filename;
matlabbatch{batch_idx}.spm.util.imcalc.outdir = {results_folder(1:end-1)};
matlabbatch{batch_idx}.spm.util.imcalc.expression = 'mean(X)';
matlabbatch{batch_idx}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{batch_idx}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{batch_idx}.spm.util.imcalc.options.mask = 0;
matlabbatch{batch_idx}.spm.util.imcalc.options.interp = 1;
matlabbatch{batch_idx}.spm.util.imcalc.options.dtype = 4;


%% 2) 2nd level concatenation
%% enter contrasts in spm
for iContrast = 1:n_con
    con_str = num2str(iContrast);
    
    % directory for concatenated contrast
    con_folder = [results_folder, list_con{iContrast}];
    mkdir(con_folder);
    
    % start second level:
    batch_idx = batch_idx + 1;
    matlabbatch{batch_idx}.spm.stats.factorial_design.dir = {con_folder};
    % list of inputs
    conlist = cell(NS,1); % 1 con per EPI-subject
    % extract MBBjune2016 contrasts
    for iSubject = 1:NS
        subject_folder = [root,subject_id{iSubject}, filesep, 'fMRI_analysis' filesep, 'functional' filesep, 'GLM',GLM_str, filesep];
        if iContrast < 10
            conlist(iSubject) = {[subject_folder,'con_000',con_str,'.nii,1']};
        elseif iContrast >= 10 && iContrast < 100
            conlist(iSubject) = {[subject_folder,'con_00',con_str,'.nii,1']};
        elseif iContrast >= 100 && iContrast < 1000
            conlist(iSubject) = {[subject_folder,'con_0',con_str,'.nii,1']};
        elseif iContrast >= 1000 && iContrast < 10000
            conlist(iSubject) = {[subject_folder,'con_',con_str,'.nii,1']};
        end
    end
    
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t1.scans = conlist;
    
    % default parameters
    matlabbatch{batch_idx}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{batch_idx}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{batch_idx}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{batch_idx}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{batch_idx}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{batch_idx}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    % model estimation
    batch_model_rtg = batch_idx;
    batch_idx = batch_idx + 1;
    matlabbatch{batch_idx}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{batch_model_rtg}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{batch_idx}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{batch_idx}.spm.stats.fmri_est.method.Classical = 1;
    % contrast
    batch_estm_rtg = batch_idx;
    batch_idx = batch_idx + 1;
    matlabbatch{batch_idx}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{batch_estm_rtg}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.name = list_con{iContrast};
    matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{batch_idx}.spm.stats.con.delete = 0;
    
end

cd(results_folder);

% display spm batch before running it
% spm_jobman('interactive',matlabbatch);
spm_jobman('run',matlabbatch);
end