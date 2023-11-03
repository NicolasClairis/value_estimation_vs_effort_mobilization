% This script will group 2 groups of t.maps together in order to perform a
% conjunction between 2 main contrasts of a GLM.
%
% See also: First_level_batch, contrasts_batch

clear all; close all; clc;

%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% which GLM?
GLM = spm_input('GLM number?',1,'e');
GLM_str = num2str(GLM);

%% working directory
root = 'enter path here';
%% load subjects
subject_id = {'enter list of subjects here'};
NS = length(subject_id);

% load GLM parameters
[GLMprm] = which_GLM_MS2(GLM);
%% load contrasts (names + order)
[list_con, ~] = MS2_load_con(GLMprm, subject_id{1});
nb_max_con = length(list_con);

% list of contrasts of interest
list_con_question = list_con{1};
for con = 2:nb_max_con
    list_con_question = [list_con_question,' | ',list_con{con}];
end
list_con_question = [list_con_question,' | If all wished contrasts selected click here'];

% select contrasts to display at the end
which_con = zeros(1,nb_max_con); % column with 0/1 to see which contrasts to display at the end
stop_con_loop = 0;
while stop_con_loop == 0 % select the contrasts you want for the conjunction
    selectedContrast = spm_input('Which contrasts do you want to include for the conjunction?',1,'m',...
        list_con_question, ...
        1:(nb_max_con+1), 0);
    if selectedContrast < nb_max_con + 1
        which_con(selectedContrast) = 1; % 1 for selected contrasts
    elseif selectedContrast == nb_max_con + 1 % stop contrast selection loop
        stop_con_loop = 1;
    end
end

conOfInterest = find(which_con == 1);
conNames = [list_con{conOfInterest(1)},'_', list_con{conOfInterest(2)}];
con_length = length(conOfInterest);
NS_str = num2str(NS);

% create results folder
results_folder = [root,'Second_level' filesep,'Conjunction_Second_level',filesep,'GLM',GLM_str,'_',NS_str,'subs_',conNames,filesep];
if exist(results_folder,'dir') ~= 7
    mkdir(results_folder);
end

taskResults_folder = [results_folder, 'GLM',GLM_str,'_',NS_str,'subs_',conNames,filesep];
% create results folder
if exist(taskResults_folder,'dir') ~= 7
    mkdir(taskResults_folder);
end

batch_idx = 0;

%% 1) take mean anatomy
% extract parameters for mean calculation
mean_filename = ['mean_anat_',NS_str,'_subjects'];
wms_anat = cell(NS,1);
% add all the anat files for all the subjects
% extract anat EPI
for iS = 1:NS
    wms_anat_name = ls([root,subject_id{iS}, filesep, 'fMRI_analysis',filesep,'anatomical',filesep 'wms*']);
    wms_anat(iS) = {[root,subject_id{iS}, filesep, 'fMRI_analysis', filesep, 'anatomical', filesep, wms_anat_name]};
end

% spm calculation of the mean
batch_idx = batch_idx + 1;
matlabbatch{batch_idx}.spm.util.imcalc.input = wms_anat;
matlabbatch{batch_idx}.spm.util.imcalc.output = mean_filename;
matlabbatch{batch_idx}.spm.util.imcalc.outdir = {taskResults_folder(1:end-1)};
matlabbatch{batch_idx}.spm.util.imcalc.expression = 'mean(X)';
matlabbatch{batch_idx}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{batch_idx}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{batch_idx}.spm.util.imcalc.options.mask = 0;
matlabbatch{batch_idx}.spm.util.imcalc.options.interp = 1;
matlabbatch{batch_idx}.spm.util.imcalc.options.dtype = 4;


%% 2nd level concatenation of EPIs
%% enter contrasts in spm
batch_idx = batch_idx + 1;

% start second level:
matlabbatch{batch_idx}.spm.stats.factorial_design.dir = {taskResults_folder};
% list of inputs
conlist = cell(NS,con_length); % 1 con per EPI-subject
% extract MBBjune2016 contrasts
for iContrast = 1:con_length
    con_str = num2str(conOfInterest(iContrast));
    for iS = 1:NS
        sub_firstLevel_folder = [root,subject_id{iS}, filesep, 'fMRI_analysis' filesep, 'functional' filesep];
        subject_folder = [sub_firstLevel_folder, 'GLM',GLM_str, filesep];
        if conOfInterest(iContrast) < 10
            conlist(iS,iContrast) = {[subject_folder,'con_000',con_str,'.nii,1']};
        elseif conOfInterest(iContrast) >= 10 && conOfInterest(iContrast) < 100
            conlist(iS,iContrast) = {[subject_folder,'con_00',con_str,'.nii,1']};
        elseif conOfInterest(iContrast) >= 100 && conOfInterest(iContrast) < 1000
            conlist(iS,iContrast) = {[subject_folder,'con_0',con_str,'.nii,1']};
        end
    end
end

%% Be careful t2.scans (for two-sample t.test)
matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.scans1 = conlist(:,1);
matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.scans2 = conlist(:,2);

matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.dept = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{batch_idx}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});

matlabbatch{batch_idx}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{batch_idx}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{batch_idx}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.glonorm = 1;

%% model estimation
batch_model_rtg = batch_idx;
batch_idx = batch_idx + 1;
matlabbatch{batch_idx}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{batch_model_rtg}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{batch_idx}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{batch_idx}.spm.stats.fmri_est.method.Classical = 1;

%% t.contrast
batch_estm_rtg = batch_idx;
batch_idx = batch_idx + 1;
matlabbatch{batch_idx}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{batch_estm_rtg}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

for iContrast = 1:con_length
    matlabbatch{batch_idx}.spm.stats.con.consess{iContrast}.tcon.name = list_con{conOfInterest(iContrast)};
    matlabbatch{batch_idx}.spm.stats.con.consess{iContrast}.tcon.weights = [zeros(1, iContrast - 1), 1, zeros(1,con_length - iContrast)];
    matlabbatch{batch_idx}.spm.stats.con.consess{iContrast}.tcon.sessrep = 'none';
    matlabbatch{batch_idx}.spm.stats.con.delete = 0;
end

%% f.contrast
fIdentityConne = eye(con_length);
batch_idx = batch_idx + 1;
matlabbatch{batch_idx}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{batch_estm_rtg}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{batch_idx}.spm.stats.con.consess{1}.fcon.name = 'Fcontrast';
matlabbatch{batch_idx}.spm.stats.con.consess{1}.fcon.weights = fIdentityConne;
matlabbatch{batch_idx}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.delete = 0;

%%
cd(results_folder);

%% display spm batch before running it
% spm_jobman('interactive',matlabbatch);
spm_jobman('run',matlabbatch);