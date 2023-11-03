function[stats_Val_DT] = MS2_RL_Val_DT_correl(model_n, SV_formula)
%[stats_Val_DT] = MS2_RL_Val_DT_correl(model_n, SV_formula)
% MS2_RL_Val_DT_correl will check at the correlation between value (Val) and
% deliberation time (DT) in the Reinforcement Learning (RL) task.
%
% INPUTS
% model_n: model number (model 6 by default)
%
% SV_formula: use QA+QB ('qAqB') or pA*QA+pB*QB ('pAqApBqB')
% (the latter will be selected by default)
%
% OUTPUTS
% stats_Val_DT: statistics for Val and DT correlation including average
% p.value, average R2, average Pearson's correlation coefficient, etc.
%
% See also MS2_RL_RT_ter for a different way of computing this correlation

%% working directory
root = 'define path here';

%% subject identification
subject_id = {'enter list of subjects here'};
NS = length(subject_id);
%% model identification
if ~exist('model_n','var') || isempty(model_n)
    model_n = 6;
end
if ~exist('SV_formula','var') || isempty(SV_formula)
    SV_formula = 'pAqApBqB';
end
%% main parameters
n_RL_runs = 3;

[stats_Val_DT.R2.perRun_perSub,...
    stats_Val_DT.pval.perRun_perSub,...
    stats_Val_DT.pearsonR.perRun_perSub] = deal(NaN(n_RL_runs,NS));
[stats_Val_DT.R2.perSub,...
    stats_Val_DT.pval.perSub,...
    stats_Val_DT.pearsonR.perSub] = deal(NaN(1,NS));

%% subject loop
for iS = 1:NS
    
    % subject identification
    sub_nm = subject_id{iS};
    if strcmp(sub_nm(3),'_')
        subid   = sub_nm(2);
    elseif strcmp(sub_nm(3),'_') == 0 && strcmp(sub_nm(4),'_')
        subid   = sub_nm(2:3);
    end
    % subject working directories
    sub_folder      = [root filesep sub_nm filesep];
    onsets_folder = [sub_folder 'fMRI_analysis' filesep];
    RL_runs_idx = MS2_task_runs_extraction('RL',sub_nm);
    
    %% run loop
    for iRun = 1:n_RL_runs
        run_nm = num2str(RL_runs_idx(iRun));
        %% load data
        load_onsetStruct = getfield(load([onsets_folder,'onsets_sub',subid,'_learning_run',run_nm,'.mat']),'learn');
        RT_GLP_run_tmp = load_onsetStruct.mod.RT_fp.GL_Pairs.main.raw;
        switch SV_formula
            case 'qAqB'
                SV_GLP_run_tmp = load_onsetStruct.mod.Q_model(model_n).raw.GL_Pairs.GL_PairsTrials.bestItem +...
                    load_onsetStruct.mod.Q_model(model_n).raw.GL_Pairs.GL_PairsTrials.worseItem;
            case 'pAqApBqB'
                SV_GLP_run_tmp = load_onsetStruct.mod.Q_model(model_n).raw.SV.GL_Pairs.GL_PairsTrials;
        end
        
        %% test correlation
        % GLM to get p.value
        [~,~,stats_tmp1] = glmfit(SV_GLP_run_tmp, RT_GLP_run_tmp, 'normal');
        stats_Val_DT.pval.perRun_perSub(iRun, iS) = stats_tmp1.p(2);
        % GLM to get R2
        stats_tmp2 = fitlm(SV_GLP_run_tmp, RT_GLP_run_tmp);
        stats_Val_DT.R2.perRun_perSub(iRun, iS) = stats_tmp2.Rsquared.Adjusted;
        % Pearson's correlation
        coeffs = corrcoef(SV_GLP_run_tmp, RT_GLP_run_tmp);
        stats_Val_DT.pearsonR.perRun_perSub(iRun, iS) = coeffs(1,2);
    end % run loop
    
    %% average correlation across runs
    stats_Val_DT.R2.perSub(iS) = mean(stats_Val_DT.R2.perRun_perSub(:, iS),1,'omitnan');
    stats_Val_DT.pval.perSub(iS) = mean(stats_Val_DT.pval.perRun_perSub(:, iS),1,'omitnan');
    stats_Val_DT.pearsonR.perSub(iS) = mean(stats_Val_DT.pearsonR.perRun_perSub(:, iS),1,'omitnan');
end % subject loop

%% average correlation across subjects
[stats_Val_DT.R2.aSubs.mean,...
    stats_Val_DT.R2.aSubs.sem] = mean_sem_sd(stats_Val_DT.R2.perSub, 2);
[stats_Val_DT.pval.aSubs.mean,...
    stats_Val_DT.pval.aSubs.sem] = mean_sem_sd(stats_Val_DT.pval.perSub, 2);
[stats_Val_DT.pearsonR.aSubs.mean,...
    stats_Val_DT.pearsonR.aSubs.sem] = mean_sem_sd(stats_Val_DT.pearsonR.perSub, 2);
end % function