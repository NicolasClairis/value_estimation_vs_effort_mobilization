function [ mChoices ] = MS2_RL_avg_choices( analysis_type )
%[ mChoices ] = MS2_RL_avg_choices( analysis_type )
% MS2_RL_avg_choices will get the average amount obtained per session for each
% subject in the reinforcement learning task.
%
% INPUT
% analysis_type: string indicating which subjects to include
%
% OUTPUT
% mChoices: structure with average proportion of choices

%% working directory
root = 'enter path here';

%% subject identification
subject_id = {'enter list of subjects here'};
NS = length(subject_id);

%% load data
n_RL_runs = 3;
[mChoices.perSub.perRun.gains, mChoices.perSub.perRun.losses] = deal(NaN(n_RL_runs, NS));
[mChoices.perSub.aRuns.gains, mChoices.perSub.aRuns.losses] = deal(NaN(1,NS));
for iS = 1:NS
    %% subject info
    sub_nm = subject_id{iS};
   if strcmp(sub_nm(3),'_')
        subid   = sub_nm(2);
    elseif ~strcmp(sub_nm(3),'_') && strcmp(sub_nm(4),'_')
        subid = sub_nm(2:3);
    end
    sub_folder = [root, sub_nm, filesep 'behavior', filesep];
    RL_runs_idx = MS2_task_runs_extraction('RL',sub_nm);
    
    %% load data
    for iRL_run = 1:n_RL_runs
        run_nb = RL_runs_idx(iRL_run);
        data_tmp = load([sub_folder, 'global_sub_',subid,'_session_',num2str(run_nb),'_learning.mat'],'response','npair');
        response_tmp = data_tmp.response;
        npair_tmp = data_tmp.npair;
        mChoices.perSub.perRun.gains(iRL_run,iS) = sum(response_tmp(npair_tmp == 1)==1)./sum(ismember(response_tmp(npair_tmp == 1),[-1,1]));
        mChoices.perSub.perRun.losses(iRL_run,iS) = sum(response_tmp(npair_tmp == 3)==1)./sum(ismember(response_tmp(npair_tmp == 1),[-1,1]));
    end
    %% average across runs
    mChoices.perSub.aRuns.gains(iS) = mean(mChoices.perSub.perRun.gains(:,iS),1,'omitnan');
    mChoices.perSub.aRuns.losses(iS) = mean(mChoices.perSub.perRun.losses(:,iS),1,'omitnan');
end

%% average across subjects
[mChoices.m_aSubs.gains,...
    mChoices.sem_aSubs.gains] = mean_sem_sd(mChoices.perSub.aRuns.gains,2);
[mChoices.m_aSubs.losses,...
    mChoices.sem_aSubs.losses] = mean_sem_sd(mChoices.perSub.aRuns.losses,2);

end % function