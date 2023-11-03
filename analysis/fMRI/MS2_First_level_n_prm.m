function [ n_prm, prm_idx ] = MS2_First_level_n_prm( GLMprm, sub_nm )
%[ n_prm, prm_idx ] = MS2_First_level_n_prm( GLMprm, sub_nm )
% estimate number and identity of regressors for each run for each task for
% each subject of Motiscan 2
%
% INPUTS
% GLMprm: structure with GLM parameters
%
% sub_nm: subject identification name
%
% OUTPUTS
% n_prm: structure containing the number of parameters for each run for
% each task for the current subject
%
% prm_idx: structure with corresponding index for each regressor
%

%% extract main relevant parameters for the script to work

% need to know add derivative to know number of betas/regressor
add_drv = GLMprm.gal.add_drv;
% infer number of parameters/regressor accordingly
switch add_drv
    case 0
        nBperR = 1;
    case 1 % normal + temporal derivative
        nBperR = 2;
    case 2 % normal + temporal derivative + spatial derivative
        nBperR = 3;
end

% reinforcement-learning task
RLprm = GLMprm.RL;
% grip task
gripRPprm = GLMprm.grip;
% stroop task
stroopRPprm = GLMprm.stroop;

% total number of runs
n_total_RL_runs     = 3;
n_total_grip_runs   = 2;
n_total_stroop_runs = 2;
n_total_runs = n_total_RL_runs + n_total_grip_runs + n_total_stroop_runs;

%% movement parameters
n_mvmt_perRun = 6*ones(1,n_total_runs);

%% parameters of interest for output
n_prm.all = 0 ;
task_names = {'RL','grip','stroop'};
n_tasks = length(task_names);
for iTask = 1:n_tasks
    task_nm = task_names{iTask};
    
    % define run index for each task
    [runs.(task_nm)] = MS2_task_runs_extraction(task_nm, sub_nm);
    
    % set initial number of parameters per task to zero
    for iRun = 1:length(runs.(task_nm))
        jRun = runs.(task_nm)(iRun);
        run_nm = ['run',num2str(jRun)];
        n_prm.(task_nm).(run_nm) = 0;
    end
end % task loop

%% check which onsets and modulators are used in the current GLM and assign
% them the corresponding index in the list of parameters
bf_idx = 0; % number of regressors before the current one
prm_idx = struct;

for iRun = 1:n_total_runs
    run_nm = ['run',num2str(iRun)];
    
    
    %% RL
    if ismember(iRun, runs.RL)
        task_nm = 'RL';
        task_id = 'L_';
        
        %% stim display
        o_stim = RLprm.o_stim;
        mod_stim = RLprm.mod_stim;
        switch o_stim
            case 1 % all pairs pooled
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim'], bf_idx, sub_nm);
                
                % modulators
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT'], bf_idx, sub_nm);
                end
                % trial number
                if ismember(mod_stim.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN'], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,1:5)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV'], bf_idx, sub_nm);
                end
                
                % dQ/dQ
                if mod_stim.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ'], bf_idx, sub_nm);
                end
                
                % p(choice=best option)
                if mod_stim.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest'], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN'], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_stim.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity'], bf_idx, sub_nm);
                end
                
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT'], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,6:10)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV'], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN'], bf_idx, sub_nm);
                end
                
            case 2 % separate depending on pair type gain/neutral/loss
                
                %% gain pair
                trialType_nm = '_gainPair';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,1:5)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % dQ
                if mod_stim.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx, sub_nm);
                end
                
                % p(choice=best option)
                if mod_stim.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_stim.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,6:10)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% neutral pair
                trialType_nm = '_ntalPair';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % no SV or dQ for neutral pair (always equal to zero)
                
                % trial number
                if ismember(mod_stim.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_stim.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% loss pair
                trialType_nm = '_lossPair';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                % trial number
                if ismember(mod_stim.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,1:5)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % dQ
                if mod_stim.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx, sub_nm);
                end
                
                % p(choice = best option)
                if mod_stim.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_stim.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,6:10)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
            case 3 % gain + loss pooled, neutral apart
                %% gain + loss pair
                trialType_nm = '_GL_Pairs';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,1:5)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % dQ
                if mod_stim.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx, sub_nm);
                end
                
                % p(choice=best option)
                if mod_stim.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_stim.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,6:10)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% neutral pair
                trialType_nm = '_ntalPair';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % no SV or dQ for neutral pair (always equal to zero)
                
                % trial number
                if ismember(mod_stim.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_stim.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
            case 4
                error('not ready yet');
            case 5 % split gain/neutral/loss and first/second half
                %% gain pair first trials
                trialType_nm = '_gainPair_first';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,1:5)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % dQ
                if mod_stim.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx, sub_nm);
                end
                
                % p(choice=best option)
                if mod_stim.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_stim.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,6:10)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% neutral pair first trials
                trialType_nm = '_ntalPair_first';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % no SV, no dQ for neutral pair
                
                % trial number
                if ismember(mod_stim.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_stim.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                
                %% loss pair first trials
                trialType_nm = '_lossPair_first';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,1:5)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % dQ
                if mod_stim.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx, sub_nm);
                end
                
                % p(choice = best option)
                if mod_stim.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_stim.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,6:10)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% gain pair last trials
                trialType_nm = '_gainPair_last';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,1:5)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % dQ
                if mod_stim.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx, sub_nm);
                end
                
                % p(choice=best option)
                if mod_stim.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_stim.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,6:10)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% neutral pair last trials
                trialType_nm = '_ntalPair_last';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % no SV, no dQ for neutral pair
                
                % trial number
                if ismember(mod_stim.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_stim.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% loss pair last trials
                trialType_nm = '_lossPair_last';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,1:5)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % dQ
                if mod_stim.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx, sub_nm);
                end
                
                % p(choice = best option)
                if mod_stim.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_stim.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,6:10)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
            case 6 % pool gain+loss/neutral apart and split first/second half trials also
                %% gain + loss pair - first trials
                trialType_nm = '_GL_Pairs_first';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,1:5)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % dQ
                if mod_stim.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx, sub_nm);
                end
                
                % p(choice=best option)
                if mod_stim.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_stim.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,6:10)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% neutral pair - first trials
                trialType_nm = '_ntalPair_first';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % no SV or dQ for neutral pair (always equal to zero)
                
                % trial number
                if ismember(mod_stim.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_stim.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% gain + loss pair - last trials
                trialType_nm = '_GL_Pairs_last';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,1:5)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % dQ
                if mod_stim.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx, sub_nm);
                end
                
                % p(choice=best option)
                if mod_stim.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_stim.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if ismember(mod_stim.SV,6:10)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% neutral pair - last trials
                trialType_nm = '_ntalPair_last';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % no SV or dQ for neutral pair (always equal to zero)
                
                % trial number
                if ismember(mod_stim.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_stim.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_stim.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx, sub_nm);
                end
        end % stimulus onset
        
        %% answer
        o_answer = RLprm.o_answer;
        %         mod_answer = RLprm.mod_answer;
        switch o_answer
            case 1
                % onset answer
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_answer'], bf_idx, sub_nm);
                
            case 2 % separate depending on pair type
                
                %% gain pair
                trialType_nm = '_gainPair';
                % onset answer
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_answer',trialType_nm], bf_idx, sub_nm);
                
                %% neutral pair
                trialType_nm = '_ntalPair';
                % onset answer
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_answer',trialType_nm], bf_idx, sub_nm);
                
                %% loss pair
                trialType_nm = '_lossPair';
                % onset answer
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_answer',trialType_nm], bf_idx, sub_nm);
                
        end % answer onset
        
        %% chosen option in red
        o_chosen = RLprm.o_chosen;
        mod_chosen = RLprm.mod_chosen;
        switch o_chosen
            case 1
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen'], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if ismember(mod_chosen.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN'], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if mod_chosen.SV ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_SV'], bf_idx, sub_nm);
                end
                
                % dQ/dQ
                if mod_chosen.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_dQ'], bf_idx, sub_nm);
                end
                
                % p(choice=best option)
                if mod_chosen.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_pBest'], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN'], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_chosen.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity'], bf_idx, sub_nm);
                end
                
                % RT
                if mod_chosen.RT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT'], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN'], bf_idx, sub_nm);
                end
                
            case 2 % separate depending on pair type
                
                %% gain pair
                trialType_nm = '_gainPair';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if ismember(mod_chosen.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if mod_chosen.SV ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % dQ
                if mod_chosen.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_dQ',trialType_nm], bf_idx, sub_nm);
                end
                
                % p(choice=best option)
                if mod_chosen.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_pBest',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_chosen.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_chosen.RT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% neutral pair
                trialType_nm = '_ntalPair';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if ismember(mod_chosen.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % no SV or dQ for neutral pair (always equal to zero)
                
                % trial number
                if ismember(mod_chosen.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_chosen.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_chosen.RT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% loss pair
                trialType_nm = '_lossPair';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if ismember(mod_chosen.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if mod_chosen.SV ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % dQ
                if mod_chosen.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_dQ',trialType_nm], bf_idx, sub_nm);
                end
                
                % p(choice = best option)
                if mod_chosen.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_pBest',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_chosen.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_chosen.RT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
            case 3
                %% gain + loss pair
                trialType_nm = '_GL_Pairs';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if ismember(mod_chosen.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if mod_chosen.SV ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % dQ
                if mod_chosen.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_dQ',trialType_nm], bf_idx, sub_nm);
                end
                
                % p(choice=best option)
                if mod_chosen.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_pBest',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_chosen.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_chosen.RT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% neutral pair
                trialType_nm = '_ntalPair';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if ismember(mod_chosen.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % no SV or dQ for neutral pair (always equal to zero)
                
                % trial number
                if ismember(mod_chosen.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_chosen.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_chosen.RT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
            case 4
                error('not ready yet');
            case 5
                %% gain pair first trials
                trialType_nm = '_gainPair_first';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if ismember(mod_chosen.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if mod_chosen.SV ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % dQ
                if mod_chosen.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_dQ',trialType_nm], bf_idx, sub_nm);
                end
                
                % p(choice=best option)
                if mod_chosen.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_pBest',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_chosen.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_chosen.RT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% neutral pair first trials
                trialType_nm = '_ntalPair_first';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if ismember(mod_chosen.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % no SV, no dQ for neutral pair
                
                % trial number
                if ismember(mod_chosen.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_chosen.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_chosen.RT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                
                %% loss pair first trials
                trialType_nm = '_lossPair_first';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if ismember(mod_chosen.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if mod_chosen.SV ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % dQ
                if mod_chosen.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_dQ',trialType_nm], bf_idx, sub_nm);
                end
                
                % p(choice = best option)
                if mod_chosen.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_pBest',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_chosen.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_chosen.RT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% gain pair last trials
                trialType_nm = '_gainPair_last';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if ismember(mod_chosen.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if mod_chosen.SV ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % dQ
                if mod_chosen.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_dQ',trialType_nm], bf_idx, sub_nm);
                end
                
                % p(choice=best option)
                if mod_chosen.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_pBest',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_chosen.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_chosen.RT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% neutral pair last trials
                trialType_nm = '_ntalPair_last';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if ismember(mod_chosen.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % no SV, no dQ for neutral pair
                
                % trial number
                if ismember(mod_chosen.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_chosen.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_chosen.RT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% loss pair last trials
                trialType_nm = '_lossPair_last';
                % onset stim
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if ismember(mod_chosen.trialN,[1,2,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % SV = pA*QA+pB*QB
                if mod_chosen.SV ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_SV',trialType_nm], bf_idx, sub_nm);
                end
                
                % dQ
                if mod_chosen.dQ ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_dQ',trialType_nm], bf_idx, sub_nm);
                end
                
                % p(choice = best option)
                if mod_chosen.pBest ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_pBest',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[4,5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_chosen.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT
                if mod_chosen.RT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number
                if ismember(mod_chosen.trialN,[7,8,9])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx, sub_nm);
                end
        end % chosen option in red
        
        %% feedback
        o_fbk = RLprm.o_fbk;
        mod_fbk = RLprm.mod_fbk;
        switch o_fbk
            case 1
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk'], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if mod_fbk.trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN'], bf_idx, sub_nm);
                end
                
                % feedback
                if mod_fbk.fbk ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk'], bf_idx, sub_nm);
                end
                
                % PE
                if mod_fbk.PE ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE'], bf_idx, sub_nm);
                end
                
                % PE
                if mod_fbk.PE_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis'], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain'], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity'], bf_idx, sub_nm);
                end
                
            case 2 % separate depending on pair type (gain/neutral/loss)
                
                %% gain pair
                trialType_nm = '_gainPair';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if mod_fbk.trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback
                if mod_fbk.fbk ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE
                if mod_fbk.PE ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE bis
                if mod_fbk.PE_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                %% neutral pair
                trialType_nm = '_ntalPair';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if mod_fbk.trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                %% loss pair
                trialType_nm = '_lossPair';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if mod_fbk.trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback
                if mod_fbk.fbk ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE
                if mod_fbk.PE ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE bis
                if mod_fbk.PE_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
            case 3 % separate depending on feedback type (gain/neutral/loss)
                %% gain feedback
                trialType_nm = '_gainFbk';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if mod_fbk.trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback
                if mod_fbk.fbk ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE
                if mod_fbk.PE ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE bis
                if mod_fbk.PE_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                %% neutral feedback
                trialType_nm = '_ntalFbk';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if mod_fbk.trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback: none always 0 for neutral feedback...
                
                % PE
                if mod_fbk.PE ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE bis
                if mod_fbk.PE_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                %% loss feedback
                trialType_nm = '_lossFbk';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if mod_fbk.trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback
                if mod_fbk.fbk ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE
                if mod_fbk.PE ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE bis
                if mod_fbk.PE_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
            case 4 % separate depending on feedback AND pair type
                %% gain pair - gain feedback
                trialType_nm = '_gainPair_gainFbk';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if mod_fbk.trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback
                if mod_fbk.fbk ~= 0
                    error('Pointless: feedback = stable for a given feedback type.');
                end
                
                % PE
                if mod_fbk.PE ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE bis
                if mod_fbk.PE_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                %% gain pair - neutral feedback
                trialType_nm = '_gainPair_ntalFbk';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if mod_fbk.trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback
                if mod_fbk.fbk ~= 0
                    error('Pointless: feedback = stable for a given feedback type.');
                end
                
                % PE
                if mod_fbk.PE ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE bis
                if mod_fbk.PE_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                %% neutral pair
                trialType_nm = '_ntalPair_ntalFbk';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if mod_fbk.trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback
                if mod_fbk.fbk ~= 0
                    error('Pointless: feedback = stable for a given feedback type.');
                end
                
                % PE: pointless should remain equal 0 always for neutral
                % pair
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                %% loss pair
                trialType_nm = '_lossPair_ntalFbk';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if mod_fbk.trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback
                if mod_fbk.fbk ~= 0
                    error('Pointless: feedback = stable for a given feedback type.');
                end
                
                % PE
                if mod_fbk.PE ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE bis
                if mod_fbk.PE_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                %% loss pair - loss feedback
                trialType_nm = '_lossPair_lossFbk';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if mod_fbk.trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback
                if mod_fbk.fbk ~= 0
                    error('Pointless: feedback = stable for a given feedback type.');
                end
                
                % PE
                if mod_fbk.PE ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE bis
                if mod_fbk.PE_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
            case 5 % pool gain+loss pairs, separate neutral apart
                
                %% gain +loss pair
                trialType_nm = '_GL_Pairs';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if mod_fbk.trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback
                if mod_fbk.fbk ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE
                if mod_fbk.PE ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE bis
                if mod_fbk.PE_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                %% neutral pair
                trialType_nm = '_ntalPair';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if mod_fbk.trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
            case 6 % pool gain+loss pairs (but split first/second half trials) and keep neutral pair apart
                
                %% gain +loss pair (first trials)
                trialType_nm = '_GL_Pairs_first';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if mod_fbk.trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback
                if mod_fbk.fbk ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE
                if mod_fbk.PE ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE bis
                if mod_fbk.PE_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                %% gain +loss pair (last trials)
                trialType_nm = '_GL_Pairs_last';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if mod_fbk.trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback
                if mod_fbk.fbk ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE
                if mod_fbk.PE ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx, sub_nm);
                end
                
                % PE bis
                if mod_fbk.PE_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                %% neutral pair
                trialType_nm = '_ntalPair';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % trial number
                if mod_fbk.trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
        end
        
        %% cross
        o_cross = RLprm.o_cross;
        switch o_cross
            case 1
                % onset cross
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_cross'], bf_idx, sub_nm);
        end
        
        %% missed trials
        o_missed_trials_stim = RLprm.o_missed_trials_stim;
        switch o_missed_trials_stim
            case 1
                % onset answer
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_missed_trials_stim'], bf_idx, sub_nm);
        end
        
        
        
        %% grip
    elseif ismember(iRun, runs.grip)
        task_nm = 'grip';
        task_id = 'G_';
        
        %% incentive display
        o_inc = gripRPprm.o_inc;
        mod_inc = gripRPprm.mod_inc;
        switch o_inc
            case 1 % all onsets grouped
                % onset incentive
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_inc'], bf_idx, sub_nm);
                
                % modulators
                % Gain/Loss condition
                if mod_inc.GL_cond ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_GL_cond'], bf_idx, sub_nm);
                end
                
                % luminance (first)
                if ismember(mod_inc.lum,[2,4,6,8])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_lum'], bf_idx, sub_nm);
                end
                
                % trial number (first)
                if ismember(mod_inc.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN'], bf_idx, sub_nm);
                end

                % RT first press
                if ismember(mod_inc.RT_fp,[2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_fp'], bf_idx, sub_nm);
                end

                % model variables
                if mod_inc.mdl_n ~= 0
                    % ressource
                    if ismember(mod_inc.E_pred,[5,6])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred'], bf_idx, sub_nm);
                    end
                end
                
                % reward displayed as piece of money or bill?
                if mod_inc.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_R_type'], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_inc.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc'], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_inc.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc_bis'], bf_idx, sub_nm);
                end

                % confidence
                if mod_inc.conf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_conf'], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_inc.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_GLM',mod_inc.ROI_activity_GLM,...
                    '_',mod_inc.ROI_activity_period,'_period_',mod_inc.ROI_activity_ROI_nm,'_activity'], bf_idx, sub_nm);
                end
                
                % effort
                if mod_inc.Effort ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_effort'], bf_idx, sub_nm);
                end
                
                % model variables
                if mod_inc.mdl_n ~= 0
                    
                    % ressource
                    if ismember(mod_inc.E_pred,[3,4])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred'], bf_idx, sub_nm);
                    end
                    
                    % cost
                    if mod_inc.mdl_cost ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_cost'], bf_idx, sub_nm);
                    end
                    
                    % benefit
                    if mod_inc.mdl_benefit ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_benefit'], bf_idx, sub_nm);
                    end
                    
                    % EV
                    if mod_inc.mdl_EV ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_EV'], bf_idx, sub_nm);
                    end
                end
                
                % trial number (middle)
                if ismember(mod_inc.trialN,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN'], bf_idx, sub_nm);
                end
                
                % RT first press
                if ismember(mod_inc.RT_fp,[1])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_fp'], bf_idx, sub_nm);
                end
                
                % total gain until now
                if mod_inc.totalGain_prev ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_totalGain_prev'], bf_idx, sub_nm);
                end
                
                % sum(performance previous trials)
                if mod_inc.sumPerfPrev ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_sumPerfPrev'], bf_idx, sub_nm);
                end
                
                % luminance (last)
                if ismember(mod_inc.lum,[1,3,5,7])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_lum'], bf_idx, sub_nm);
                end
                
                % trial number (last)
                if ismember(mod_inc.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN'], bf_idx, sub_nm);
                end
                
                % incentive*trial number
                if mod_inc.inc_x_trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc_x_trialN'], bf_idx, sub_nm);
                end
                
                % ressource
                if ismember(mod_inc.E_pred,[1,2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred'], bf_idx, sub_nm);
                end
                
                % performance
                if mod_inc.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_perf'], bf_idx, sub_nm);
                end
                
            case 2 % to gain/to lose trials split
                %% to gain trials
                trialType_nm = '_toGain';
                % onset incentive
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_inc',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % luminance (first)
                if ismember(mod_inc.lum,[2,4,6,8])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (first)
                if ismember(mod_inc.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % RT first press
                if ismember(mod_inc.RT_fp,[2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_fp',trialType_nm], bf_idx, sub_nm);
                end

                % model variables
                if mod_inc.mdl_n ~= 0
                    % ressource
                    if ismember(mod_inc.E_pred,[5,6])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % reward displayed as piece of money or bill?
                if mod_inc.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_inc.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_inc.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc_bis',trialType_nm], bf_idx, sub_nm);
                end

                % confidence
                if mod_inc.conf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_conf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_inc.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_GLM',mod_inc.ROI_activity_GLM,...
                    '_',mod_inc.ROI_activity_period,'_period_',mod_inc.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % effort
                if mod_inc.Effort ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_effort',trialType_nm], bf_idx, sub_nm);
                end
                
                % model variables
                if mod_inc.mdl_n ~= 0
                    
                    % ressource
                    if ismember(mod_inc.E_pred,[3,4])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % cost
                    if mod_inc.mdl_cost ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_cost',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % benefit
                    if mod_inc.mdl_benefit ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_benefit',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % EV
                    if mod_inc.mdl_EV ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_EV',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % trial number (middle)
                if ismember(mod_inc.trialN,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT first press
                if ismember(mod_inc.RT_fp,[1])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_fp',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain until now
                if mod_inc.totalGain_prev ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_totalGain_prev',trialType_nm], bf_idx, sub_nm);
                end
                
                % sum(performance previous trials)
                if mod_inc.sumPerfPrev ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_sumPerfPrev',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                if ismember(mod_inc.lum,[1,3,5,7])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (last)
                if ismember(mod_inc.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive*trial number
                if mod_inc.inc_x_trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc_x_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if ismember(mod_inc.E_pred,[1,2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_inc.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_perf',trialType_nm], bf_idx, sub_nm);
                end
                
                %% to lose trials
                trialType_nm = '_toLose';
                % onset incentive
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_inc',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % luminance (first)
                if ismember(mod_inc.lum,[2,4,6,8])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (first)
                if ismember(mod_inc.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % RT first press
                if ismember(mod_inc.RT_fp,[2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_fp',trialType_nm], bf_idx, sub_nm);
                end

                % model variables
                if mod_inc.mdl_n ~= 0
                    % ressource
                    if ismember(mod_inc.E_pred,[5,6])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % reward displayed as piece of money or bill?
                if mod_inc.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_inc.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_inc.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc_bis',trialType_nm], bf_idx, sub_nm);
                end

                % confidence
                if mod_inc.conf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_conf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_inc.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_GLM',mod_inc.ROI_activity_GLM,...
                    '_',mod_inc.ROI_activity_period,'_period_',mod_inc.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % effort
                if mod_inc.Effort ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_effort',trialType_nm], bf_idx, sub_nm);
                end
                
                % model variables
                if mod_inc.mdl_n ~= 0
                    
                    % ressource
                    if ismember(mod_inc.E_pred,[3,4])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % cost
                    if mod_inc.mdl_cost ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_cost',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % benefit
                    if mod_inc.mdl_benefit ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_benefit',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % EV
                    if mod_inc.mdl_EV ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_EV',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % trial number (middle)
                if ismember(mod_inc.trialN,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT first press
                if ismember(mod_inc.RT_fp,[1])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_fp',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain until now
                if mod_inc.totalGain_prev ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_totalGain_prev',trialType_nm], bf_idx, sub_nm);
                end
                
                % sum(performance previous trials)
                if mod_inc.sumPerfPrev ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_sumPerfPrev',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                if ismember(mod_inc.lum,[1,3,5,7])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (last)
                if ismember(mod_inc.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive*trial number
                if mod_inc.inc_x_trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc_x_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if ismember(mod_inc.E_pred,[1,2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_inc.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_perf',trialType_nm], bf_idx, sub_nm);
                end
                
        end % incentive
        
        %% effort scale display
        o_dispE = gripRPprm.o_dispE;
        mod_dispE = gripRPprm.mod_dispE;
        switch o_dispE
            case 1
                % onset effort scale
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_dispE'], bf_idx, sub_nm);
                
                % modulators
                % Gain/Loss condition
                if mod_dispE.GL_cond ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_GL_cond'], bf_idx, sub_nm);
                end
                
                % luminance (first)
                if ismember(mod_dispE.lum,[2,4,6,8])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_lum'], bf_idx, sub_nm);
                end
                
                % trial number (first)
                if ismember(mod_dispE.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN'], bf_idx, sub_nm);
                end

                % RT first press
                if ismember(mod_dispE.RT_fp,[2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_fp'], bf_idx, sub_nm);
                end

                % model variables
                if mod_dispE.mdl_n ~= 0
                    % ressource
                    if ismember(mod_dispE.E_pred,[5,6])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred'], bf_idx, sub_nm);
                    end
                end
                
                % reward displayed as piece of money or bill?
                if mod_dispE.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_R_type'], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_dispE.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc'], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_dispE.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc_bis'], bf_idx, sub_nm);
                end

                % confidence
                if mod_dispE.conf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_conf'], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_dispE.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_GLM',mod_dispE.ROI_activity_GLM,...
                    '_',mod_dispE.ROI_activity_period,'_period_',mod_dispE.ROI_activity_ROI_nm,'_activity'], bf_idx, sub_nm);
                end
                
                % force exerted
                if mod_dispE.Effort ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_Effort'], bf_idx, sub_nm);
                end
                
                % model variables
                if mod_dispE.mdl_n ~= 0
                    
                    % ressource
                    if ismember(mod_dispE.E_pred,[3,4])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred'], bf_idx, sub_nm);
                    end
                    
                    % cost
                    if mod_dispE.mdl_cost ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_cost'], bf_idx, sub_nm);
                    end
                    
                    % benefit
                    if mod_dispE.mdl_benefit ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_benefit'], bf_idx, sub_nm);
                    end
                    
                    % EV
                    if mod_dispE.mdl_EV ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_EV'], bf_idx, sub_nm);
                    end
                end
                
                % trial number (middle)
                if ismember(mod_dispE.trialN,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN'], bf_idx, sub_nm);
                end
                
                % RT first press
                if ismember(mod_dispE.RT_fp,[1])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_fp'], bf_idx, sub_nm);
                end
                
                % luminance (last)
                if ismember(mod_dispE.lum,[1,3,5,7])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_lum'], bf_idx, sub_nm);
                end
                
                % trial number (last)
                if ismember(mod_dispE.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN'], bf_idx, sub_nm);
                end
                
                % incentive*trial number
                if mod_dispE.inc_x_trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc_x_trialN'], bf_idx, sub_nm);
                end
                
                % ressource
                if ismember(mod_dispE.E_pred,[1,2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred'], bf_idx, sub_nm);
                end
                
                % performance
                if mod_dispE.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_perf'], bf_idx, sub_nm);
                end
                
            case 2
                %% to gain trials
                trialType_nm = '_toGain';
                % onset effort scale
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_dispE',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % luminance (first)
                if ismember(mod_dispE.lum,[2,4,6,8])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (first)
                if ismember(mod_dispE.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % RT first press
                if ismember(mod_dispE.RT_fp,[2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_fp',trialType_nm], bf_idx, sub_nm);
                end

                % model variables
                if mod_dispE.mdl_n ~= 0
                    % ressource
                    if ismember(mod_dispE.E_pred,[5,6])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % reward displayed as piece of money or bill?
                if mod_dispE.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_dispE.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_dispE.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc_bis',trialType_nm], bf_idx, sub_nm);
                end

                % confidence
                if mod_dispE.conf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_conf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_dispE.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_GLM',mod_dispE.ROI_activity_GLM,...
                    '_',mod_dispE.ROI_activity_period,'_period_',mod_dispE.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % force exerted
                if mod_dispE.Effort ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_Effort',trialType_nm], bf_idx, sub_nm);
                end
                
                % model variables
                if mod_dispE.mdl_n ~= 0
                    
                    % ressource
                    if ismember(mod_dispE.E_pred,[3,4])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                
                    % cost
                    if mod_dispE.mdl_cost ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_cost',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % benefit
                    if mod_dispE.mdl_benefit ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_benefit',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % EV
                    if mod_dispE.mdl_EV ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_EV',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % trial number (middle)
                if ismember(mod_dispE.trialN,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT first press
                if ismember(mod_dispE.RT_fp,[1])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_fp',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                if ismember(mod_dispE.lum,[1,3,5,7])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (last)
                if ismember(mod_dispE.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive*trial number
                if mod_dispE.inc_x_trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc_x_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if ismember(mod_dispE.E_pred,[1,2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_dispE.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_perf',trialType_nm], bf_idx, sub_nm);
                end
                
                %% to lose trials
                trialType_nm = '_toLose';
                % onset effort scale
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_dispE',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % luminance (first)
                if ismember(mod_dispE.lum,[2,4,6,8])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (first)
                if ismember(mod_dispE.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % RT first press
                if ismember(mod_dispE.RT_fp,[2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_fp',trialType_nm], bf_idx, sub_nm);
                end

                % model variables
                if mod_dispE.mdl_n ~= 0
                    % ressource
                    if ismember(mod_dispE.E_pred,[5,6])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % reward displayed as piece of money or bill?
                if mod_dispE.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_dispE.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_dispE.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc_bis',trialType_nm], bf_idx, sub_nm);
                end

                % confidence
                if mod_dispE.conf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_conf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_dispE.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_GLM',mod_dispE.ROI_activity_GLM,...
                    '_',mod_dispE.ROI_activity_period,'_period_',mod_dispE.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % force exerted
                if mod_dispE.Effort ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_Effort',trialType_nm], bf_idx, sub_nm);
                end
                
                % model variables
                if mod_dispE.mdl_n ~= 0
                    
                    % ressource
                    if ismember(mod_dispE.E_pred,[3,4])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                
                    % cost
                    if mod_dispE.mdl_cost ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_cost',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % benefit
                    if mod_dispE.mdl_benefit ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_benefit',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % EV
                    if mod_dispE.mdl_EV ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_EV',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % trial number (middle)
                if ismember(mod_dispE.trialN,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT first press
                if ismember(mod_dispE.RT_fp,[1])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_fp',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                if ismember(mod_dispE.lum,[1,3,5,7])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (last)
                if ismember(mod_dispE.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive*trial number
                if mod_dispE.inc_x_trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc_x_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if ismember(mod_dispE.E_pred,[1,2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_dispE.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_perf',trialType_nm], bf_idx, sub_nm);
                end
                
        end
        
        %% start of effort
        o_perfE = gripRPprm.o_perfE;
        mod_perfE = gripRPprm.mod_perfE;
        switch o_perfE
            case 1
                % onset effort perf
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_perfE'], bf_idx, sub_nm);
                
                % modulators
                % Gain/Loss condition
                if mod_perfE.GL_cond ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_GL_cond'], bf_idx, sub_nm);
                end
                
                % luminance (first)
                
                % trial number (first)
                if ismember(mod_perfE.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_trialN'], bf_idx, sub_nm);
                end
                
                % reward displayed as piece of money or bill?
                if mod_perfE.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_R_type'], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_perfE.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_inc'], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_perfE.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_inc_bis'], bf_idx, sub_nm);
                end
                
                % force exerted
                if mod_perfE.Effort ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_Effort'], bf_idx, sub_nm);
                end
                
                % luminance (last)
                
                % trial number (last)
                if ismember(mod_perfE.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_trialN'], bf_idx, sub_nm);
                end
                
            case 2
                %% to gain trials
                trialType_nm = '_toGain';
                % onset effort perf
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_perfE',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % luminance (first)
                
                % trial number (first)
                if ismember(mod_perfE.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % reward displayed as piece of money or bill?
                if mod_perfE.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_perfE.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_perfE.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_inc_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % force exerted
                if mod_perfE.Effort ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_Effort',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                
                % trial number (last)
                if ismember(mod_perfE.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% to lose trials
                trialType_nm = '_toLose';
                % onset effort perf
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_perfE',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % luminance (first)
                
                % trial number (first)
                if ismember(mod_perfE.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % reward displayed as piece of money or bill?
                if mod_perfE.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_perfE.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_perfE.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_inc_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % force exerted
                if mod_perfE.Effort ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_Effort',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                
                % trial number (last)
                if ismember(mod_perfE.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
        end % effort performance
        
        %% feedback
        o_fbk = gripRPprm.o_fbk;
        mod_fbk = gripRPprm.mod_fbk;
        switch o_fbk
            case 1
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk'], bf_idx, sub_nm);
                
                % modulators
                % Gain/Loss condition
                if mod_fbk.GL_cond ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GL_cond'], bf_idx, sub_nm);
                end
                
                % luminance (first)
                
                % trial number (first)
                if ismember(mod_fbk.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN'], bf_idx, sub_nm);
                end
                
                % ressource
                if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_E_pred'], bf_idx, sub_nm);
                end

                % feedback current trial
                if mod_fbk.fbk ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk'], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain'], bf_idx, sub_nm);
                end
                
                % ressource
                if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[3,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_E_pred'], bf_idx, sub_nm);
                end
                
                % performance
                if mod_fbk.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_perf'], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity'], bf_idx, sub_nm);
                end
                
                % luminance (last)
                
                % trial number (last)
                if ismember(mod_fbk.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN'], bf_idx, sub_nm);
                end
                
            case 2
                %% to gain trials
                trialType_nm = '_toGain';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % luminance (first)
                
                % trial number (first)
                if ismember(mod_fbk.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % ressource
                if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback current trial
                if mod_fbk.fbk ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[3,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_fbk.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_perf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                
                % trial number (last)
                if ismember(mod_fbk.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% to lose trials
                trialType_nm = '_toLose';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % luminance (first)
                
                % trial number (first)
                if ismember(mod_fbk.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % ressource
                if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback current trial
                if mod_fbk.fbk ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[3,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_fbk.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_perf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                
                % trial number (last)
                if ismember(mod_fbk.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
        end % feedback
        
        %% cross
        o_cross = gripRPprm.o_cross;
        switch o_cross
            case 1
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_cross'], bf_idx, sub_nm);
        end
        
        %% missed trials
        o_missed_trials_dispE = gripRPprm.o_missed_trials_dispE;
        switch o_missed_trials_dispE
            case 1
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_missed_trials_dispE'], bf_idx, sub_nm);
        end
        
        
        
        %% stroop
    elseif ismember(iRun, runs.stroop)
        task_nm = 'stroop';
        task_id = 'S_';
        
        %% incentive display
        o_inc = stroopRPprm.o_inc;
        mod_inc = stroopRPprm.mod_inc;
        switch o_inc
            case 1 % all onsets grouped
                % onset incentive
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_inc'], bf_idx, sub_nm);
                
                % modulators
                % Gain/Loss condition
                if mod_inc.GL_cond ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_GL_cond'], bf_idx, sub_nm);
                end
                
                % luminance (first)
                if ismember(mod_inc.lum,[2,4,6,8])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_lum'], bf_idx, sub_nm);
                end
                
                % trial number (first)
                if ismember(mod_inc.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN'], bf_idx, sub_nm);
                end

                % RT first pair
                if ismember(mod_inc.RT_fp,[2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_fp'], bf_idx, sub_nm);
                end

                % model variables
                if mod_inc.mdl_n ~= 0
                    % ressource
                    if ismember(mod_inc.E_pred,[5,6])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred'], bf_idx, sub_nm);
                    end
                end
                
                % reward displayed as piece of money or bill?
                if mod_inc.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_R_type'], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_inc.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc'], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_inc.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc_bis'], bf_idx, sub_nm);
                end

                % incentive bis
                if mod_inc.conf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_conf'], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_inc.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_GLM',mod_inc.ROI_activity_GLM,...
                    '_',mod_inc.ROI_activity_period,'_period_',mod_inc.ROI_activity_ROI_nm,'_activity'], bf_idx, sub_nm);
                end
                
                % number of pairs solved
                if mod_inc.n_pairs_solved ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_n_pairs_solved'], bf_idx, sub_nm);
                end

                % number of incongruent pairs solved
                if mod_inc.n_incong ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_n_incong'], bf_idx, sub_nm);
                end
                
                % were errors made?
                if mod_inc.n_errors ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_n_errors'], bf_idx, sub_nm);
                end
                
                % model variables
                if mod_inc.mdl_n ~= 0
                    
                    % ressource
                    if ismember(mod_inc.E_pred,[3,4])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred'], bf_idx, sub_nm);
                    end
                
                    % cost
                    if mod_inc.mdl_cost ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_cost'], bf_idx, sub_nm);
                    end
                    
                    % benefit
                    if mod_inc.mdl_benefit ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_benefit'], bf_idx, sub_nm);
                    end
                    
                    % EV
                    if mod_inc.mdl_EV ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_EV'], bf_idx, sub_nm);
                    end
                end
                
                % trial number (middle)
                if ismember(mod_inc.trialN,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN'], bf_idx, sub_nm);
                end
                
                % RT first pair
                if ismember(mod_inc.RT_fp,[1])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_fp'], bf_idx, sub_nm);
                end
                
                % mean RT across pairs (except first pair)
                if mod_inc.RT_mRT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_mRT'], bf_idx, sub_nm);
                end
                
                % total gain until now
                if mod_inc.totalGain_prev ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_totalGain_prev'], bf_idx, sub_nm);
                end
                
                % sum(performance previous trials)
                if mod_inc.sumPerfPrev ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_sumPerfPrev'], bf_idx, sub_nm);
                end
                
                % luminance (last)
                if ismember(mod_inc.lum,[1,3,5,7])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_lum'], bf_idx, sub_nm);
                end
                
                % trial number (last)
                if ismember(mod_inc.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN'], bf_idx, sub_nm);
                end
                
                % incentive*trial number
                if mod_inc.inc_x_trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc_x_trialN'], bf_idx, sub_nm);
                end
                
                % ressource
                if ismember(mod_inc.E_pred,[1,2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred'], bf_idx, sub_nm);
                end
                
                % performance
                if mod_inc.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_perf'], bf_idx, sub_nm);
                end
                
            case 2 % to gain/to lose trials split
                %% to gain trials
                trialType_nm = '_toGain';
                % onset incentive
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_inc',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % luminance (first)
                if ismember(mod_inc.lum,[2,4,6,8])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (first)
                if ismember(mod_inc.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % RT first pair
                if ismember(mod_inc.RT_fp,[2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_fp',trialType_nm], bf_idx, sub_nm);
                end

                % model variables
                if mod_inc.mdl_n ~= 0
                    % ressource
                    if ismember(mod_inc.E_pred,[5,6])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % reward displayed as piece of money or bill?
                if mod_inc.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_inc.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_inc.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc_bis',trialType_nm], bf_idx, sub_nm);
                end

                % confidence
                if mod_inc.conf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_conf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_inc.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_GLM',mod_inc.ROI_activity_GLM,...
                    '_',mod_inc.ROI_activity_period,'_period_',mod_inc.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % number of pairs solved
                if mod_inc.n_pairs_solved ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_n_pairs_solved',trialType_nm], bf_idx, sub_nm);
                end

                % number of incongruent pairs solved
                if mod_inc.n_incong ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_n_incong',trialType_nm], bf_idx, sub_nm);
                end
                
                % were errors made?
                if mod_inc.n_errors ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_n_errors',trialType_nm], bf_idx, sub_nm);
                end
                
                % model variables
                if mod_inc.mdl_n ~= 0
                    
                    % ressource
                    if ismember(mod_inc.E_pred,[3,4])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % cost
                    if mod_inc.mdl_cost ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_cost',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % benefit
                    if mod_inc.mdl_benefit ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_benefit',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % EV
                    if mod_inc.mdl_EV ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_EV',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % trial number (middle)
                if ismember(mod_inc.trialN,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT first pair
                if ismember(mod_inc.RT_fp,[1])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_fp',trialType_nm], bf_idx, sub_nm);
                end
                
                % mean RT across pairs (except first pair)
                if mod_inc.RT_mRT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_mRT',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain until now
                if mod_inc.totalGain_prev ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_totalGain_prev',trialType_nm], bf_idx, sub_nm);
                end
                
                % sum(performance previous trials)
                if mod_inc.sumPerfPrev ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_sumPerfPrev',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                if ismember(mod_inc.lum,[1,3,5,7])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (last)
                if ismember(mod_inc.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive*trial number
                if mod_inc.inc_x_trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc_x_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if ismember(mod_inc.E_pred,[1,2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_inc.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_perf',trialType_nm], bf_idx, sub_nm);
                end
                
                %% to lose trials
                trialType_nm = '_toLose';
                % onset incentive
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_inc',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % luminance (first)
                if ismember(mod_inc.lum,[2,4,6,8])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (first)
                if ismember(mod_inc.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % RT first pair (first)
                if ismember(mod_inc.RT_fp,[2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_fp',trialType_nm], bf_idx, sub_nm);
                end

                % model variables
                if mod_inc.mdl_n ~= 0
                    % ressource
                    if ismember(mod_inc.E_pred,[5,6])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % reward displayed as piece of money or bill?
                if mod_inc.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_inc.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_inc.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc_bis',trialType_nm], bf_idx, sub_nm);
                end

                % confidence
                if mod_inc.conf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_conf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_inc.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_GLM',mod_inc.ROI_activity_GLM,...
                    '_',mod_inc.ROI_activity_period,'_period_',mod_inc.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % number of pairs solved
                if mod_inc.n_pairs_solved ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_n_pairs_solved',trialType_nm], bf_idx, sub_nm);
                end

                % number of incongruent pairs solved
                if mod_inc.n_incong ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_n_incong',trialType_nm], bf_idx, sub_nm);
                end
                
                % were errors made?
                if mod_inc.n_errors ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_n_errors',trialType_nm], bf_idx, sub_nm);
                end
                
                % model variables
                if mod_inc.mdl_n ~= 0
                    
                    % ressource
                    if ismember(mod_inc.E_pred,[3,4])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % cost
                    if mod_inc.mdl_cost ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_cost',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % benefit
                    if mod_inc.mdl_benefit ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_benefit',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % EV
                    if mod_inc.mdl_EV ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_EV',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % trial number (middle)
                if ismember(mod_inc.trialN,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT first pair
                if ismember(mod_inc.RT_fp,[1])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_fp',trialType_nm], bf_idx, sub_nm);
                end
                
                % mean RT across pairs (except first pair)
                if mod_inc.RT_mRT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_mRT',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain until now
                if mod_inc.totalGain_prev ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_totalGain_prev',trialType_nm], bf_idx, sub_nm);
                end
                
                % sum(performance previous trials)
                if mod_inc.sumPerfPrev ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_sumPerfPrev',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                if ismember(mod_inc.lum,[1,3,5,7])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (last)
                if ismember(mod_inc.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive*trial number
                if mod_inc.inc_x_trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc_x_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if ismember(mod_inc.E_pred,[1,2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_inc.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_perf',trialType_nm], bf_idx, sub_nm);
                end
                
            case 3 % no-error/error trials split
                %% no-error trials
                trialType_nm = '_noErrorTrials';
                % onset incentive
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_inc',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % Gain/Loss condition
                if mod_inc.GL_cond ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_GL_cond',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (first)
                if ismember(mod_inc.lum,[2,4,6,8])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (first)
                if ismember(mod_inc.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % RT first pair (first)
                if ismember(mod_inc.RT_fp,[2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_fp',trialType_nm], bf_idx, sub_nm);
                end

                % model variables
                if mod_inc.mdl_n ~= 0
                    % ressource
                    if ismember(mod_inc.E_pred,[5,6])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % reward displayed as piece of money or bill?
                if mod_inc.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_inc.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_inc.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc_bis',trialType_nm], bf_idx, sub_nm);
                end

                % confidence
                if mod_inc.conf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_conf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_inc.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_GLM',mod_inc.ROI_activity_GLM,...
                    '_',mod_inc.ROI_activity_period,'_period_',mod_inc.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % number of pairs solved
                if mod_inc.n_pairs_solved ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_n_pairs_solved',trialType_nm], bf_idx, sub_nm);
                end

                % number of incongruent pairs solved
                if mod_inc.n_incong ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_n_incong',trialType_nm], bf_idx, sub_nm);
                end
                
                % no error for no-error trials condition
                
                % model variables
                if mod_inc.mdl_n ~= 0
                    
                    % ressource
                    if ismember(mod_inc.E_pred,[3,4])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % cost
                    if mod_inc.mdl_cost ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_cost',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % benefit
                    if mod_inc.mdl_benefit ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_benefit',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % EV
                    if mod_inc.mdl_EV ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_EV',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % trial number (middle)
                if ismember(mod_inc.trialN,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT first pair
                if ismember(mod_inc.RT_fp,[1])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_fp',trialType_nm], bf_idx, sub_nm);
                end
                
                % mean RT across pairs (except first pair)
                if mod_inc.RT_mRT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_mRT',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain until now
                if mod_inc.totalGain_prev ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_totalGain_prev',trialType_nm], bf_idx, sub_nm);
                end
                
                % sum(performance previous trials)
                if mod_inc.sumPerfPrev ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_sumPerfPrev',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                if ismember(mod_inc.lum,[1,3,5,7])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (last)
                if ismember(mod_inc.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive*trial number
                if mod_inc.inc_x_trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc_x_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if ismember(mod_inc.E_pred,[1,2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_inc.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_perf',trialType_nm], bf_idx, sub_nm);
                end
                
                %% error trials
                trialType_nm = '_ErrorTrials';
                % onset incentive
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_inc',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % Gain/Loss condition
                if mod_inc.GL_cond ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_GL_cond',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (first)
                if ismember(mod_inc.lum,[2,4,6,8])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (first)
                if ismember(mod_inc.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % RT first pair (first)
                if ismember(mod_inc.RT_fp,[2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_fp',trialType_nm], bf_idx, sub_nm);
                end

                % model variables
                if mod_inc.mdl_n ~= 0
                    % ressource
                    if ismember(mod_inc.E_pred,[5,6])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % reward displayed as piece of money or bill?
                if mod_inc.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_inc.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_inc.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc_bis',trialType_nm], bf_idx, sub_nm);
                end

                % confidence
                if mod_inc.conf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_conf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_inc.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_GLM',mod_inc.ROI_activity_GLM,...
                    '_',mod_inc.ROI_activity_period,'_period_',mod_inc.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % number of pairs solved
                if mod_inc.n_pairs_solved ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_n_pairs_solved',trialType_nm], bf_idx, sub_nm);
                end

                % number of incongruent pairs solved
                if mod_inc.n_incong ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_n_incong',trialType_nm], bf_idx, sub_nm);
                end
                
                % were errors made?
                if mod_inc.n_errors ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_n_errors',trialType_nm], bf_idx, sub_nm);
                end
                
                % model variables
                if mod_inc.mdl_n ~= 0
                    
                    % ressource
                    if ismember(mod_inc.E_pred,[3,4])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                
                    % cost
                    if mod_inc.mdl_cost ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_cost',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % benefit
                    if mod_inc.mdl_benefit ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_benefit',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % EV
                    if mod_inc.mdl_EV ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_model',num2str(mod_inc.mdl_n),'_EV',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % trial number (middle)
                if ismember(mod_inc.trialN,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % RT first pair
                if ismember(mod_inc.RT_fp,[1])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_fp',trialType_nm], bf_idx, sub_nm);
                end
                
                % mean RT across pairs (except first pair)
                if mod_inc.RT_mRT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_RT_mRT',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain until now
                if mod_inc.totalGain_prev ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_totalGain_prev',trialType_nm], bf_idx, sub_nm);
                end
                
                % sum(performance previous trials)
                if mod_inc.sumPerfPrev ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_sumPerfPrev',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                if ismember(mod_inc.lum,[1,3,5,7])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (last)
                if ismember(mod_inc.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive*trial number
                if mod_inc.inc_x_trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_inc_x_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if ismember(mod_inc.E_pred,[1,2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_inc.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_inc_perf',trialType_nm], bf_idx, sub_nm);
                end
                
        end % incentive
        
        %% effort scale display
        o_dispE = stroopRPprm.o_dispE;
        mod_dispE = stroopRPprm.mod_dispE;
        switch o_dispE
            case 1
                % onset effort scale
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_dispE'], bf_idx, sub_nm);
                
                % modulators
                % Gain/Loss condition
                if mod_dispE.GL_cond ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_GL_cond'], bf_idx, sub_nm);
                end
                
                % luminance (first)
                if ismember(mod_dispE.lum,[2,4,6,8])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_lum'], bf_idx, sub_nm);
                end
                
                % trial number (first)
                if ismember(mod_dispE.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN'], bf_idx, sub_nm);
                end

                % RT first pair
                if ismember(mod_dispE.RT_fp,[2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_fp'], bf_idx, sub_nm);
                end

                % model variables
                if mod_dispE.mdl_n ~= 0
                    % ressource
                    if ismember(mod_dispE.E_pred,[5,6])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred'], bf_idx, sub_nm);
                    end
                end
                
                % reward displayed as piece of money or bill?
                if mod_dispE.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_R_type'], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_dispE.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc'], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_dispE.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc_bis'], bf_idx, sub_nm);
                end

                % confidence
                if mod_dispE.conf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_conf'], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_dispE.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_GLM',mod_dispE.ROI_activity_GLM,...
                    '_',mod_dispE.ROI_activity_period,'_period_',mod_dispE.ROI_activity_ROI_nm,'_activity'], bf_idx, sub_nm);
                end
                
                % number of pairs solved
                if mod_dispE.n_pairs_solved ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_pairs_solved'], bf_idx, sub_nm);
                end
                
                % number/percentage incongruent pairs solved
                if mod_dispE.n_incong ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_incong'], bf_idx, sub_nm);
                end
                
                % distance between numbers of the solved pairs
                if mod_dispE.n_dist ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_dist'], bf_idx, sub_nm);
                end
                
                % were errors made?
                if mod_dispE.n_errors ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_errors'], bf_idx, sub_nm);
                end
                
                % model variables
                if mod_dispE.mdl_n ~= 0
                    
                    % ressource
                    if ismember(mod_dispE.E_pred,[3,4])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred'], bf_idx, sub_nm);
                    end
                
                    % cost
                    if mod_dispE.mdl_cost ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_cost'], bf_idx, sub_nm);
                    end
                    
                    % benefit
                    if mod_dispE.mdl_benefit ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_benefit'], bf_idx, sub_nm);
                    end
                    
                    % EV
                    if mod_dispE.mdl_EV ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_EV'], bf_idx, sub_nm);
                    end
                end
                
                % RT first pair
                if ismember(mod_dispE.RT_fp,[1])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_fp'], bf_idx, sub_nm);
                end
                
                % mean RT across pairs (except first pair)
                if mod_dispE.RT_mRT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_mRT'], bf_idx, sub_nm);
                end
                
                % luminance (last)
                if ismember(mod_dispE.lum,[1,3,5,7])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_lum'], bf_idx, sub_nm);
                end
                
                % trial number (last)
                if ismember(mod_dispE.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN'], bf_idx, sub_nm);
                end
                
                % incentive*trial number
                if mod_dispE.inc_x_trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc_x_trialN'], bf_idx, sub_nm);
                end
                
                % ressource
                if ismember(mod_dispE.E_pred,[1,2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred'], bf_idx, sub_nm);
                end
                
                % performance
                if mod_dispE.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_perf'], bf_idx, sub_nm);
                end
                
            case 2 % to gain/to lose
                %% to gain trials
                trialType_nm = '_toGain';
                % onset effort scale
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_dispE',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % luminance (first)
                if ismember(mod_dispE.lum,[2,4,6,8])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (first)
                if ismember(mod_dispE.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % RT first pair (first)
                if ismember(mod_dispE.RT_fp,[2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_fp',trialType_nm], bf_idx, sub_nm);
                end
                
                % model variables
                if mod_dispE.mdl_n ~= 0
                    % ressource
                    if ismember(mod_dispE.E_pred,[5,6])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                end

                % reward displayed as piece of money or bill?
                if mod_dispE.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_dispE.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_dispE.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc_bis',trialType_nm], bf_idx, sub_nm);
                end

                % confidence
                if mod_dispE.conf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_conf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_dispE.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_GLM',mod_dispE.ROI_activity_GLM,...
                    '_',mod_dispE.ROI_activity_period,'_period_',mod_dispE.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % number of pairs solved
                if mod_dispE.n_pairs_solved ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_pairs_solved',trialType_nm], bf_idx, sub_nm);
                end
                
                % number/percentage incongruent pairs solved
                if mod_dispE.n_incong ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_incong',trialType_nm], bf_idx, sub_nm);
                end
                
                % distance between numbers of the solved pairs
                if mod_dispE.n_dist ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_dist',trialType_nm], bf_idx, sub_nm);
                end
                
                % were errors made?
                if mod_dispE.n_errors ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_errors',trialType_nm], bf_idx, sub_nm);
                end
                
                % model variables
                if mod_dispE.mdl_n ~= 0
                    
                    % ressource
                    if ismember(mod_dispE.E_pred,[3,4])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % cost
                    if mod_dispE.mdl_cost ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_cost',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % benefit
                    if mod_dispE.mdl_benefit ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_benefit',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % EV
                    if mod_dispE.mdl_EV ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_EV',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % RT first pair
                if ismember(mod_dispE.RT_fp,[1])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_fp',trialType_nm], bf_idx, sub_nm);
                end
                
                % mean RT across pairs (except first pair)
                if mod_dispE.RT_mRT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_mRT',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                if ismember(mod_dispE.lum,[1,3,5,7])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (last)
                if ismember(mod_dispE.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive*trial number
                if mod_dispE.inc_x_trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc_x_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if ismember(mod_dispE.E_pred,[1,2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_dispE.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_perf',trialType_nm], bf_idx, sub_nm);
                end
                
                %% to lose trials
                trialType_nm = '_toLose';
                % onset effort scale
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_dispE',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % luminance (first)
                if ismember(mod_dispE.lum,[2,4,6,8])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (first)
                if ismember(mod_dispE.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % RT first pair (first)
                if ismember(mod_dispE.RT_fp,[2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_fp',trialType_nm], bf_idx, sub_nm);
                end

                % model variables
                if mod_dispE.mdl_n ~= 0
                    % ressource
                    if ismember(mod_dispE.E_pred,[5,6])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % reward displayed as piece of money or bill?
                if mod_dispE.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_dispE.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_dispE.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc_bis',trialType_nm], bf_idx, sub_nm);
                end

                % confidence
                if mod_dispE.conf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_conf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_dispE.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_GLM',mod_dispE.ROI_activity_GLM,...
                    '_',mod_dispE.ROI_activity_period,'_period_',mod_dispE.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % number of pairs solved
                if mod_dispE.n_pairs_solved ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_pairs_solved',trialType_nm], bf_idx, sub_nm);
                end
                
                % number/percentage incongruent pairs solved
                if mod_dispE.n_incong ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_incong',trialType_nm], bf_idx, sub_nm);
                end
                
                % distance between numbers of the solved pairs
                if mod_dispE.n_dist ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_dist',trialType_nm], bf_idx, sub_nm);
                end
                
                % were errors made?
                if mod_dispE.n_errors ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_errors',trialType_nm], bf_idx, sub_nm);
                end
                
                % model variables
                if mod_dispE.mdl_n ~= 0
                    
                    % ressource
                    if ismember(mod_dispE.E_pred,[3,4])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % cost
                    if mod_dispE.mdl_cost ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_cost',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % benefit
                    if mod_dispE.mdl_benefit ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_benefit',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % EV
                    if mod_dispE.mdl_EV ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_EV',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % RT first pair
                if ismember(mod_dispE.RT_fp,[1])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_fp',trialType_nm], bf_idx, sub_nm);
                end
                
                % mean RT across pairs (except first pair)
                if mod_dispE.RT_mRT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_mRT',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                if ismember(mod_dispE.lum,[1,3,5,7])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (last)
                if ismember(mod_dispE.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive*trial number
                if mod_dispE.inc_x_trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc_x_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if ismember(mod_dispE.E_pred,[1,2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_dispE.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_perf',trialType_nm], bf_idx, sub_nm);
                end
                
            case 3 % no-error/error
                %% no-error trials
                trialType_nm = '_noErrorTrials';
                % onset effort scale
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_dispE',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % Gain/Loss condition
                if mod_dispE.GL_cond ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_GL_cond',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (first)
                if ismember(mod_dispE.lum,[2,4,6,8])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (first)
                if ismember(mod_dispE.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % RT first pair (first)
                if ismember(mod_dispE.RT_fp,[2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_fp',trialType_nm], bf_idx, sub_nm);
                end
                
                % model variables
                if mod_dispE.mdl_n ~= 0
                    % ressource
                    if ismember(mod_dispE.E_pred,[5,6])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                end

                % reward displayed as piece of money or bill?
                if mod_dispE.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_dispE.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_dispE.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc_bis',trialType_nm], bf_idx, sub_nm);
                end

                % confidence
                if mod_dispE.conf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_conf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_dispE.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_GLM',mod_dispE.ROI_activity_GLM,...
                    '_',mod_dispE.ROI_activity_period,'_period_',mod_dispE.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % number of pairs solved
                if mod_dispE.n_pairs_solved ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_pairs_solved',trialType_nm], bf_idx, sub_nm);
                end
                
                % number/percentage incongruent pairs solved
                if mod_dispE.n_incong ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_incong',trialType_nm], bf_idx, sub_nm);
                end
                
                % distance between numbers of the solved pairs
                if mod_dispE.n_dist ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_dist',trialType_nm], bf_idx, sub_nm);
                end
                
                % no error regressor for no-error trials condition
                
                % model variables
                if mod_dispE.mdl_n ~= 0
                    
                    % ressource
                    if ismember(mod_dispE.E_pred,[3,4])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % cost
                    if mod_dispE.mdl_cost ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_cost',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % benefit
                    if mod_dispE.mdl_benefit ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_benefit',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % EV
                    if mod_dispE.mdl_EV ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_EV',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % RT first pair
                if ismember(mod_dispE.RT_fp,[1])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_fp',trialType_nm], bf_idx, sub_nm);
                end
                
                % mean RT across pairs (except first pair)
                if mod_dispE.RT_mRT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_mRT',trialType_nm], bf_idx, sub_nm);
                endx_
                
                % luminance (last)
                if ismember(mod_dispE.lum,[1,3,5,7])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (last)
                if ismember(mod_dispE.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive*trial number
                if mod_dispE.inc_x_trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc_x_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if ismember(mod_dispE.E_pred,[1,2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_dispE.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_perf',trialType_nm], bf_idx, sub_nm);
                end
                
                %% error trials
                trialType_nm = '_ErrorTrials';
                % onset effort scale
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_dispE',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % Gain/Loss condition
                if mod_dispE.GL_cond ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_GL_cond',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (first)
                if ismember(mod_dispE.lum,[2,4,6,8])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (first)
                if ismember(mod_dispE.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % RT first pair (first)
                if ismember(mod_dispE.RT_fp,[2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_fp',trialType_nm], bf_idx, sub_nm);
                end
                
                % model variables
                if mod_dispE.mdl_n ~= 0
                    % ressource
                    if ismember(mod_dispE.E_pred,[5,6])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                end

                % reward displayed as piece of money or bill?
                if mod_dispE.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_dispE.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_dispE.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc_bis',trialType_nm], bf_idx, sub_nm);
                end

                % confidence
                if mod_dispE.conf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_conf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_dispE.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_GLM',mod_dispE.ROI_activity_GLM,...
                    '_',mod_dispE.ROI_activity_period,'_period_',mod_dispE.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % number of pairs solved
                if mod_dispE.n_pairs_solved ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_pairs_solved',trialType_nm], bf_idx, sub_nm);
                end
                
                % number/percentage incongruent pairs solved
                if mod_dispE.n_incong ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_incong',trialType_nm], bf_idx, sub_nm);
                end
                
                % distance between numbers of the solved pairs
                if mod_dispE.n_dist ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_dist',trialType_nm], bf_idx, sub_nm);
                end
                
                % were errors made?
                if mod_dispE.n_errors ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_n_errors',trialType_nm], bf_idx, sub_nm);
                end
                
                % model variables
                if mod_dispE.mdl_n ~= 0
                    
                    % ressource
                    if ismember(mod_dispE.E_pred,[3,4])
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % cost
                    if mod_dispE.mdl_cost ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_cost',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % benefit
                    if mod_dispE.mdl_benefit ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_benefit',trialType_nm], bf_idx, sub_nm);
                    end
                    
                    % EV
                    if mod_dispE.mdl_EV ~= 0
                        [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_model',num2str(mod_dispE.mdl_n),'_EV',trialType_nm], bf_idx, sub_nm);
                    end
                end
                
                % RT first pair
                if ismember(mod_dispE.RT_fp,[1])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_fp',trialType_nm], bf_idx, sub_nm);
                end
                
                % mean RT across pairs (except first pair)
                if mod_dispE.RT_mRT ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_RT_mRT',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                if ismember(mod_dispE.lum,[1,3,5,7])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_lum',trialType_nm], bf_idx, sub_nm);
                end
                
                % trial number (last)
                if ismember(mod_dispE.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive*trial number
                if mod_dispE.inc_x_trialN ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_inc_x_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if ismember(mod_dispE.E_pred,[1,2])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_dispE.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_dispE_perf',trialType_nm], bf_idx, sub_nm);
                end
        end
        
        %% start of effort
        o_perfE = stroopRPprm.o_perfE;
        mod_perfE = stroopRPprm.mod_perfE;
        switch o_perfE
            case 1
                % onset effort perf
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_perfE'], bf_idx, sub_nm);
                
                % modulators
                % Gain/Loss condition
                if mod_perfE.GL_cond ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_GL_cond'], bf_idx, sub_nm);
                end
                
                % luminance (first)
                
                % trial number (first)
                if ismember(mod_perfE.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_trialN'], bf_idx, sub_nm);
                end
                
                % reward displayed as piece of money or bill?
                if mod_perfE.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_R_type'], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_perfE.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_inc'], bf_idx, sub_nm);
                end
                
                % number of pairs solved
                if mod_perfE.n_pairs_solved ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_pairs_solved'], bf_idx, sub_nm);
                end
                
                % number/percentage incongruent pairs solved
                if mod_perfE.n_incong ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_incong'], bf_idx, sub_nm);
                end
                
                % distance between numbers of the solved pairs
                if mod_perfE.n_dist ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_dist'], bf_idx, sub_nm);
                end
                
                % were errors made?
                if mod_perfE.n_errors ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_errors'], bf_idx, sub_nm);
                end
                
                % luminance (last)
                
                % trial number (last)
                if ismember(mod_perfE.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_trialN'], bf_idx, sub_nm);
                end
                
            case 2 % to gain/to lose
                %% to gain trials
                trialType_nm = '_toGain';
                % onset effort perf
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_perfE',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % luminance (first)
                
                % trial number (first)
                if ismember(mod_perfE.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % reward displayed as piece of money or bill?
                if mod_perfE.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_perfE.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_perfE.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_inc_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % number of pairs solved
                if mod_perfE.n_pairs_solved ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_pairs_solved',trialType_nm], bf_idx, sub_nm);
                end
                
                % number/percentage incongruent pairs solved
                if mod_perfE.n_incong ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_incong',trialType_nm], bf_idx, sub_nm);
                end
                
                % distance between numbers of the solved pairs
                if mod_perfE.n_dist ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_dist',trialType_nm], bf_idx, sub_nm);
                end
                
                % were errors made?
                if mod_perfE.n_errors ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_errors',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                
                % trial number (last)
                if ismember(mod_perfE.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% to lose trials
                trialType_nm = '_toLose';
                % onset effort perf
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_perfE',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % luminance (first)
                
                % trial number (first)
                if ismember(mod_perfE.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % reward displayed as piece of money or bill?
                if mod_perfE.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_perfE.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_perfE.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_inc_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % number of pairs solved
                if mod_perfE.n_pairs_solved ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_pairs_solved',trialType_nm], bf_idx, sub_nm);
                end
                
                % number/percentage incongruent pairs solved
                if mod_perfE.n_incong ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_incong',trialType_nm], bf_idx, sub_nm);
                end
                
                % distance between numbers of the solved pairs
                if mod_perfE.n_dist ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_dist',trialType_nm], bf_idx, sub_nm);
                end
                
                % were errors made?
                if mod_perfE.n_errors ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_errors',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                
                % trial number (last)
                if ismember(mod_perfE.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
            case 3 % no-error/error trials
                %% no-error trials
                trialType_nm = '_noErrorTrials';
                % onset effort perf
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_perfE',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % Gain/Loss condition
                if mod_perfE.GL_cond ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_GL_cond',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (first)
                
                % trial number (first)
                if ismember(mod_perfE.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % reward displayed as piece of money or bill?
                if mod_perfE.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_perfE.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_perfE.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_inc_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % number of pairs solved
                if mod_perfE.n_pairs_solved ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_pairs_solved',trialType_nm], bf_idx, sub_nm);
                end
                
                % number/percentage incongruent pairs solved
                if mod_perfE.n_incong ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_incong',trialType_nm], bf_idx, sub_nm);
                end
                
                % distance between numbers of the solved pairs
                if mod_perfE.n_dist ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_dist',trialType_nm], bf_idx, sub_nm);
                end
                
                % no error regression in no-error trials condition
                
                % luminance (last)
                
                % trial number (last)
                if ismember(mod_perfE.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% error trials
                trialType_nm = '_ErrorTrials';
                % onset effort perf
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_perfE',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % Gain/Loss condition
                if mod_perfE.GL_cond ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_GL_cond',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (first)
                
                % trial number (first)
                if ismember(mod_perfE.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                % reward displayed as piece of money or bill?
                if mod_perfE.R_type ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_R_type',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive
                if mod_perfE.inc ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_inc',trialType_nm], bf_idx, sub_nm);
                end
                
                % incentive bis
                if mod_perfE.inc_bis ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_inc_bis',trialType_nm], bf_idx, sub_nm);
                end
                
                % number of pairs solved
                if mod_perfE.n_pairs_solved ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_pairs_solved',trialType_nm], bf_idx, sub_nm);
                end
                
                % number/percentage incongruent pairs solved
                if mod_perfE.n_incong ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_incong',trialType_nm], bf_idx, sub_nm);
                end
                
                % distance between numbers of the solved pairs
                if mod_perfE.n_dist ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_dist',trialType_nm], bf_idx, sub_nm);
                end
                
                % were errors made?
                if mod_perfE.n_errors ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_n_errors',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                
                % trial number (last)
                if ismember(mod_perfE.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_perfE_trialN',trialType_nm], bf_idx, sub_nm);
                end
        end
        
        %% feedback
        o_fbk = stroopRPprm.o_fbk;
        mod_fbk = stroopRPprm.mod_fbk;
        switch o_fbk
            case 1
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk'], bf_idx, sub_nm);
                
                % modulators
                % Gain/Loss condition
                if mod_fbk.GL_cond ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GL_cond'], bf_idx, sub_nm);
                end
                
                % luminance (first)
                
                % trial number (first)
                if ismember(mod_fbk.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN'], bf_idx, sub_nm);
                end

                % ressource
                if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_E_pred'], bf_idx, sub_nm);
                end
                
                % feedback current trial
                if mod_fbk.fbk ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk'], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain'], bf_idx, sub_nm);
                end
                
                % ressource
                if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[1,2,3,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_E_pred'], bf_idx, sub_nm);
                end
                
                % performance
                if mod_fbk.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_perf'], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity'], bf_idx, sub_nm);
                end
                
                % luminance (last)
                
                % trial number (last)
                if ismember(mod_fbk.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN'], bf_idx, sub_nm);
                end
                
            case 2 % to gain/to lose
                %% to gain trials
                trialType_nm = '_toGain';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % luminance (first)
                
                % trial number (first)
                if ismember(mod_fbk.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % ressource
                if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback current trial
                if mod_fbk.fbk ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[1,2,3,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_fbk.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_perf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                
                % trial number (last)
                if ismember(mod_fbk.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% to lose trials
                trialType_nm = '_toLose';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % luminance (first)
                
                % trial number (first)
                if ismember(mod_fbk.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % ressource
                if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback current trial
                if mod_fbk.fbk ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[1,2,3,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_fbk.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_perf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                
                % trial number (last)
                if ismember(mod_fbk.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
            case 3 % no-error/error
                %% no-error trials
                trialType_nm = '_noErrorTrials';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % Gain/Loss condition
                if mod_fbk.GL_cond ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GL_cond',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (first)
                
                % trial number (first)
                if ismember(mod_fbk.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % ressource
                if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback current trial
                if mod_fbk.fbk ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[1,2,3,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_fbk.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_perf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                
                % trial number (last)
                if ismember(mod_fbk.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
                
                %% error trials
                trialType_nm = '_ErrorTrials';
                % onset feedback
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx, sub_nm);
                
                % modulators
                % Gain/Loss condition
                if mod_fbk.GL_cond ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GL_cond',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (first)
                
                % trial number (first)
                if ismember(mod_fbk.trialN,[2,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end

                % ressource
                if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[5,6])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % feedback current trial
                if mod_fbk.fbk ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx, sub_nm);
                end
                
                % total gain
                if mod_fbk.totalGain ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx, sub_nm);
                end
                
                % ressource
                if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[1,2,3,4])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_E_pred',trialType_nm], bf_idx, sub_nm);
                end
                
                % performance
                if mod_fbk.perf ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_perf',trialType_nm], bf_idx, sub_nm);
                end
                
                % ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx, sub_nm);
                end
                
                % luminance (last)
                
                % trial number (last)
                if ismember(mod_fbk.trialN,[1,3])
                    [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx, sub_nm);
                end
        end % feedback
        
        %% cross
        o_cross = stroopRPprm.o_cross;
        switch o_cross
            case 1
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_cross'], bf_idx, sub_nm);
        end
        
        %% missed trials
        o_missed_trials_dispE = stroopRPprm.o_missed_trials_dispE;
        switch o_missed_trials_dispE
            case 1
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_missed_trials_dispE'], bf_idx, sub_nm);
        end
        
    else
        error('should not have happened, check runs');
    end % task filter
    
    %% add movement parameters at the end of each run
    n_prm.all = n_prm.all + n_mvmt_perRun(iRun);
    n_prm.(task_nm).(run_nm) = n_prm.(task_nm).(run_nm) + n_mvmt_perRun(iRun);
    bf_idx = bf_idx + n_mvmt_perRun(iRun);
    
end % run loop

%% add one constant per task at the end
n_prm.all = n_prm.all + n_total_runs;

end % function end


function[n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, curr_reg_nm, bf_idx, sub_nm)
%[n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, curr_reg_nm, bf_idx, sub_nm)
% sub-function to update each regressor with one single line

n_prm.all = n_prm.all + nBperR; % update total number of parameters

n_prm.(task_nm).(run_nm) = n_prm.(task_nm).(run_nm) + nBperR; % update RL number of parameters for the current run

% extract index
% if ( strcmp(sub_nm,'s1_030317') && ismember(run_nm,{'run3','run5','run6'}) ) ||...
%         ( strcmp(sub_nm,'s2_030317') && ismember(run_nm,{'run5'}) ) ||...
%         ( strcmp(sub_nm,'s3_040317') && ismember(run_nm,{'run3','run6'}) ) ||...
%         ( strcmp(sub_nm,'s4_040317') && ismember(run_nm,{'run2','run3','run4','run5','run6','run7'}) ) ||...
%         ( strcmp(sub_nm,'s6_050317') && ismember(run_nm,{'run5','run7'}) ) ||...
%         ( strcmp(sub_nm,'s9_080317') && ismember(run_nm,{'run3','run5','run6'}) ) ||...
%         ( strcmp(sub_nm,'s22_080317') && ismember(run_nm,{'run5'}) ) % runs with too much movement => not too include in contrast
% else
    curr_reg_idx = bf_idx + 1;
    if ~isfield(prm_idx, curr_reg_nm)
        prm_idx.(curr_reg_nm) = curr_reg_idx; % extract index for this regressor
    else % means that the current regressor has already been extracted in a previous run for example
        % => extract all the index when it is used inside the same variable
        prm_idx.(curr_reg_nm) = [prm_idx.(curr_reg_nm), curr_reg_idx];
    end
% end

bf_idx = bf_idx + nBperR; % update index (need to be done after extracting current regressor index)

end
