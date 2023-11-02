function [  ] = MS2_pupil_RL_Val_Conf_DT(  )
% MS2_pupil_RL_Val_Conf_DT will compute the correlation between value (Val),
% confidence (Conf) and deliberation time (DT) and the pupil diameter 
% in the learning (RL) task. The pupil data needs to have been preprocessed
% beforehand to be able to launch this script.

%% set main directory and list of subjects
subject_id = {'enter list of subjects here'};
NS = length(subject_id);

%% working directory
root = 'enter path where main data is here';
% add path where luminance data lies
lumPath = 'enter path where luminance was extracted';
addpath(lumPath);

%% check if VBA toolbox is in the path already or not
% if the VBA toolbox is not installed yet, please do it otherwise the
% script may not apply the RFT
wherever = pwd;
if ~exist('RFT_GLM_contrast.m','file')
    error(['Please install the VBA toolbox at ',...
        'https://mbb-team.github.io/VBA-toolbox/download/ to be able to ',...
        'apply the RFT on the pupil data.']);
end

%% initialize variables of interest
nTrials_per_run = 60;
nRuns                   = 3; % 3 RL runs
nTotal_trials           = nTrials_per_run*nRuns;

% time to check
% onset-lock
time_bf_onset   = 500; % start checking how many ms before the onset
time_af_RT      = 500; % end checking how many ms after the reaction time
n_effort_time   = 3000;
nTotalTimeChecking = n_effort_time + time_bf_onset + time_af_RT;
% choice-lock
time_bf_choice = 3000;
time_af_choice = 3000;
nTotalTimeChecking_choice = time_bf_choice + time_af_choice;

eye_data_type = {'Xcoord','Ycoord','pupilD'};
n_eye_data_type = length(eye_data_type);
for iEye_type = 1:n_eye_data_type
    eye_type_nm = eye_data_type{iEye_type};
    [eyeXYpupil.onset_lock.(eye_type_nm)] = deal(NaN(nTotal_trials, nTotalTimeChecking, NS));
    [eyeXYpupil.choice_lock.(eye_type_nm)] = deal(NaN(nTotal_trials, nTotalTimeChecking_choice, NS));
    % median split data
    [eyeXYpupil.medSplit.onset_lock.Val.low.(eye_type_nm),...
        eyeXYpupil.medSplit.onset_lock.Val.high.(eye_type_nm),...
        eyeXYpupil.medSplit.onset_lock.Conf.low.(eye_type_nm),...
        eyeXYpupil.medSplit.onset_lock.Conf.high.(eye_type_nm),...
        eyeXYpupil.medSplit.onset_lock.RT.low.(eye_type_nm),...
        eyeXYpupil.medSplit.onset_lock.RT.high.(eye_type_nm)] = deal(NaN(nTotal_trials/2, nTotalTimeChecking, NS));
    [eyeXYpupil.medSplit.choice_lock.Val.low.(eye_type_nm),...
        eyeXYpupil.medSplit.choice_lock.Val.high.(eye_type_nm),...
        eyeXYpupil.medSplit.choice_lock.Conf.low.(eye_type_nm),...
        eyeXYpupil.medSplit.choice_lock.Conf.high.(eye_type_nm),...
        eyeXYpupil.medSplit.choice_lock.RT.low.(eye_type_nm),...
        eyeXYpupil.medSplit.choice_lock.RT.high.(eye_type_nm)] = deal(NaN(nTotal_trials/2, nTotalTimeChecking_choice, NS));
end
regs = {'run1_cstt','run2_cstt','run3_cstt','Lum','Val','Conf','DT'};
nReg = length(regs);
for iReg = 1:nReg
    reg_nm = regs{iReg};
    [RL_regs.raw.(reg_nm),...
        RL_regs.zscore_acrossRuns.(reg_nm),...
        RL_reg_to_use.(reg_nm)] = deal( NaN(nTotal_trials, NS) );
end

% betas for the incentive
betas.onset_lock = NaN(NS, nReg, nTotalTimeChecking);
betas.choice_lock= NaN(NS, nReg, nTotalTimeChecking_choice); % betas regressors of interest

%% correct for baseline or not?
corr_base = true;
switch corr_base
    case true
        corr_base_nm = '_baseCorr';
    case false
        corr_base_nm = '';
end

%% use raw pupil data (just excluding out of screen gaze) (0) or preprocessed
% eye-data (1) ?
raw_or_preproc_data = true;
switch raw_or_preproc_data
    case false % add a name to identify data based on preprocessed files versus not
        raw_or_ppc = '_raw';
    case true
        raw_or_ppc = '_ppc_data';
end

%% use raw values of regressors or zscore across runs
raw_or_z_regs = true;
switch raw_or_z_regs
    case false
        RL_type_to_use = 'raw';
    case true
        RL_type_to_use = 'zscore_acrossRuns';
end

%% trials to focus on
curr_modStim_cond2 = 'GL_Pairs';

%% Q.learning model to use
mdl_type_toUse = 'Nico_VBA_models';
n_Qmodel_toUse = 6;

%% load luminance
lum = getfield(load([lumPath,filesep,...
    'MS2_RL_luminance_perStim.mat']),'lum');

%% load excluded runs
list_skipped_files = getfield( load([root 'behavior_summary',...
    filesep 'eye' filesep 'list_skipped_runs.mat']),'list_skipped_files');

%% loop through subjects
for iSub = 1:NS
    %% subject id and paths
    sub_nm = subject_id{iSub};
    if strcmp(sub_nm(3),'_')
        subid   = sub_nm(2);
    elseif ~strcmp(sub_nm(3),'_') && strcmp(sub_nm(4),'_')
        subid = sub_nm(2:3);
    end
    subj_folder                 = fullfile(root, sub_nm);
    subj_behavior_folder        = [fullfile(subj_folder, 'behavior'), filesep];
    subj_fMRI_analysis_folder   = [fullfile(subj_folder,'fMRI_analysis'), filesep];
    subj_eye_folder             = [fullfile(subj_folder, 'eye_analysis','preprocessed_files'), filesep];
    
    %% extract RL runs
    [RL_runs]     = task_runs_extraction('RL', sub_nm);
    
    %% check if some runs need to be excluded
    RL_runs_ok = [];
    for iRLRun = 1:nRuns
        eye_nm = ['MS2s',subid,'r',num2str(RL_runs(iRLRun)),'.asc'];
        if ~ismember(eye_nm, list_skipped_files)
            RL_runs_ok = [RL_runs_ok, iRLRun];
        end
    end
    
    if ~isempty(RL_runs_ok) % ignore subjects where no RL run is usable
        
        %% loop through RL runs
        for iRLRun = RL_runs_ok
            
            run_nm = num2str(RL_runs(iRLRun));
            sess_nm = ['session',num2str(iRLRun)];
            
            %% extract relevant variables
            RL_behav_loadStruct.(sess_nm)  = load([subj_behavior_folder,'global_sub_',subid,'_session_',run_nm,'_learning.mat'],...
                'onset','T0','npair','lottery','choice','response',...
                'feedback','gain','rt_fp','pairForThisRun','nstim');
            switch raw_or_preproc_data
                case false
                    error('raw data not ready yet');
                case true
                    RL_eye_loadStruct.(sess_nm)    = load([subj_eye_folder,'eye_preproc_',sub_nm,'_RL_run',run_nm,'.mat'],...
                        'X_eyeCoord','Y_eyeCoord','z_pupil_diam','pupil_diam','eyeTime');
            end

            % load Q.values
            RL_model_loadStruct.(sess_nm) = getfield(getfield( load([subj_fMRI_analysis_folder,...
                'onsets_sub',subid,'_learning_run',run_nm],'learn'), 'learn'),'mod');
            % loadpA*QA+pB*QB
            RL.(sess_nm).SV = RL_model_loadStruct.(sess_nm).(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.SV.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']);
            % load confidence
            RL.(sess_nm).Conf = (RL_model_loadStruct.(sess_nm).(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.pChoice.best.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']) - 0.5).^2;
            % trial number
            trialN = RL_model_loadStruct.(sess_nm).trialN.raw.(curr_modStim_cond2).main;

            % RL data
            % onsets
            RL.(sess_nm).onset     = RL_behav_loadStruct.(sess_nm).onset;
            % T0
            RL.(sess_nm).T0        = RL_behav_loadStruct.(sess_nm).T0;
            % X coordinates
            RL.(sess_nm).Xcoord    = RL_eye_loadStruct.(sess_nm).X_eyeCoord;
            % Y coordinates
            RL.(sess_nm).Ycoord    = RL_eye_loadStruct.(sess_nm).Y_eyeCoord;
            % pupil diameter
            %             RL.(sess_nm).pupilD    = RL_eye_loadStruct.(sess_nm).pupil_diam;
            RL.(sess_nm).pupilD    = RL_eye_loadStruct.(sess_nm).z_pupil_diam;
            % eye timing
            RL.(sess_nm).eyeTime   = RL_eye_loadStruct.(sess_nm).eyeTime;
            
            % behavior
            npair          = RL_behav_loadStruct.(sess_nm).npair;
            pairForThisRun = RL_behav_loadStruct.(sess_nm).pairForThisRun;
            nstim          = RL_behav_loadStruct.(sess_nm).nstim;
            RL.(sess_nm).lottery        = RL_behav_loadStruct.(sess_nm).lottery;
            RL.(sess_nm).choice         = RL_behav_loadStruct.(sess_nm).choice;
            RL.(sess_nm).response       = RL_behav_loadStruct.(sess_nm).response;
            RL.(sess_nm).feedback       = RL_behav_loadStruct.(sess_nm).feedback;
            RL.(sess_nm).gain           = RL_behav_loadStruct.(sess_nm).gain;
            RL.(sess_nm).DT             = RL_behav_loadStruct.(sess_nm).rt_fp;
            RL.(sess_nm).trial_nb       = 1:nTrials_per_run;

            % define luminance for each trial
            gain_lum_tmp = lum.(['Stim',num2str(pairForThisRun),...
                num2str(nstim(1))]);
            ntal_lum_tmp = lum.(['Stim',num2str(pairForThisRun),...
                num2str(nstim(2))]);
            loss_lum_tmp = lum.(['Stim',num2str(pairForThisRun),...
                num2str(nstim(3))]);
            RL.(sess_nm).Lum = gain_lum_tmp.*(npair == 1) +...
                ntal_lum_tmp.*(npair == 2) +...
                loss_lum_tmp.*(npair == 3);
            
            %% RL task
            T0 = RL.(sess_nm).T0;
            onset_inc   = RL.(sess_nm).onset.displayOptions - T0; % in seconds
            choice_inc  = RL.(sess_nm).onset.choice - T0; % in seconds
            X_eyeCoord  = RL.(sess_nm).Xcoord;
            Y_eyeCoord  = RL.(sess_nm).Ycoord;
            pupilD      = RL.(sess_nm).pupilD;
            eyeTime     = RL.(sess_nm).eyeTime;
            % dur_inc = duration.displayIncentive;
            
            for iTrial = 1:nTrials_per_run
                jTrial = nTrials_per_run*(iRLRun - 1) + iTrial;
                
                if ~isnan( onset_inc(iTrial) ) &&...
                        (npair(iTrial) ~= 2) &&...
                        (ismember(iTrial,trialN))
                    onset_time = floor(onset_inc(iTrial)*1000); % ms
                    choice_time = floor(choice_inc(iTrial)*1000);
                    start_extraction_onset = onset_time - time_bf_onset;
                    start_extraction_choice = choice_time - time_bf_choice;
                    
                    % onset-lock
                    for iTime = 1:nTotalTimeChecking
                        checkTime_onset = floor((start_extraction_onset + iTime)/1000); % ms to seconds
                        eyeXYpupil.onset_lock.Xcoord(jTrial,iTime,iSub) = X_eyeCoord(eyeTime == checkTime_onset);
                        eyeXYpupil.onset_lock.Ycoord(jTrial,iTime,iSub) = Y_eyeCoord(eyeTime == checkTime_onset);
                        eyeXYpupil.onset_lock.pupilD(jTrial,iTime,iSub) = pupilD(eyeTime == checkTime_onset);
                    end
                    
                    % choice-lock
                    for iTime = 1:nTotalTimeChecking_choice
                        checkTime_choice = floor((start_extraction_choice + iTime)/1000); % ms to seconds
                        eyeXYpupil.choice_lock.Xcoord(jTrial,iTime,iSub) = X_eyeCoord(eyeTime == checkTime_choice);
                        eyeXYpupil.choice_lock.Ycoord(jTrial,iTime,iSub) = Y_eyeCoord(eyeTime == checkTime_choice);
                        eyeXYpupil.choice_lock.pupilD(jTrial,iTime,iSub) = pupilD(eyeTime == checkTime_choice);
                    end
                    
                    % extract regressors across runs
                    switch iRLRun
                        case 1
                            RL_regs.raw.run1_cstt(jTrial, iSub) = 1;
                            RL_regs.raw.run2_cstt(jTrial, iSub) = 0;
                            RL_regs.raw.run3_cstt(jTrial, iSub) = 0;
                        case 2
                            RL_regs.raw.run1_cstt(jTrial, iSub) = 0;
                            RL_regs.raw.run2_cstt(jTrial, iSub) = 1;
                            RL_regs.raw.run3_cstt(jTrial, iSub) = 0;
                        case 3
                            RL_regs.raw.run1_cstt(jTrial, iSub) = 0;
                            RL_regs.raw.run2_cstt(jTrial, iSub) = 0;
                            RL_regs.raw.run3_cstt(jTrial, iSub) = 1;
                    end
                    RL_regs.raw.Lum(jTrial, iSub)       = RL.(sess_nm).Lum(iTrial);
                    RL_regs.raw.trial_nb(jTrial, iSub)  = RL.(sess_nm).trial_nb(iTrial);
                    RL_regs.raw.Val(jTrial, iSub)       = RL.(sess_nm).SV(trialN == iTrial);
                    RL_regs.raw.Conf(jTrial, iSub)      = RL.(sess_nm).Conf(trialN == iTrial);
                    RL_regs.raw.DT(jTrial, iSub)        = RL.(sess_nm).DT(iTrial);
                    
                    %% correct for baseline
                    if corr_base == true
                        % onset-lock
                        baseline_onset_idx_start  = find(eyeTime == floor(start_extraction_onset/1000));
                        baseline_onset_idx_end    = find(eyeTime == floor(onset_time/1000));
                        baseline_onset_pupil_for_this_trial = mean( pupilD(baseline_onset_idx_start:baseline_onset_idx_end) ,'omitnan'); % extract baseline for this trial
                        eyeXYpupil.onset_lock.pupilD(jTrial,:,iSub) = eyeXYpupil.onset_lock.pupilD(jTrial,:,iSub) - baseline_onset_pupil_for_this_trial; % correct the data for the baseline
                        
                        % choice-lock baseline
                        baseline_choice_idx_start  = find(eyeTime == floor(start_extraction_choice/1000));
                        baseline_choice_idx_end    = find(eyeTime == floor(choice_time/1000));
                        baseline_choice_pupil_for_this_trial = mean( pupilD(baseline_choice_idx_start:baseline_choice_idx_end) ,'omitnan'); % extract baseline for this trial
                        eyeXYpupil.choice_lock.pupilD(jTrial,:,iSub) = eyeXYpupil.choice_lock.pupilD(jTrial,:,iSub) - baseline_choice_pupil_for_this_trial; % correct the data for the baseline
                    end
                    
                end % trial filter (avoid too slow trials)
                
            end % trial loop
            
        end % run loop
        
        for iReg = 1:nReg
            reg_nm = regs{iReg};


            % ignore run constants from this procedure (zscoring wouldn't make sense)
            if ~ismember(reg_nm,{'run1_cstt','run2_cstt','run3_cstt'})
                %% zscore regressors of interest
                RL_regs.zscore_acrossRuns.(reg_nm)(:,iSub)    = nanzscore( RL_regs.raw.(reg_nm)(:,iSub) );

                %% select between raw or zscored variable which to use
                RL_reg_to_use.(reg_nm)(:,iSub) = RL_regs.(RL_type_to_use).(reg_nm)(:,iSub);
            else
                RL_reg_to_use.(reg_nm)(:,iSub) = RL_regs.raw.(reg_nm)(:,iSub);
            end
        end
        
        %% correlation pupil and regressors of interest
        
        % exclude NaN trials from the correlation
        ok_trials = ( (~isnan( RL_reg_to_use.run1_cstt(:,iSub) )).*...
            (~isnan( RL_reg_to_use.run2_cstt(:,iSub) ) ).*...
            (~isnan( RL_reg_to_use.run3_cstt(:,iSub) ) )) == 1;
        % correlation
        %  perform GLM
        x0_vars  = [RL_reg_to_use.run1_cstt(ok_trials,iSub),...
            RL_reg_to_use.run2_cstt(ok_trials,iSub),...
            RL_reg_to_use.run3_cstt(ok_trials,iSub),...
            RL_reg_to_use.Lum(ok_trials,iSub),...
            RL_reg_to_use.Val(ok_trials,iSub),...
            RL_reg_to_use.Conf(ok_trials,iSub),...
            RL_reg_to_use.DT(ok_trials,iSub)];
        % orthogonalize the variables in X_grip
        x_regs                  = mtrx_orthog(x0_vars);
        Y_pupilD_onset          = eyeXYpupil.onset_lock.pupilD(ok_trials,:,iSub); % trials in lines and time in columns
        Y_pupilD_choice         = eyeXYpupil.choice_lock.pupilD(ok_trials,:,iSub); % trials in lines and time in columns
        for iTimePoint = 1:nTotalTimeChecking
            [betas.onset_lock(iSub,:, iTimePoint), ~, stats_onset] =...
                glmfit(x_regs, Y_pupilD_onset(:,iTimePoint), 'normal',...
                'constant','off');
        end
        for iTimePoint_Choice = 1:nTotalTimeChecking_choice
            [betas.choice_lock(iSub,:, iTimePoint_Choice), ~, stats_choice] =...
                glmfit(x_regs, Y_pupilD_choice(:,iTimePoint_Choice), 'normal',...
                'constant','off');
        end
        
        %% median split based on each regressor (ignoring run constants)
        for iReg = 4:nReg
            reg_nm = regs{iReg};
            
            % median
            med_reg = median( RL_reg_to_use.(reg_nm)(:,iSub) ,'omitnan');
            % low trials
            low_trials = RL_reg_to_use.(reg_nm)(:,iSub) <= med_reg;
            n_low_trials = sum( low_trials );
            low_reg_trials = 1:n_low_trials;
            % high trials
            high_trials = RL_reg_to_use.(reg_nm)(:,iSub) >= med_reg;
            n_high_trials = sum( high_trials );
            high_reg_trials = 1:n_high_trials;
            
            for iEye_type = 1:n_eye_data_type
                eye_type_nm = eye_data_type{iEye_type};
                
                % low
                % onset_lock
                eyeXYpupil.onset_lock.(reg_nm).low.(eye_type_nm)(low_reg_trials,:,iSub) = eyeXYpupil.onset_lock.(eye_type_nm)(low_trials,:,iSub);
                % choice-lock
                eyeXYpupil.choice_lock.(reg_nm).low.(eye_type_nm)(low_reg_trials,:,iSub) = eyeXYpupil.choice_lock.(eye_type_nm)(low_trials,:,iSub);
                
                % high
                % onset-lock
                eyeXYpupil.onset_lock.(reg_nm).high.(eye_type_nm)(high_reg_trials,:,iSub) = eyeXYpupil.onset_lock.(eye_type_nm)(high_trials,:,iSub);
                % choice-lock
                eyeXYpupil.choice_lock.(reg_nm).high.(eye_type_nm)(high_reg_trials,:,iSub) = eyeXYpupil.choice_lock.(eye_type_nm)(high_trials,:,iSub);
            end % eye data type X/Y/pupilD
        end % regressor loop
        
        %% keep track of how many subjects done
        disp([sub_nm,' extracted']);
        
    else
        warning(['Subject ',sub_nm, 'skipped because no RL run was usable.']);
    end % ignore subjects where no run is good
    
end % subject loop

med_splt_cat = {'low','high'};
eye_lock_type = {'onset_lock','choice_lock'};
n_eye_lock_type = length(eye_lock_type);

%% mean +/-SEM across trials

% eye data
for iEye_type = 1:n_eye_data_type
    eye_type_nm = eye_data_type{iEye_type};
    eyeXYpupil.mean_per_sub.onset_lock.(eye_type_nm)    = mean(eyeXYpupil.onset_lock.(eye_type_nm),1,'omitnan');
    eyeXYpupil.mean_per_sub.choice_lock.(eye_type_nm)   = mean(eyeXYpupil.choice_lock.(eye_type_nm),1,'omitnan');
    
    % median split data
    for iMed_splt = 1:length(med_splt_cat)
        med_splt_nm = med_splt_cat{iMed_splt};
        
        for iReg = 4:nReg
            reg_nm = regs{iReg};
            eyeXYpupil.mean_per_sub.onset_lock.(reg_nm).(med_splt_nm).(eye_type_nm)     = mean(eyeXYpupil.onset_lock.(reg_nm).(med_splt_nm).(eye_type_nm),1,'omitnan');
            eyeXYpupil.mean_per_sub.choice_lock.(reg_nm).(med_splt_nm).(eye_type_nm)    = mean(eyeXYpupil.choice_lock.(reg_nm).(med_splt_nm).(eye_type_nm),1,'omitnan');
        end
    end
end

%% mean +/-SEM across subjects
% eye data
for iEye_type = 1:n_eye_data_type
    eye_type_nm = eye_data_type{iEye_type};
    
    % mean eye coords and pupil
    eyeXYpupil.mean.onset_lock.(eye_type_nm)(1,:)    = mean(eyeXYpupil.mean_per_sub.onset_lock.(eye_type_nm),3,'omitnan');
    eyeXYpupil.mean.choice_lock.(eye_type_nm)(1,:)    = mean(eyeXYpupil.mean_per_sub.choice_lock.(eye_type_nm),3,'omitnan');
    % SEM eye coords and pupil
    eyeXYpupil.sem.onset_lock.(eye_type_nm)(1,:)     = sem(eyeXYpupil.mean_per_sub.onset_lock.(eye_type_nm),3);
    eyeXYpupil.sem.choice_lock.(eye_type_nm)(1,:)     = sem(eyeXYpupil.mean_per_sub.choice_lock.(eye_type_nm),3);
    
    % median split data
    for iMed_splt = 1:length(med_splt_cat)
        med_splt_nm = med_splt_cat{iMed_splt};
        
        for iReg = 4:nReg
            reg_nm = regs{iReg};
            
            % mean
            eyeXYpupil.mean.onset_lock.(reg_nm).(med_splt_nm).(eye_type_nm)(1,:)   = mean(eyeXYpupil.mean_per_sub.onset_lock.(reg_nm).(med_splt_nm).(eye_type_nm),3,'omitnan');
            eyeXYpupil.mean.choice_lock.(reg_nm).(med_splt_nm).(eye_type_nm)(1,:)  = mean(eyeXYpupil.mean_per_sub.choice_lock.(reg_nm).(med_splt_nm).(eye_type_nm),3,'omitnan');
            % SEM
            eyeXYpupil.sem.onset_lock.(reg_nm).(med_splt_nm).(eye_type_nm)(1,:)   = sem(eyeXYpupil.mean_per_sub.onset_lock.(reg_nm).(med_splt_nm).(eye_type_nm),3);
            eyeXYpupil.sem.choice_lock.(reg_nm).(med_splt_nm).(eye_type_nm)(1,:)   = sem(eyeXYpupil.mean_per_sub.choice_lock.(reg_nm).(med_splt_nm).(eye_type_nm),3);
        end
    end
end

%% slight smooth of the betas to improve RFT application
smooth_kernel = 100; % short smoothing for betas
for iLock_type = 1:n_eye_lock_type
    eye_lock_nm = eye_lock_type{iLock_type};
    
    switch eye_lock_nm
        case 'onset_lock'
            nTime = nTotalTimeChecking;
        case 'choice_lock'
            nTime = nTotalTimeChecking_choice;
    end
    
    conv_betas.(eye_lock_nm) = NaN(NS, nReg, nTime);
    for iReg = 1:nReg
        for iS = 1:NS
            beta_for_smooth = NaN(1,nTime);
            beta_for_smooth(1,:) = betas.(eye_lock_nm)(iS,iReg,:);
             smoothed_beta = smooth_NC(beta_for_smooth, smooth_kernel);
             for iTime = 1:nTime
                 conv_betas.(eye_lock_nm)(iS,iReg,iTime) = smoothed_beta(iTime);
             end
        end
    end
end

%% extract mean betas after smoothing them
for iReg = 4:nReg % ignore constants
    reg_nm = regs{iReg};
    
    beta_to_check.onset_lock(1:NS, 1:nTotalTimeChecking) = conv_betas.onset_lock(:,iReg,:);
    beta_to_check.choice_lock(1:NS, 1:nTotalTimeChecking_choice) = conv_betas.choice_lock(:,iReg,:);
    
    for iLock = 1:n_eye_lock_type
        eye_lock_nm = eye_lock_type{iLock};
        
        % mean and SEM
        mBetas.(eye_lock_nm).(reg_nm)      = mean( beta_to_check.(eye_lock_nm), 1,'omitnan');
        semBetas.(eye_lock_nm).(reg_nm)    = sem( beta_to_check.(eye_lock_nm), 1);
    end % locking type loop
end % regressor loop

%% RFT to check for statistical differences
RFT_prm.Xvector = ones(NS,1);
RFT_prm.verbose = 0;
RFT_prm.con = 1;
RFT_prm.conType = 'F';
%         RFT_prm.time_limits = [];
RFT_prm.pValCorr_threshold  = 0.05;

for iReg = 4:nReg % ignore run constants
    reg_nm = regs{iReg};
    
    for iLock = 1:n_eye_lock_type
        eye_lock_nm = eye_lock_type{iLock};
        switch eye_lock_nm
            case 'onset_lock'
                nTime = nTotalTimeChecking;
            case 'choice_lock'
                nTime = nTotalTimeChecking_choice;
        end
        RFT_prm.subfield_nm = reg_nm;
        
         % betas
        data_to_test.(reg_nm) = NaN(NS, nTime);
        data_to_test.(reg_nm)(1:NS, 1:nTime)    = conv_betas.(eye_lock_nm)(1:NS, iReg, 1:nTime);
        [ pVal_X_tmp, ~ ] = RFT_extraction( data_to_test, RFT_prm );
        pval_X_signif_clusters.(eye_lock_nm).(reg_nm) = pVal_X_tmp.(reg_nm);
        % extract corresponding Y values for display on the mean
        pval_Y_signif_clusters.(eye_lock_nm).(reg_nm) = mBetas.(eye_lock_nm).(reg_nm)( pval_X_signif_clusters.(eye_lock_nm).(reg_nm) );
        
        % median split
        [data_to_test_mSplit.(reg_nm),...
            sm_data_to_test_mSplit.(reg_nm)] = deal(NaN(NS, nTime));
        for iS = 1:NS
            for iTime = 1:nTime
                data_to_test_mSplit.(reg_nm)(iS, iTime) = eyeXYpupil.mean_per_sub.(eye_lock_nm).(reg_nm).low.pupilD(1,iTime,iS) - eyeXYpupil.mean_per_sub.(eye_lock_nm).(reg_nm).high.pupilD(1,iTime,iS);
            end % time loop
            sm_data_to_test_mSplit.(reg_nm)(iS, :) = smooth_NC(data_to_test_mSplit.(reg_nm)(iS, :), smooth_kernel);
        end % subject loop
        [ pVal_X_mSplit_tmp, ~ ] = RFT_extraction( sm_data_to_test_mSplit, RFT_prm );
        pval_X_signif_mSplit.(eye_lock_nm).(reg_nm) = pVal_X_mSplit_tmp.(reg_nm);
        % extract corresponding Y values for display on the mean
        pval_Y_signif_mSplit.(eye_lock_nm).(reg_nm).low   = eyeXYpupil.mean.(eye_lock_nm).(reg_nm).low.pupilD( pval_X_signif_mSplit.(eye_lock_nm).(reg_nm) );
        pval_Y_signif_mSplit.(eye_lock_nm).(reg_nm).high  = eyeXYpupil.mean.(eye_lock_nm).(reg_nm).high.pupilD( pval_X_signif_mSplit.(eye_lock_nm).(reg_nm) );
        
    end % lock type
end % regressors

%% graphs
stimEnd_onset_lock     = 3*1000; % ms
colours = {'k','r','b','g','k','y'};
x_eye.onset_lock    = (-time_bf_onset):(nTotalTimeChecking - time_bf_onset - 1);
x_eye.choice_lock   = (-time_bf_choice):(nTotalTimeChecking_choice - time_bf_choice - 1);
pSize = 30; % size for the legends, xlabel, ylabel, etc.
lWidthTraits = 3;

%% display average pupil shape
fig();

% pupil shape
for iLock = 1:n_eye_lock_type
    eye_lock_nm = eye_lock_type{iLock};
    subplot(1,n_eye_lock_type,iLock)
    jbfill(x_eye.(eye_lock_nm),...
        eyeXYpupil.mean.(eye_lock_nm).pupilD + eyeXYpupil.sem.(eye_lock_nm).pupilD,...
        eyeXYpupil.mean.(eye_lock_nm).pupilD - eyeXYpupil.sem.(eye_lock_nm).pupilD,...
        eyeXYpupil.mean.(eye_lock_nm).pupilD,...
        'b');
    hold on;
    line([stimEnd_onset_lock stimEnd_onset_lock],ylim(),...
        'Color','k','LineWidth',lWidthTraits);
    line(xlim,[0 0],'Color','k','LineWidth',lWidthTraits);
    line([0 0],ylim(),'Color','k','LineWidth',lWidthTraits);
    xlim([x_eye.(eye_lock_nm)(1),x_eye.(eye_lock_nm)(end)]);
    switch eye_lock_nm
        case 'onset_lock'
            xlabel('Time after stimulus onset (ms)');
        case 'choice_lock'
            xlabel('Time after choice (ms)');
    end
    ylabel('Pupil diameter');
    
    legend_size(pSize);
    line([0 0],ylim(),'Color','k','LineWidth',lWidthTraits);
end


%% display betas
for iLock = 1:n_eye_lock_type
    eye_lock_nm = eye_lock_type{iLock};
    
    fig_hdl = fig();
    
    % plot betas
    legend_list = [];
    for iReg = 4:nReg
        jRegCol = iReg - 3;
        reg_nm = regs{iReg};
        % show average beta
        handle_curve.(reg_nm) = jbfill(x_eye.(eye_lock_nm),...
            mBetas.(eye_lock_nm).(reg_nm) + semBetas.(eye_lock_nm).(reg_nm),...
            mBetas.(eye_lock_nm).(reg_nm) - semBetas.(eye_lock_nm).(reg_nm),...
            mBetas.(eye_lock_nm).(reg_nm),...
            colours{jRegCol});
        hold on;
        legend_list = [legend_list, handle_curve.(reg_nm)];
        
        % mark significant clusters from the RFT analysis
        x_signif = x_eye.(eye_lock_nm)( pval_X_signif_clusters.(eye_lock_nm).(reg_nm) );
        % mark the mean
        plot(x_signif,...
            pval_Y_signif_clusters.(eye_lock_nm).(reg_nm),...
            [colours{jRegCol},' *']);
        % mark significant clusters on the top of the graph
        ypos_perc = (1/100)*iReg;
        place_signif_curve( fig_hdl, x_signif, colours{jRegCol}, ypos_perc  );
    end
    
    xlim([x_eye.(eye_lock_nm)(1),x_eye.(eye_lock_nm)(end)]);
    line([0 0],ylim(),'Color','k','LineWidth',3);
    if strcmp(eye_lock_nm,'onset_lock')
        line([stimEnd_onset_lock stimEnd_onset_lock],ylim(),...
            'Color','k','LineWidth',3);
    end
    line(xlim(),[0 0]);
    switch eye_lock_nm
        case {'onset_lock'}
            xlabel('Time after stimulus onset (ms)');
        case {'choice_lock'}
            xlabel('Time after choice (ms)');
    end
    ylabel('Parameter estimates (a.u.)');
    legend_size(pSize);
    % legend
    legend(legend_list,regs(4:end),...
        'Location','southwest','fontsize',16,'FontWeight','bold');
    legend('boxoff');
end


%% median splits
for iReg = 4:nReg
    reg_nm = regs{iReg};
    
    fig_hdl = fig();
    
    for iLock = 1:n_eye_lock_type
        eye_lock_nm = eye_lock_type{iLock};
        
        subplot(1,n_eye_lock_type,iLock);
        % low
        low_reg_curve = jbfill(x_eye.(eye_lock_nm),...
            eyeXYpupil.mean.(eye_lock_nm).(reg_nm).low.pupilD + eyeXYpupil.sem.(eye_lock_nm).(reg_nm).low.pupilD,...
            eyeXYpupil.mean.(eye_lock_nm).(reg_nm).low.pupilD - eyeXYpupil.sem.(eye_lock_nm).(reg_nm).low.pupilD,...
            eyeXYpupil.mean.(eye_lock_nm).(reg_nm).low.pupilD,...
            'b');
        % high
        high_reg_curve = jbfill(x_eye.(eye_lock_nm),...
            eyeXYpupil.mean.(eye_lock_nm).(reg_nm).high.pupilD + eyeXYpupil.sem.(eye_lock_nm).(reg_nm).high.pupilD,...
            eyeXYpupil.mean.(eye_lock_nm).(reg_nm).high.pupilD - eyeXYpupil.sem.(eye_lock_nm).(reg_nm).high.pupilD,...
            eyeXYpupil.mean.(eye_lock_nm).(reg_nm).high.pupilD,...
            'r');
        
        xlim([x_eye.(eye_lock_nm)(1), x_eye.(eye_lock_nm)(end)]);
        
        % mark significant clusters from the RFT analysis
        x_signif = x_eye.(eye_lock_nm)( pval_X_signif_mSplit.(eye_lock_nm).(reg_nm) );
        % mark the mean
        plot(x_signif,...
            pval_Y_signif_mSplit.(eye_lock_nm).(reg_nm).low,['b',' *']);
        plot(x_signif,...
            pval_Y_signif_mSplit.(eye_lock_nm).(reg_nm).high,['r',' *']);
        % mark significant clusters on the top of the graph
        ypos_perc = 1/100;
        place_signif_curve( fig_hdl, x_signif, 'g', ypos_perc  );
        
        line([0 0],ylim(),'Color','k','LineWidth',3);
        line(xlim(),[0 0],'Color','b','LineWidth',3);
        if strcmp(eye_lock_nm,'onset_lock')
            line([stimEnd_onset_lock stimEnd_onset_lock],ylim(),'Color','k','LineWidth',3);
        end
        legend([high_reg_curve, low_reg_curve],{['high ',reg_nm],['low ',reg_nm]},...
            'Location','southwest','fontsize',16,'FontWeight','bold');
        legend('boxoff');
        switch eye_lock_nm
            case {'onset_lock'}
                xlabel('Time after stimulus onset (ms)');
            case {'choice_lock'}
                xlabel('Time after choice (ms)');
        end
        ylabel('Pupil diameter');
        legend_size(pSize);
    end % lock type loop
end % regressor loop

end % function end