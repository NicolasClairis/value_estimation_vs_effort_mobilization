function[eyeXYpupil, eyeXYpupil_mSplit] = MS2_pupil_GS_X_I_RT(task_id)
%[eyeXYpupil, eyeXYpupil_mSplit] = MS2_eye_GS_GLM_bis(task_id)
% MS2_eye_GS_GLM_bis performs a GLM and a median split on the pupil diameter
% according to several regressors of interest in the grip ('G') or the
% stroop ('S') tasks.
% Test performed with data locked on the onset of incentives, of the effort
% scale onset and of the effort performance onset.
%
%
% INPUTS
% task_id: 'G'/'S' = grip or stroop task
%
% OUTPUTS
% eyeXYpupil: structure with pupil data organized by trials and subjects
% and averaged
%
% eyeXYpupil_mSplit: structure with pupil data median split by each
% regressor of interest (before and after correcting pupil for brightness
% effects)

%% task id
switch task_id
    case 'G'
        task_nm = 'grip';
    case 'S'
        task_nm = 'stroop';
end

%% working directories
root = 'enter path here';
% add path where luminance data lies
addpath(fullfile('enter luminance path here'));

%% subject identification
subject_id = {'enter list of subjects here'};
NS = length(subject_id);
%%

%% check if VBA toolbox is in the path already or not
% if the VBA toolbox is not installed yet, installs it so that the RFT can
% be made
if ~exist('RFT_GLM_contrast.m','file')
    error(['Please install the VBA toolbox at ',...
        'https://mbb-team.github.io/VBA-toolbox/download/ to be able to ',...
        'apply the RFT on the pupil data.']);
end

%% initialize variables of interest
regs_names = {'run1_cstt',...
    'run2_cstt',...
    'Lum',...
    'RT','E'};
nReg = length(regs_names);

% eye splitting types
med_splt_cat = {'low','high'};
n_mSplit = length(med_splt_cat);
eye_data_type = {'pupilD','pupilD_noLum'};
n_eyeDataType = length(eye_data_type);

% number of trials
nTrials_per_run         = 60;
nRuns                   = 2;
nTotal_trials           = nTrials_per_run*nRuns;
% time to check
time_bf_onset   = 2000; % start checking how many ms before the onset
time_af_RT      = 2000; % end checking how many ms after the reaction time
switch task_id
    case 'G'
        n_effort_time   = 5000; % fixed time during which effort is performed (for dispE) and use same for inc period
    case 'S'
        n_effort_time   = 8000; % fixed time during which effort is performed (for dispE) and use same for inc period
end
nTotalTimeChecking      = n_effort_time + time_bf_onset + time_af_RT;
for iEye_type = 1:n_eyeDataType
    eye_type_nm = eye_data_type{iEye_type};
    [eyeXYpupil.inc.(eye_type_nm),...
        eyeXYpupil.dispE.(eye_type_nm),...
        eyeXYpupil.Eperf.(eye_type_nm)] = deal( NaN(nTotal_trials, nTotalTimeChecking, NS));
end

% median split eye data
for iMsplit_type = 1:n_mSplit
    mSplit_nm = med_splt_cat{iMsplit_type};
    for iEye_type = 1:n_eyeDataType
        eye_type_nm = eye_data_type{iEye_type};
        for iReg = 1:nReg
            reg_nm = regs_names{iReg};
            eyeXYpupil_mSplit.inc.(reg_nm).(mSplit_nm).(eye_type_nm) = NaN(nTotal_trials/2, nTotalTimeChecking, NS);
            eyeXYpupil_mSplit.dispE.(reg_nm).(mSplit_nm).(eye_type_nm) = NaN(nTotal_trials/2, nTotalTimeChecking, NS);
            eyeXYpupil_mSplit.Eperf.(reg_nm).(mSplit_nm).(eye_type_nm) = NaN(nTotal_trials/2, nTotalTimeChecking, NS);
            eyeXYpupil_mSplit.inc.(reg_nm).mean_per_sub.(mSplit_nm).(eye_type_nm) = NaN(nTotalTimeChecking, NS);
            eyeXYpupil_mSplit.dispE.(reg_nm).mean_per_sub.(mSplit_nm).(eye_type_nm) = NaN(nTotalTimeChecking, NS);
            eyeXYpupil_mSplit.Eperf.(reg_nm).mean_per_sub.(mSplit_nm).(eye_type_nm) = NaN(nTotalTimeChecking, NS);
        end
    end
end

% betas
for iReg = 1:nReg
    reg_nm = regs_names{iReg};
    [betas.inc.(reg_nm),...
        betas.dispE.(reg_nm),...
        betas.Eperf.(reg_nm),...
        betas_mSplit.inc.(reg_nm),...
        betas_mSplit.dispE.(reg_nm),...
        betas_mSplit.Eperf.(reg_nm)]    = deal( NaN(nTotalTimeChecking, NS) );
end

% actual live effort (force for grip and RT for response for Stroop)
switch task_id
    case 'G'
        n_grip_samples_per_trial = 300;
        [gripForce.dispE,...
            gripForceDisplay.dispE,...
            gripForceTime.inc,...
            gripForceTime.dispE,...
            gripForceTime.Eperf] = deal(NaN(nTotal_trials, n_grip_samples_per_trial, NS));
        [gripForce_perSub.dispE,...
            gripForceDisplay_perSub.dispE,...
            gripForceTime_perSub.inc,...
            gripForceTime_perSub.dispE,...
            gripForceTime_perSub.Eperf] = deal(NaN(n_grip_samples_per_trial, NS));
    case 'S'
        nMaxStroop = 10;
        stroopRT = NaN(nTotal_trials, nMaxStroop, NS);
        [inc_dur, E_perf_onset] = deal(NaN(nTotal_trials, NS));
        stroopRT_perSub = NaN(nMaxStroop, NS);
        [mTime_perPair_perSub_perTrial.inc,...
            mTime_perPair_perSub_perTrial.dispE,...
            mTime_perPair_perSub_perTrial.Eperf] = deal(NaN(nTotal_trials, nMaxStroop, NS));
        [mTime_perPair_perSub.inc,...
            mTime_perPair_perSub.dispE,...
            mTime_perPair_perSub.Eperf] = deal(NaN(nMaxStroop,NS));
end

%% correct for baseline or not?
corr_base = true;
switch corr_base
    case true
        corr_base_nm = '_baseCorr';
    case false
        corr_base_nm = '';
end

%% use raw data (just excluding out of screen gaze) (0) or preprocessed
% eye-data (1) ?
raw_or_preproc_data = 1;
switch raw_or_preproc_data
    case 0 % add a name to identify data based on preprocessed files versus not
        raw_or_ppc = '_raw';
    case 1
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

%% which model to use
n_model_toUse = 96;
model_nm = (['model_',num2str(n_model_toUse)]);

%% load excluded runs
list_skipped_files = getfield( load([root 'behavior_summary' filesep,...
    'eye' filesep 'list_skipped_runs.mat']),'list_skipped_files');

%% loop through subjects
for iS = 1:NS
    %% subject id and paths
    sub_nm = subject_id{iS};
    if strcmp(sub_nm(3),'_')
        subid   = sub_nm(2);
    elseif ~strcmp(sub_nm(3),'_') && strcmp(sub_nm(4),'_')
        subid = sub_nm(2:3);
    end
    subj_folder             = fullfile(root, sub_nm);
    sub_onsets_folder       = [fullfile(subj_folder, 'fMRI_analysis'), filesep];
    subj_eye_folder         = [fullfile(subj_folder, 'eye_analysis','preprocessed_files'), filesep];
    sub_bhv_folder          = [fullfile(subj_folder, 'behavior'), filesep];
    
    %% extract runs
    [task_runs]     = MS2_task_runs_extraction(task_nm, sub_nm);
    
    %% check if some runs need to be excluded
    task_runs_ok = [];
    for iRun = task_runs
        eye_nm = ['MS2s',subid,'r',num2str(iRun),'.asc'];
        if ~ismember(eye_nm, list_skipped_files)
            task_runs_ok = [task_runs_ok, iRun];
        end
    end
    
    %% initialize vars of interest
    [run1_cstt_sub,...
        run2_cstt_sub,...
        lum_inc_sub,...
        trialN_sub,...
        onset_inc_sub,...
        onset_dispE_sub,...
        onset_Eperf_sub,...
        X_sub,...
        RT_sub] = deal( NaN(nTotal_trials,1) );
    %% loop through runs
    jRun = 0;
    if ~isempty(task_runs_ok)
        for iRun = task_runs_ok
            jRun = jRun + 1;
            run_nm = num2str(iRun);
            
            %% load data
            % load behavioral vars
            behavioralData = getfield( load([sub_onsets_folder,'onsets_sub',subid,'_',task_nm,'_run',run_nm,'.mat'],task_nm),task_nm);
            % load model data
            loadFit = load([sub_onsets_folder,'GS_model_sub',subid,'_',task_nm,'_run',run_nm,'.mat']);
            X_pred_tmp = loadFit.X_pred.(model_nm);
            % load eye data
            switch raw_or_preproc_data
                case 0
                    error('raw data not ready yet');
                case 1
                    eye_loadStruct    = load([subj_eye_folder,'eye_preproc_',sub_nm,'_',task_nm,'_run',run_nm,'.mat'],...
                        'X_eyeCoord','Y_eyeCoord','z_pupil_diam','pupil_diam','eyeTime');
            end
            % load data to extract force and stroop RT
            switch task_id
                case 'G'
                    grip_loadStruct_tmp  = load([sub_bhv_folder,'global_sub_',subid,'_session_',run_nm,'_gripRP.mat'],...
                        'onset','T0','gripData','Fmax');
                    gripB_tmp     = grip_loadStruct_tmp.gripData.gripB;
                    gripDisplay_tmp     = grip_loadStruct_tmp.gripData.level;
                    grip_time_tmp = grip_loadStruct_tmp.gripData.time;
                    grip_onsets_time_tmp = grip_loadStruct_tmp.onset;
                    Fmax_sub = max(grip_loadStruct_tmp.Fmax,[],2,'omitnan');
                    T0 = grip_loadStruct_tmp.T0;
                case 'S'
                    inc_dur_tmp = behavioralData.duration.all.incentive;
                    RT_perPair_tmp = behavioralData.mod.all.RT_per_pair;
            end
            
            %% extract relevant variables
            trialN_tmp      = behavioralData.mod.all.trialN;
            % behavior data
            onset_inc_tmp   = behavioralData.onset.all.incentive;
            onset_dispE_tmp = behavioralData.onset.all.effortScale;
            onset_Eperf_tmp = behavioralData.onset.all.firstPress;
            lum_inc_tmp     = behavioralData.mod.all.lum_inc;
            %             absInc_tmp         = behavioralData.mod.all.absIncentive;
            inc_tmp         = behavioralData.mod.all.incentive;
            %             trialVal_tmp    = behavioralData.mod.all.trialValence;
            RT_fp_tmp       = behavioralData.mod.all.RT_fp;
            %             perf_tmp        = behavioralData.mod.all.perf;
            % eye data
            pupilD_tmp      = eye_loadStruct.z_pupil_diam;
            eyeTime_tmp     = eye_loadStruct.eyeTime;
            
            %% loop through trials
            for iTrial = 1:length(trialN_tmp)
                trial_idx = trialN_tmp(iTrial);
                jTrial = trial_idx + nTrials_per_run*(jRun - 1);
                
                %% store explicative variables across runs
                if ismember(iRun,[2,3])
                    run1_cstt_sub(jTrial)       = 1;
                    run2_cstt_sub(jTrial)       = 0;
                elseif ismember(iRun,[5,6])
                    run1_cstt_sub(jTrial)       = 0;
                    run2_cstt_sub(jTrial)       = 1;
                end
                lum_inc_sub(jTrial)         = lum_inc_tmp(iTrial);
                trialN_sub(jTrial)          = trial_idx;
                onset_inc_sub(jTrial)       = onset_inc_tmp(iTrial);
                onset_dispE_sub(jTrial)     = onset_dispE_tmp(iTrial);
                onset_Eperf_sub(jTrial)     = onset_Eperf_tmp(iTrial);
                RT_sub(jTrial) = RT_fp_tmp(iTrial);
                X_sub(jTrial) = X_pred_tmp(iTrial);
                switch task_id
                    case 'G'
                        % extract force normalized by Fmax (so force will
                        % be actual force exerted, not force displayed on
                        % screen!)
                        gripForce.dispE(jTrial, :, iS) = gripB_tmp{1,trial_idx}(1:n_grip_samples_per_trial)./Fmax_sub;
                        gripForceDisplay.dispE(jTrial, :, iS) = gripDisplay_tmp{1,trial_idx}(1:n_grip_samples_per_trial);
                        % extract timing and convert from s to ms
                        gripForceTime.inc(jTrial, :, iS) = (grip_time_tmp{1,trial_idx}(1:n_grip_samples_per_trial) - grip_onsets_time_tmp.displayIncentive(trial_idx)).*1000;
                        gripForceTime.dispE(jTrial, :, iS) = (grip_time_tmp{1,trial_idx}(1:n_grip_samples_per_trial) - grip_onsets_time_tmp.displayEffortScale(trial_idx)).*1000;
                        gripForceTime.Eperf(jTrial, :, iS) = (grip_time_tmp{1,trial_idx}(1:n_grip_samples_per_trial)- T0 - onset_Eperf_tmp(iTrial)).*1000;
                    case 'S'
                        inc_dur(jTrial, iS) = inc_dur_tmp(iTrial);
                        E_perf_onset(jTrial,iS) = onset_Eperf_sub(jTrial) - onset_dispE_sub(jTrial);
                        % extract timings (+ convert time from s to ms (*1000))
                        for iPair = 1:nMaxStroop
                            stroopRT(jTrial,iPair,iS) = RT_perPair_tmp(iPair,iTrial);
                            mTime_perPair_perSub_perTrial.inc(jTrial, iPair, iS) = (sum(RT_perPair_tmp(1:iPair,iTrial)) + inc_dur_tmp(iTrial)).*1000;
                            mTime_perPair_perSub_perTrial.dispE(jTrial, iPair, iS) = (sum(RT_perPair_tmp(1:iPair,iTrial))).*1000;
                            mTime_perPair_perSub_perTrial.Eperf(jTrial, iPair, iS) = (sum(RT_perPair_tmp(1:iPair,iTrial)) - E_perf_onset(jTrial,iS)).*1000;
                        end
                end
                
                %% onset incentive extraction
                onset_inc_time = floor(onset_inc_tmp(iTrial)*1000); % ms
                start_inc_extraction    = onset_inc_time - time_bf_onset; % ms
                %             end_extraction      = onset_time + n_effort_time + time_af_RT; % ms
                
                for iTime = 1:nTotalTimeChecking
                    inc_checkTime = floor((start_inc_extraction + iTime)/1000); % ms to seconds
                    %                     eyeXYpupil.inc.Xcoord(jTrial,iTime,iSub) = X_eyeCoord(eyeTime == checkTime);
                    %                     eyeXYpupil.inc.Ycoord(jTrial,iTime,iSub) = Y_eyeCoord(eyeTime == checkTime);
                    eyeXYpupil.inc.pupilD(jTrial,iTime,iS) = pupilD_tmp(eyeTime_tmp == inc_checkTime);
                end
                
                %% onset effort scale extraction
                onset_dispE_time = floor(onset_dispE_tmp(iTrial)*1000); % ms
                start_dispE_extraction    = onset_dispE_time - time_bf_onset; % ms
                for iTime = 1:nTotalTimeChecking
                    dispE_checkTime = floor((start_dispE_extraction + iTime)/1000); % ms to seconds
                    %                     eyeXYpupil.inc.Xcoord(jTrial,iTime,iSub) = X_eyeCoord(eyeTime == checkTime);
                    %                     eyeXYpupil.inc.Ycoord(jTrial,iTime,iSub) = Y_eyeCoord(eyeTime == checkTime);
                    eyeXYpupil.dispE.pupilD(jTrial,iTime,iS) = pupilD_tmp(eyeTime_tmp == dispE_checkTime);
                end
                
                %% onset effort performance extraction
                onset_Eperf_time = floor(onset_Eperf_tmp(iTrial)*1000); % ms
                start_Eperf_extraction    = onset_Eperf_time - time_bf_onset; % ms
                for iTime = 1:nTotalTimeChecking
                    Eperf_checkTime = floor((start_Eperf_extraction + iTime)/1000); % ms to seconds
                    %                     eyeXYpupil.inc.Xcoord(jTrial,iTime,iSub) = X_eyeCoord(eyeTime == checkTime);
                    %                     eyeXYpupil.inc.Ycoord(jTrial,iTime,iSub) = Y_eyeCoord(eyeTime == checkTime);
                    eyeXYpupil.Eperf.pupilD(jTrial,iTime,iS) = pupilD_tmp(eyeTime_tmp == Eperf_checkTime);
                end
                
                %% correct for baseline
                if corr_base == true
                    % incentive
                    baseline_inc_idx_start  = find(eyeTime_tmp == floor(start_inc_extraction/1000));
                    baseline_inc_idx_end    = find(eyeTime_tmp == floor(onset_inc_time/1000));
                    baseline_inc_pupil_tmp = mean( pupilD_tmp(baseline_inc_idx_start:baseline_inc_idx_end) ,'omitnan'); % extract baseline for this trial
                    eyeXYpupil.inc.pupilD(jTrial,:,iS) = eyeXYpupil.inc.pupilD(jTrial,:,iS) - baseline_inc_pupil_tmp; % correct the data for the baseline
                    
                    % effort scale
                    baseline_dispE_idx_start  = find(eyeTime_tmp == floor(start_dispE_extraction/1000));
                    baseline_dispE_idx_end    = find(eyeTime_tmp == floor(onset_dispE_time/1000));
                    baseline_dispE_pupil_tmp = mean( pupilD_tmp(baseline_dispE_idx_start:baseline_dispE_idx_end) ,'omitnan'); % extract baseline for this trial
                    eyeXYpupil.dispE.pupilD(jTrial,:,iS) = eyeXYpupil.dispE.pupilD(jTrial,:,iS) - baseline_dispE_pupil_tmp; % correct the data for the baseline
                    
                    % effort performance
                    baseline_Eperf_idx_start  = find(eyeTime_tmp == floor(start_Eperf_extraction/1000));
                    baseline_Eperf_idx_end    = find(eyeTime_tmp == floor(onset_Eperf_time/1000));
                    baseline_Eperf_pupil_tmp = mean( pupilD_tmp(baseline_Eperf_idx_start:baseline_Eperf_idx_end) ,'omitnan'); % extract baseline for this trial
                    eyeXYpupil.Eperf.pupilD(jTrial,:,iS) = eyeXYpupil.Eperf.pupilD(jTrial,:,iS) - baseline_Eperf_pupil_tmp; % correct the data for the baseline
                end
                
            end % trial loop
            
        end % run loop
        
        %% stats
        if raw_or_z_regs == true
            % zscore all (but binary) variables
            lum_inc_sub = nanzscore(lum_inc_sub);
            RT_sub      = nanzscore(RT_sub);
            X_sub       = nanzscore(X_sub);
        end
        % exclude NaN trials from the correlation
        ok_trials = ~isnan(trialN_sub);
        % correlation
        x0_regs = [run1_cstt_sub,...
            run2_cstt_sub,...
            lum_inc_sub,...
            RT_sub,...
            X_sub];
        x_regs = x0_regs(ok_trials,:);
%         x_regs = mtrx_orthog(x0_regs(ok_trials,:)); % orthogonalize all vars before modelling
        [betas_inc_tmp, betas_dispE_tmp, betas_Eperf_tmp] = deal( NaN(nReg,nTotalTimeChecking));
        for iTime = 1:nTotalTimeChecking
            % incentive period GLM
            betas_inc_tmp(:,iTime) = glmfit(x_regs,...
                eyeXYpupil.inc.pupilD(ok_trials,iTime,iS),...
                'normal','constant','off'); % no global constant (taken in run constant)
            % effort scale period
            betas_dispE_tmp(:,iTime) = glmfit(x_regs,...
                eyeXYpupil.dispE.pupilD(ok_trials,iTime,iS),...
                'normal','constant','off'); % no global constant (taken in run constant)
            % effort performance period
            betas_Eperf_tmp(:,iTime) = glmfit(x_regs,...
                eyeXYpupil.Eperf.pupilD(ok_trials,iTime,iS),...
                'normal','constant','off'); % no global constant (taken in run constant)
        end
        % store betas
        for iBeta = 1:nReg
            reg_nm = regs_names{iBeta};
            betas.inc.(reg_nm)(:,iS) = betas_inc_tmp(iBeta, :);
            betas.dispE.(reg_nm)(:,iS) = betas_dispE_tmp(iBeta, :);
            betas.Eperf.(reg_nm)(:,iS) = betas_Eperf_tmp(iBeta, :);
        end % beta loop
        
        %% median split (after removing luminance effect)
        
        % remove luminance effect from pupil
        for iTime = 1:nTotalTimeChecking
            eyeXYpupil.inc.pupilD_noLum(:,iTime,iS)     = eyeXYpupil.inc.pupilD(:,iTime,iS)     - betas.inc.Lum(iTime,iS).*lum_inc_sub;
            eyeXYpupil.dispE.pupilD_noLum(:,iTime,iS)   = eyeXYpupil.dispE.pupilD(:,iTime,iS)   - betas.dispE.Lum(iTime,iS).*lum_inc_sub;
            eyeXYpupil.Eperf.pupilD_noLum(:,iTime,iS)   = eyeXYpupil.Eperf.pupilD(:,iTime,iS)   - betas.Eperf.Lum(iTime,iS).*lum_inc_sub;
        end
        
        % median split
        for iReg = 3:nReg % avoid run constants
            reg_nm = regs_names{iReg};
            switch reg_nm
                case 'Lum'
                    reg_vals = lum_inc_sub;
                case 'RT'
                    reg_vals = RT_sub;
                case 'E'
                    reg_vals = X_sub;
%                 case 'Inc'
%                     reg_vals = inc_sub;
            end
            med = median(reg_vals, 'omitnan');
            low_trials = reg_vals <= med;
            high_trials = reg_vals >= med;
            
            
            trialN_mSplit_low = 1:sum(low_trials ~= 0);
            trialN_mSplit_high = 1:sum(high_trials ~= 0); % different number btw low and high trials could happen for run split
            for iEye_type = 1:n_eyeDataType
                eye_type_nm = eye_data_type{iEye_type};
                if sum(low_trials) > 0 % for run median split especially
                    eyeXYpupil_mSplit.inc.(reg_nm).low.(eye_type_nm)(trialN_mSplit_low,:,iS)    = eyeXYpupil.inc.(eye_type_nm)(low_trials,:,iS);
                    eyeXYpupil_mSplit.dispE.(reg_nm).low.(eye_type_nm)(trialN_mSplit_low,:,iS)  = eyeXYpupil.dispE.(eye_type_nm)(low_trials,:,iS);
                    eyeXYpupil_mSplit.Eperf.(reg_nm).low.(eye_type_nm)(trialN_mSplit_low,:,iS)  = eyeXYpupil.Eperf.(eye_type_nm)(low_trials,:,iS);
                end
                if sum(high_trials) > 0 % for run median split especially
                    eyeXYpupil_mSplit.inc.(reg_nm).high.(eye_type_nm)(trialN_mSplit_high,:,iS)   = eyeXYpupil.inc.(eye_type_nm)(high_trials,:,iS);
                    eyeXYpupil_mSplit.dispE.(reg_nm).high.(eye_type_nm)(trialN_mSplit_high,:,iS) = eyeXYpupil.dispE.(eye_type_nm)(high_trials,:,iS);
                    eyeXYpupil_mSplit.Eperf.(reg_nm).high.(eye_type_nm)(trialN_mSplit_high,:,iS)   = eyeXYpupil.Eperf.(eye_type_nm)(high_trials,:,iS);
                end
            end
        end % regressor loop
        
    end % run filter
    
    %% average per subject across trials
    for iReg = 1:nReg
        reg_nm = regs_names{iReg};
        for iEye_type = 1:n_eyeDataType
            eye_type_nm = eye_data_type{iEye_type};
            eyeXYpupil_mSplit.inc.mean_per_sub.(reg_nm).low.(eye_type_nm)(:,iS)       = mean(eyeXYpupil_mSplit.inc.(reg_nm).low.(eye_type_nm)(:,:,iS), 1,'omitnan');
            eyeXYpupil_mSplit.inc.mean_per_sub.(reg_nm).high.(eye_type_nm)(:,iS)      = mean(eyeXYpupil_mSplit.inc.(reg_nm).high.(eye_type_nm)(:,:,iS), 1,'omitnan');
            eyeXYpupil_mSplit.dispE.mean_per_sub.(reg_nm).low.(eye_type_nm)(:,iS)     = mean(eyeXYpupil_mSplit.dispE.(reg_nm).low.(eye_type_nm)(:,:,iS), 1,'omitnan');
            eyeXYpupil_mSplit.dispE.mean_per_sub.(reg_nm).high.(eye_type_nm)(:,iS)    = mean(eyeXYpupil_mSplit.dispE.(reg_nm).high.(eye_type_nm)(:,:,iS), 1,'omitnan');
            eyeXYpupil_mSplit.Eperf.mean_per_sub.(reg_nm).low.(eye_type_nm)(:,iS)     = mean(eyeXYpupil_mSplit.Eperf.(reg_nm).low.(eye_type_nm)(:,:,iS), 1,'omitnan');
            eyeXYpupil_mSplit.Eperf.mean_per_sub.(reg_nm).high.(eye_type_nm)(:,iS)    = mean(eyeXYpupil_mSplit.Eperf.(reg_nm).high.(eye_type_nm)(:,:,iS), 1,'omitnan');
        end
    end % regressor loop
    
    % same for force/RT
    switch task_id
        case 'G'
            gripForce_perSub.dispE(:,iS) = mean(gripForce.dispE(:,:,iS),1,'omitnan');
            gripForceDisplay_perSub.dispE(:,iS) = mean(gripForceDisplay.dispE(:,:,iS),1,'omitnan');
            gripForceTime_perSub.inc(:,iS) = mean(gripForceTime.inc(:,:,iS),1,'omitnan');
            gripForceTime_perSub.dispE(:,iS) = mean(gripForceTime.dispE(:,:,iS),1,'omitnan');
            gripForceTime_perSub.Eperf(:,iS) = mean(gripForceTime.Eperf(:,:,iS),1,'omitnan');
        case 'S'
            for iPair = 1:nMaxStroop
                stroopRT_perSub(iPair,iS) = mean(stroopRT(:,iPair,iS),1,'omitnan');
                
                 % convert all timings from s to ms (*1000)
                 % extract timings per phase of the trial
                mTime_perPair_perSub.inc(iPair,iS) = mean(mTime_perPair_perSub_perTrial.inc(:,iPair,iS),1,'omitnan');
                mTime_perPair_perSub.dispE(iPair,iS) = mean(mTime_perPair_perSub_perTrial.dispE(:,iPair,iS),1,'omitnan');
                mTime_perPair_perSub.Eperf(iPair,iS) = mean(mTime_perPair_perSub_perTrial.Eperf(:,iPair,iS),1,'omitnan');
            end
    end
    
    %% signal subject done
    disp(['subject ',num2str(iS),'/',num2str(NS),' done']);
end % subject loop

%% slight smooth of the data to improve RFT application and graph visualization
smooth_kernel = 200; % short smoothing
for iReg = 1:nReg
    reg_nm = regs_names{iReg};
    
    for iS = 1:NS
        % betas
        betas.inc.(reg_nm)(:,iS) = smooth_NC(betas.inc.(reg_nm)(:,iS)', smooth_kernel);
        betas.dispE.(reg_nm)(:,iS) = smooth_NC(betas.dispE.(reg_nm)(:,iS)', smooth_kernel);
        betas.Eperf.(reg_nm)(:,iS) = smooth_NC(betas.Eperf.(reg_nm)(:,iS)', smooth_kernel);
        
        % eye data
        for iEye_type = 1:n_eyeDataType
            eye_type_nm = eye_data_type{iEye_type};
            eyeXYpupil_mSplit.inc.mean_per_sub.(reg_nm).low.(eye_type_nm)(:,iS)     = smooth_NC(eyeXYpupil_mSplit.inc.mean_per_sub.(reg_nm).low.(eye_type_nm)(:,iS)', smooth_kernel);
            eyeXYpupil_mSplit.inc.mean_per_sub.(reg_nm).high.(eye_type_nm)(:,iS)    = smooth_NC(eyeXYpupil_mSplit.inc.mean_per_sub.(reg_nm).high.(eye_type_nm)(:,iS)', smooth_kernel);
            eyeXYpupil_mSplit.dispE.mean_per_sub.(reg_nm).low.(eye_type_nm)(:,iS)   = smooth_NC(eyeXYpupil_mSplit.dispE.mean_per_sub.(reg_nm).low.(eye_type_nm)(:,iS)', smooth_kernel);
            eyeXYpupil_mSplit.dispE.mean_per_sub.(reg_nm).high.(eye_type_nm)(:,iS)  = smooth_NC(eyeXYpupil_mSplit.dispE.mean_per_sub.(reg_nm).high.(eye_type_nm)(:,iS)', smooth_kernel);
            eyeXYpupil_mSplit.Eperf.mean_per_sub.(reg_nm).low.(eye_type_nm)(:,iS)   = smooth_NC(eyeXYpupil_mSplit.Eperf.mean_per_sub.(reg_nm).low.(eye_type_nm)(:,iS)', smooth_kernel);
            eyeXYpupil_mSplit.Eperf.mean_per_sub.(reg_nm).high.(eye_type_nm)(:,iS)  = smooth_NC(eyeXYpupil_mSplit.Eperf.mean_per_sub.(reg_nm).high.(eye_type_nm)(:,iS)', smooth_kernel);
        end
    end
end

%% average + SEM across subs
for iReg = 1:nReg
    reg_nm = regs_names{iReg};
    
    % betas
    % mean
    betas.inc.mean_aSubs.(reg_nm)   = mean(betas.inc.(reg_nm), 2, 'omitnan');
    betas.dispE.mean_aSubs.(reg_nm) = mean(betas.dispE.(reg_nm), 2, 'omitnan');
    betas.Eperf.mean_aSubs.(reg_nm) = mean(betas.Eperf.(reg_nm), 2, 'omitnan');
    % SEM
    betas.inc.sem_aSubs.(reg_nm)   = sem(betas.inc.(reg_nm), 2);
    betas.dispE.sem_aSubs.(reg_nm) = sem(betas.dispE.(reg_nm), 2);
    betas.Eperf.sem_aSubs.(reg_nm) = sem(betas.Eperf.(reg_nm), 2);
    
    % eye data
    for iEye_type = 1:n_eyeDataType
        eye_type_nm = eye_data_type{iEye_type};
        % mean
        eyeXYpupil_mSplit.inc.mean_aSubs.(reg_nm).low.(eye_type_nm)     = mean(eyeXYpupil_mSplit.inc.mean_per_sub.(reg_nm).low.(eye_type_nm), 2, 'omitnan');
        eyeXYpupil_mSplit.inc.mean_aSubs.(reg_nm).high.(eye_type_nm)    = mean(eyeXYpupil_mSplit.inc.mean_per_sub.(reg_nm).high.(eye_type_nm), 2, 'omitnan');
        eyeXYpupil_mSplit.dispE.mean_aSubs.(reg_nm).low.(eye_type_nm)   = mean(eyeXYpupil_mSplit.dispE.mean_per_sub.(reg_nm).low.(eye_type_nm), 2, 'omitnan');
        eyeXYpupil_mSplit.dispE.mean_aSubs.(reg_nm).high.(eye_type_nm)  = mean(eyeXYpupil_mSplit.dispE.mean_per_sub.(reg_nm).high.(eye_type_nm), 2, 'omitnan');
        eyeXYpupil_mSplit.Eperf.mean_aSubs.(reg_nm).low.(eye_type_nm)     = mean(eyeXYpupil_mSplit.Eperf.mean_per_sub.(reg_nm).low.(eye_type_nm), 2, 'omitnan');
        eyeXYpupil_mSplit.Eperf.mean_aSubs.(reg_nm).high.(eye_type_nm)    = mean(eyeXYpupil_mSplit.Eperf.mean_per_sub.(reg_nm).high.(eye_type_nm), 2, 'omitnan');
        % sem
        eyeXYpupil_mSplit.inc.sem_aSubs.(reg_nm).low.(eye_type_nm)     = sem(eyeXYpupil_mSplit.inc.mean_per_sub.(reg_nm).low.(eye_type_nm), 2);
        eyeXYpupil_mSplit.inc.sem_aSubs.(reg_nm).high.(eye_type_nm)    = sem(eyeXYpupil_mSplit.inc.mean_per_sub.(reg_nm).high.(eye_type_nm), 2);
        eyeXYpupil_mSplit.dispE.sem_aSubs.(reg_nm).low.(eye_type_nm)   = sem(eyeXYpupil_mSplit.dispE.mean_per_sub.(reg_nm).low.(eye_type_nm), 2);
        eyeXYpupil_mSplit.dispE.sem_aSubs.(reg_nm).high.(eye_type_nm)  = sem(eyeXYpupil_mSplit.dispE.mean_per_sub.(reg_nm).high.(eye_type_nm), 2);
        eyeXYpupil_mSplit.Eperf.sem_aSubs.(reg_nm).low.(eye_type_nm)   = sem(eyeXYpupil_mSplit.Eperf.mean_per_sub.(reg_nm).low.(eye_type_nm), 2);
        eyeXYpupil_mSplit.Eperf.sem_aSubs.(reg_nm).high.(eye_type_nm)  = sem(eyeXYpupil_mSplit.Eperf.mean_per_sub.(reg_nm).high.(eye_type_nm), 2);
    end % luminance-corrected or not
end % regressor loop
% average force/RT across subjects
switch task_id
    case 'G'
        [gripForce_m_aSubs.dispE,...
            gripForce_sem_aSubs.dispE]  = mean_sem_sd(gripForce_perSub.dispE.*100,2);
        [gripForceDisplay_m_aSubs.dispE,...
            gripForceDisplay_sem_aSubs.dispE]  = mean_sem_sd(gripForceDisplay_perSub.dispE,2);
        gripForceTime_m_aSubs.inc = mean(gripForceTime_perSub.inc,2,'omitnan');
        gripForceTime_m_aSubs.dispE = mean(gripForceTime_perSub.dispE,2,'omitnan');
        gripForceTime_m_aSubs.Eperf = mean(gripForceTime_perSub.Eperf,2,'omitnan');
    case 'S'
        [stroopRT_m_aSubs,...
            stroopRT_sem_aSubs] = mean_sem_sd(stroopRT_perSub,2);
        [mTime_perPair_m_aSubs.inc,...
            mTime_perPair_sem_aSubs.inc]  = mean_sem_sd(mTime_perPair_perSub.inc,2);
        [mTime_perPair_m_aSubs.dispE,...
            mTime_perPair_sem_aSubs.dispE]  = mean_sem_sd(mTime_perPair_perSub.dispE,2);
        [mTime_perPair_m_aSubs.Eperf,...
            mTime_perPair_sem_aSubs.Eperf]  = mean_sem_sd(mTime_perPair_perSub.Eperf,2);
end

%% see significant parts
RFT_prm.Xvector = ones(NS,1);
RFT_prm.verbose = 0;
RFT_prm.con = 1;
RFT_prm.conType = 'F';
RFT_prm.pValCorr_threshold  = 0.05;
for iReg = 3:nReg
    reg_nm = regs_names{iReg};
    
    RFT_prm.subfield_nm = reg_nm;
    
    %% perform RFT
    % beta inc
    data_to_test_inc.(reg_nm) = NaN(NS, nTotalTimeChecking);
    data_to_test_inc.(reg_nm)(1:NS,1:nTotalTimeChecking) = betas.inc.(reg_nm)';
    [ pVal_X_inc_tmp, ~ ] = RFT_extraction( data_to_test_inc, RFT_prm );
    pval_X_signif_clusters_inc.(reg_nm) = pVal_X_inc_tmp.(reg_nm);
    % extract corresponding Y values for display on the mean
    pval_Y_signif_clusters_inc.(reg_nm) = betas.inc.mean_aSubs.(reg_nm)( pval_X_signif_clusters_inc.(reg_nm) );
    
    % beta dispE
    data_to_test_dispE.(reg_nm) = NaN(NS, nTotalTimeChecking);
    data_to_test_dispE.(reg_nm)(1:NS,1:nTotalTimeChecking) = betas.dispE.(reg_nm)';
    [ pVal_X_dispE_tmp, ~ ] = RFT_extraction( data_to_test_dispE, RFT_prm );
    pval_X_signif_clusters_dispE.(reg_nm) = pVal_X_dispE_tmp.(reg_nm);
    % extract corresponding Y values for display on the mean
    pval_Y_signif_clusters_dispE.(reg_nm) = betas.dispE.mean_aSubs.(reg_nm)( pval_X_signif_clusters_dispE.(reg_nm) );
    
    % beta Eperf
    data_to_test_Eperf.(reg_nm) = NaN(NS, nTotalTimeChecking);
    data_to_test_Eperf.(reg_nm)(1:NS,1:nTotalTimeChecking) = betas.Eperf.(reg_nm)';
    [ pVal_X_Eperf_tmp, ~ ] = RFT_extraction( data_to_test_Eperf, RFT_prm );
    pval_X_signif_clusters_Eperf.(reg_nm) = pVal_X_Eperf_tmp.(reg_nm);
    % extract corresponding Y values for display on the mean
    pval_Y_signif_clusters_Eperf.(reg_nm) = betas.Eperf.mean_aSubs.(reg_nm)( pval_X_signif_clusters_Eperf.(reg_nm) );
    
    %% median split
    % median split inc
    [data_to_test_inc_mSplit_pupilD.(reg_nm),...
        data_to_test_inc_mSplit_pupilD_noLum.(reg_nm)] = deal( NaN(NS, nTotalTimeChecking) );
    for iS = 1:NS
        for iTime = 1:nTotalTimeChecking
            data_to_test_inc_mSplit_pupilD.(reg_nm)(iS, iTime)          = eyeXYpupil_mSplit.inc.mean_per_sub.(reg_nm).low.pupilD(iTime,iS)          - eyeXYpupil_mSplit.inc.mean_per_sub.(reg_nm).high.pupilD(iTime,iS);
            data_to_test_inc_mSplit_pupilD_noLum.(reg_nm)(iS, iTime)    = eyeXYpupil_mSplit.inc.mean_per_sub.(reg_nm).low.pupilD_noLum(iTime,iS)    - eyeXYpupil_mSplit.inc.mean_per_sub.(reg_nm).high.pupilD_noLum(iTime,iS);
        end
    end
    % exclude NaN subjects (for run constant especially)
    ok_subs_mSplit_inc     = ~isnan(data_to_test_inc_mSplit_pupilD.(reg_nm)(:, 1));
    RFT_prm_mSplit_inc = RFT_prm;
    if sum(ok_subs_mSplit_inc) < NS
        RFT_prm_mSplit_inc.Xvector = RFT_prm_mSplit_inc.Xvector(ok_subs_mSplit_inc);
        data_to_test_inc_mSplit_pupilD.(reg_nm)         = data_to_test_inc_mSplit_pupilD.(reg_nm)(ok_subs_mSplit_inc,:);
        data_to_test_inc_mSplit_pupilD_noLum.(reg_nm)   = data_to_test_inc_mSplit_pupilD_noLum.(reg_nm)(ok_subs_mSplit_inc,:);
    end
    [ pVal_X_inc_mSplit_pupilD_tmp, ~ ]         = RFT_extraction( data_to_test_inc_mSplit_pupilD, RFT_prm_mSplit_inc );
    [ pVal_X_inc_mSplit_pupilD_noLum_tmp, ~ ]   = RFT_extraction( data_to_test_inc_mSplit_pupilD_noLum, RFT_prm_mSplit_inc );
    pval_X_signif_inc_mSplit_pupilD.(reg_nm)         = pVal_X_inc_mSplit_pupilD_tmp.(reg_nm);
    pval_X_signif_inc_mSplit_pupilD_noLum.(reg_nm)   = pVal_X_inc_mSplit_pupilD_noLum_tmp.(reg_nm);
    % extract corresponding Y values for display on the mean
    pval_Y_signif_inc_mSplit_pupilD.(reg_nm).low        = eyeXYpupil_mSplit.inc.mean_aSubs.(reg_nm).low.pupilD( pval_X_signif_inc_mSplit_pupilD.(reg_nm) );
    pval_Y_signif_inc_mSplit_pupilD.(reg_nm).high       = eyeXYpupil_mSplit.inc.mean_aSubs.(reg_nm).high.pupilD( pval_X_signif_inc_mSplit_pupilD.(reg_nm) );
    pval_Y_signif_inc_mSplit_pupilD_noLum.(reg_nm).low  = eyeXYpupil_mSplit.inc.mean_aSubs.(reg_nm).low.pupilD_noLum( pval_X_signif_inc_mSplit_pupilD_noLum.(reg_nm) );
    pval_Y_signif_inc_mSplit_pupilD_noLum.(reg_nm).high = eyeXYpupil_mSplit.inc.mean_aSubs.(reg_nm).high.pupilD_noLum( pval_X_signif_inc_mSplit_pupilD_noLum.(reg_nm) );
    
    % median split dispE
    [data_to_test_dispE_mSplit_pupilD.(reg_nm),...
        data_to_test_dispE_mSplit_pupilD_noLum.(reg_nm)] = deal( NaN(NS, nTotalTimeChecking) );
    for iS = 1:NS
        for iTime = 1:nTotalTimeChecking
            data_to_test_dispE_mSplit_pupilD.(reg_nm)(iS, iTime)        = eyeXYpupil_mSplit.dispE.mean_per_sub.(reg_nm).low.pupilD(iTime,iS)        - eyeXYpupil_mSplit.dispE.mean_per_sub.(reg_nm).high.pupilD(iTime,iS);
            data_to_test_dispE_mSplit_pupilD_noLum.(reg_nm)(iS, iTime)  = eyeXYpupil_mSplit.dispE.mean_per_sub.(reg_nm).low.pupilD_noLum(iTime,iS)  - eyeXYpupil_mSplit.dispE.mean_per_sub.(reg_nm).high.pupilD_noLum(iTime,iS);
        end
    end
    % exclude NaN subjects (for run constant especially)
    ok_subs_mSplit_dispE     = ~isnan(data_to_test_dispE_mSplit_pupilD.(reg_nm)(:, 1));
    RFT_prm_mSplit_dispE = RFT_prm;
    if sum(ok_subs_mSplit_dispE) < NS
        RFT_prm_mSplit_dispE.Xvector = RFT_prm_mSplit_dispE.Xvector(ok_subs_mSplit_dispE);
        data_to_test_dispE_mSplit_pupilD.(reg_nm)         = data_to_test_dispE_mSplit_pupilD.(reg_nm)(ok_subs_mSplit_dispE,:);
        data_to_test_dispE_mSplit_pupilD_noLum.(reg_nm)   = data_to_test_dispE_mSplit_pupilD_noLum.(reg_nm)(ok_subs_mSplit_dispE,:);
    end
    [ pVal_X_dispE_mSplit_pupilD_tmp, ~ ]         = RFT_extraction( data_to_test_dispE_mSplit_pupilD, RFT_prm_mSplit_dispE );
    [ pVal_X_dispE_mSplit_pupilD_noLum_tmp, ~ ]   = RFT_extraction( data_to_test_dispE_mSplit_pupilD_noLum, RFT_prm_mSplit_dispE );
    pval_X_signif_dispE_mSplit_pupilD.(reg_nm)         = pVal_X_dispE_mSplit_pupilD_tmp.(reg_nm);
    pval_X_signif_dispE_mSplit_pupilD_noLum.(reg_nm)   = pVal_X_dispE_mSplit_pupilD_noLum_tmp.(reg_nm);
    % extract corresponding Y values for display on the mean
    pval_Y_signif_dispE_mSplit_pupilD.(reg_nm).low = eyeXYpupil_mSplit.dispE.mean_aSubs.(reg_nm).low.pupilD( pval_X_signif_dispE_mSplit_pupilD.(reg_nm) );
    pval_Y_signif_dispE_mSplit_pupilD.(reg_nm).high = eyeXYpupil_mSplit.dispE.mean_aSubs.(reg_nm).high.pupilD(  pval_X_signif_dispE_mSplit_pupilD.(reg_nm) );
    pval_Y_signif_dispE_mSplit_pupilD_noLum.(reg_nm).low = eyeXYpupil_mSplit.dispE.mean_aSubs.(reg_nm).low.pupilD_noLum( pval_X_signif_dispE_mSplit_pupilD_noLum.(reg_nm) );
    pval_Y_signif_dispE_mSplit_pupilD_noLum.(reg_nm).high = eyeXYpupil_mSplit.dispE.mean_aSubs.(reg_nm).high.pupilD_noLum( pval_X_signif_dispE_mSplit_pupilD_noLum.(reg_nm) );
    
    % median split Effort performance
    [data_to_test_Eperf_mSplit_pupilD.(reg_nm),...
        data_to_test_Eperf_mSplit_pupilD_noLum.(reg_nm)] = deal( NaN(NS, nTotalTimeChecking) );
    for iS = 1:NS
        for iTime = 1:nTotalTimeChecking
            data_to_test_Eperf_mSplit_pupilD.(reg_nm)(iS, iTime)          = eyeXYpupil_mSplit.Eperf.mean_per_sub.(reg_nm).low.pupilD(iTime,iS)          - eyeXYpupil_mSplit.Eperf.mean_per_sub.(reg_nm).high.pupilD(iTime,iS);
            data_to_test_Eperf_mSplit_pupilD_noLum.(reg_nm)(iS, iTime)    = eyeXYpupil_mSplit.Eperf.mean_per_sub.(reg_nm).low.pupilD_noLum(iTime,iS)    - eyeXYpupil_mSplit.Eperf.mean_per_sub.(reg_nm).high.pupilD_noLum(iTime,iS);
        end
    end
    % exclude NaN subjects (for run constant especially)
    ok_subs_mSplit_Eperf     = ~isnan(data_to_test_Eperf_mSplit_pupilD.(reg_nm)(:, 1));
    RFT_prm_mSplit_Eperf = RFT_prm;
    if sum(ok_subs_mSplit_Eperf) < NS
        RFT_prm_mSplit_Eperf.Xvector = RFT_prm_mSplit_Eperf.Xvector(ok_subs_mSplit_Eperf);
        data_to_test_Eperf_mSplit_pupilD.(reg_nm)         = data_to_test_Eperf_mSplit_pupilD.(reg_nm)(ok_subs_mSplit_Eperf,:);
        data_to_test_Eperf_mSplit_pupilD_noLum.(reg_nm)   = data_to_test_Eperf_mSplit_pupilD_noLum.(reg_nm)(ok_subs_mSplit_Eperf,:);
    end
    [ pVal_X_Eperf_mSplit_pupilD_tmp, ~ ]         = RFT_extraction( data_to_test_Eperf_mSplit_pupilD, RFT_prm_mSplit_Eperf );
    [ pVal_X_Eperf_mSplit_pupilD_noLum_tmp, ~ ]   = RFT_extraction( data_to_test_Eperf_mSplit_pupilD_noLum, RFT_prm_mSplit_Eperf );
    pval_X_signif_Eperf_mSplit_pupilD.(reg_nm)         = pVal_X_Eperf_mSplit_pupilD_tmp.(reg_nm);
    pval_X_signif_Eperf_mSplit_pupilD_noLum.(reg_nm)   = pVal_X_Eperf_mSplit_pupilD_noLum_tmp.(reg_nm);
    % extract corresponding Y values for display on the mean
    pval_Y_signif_Eperf_mSplit_pupilD.(reg_nm).low        = eyeXYpupil_mSplit.Eperf.mean_aSubs.(reg_nm).low.pupilD( pval_X_signif_Eperf_mSplit_pupilD.(reg_nm) );
    pval_Y_signif_Eperf_mSplit_pupilD.(reg_nm).high       = eyeXYpupil_mSplit.Eperf.mean_aSubs.(reg_nm).high.pupilD( pval_X_signif_Eperf_mSplit_pupilD.(reg_nm) );
    pval_Y_signif_Eperf_mSplit_pupilD_noLum.(reg_nm).low  = eyeXYpupil_mSplit.Eperf.mean_aSubs.(reg_nm).low.pupilD_noLum( pval_X_signif_Eperf_mSplit_pupilD_noLum.(reg_nm) );
    pval_Y_signif_Eperf_mSplit_pupilD_noLum.(reg_nm).high = eyeXYpupil_mSplit.Eperf.mean_aSubs.(reg_nm).high.pupilD_noLum( pval_X_signif_Eperf_mSplit_pupilD_noLum.(reg_nm) );
    
end % regressor loop

%% graph
x_inc_eye   = (-time_bf_onset):(nTotalTimeChecking - time_bf_onset - 1);
x_dispE_eye  = (-time_bf_onset):(nTotalTimeChecking - time_bf_onset - 1);
x_Eperf_eye   = (-time_bf_onset):(nTotalTimeChecking - time_bf_onset - 1);
pSize = 25;
colours = {'k','k',...
    'k','m','g','m','m','k','b','c','r'};
lWidth = 3;

%% betas figure
betas_fig = fig();

% incentive period
subplot(1,3,1);
switch task_id
    case 'G'
        colororder({'k','b'});
        yyaxis left;
end
for iReg = 4:nReg
    reg_nm = regs_names{iReg};
    col_nm = colours{iReg};
    curve_inc_hdl.(reg_nm) = jbfill(x_inc_eye,...
        betas.inc.mean_aSubs.(reg_nm)' + betas.inc.sem_aSubs.(reg_nm)',...
        betas.inc.mean_aSubs.(reg_nm)' - betas.inc.sem_aSubs.(reg_nm)',...
        betas.inc.mean_aSubs.(reg_nm)',...
        col_nm);
    hold on;
    
    % mark significant clusters from the RFT analysis
    x_signif = x_inc_eye( pval_X_signif_clusters_inc.(reg_nm) );
    % mark the mean
    plot(x_signif,...
        pval_Y_signif_clusters_inc.(reg_nm),[col_nm,' *']);
    % mark significant clusters on the top of the graph
    ypos_perc = (1/100)*iReg;
    place_signif_curve( betas_fig, x_signif, col_nm, ypos_perc  );
end % regresssor loop
line([0 0],ylim(),'Color','k','LineWidth',lWidth); % mark onset of incentive
line(xlim,[0 0],'Color','k','LineWidth',lWidth); % mark beta at zero
xlim([x_inc_eye(1), x_inc_eye(end)]);
xlabel('Time (ms) (incentive period)');
ylabel('Parameter estimates (a.u.)');
% add force (grip)/RT (stroop) overlap
switch task_id
    case 'G'
        yyaxis right;
        jbfill(gripForceTime_m_aSubs.inc',...
            gripForce_m_aSubs.dispE' + gripForce_sem_aSubs.dispE',...
            gripForce_m_aSubs.dispE' - gripForce_sem_aSubs.dispE',...
            gripForce_m_aSubs.dispE',...
            'b');
        ylabel('Force (%)');
        ylim([0 50]);
    case 'S'
        y_size = ylim();
        for iPair = 1:nMaxStroop
            RT_pair_inc = mTime_perPair_m_aSubs.inc(iPair);
            RT_pair_inc_top = mTime_perPair_m_aSubs.inc(iPair) + mTime_perPair_sem_aSubs.inc(iPair);
            RT_pair_inc_low = mTime_perPair_m_aSubs.inc(iPair) - mTime_perPair_sem_aSubs.inc(iPair);
            jbfill_x([RT_pair_inc, RT_pair_inc],...
                [RT_pair_inc_top RT_pair_inc_top],...
                [RT_pair_inc_low RT_pair_inc_low],...
                y_size,...
                'r ');
        end
end
% legend([curve_inc_hdl.Lum,...
%     curve_inc_hdl.RT,...
%     curve_inc_hdl.X],...
%     regs_names(3:end),...
%     'Location','southeast');
legend([curve_inc_hdl.RT,...
    curve_inc_hdl.X],...
    {'RT','E'},...
    'Location','southeast');
legend('boxoff');
legend_size(pSize);

% effort scale period
subplot(1,3,2);
switch task_id
    case 'G'
        colororder({'k','b'});
        yyaxis left;
end
for iReg = 4:nReg
    reg_nm = regs_names{iReg};
    col_nm = colours{iReg};
    curve_dispE_hdl.(reg_nm) = jbfill(x_dispE_eye,...
        betas.dispE.mean_aSubs.(reg_nm)' + betas.dispE.sem_aSubs.(reg_nm)',...
        betas.dispE.mean_aSubs.(reg_nm)' - betas.dispE.sem_aSubs.(reg_nm)',...
        betas.dispE.mean_aSubs.(reg_nm)',...
        col_nm);
    
    % mark significant clusters from the RFT analysis
    x_signif = x_dispE_eye( pval_X_signif_clusters_dispE.(reg_nm) );
    % mark the mean
    plot(x_signif,...
        pval_Y_signif_clusters_dispE.(reg_nm),[col_nm,' *']);
    % mark significant clusters on the top of the graph
    ypos_perc = (1/100)*iReg;
    place_signif_curve( betas_fig, x_signif, col_nm, ypos_perc  );
end % regresssor loop
line([0 0],ylim(),'Color','k','LineWidth',lWidth); % mark onset of effort scale
line(xlim,[0 0],'Color','k','LineWidth',lWidth); % mark beta at zero
if strcmp(task_id,'G')
    line([5000 5000],ylim(),'Color','k','LineWidth',lWidth); % mark end of effort period
end
xlim([x_dispE_eye(1), x_dispE_eye(end)]);
xlabel('Time (ms) (effort scale period)');
ylabel('Parameter estimates (a.u.)');
% add force (grip)/RT (stroop) overlap
switch task_id
    case 'G'
        yyaxis right;
        jbfill(gripForceTime_m_aSubs.dispE',...
            gripForce_m_aSubs.dispE' + gripForce_sem_aSubs.dispE',...
            gripForce_m_aSubs.dispE' - gripForce_sem_aSubs.dispE',...
            gripForce_m_aSubs.dispE',...
            'b');
        ylabel('Force (%)');
        ylim([0 50]);
    case 'S'
        y_size = ylim();
        for iPair = 1:nMaxStroop
            RT_pair_dispE = mTime_perPair_m_aSubs.dispE(iPair);
            RT_pair_dispE_top = mTime_perPair_m_aSubs.dispE(iPair) + mTime_perPair_sem_aSubs.dispE(iPair);
            RT_pair_dispE_low = mTime_perPair_m_aSubs.dispE(iPair) - mTime_perPair_sem_aSubs.dispE(iPair);
            jbfill_x([RT_pair_dispE, RT_pair_dispE],...
                [RT_pair_dispE_top RT_pair_dispE_top],...
                [RT_pair_dispE_low RT_pair_dispE_low],...
                y_size,...
                'r ');
        end
end
% legend([curve_dispE_hdl.Lum,...
%     curve_dispE_hdl.RT,...
%     curve_dispE_hdl.X],...
%     regs_names(3:end),...
%     'Location','northwest');
legend([curve_dispE_hdl.RT,...
    curve_dispE_hdl.X],...
    {'RT','E'},...
    'Location','northwest');
legend('boxoff');
legend_size(pSize);

% effort performance
subplot(1,3,3);
switch task_id
    case 'G'
        colororder({'k','b'});
        yyaxis left;
end
for iReg = 4:nReg
    reg_nm = regs_names{iReg};
    col_nm = colours{iReg};
    curve_Eperf_hdl.(reg_nm) = jbfill(x_Eperf_eye,...
        betas.Eperf.mean_aSubs.(reg_nm)' + betas.Eperf.sem_aSubs.(reg_nm)',...
        betas.Eperf.mean_aSubs.(reg_nm)' - betas.Eperf.sem_aSubs.(reg_nm)',...
        betas.Eperf.mean_aSubs.(reg_nm)',...
        col_nm);
    hold on;
    
    % mark significant clusters from the RFT analysis
    x_Eperf_signif = x_Eperf_eye( pval_X_signif_clusters_Eperf.(reg_nm) );
    % mark the mean
    plot(x_Eperf_signif,...
        pval_Y_signif_clusters_Eperf.(reg_nm),[col_nm,' *']);
    % mark significant clusters on the top of the graph
    ypos_perc = (1/100)*iReg;
    place_signif_curve( betas_fig, x_Eperf_signif, col_nm, ypos_perc  );
end % regresssor loop
line([0 0],ylim(),'Color','k','LineWidth',lWidth); % mark onset of incentive
line(xlim,[0 0],'Color','k','LineWidth',lWidth); % mark beta at zero
xlim([x_Eperf_eye(1), x_Eperf_eye(end)]);
xlabel('Time after start of effort (ms)');
ylabel('Parameter estimates (a.u.)');
% add force (grip)/RT (stroop) overlap
switch task_id
    case 'G'
        yyaxis right;
        jbfill(gripForceTime_m_aSubs.Eperf',...
            gripForce_m_aSubs.dispE' + gripForce_sem_aSubs.dispE',...
            gripForce_m_aSubs.dispE' - gripForce_sem_aSubs.dispE',...
            gripForce_m_aSubs.dispE',...
            'b');
        ylabel('Force (%)');
        ylim([0 50]);
    case 'S'
        y_size = ylim();
        for iPair = 1:nMaxStroop
            RT_pair_Eperf = mTime_perPair_m_aSubs.Eperf(iPair);
            RT_pair_Eperf_top = mTime_perPair_m_aSubs.Eperf(iPair) + mTime_perPair_sem_aSubs.Eperf(iPair);
            RT_pair_Eperf_low = mTime_perPair_m_aSubs.Eperf(iPair) - mTime_perPair_sem_aSubs.Eperf(iPair);
            jbfill_x([RT_pair_Eperf, RT_pair_Eperf],...
                [RT_pair_Eperf_top RT_pair_Eperf_top],...
                [RT_pair_Eperf_low RT_pair_Eperf_low],...
                y_size,...
                'r ');
        end
end
% legend([curve_Eperf_hdl.Lum,...
%     curve_Eperf_hdl.RT,...
%     curve_Eperf_hdl.X],...
%     regs_names(3:end),...
%     'Location','southeast');
legend([curve_Eperf_hdl.RT,...
    curve_Eperf_hdl.X],...
    {'RT','E'},...
    'Location','southeast');
legend('boxoff');
legend_size(pSize);

% % save
% img_name = [saveFolder task_nm '_pupil',corr_base_nm,'_betas_' num2str(NS),'subs.png'];
% save_fig(betas_fig, saveFolder, img_name);

%% median split figures
for iEye_type = 1:n_eyeDataType
    eye_type_nm = eye_data_type{iEye_type};
    for iReg = 3:nReg % ignore run constants
        reg_nm = regs_names{iReg};
        
        mSplit_fig.(reg_nm) = fig();
        
        %% incentive period
        subplot(1,3,1);
        low_mSplit_inc_curve = jbfill(x_inc_eye,...
            eyeXYpupil_mSplit.inc.mean_aSubs.(reg_nm).low.(eye_type_nm)' + eyeXYpupil_mSplit.inc.sem_aSubs.(reg_nm).low.(eye_type_nm)',...
            eyeXYpupil_mSplit.inc.mean_aSubs.(reg_nm).low.(eye_type_nm)' - eyeXYpupil_mSplit.inc.sem_aSubs.(reg_nm).low.(eye_type_nm)',...
            eyeXYpupil_mSplit.inc.mean_aSubs.(reg_nm).low.(eye_type_nm)',...
            'b');
        hold on;
        high_mSplit_inc_curve = jbfill(x_inc_eye,...
            eyeXYpupil_mSplit.inc.mean_aSubs.(reg_nm).high.(eye_type_nm)' + eyeXYpupil_mSplit.inc.sem_aSubs.(reg_nm).high.(eye_type_nm)',...
            eyeXYpupil_mSplit.inc.mean_aSubs.(reg_nm).high.(eye_type_nm)' - eyeXYpupil_mSplit.inc.sem_aSubs.(reg_nm).high.(eye_type_nm)',...
            eyeXYpupil_mSplit.inc.mean_aSubs.(reg_nm).high.(eye_type_nm)',...
            'r');
        switch eye_type_nm
            case 'pupilD'
                % mark significant clusters from the RFT analysis
                x_signif = x_inc_eye( pval_X_signif_inc_mSplit_pupilD.(reg_nm) );
                % mark the mean
                plot(x_signif,...
                    pval_Y_signif_inc_mSplit_pupilD.(reg_nm).low,['b',' *']);
                plot(x_signif,...
                    pval_Y_signif_inc_mSplit_pupilD.(reg_nm).high,['r',' *']);
            case 'pupilD_noLum'
                % mark significant clusters from the RFT analysis
                x_signif = x_inc_eye( pval_X_signif_inc_mSplit_pupilD_noLum.(reg_nm) );
                % mark the mean
                plot(x_signif,...
                    pval_Y_signif_inc_mSplit_pupilD_noLum.(reg_nm).low,['b',' *']);
                plot(x_signif,...
                    pval_Y_signif_inc_mSplit_pupilD_noLum.(reg_nm).high,['r',' *']);
        end
        % mark significant clusters on the top of the graph
        ypos_perc = 1/100;
        place_signif_curve( mSplit_fig.(reg_nm), x_signif, 'g', ypos_perc  );
        
        xlim([x_inc_eye(1), x_inc_eye(end)]);
        line([0 0],ylim(),'Color','k','LineWidth',lWidth); % mark onset of incentive
        line(xlim,[0 0],'Color','k','LineWidth',lWidth); % mark pupil at zero
        legend([high_mSplit_inc_curve, low_mSplit_inc_curve],...
            {['high ',reg_nm],['low ',reg_nm]},...
            'Location','southeast');
        legend('boxoff');
        xlabel('Time (ms) (incentive period)');
        switch eye_type_nm
            case 'pupilD'
                ylabel('Pupil diameter');
            case 'pupilD_noLum'
                ylabel('Pupil diameter (corrected for luminance)');
        end
        legend_size(pSize);
        
        %% effort scale period
        subplot(1,3,2);
        low_mSplit_dispE_curve = jbfill(x_dispE_eye,...
            eyeXYpupil_mSplit.dispE.mean_aSubs.(reg_nm).low.(eye_type_nm)' + eyeXYpupil_mSplit.dispE.sem_aSubs.(reg_nm).low.(eye_type_nm)',...
            eyeXYpupil_mSplit.dispE.mean_aSubs.(reg_nm).low.(eye_type_nm)' - eyeXYpupil_mSplit.dispE.sem_aSubs.(reg_nm).low.(eye_type_nm)',...
            eyeXYpupil_mSplit.dispE.mean_aSubs.(reg_nm).low.(eye_type_nm)',...
            'b');
        hold on;
        high_mSplit_dispE_curve = jbfill(x_dispE_eye,...
            eyeXYpupil_mSplit.dispE.mean_aSubs.(reg_nm).high.(eye_type_nm)' + eyeXYpupil_mSplit.dispE.sem_aSubs.(reg_nm).high.(eye_type_nm)',...
            eyeXYpupil_mSplit.dispE.mean_aSubs.(reg_nm).high.(eye_type_nm)' - eyeXYpupil_mSplit.dispE.sem_aSubs.(reg_nm).high.(eye_type_nm)',...
            eyeXYpupil_mSplit.dispE.mean_aSubs.(reg_nm).high.(eye_type_nm)',...
            'r');
        switch eye_type_nm
            case 'pupilD'
                % mark significant clusters from the RFT analysis
                x_signif = x_dispE_eye( pval_X_signif_dispE_mSplit_pupilD.(reg_nm) );
                % mark the mean
                plot(x_signif,...
                    pval_Y_signif_dispE_mSplit_pupilD.(reg_nm).low,['b',' *']);
                plot(x_signif,...
                    pval_Y_signif_dispE_mSplit_pupilD.(reg_nm).high,['r',' *']);
            case 'pupilD_noLum'
                % mark significant clusters from the RFT analysis
                x_signif = x_dispE_eye( pval_X_signif_dispE_mSplit_pupilD_noLum.(reg_nm) );
                % mark the mean
                plot(x_signif,...
                    pval_Y_signif_dispE_mSplit_pupilD_noLum.(reg_nm).low,['b',' *']);
                plot(x_signif,...
                    pval_Y_signif_dispE_mSplit_pupilD_noLum.(reg_nm).high,['r',' *']);
        end
        % mark significant clusters on the top of the graph
        ypos_perc = 1/100;
        place_signif_curve( mSplit_fig.(reg_nm), x_signif, 'g', ypos_perc  );
        
        xlim([x_dispE_eye(1), x_dispE_eye(end)]);
        line([0 0],ylim(),'Color','k','LineWidth',3); % mark onset of effort scale
        if strcmp(task_id,'G')
            line([5000 5000],ylim(),'Color','k','LineWidth',3); % mark end of effort period
        end
        line(xlim,[0 0],'Color','k','LineWidth',2); % mark pupil at zero
        legend([high_mSplit_dispE_curve, low_mSplit_dispE_curve],...
            {['high ',reg_nm],['low ',reg_nm]},...
            'Location','southeast');
        legend('boxoff');
        xlabel('Time (ms) (effort scale period)');
        switch eye_type_nm
            case 'pupilD'
                ylabel('Pupil diameter');
            case 'pupilD_noLum'
                ylabel('Pupil diameter (corrected for luminance)');
        end
        legend_size(pSize);
        
        %% time locked to effort performance
        subplot(1,3,3);
        low_mSplit_Eperf_curve = jbfill(x_Eperf_eye,...
            eyeXYpupil_mSplit.Eperf.mean_aSubs.(reg_nm).low.(eye_type_nm)' + eyeXYpupil_mSplit.Eperf.sem_aSubs.(reg_nm).low.(eye_type_nm)',...
            eyeXYpupil_mSplit.Eperf.mean_aSubs.(reg_nm).low.(eye_type_nm)' - eyeXYpupil_mSplit.Eperf.sem_aSubs.(reg_nm).low.(eye_type_nm)',...
            eyeXYpupil_mSplit.Eperf.mean_aSubs.(reg_nm).low.(eye_type_nm)',...
            'b');
        hold on;
        high_mSplit_Eperf_curve = jbfill(x_Eperf_eye,...
            eyeXYpupil_mSplit.Eperf.mean_aSubs.(reg_nm).high.(eye_type_nm)' + eyeXYpupil_mSplit.Eperf.sem_aSubs.(reg_nm).high.(eye_type_nm)',...
            eyeXYpupil_mSplit.Eperf.mean_aSubs.(reg_nm).high.(eye_type_nm)' - eyeXYpupil_mSplit.Eperf.sem_aSubs.(reg_nm).high.(eye_type_nm)',...
            eyeXYpupil_mSplit.Eperf.mean_aSubs.(reg_nm).high.(eye_type_nm)',...
            'r');
        switch eye_type_nm
            case 'pupilD'
                % mark significant clusters from the RFT analysis
                x_signif = x_Eperf_eye( pval_X_signif_Eperf_mSplit_pupilD.(reg_nm) );
                % mark the mean
                plot(x_signif,...
                    pval_Y_signif_Eperf_mSplit_pupilD.(reg_nm).low,['b',' *']);
                plot(x_signif,...
                    pval_Y_signif_Eperf_mSplit_pupilD.(reg_nm).high,['r',' *']);
            case 'pupilD_noLum'
                % mark significant clusters from the RFT analysis
                x_signif = x_Eperf_eye( pval_X_signif_Eperf_mSplit_pupilD_noLum.(reg_nm) );
                % mark the mean
                plot(x_signif,...
                    pval_Y_signif_Eperf_mSplit_pupilD_noLum.(reg_nm).low,['b',' *']);
                plot(x_signif,...
                    pval_Y_signif_Eperf_mSplit_pupilD_noLum.(reg_nm).high,['r',' *']);
        end
        % mark significant clusters on the top of the graph
        ypos_perc = 1/100;
        place_signif_curve( mSplit_fig.(reg_nm), x_signif, 'g', ypos_perc  );
        
        xlim([x_Eperf_eye(1), x_Eperf_eye(end)]);
        line([0 0],ylim(),'Color','k','LineWidth',lWidth); % mark onset of incentive
        line(xlim,[0 0],'Color','k','LineWidth',lWidth); % mark pupil at zero
        legend([high_mSplit_Eperf_curve, low_mSplit_Eperf_curve],...
            {['high ',reg_nm],['low ',reg_nm]},...
            'Location','southeast');
        legend('boxoff');
        xlabel('Time locked to effort (ms)');
        switch eye_type_nm
            case 'pupilD'
                ylabel('Pupil diameter');
            case 'pupilD_noLum'
                ylabel('Pupil diameter (corrected for luminance)');
        end
        legend_size(pSize);
        
        %         %% save
        %         img_name = [saveFolder task_nm,'_',eye_type_nm,corr_base_nm,'_' reg_nm '_mSplit_' num2str(NS),'subs.png'];
        %         save_fig(mSplit_fig.(reg_nm), saveFolder, img_name);
    end % regressor
end % eye data type


%% initialize figure
figure;
black = [0 0 0];
% X_col = [178 255 178]./255;
X_col = [124 255 124]./255;
% perf_col = [171 217 233]./255;
perf_col = [158 202 225]./255;
% RT_col = [178 178 178]./255;
RT_col = [124 124 124]./255;
pSize = 40;
switch task_id
    case 'G'
        colororder([black;perf_col]);
end
ylim_vals = [-0.25 0.25];
lWidth_graphAxis = 1;
lWidth_mean = 3;
%% add force (grip)/RT (stroop) overlap on the right of the graph
% start with this so that it's in the font and the more relevant curves are
% displayed above it
switch task_id
    case 'G'
        yyaxis right;
        jbfill(gripForceTime_m_aSubs.dispE',...
            gripForce_m_aSubs.dispE' + gripForce_sem_aSubs.dispE',...
            gripForce_m_aSubs.dispE' - gripForce_sem_aSubs.dispE',...
            gripForce_m_aSubs.dispE',...
            perf_col);
        ylabel('Force (%)');
        ylim([0 50]);
    case 'S'
        y_size = ylim_vals;
        for iPair = 1:nMaxStroop
            RT_pair_dispE = mTime_perPair_m_aSubs.dispE(iPair);
            RT_pair_dispE_top = mTime_perPair_m_aSubs.dispE(iPair) + mTime_perPair_sem_aSubs.dispE(iPair);
            RT_pair_dispE_low = mTime_perPair_m_aSubs.dispE(iPair) - mTime_perPair_sem_aSubs.dispE(iPair);
            jbfill_x([RT_pair_dispE, RT_pair_dispE],...
                [RT_pair_dispE_top RT_pair_dispE_top],...
                [RT_pair_dispE_low RT_pair_dispE_low],...
                y_size,...
                perf_col);
        end
end

%% add main result with betas
switch task_id
    case 'G'
        yyaxis left;
end
ylim(ylim_vals);
for iReg = 4:nReg
    reg_nm = regs_names{iReg};
    col_nm = colours{iReg};
    switch reg_nm
        case 'X'
            col_tmp = X_col;
        case 'RT'
            col_tmp = RT_col;
    end
    curve_dispE_hdl.(reg_nm) = jbfill(x_dispE_eye,...
        betas.dispE.mean_aSubs.(reg_nm)' + betas.dispE.sem_aSubs.(reg_nm)',...
        betas.dispE.mean_aSubs.(reg_nm)' - betas.dispE.sem_aSubs.(reg_nm)',...
        betas.dispE.mean_aSubs.(reg_nm)',...
        col_tmp);
    % mark significant clusters from the RFT analysis
    x_signif = x_dispE_eye( pval_X_signif_clusters_dispE.(reg_nm) );
    % mark the mean
    plot(x_signif,...
        pval_Y_signif_clusters_dispE.(reg_nm),...
        ' *','Color',col_tmp,'LineWidth',lWidth_mean);
%     % mark significant clusters on the top of the graph
%     ypos_perc = (1/100)*iReg;
% %     ypos_perc = ylim_vals(2) - (iReg/100)*abs(ylim_vals(2) - ylim_vals(1));
%     place_signif_curve( betas_fig, x_signif, col_tmp, ypos_perc  );
end % regresssor loop
line([0 0],ylim_vals,'Color','k','LineWidth',lWidth_graphAxis); % mark onset of effort scale
line(xlim,[0 0],'Color','k','LineWidth',lWidth_graphAxis); % mark beta at zero
if strcmp(task_id,'G')
    line([5000 5000],ylim_vals,'Color','k','LineWidth',lWidth_graphAxis); % mark end of effort period
end
ylim(ylim_vals);
xlabel('Time after scale onset (ms)');
ylabel('Regression estimate');
legend([curve_dispE_hdl.RT,...
    curve_dispE_hdl.X],...
    {'RT','E*'},...
    'Location','northwest');
legend('boxoff');
xlim([-500 6000]);
legend_size(pSize);

%% for grip task, display also graph with performance instead of force overlap
if strcmp(task_id,'G')
    betas_fig = fig();
    switch task_id
        case 'G'
            colororder([black;perf_col]);
    end
    %% add performance (grip)/RT (stroop) overlap
    switch task_id
        case 'G'
            yyaxis right;
            jbfill(gripForceTime_m_aSubs.dispE',...
                gripForceDisplay_m_aSubs.dispE' + gripForceDisplay_sem_aSubs.dispE',...
                gripForceDisplay_m_aSubs.dispE' - gripForceDisplay_sem_aSubs.dispE',...
                gripForceDisplay_m_aSubs.dispE',...
                perf_col);
            ylabel('Performance (%)');
            ylim([0 50]);
        case 'S'
            y_size = ylim();
            for iPair = 1:nMaxStroop
                RT_pair_dispE = mTime_perPair_m_aSubs.dispE(iPair);
                RT_pair_dispE_top = mTime_perPair_m_aSubs.dispE(iPair) + mTime_perPair_sem_aSubs.dispE(iPair);
                RT_pair_dispE_low = mTime_perPair_m_aSubs.dispE(iPair) - mTime_perPair_sem_aSubs.dispE(iPair);
                jbfill_x([RT_pair_dispE, RT_pair_dispE],...
                    [RT_pair_dispE_top RT_pair_dispE_top],...
                    [RT_pair_dispE_low RT_pair_dispE_low],...
                    y_size,...
                    perf_col);
            end
    end
    
    %% add main result with betas
    switch task_id
        case 'G'
            yyaxis left;
            ylim(ylim_vals);
    end
    for iReg = 4:nReg
        reg_nm = regs_names{iReg};
        switch reg_nm
            case 'X'
                col_tmp = X_col;
            case 'RT'
                col_tmp = RT_col;
        end
        curve_dispE_hdl.(reg_nm) = jbfill(x_dispE_eye,...
            betas.dispE.mean_aSubs.(reg_nm)' + betas.dispE.sem_aSubs.(reg_nm)',...
            betas.dispE.mean_aSubs.(reg_nm)' - betas.dispE.sem_aSubs.(reg_nm)',...
            betas.dispE.mean_aSubs.(reg_nm)',...
            col_tmp);
        % mark significant clusters from the RFT analysis
        x_signif = x_dispE_eye( pval_X_signif_clusters_dispE.(reg_nm) );
        % mark the mean
        plot(x_signif,...
            pval_Y_signif_clusters_dispE.(reg_nm),...
            ' *','Color',col_tmp,'LineWidth',lWidth_mean);
%         % mark significant clusters on the top of the graph
%         ypos_perc = (1/100)*iReg;
%         place_signif_curve( betas_fig, x_signif, col_tmp, ypos_perc  );
    end % regresssor loop
    line([0 0],ylim_vals,'Color','k','LineWidth',lWidth_graphAxis); % mark onset of effort scale
    line(xlim,[0 0],'Color','k','LineWidth',lWidth_graphAxis); % mark beta at zero
    if strcmp(task_id,'G')
        line([5000 5000],ylim_vals,'Color','k','LineWidth',lWidth_graphAxis); % mark end of effort period
    end
    ylim(ylim_vals);
    xlabel('Time after scale onset (ms)');
    ylabel('Regression estimate');
    
    legend([curve_dispE_hdl.RT,...
        curve_dispE_hdl.X],...
        {'RT','E*'},...
        'Location','northwest');
    legend('boxoff');
    xlim([-500 6000]);
    legend_size(pSize);
end % grip filter

%% perform t.test on the average beta across the whole relevant period
switch task_id
    case 'G'
        effort_dur = (1:n_effort_time) + time_bf_onset;
        m_betas_to_test.RT = mean(betas.dispE.RT(effort_dur,:),1,'omitnan');
        m_betas_to_test.X = mean(betas.dispE.X(effort_dur,:),1,'omitnan');
        [~,pval.RT] = ttest(m_betas_to_test.RT);
        [~,pval.X] = ttest(m_betas_to_test.X);
        mean(m_betas_to_test.X,2)
        sem(m_betas_to_test.X,2)
        mean(m_betas_to_test.RT,2)
        sem(m_betas_to_test.RT,2)
    case 'S'
        % average across 8s
        effort_dur = (1:n_effort_time) + time_bf_onset;
        m_betas_to_test.RT = mean(betas.dispE.RT(effort_dur,:),1,'omitnan');
        m_betas_to_test.X = mean(betas.dispE.X(effort_dur,:),1,'omitnan');
        [~,pval.RT] = ttest(m_betas_to_test.RT);
        [~,pval.X] = ttest(m_betas_to_test.X);
        mean(m_betas_to_test.X,2)
        sem(m_betas_to_test.X,2)
        mean(m_betas_to_test.RT,2)
        sem(m_betas_to_test.RT,2)
        
        % average across 5s
        effort_dur_bis = (1:5000) + time_bf_onset;
        m_betas_to_test_bis.RT = mean(betas.dispE.RT(effort_dur_bis,:),1,'omitnan');
        m_betas_to_test_bis.X = mean(betas.dispE.X(effort_dur_bis,:),1,'omitnan');
        [~,pval_bis.RT] = ttest(m_betas_to_test_bis.RT);
        [~,pval_bis.X] = ttest(m_betas_to_test_bis.X);
        mean(m_betas_to_test_bis.X,2)
        sem(m_betas_to_test_bis.X,2)
        mean(m_betas_to_test_bis.RT,2)
        sem(m_betas_to_test_bis.RT,2)
end

% uncorrected stats
data_to_test=data_to_test_dispE.X;
nSamples = size(data_to_test,2);
pval_uncorrected = NaN(1,nSamples);
for iSample = 1:nSamples
    [~,pval_uncorrected(iSample)] = ttest(data_to_test(:,iSample));
end % sample loop

end % function