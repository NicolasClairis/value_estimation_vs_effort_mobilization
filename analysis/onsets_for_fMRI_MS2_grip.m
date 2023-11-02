function [] = onsets_for_fMRI_MS2_grip(behaviorDir, fMRI_Dir, subid)
%[] = onsets_for_fMRI_MS2_grip(behaviorDir, fMRI_Dir, subid)
% onsets_for_fMRI_MS2_grip extracts the onsets and relevant variables for
% performing the analysis of the grip task in the Motiscan2 2017 march
% study.
%
% INPUTS
% behaviorDir: path to subject behavioral data
%
% fMRI_Dir: path to subject fMRI data
%
% subid: subject identification name (string)
%
% OUTPUTS
% (grip: : structure containing all the relevant infos saved in each subject
% folder)
%
% See also onsets_for_fMRI_MS2, onsets_for_fMRI_MS2_stroop,
% onsets_for_fMRI_MS2_RL, First_level_MS2_megaconcatenation_NicoC_batch &
% First_level_MS2_gripRP_prm

%% main infos
n_grip_runs = 2; % number of grip runs
% 2 stars required in old matlab versions but now 1 star should be enough
% grip_runnames = ls([behaviorDir, filesep, 'MBB_battery_*','*gripRP_onsets_sub',subid,'_sess*','*.mat']);
grip_runnames = ls([behaviorDir, 'MBB_battery_*gripRP_onsets_sub',subid,'_sess*.mat']);
grip_runs = nan(1,n_grip_runs);
for iGripRun = 1:n_grip_runs
    grip_runs(iGripRun) = str2double(grip_runnames(iGripRun,end-4));
end
n_totalTrials = 60;

onset_types = {'ITI','incentive','effortScale','feedback','firstPress'};
n_onsets = length( onset_types );
dur_types   = {'ITI','incentive','effortScale','feedback','perfE','perfE_until_fbk'};
n_dur = length(dur_types);
% modulators
mod_types = {'incentive','incentiveRank','absIncentive','absIncentiveRank',...
    'inc_type',...
    'trialValence',...
    'task_inc',...
    'RT_fp',...
    'perf','perf_orth','peakForce','intForce','intForcePrev','normIntForce','Fmax',...
    'sumPerf_prev',...
    'trialN',...
    'gain','totalFeedback','totalGain_prev',...
    'lum_inc'};
n_mod = length(mod_types);

% trial valence
trialTypes_valence = {'toGain','toLose'};
n_trialTypes = length(trialTypes_valence);

%% load luminance (at incentive)
lum_path = fullfile('C:','Users','clairis','Desktop','GitHub','Bat-Motiv',...
    'Bat-analyse','various','NicoC_analysis_behavior',...
    'eye_tracking_functions','evaluate_luminance','MS2');
lum_inc_vals = getfield( load([lum_path,filesep,'MS2_luminance_incentive.mat'],'lum_inc'),'lum_inc');

%% fit performance 
[ ~, perf_orth_r1, perf_orth_r2 ] = MS2_GS_perf_orth( behaviorDir, subid, grip_runs, 'grip' );

%% loop through grip runs
jGripRun = 0;
for iGripRun = grip_runs
    sessnber = num2str(iGripRun);
    
    % load onsets and duration
    jGripRun = jGripRun + 1;
    Grip_onset_loadStruct = load([behaviorDir grip_runnames(jGripRun,:)],'onset','duration');
    onset       = Grip_onset_loadStruct.onset;
    duration    = Grip_onset_loadStruct.duration;
    % load T0 for this run
    T0 = getfield(load([behaviorDir 'TTL_sub' subid '_sess' sessnber '.mat'],'T0'),'T0');
    % load behavioral modulators: Fmax used for each trial, performance for
    % the trial (in percentage), actual performance (raw force value),
    % values of incentive
    Grip_behav_loadStruct = load([behaviorDir 'global_sub_' subid '_session_' sessnber '_gripRP.mat'],...
        'peakForce','Fmax','perf','gripData',...
        'outcomeValues','task_inc','trialValence',...
        'gain','feedback');
    peakForce       = Grip_behav_loadStruct.peakForce;
    Fmax            = Grip_behav_loadStruct.Fmax;
    perf            = Grip_behav_loadStruct.perf;
    gripData        = Grip_behav_loadStruct.gripData;
    outcomeValues   = Grip_behav_loadStruct.outcomeValues;
    task_inc        = Grip_behav_loadStruct.task_inc;
    trialValence    = Grip_behav_loadStruct.trialValence;
    gain            = Grip_behav_loadStruct.gain;
    totalGain       = Grip_behav_loadStruct.feedback;
    sumPerf_prev = NaN(1,n_totalTrials);
    sumPerf_prev(1) = 0;
    for iT = 2:n_totalTrials
        sumPerf_prev(iT) = sum(perf(1:(iT-1))./100);
    end % trial loop
    if ismember(sessnber,{'2','3'})
        init_total = 0; % start at 0 euros
    elseif ismember(sessnber,{'5','6'}) % gain cumulated from previous sessions
        init_total = totalGain(1) - gain(1);
    end
    totalGain_prev  = [init_total, totalGain(1:(end-1))];
    [intForce, intForcePrev, normIntForce] = deal( NaN(1,n_totalTrials) );
    
    % extract luminance for incentive period
    lum_inc = NaN(1,n_totalTrials);
    for iTrial = 1:n_totalTrials
        switch trialValence(iTrial)
            case -1
                lum_inc(iTrial) = lum_inc_vals.loss(task_inc(iTrial));
            case 1
                lum_inc(iTrial) = lum_inc_vals.gain(task_inc(iTrial));
        end
    end
    
    % extract integral force
    runFmax = nanmax(Fmax);
    for iTrial = 1:n_totalTrials
        intForce(iTrial) = nansum( gripData.gripB{1,iTrial} );
        
        n_trial_samples = length(gripData.gripB{1,iTrial});
        normIntForce(iTrial) = intForce(iTrial)/(nansum( n_trial_samples*runFmax) ); % normalize by as if the subject did his maximal force for the whole trial
        
        if iTrial == 1
            intForcePrev(1) = 0;
        else
            intForcePrev(iTrial) = nansum( intForce(1:(iTrial - 1)) );
        end
    end
    
    %% extract T0 from onsets of interest
    grip.onset.all.ITI          = onset.cross_ITI           - T0;
    grip.onset.all.incentive    = onset.displayIncentive    - T0;
    grip.onset.all.effortScale  = onset.displayEffortScale  - T0;
    grip.onset.all.feedback     = onset.feedback            - T0;
    
    %% duration
    grip.duration.all.ITI          = duration.cross_ITI;
    grip.duration.all.incentive    = duration.displayIncentive;
    grip.duration.all.effortScale  = duration.displayEffortScale;
    grip.duration.all.feedback     = duration.feedback;
    
    %% modulators
    % extract RT
    force_min_threshold_for_RT_fp = 1;
    [grip_mod_all_RT_fp, grip.onset.all.firstPress,...
        grip_mod_all_RT_lp, grip.onset.all.endPress] = deal( NaN(1,n_totalTrials) );
    for iTrial = 1:n_totalTrials
        iGrip = 1;
        trialGripSamples = length(gripData.gripB{iTrial});
        while gripData.level{iTrial}(iGrip) <= force_min_threshold_for_RT_fp && iGrip < trialGripSamples % try to find when pressure is exerted on the grip
            iGrip = iGrip + 1;
        end
        if iGrip < trialGripSamples % if no force exerted at all in the trial => keeps NaN value, otherwise extract time of onset and RT
            grip.onset.all.firstPress(iTrial)   = gripData.time{iTrial}(iGrip) - T0;
            grip_mod_all_RT_fp(iTrial)          = gripData.time{iTrial}(iGrip) - onset.displayEffortScale(iTrial); % grip.onset.grip(ntrial) - grip.onset.effortScale(ntrial);
            
            % extract end of effort (if an effort was performed at all
            % during the current trial)
            jGrip = iGrip;
            while gripData.level{iTrial}(jGrip) > force_min_threshold_for_RT_fp && jGrip < trialGripSamples
                jGrip = jGrip + 1;
            end
            grip.onset.all.endPress(iTrial) = gripData.time{iTrial}(jGrip) - T0;
            grip_mod_all_RT_lp(iTrial)      = gripData.time{iTrial}(jGrip) - onset.displayEffortScale(iTrial);
        end
    end % trial loop
    
    % corresponding duration
    grip.duration.all.perfE             = grip.onset.all.endPress - grip.onset.all.firstPress;
    grip.duration.all.perfE_until_fbk   = grip.onset.all.feedback - grip.onset.all.firstPress;
    
    % extract incentive
    grip.mod.all.incentive          = outcomeValues(task_inc).*trialValence; % incentive value (-0.01 to -20 and +0.01 to +20)
    grip.mod.all.incentiveRank      = task_inc.*trialValence; % incentive rank (-1 to -6 and +1 to +6)
    grip.mod.all.absIncentive       = abs(grip.mod.all.incentive); % incentive value (+0.01 to +20)
    grip.mod.all.absIncentiveRank   = abs(grip.mod.all.incentiveRank); % incentive rank (+1 to +6)
    grip.mod.all.inc_type           = ismember( abs(grip.mod.all.incentiveRank), [5,6]); % 1 for bills (5€ and 20€), 0 for pieces of money
    
    % trial number
    trialN = 1:n_totalTrials;
    
    %% if no force exerted at all => exclude trial from analysis
    Grip_nonExcludedTrials_idx       = find(isnan(grip_mod_all_RT_fp) == 0); % extract all trials where some force was exerted (=non-NaN trials)
    % onsets
    for iOnset = 1:n_onsets
        curr_onset_nm = onset_types{iOnset};
        grip_onset_all.(curr_onset_nm) = grip.onset.all.(curr_onset_nm);
        grip.onset.all.(curr_onset_nm) = grip_onset_all.(curr_onset_nm)(Grip_nonExcludedTrials_idx);
    end
    % duration
    for iDur = 1:n_dur
        curr_dur_nm = dur_types{iDur};
        grip_dur_all.(curr_dur_nm) = grip.duration.all.(curr_dur_nm);
        grip.duration.all.(curr_dur_nm) = grip_dur_all.(curr_dur_nm)(Grip_nonExcludedTrials_idx);
    end
    % modulators
    grip.mod.all.RT_fp              = grip_mod_all_RT_fp(Grip_nonExcludedTrials_idx);
    grip.mod.all.incentive          = grip.mod.all.incentive(Grip_nonExcludedTrials_idx); % incentive value (-0.01 to -20 and +0.01 to +20)
    grip.mod.all.incentiveRank      = grip.mod.all.incentiveRank(Grip_nonExcludedTrials_idx); % incentive rank (-1 to -6 and +1 to +6)
    grip.mod.all.absIncentive       = grip.mod.all.absIncentive(Grip_nonExcludedTrials_idx); % incentive value (+0.01 to +20)
    grip.mod.all.absIncentiveRank   = grip.mod.all.absIncentiveRank(Grip_nonExcludedTrials_idx); % incentive rank (+1 to +6)
    grip.mod.all.inc_type           = grip.mod.all.inc_type(Grip_nonExcludedTrials_idx); % 0/1 pieces/bills
    grip.mod.all.peakForce          = peakForce(Grip_nonExcludedTrials_idx);
    grip.mod.all.intForce           = intForce(Grip_nonExcludedTrials_idx);
    grip.mod.all.intForcePrev       = intForcePrev(Grip_nonExcludedTrials_idx);
    grip.mod.all.normIntForce       = normIntForce(Grip_nonExcludedTrials_idx);
    grip.mod.all.Fmax               = Fmax(Grip_nonExcludedTrials_idx);
    grip.mod.all.perf               = perf(Grip_nonExcludedTrials_idx);
    grip.mod.all.sumPerf_prev       = sumPerf_prev(Grip_nonExcludedTrials_idx);
    switch jGripRun
        case 1
            grip.mod.all.perf_orth  = perf_orth_r1(Grip_nonExcludedTrials_idx);
        case 2
            grip.mod.all.perf_orth  = perf_orth_r2(Grip_nonExcludedTrials_idx);
    end
    grip.mod.all.task_inc           = task_inc(Grip_nonExcludedTrials_idx);
    grip.mod.all.trialValence       = trialValence(Grip_nonExcludedTrials_idx);
    grip.mod.all.gain               = gain(Grip_nonExcludedTrials_idx);
    grip.mod.all.totalFeedback      = totalGain(Grip_nonExcludedTrials_idx);
    grip.mod.all.totalGain_prev     = totalGain_prev(Grip_nonExcludedTrials_idx);
    grip.mod.all.trialN             = trialN(Grip_nonExcludedTrials_idx);
    grip.mod.all.lum_inc            = lum_inc(Grip_nonExcludedTrials_idx);
    
    
    %% extract missed trials
    Grip_excludedTrials_idx = find(isnan(grip_mod_all_RT_fp)); % if no answer was provided => exclude the trial
    if ~isempty(Grip_excludedTrials_idx)
        % onsets
        for iOnset = 1:n_onsets
            curr_onset_nm1 = onset_types{iOnset};
            curr_onset_nm2 = [curr_onset_nm1,'_missedTrial'];
            grip.onset.all.(curr_onset_nm2)      = grip_onset_all.(curr_onset_nm1)( Grip_excludedTrials_idx );
        end % onsets loop
        
        % duration
        for iDur = 1:n_dur
            curr_dur_nm1 = dur_types{iDur};
            curr_dur_nm2 = [curr_dur_nm1,'_missedTrial'];
            grip.duration.all.(curr_dur_nm2)   = grip_dur_all.(curr_dur_nm1)( Grip_excludedTrials_idx );
        end % onsets loop
    end
    
    %% separate between gain and loss trials
    % gain/loss trials : be careful that index depends on whether you
    % already excluded the NaN trials or not => be very careful on how you
    % use this to avoid any weird data point
    toGainTrials = (grip.mod.all.trialValence == 1);
    toLoseTrials = (grip.mod.all.trialValence == -1);
    
    % onsets
    for iOnset = 1:n_onsets
        curr_onset_nm = onset_types{iOnset};
        % onsets gain trials
        grip.onset.toGain.(curr_onset_nm)       = grip.onset.all.(curr_onset_nm)(toGainTrials);
        % onsets loss trials
        grip.onset.toLose.(curr_onset_nm)       = grip.onset.all.(curr_onset_nm)(toLoseTrials);
    end
    
    % duration
    for iDur = 1:n_dur
        curr_dur_nm = dur_types{iDur};
        % duration gain trials
        grip.duration.toGain.(curr_dur_nm)    = grip.duration.all.(curr_dur_nm)(toGainTrials);
        % duration loss trials
        grip.duration.toLose.(curr_dur_nm)    = grip.duration.all.(curr_dur_nm)(toLoseTrials);
    end
    
    % modulators
    for iMod = 1:n_mod
        curr_mod_nm = mod_types{iMod};
        % gain trials
        grip.mod.toGain.(curr_mod_nm) = grip.mod.all.(curr_mod_nm)(toGainTrials);
        % loss trials
        grip.mod.toLose.(curr_mod_nm) = grip.mod.all.(curr_mod_nm)(toLoseTrials);
    end % modulators loop
    
    %% extract missed trials
    Grip_excludedTrials_idx_toGain = ( isnan(grip_mod_all_RT_fp).*(trialValence == 1) ) == 1; % if no trial correctly solved at all => exclude the trial (even if a bad answer was provided)
    Grip_excludedTrials_idx_toLose = ( isnan(grip_mod_all_RT_fp).*( trialValence == -1) ) == 1; % if no trial correctly solved at all => exclude the trial (even if a bad answer was provided)
    % loop through onsets
    % to gain
    if ~isempty(Grip_excludedTrials_idx_toGain)
        % onsets
        for iOnset = 1:n_onsets
            curr_onset_nm1 = [onset_types{iOnset}];
            curr_onset_nm2 = [curr_onset_nm1,'_missedTrial'];
            grip.onset.toGain.(curr_onset_nm2)       = grip_onset_all.(curr_onset_nm1)( Grip_excludedTrials_idx_toGain );
        end % onsets loop
        
        % duration
        for iDur = 1:n_dur
            curr_dur_nm1 = dur_types{iDur};
            curr_dur_nm2 = [curr_dur_nm1,'_missedTrial'];
            grip.duration.toGain.(curr_dur_nm2)    = grip_dur_all.(curr_dur_nm1)( Grip_excludedTrials_idx_toGain );
        end % onsets loop
    end
    % to lose
    if ~isempty(Grip_excludedTrials_idx_toLose)
        % onsets
        for iOnset = 1:n_onsets
            curr_onset_nm1 = onset_types{iOnset};
            curr_onset_nm2 = [curr_onset_nm1,'_missedTrial'];
            grip.onset.toLose.(curr_onset_nm2)       = grip_onset_all.(curr_onset_nm1)( Grip_excludedTrials_idx_toLose );
        end % onsets loop
        
        % duration
        for iDur = 1:n_dur
            curr_dur_nm1 = dur_types{iDur};
            curr_dur_nm2 = [curr_dur_nm1,'_missedTrial'];
            grip.duration.toLose.(curr_dur_nm2)    = grip_dur_all.(curr_dur_nm1)( Grip_excludedTrials_idx_toLose );
        end % duration loop
    end
    
    %% pool all data which is relevant for the GLM together in a new file
    save([fMRI_Dir 'onsets_sub',subid,'_grip_run',sessnber,'.mat'],'grip');
    
    %% clear grip structure before next run to avoid interferences
    clear('grip');
end % runs

end % function