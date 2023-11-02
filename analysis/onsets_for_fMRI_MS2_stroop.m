function [] = onsets_for_fMRI_MS2_stroop(behaviorDir, fMRI_Dir,...
    subid)
%[ stroop ] = onsets_for_fMRI_MS2_stroop(behaviorDir, fMRI_Dir, subid)
% onsets_for_fMRI_MS2_stroop extracts the onsets and relevant variables for
% performing the analysis of the stroop tasks in the Motiscan2 2017 march
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
% (stroop: structure containing all the relevant infos saved in each subject
% folder)
%
% See also onsets_for_fMRI_MS2, onsets_for_fMRI_MS2_grip,
% onsets_for_fMRI_MS2_RL, First_level_MS2_megaconcatenation_NicoC_batch &
% First_level_MS2_mentalRP_prm

%% main infos
n_stroop_runs = 2; % number of stroop runs
% 2 stars required in old matlab versions but now 1 star should be enough
% stroop_runnames = ls([behaviorDir, filesep, 'MBB_battery_*','*mentalRP_onsets_sub',subid,'_sess*','*.mat']);
stroop_runnames = ls([behaviorDir, 'MBB_battery_*mentalRP_onsets_sub',subid,'_sess*.mat']);
stroop_runs = nan(1,n_stroop_runs);
for iStroopRun = 1:n_stroop_runs
    stroop_runs(iStroopRun) = str2double(stroop_runnames(iStroopRun,end-4));
end
n_totalTrials = 60;
n_totalPairs_per_trial = 10;

%% main parameters of interest
% onsets
onset_types = {'ITI','incentive','effortScale','firstPress','feedback'};
n_onsets = length(onset_types);
% modulators
mod_types = {'incentive','incentiveRank','absIncentive','absIncentiveRank',...
    'inc_type',...
    'trialValence',...
    'task_inc',...
    'RT_fp','mRT',...
    'mRT_without_fPair','mRT_without_Errors','mRT_without_fPair_without_Errors',...
    'RT_per_pair','RT_per_pair_first_answer',...
    'diff_per_pair',...
    'perf','perf_orth','n_errorTrials','n_errors_per_trial','perc_errors_per_trial',...
    'sumPerf_prev',...
    'first_answer_per_pair_error_or_correct',...
    'pairs_with_errors','n_pairs_solved_with_errors',...
    'n_pairs_solved','n_cong','n_incong','perc_incong',...
    'n_cong_with_error','n_cong_without_error','n_incong_with_error','n_incong_without_error',...
    'congIncong_per_pair','nDist_per_pair','cumul_n_incong_pairs_solved',...
    'm_nDist','m_inverse_nDist',...
    'trialN',...
    'gain','totalFeedback','totalGain_prev',...
    'lum_inc'};
n_mod = length(mod_types);

% trial valence
trialTypes_valence = {'toGain','toLose'};
n_trialTypes = length(trialTypes_valence);

n_min_solved_pairs = 1; % how many pairs should have been solved at minimum so that the trial is considered?
% 0 = trial not even tried at all, 1= not more than one pair solved,
% necessary if you want to consider mean(RT) of all pairs but first

%% load luminance (at incentive)
lum_path = fullfile('C:','Users','clairis','Desktop','GitHub','Bat-Motiv',...
    'Bat-analyse','various','NicoC_analysis_behavior',...
    'eye_tracking_functions','evaluate_luminance','MS2');
lum_inc_vals = getfield( load([lum_path,filesep,'MS2_luminance_incentive.mat'],'lum_inc'),'lum_inc');

%% fit performance 
[ ~, perf_orth_r1, perf_orth_r2 ] = MS2_GS_perf_orth( behaviorDir, subid, stroop_runs, 'mental' );

%% loop through stroop runs
jStroopRun = 0;
for iStroopRun = stroop_runs
    sessnber = num2str(iStroopRun);
    
    % load onsets and duration
    jStroopRun = jStroopRun + 1;
    Stroop_onset_loadStruct = load([behaviorDir stroop_runnames(jStroopRun,:)],'onset','duration');
    onset = Stroop_onset_loadStruct.onset;
    duration = Stroop_onset_loadStruct.duration;
    % load T0 for this run
    T0 = getfield(load([behaviorDir 'TTL_sub' subid '_sess' sessnber '.mat'],'T0'),'T0');
    
    % load behavioral modulators: number of pairs solved in each trial,
    % values of incentive
    Stroop_behav_loadStruct = load([behaviorDir 'global_sub_' subid '_session_' sessnber '_mentalRP.mat'],...
        'taskPerf','totalAnswers','goodAnswers','taskNumbers',...
        'data_Task',...
        'outcomeValues','task_inc','trialValence',...
        'gain','feedback');
    perf            = Stroop_behav_loadStruct.taskPerf;
    totalAnswers    = Stroop_behav_loadStruct.totalAnswers;
    goodAnswers     = Stroop_behav_loadStruct.goodAnswers;
    data_Task       = Stroop_behav_loadStruct.data_Task;
    taskNumbers     = Stroop_behav_loadStruct.taskNumbers;
    outcomeValues   = Stroop_behav_loadStruct.outcomeValues;
    task_inc        = Stroop_behav_loadStruct.task_inc;
    trialValence    = Stroop_behav_loadStruct.trialValence;
    gain            = Stroop_behav_loadStruct.gain;
    totalGain       = Stroop_behav_loadStruct.feedback;
    if ismember(sessnber,{'2','3'})
        init_total = 0; % start at 0 euros
    elseif ismember(sessnber,{'5','6'}) % gain cumulated from previous sessions
        init_total = totalGain(1) - gain(1);
    end
    totalGain_prev  = [init_total, totalGain(1:(end-1))];
    
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
    
    %% extract T0 from onsets of interest
    stroop.onset.all.ITI                = onset.cross_ITI - T0;
    stroop.onset.all.incentive          = onset.displayIncentive - T0;
    stroop.onset.all.effortScale        = onset.displayEffortScale - T0;
    stroop.onset.all.feedback           = onset.feedback - T0;
    
    %% duration
    stroop.duration.all.ITI                = duration.cross_ITI;
    stroop.duration.all.incentive          = duration.displayIncentive;
    stroop.duration.all.effortScale        = duration.displayEffortScale;
    stroop.duration.all.feedback           = duration.feedback;
    
    %% modulators
    % extract RT first press and mean RT per trial across pairs
    stroop.onset.all.firstPress = NaN(n_totalTrials,1);
    [stroop.mod.all.RT_fp,...
        stroop.mod.all.mRT,...
        stroop.mod.all.mRT_without_fPair,...
        stroop.mod.all.mRT_without_Errors,...
        stroop.mod.all.mRT_without_fPair_without_Errors] = deal( NaN(1,n_totalTrials) );
    [stroop.mod.all.RT_per_pair,... % this will consider the RT for the final answer
        stroop.mod.all.RT_per_pair_first_answer,... % this will consider the first RT for each pair (including errors)
        stroop.mod.all.first_answer_per_pair_error_or_correct,... % is the first answer correct (1) or wrong(0)?
        stroop.mod.all.congIncong_per_pair,...
        stroop.mod.all.nDist_per_pair,...
        stroop.mod.all.diff_per_pair] = deal(NaN(n_totalPairs_per_trial, n_totalTrials));
    firstPairPress = 1;
    for iTrial = 1:n_totalTrials
        if goodAnswers(iTrial) > n_min_solved_pairs % if no pair solved at all for this trial = exclude the trial (even if the subject answered with a wrong answer)
            
            % first press onset and RT
            stroop.onset.all.firstPress(iTrial) = data_Task.answertime{iTrial}(firstPairPress) - T0;
            curr_RT_fp = data_Task.answertime{iTrial}(firstPairPress) - onset.displayEffortScale(iTrial);
            stroop.mod.all.RT_fp(iTrial) = curr_RT_fp;
            
            % RT per pair (for correct answer)
            curr_RT = diff(data_Task.answertime{iTrial});
            curr_RT_bis = [curr_RT_fp, curr_RT]; % add RT first pair
            curr_correctAnswers = data_Task.answer{iTrial} == 1; % correct answers
            curr_RT_ter = curr_RT_bis( curr_correctAnswers ); % filter errors
            stroop.mod.all.RT_per_pair(1:length(curr_RT_ter), iTrial) = curr_RT_ter;
            
            jPair = 0;
            for iA = 1:length(curr_RT_bis)
                if (iA == 1) ||... % first answer
                        (curr_correctAnswers(iA - 1) == 1) % if previous answer was correct, then this is the first answer for the next pair (independently of correct or error)
                    jPair = jPair + 1;
                    curr_RT_first_answer = curr_RT_bis(iA);
                    stroop.mod.all.RT_per_pair_first_answer(jPair, iTrial) = curr_RT_first_answer;
                    
                    stroop.mod.all.first_answer_per_pair_error_or_correct(jPair, iTrial) = curr_correctAnswers(iA); % 0 if error, 1 if correct
                end
            end
            
            % mean RT across pairs (except first pair, including errors)
            curr_mRT = nanmean(curr_RT);
            stroop.mod.all.mRT(iTrial) = curr_mRT;
            
            % mean RT across pairs (including first pair)
            curr_mRT_bis = nanmean(curr_RT_bis);
            stroop.mod.all.mRT_without_fPair(iTrial) = curr_mRT_bis;
            
             % mean RT across pairs (including first pair but without
            % errors)
            curr_mRT_ter = nanmean(curr_RT_ter);
            stroop.mod.all.mRT_without_Errors(iTrial) = curr_mRT_ter;
            
            % mean RT across pairs (except first pair AND errors)
            curr_mRT_quater = nanmean(curr_RT_ter(2:end));
            stroop.mod.all.mRT_without_fPair_without_Errors(iTrial) = curr_mRT_quater;
            
        end % filter trials
    end % trial loop
    
    % duration from first pair RT to feedback
    stroop.duration.all.firstPress = stroop.onset.all.feedback - stroop.onset.all.firstPress;
    
    % extract error pairs
    errorPairs = NaN(n_totalPairs_per_trial, n_totalTrials);% identify for which pairs of each trial errors were made if so
    for iTrial = 1:n_totalTrials
        trial_answers = data_Task.answer{iTrial};
        jPair = 1;
        for iAnswer = 1:length(trial_answers)
            switch trial_answers(iAnswer) 
                case 0 % error
                    errorPairs(jPair, iTrial) = 1;
                case 1 % correct answer
                    if isnan(errorPairs(jPair, iTrial))
                        errorPairs(jPair, iTrial) = 0; % pair correctly answered without error
                    end
                    jPair = jPair + 1;
            end
        end % answer loop
    end % trial loop
    stroop.mod.all.pairs_with_errors = errorPairs;
    stroop.mod.all.n_pairs_solved_with_errors = nansum(errorPairs, 1); % sum to see how many pairs contain at least 1 error/trial
    
    % incentives
    % incentive value and incentive rank
    stroop.mod.all.incentive      = outcomeValues(task_inc).*trialValence;
    stroop.mod.all.incentiveRank  = task_inc.*trialValence;
    stroop.mod.all.absIncentive      = abs(stroop.mod.all.incentive);
    stroop.mod.all.absIncentiveRank  = abs(stroop.mod.all.incentiveRank);
    stroop.mod.all.inc_type       = ismember( abs(stroop.mod.all.incentiveRank), [5,6]);
    
    % performance
    stroop.mod.all.perf           = perf;
    % sum performance previous trials
    sumPerf_prev = NaN(1,n_totalTrials);
    sumPerf_prev(1) = 0;
    for iT = 2:n_totalTrials
        sumPerf_prev(iT) = sum(perf(1:(iT-1))./100);
    end % trial loop
    stroop.mod.all.sumPerf_prev = sumPerf_prev;
    % performance orthogonalized
    switch jStroopRun
        case 1
            stroop.mod.all.perf_orth  = perf_orth_r1';
        case 2
            stroop.mod.all.perf_orth  = perf_orth_r2';
    end
    
    % trial type
    stroop.mod.all.trialValence   = trialValence;
    stroop.mod.all.task_inc       = task_inc;
    
    % gains
    stroop.mod.all.gain           = gain;
    stroop.mod.all.totalFeedback  = totalGain;
    stroop.mod.all.totalGain_prev  = totalGain_prev;
    
    % number of pairs solved
    stroop.mod.all.n_pairs_solved = goodAnswers;
    
    % luminance
    stroop.mod.all.lum_inc = lum_inc;
    
    % trial number
    stroop.mod.all.trialN = 1:n_totalTrials;
    
    % extract list of numbers and fontsize for each pair for each trial
    [n_cong_pairs_solved, n_incong_pairs_solved,...
        n_cong_pairs_solved_without_any_error, n_incong_pairs_solved_without_any_error,...
        n_incong_pairs_solved_with_an_error, n_cong_pairs_solved_with_an_error,...
        cumul_n_incong_pairs_solved] = deal( zeros(1,n_totalTrials) );
    [mean_dNbers_pairs_solved, mean_dNbers_pairs_solved_bis] = deal( NaN(1,n_totalTrials) );
    for iTrial = 1:n_totalTrials
        n_currPairs_solved = goodAnswers(iTrial);
        curr_trials_solved = 1:n_currPairs_solved;
        % extract fontsize and list of numbers for the current trial
        fontSize_curr_trial     = taskNumbers{1,iTrial}.fontSize(:,curr_trials_solved);
        list_nbers_curr_trial   = taskNumbers{1,iTrial}.list(:,curr_trials_solved);
        
        dNbers_curr_trial = NaN(1,n_currPairs_solved);
        
        for iPair = 1:n_currPairs_solved % loop through pairs solved
            
            % sort congruent/incongruent pairs
            if (fontSize_curr_trial(1,iPair) < fontSize_curr_trial(2,iPair) &&...
                    list_nbers_curr_trial(1,iPair) < list_nbers_curr_trial(2,iPair) ) ||...
                    (fontSize_curr_trial(1,iPair) > fontSize_curr_trial(2,iPair) &&...
                    list_nbers_curr_trial(1,iPair) > list_nbers_curr_trial(2,iPair) ) % congruent pair
                
                cong_curr_pair   = 1;
                switch errorPairs(iPair, iTrial)
                    case 0
                        cong_error_currPair = 0;
                        incong_error_currPair = 0;
                    case 1
                        cong_error_currPair = 1;
                        incong_error_currPair = 0;
                end
                
            elseif (fontSize_curr_trial(1,iPair) < fontSize_curr_trial(2,iPair) &&...
                    list_nbers_curr_trial(1,iPair) > list_nbers_curr_trial(2,iPair) ) ||...
                    (fontSize_curr_trial(1,iPair) > fontSize_curr_trial(2,iPair) &&...
                    list_nbers_curr_trial(1,iPair) < list_nbers_curr_trial(2,iPair) )% incongruent pair
                
                cong_curr_pair   = 0;
                switch errorPairs(iPair, iTrial)
                    case 0
                        cong_error_currPair = 0;
                        incong_error_currPair = 0;
                    case 1
                        cong_error_currPair = 0;
                        incong_error_currPair = 1;
                end
            end % congruent/incongruent trials sorting
            incong_curr_pair = 1 - cong_curr_pair;
            n_cong_pairs_solved(iTrial)     = n_cong_pairs_solved(iTrial)   + cong_curr_pair;
            n_incong_pairs_solved(iTrial) 	= n_incong_pairs_solved(iTrial) + incong_curr_pair;
            
            n_cong_pairs_solved_without_any_error(iTrial)   = n_cong_pairs_solved_without_any_error(iTrial) + cong_curr_pair - cong_error_currPair;
            n_cong_pairs_solved_with_an_error(iTrial)       = n_cong_pairs_solved_with_an_error(iTrial) + cong_error_currPair;
            n_incong_pairs_solved_without_any_error(iTrial) = n_incong_pairs_solved_without_any_error(iTrial) + incong_curr_pair - incong_error_currPair;
            n_incong_pairs_solved_with_an_error(iTrial)     = n_incong_pairs_solved_with_an_error(iTrial) + incong_error_currPair;
            
            stroop.mod.all.congIncong_per_pair(iPair, iTrial) = cong_curr_pair; % 1 if congruent, 0 if incongruent
            
            % extract distance between the numbers
            dN_currPair = abs( list_nbers_curr_trial(1,iPair) - list_nbers_curr_trial(2,iPair) );
            dNbers_curr_trial(iPair) = dN_currPair; % absolute number difference between two options
            stroop.mod.all.nDist_per_pair(iPair, iTrial) = dN_currPair;
            
            % give 1 difficulty level per pair (according to
            % congruent/incongruent & to distance between numbers)
            switch cong_curr_pair % congruent pairs
                case 1
                    switch dN_currPair
                        case 5
                            stroop.mod.all.diff_per_pair(iPair, iTrial) = 1;
                        case 4
                            stroop.mod.all.diff_per_pair(iPair, iTrial) = 2;
                        case 3
                            stroop.mod.all.diff_per_pair(iPair, iTrial) = 3;
                        case 2
                            stroop.mod.all.diff_per_pair(iPair, iTrial) = 4;
                        case 1
                            stroop.mod.all.diff_per_pair(iPair, iTrial) = 5;
                    end
                case 0 % incongruent pairs
                    switch dN_currPair
                        case 5
                            stroop.mod.all.diff_per_pair(iPair, iTrial) = 6;
                        case 4
                            stroop.mod.all.diff_per_pair(iPair, iTrial) = 7;
                        case 3
                            stroop.mod.all.diff_per_pair(iPair, iTrial) = 8;
                        case 2
                            stroop.mod.all.diff_per_pair(iPair, iTrial) = 9;
                        case 1
                            stroop.mod.all.diff_per_pair(iPair, iTrial) = 10;
                    end
            end
        end % pairs solved loop
        
        % sanity check
        if ((n_cong_pairs_solved_without_any_error(iTrial) + n_cong_pairs_solved_with_an_error(iTrial)) ~= n_cong_pairs_solved(iTrial)) ||...
                ((n_incong_pairs_solved_without_any_error(iTrial) + n_incong_pairs_solved_with_an_error(iTrial)) ~= n_incong_pairs_solved(iTrial)) ||...
                (n_cong_pairs_solved(iTrial) + n_incong_pairs_solved(iTrial) ~= n_currPairs_solved)
            error(['problem in extraction of congruent/incongruent pairs solved - subject ',subid,' trial ',num2str(iTrial)]);
        end
        
        % cumulated number of incongruent pairs solved until this trial
        if iTrial == 1
            cumul_n_incong_pairs_solved(iTrial) = 0;
        else
            cumul_n_incong_pairs_solved(iTrial) = nansum(n_incong_pairs_solved(1:(iTrial - 1)));
        end
        
        % average distance between numbers for each trial
        mean_dNbers_pairs_solved(iTrial)        = nanmean( dNbers_curr_trial );
        mean_dNbers_pairs_solved_bis(iTrial)    = nanmean( 1./dNbers_curr_trial );
        
    end % trial loop
    
    % check up
    if any(goodAnswers ~= n_cong_pairs_solved + n_incong_pairs_solved)
        error('problem with number of pairs solved and congruent/incongruent pairs solved, numbers don''t match.');
    end
    
    % number of errors
    % binary variable to see whether at least one error was made during
    % a given trial or not
    errorTrials_idx = (totalAnswers - goodAnswers > 0); % when equal zero = no errors, when higher than zero = at least one error
    stroop.mod.all.n_errorTrials = double(errorTrials_idx); % convert from logical to numeric value (otherwise will create bugs in 1st level
    % number of errors made in each trial across pairs
    n_errors_per_trial = totalAnswers - goodAnswers;
    stroop.mod.all.n_errors_per_trial = n_errors_per_trial;
    % percentage of errors made in each trial
    stroop.mod.all.perc_errors_per_trial = n_errors_per_trial./totalAnswers;
    
    % number of congruent pairs solved
    stroop.mod.all.n_cong               = n_cong_pairs_solved;
    stroop.mod.all.n_cong_with_error    = n_cong_pairs_solved_with_an_error;
    stroop.mod.all.n_cong_without_error = n_cong_pairs_solved_without_any_error;
    % number of incongruent pairs solved
    stroop.mod.all.n_incong                 = n_incong_pairs_solved;
    stroop.mod.all.n_incong_with_error      = n_incong_pairs_solved_with_an_error;
    stroop.mod.all.n_incong_without_error   = n_incong_pairs_solved_without_any_error;
    % incongruent pairs percentage among pairs solved
    stroop.mod.all.perc_incong          = n_incong_pairs_solved./goodAnswers;
    % cumulated number of incongruent pairs solved before each trial
    stroop.mod.all.cumul_n_incong_pairs_solved = cumul_n_incong_pairs_solved;
    
    % distance between numbers
    stroop.mod.all.m_nDist              = mean_dNbers_pairs_solved;
    stroop.mod.all.m_inverse_nDist      = mean_dNbers_pairs_solved_bis;
    
    %% split trials according to if contain or not at least one error
    noErrorTrials_idx = (totalAnswers - goodAnswers == 0);
    
    % loop through onsets & durations
    for iOnset = 1:n_onsets
        curr_onset_nm = onset_types{iOnset};
        if sum(errorTrials_idx) > 0
            stroop.onset.all_ErrorTrials.(curr_onset_nm)    = stroop.onset.all.(curr_onset_nm)( errorTrials_idx );
            stroop.duration.all_ErrorTrials.(curr_onset_nm) = stroop.duration.all.(curr_onset_nm)( errorTrials_idx );
        else
            stroop.onset.all_ErrorTrials.(curr_onset_nm)    = [];
            stroop.duration.all_ErrorTrials.(curr_onset_nm) = [];
        end
        if sum(noErrorTrials_idx) > 0
            stroop.onset.all_noErrorTrials.(curr_onset_nm)      = stroop.onset.all.(curr_onset_nm)( noErrorTrials_idx );
            stroop.duration.all_noErrorTrials.(curr_onset_nm)   = stroop.duration.all.(curr_onset_nm)( noErrorTrials_idx );
        else
            stroop.onset.all_noErrorTrials.(curr_onset_nm)      = [];
            stroop.duration.all_noErrorTrials.(curr_onset_nm)   = [];
        end
    end % onsets loop
    
    % loop through modulators
    for iMod = 1:n_mod
        curr_mod_nm = mod_types{iMod};
        if sum(errorTrials_idx) > 0
            stroop.mod.all_ErrorTrials.(curr_mod_nm) = stroop.mod.all.(curr_mod_nm)(:, errorTrials_idx );
        else
            stroop.mod.all_ErrorTrials.(curr_mod_nm) = [];
        end
        if sum(noErrorTrials_idx) > 0
            stroop.mod.all_noErrorTrials.(curr_mod_nm) = stroop.mod.all.(curr_mod_nm)(:, noErrorTrials_idx );
        else
            stroop.mod.all_noErrorTrials.(curr_mod_nm) = [];
        end
    end % modulators loop
    
    %% exclude trials with no answer at all
    Stroop_nonExcludedTrials_idx = (goodAnswers > n_min_solved_pairs);
    % if no trial correctly solved at all => exclude the trial (even if a bad answer was provided)
    % if only one pair solved, exclude the trial as well because will
    % create NaN values for some regressors
    
    % same but for no error/at least one error case
    goodAnswers_ErrorTrials = goodAnswers(errorTrials_idx);
    goodAnswers_noErrorTrials = goodAnswers(noErrorTrials_idx);
    Stroop_nonExcludedTrials_errorTrials_idx = (goodAnswers_ErrorTrials > n_min_solved_pairs);
    Stroop_nonExcludedTrials_noErrorTrials_idx = (goodAnswers_noErrorTrials > n_min_solved_pairs);
    
    % loop through onsets & durations
    for iOnset = 1:n_onsets
        curr_onset_nm = onset_types{iOnset};
        
        %% keep data temporarily to get the values for the missed trials
        % onsets
        stroop_onset_all.(curr_onset_nm)            = stroop.onset.all.(curr_onset_nm);
        stroop_onset_toGain.(curr_onset_nm)         = stroop_onset_all.(curr_onset_nm)(trialValence == 1);
        stroop_onset_toLose.(curr_onset_nm)         = stroop_onset_all.(curr_onset_nm)(trialValence == -1);
        stroop_onset_ErrorTrials.(curr_onset_nm)    = stroop.onset.all_ErrorTrials.(curr_onset_nm);
        stroop_onset_noErrorTrials.(curr_onset_nm)  = stroop.onset.all_noErrorTrials.(curr_onset_nm);
        % durations
        stroop_dur_all.(curr_onset_nm)              = stroop.duration.all.(curr_onset_nm);
        stroop_dur_toGain.(curr_onset_nm)           = stroop_dur_all.(curr_onset_nm)(trialValence == 1);
        stroop_dur_toLose.(curr_onset_nm)           = stroop_dur_all.(curr_onset_nm)(trialValence == -1);
        stroop_dur_ErrorTrials.(curr_onset_nm)      = stroop.duration.all_ErrorTrials.(curr_onset_nm);
        stroop_dur_noErrorTrials.(curr_onset_nm)    = stroop.duration.all_noErrorTrials.(curr_onset_nm);
        
        %% exclude bad trials
        stroop.onset.all.(curr_onset_nm)        = stroop_onset_all.(curr_onset_nm)( Stroop_nonExcludedTrials_idx );
        stroop.duration.all.(curr_onset_nm)     = stroop_dur_all.(curr_onset_nm)( Stroop_nonExcludedTrials_idx );
        
        if sum(Stroop_nonExcludedTrials_errorTrials_idx) > 0
            stroop.onset.all_ErrorTrials.(curr_onset_nm)    = stroop_onset_ErrorTrials.(curr_onset_nm)( Stroop_nonExcludedTrials_errorTrials_idx );
            stroop.duration.all_ErrorTrials.(curr_onset_nm) = stroop_dur_ErrorTrials.(curr_onset_nm)( Stroop_nonExcludedTrials_errorTrials_idx );
        else
            stroop.onset.all_ErrorTrials.(curr_onset_nm) = [];
            stroop.duration.all_ErrorTrials.(curr_onset_nm) = [];
        end
        if sum(Stroop_nonExcludedTrials_noErrorTrials_idx) > 0
            stroop.onset.all_noErrorTrials.(curr_onset_nm)      = stroop_onset_noErrorTrials.(curr_onset_nm)( Stroop_nonExcludedTrials_noErrorTrials_idx );
            stroop.duration.all_noErrorTrials.(curr_onset_nm)   = stroop_dur_noErrorTrials.(curr_onset_nm)( Stroop_nonExcludedTrials_noErrorTrials_idx );
        else
            stroop.onset.all_noErrorTrials.(curr_onset_nm) = [];
            stroop.duration.all_noErrorTrials.(curr_onset_nm) = [];
        end
    end % onsets loop
    
    % loop through modulators
    for iMod = 1:n_mod
        curr_mod_nm = mod_types{iMod};
%         stroop_mod_all.(curr_mod_nm) = stroop.mod.all.(curr_mod_nm);
        stroop.mod.all.(curr_mod_nm) = stroop.mod.all.(curr_mod_nm)(:, Stroop_nonExcludedTrials_idx );
        
        if sum(Stroop_nonExcludedTrials_errorTrials_idx) > 0
            stroop.mod.all_ErrorTrials.(curr_mod_nm) = stroop.mod.all_ErrorTrials.(curr_mod_nm)(:, Stroop_nonExcludedTrials_errorTrials_idx );
        else
            stroop.mod.all_ErrorTrials.(curr_mod_nm) = [];
        end
        if sum(Stroop_nonExcludedTrials_noErrorTrials_idx) > 0
            stroop.mod.all_noErrorTrials.(curr_mod_nm) = stroop.mod.all_noErrorTrials.(curr_mod_nm)(:, Stroop_nonExcludedTrials_noErrorTrials_idx );
        else
            stroop.mod.all_noErrorTrials.(curr_mod_nm) = [];
        end
    end % modulators loop
    
    %% extract missed trials
    Stroop_excludedTrials_idx = (goodAnswers <= n_min_solved_pairs); % if no pair correctly solved at all => exclude the trial (even if a bad answer was provided)
    % loop through onsets
    if ~isempty(Stroop_excludedTrials_idx)
        for iOnset = 1:n_onsets
            curr_onset_nm = onset_types{iOnset};
            curr_onset_missed_nm = [onset_types{iOnset},'_missedTrial'];
            stroop.onset.all.(curr_onset_missed_nm)     = stroop_onset_all.(curr_onset_nm)( Stroop_excludedTrials_idx );
            stroop.duration.all.(curr_onset_missed_nm)  = stroop_dur_all.(curr_onset_nm)( Stroop_excludedTrials_idx );
        end % onsets loop
    end
    
    % same for error trials
    Stroop_errorTrials_excludedTrials_idx = (goodAnswers_ErrorTrials <= n_min_solved_pairs); % if no pair correctly solved at all => exclude the trial (even if a bad answer was provided)
    % loop through onsets
    if ~isempty(Stroop_errorTrials_excludedTrials_idx)
        for iOnset = 1:n_onsets
            curr_onset_nm = onset_types{iOnset};
            curr_onset_missed_nm = [onset_types{iOnset},'_missedTrial'];
            stroop.onset.all_ErrorTrials.(curr_onset_missed_nm)     = stroop_onset_ErrorTrials.(curr_onset_nm)( Stroop_errorTrials_excludedTrials_idx );
            stroop.duration.all_ErrorTrials.(curr_onset_missed_nm)  = stroop_dur_ErrorTrials.(curr_onset_nm)( Stroop_errorTrials_excludedTrials_idx );
        end % onsets loop
    end
    
    % same for no-error trials
    Stroop_noErrorTrials_excludedTrials_idx = (goodAnswers_noErrorTrials <= n_min_solved_pairs); % if no pair correctly solved at all => exclude the trial (even if a bad answer was provided)
    % loop through onsets
    if ~isempty(Stroop_noErrorTrials_excludedTrials_idx)
        for iOnset = 1:n_onsets
            curr_onset_nm = onset_types{iOnset};
            curr_onset_missed_nm = [onset_types{iOnset},'_missedTrial'];
            stroop.onset.all_noErrorTrials.(curr_onset_missed_nm)       = stroop_onset_noErrorTrials.(curr_onset_nm)( Stroop_noErrorTrials_excludedTrials_idx );
            stroop.duration.all_noErrorTrials.(curr_onset_missed_nm)    = stroop_dur_noErrorTrials.(curr_onset_nm)( Stroop_noErrorTrials_excludedTrials_idx );
        end % onsets loop
    end
    
    %% separate to-gain and to-lose trials
    
    % gain/loss trials
    toGainTrials = (stroop.mod.all.trialValence == 1);
    toLoseTrials = (stroop.mod.all.trialValence == -1);
    
    % same for error/no-error trials
    trialValence_ErrorTrials = stroop.mod.all_ErrorTrials.trialValence;
    trialValence_noErrorTrials = stroop.mod.all_noErrorTrials.trialValence;
    toGain_ErrorTrials      = (trialValence_ErrorTrials == 1);
    toLose_ErrorTrials      = (trialValence_ErrorTrials == -1);
    toGain_noErrorTrials    = (trialValence_noErrorTrials == 1);
    toLose_noErrorTrials    = (trialValence_noErrorTrials == -1);
    
    for iTrialType = 1:n_trialTypes
        
        % select trial type and corresponding trials index
        curr_trial_type             = trialTypes_valence{iTrialType};
        curr_trial_type_errors      = [curr_trial_type,'_ErrorTrials'];
        curr_trial_type_noErrors    = [curr_trial_type,'_noErrorTrials'];
        switch curr_trial_type
            case 'toGain'
                trials_idx          = toGainTrials;
                valence_errorTrials_idx     = toGain_ErrorTrials;
                valence_noErrorTrials_idx   = toGain_noErrorTrials;
            case 'toLose'
                trials_idx          = toLoseTrials;
                valence_errorTrials_idx     = toLose_ErrorTrials;
                valence_noErrorTrials_idx   = toLose_noErrorTrials;
        end
        
        % loop through onsets
        for iOnset = 1:n_onsets
            curr_onset_nm = onset_types{iOnset};
            stroop.onset.(curr_trial_type).(curr_onset_nm)                  = stroop.onset.all.(curr_onset_nm)( trials_idx );
            stroop.duration.(curr_trial_type).(curr_onset_nm)               = stroop.duration.all.(curr_onset_nm)( trials_idx );
            
            if ~isempty(valence_errorTrials_idx)
                stroop.onset.(curr_trial_type_errors).(curr_onset_nm)       = stroop.onset.all_ErrorTrials.(curr_onset_nm)( valence_errorTrials_idx );
                stroop.duration.(curr_trial_type_errors).(curr_onset_nm)    = stroop.duration.all_ErrorTrials.(curr_onset_nm)( valence_errorTrials_idx );
            else
                stroop.onset.(curr_trial_type_errors).(curr_onset_nm)       = [];
                stroop.duration.(curr_trial_type_errors).(curr_onset_nm)    = [];
            end
            if ~isempty(valence_noErrorTrials_idx)
                stroop.onset.(curr_trial_type_noErrors).(curr_onset_nm)     = stroop.onset.all_noErrorTrials.(curr_onset_nm)( valence_noErrorTrials_idx );
                stroop.duration.(curr_trial_type_noErrors).(curr_onset_nm)  = stroop.duration.all_noErrorTrials.(curr_onset_nm)( valence_noErrorTrials_idx );
            else
                stroop.onset.(curr_trial_type_noErrors).(curr_onset_nm)     = [];
                stroop.duration.(curr_trial_type_noErrors).(curr_onset_nm)  = [];
            end
        end % onsets loop
        
        % loop through modulators
        for iMod = 1:n_mod
            curr_mod_nm = mod_types{iMod};
            stroop.mod.(curr_trial_type).(curr_mod_nm) = stroop.mod.all.(curr_mod_nm)(:, trials_idx );
            
            if ~isempty(valence_errorTrials_idx)
                stroop.mod.(curr_trial_type_errors).(curr_mod_nm) = stroop.mod.all_ErrorTrials.(curr_mod_nm)(:, valence_errorTrials_idx );
            else
                stroop.mod.(curr_trial_type_errors).(curr_mod_nm)   = [];
            end
            if ~isempty(valence_noErrorTrials_idx)
                stroop.mod.(curr_trial_type_noErrors).(curr_mod_nm) = stroop.mod.all_noErrorTrials.(curr_mod_nm)(:, valence_noErrorTrials_idx );
            else
                stroop.mod.(curr_trial_type_noErrors).(curr_mod_nm) = [];
            end
        end % modulators loop
        
    end % trial type loop
    
    %% extract missed trials
    
    noAnswerTrials          = goodAnswers <= n_min_solved_pairs;
    noAnswerErrorTrials     = goodAnswers_ErrorTrials <= n_min_solved_pairs;
    noAnswerNoErrorTrials   = goodAnswers_noErrorTrials <= n_min_solved_pairs;
    % to gain
    goodAnswers_toGain = noAnswerTrials(trialValence == 1);
    Stroop_excludedTrials_toGain_idx = (goodAnswers_toGain <= n_min_solved_pairs); % if no pair correctly solved at all => exclude the trial (even if a bad answer was provided)
    if ~isempty(Stroop_excludedTrials_toGain_idx)
        for iOnset = 1:n_onsets % loop through onsets
            curr_onset_nm = onset_types{iOnset};
            curr_onset_missed_nm = [curr_onset_nm,'_missedTrial'];
            stroop.onset.toGain.(curr_onset_missed_nm)      = stroop_onset_toGain.(curr_onset_nm)( Stroop_excludedTrials_toGain_idx );
            stroop.duration.toGain.(curr_onset_missed_nm)   = stroop_dur_toGain.(curr_onset_nm)( Stroop_excludedTrials_toGain_idx );
        end % onsets loop
    end
    % error trials to gain
    goodAnswers_errorTrials_toGain = noAnswerErrorTrials(toGain_ErrorTrials);
    Stroop_excludedTrials_errorTrials_toGain_idx = (goodAnswers_errorTrials_toGain <= n_min_solved_pairs); % if no pair correctly solved at all => exclude the trial (even if a bad answer was provided)
    if ~isempty(Stroop_excludedTrials_errorTrials_toGain_idx)
        for iOnset = 1:n_onsets % loop through onsets
            curr_onset_nm = onset_types{iOnset};
            curr_onset_errorT_missed_nm = [curr_onset_nm,'_errorTrials','_missedTrial'];
            stroop.onset.toGain_errorTrials.(curr_onset_errorT_missed_nm)       = stroop.onset.toGain_ErrorTrials.(curr_onset_nm)( Stroop_excludedTrials_errorTrials_toGain_idx );
            stroop.duration.toGain_errorTrials.(curr_onset_errorT_missed_nm)    = stroop.duration.toGain_ErrorTrials.(curr_onset_nm)( Stroop_excludedTrials_errorTrials_toGain_idx );
        end % onsets loop
    end
    % no error trials to gain
    goodAnswers_noErrorTrials_toGain = noAnswerNoErrorTrials(toGain_noErrorTrials);
    Stroop_excludedTrials_noErrorTrials_toGain_idx = (goodAnswers_noErrorTrials_toGain <= n_min_solved_pairs); % if no pair correctly solved at all => exclude the trial (even if a bad answer was provided)
    if ~isempty(Stroop_excludedTrials_noErrorTrials_toGain_idx)
        for iOnset = 1:n_onsets % loop through onsets
            curr_onset_nm = onset_types{iOnset};
            curr_onset_noErrorT_missed_nm = [curr_onset_nm,'_noErrorTrials','_missedTrial'];
            stroop.onset.toGain_errorTrials.(curr_onset_noErrorT_missed_nm)     = stroop.onset.toGain_noErrorTrials.(curr_onset_nm)( Stroop_excludedTrials_noErrorTrials_toGain_idx );
            stroop.duration.toGain_errorTrials.(curr_onset_noErrorT_missed_nm)  = stroop.duration.toGain_noErrorTrials.(curr_onset_nm)( Stroop_excludedTrials_noErrorTrials_toGain_idx );
        end % onsets loop
    end
    
    % to lose
    goodAnswers_toLose = noAnswerTrials(trialValence == -1);
    Stroop_excludedTrials_toLose_idx = (goodAnswers_toLose <= n_min_solved_pairs); % if no pair correctly solved at all => exclude the trial (even if a bad answer was provided)
    if ~isempty(Stroop_excludedTrials_toLose_idx)
        for iOnset = 1:n_onsets % loop through onsets
            curr_onset_nm = onset_types{iOnset};
            curr_onset_missed_nm = [onset_types{iOnset},'_missedTrial'];
            stroop.onset.toLose_ErrorTrials.(curr_onset_missed_nm)      = stroop_onset_toLose.(curr_onset_nm)( Stroop_excludedTrials_toLose_idx );
            stroop.duration.toLose_ErrorTrials.(curr_onset_missed_nm)   = stroop_dur_toLose.(curr_onset_nm)( Stroop_excludedTrials_toLose_idx );
        end % onsets loop
    end
    % error trials to lose
    goodAnswers_errorTrials_toLose = noAnswerErrorTrials(toLose_ErrorTrials);
    Stroop_excludedTrials_errorTrials_toLose_idx = (goodAnswers_errorTrials_toLose <= n_min_solved_pairs); % if no pair correctly solved at all => exclude the trial (even if a bad answer was provided)
    if ~isempty(Stroop_excludedTrials_errorTrials_toLose_idx)
        for iOnset = 1:n_onsets % loop through onsets
            curr_onset_nm = onset_types{iOnset};
            curr_onset_errorT_missed_nm = [curr_onset_nm,'_errorTrials','_missedTrial'];
            stroop.onset.toLose_ErrorTrials.(curr_onset_errorT_missed_nm)       = stroop.onset.toLose_ErrorTrials.(curr_onset_nm)( Stroop_excludedTrials_errorTrials_toLose_idx );
            stroop.duration.toLose_ErrorTrials.(curr_onset_errorT_missed_nm)    = stroop.duration.toLose_ErrorTrials.(curr_onset_nm)( Stroop_excludedTrials_errorTrials_toLose_idx );
        end % onsets loop
    end
    % no error trials to lose
    goodAnswers_noErrorTrials_toLose = noAnswerNoErrorTrials(toLose_noErrorTrials);
    Stroop_excludedTrials_noErrorTrials_toLose_idx = (goodAnswers_noErrorTrials_toLose <= n_min_solved_pairs); % if no pair correctly solved at all => exclude the trial (even if a bad answer was provided)
    if ~isempty(Stroop_excludedTrials_noErrorTrials_toLose_idx)
        for iOnset = 1:n_onsets % loop through onsets
            curr_onset_nm = onset_types{iOnset};
            curr_onset_noErrorT_missed_nm = [curr_onset_nm,'_noErrorTrials','_missedTrial'];
            stroop.onset.toLose_ErrorTrials.(curr_onset_noErrorT_missed_nm)     = stroop.onset.toLose_noErrorTrials.(curr_onset_nm)( Stroop_excludedTrials_noErrorTrials_toLose_idx );
            stroop.duration.toLose_ErrorTrials.(curr_onset_noErrorT_missed_nm)  = stroop.duration.toLose_noErrorTrials.(curr_onset_nm)( Stroop_excludedTrials_noErrorTrials_toLose_idx );
        end % onsets loop
    end
    
    %% quick check that no onset or modulator has been forgotten in the initial list
    % (mostly relevant to be sure that when you extract a new variable to
    % extract inside the script, you do not forget to also add it to the
    % initial list of variables)
    
    % onsets
    stroop_onset_fields = fieldnames(stroop.onset.all);
    n_total_possible_onsets = n_onsets*2; % if you include missed trials
    if size(stroop_onset_fields,1) > n_total_possible_onsets
        error('Number of fields in onsets does not match with onset_types. You must have forgotten some. Please fix it');
    end
    
    % modulators
    stroop_mod_fields = fieldnames(stroop.mod.all);
    n_total_possible_mods = n_mod*2; % if you include missed trials
    if size(stroop_mod_fields,1) > n_total_possible_mods
        error('Number of fields in modulators does not match with mod_types. You must have forgotten some. Please fix it');
    end
    
    %% pool all data which is relevant for the GLM together in a new file
    save([fMRI_Dir 'onsets_sub',subid,'_stroop_run',sessnber,'.mat'],'stroop');
    
    %% clear stroop structure before next run to avoid interferences
    clear('stroop');
end % run

end