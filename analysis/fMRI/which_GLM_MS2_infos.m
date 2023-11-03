function [ ] = which_GLM_MS2_infos( GLMprm )
%which_GLM_MS2_infos( GLMprm ) outputs the informations about the GLM based
% on the list of parameters designated inside GLMprm.
%
% INPUTS
% GLMprm: structure with GLM parameters for each task extracted with
% which_GLM_MS2.
%
% See also which_GLM_MS2

gal     = GLMprm.gal;
RL      = GLMprm.RL;
grip    = GLMprm.grip;
stroop  = GLMprm.stroop;

%% main parameters info
disp('***main parameters');

switch gal.grey_mask
    case 1
        disp('Use of 1st level probability grey mask for each subject.');
end

switch gal.add_drv
    case 1
        disp('Use of temporal derivative => be aware that each regressor will be doubled in the matrix list.');
    case 2
        disp('Use of temporal and spatial derivative => be aware that each regressor will be tripled in the matrix list');
end

%% reinforcement-learning parameters info
disp('***RL parameters');

% stimulus display
o_stim = RL.o_stim;
dur_stim = RL.dur_stim;
mod_stim = RL.mod_stim;
if o_stim ~= 0
    switch o_stim
        case 1 % all pairs grouped
            fprintf('**Stimulus display (all pairs): ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset to end of trial');
                case 2
                    fprintf('boxcar from onset to RT');
            end
            fprintf('\n');
            
            switch mod_stim.trial_nber
                case 1
                    fprintf('trial number');
                case 2
                    fprintf('zPerRun trial number');
                case 3
                    fprintf('zAcrossRuns trial number');
            end
            fprintf('/');
            switch mod_stim.dQ
                case 1
                    fprintf('Q best cue - Q non-best cue');
                case 2
                    fprintf('Q best cue - Q non-best cue (zPerRun)');
                case 3
                    fprintf('Q best cue - Q non-best cue (zAcrossRuns)');
                case 4
                    fprintf('Q best cue - Q non-best cue (zPerRun/pair)');
                case 5
                    fprintf('Q best cue - Q non-best cue (zAcrossRuns per pair)');
            end
            fprintf('/');
            if mod_stim.RT ~= 0
                switch mod_stim.RT
                    case {1,2,3,4}
                        RTcorr_nm = '';
                    case {5,6,7,8}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {1,5}
                        fprintf(['RT',RTcorr_nm]);
                    case {2,6}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {3,7}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {4,8}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                end
                fprintf('/');
            end
            fprintf('\n');
            
        case 2 % split gain/neutral/loss pairs
            % gain pairs
            fprintf('**Stimulus display (gain): ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset to end of trial');
                case 2
                    fprintf('boxcar from onset to RT');
            end
            fprintf('\n');
            switch mod_stim.trial_nber
                case 1
                    fprintf('trial number');
                case 2
                    fprintf('zPerRun trial number');
                case 3
                    fprintf('zAcrossRuns trial number');
            end
            fprintf('/');
            switch mod_stim.dQ
                case 1
                    fprintf('Q gain cue - Q neutral cue');
                case 2
                    fprintf('Q gain cue - Q neutral cue (zPerRun)');
                case 3
                    fprintf('Q gain cue - Q neutral cue (zAcrossRuns)');
                case 4
                    fprintf('Q gain cue - Q neutral cue (zPerRun/pair)');
                case 5
                    fprintf('Q gain cue - Q neutral cue (zAcrossRuns per pair)');
            end
            fprintf('/');
            if mod_stim.RT ~= 0
                switch mod_stim.RT
                    case {1,2,3,4}
                        RTcorr_nm = '';
                    case {5,6,7,8}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {1,5}
                        fprintf(['RT',RTcorr_nm]);
                    case {2,6}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {3,7}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {4,8}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                end
                fprintf('/');
            end
            fprintf('\n');
            
            
            % neutral
            fprintf('**Stimulus display (neutral): ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset to end of trial');
                case 2
                    fprintf('boxcar from onset to RT');
            end
            fprintf('\n');
            switch mod_stim.trial_nber
                case 1
                    fprintf('trial number');
                case 2
                    fprintf('zPerRun trial number');
                case 3
                    fprintf('zAcrossRuns trial number');
            end
            fprintf('/');
            % no dQ for neutral pair
            if mod_stim.RT ~= 0
                switch mod_stim.RT
                    case {1,2,3,4}
                        RTcorr_nm = '';
                    case {5,6,7,8}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {1,5}
                        fprintf(['RT',RTcorr_nm]);
                    case {2,6}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {3,7}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {4,8}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                end
                fprintf('/');
            end
            fprintf('\n');
            
            
            % loss pairs
            fprintf('**Stimulus display (loss): ');
            switch dur_stim
                case 0
                    fprintf('stick function \n');
                case 1
                    fprintf('boxcar from onset to end of trial \n');
                case 2
                    fprintf('boxcar from onset to RT \n');
            end
            fprintf('\n');
            
            switch mod_stim.trial_nber
                case 1
                    fprintf('trial number');
                case 2
                    fprintf('zPerRun trial number');
                case 3
                    fprintf('zAcrossRuns trial number');
            end
            fprintf('/');
            switch mod_stim.dQ
                case 1
                    fprintf('Q neutral cue - Q loss cue');
                case 2
                    fprintf('Q neutral cue - Q loss cue (zPerRun)');
                case 3
                    fprintf('Q neutral cue - Q loss cue (zAcrossRuns)');
                case 4
                    fprintf('Q neutral cue - Q loss cue (zPerRun/pair)');
                case 5
                    fprintf('Q neutral cue - Q loss cue (zAcrossRuns per pair)');
            end
            fprintf('/');
            if mod_stim.RT ~= 0
                switch mod_stim.RT
                    case {1,2,3,4}
                        RTcorr_nm = '';
                    case {5,6,7,8}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {1,5}
                        fprintf(['RT',RTcorr_nm]);
                    case {2,6}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {3,7}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {4,8}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                end
                fprintf('/');
            end
            fprintf('\n');
    end
end

% RT
o_answer = RL.o_answer;
dur_answer = RL.dur_answer;
% mod_answer = RL.mod_answer;
if o_answer ~= 0
    switch o_answer
        case 1 % all pairs grouped
            fprintf('**Answer (all pairs): ');
            switch dur_answer
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from answer to end of trial');
            end
            fprintf('\n');
        case 2 % split gain/neutral/loss pair
            fprintf('**Answer (gain pair): ');
            switch dur_answer
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from answer to end of trial');
            end
            fprintf('\n');
            
            fprintf('**Answer (neutral pair): ');
            switch dur_answer
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from answer to end of trial');
            end
            fprintf('/');
            fprintf('\n');
            
            fprintf('**Answer (loss pair): ');
            switch dur_answer
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from answer to end of trial');
            end
            fprintf('/');
            fprintf('\n');
    end
end


% feedback
o_fbk = RL.o_fbk;
dur_fbk = RL.dur_fbk;
mod_fbk = RL.mod_fbk;
if o_fbk ~= 0
    switch o_fbk
        case 1 % group all feedbacks
            fprintf('**Feedback (all pairs): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
            end
            fprintf('\n');
            
            if mod_fbk.trial_nber ~= 0
                switch mod_stim.trial_nber
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            if mod_fbk.fbk ~= 0
                switch mod_fbk.fbk
                    case 1
                        fprintf('feedback (current trial)');
                    case 2
                        fprintf('Prediction Error');
                end
                fprintf('/');
            end
            
            if mod_fbk.totalGain ~= 0
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
                fprintf('/');
            end
            
        case 2 % slit by pair type gain pair/neutral pair/loss pair
            fprintf('**Feedback (gain pair): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
            end
            fprintf('\n');
            
            fprintf('**Feedback (neutral pair): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
            end
            fprintf('\n');
            
            fprintf('**Feedback (loss pair): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
            end
            fprintf('\n');
            
        case 3 % split by feedback type gain feedback/neutral feedback/loss feedback
            fprintf('**Feedback (gain feedback): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
            end
            fprintf('\n');
            
            fprintf('**Feedback (neutral feedback): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
            end
            fprintf('\n');
            
            fprintf('**Feedback (loss feedback): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
            end
            fprintf('\n');
            
        case 4 % split by gain and feedback type gain pair gain feedback/gain pair neutral feedback/neutral pair & feedback/loss pair neutral feedback/loss pair loss feedback
            fprintf('**Feedback (gain pair - gain feedback): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
            end
            fprintf('\n');
            
            fprintf('**Feedback (gain pair - neutral feedback): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
            end
            fprintf('\n');
            
            fprintf('**Feedback (neutral pair & feedback): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
            end
            fprintf('\n');
            
            fprintf('**Feedback (loss pair - neutral feedback): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
            end
            fprintf('\n');
            
            fprintf('**Feedback (loss pair - loss feedback): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
            end
            fprintf('\n');
    end
end


% cross
o_cross = RL.o_cross;
dur_cross = RL.dur_cross;
if o_cross ~= 0
    switch o_cross
        case 1
            fprintf('**Cross (all pairs): ');
            switch dur_cross
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
            end
            fprintf('\n');
    end
end


% missed trials
o_missed_trials_stim = RL.o_missed_trials_stim;
dur_missed_trials_stim = RL.dur_missed_trials_stim;
switch o_missed_trials_stim
    case 1
        fprintf('**Missed Trials (all pairs): ');
        switch dur_missed_trials_stim
            case 0
                fprintf('stick function');
            case 1
                fprintf('boxcar from display of stimuli to end of trial');
        end
end

%% grip parameters info
disp('***Grip parameters');

%% stroop parameters info
disp('***Stroop parameters');

end