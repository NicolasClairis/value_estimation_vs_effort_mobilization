function [ n_prm_RL ] = which_GLM_MS2_infos_RL( GLMprm )
%which_GLM_MS2_infos_RL( GLMprm ) outputs the informations about the GLM based
% on the list of parameters designated inside GLMprm for the
% reinforcement-learning task.
%
% INPUTS
% GLMprm: structure with GLM parameters for each task extracted with
% which_GLM_MS2.
%
% OUTPUTS
% n_prm_RL: number of parameters in RL task (includes onsets and
% regressors) for a given run (without taking into account derivative)
%
% See also which_GLM_MS2

%% reinforcement-learning parameters info
RL      = GLMprm.RL;
disp('***RL parameters');
n_prm_RL = 0;

%% Chosen stimulus display
o_stim = RL.o_stim;
dur_stim = RL.dur_stim;
mod_stim = RL.mod_stim;
if o_stim ~= 0
    switch o_stim
        case 1 % all pairs grouped
            
            %% all pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Stimulus display (all pairs): ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until choice in red');
                case 2
                    fprintf('boxcar from onset until RT');
                case 3
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            %% RT
            if ismember(mod_stim.RT,9:16)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.RT
                    case {9,10,11,12}
                        RTcorr_nm = '';
                    case {13,14,15,16}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {9,13}
                        fprintf(['RT',RTcorr_nm]);
                    case {10,14}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {11,15}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {12,16}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number
            if ismember(mod_stim.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            
            %% model used for SV, dQ and p(choice=best)
            if ~isempty(mod_stim.mdl_type) && mod_stim.mdl_n ~= 0
                fprintf([mod_stim.mdl_type,' model ',num2str(mod_stim.mdl_n),' ']);
            end
            %% SV
            if ismember(mod_stim.SV,1:5)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 1
                        fprintf('pA*QA+pB*QB');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% dQ
            if mod_stim.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.dQ
                    case {1}
                        fprintf('Q best cue - Q non-best cue');
                    case {2}
                        fprintf('Q best cue - Q non-best cue (zPerRun)');
                    case {3}
                        fprintf('Q best cue - Q non-best cue (zAcrossRuns)');
                    case {4}
                        fprintf('Q best cue - Q non-best cue (zPerRun/pair)');
                    case {5}
                        fprintf('Q best cue - Q non-best cue (zAcrossRuns per pair)');
                    case {6}
                        fprintf('(Q best cue - Q non-best cue)^2');
                    case {7}
                        fprintf('|Q best cue - Q non-best cue|/sqrt((Qbest+0.5)^2+(QnonBest+0.5)^2)');
                    case {8}
                        fprintf('|Q best cue - Q non-best cue|/sqrt((Qbest+0.5)^2+(QnonBest+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice=best option)
            if mod_stim.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.pBest
                    case 1
                        fprintf('[p(choice=best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number
            if ismember(mod_stim.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor
            if mod_stim.ROI_activity_yn ~= 0
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_stim.ROI_activity_GLM,...
                    ' ',mod_stim.ROI_activity_period,' period ',...
                    mod_stim.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT
            if ismember(mod_stim.RT,1:8)
                n_prm_RL = n_prm_RL + 1;
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
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% SV
            if ismember(mod_stim.SV,6:10)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 6
                        fprintf('pA*QA+pB*QB');
                    case 7
                        fprintf('QA+QB');
                    case 8
                        fprintf('Qch-Qunch');
                    case 9
                        fprintf('Qch');
                    case 10
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number
            if ismember(mod_stim.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
        case 2 % split gain/neutral/loss pairs
            %% gain pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Stimulus display (gain): ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until choice in red');
                case 2
                    fprintf('boxcar from onset until RT');
                case 3
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% RT (gain pair)
            if ismember(mod_stim.RT,9:16)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.RT
                    case {9,10,11,12}
                        RTcorr_nm = '';
                    case {13,14,15,16}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {9,13}
                        fprintf(['RT',RTcorr_nm]);
                    case {10,14}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {11,15}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {12,16}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain pair)
            if ismember(mod_stim.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% model used for SV, dQ and p(choice=best)
            if ~isempty(mod_stim.mdl_type) && mod_stim.mdl_n ~= 0
                fprintf([mod_stim.mdl_type,' model ',num2str(mod_stim.mdl_n),' ']);
            end
            %% SV (gain pair)
            if ismember(mod_stim.SV,1:5)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 1
                        fprintf('pGainItem*QgainItem+pNeutralItem*QneutralItem');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% dQ (gain pair)
            if mod_stim.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
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
                    case 6
                        fprintf('(Q gain cue - Q neutral cue)^2');
                    case 7
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2)');
                    case 8
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice=best option) (gain pair)
            if mod_stim.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.pBest
                    case 1
                        fprintf('[p(choice=best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain pair)
            if ismember(mod_stim.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (gain pair)
            if mod_stim.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_stim.ROI_activity_GLM,...
                    ' ',mod_stim.ROI_activity_period,' period ',...
                    mod_stim.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (gain pair)
            if ismember(mod_stim.RT,1:8)
                n_prm_RL = n_prm_RL + 1;
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
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% SV (gain pair)
            if ismember(mod_stim.SV,6:10)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 6
                        fprintf('pGainItem*QgainItem+pNeutralItem*QneutralItem');
                    case 7
                        fprintf('QA+QB');
                    case 8
                        fprintf('Qch-Qunch');
                    case 9
                        fprintf('Qch');
                    case 10
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain pair)
            if ismember(mod_stim.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            
            %% neutral pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Stimulus display (neutral): ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until choice in red');
                case 2
                    fprintf('boxcar from onset until RT');
                case 3
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% RT (neutral pair)
            if ismember(mod_stim.RT,9:16)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.RT
                    case {9,10,11,12}
                        RTcorr_nm = '';
                    case {13,14,15,16}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {9,13}
                        fprintf(['RT',RTcorr_nm]);
                    case {10,14}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {11,15}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {12,16}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (neutral pair)
            if ismember(mod_stim.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% no SV for neutral pair
            %% no dQ for neutral pair
            
            %% trial number (neutral pair)
            if ismember(mod_stim.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (neutral pair)
            if mod_stim.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_stim.ROI_activity_GLM,...
                    ' ',mod_stim.ROI_activity_period,' period ',...
                    mod_stim.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (neutral pair)
            if ismember(mod_stim.RT,1:8)
                n_prm_RL = n_prm_RL + 1;
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
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (neutral pair)
            if ismember(mod_stim.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            
            %% loss pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Stimulus display (loss): ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until choice in red');
                case 2
                    fprintf('boxcar from onset until RT');
                case 3
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% RT (loss pair)
            if ismember(mod_stim.RT,9:16)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.RT
                    case {9,10,11,12}
                        RTcorr_nm = '';
                    case {13,14,15,16}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {9,13}
                        fprintf(['RT',RTcorr_nm]);
                    case {10,14}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {11,15}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {12,16}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (loss pair)
            if ismember(mod_stim.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% model used for SV, dQ and p(choice=best)
            if ~isempty(mod_stim.mdl_type) && mod_stim.mdl_n ~= 0
                fprintf([mod_stim.mdl_type,' model ',num2str(mod_stim.mdl_n),' ']);
            end
            %% SV (loss pair)
            if ismember(mod_stim.SV,1:5)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 1
                        fprintf('pLossItem*QlossItem+pNeutralItem*QneutralItem');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% dQ (loss pair)
            if mod_stim.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
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
                    case 6
                        fprintf('(Q neutral cue - Q loss cue)^2');
                    case 7
                        fprintf('|Q neutral cue - Q loss cue|/sqrt((QneutralCue+0.5)^2+(QlossCue+0.5)^2)');
                    case 8
                        fprintf('|Q neutral cue - Q loss cue|/sqrt((QneutralCue+0.5)^2+(QlossCue+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice = best option)
            if mod_stim.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.pBest
                    case 1
                        fprintf('[p(choice = best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (loss pair)
            if ismember(mod_stim.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (loss pair)
            if mod_stim.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_stim.ROI_activity_GLM,...
                    ' ',mod_stim.ROI_activity_period,' period ',...
                    mod_stim.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (loss pair)
            if ismember(mod_stim.RT,1:8)
                n_prm_RL = n_prm_RL + 1;
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
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% SV (loss pair)
            if ismember(mod_stim.SV,6:10)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 6
                        fprintf('pLossItem*QlossItem+pNeutralItem*QneutralItem');
                    case 7
                        fprintf('QA+QB');
                    case 8
                        fprintf('Qch-Qunch');
                    case 9
                        fprintf('Qch');
                    case 10
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (loss pair)
            if ismember(mod_stim.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
        case 3 % pool gain and loss together/neutral pair apart
            %% gain + loss pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Stimulus display (gain + loss): ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until choice in red');
                case 2
                    fprintf('boxcar from onset until RT');
                case 3
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% RT (gain+loss pair)
            if ismember(mod_stim.RT,9:16)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.RT
                    case {9,10,11,12}
                        RTcorr_nm = '';
                    case {13,14,15,16}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {9,13}
                        fprintf(['RT',RTcorr_nm]);
                    case {10,14}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {11,15}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {12,16}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain + loss pair)
            if ismember(mod_stim.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% model used for SV, dQ and p(choice=best)
            if ~isempty(mod_stim.mdl_type) && mod_stim.mdl_n ~= 0
                fprintf([mod_stim.mdl_type,' model ',num2str(mod_stim.mdl_n),' ']);
            end
            %% SV (gain + loss pair)
            if ismember(mod_stim.SV,1:5)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 1
                        fprintf('pBestItem*QbestItem+pWorseItem*QworseItem');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% dQ (gain + loss pair)
            if mod_stim.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
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
                    case 6
                        fprintf('(Q gain cue - Q neutral cue)^2');
                    case 7
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2)');
                    case 8
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice=best option) (gain + loss pair)
            if mod_stim.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.pBest
                    case 1
                        fprintf('[p(choice=best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain + loss pair)
            if ismember(mod_stim.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (gain + loss pair)
            if mod_stim.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_stim.ROI_activity_GLM,...
                    ' ',mod_stim.ROI_activity_period,' period ',...
                    mod_stim.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (gain + loss pair)
            if ismember(mod_stim.RT,1:8)
                n_prm_RL = n_prm_RL + 1;
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
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% SV (gain + loss pair)
            if ismember(mod_stim.SV,6:10)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 6
                        fprintf('pBestItem*QbestItem+pWorseItem*QworseItem');
                    case 7
                        fprintf('QA+QB');
                    case 8
                        fprintf('Qch-Qunch');
                    case 9
                        fprintf('Qch');
                    case 10
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain + loss pair)
            if ismember(mod_stim.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            
            %% neutral pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Stimulus display (neutral): ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until choice in red');
                case 2
                    fprintf('boxcar from onset until RT');
                case 3
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% RT (neutral pair)
            if ismember(mod_stim.RT,9:16)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.RT
                    case {9,10,11,12}
                        RTcorr_nm = '';
                    case {13,14,15,16}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {9,13}
                        fprintf(['RT',RTcorr_nm]);
                    case {10,14}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {11,15}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {12,16}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (neutral pair)
            if ismember(mod_stim.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% no SV for neutral pair
            %% no dQ for neutral pair
            
            %% trial number (neutral pair)
            if ismember(mod_stim.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (neutral pair)
            if mod_stim.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_stim.ROI_activity_GLM,...
                    ' ',mod_stim.ROI_activity_period,' period ',...
                    mod_stim.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (neutral pair)
            if ismember(mod_stim.RT,1:8)
                n_prm_RL = n_prm_RL + 1;
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
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (neutral pair)
            if ismember(mod_stim.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
        case 4
            error('not ready yet');
        case 5
            %% gain pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Stimulus display (gain first trials): ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until choice in red');
                case 2
                    fprintf('boxcar from onset until RT');
                case 3
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% RT (gain pair first trials)
            if ismember(mod_stim.RT,9:16)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.RT
                    case {9,10,11,12}
                        RTcorr_nm = '';
                    case {13,14,15,16}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {9,13}
                        fprintf(['RT',RTcorr_nm]);
                    case {10,14}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {11,15}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {12,16}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain pair first trials)
            if ismember(mod_stim.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% model used for SV, dQ and p(choice=best)
            if ~isempty(mod_stim.mdl_type) && mod_stim.mdl_n ~= 0
                fprintf([mod_stim.mdl_type,' model ',num2str(mod_stim.mdl_n),' ']);
            end
            %% SV (gain pair first trials)
            if ismember(mod_stim.SV,1:5)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 1
                        fprintf('pGainItem*QgainItem+pNeutralItem*QneutralItem');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% dQ (gain pair first trials)
            if mod_stim.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
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
                    case 6
                        fprintf('(Q gain cue - Q neutral cue)^2');
                    case 7
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2)');
                    case 8
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice=best option) (gain pair first trials)
            if mod_stim.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.pBest
                    case 1
                        fprintf('[p(choice=best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain pair first trials)
            if ismember(mod_stim.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (gain pair first trials)
            if mod_stim.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_stim.ROI_activity_GLM,...
                    ' ',mod_stim.ROI_activity_period,' period ',...
                    mod_stim.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (gain pair first trials)
            if ismember(mod_stim.RT,1:8)
                n_prm_RL = n_prm_RL + 1;
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
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% SV (gain pair first trials)
            if ismember(mod_stim.SV,6:10)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 6
                        fprintf('pGainItem*QgainItem+pNeutralItem*QneutralItem');
                    case 7
                        fprintf('QA+QB');
                    case 8
                        fprintf('Qch-Qunch');
                    case 9
                        fprintf('Qch');
                    case 10
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain pair first trials)
            if ismember(mod_stim.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            
            %% neutral pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Stimulus display (neutral first trials): ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until choice in red');
                case 2
                    fprintf('boxcar from onset until RT');
                case 3
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% RT (neutral pair first trials)
            if ismember(mod_stim.RT,9:16)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.RT
                    case {9,10,11,12}
                        RTcorr_nm = '';
                    case {13,14,15,16}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {9,13}
                        fprintf(['RT',RTcorr_nm]);
                    case {10,14}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {11,15}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {12,16}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (neutral pair first trials)
            if ismember(mod_stim.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% no SV for neutral pair
            %% no dQ for neutral pair
            
            %% trial number (neutral pair, first trials)
            if ismember(mod_stim.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (neutral pair first trials)
            if mod_stim.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_stim.ROI_activity_GLM,...
                    ' ',mod_stim.ROI_activity_period,' period ',...
                    mod_stim.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (neutral pair first trials)
            if ismember(mod_stim.RT,1:8)
                n_prm_RL = n_prm_RL + 1;
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
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (neutral pair first trials)
            if ismember(mod_stim.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            
            %% loss pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Stimulus display (loss first trials): ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until choice in red');
                case 2
                    fprintf('boxcar from onset until RT');
                case 3
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% RT (loss pair first trials)
            if ismember(mod_stim.RT,9:16)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.RT
                    case {9,10,11,12}
                        RTcorr_nm = '';
                    case {13,14,15,16}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {9,13}
                        fprintf(['RT',RTcorr_nm]);
                    case {10,14}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {11,15}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {12,16}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (loss pair first trials)
            if ismember(mod_stim.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% model used for SV, dQ and p(choice=best)
            if ~isempty(mod_stim.mdl_type) && mod_stim.mdl_n ~= 0
                fprintf([mod_stim.mdl_type,' model ',num2str(mod_stim.mdl_n),' ']);
            end
            %% SV (loss pair first trials)
            if ismember(mod_stim.SV,1:5)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 1
                        fprintf('pLossItem*QlossItem+pNeutralItem*QneutralItem');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% dQ (loss pair first trials)
            if mod_stim.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
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
                    case 6
                        fprintf('(Q neutral cue - Q loss cue)^2');
                    case 7
                        fprintf('|Q neutral cue - Q loss cue|/sqrt((QneutralCue+0.5)^2+(QlossCue+0.5)^2)');
                    case 8
                        fprintf('|Q neutral cue - Q loss cue|/sqrt((QneutralCue+0.5)^2+(QlossCue+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice = best option) (loss pair first trials)
            if mod_stim.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.pBest
                    case 1
                        fprintf('[p(choice = best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (loss pair first trials)
            if ismember(mod_stim.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (loss pair first trials)
            if mod_stim.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_stim.ROI_activity_GLM,...
                    ' ',mod_stim.ROI_activity_period,' period ',...
                    mod_stim.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (loss pair first trials)
            if ismember(mod_stim.RT,1:8)
                n_prm_RL = n_prm_RL + 1;
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
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% SV (loss pair first trials)
            if ismember(mod_stim.SV,6:10)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 6
                        fprintf('pLossItem*QlossItem+pNeutralItem*QneutralItem');
                    case 7
                        fprintf('QA+QB');
                    case 8
                        fprintf('Qch-Qunch');
                    case 9
                        fprintf('Qch');
                    case 10
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (loss pair first trials)
            if ismember(mod_stim.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            %% gain pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Stimulus display (gain last trials): ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until choice in red');
                case 2
                    fprintf('boxcar from onset until RT');
                case 3
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% RT (gain pair last trials)
            if ismember(mod_stim.RT,9:16)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.RT
                    case {9,10,11,12}
                        RTcorr_nm = '';
                    case {13,14,15,16}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {9,13}
                        fprintf(['RT',RTcorr_nm]);
                    case {10,14}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {11,15}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {12,16}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain pair last trials)
            if ismember(mod_stim.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% model used for SV, dQ and p(choice=best)
            if ~isempty(mod_stim.mdl_type) && mod_stim.mdl_n ~= 0
                fprintf([mod_stim.mdl_type,' model ',num2str(mod_stim.mdl_n),' ']);
            end
            %% SV (gain pair last trials)
            if ismember(mod_stim.SV,1:5)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 1
                        fprintf('pGainItem*QgainItem+pNeutralItem*QneutralItem');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% dQ (gain pair last trials)
            if mod_stim.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
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
                    case 6
                        fprintf('(Q gain cue - Q neutral cue)^2');
                    case 7
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2)');
                    case 8
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice=best option) (gain pair last trials)
            if mod_stim.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.pBest
                    case 1
                        fprintf('[p(choice=best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain pair last trials)
            if ismember(mod_stim.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (gain pair last trials)
            if mod_stim.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_stim.ROI_activity_GLM,...
                    ' ',mod_stim.ROI_activity_period,' period ',...
                    mod_stim.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (gain pair last trials)
            if ismember(mod_stim.RT,1:8)
                n_prm_RL = n_prm_RL + 1;
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
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% SV (gain pair last trials)
            if ismember(mod_stim.SV,6:10)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 6
                        fprintf('pGainItem*QgainItem+pNeutralItem*QneutralItem');
                    case 7
                        fprintf('QA+QB');
                    case 8
                        fprintf('Qch-Qunch');
                    case 9
                        fprintf('Qch');
                    case 10
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain pair last trials)
            if ismember(mod_stim.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            
            %% neutral pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Stimulus display (neutral last trials): ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until choice in red');
                case 2
                    fprintf('boxcar from onset until RT');
                case 3
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% RT (neutral pair last trials)
            if ismember(mod_stim.RT,9:16)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.RT
                    case {9,10,11,12}
                        RTcorr_nm = '';
                    case {13,14,15,16}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {9,13}
                        fprintf(['RT',RTcorr_nm]);
                    case {10,14}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {11,15}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {12,16}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (neutral pair last trials)
            if ismember(mod_stim.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% no SV for neutral pair
            %% no dQ for neutral pair
            
            %% trial number (neutral pair last trials)
            if ismember(mod_stim.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (neutral pair last trials)
            if mod_stim.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_stim.ROI_activity_GLM,...
                    ' ',mod_stim.ROI_activity_period,' period ',...
                    mod_stim.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (neutral pair last trials)
            if ismember(mod_stim.RT,1:8)
                n_prm_RL = n_prm_RL + 1;
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
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (neutral pair last trials)
            if ismember(mod_stim.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            
            %% loss pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Stimulus display (loss last trials): ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until choice in red');
                case 2
                    fprintf('boxcar from onset until RT');
                case 3
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% RT (loss pair last trials)
            if ismember(mod_stim.RT,9:16)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.RT
                    case {9,10,11,12}
                        RTcorr_nm = '';
                    case {13,14,15,16}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {9,13}
                        fprintf(['RT',RTcorr_nm]);
                    case {10,14}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {11,15}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {12,16}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (loss pair last trials)
            if ismember(mod_stim.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% model used for SV, dQ and p(choice=best)
            if ~isempty(mod_stim.mdl_type) && mod_stim.mdl_n ~= 0
                fprintf([mod_stim.mdl_type,' model ',num2str(mod_stim.mdl_n),' ']);
            end
            %% SV (loss pair last trials)
            if ismember(mod_stim.SV,1:5)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 1
                        fprintf('pLossItem*QlossItem+pNeutralItem*QneutralItem');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% dQ (loss pair last trials)
            if mod_stim.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
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
                    case 6
                        fprintf('(Q neutral cue - Q loss cue)^2');
                    case 7
                        fprintf('|Q neutral cue - Q loss cue|/sqrt((QneutralCue+0.5)^2+(QlossCue+0.5)^2)');
                    case 8
                        fprintf('|Q neutral cue - Q loss cue|/sqrt((QneutralCue+0.5)^2+(QlossCue+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice = best option) (loss pair last trials)
            if mod_stim.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.pBest
                    case 1
                        fprintf('[p(choice = best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (loss pair last trials)
            if ismember(mod_stim.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (loss pair last trials)
            if mod_stim.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_stim.ROI_activity_GLM,...
                    ' ',mod_stim.ROI_activity_period,' period ',...
                    mod_stim.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (loss pair last trials)
            if ismember(mod_stim.RT,1:8)
                n_prm_RL = n_prm_RL + 1;
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
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% SV (loss pair last trials)
            if ismember(mod_stim.SV,6:10)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 6
                        fprintf('pLossItem*QlossItem+pNeutralItem*QneutralItem');
                    case 7
                        fprintf('QA+QB');
                    case 8
                        fprintf('Qch-Qunch');
                    case 9
                        fprintf('Qch');
                    case 10
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (loss pair last trials)
            if ismember(mod_stim.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
        case 6 % pool gain+loss/ neutral apart + split first/second half
            %% gain + loss pairs onset + duration - first trials
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Stimulus display (gain + loss) - first trials: ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until choice in red');
                case 2
                    fprintf('boxcar from onset until RT');
                case 3
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% RT (gain+loss pair) - first trials
            if ismember(mod_stim.RT,9:16)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.RT
                    case {9,10,11,12}
                        RTcorr_nm = '';
                    case {13,14,15,16}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {9,13}
                        fprintf(['RT',RTcorr_nm]);
                    case {10,14}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {11,15}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {12,16}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain + loss pair) - first trials
            if ismember(mod_stim.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% model used for SV, dQ and p(choice=best) - first trials
            if ~isempty(mod_stim.mdl_type) && mod_stim.mdl_n ~= 0
                fprintf([mod_stim.mdl_type,' model ',num2str(mod_stim.mdl_n),' ']);
            end
            %% SV (gain + loss pair) - first trials
            if ismember(mod_stim.SV,1:5)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 1
                        fprintf('pBestItem*QbestItem+pWorseItem*QworseItem');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% dQ (gain + loss pair) - first trials
            if mod_stim.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
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
                    case 6
                        fprintf('(Q gain cue - Q neutral cue)^2');
                    case 7
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2)');
                    case 8
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice=best option) (gain + loss pair) - first trials
            if mod_stim.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.pBest
                    case 1
                        fprintf('[p(choice=best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain + loss pair) - first trials
            if ismember(mod_stim.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (gain + loss pair) - first trials
            if mod_stim.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_stim.ROI_activity_GLM,...
                    ' ',mod_stim.ROI_activity_period,' period ',...
                    mod_stim.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (gain + loss pair) - first trials
            if ismember(mod_stim.RT,1:8)
                n_prm_RL = n_prm_RL + 1;
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
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% SV (gain + loss pair) - first trials
            if ismember(mod_stim.SV,6:10)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 6
                        fprintf('pBestItem*QbestItem+pWorseItem*QworseItem');
                    case 7
                        fprintf('QA+QB');
                    case 8
                        fprintf('Qch-Qunch');
                    case 9
                        fprintf('Qch');
                    case 10
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain + loss pair) - first trials
            if ismember(mod_stim.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            
            %% neutral pairs onset + duration - first trials
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Stimulus display (neutral) - first trials: ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until choice in red');
                case 2
                    fprintf('boxcar from onset until RT');
                case 3
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% RT (neutral pair) - first trials
            if ismember(mod_stim.RT,9:16)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.RT
                    case {9,10,11,12}
                        RTcorr_nm = '';
                    case {13,14,15,16}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {9,13}
                        fprintf(['RT',RTcorr_nm]);
                    case {10,14}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {11,15}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {12,16}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (neutral pair) - first trials
            if ismember(mod_stim.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% no SV for neutral pair
            %% no dQ for neutral pair
            
            %% trial number (neutral pair) - first trials
            if ismember(mod_stim.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (neutral pair) - first trials
            if mod_stim.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_stim.ROI_activity_GLM,...
                    ' ',mod_stim.ROI_activity_period,' period ',...
                    mod_stim.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (neutral pair) - first trials
            if ismember(mod_stim.RT,1:8)
                n_prm_RL = n_prm_RL + 1;
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
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (neutral pair) - first trials
            if ismember(mod_stim.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            %% gain + loss pairs onset + duration - last trials
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Stimulus display (gain + loss) - last trials: ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until choice in red');
                case 2
                    fprintf('boxcar from onset until RT');
                case 3
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% RT (gain+loss pair) - last trials
            if ismember(mod_stim.RT,9:16)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.RT
                    case {9,10,11,12}
                        RTcorr_nm = '';
                    case {13,14,15,16}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {9,13}
                        fprintf(['RT',RTcorr_nm]);
                    case {10,14}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {11,15}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {12,16}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain + loss pair) - last trials
            if ismember(mod_stim.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% model used for SV, dQ and p(choice=best) - last trials
            if ~isempty(mod_stim.mdl_type) && mod_stim.mdl_n ~= 0
                fprintf([mod_stim.mdl_type,' model ',num2str(mod_stim.mdl_n),' ']);
            end
            %% SV (gain + loss pair) - last trials
            if ismember(mod_stim.SV,1:5)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 1
                        fprintf('pBestItem*QbestItem+pWorseItem*QworseItem');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% dQ (gain + loss pair) - last trials
            if mod_stim.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
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
                    case 6
                        fprintf('(Q gain cue - Q neutral cue)^2');
                    case 7
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2)');
                    case 8
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice=best option) (gain + loss pair) - last trials
            if mod_stim.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.pBest
                    case 1
                        fprintf('[p(choice=best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain + loss pair) - last trials
            if ismember(mod_stim.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (gain + loss pair) - last trials
            if mod_stim.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_stim.ROI_activity_GLM,...
                    ' ',mod_stim.ROI_activity_period,' period ',...
                    mod_stim.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (gain + loss pair) - last trials
            if ismember(mod_stim.RT,1:8)
                n_prm_RL = n_prm_RL + 1;
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
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% SV (gain + loss pair) - last trials
            if ismember(mod_stim.SV,6:10)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.SV
                    case 6
                        fprintf('pBestItem*QbestItem+pWorseItem*QworseItem');
                    case 7
                        fprintf('QA+QB');
                    case 8
                        fprintf('Qch-Qunch');
                    case 9
                        fprintf('Qch');
                    case 10
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain + loss pair) - last trials
            if ismember(mod_stim.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            
            %% neutral pairs onset + duration - last trials
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Stimulus display (neutral) - last trials: ');
            switch dur_stim
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until choice in red');
                case 2
                    fprintf('boxcar from onset until RT');
                case 3
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% RT (neutral pair) - last trials
            if ismember(mod_stim.RT,9:16)
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.RT
                    case {9,10,11,12}
                        RTcorr_nm = '';
                    case {13,14,15,16}
                        RTcorr_nm = ' corrected';
                end
                switch mod_stim.RT
                    case {9,13}
                        fprintf(['RT',RTcorr_nm]);
                    case {10,14}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {11,15}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {12,16}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (neutral pair) - last trials
            if ismember(mod_stim.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% no SV for neutral pair
            %% no dQ for neutral pair
            
            %% trial number (neutral pair) - last trials
            if ismember(mod_stim.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (neutral pair) - last trials
            if mod_stim.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_stim.ROI_activity_GLM,...
                    ' ',mod_stim.ROI_activity_period,' period ',...
                    mod_stim.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (neutral pair) - last trials
            if ismember(mod_stim.RT,1:8)
                n_prm_RL = n_prm_RL + 1;
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
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (neutral pair) - last trials
            if ismember(mod_stim.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_stim.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
    end
end

%% onset at the moment when subject answers (RT) + duration
o_answer = RL.o_answer;
dur_answer = RL.dur_answer;
% mod_answer = RL.mod_answer;
if o_answer ~= 0
    switch o_answer
        case 1 % all pairs grouped
            %% RT onset + duration (all pairs)
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Answer (all pairs): ');
            switch dur_answer
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from RT until choice in red');
                case 2
                    fprintf('boxcar from RT until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
        case 2 % split gain/neutral/loss pair
            %% gain pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Answer (gain pair): ');
            switch dur_answer
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from RT until choice in red');
                case 2
                    fprintf('boxcar from RT until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% neutral pair onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Answer (neutral pair): ');
            switch dur_answer
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from RT until choice in red');
                case 2
                    fprintf('boxcar from RT until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% loss pair onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Answer (loss pair): ');
            switch dur_answer
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from RT until choice in red');
                case 2
                    fprintf('boxcar from RT until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
    end
end

%% chosen option in red
o_chosen    = RL.o_chosen;
dur_chosen  = RL.dur_chosen;
mod_chosen  = RL.mod_chosen;
if o_chosen ~= 0
    switch o_chosen
        case 1 % all pairs grouped
            
            %% all pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Chosen stimulus display (all pairs): ');
            switch dur_chosen
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from choice in red until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% trial number
            if ismember(mod_chosen.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% model used for SV, dQ and p(choice=best)
            if ~isempty(mod_chosen.mdl_type) && mod_chosen.mdl_n ~= 0
                fprintf([mod_chosen.mdl_type,' model ',num2str(mod_chosen.mdl_n),' ']);
            end
            %% SV
            if mod_chosen.SV ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.SV
                    case 1
                        fprintf('pA*QA+pB*QB');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% dQ
            if mod_chosen.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.dQ
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
                    case 6
                        fprintf('(Q best cue - Q non-best cue)^2');
                    case 7
                        fprintf('|Q best cue - Q non-best cue|/sqrt((Qbest+0.5)^2+(QnonBest+0.5)^2)');
                    case 8
                        fprintf('|Q best cue - Q non-best cue|/sqrt((Qbest+0.5)^2+(QnonBest+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice=best option)
            if mod_chosen.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.pBest
                    case 1
                        fprintf('[p(choice=best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number
            if ismember(mod_chosen.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor
            if mod_chosen.ROI_activity_yn ~= 0
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_chosen.ROI_activity_GLM,...
                    ' ',mod_chosen.ROI_activity_period,' period ',...
                    mod_chosen.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT
            if mod_chosen.RT ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.RT
                    case {1,2,3,4}
                        RTcorr_nm = '';
                    case {5,6,7,8}
                        RTcorr_nm = ' corrected';
                end
                switch mod_chosen.RT
                    case {1,5}
                        fprintf(['RT',RTcorr_nm]);
                    case {2,6}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {3,7}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {4,8}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number
            if ismember(mod_chosen.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
        case 2 % split gain/neutral/loss pairs
            %% gain pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Chosen stimulus display (gain): ');
            switch dur_chosen
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from choice in red until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% trial number (gain pair)
            if ismember(mod_chosen.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% model used for SV, dQ and p(choice=best)
            if ~isempty(mod_chosen.mdl_type) && mod_chosen.mdl_n ~= 0
                fprintf([mod_chosen.mdl_type,' model ',num2str(mod_chosen.mdl_n),' ']);
            end
            %% SV (gain pair)
            if mod_chosen.SV ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.SV
                    case 1
                        fprintf('pGainItem*QgainItem+pNeutralItem*QneutralItem');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            
            %% dQ (gain pair)
            if mod_chosen.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.dQ
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
                    case 6
                        fprintf('(Q gain cue - Q neutral cue)^2');
                    case 7
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2)');
                    case 8
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice=best option) (gain pair)
            if mod_chosen.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.pBest
                    case 1
                        fprintf('[p(choice=best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain pair)
            if ismember(mod_chosen.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (gain pair)
            if mod_chosen.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_chosen.ROI_activity_GLM,...
                    ' ',mod_chosen.ROI_activity_period,' period ',...
                    mod_chosen.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (gain pair)
            if mod_chosen.RT ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.RT
                    case {1,2,3,4}
                        RTcorr_nm = '';
                    case {5,6,7,8}
                        RTcorr_nm = ' corrected';
                end
                switch mod_chosen.RT
                    case {1,5}
                        fprintf(['RT',RTcorr_nm]);
                    case {2,6}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {3,7}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {4,8}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain pair)
            if ismember(mod_chosen.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            
            %% neutral pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Chosen stimulus display (neutral): ');
            switch dur_chosen
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from choice in red until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% trial number (neutral pair)
            if ismember(mod_chosen.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% no SV for neutral pair
            %% no dQ for neutral pair
            
            %% trial number (neutral pair)
            if ismember(mod_chosen.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (neutral pair)
            if mod_chosen.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_chosen.ROI_activity_GLM,...
                    ' ',mod_chosen.ROI_activity_period,' period ',...
                    mod_chosen.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (neutral pair)
            if mod_chosen.RT ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.RT
                    case {1,2,3,4}
                        RTcorr_nm = '';
                    case {5,6,7,8}
                        RTcorr_nm = ' corrected';
                end
                switch mod_chosen.RT
                    case {1,5}
                        fprintf(['RT',RTcorr_nm]);
                    case {2,6}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {3,7}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {4,8}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (neutral pair)
            if ismember(mod_chosen.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            
            %% loss pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Chosen stimulus display (loss): ');
            switch dur_chosen
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from choice in red until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% trial number (loss pair)
            if ismember(mod_chosen.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% model used for SV, dQ and p(choice=best)
            if ~isempty(mod_chosen.mdl_type) && mod_chosen.mdl_n ~= 0
                fprintf([mod_chosen.mdl_type,' model ',num2str(mod_chosen.mdl_n),' ']);
            end
            %% SV (loss pair)
            if mod_chosen.SV ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.SV
                    case 1
                        fprintf('pLossItem*QlossItem+pNeutralItem*QneutralItem');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% dQ (loss pair)
            if mod_chosen.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.dQ
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
                    case 6
                        fprintf('(Q neutral cue - Q loss cue)^2');
                    case 7
                        fprintf('|Q neutral cue - Q loss cue|/sqrt((QneutralCue+0.5)^2+(QlossCue+0.5)^2)');
                    case 8
                        fprintf('|Q neutral cue - Q loss cue|/sqrt((QneutralCue+0.5)^2+(QlossCue+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice = best option)
            if mod_chosen.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.pBest
                    case 1
                        fprintf('[p(choice = best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (loss pair)
            if ismember(mod_chosen.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (loss pair)
            if mod_chosen.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_chosen.ROI_activity_GLM,...
                    ' ',mod_chosen.ROI_activity_period,' period ',...
                    mod_chosen.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (loss pair)
            if mod_chosen.RT ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.RT
                    case {1,2,3,4}
                        RTcorr_nm = '';
                    case {5,6,7,8}
                        RTcorr_nm = ' corrected';
                end
                switch mod_chosen.RT
                    case {1,5}
                        fprintf(['RT',RTcorr_nm]);
                    case {2,6}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {3,7}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {4,8}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (loss pair)
            if ismember(mod_chosen.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
        case 3
            %% gain + loss pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Chosen stimulus display (gain + loss): ');
            switch dur_chosen
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from choice in red until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% trial number (gain + loss pair)
            if ismember(mod_chosen.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% model used for SV, dQ and p(choice=best)
            if ~isempty(mod_chosen.mdl_type) && mod_chosen.mdl_n ~= 0
                fprintf([mod_chosen.mdl_type,' model ',num2str(mod_chosen.mdl_n),' ']);
            end
            %% SV (gain + loss pair)
            if mod_chosen.SV ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.SV
                    case 1
                        fprintf('pBestItem*QbestItem+pWorseItem*QworseItem');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            
            %% dQ (gain + loss pair)
            if mod_chosen.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.dQ
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
                    case 6
                        fprintf('(Q gain cue - Q neutral cue)^2');
                    case 7
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2)');
                    case 8
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice=best option) (gain + loss pair)
            if mod_chosen.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.pBest
                    case 1
                        fprintf('[p(choice=best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain + loss pair)
            if ismember(mod_chosen.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (gain + loss pair)
            if mod_chosen.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_chosen.ROI_activity_GLM,...
                    ' ',mod_chosen.ROI_activity_period,' period ',...
                    mod_chosen.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (gain + loss pair)
            if mod_chosen.RT ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.RT
                    case {1,2,3,4}
                        RTcorr_nm = '';
                    case {5,6,7,8}
                        RTcorr_nm = ' corrected';
                end
                switch mod_chosen.RT
                    case {1,5}
                        fprintf(['RT',RTcorr_nm]);
                    case {2,6}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {3,7}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {4,8}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain + loss pair)
            if ismember(mod_chosen.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            
            %% neutral pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Chosen stimulus display (neutral): ');
            switch dur_chosen
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from choice in red until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% trial number (neutral pair)
            if ismember(mod_chosen.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% no SV for neutral pair
            %% no dQ for neutral pair
            
            %% trial number (neutral pair)
            if ismember(mod_chosen.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (neutral pair)
            if mod_chosen.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_chosen.ROI_activity_GLM,...
                    ' ',mod_chosen.ROI_activity_period,' period ',...
                    mod_chosen.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (neutral pair)
            if mod_chosen.RT ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.RT
                    case {1,2,3,4}
                        RTcorr_nm = '';
                    case {5,6,7,8}
                        RTcorr_nm = ' corrected';
                end
                switch mod_chosen.RT
                    case {1,5}
                        fprintf(['RT',RTcorr_nm]);
                    case {2,6}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {3,7}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {4,8}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (neutral pair)
            if ismember(mod_chosen.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
        case 4
            error('not ready yet');
        case 5
            %% gain pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Chosen stimulus display (gain first trials): ');
            switch dur_chosen
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from choice in red until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% trial number (gain pair first trials)
            if ismember(mod_chosen.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% model used for SV, dQ and p(choice=best)
            if ~isempty(mod_chosen.mdl_type) && mod_chosen.mdl_n ~= 0
                fprintf([mod_chosen.mdl_type,' model ',num2str(mod_chosen.mdl_n),' ']);
            end
            %% SV (gain pair first trials)
            if mod_chosen.SV ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.SV
                    case 1
                        fprintf('pGainItem*QgainItem+pNeutralItem*QneutralItem');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            
            %% dQ (gain pair first trials)
            if mod_chosen.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.dQ
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
                    case 6
                        fprintf('(Q gain cue - Q neutral cue)^2');
                    case 7
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2)');
                    case 8
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice=best option) (gain pair first trials)
            if mod_chosen.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.pBest
                    case 1
                        fprintf('[p(choice=best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain pair first trials)
            if ismember(mod_chosen.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (gain pair first trials)
            if mod_chosen.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_chosen.ROI_activity_GLM,...
                    ' ',mod_chosen.ROI_activity_period,' period ',...
                    mod_chosen.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (gain pair first trials)
            if mod_chosen.RT ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.RT
                    case {1,2,3,4}
                        RTcorr_nm = '';
                    case {5,6,7,8}
                        RTcorr_nm = ' corrected';
                end
                switch mod_chosen.RT
                    case {1,5}
                        fprintf(['RT',RTcorr_nm]);
                    case {2,6}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {3,7}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {4,8}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain pair first trials)
            if ismember(mod_chosen.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            
            %% neutral pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Chosen stimulus display (neutral first trials): ');
            switch dur_chosen
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from choice in red until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% trial number (neutral pair first trials)
            if ismember(mod_chosen.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% no SV for neutral pair
            %% no dQ for neutral pair
            
            %% trial number (neutral pair, first trials)
            if ismember(mod_chosen.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (neutral pair first trials)
            if mod_chosen.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_chosen.ROI_activity_GLM,...
                    ' ',mod_chosen.ROI_activity_period,' period ',...
                    mod_chosen.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (neutral pair first trials)
            if mod_chosen.RT ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.RT
                    case {1,2,3,4}
                        RTcorr_nm = '';
                    case {5,6,7,8}
                        RTcorr_nm = ' corrected';
                end
                switch mod_chosen.RT
                    case {1,5}
                        fprintf(['RT',RTcorr_nm]);
                    case {2,6}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {3,7}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {4,8}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (neutral pair first trials)
            if ismember(mod_chosen.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            
            %% loss pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Chosen stimulus display (loss first trials): ');
            switch dur_chosen
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from choice in red until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% trial number (loss pair first trials)
            if ismember(mod_chosen.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% model used for SV, dQ and p(choice=best)
            if ~isempty(mod_chosen.mdl_type) && mod_chosen.mdl_n ~= 0
                fprintf([mod_chosen.mdl_type,' model ',num2str(mod_chosen.mdl_n),' ']);
            end
            %% SV (loss pair first trials)
            if mod_chosen.SV ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.SV
                    case 1
                        fprintf('pLossItem*QlossItem+pNeutralItem*QneutralItem');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% dQ (loss pair first trials)
            if mod_chosen.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.dQ
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
                    case 6
                        fprintf('(Q neutral cue - Q loss cue)^2');
                    case 7
                        fprintf('|Q neutral cue - Q loss cue|/sqrt((QneutralCue+0.5)^2+(QlossCue+0.5)^2)');
                    case 8
                        fprintf('|Q neutral cue - Q loss cue|/sqrt((QneutralCue+0.5)^2+(QlossCue+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice = best option) (loss pair first trials)
            if mod_chosen.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.pBest
                    case 1
                        fprintf('[p(choice = best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (loss pair first trials)
            if ismember(mod_chosen.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (loss pair first trials)
            if mod_chosen.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_chosen.ROI_activity_GLM,...
                    ' ',mod_chosen.ROI_activity_period,' period ',...
                    mod_chosen.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (loss pair first trials)
            if mod_chosen.RT ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.RT
                    case {1,2,3,4}
                        RTcorr_nm = '';
                    case {5,6,7,8}
                        RTcorr_nm = ' corrected';
                end
                switch mod_chosen.RT
                    case {1,5}
                        fprintf(['RT',RTcorr_nm]);
                    case {2,6}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {3,7}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {4,8}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (loss pair first trials)
            if ismember(mod_chosen.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            %% gain pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Chosen stimulus display (gain last trials): ');
            switch dur_chosen
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from choice in red until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% trial number (gain pair last trials)
            if ismember(mod_chosen.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% model used for SV, dQ and p(choice=best)
            if ~isempty(mod_chosen.mdl_type) && mod_chosen.mdl_n ~= 0
                fprintf([mod_chosen.mdl_type,' model ',num2str(mod_chosen.mdl_n),' ']);
            end
            %% SV (gain pair last trials)
            if mod_chosen.SV ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.SV
                    case 1
                        fprintf('pGainItem*QgainItem+pNeutralItem*QneutralItem');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            
            %% dQ (gain pair last trials)
            if mod_chosen.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.dQ
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
                    case 6
                        fprintf('(Q gain cue - Q neutral cue)^2');
                    case 7
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2)');
                    case 8
                        fprintf('|Q gain cue - Q neutral cue|/sqrt((QgainCue+0.5)^2+(QneutralCue+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice=best option) (gain pair last trials)
            if mod_chosen.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.pBest
                    case 1
                        fprintf('[p(choice=best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain pair last trials)
            if ismember(mod_chosen.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (gain pair last trials)
            if mod_chosen.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_chosen.ROI_activity_GLM,...
                    ' ',mod_chosen.ROI_activity_period,' period ',...
                    mod_chosen.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (gain pair last trials)
            if mod_chosen.RT ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.RT
                    case {1,2,3,4}
                        RTcorr_nm = '';
                    case {5,6,7,8}
                        RTcorr_nm = ' corrected';
                end
                switch mod_chosen.RT
                    case {1,5}
                        fprintf(['RT',RTcorr_nm]);
                    case {2,6}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {3,7}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {4,8}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (gain pair last trials)
            if ismember(mod_chosen.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            
            %% neutral pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Chosen stimulus display (neutral last trials): ');
            switch dur_chosen
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from choice in red until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% trial number (neutral pair last trials)
            if ismember(mod_chosen.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% no SV for neutral pair
            %% no dQ for neutral pair
            
            %% trial number (neutral pair last trials)
            if ismember(mod_chosen.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (neutral pair last trials)
            if mod_chosen.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_chosen.ROI_activity_GLM,...
                    ' ',mod_chosen.ROI_activity_period,' period ',...
                    mod_chosen.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (neutral pair last trials)
            if mod_chosen.RT ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.RT
                    case {1,2,3,4}
                        RTcorr_nm = '';
                    case {5,6,7,8}
                        RTcorr_nm = ' corrected';
                end
                switch mod_chosen.RT
                    case {1,5}
                        fprintf(['RT',RTcorr_nm]);
                    case {2,6}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {3,7}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {4,8}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (neutral pair last trials)
            if ismember(mod_chosen.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            
            %% loss pairs onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Chosen stimulus display (loss last trials): ');
            switch dur_chosen
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from choice in red until feedback');
                otherwise
                    error('case not ready yet');
            end
            fprintf('\n');
            
            %% trial number (loss pair last trials)
            if ismember(mod_chosen.trialN,[1,2,3])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% model used for SV, dQ and p(choice=best)
            if ~isempty(mod_chosen.mdl_type) && mod_chosen.mdl_n ~= 0
                fprintf([mod_chosen.mdl_type,' model ',num2str(mod_chosen.mdl_n),' ']);
            end
            %% SV (loss pair last trials)
            if mod_chosen.SV ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.SV
                    case 1
                        fprintf('pLossItem*QlossItem+pNeutralItem*QneutralItem');
                    case 2
                        fprintf('QA+QB');
                    case 3
                        fprintf('Qch-Qunch');
                    case 4
                        fprintf('Qch');
                    case 5
                        fprintf('Qch./(QA+QB)');
                    otherwise
                        warning('SV formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% dQ (loss pair last trials)
            if mod_chosen.dQ ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.dQ
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
                    case 6
                        fprintf('(Q neutral cue - Q loss cue)^2');
                    case 7
                        fprintf('|Q neutral cue - Q loss cue|/sqrt((QneutralCue+0.5)^2+(QlossCue+0.5)^2)');
                    case 8
                        fprintf('|Q neutral cue - Q loss cue|/sqrt((QneutralCue+0.5)^2+(QlossCue+0.5)^2) (zPerRun)');
                    otherwise
                        warning('dQ formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% p(choice = best option) (loss pair last trials)
            if mod_chosen.pBest ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.pBest
                    case 1
                        fprintf('[p(choice = best option) - 0.5]^2');
                    case 2
                        fprintf('[p(choice=left option) - 0.5]^2');
                    case 3
                        fprintf('p(choice=chosen option)');
                    otherwise
                        warning('pBest formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (loss pair last trials)
            if ismember(mod_chosen.trialN,[4,5,6])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 4
                        fprintf('trial number');
                    case 5
                        fprintf('zPerRun trial number');
                    case 6
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (loss pair last trials)
            if mod_chosen.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_chosen.ROI_activity_GLM,...
                    ' ',mod_chosen.ROI_activity_period,' period ',...
                    mod_chosen.ROI_activity_ROI_nm,' activity ']);
            end
            %% RT (loss pair last trials)
            if mod_chosen.RT ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.RT
                    case {1,2,3,4}
                        RTcorr_nm = '';
                    case {5,6,7,8}
                        RTcorr_nm = ' corrected';
                end
                switch mod_chosen.RT
                    case {1,5}
                        fprintf(['RT',RTcorr_nm]);
                    case {2,6}
                        fprintf(['RT',RTcorr_nm, ' zPerRun']);
                    case {3,7}
                        fprintf(['RT',RTcorr_nm,' zAcrossRuns']);
                    case {4,8}
                        fprintf(['RT',RTcorr_nm,'zAcrossRuns and zPerRun']);
                    otherwise
                        warning('RT formula not determined in which_GLM_MS2_infos_RL.m');
                end
                fprintf('/');
            end
            %% trial number (loss pair last trials)
            if ismember(mod_chosen.trialN,[7,8,9])
                n_prm_RL = n_prm_RL + 1;
                switch mod_chosen.trialN
                    case 7
                        fprintf('trial number');
                    case 8
                        fprintf('zPerRun trial number');
                    case 9
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
    end
end

%% feedback
o_fbk = RL.o_fbk;
dur_fbk = RL.dur_fbk;
mod_fbk = RL.mod_fbk;
if o_fbk ~= 0
    switch o_fbk
        case 1 % group all feedbacks
            %% feedback onset + duration (all pairs)
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Feedback (all pairs): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from feedback onset to cross display');
            end
            fprintf('\n');
            
            %% trial number (all pairs feedback)
            if mod_fbk.trialN ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% feedback (all pairs feedback)
            if mod_fbk.fbk ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.fbk
                    case 1
                        fprintf('feedback (current trial)');
                end
                fprintf('/');
            end
            %% prediction error (all pairs feedback)
            if mod_fbk.PE ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% prediction error bis (all pairs feedback)
            if mod_fbk.PE_bis ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE_bis
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% total gain (all pairs feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (all pairs, feedback)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            
        case 2 % slit by pair type gain pair/neutral pair/loss pair
            
            %% feedback onset + duration (gain pair)
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Feedback (gain pair): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from feedback onset to cross display');
            end
            fprintf('\n');
            
            %% trial number (gain pair feedback)
            if mod_fbk.trialN ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% feedback (gain pair feedback)
            if mod_fbk.fbk ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.fbk
                    case 1
                        fprintf('feedback (current trial)');
                end
                fprintf('/');
            end
            %% prediction error (gain pair feedback)
            if mod_fbk.PE ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% prediction error bis (gain pair feedback)
            if mod_fbk.PE_bis ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE_bis
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% total gain (gain pair feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (gain pair, feedback)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %%
            fprintf('\n');
            
            %% Feedback neutral pair onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Feedback (neutral pair): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from feedback onset to cross display');
            end
            fprintf('\n');
            %% trial number (neutral pair feedback)
            if mod_fbk.trialN ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% feedback and PE pointless for neutral pair => always zero
            %% total gain (neutral pair feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (neutral pair, feedback)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %%
            fprintf('\n');
            
            %% Feedback loss pair onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Feedback (loss pair): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from feedback onset to cross display');
            end
            fprintf('\n');
            
            %% trial number (loss pair feedback)
            if mod_fbk.trialN ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% feedback (loss pair feedback)
            if mod_fbk.fbk ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.fbk
                    case 1
                        fprintf('feedback (current trial)');
                end
                fprintf('/');
            end
            %% prediction error (loss pair feedback)
            if mod_fbk.PE ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% prediction error bis (loss pair feedback)
            if mod_fbk.PE_bis ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE_bis
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% total gain (loss pair feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
                fprintf('/');
            end
            fprintf('\n');
            %% ROI activity as regressor (loss pair, feedback)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            
        case 3 % split by feedback type gain feedback/neutral feedback/loss feedback
            %% feedback onset + duration (gain feedback)
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Feedback (gain feedback): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from feedback onset to cross display');
            end
            fprintf('\n');
            %% trial number (gain feedback)
            if mod_fbk.trialN ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% feedback (gain feedback)
            if mod_fbk.fbk ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.fbk
                    case 1
                        fprintf('feedback (current trial)');
                end
                fprintf('/');
            end
            %% prediction error (gain feedback)
            if mod_fbk.PE ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% prediction error bis (gain feedback)
            if mod_fbk.PE_bis ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE_bis
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% totalGain (gain feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (gain feedback, feedback)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %%
            fprintf('\n');
            
            %%  feedback (neutral feedback)
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Feedback (neutral feedback): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from feedback onset to cross display');
            end
            fprintf('\n');
            %% trial number (neutral feedback)
            if mod_fbk.trialN ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% feedback (neutral feedback)
            
            %% prediction error (neutral feedback)
            if mod_fbk.PE ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% prediction error bis (neutral feedback)
            if mod_fbk.PE_bis ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE_bis
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            
            %% total gain (neutral feedback)
            
            %% ROI activity as regressor (neutral feedback, feedback)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %%
            fprintf('\n');
            
            %% feedback onset + duration (loss feedback)
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Feedback (loss feedback): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from feedback onset to cross display');
            end
            fprintf('\n');
            %% trial number (loss feedback)
            if mod_fbk.trialN ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% feedback (loss feedback)
            if mod_fbk.fbk ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.fbk
                    case 1
                        fprintf('feedback (current trial)');
                end
                fprintf('/');
            end
            %% prediction error (loss feedback)
            if mod_fbk.PE ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% prediction error bis (loss feedback)
            if mod_fbk.PE_bis ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE_bis
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% total gain (loss feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (loss feedback, feedback)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %%
            fprintf('\n');
            
        case 4 % split by gain and feedback type gain pair gain feedback/gain pair neutral feedback/neutral pair & feedback/loss pair neutral feedback/loss pair loss feedback
            %% feedback onset + duration (gain pair, gain feedback)
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Feedback (gain pair - gain feedback): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from feedback onset to cross display');
            end
            fprintf('\n');
            %% trial number (gain pair, gain feedback)
            if mod_fbk.trialN ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% feedback (gain pair, gain feedback)
            if mod_fbk.fbk ~= 0
                error('Pointless: feedback = stable for a given feedback type.');
            end
            %% prediction error (gain pair, gain feedback)
            if mod_fbk.PE ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% prediction error bis (gain pair, gain feedback)
            if mod_fbk.PE_bis ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE_bis
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% total gain (gain pair, gain feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (gain pair, gain feedback)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %%
            fprintf('\n');
            
            %% feedback onset + duration (gain pair, neutral feedback)
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Feedback (gain pair - neutral feedback): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from feedback onset to cross display');
            end
            fprintf('\n');
            %% trial number (gain pair, neutral feedback)
            if mod_fbk.trialN ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% feedback (gain pair, neutral feedback) : always zero
            if mod_fbk.fbk ~= 0
                error('Pointless: feedback = stable for a given feedback type.');
            end
            %% prediction error (gain pair, neutral feedback): no change
            
            %% prediction error bis (gain pair, neutral feedback): no change
            
            %% total gain (gain pair, neutral feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (gain pair, neutral feedback)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %%
            fprintf('\n');
            
            %% feedback onset + duration (neutral pair, neutral feedback)
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Feedback (neutral pair & feedback): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from feedback onset to cross display');
            end
            fprintf('\n');
            %% trial number (neutral pair, neutral feedback)
            if mod_fbk.trialN ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% neutral pair neutral feedback => always zero
            %% PE: no change
            %% PE bis: no change
            %% total gain (neutral pair, neutral feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (neutral pair, neutral feedback)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %%
            fprintf('\n');
            
            %% feedback onset + duration (loss pair, neutral feedback)
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Feedback (loss pair - neutral feedback): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from feedback onset to cross display');
            end
            fprintf('\n');
            %% trial number (loss pair, neutral feedback)
            if mod_fbk.trialN ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% feedback (loss pair, neutral feedback)
            if mod_fbk.fbk ~= 0
                error('Pointless: feedback = stable for a given feedback type.');
            end
            %% prediction error (loss pair, neutral feedback)
            if mod_fbk.PE ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% prediction error bis (loss pair, neutral feedback)
            if mod_fbk.PE_bis ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE_bis
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% total gain (loss pair, neutral feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (loss pair, neutral feedback)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %%
            fprintf('\n');
            
            %% feedback (loss pair, loss feedback)
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Feedback (loss pair - loss feedback): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from feedback onset to cross display');
            end
            fprintf('\n');
            %% trial number (loss pair, loss feedback)
            if mod_fbk.trialN ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% feedback (loss pair, loss feedback)
            if mod_fbk.fbk ~= 0
                error('Pointless: feedback = stable for a given feedback type.');
            end
            %% prediction error (loss pair, loss feedback)
            if mod_fbk.PE ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% prediction error bis (loss pair, loss feedback)
            if mod_fbk.PE_bis ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE_bis
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% total gain (loss pair, loss feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (loss pair, loss feedback)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %%
            fprintf('\n');
            
        case 5 % pool gain+loss/leave neutral apart
            
            %% feedback onset + duration (gain + loss pairs)
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Feedback (gain + loss pair): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from feedback onset to cross display');
            end
            fprintf('\n');
            
            %% trial number (gain + loss pairs feedback)
            if mod_fbk.trialN ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% feedback (gain + loss pairs feedback)
            if mod_fbk.fbk ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.fbk
                    case 1
                        fprintf('feedback (current trial)');
                end
                fprintf('/');
            end
            %% prediction error (gain + loss pairs feedback)
            if mod_fbk.PE ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% prediction error bis (gain + loss pairs feedback)
            if mod_fbk.PE_bis ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE_bis
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% total gain (gain + loss pairs feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (gain + loss pairs, feedback)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %%
            fprintf('\n');
            
            %% Feedback neutral pair onset + duration
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Feedback (neutral pair): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from feedback onset to cross display');
            end
            fprintf('\n');
            %% trial number (neutral pair feedback)
            if mod_fbk.trialN ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% feedback + PE pointless for neutral pair => always zero
            %% total gain (neutral pair feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (neutral pair, feedback)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %%
            fprintf('\n');
            
        case 6 % pool gain+loss/leave neutral apart and also split first/second half of the trials
            
            %% feedback onset + duration (gain + loss pairs, first trials)
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Feedback (gain + loss pair - first trials): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from feedback onset to cross display');
            end
            fprintf('\n');
            
            %% trial number (gain + loss pairs feedback, first trials)
            if mod_fbk.trialN ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% feedback (gain + loss pairs feedback, first trials)
            if mod_fbk.fbk ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.fbk
                    case 1
                        fprintf('feedback (current trial)');
                end
                fprintf('/');
            end
            %% prediction error (gain + loss pairs feedback, first trials)
            if mod_fbk.PE ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% prediction error bis (gain + loss pairs feedback, first trials)
            if mod_fbk.PE_bis ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE_bis
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% total gain (gain + loss pairs feedback, first trials)
            if mod_fbk.totalGain ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (gain + loss pairs, feedback, first trials)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %%
            fprintf('\n');
            
            %% feedback onset + duration (gain + loss pairs, last trials)
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Feedback (gain + loss pair - last trials): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from feedback onset to cross display');
            end
            fprintf('\n');
            
            %% trial number (gain + loss pairs feedback, last trials)
            if mod_fbk.trialN ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% feedback (gain + loss pairs feedback, last trials)
            if mod_fbk.fbk ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.fbk
                    case 1
                        fprintf('feedback (current trial)');
                end
                fprintf('/');
            end
            %% prediction error (gain + loss pairs feedback, last trials)
            if mod_fbk.PE ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% prediction error bis (gain + loss pairs feedback, last trials)
            if mod_fbk.PE_bis ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.PE_bis
                    case 1
                        fprintf('Prediction Error');
                    case 2
                        fprintf('|Prediction Error|');
                end
                fprintf('/');
            end
            %% total gain (gain + loss pairs feedback, last trials)
            if mod_fbk.totalGain ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (gain + loss pairs, feedback, last trials)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %%
            fprintf('\n');
            
            %% Feedback neutral pair onset + duration, last trials
            n_prm_RL = n_prm_RL + 1;
            fprintf('**Feedback (neutral pair): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from feedback onset to cross display');
            end
            fprintf('\n');
            %% trial number (neutral pair feedback)
            if mod_fbk.trialN ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trial number');
                    case 2
                        fprintf('zPerRun trial number');
                    case 3
                        fprintf('zAcrossRuns trial number');
                end
                fprintf('/');
            end
            %% feedback + PE pointless for neutral pair => always zero
            %% total gain (neutral pair feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_RL = n_prm_RL + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
                fprintf('/');
            end
            %% ROI activity as regressor (neutral pair, feedback)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_RL = n_prm_RL + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %%
            fprintf('\n');
    end % o_fbk
end % ofbk # 0


%% cross
o_cross = RL.o_cross;
dur_cross = RL.dur_cross;
if o_cross ~= 0
    switch o_cross
        case 1
            %% cross onset + duration (all pairs)
            n_prm_RL = n_prm_RL + 1;
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


%% missed trials
o_missed_trials_stim = RL.o_missed_trials_stim;
dur_missed_trials_stim = RL.dur_missed_trials_stim;
switch o_missed_trials_stim
    case 1
        %% missed trials onset + duration (all trial types)
        %         n_prm_RL = n_prm_RL + 1; % depends on the subjects and runs
        fprintf('**Missed Trials (all pairs): ');
        switch dur_missed_trials_stim
            case 0
                fprintf('stick function');
            case 1
                fprintf('boxcar from display of stimuli to end of trial');
        end
end

end % function