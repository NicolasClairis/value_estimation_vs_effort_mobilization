function [ n_prm_stroop ] = which_GLM_MS2_infos_stroop( GLMprm )
%which_GLM_MS2_infos_stroop( GLMprm ) outputs the informations about the GLM based
% on the list of parameters designated inside GLMprm for the
% stroop task.
%
% INPUTS
% GLMprm: structure with GLM parameters for each task extracted with
% which_GLM_MS2.
%
% OUTPUTS
% n_prm_stroop: number of parameters in stroop task (includes onsets and
% regressors) for a given run (without taking into account derivative)
%
% See also which_GLM_MS2

%% stroop parameters info
stroop  = GLMprm.stroop;
disp('***Stroop parameters');
n_prm_stroop = 0;

%% incentive display
o_inc = stroop.o_inc;
dur_inc = stroop.dur_inc;
mod_inc = stroop.mod_inc;
if o_inc ~= 0
    switch o_inc
        case 1 % all trials grouped
            %% incentive onset + duration
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Incentive display (all): ');
            switch dur_inc
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset to effort scale');
                case 2
                    fprintf('boxcar from onset to feedback');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% Gain/Loss condition (incentive, all trials)
            if mod_inc.GL_cond ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.GL_cond
                    case 1
                        fprintf('Gain or Loss condition (1/0)');
                    case 2
                        fprintf('Gain or Loss condition (1/-1)');
                end
                fprintf('/');
            end
            %% luminance (first regressor) (incentive, all trials)
            if mod_inc.lum ~= 0 && ismember(mod_inc.lum,[2,4,6,8])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.lum
                    case 2
                        fprintf('luminance');
                    case 4
                        fprintf('luminance zPerRun');
                    case 6
                        fprintf('luminance zAcrossRuns');
                    case 8
                        fprintf('luminance zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (first regressor) (incentive, all trials)
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first pair (incentive, all trials)
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[1])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.RT_fp
                    case 2
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% model variables
            if mod_inc.mdl_n ~= 0
                %% ressource (incentive, all trials)
                if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[5,6])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.E_pred
                        case 5
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                        case 6
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% piece/bill condition (incentive, all trials)
            if mod_inc.R_type ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.R_type
                    case 1
                        fprintf('reward type: piece of money or bill (0/1)');
                    otherwise
                        error('not ready');
                end
                fprintf('/');
            end
            %% incentive (incentive, all trials)
            if mod_inc.inc ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.inc
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivational value');
                    case 4
                        fprintf('incentive motivational value');
                end
                fprintf('/');
            end
            %% incentive bis (incentive, all trials)
            if mod_inc.inc_bis ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.inc_bis
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivational value');
                    case 4
                        fprintf('incentive motivational value');
                end
                fprintf('/');
            end
            %% confidence proxy (incentive, all trials)
            if mod_inc.conf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.conf
                    case 1
                        fprintf('confidence (incentive rank centered)');
                    case 2
                        fprintf('confidence (incentive money centered)');
                    otherwise
                        error(['mod_inc.conf = ',num2str(mod_inc.conf),' not ready yet']);
                end
                fprintf('/');
            end
            %% ROI activity as regressor (incentive, all trials)
            if mod_inc.ROI_activity_yn == 1
                n_prm_stroop = n_prm_stroop + 1;
                fprintf(['GLM',mod_inc.ROI_activity_GLM,...
                    ' ',mod_inc.ROI_activity_period,' period ',...
                    mod_inc.ROI_activity_ROI_nm,' activity ']);
            end
            %% number of pairs solved (effort #1) (incentive, all trials)
            if mod_inc.n_pairs_solved ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.n_pairs_solved
                    case 1
                        fprintf('Number of pairs solved per trial');
                end
                fprintf('/');
            end
            %% percentage of incongruent pairs (effort #2) (incentive, all trials)
            if mod_inc.n_incong ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.n_incong
                    case 1
                        fprintf('Number of incongruent pairs solved');
                    case 2
                        fprintf('Percentage of incongruent pairs solved per trial');
                end
                fprintf('/');
            end
            %% error numbers (incentive, all trials)
            if mod_inc.n_errors ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.n_errors
                    case 1
                        fprintf('Trial with at least one error (1) or no error at all (0)');
                    case 2
                        fprintf('Number of errors made per trial across pairs'); % includes errors repeated on the same pair
                    case 3
                        fprintf('Percentage of errors made per trial');
                end
                fprintf('/');
            end
            %% model variables
            if mod_inc.mdl_n ~= 0
                %% ressource (incentive, all trials)
                if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[3,4])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.E_pred
                        case 3
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                        case 4
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
                %% cost (incentive, all trials)
                if mod_inc.mdl_cost ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.mdl_cost
                        case 1
                            fprintf(['model ',num2str(mod_inc.mdl_n),' cost based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_inc.mdl_n),' cost based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_inc.mdl_n),' cost (before perf/X scaling)']);
                    end
                    fprintf('/');
                end
                %% benefit (incentive, all trials)
                if mod_inc.mdl_benefit ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.mdl_benefit
                        case 1
                            fprintf(['model ',num2str(mod_inc.mdl_n),' benefit based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_inc.mdl_n),' benefit based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_inc.mdl_n),' benefit (before perf/X scaling)']);
                        case 4
                            fprintf(['model ',num2str(mod_inc.mdl_n),' expected monetary payoff']);
                    end
                    fprintf('/');
                end
                %% EV (incentive, all trials)
                if mod_inc.mdl_EV ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.mdl_EV
                        case 1
                            fprintf(['model ',num2str(mod_inc.mdl_n),' EV based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_inc.mdl_n),' EV based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% trial number (middle regressor) (incentive, all trials)
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[5,6])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.trialN
                    case 5
                        fprintf('trialN');
                    case 6
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first pair (incentive, all trials)
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.RT_fp
                    case 1
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% mean RT all but first pair (incentive, all trials)
            if mod_inc.RT_mRT ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.RT_mRT
                    case 1
                        fprintf('mean(RT) (across pairs except 1st pair)');
                    case 2
                        fprintf('mean(RT) (across pairs except 1st pair and errors)');
                end
                fprintf('/');
            end
            %% cumulated total gain until now (incentive, all trials)
            if mod_inc.totalGain_prev ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.totalGain_prev
                    case 1
                        fprintf('total gain (until this trial)');
                end
                fprintf('/');
            end
            %% sum(performance previous trials) until now (incentive, all trials)
            if mod_inc.sumPerfPrev ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.sumPerfPrev
                    case 1
                        fprintf('sum(performance previous trials) (raw values)');
                    case 2
                        fprintf('sum(performance previous trials) (normalized values)');
                end
                fprintf('/');
            end
            %% luminance (last regressor) (incentive, all trials)
            if mod_inc.lum ~= 0 && ismember(mod_inc.lum,[1,3,5,7])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.lum
                    case 1
                        fprintf('luminance');
                    case 3
                        fprintf('luminance zPerRun');
                    case 5
                        fprintf('luminance zAcrossRuns');
                    case 7
                        fprintf('luminance zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (last regressor) (incentive, all trials)
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% incentive*trial number (incentive, all trials)
            if mod_inc.inc_x_trialN ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.inc_x_trialN
                    case 1
                        fprintf('inc_x_trialN');
                end
                fprintf('/');
            end
            %% ressource (incentive, all trials)
            if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[1,2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.E_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (incentive, all trials)
            if mod_inc.perf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.perf
                    case 1
                        fprintf('performance (1-100)');
                    case 2
                        fprintf('performance (0-1)');
                    case 3
                        fprintf('performance (0-1) corrected');
                    case 4
                        fprintf(['performance (0-1) predicted by model ',num2str(mod_inc.mdl_n)]);
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
        case 2 % split by to gain/to lose trial type
            %% incentive onset + duration (to gain trials)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Incentive display (to gain): ');
            switch dur_inc
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset to effort scale');
                case 2
                    fprintf('boxcar from onset to feedback');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% Luminance (first regressor) (incentive, to gain trials)
            if mod_inc.lum ~= 0 && ismember(mod_inc.lum,[2,4,6,8])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.lum
                    case 2
                        fprintf('luminance');
                    case 4
                        fprintf('luminance zPerRun');
                    case 6
                        fprintf('luminance zAcrossRuns');
                    case 8
                        fprintf('luminance zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (first regressor) (incentive, to gain trials)
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first pair (incentive, to gain trials)
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[1])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.RT_fp
                    case 2
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% model variables
            if mod_inc.mdl_n ~= 0
                %% ressource (incentive, to gain trials)
                if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[5,6])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.E_pred
                        case 5
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                        case 6
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% piece/bill condition (incentive, to gain trials)
            if mod_inc.R_type ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.R_type
                    case 1
                        fprintf('reward type: piece of money or bill (0/1)');
                    otherwise
                        error('not ready');
                end
                fprintf('/');
            end
            %% incentive (incentive, to gain trials)
            if mod_inc.inc ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.inc
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivational value');
                    case 4
                        fprintf('incentive motivational value');
                end
                fprintf('/');
            end
            %% incentive bis (incentive, to gain trials)
            if mod_inc.inc_bis ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.inc_bis
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivational value');
                    case 4
                        fprintf('incentive motivational value');
                end
                fprintf('/');
            end
            %% confidence proxy (incentive, to gain trials)
            if mod_inc.conf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.conf
                    case 1
                        fprintf('confidence (incentive rank centered)');
                    case 2
                        fprintf('confidence (incentive money centered)');
                    otherwise
                        error(['mod_inc.conf = ',num2str(mod_inc.conf),' not ready yet']);
                end
                fprintf('/');
            end
            %% ROI activity as regressor (incentive, to gain trials)
            if mod_inc.ROI_activity_yn == 1
                n_prm_stroop = n_prm_stroop + 1;
                fprintf(['GLM',mod_inc.ROI_activity_GLM,...
                    ' ',mod_inc.ROI_activity_period,' period ',...
                    mod_inc.ROI_activity_ROI_nm,' activity ']);
            end
            %% number of pairs solved (effort #1) (incentive, to gain trials)
            if mod_inc.n_pairs_solved ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.n_pairs_solved
                    case 1
                        fprintf('Number of pairs solved per trial');
                end
                fprintf('/');
            end
            %% percentage of incongruent pairs (effort #2) (incentive, to gain trials)
            if mod_inc.n_incong ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.n_incong
                    case 1
                        fprintf('Number of incongruent pairs solved');
                    case 2
                        fprintf('Percentage of incongruent pairs solved per trial');
                end
                fprintf('/');
            end
            %% error numbers (incentive, to gain trials)
            if mod_inc.n_errors ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.n_errors
                    case 1
                        fprintf('Trial with at least one error (1) or no error at all (0)');
                    case 2
                        fprintf('Number of errors made per trial across pairs'); % includes errors repeated on the same pair
                    case 3
                        fprintf('Percentage of errors made per trial');
                end
                fprintf('/');
            end
            %% model variables
            if mod_inc.mdl_n ~= 0
                %% ressource (incentive, to gain trials)
                if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[3,4])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.E_pred
                        case 3
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                        case 4
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
                %% cost (incentive, to gain trials)
                if mod_inc.mdl_cost ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.mdl_cost
                        case 1
                            fprintf(['model ',num2str(mod_inc.mdl_n),' cost based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_inc.mdl_n),' cost based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_inc.mdl_n),' cost (before perf/X scaling)']);
                    end
                    fprintf('/');
                end
                %% benefit (incentive, to gain trials)
                if mod_inc.mdl_benefit ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.mdl_benefit
                        case 1
                            fprintf(['model ',num2str(mod_inc.mdl_n),' benefit based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_inc.mdl_n),' benefit based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_inc.mdl_n),' benefit (before perf/X scaling)']);
                        case 4
                            fprintf(['model ',num2str(mod_inc.mdl_n),' expected monetary payoff']);
                    end
                    fprintf('/');
                end
                %% EV (incentive, to gain trials)
                if mod_inc.mdl_EV ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.mdl_EV
                        case 1
                            fprintf(['model ',num2str(mod_inc.mdl_n),' EV based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_inc.mdl_n),' EV based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% trial number (middle regressor) (incentive, to gain trials)
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[5,6])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.trialN
                    case 5
                        fprintf('trialN');
                    case 6
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first pair (incentive, to gain trials)
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.RT_fp
                    case 1
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% mean RT all but first pair (incentive, to gain trials)
            if mod_inc.RT_mRT ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.RT_mRT
                    case 1
                        fprintf('mean(RT) (across pairs except 1st pair)');
                end
                fprintf('/');
            end
            %% cumulated total gain until now (incentive, to gain trials)
            if mod_inc.totalGain_prev ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.totalGain_prev
                    case 1
                        fprintf('total gain (until this trial)');
                end
                fprintf('/');
            end
            %% sum(performance previous trials) until now (incentive, to gain trials)
            if mod_inc.sumPerfPrev ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.sumPerfPrev
                    case 1
                        fprintf('sum(performance previous trials) (raw values)');
                    case 2
                        fprintf('sum(performance previous trials) (normalized values)');
                end
                fprintf('/');
            end
            %% luminance (last regressor) (incentive, to gain trials)
            if mod_inc.lum ~= 0 && ismember(mod_inc.lum,[1,3,5,7])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.lum
                    case 1
                        fprintf('luminance');
                    case 3
                        fprintf('luminance zPerRun');
                    case 5
                        fprintf('luminance zAcrossRuns');
                    case 7
                        fprintf('luminance zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (last regressor) (incentive, to gain trials)
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% incentive*trial number (incentive, to gain trials)
            if mod_inc.inc_x_trialN ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.inc_x_trialN
                    case 1
                        fprintf('inc_x_trialN');
                end
                fprintf('/');
            end
            %% ressource (incentive, to gain trials)
            if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[1,2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.E_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (incentive, to gain trials)
            if mod_inc.perf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.perf
                    case 1
                        fprintf('performance (1-100)');
                    case 2
                        fprintf('performance (0-1)');
                    case 3
                        fprintf('performance (0-1) corrected');
                    case 4
                        fprintf(['performance (0-1) predicted by model ',num2str(mod_inc.mdl_n)]);
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            %% incentive onset + duration (incentive, to lose trials)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Incentive display (to lose): ');
            switch dur_inc
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset to effort scale');
                case 2
                    fprintf('boxcar from onset to feedback');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% Luminance (first regressor) (incentive, to lose trials)
            if mod_inc.lum ~= 0 && ismember(mod_inc.lum,[2,4,6,8])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.lum
                    case 2
                        fprintf('luminance');
                    case 4
                        fprintf('luminance zPerRun');
                    case 6
                        fprintf('luminance zAcrossRuns');
                    case 8
                        fprintf('luminance zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (first regressor) (incentive, to lose trials)
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first pair (incentive, to lose trials)
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[1])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.RT_fp
                    case 2
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% model variables
            if mod_inc.mdl_n ~= 0
                %% ressource (incentive, to lose trials)
                if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[5,6])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.E_pred
                        case 5
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                        case 6
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% piece/bill condition (incentive, to lose trials)
            if mod_inc.R_type ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.R_type
                    case 1
                        fprintf('reward type: piece of money or bill (0/1)');
                    otherwise
                        error('not ready');
                end
                fprintf('/');
            end
            %% incentive (incentive, to lose trials)
            if mod_inc.inc ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.inc
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivational value');
                    case 4
                        fprintf('incentive motivational value');
                end
                fprintf('/');
            end
            %% incentive bis (incentive, to lose trials)
            if mod_inc.inc_bis ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.inc_bis
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivational value');
                    case 4
                        fprintf('incentive motivational value');
                end
                fprintf('/');
            end
            %% confidence proxy (incentive, to lose trials)
            if mod_inc.conf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.conf
                    case 1
                        fprintf('confidence (incentive rank centered)');
                    case 2
                        fprintf('confidence (incentive money centered)');
                    otherwise
                        error(['mod_inc.conf = ',num2str(mod_inc.conf),' not ready yet']);
                end
                fprintf('/');
            end
            %% ROI activity as regressor (incentive, to lose trials)
            if mod_inc.ROI_activity_yn == 1
                n_prm_stroop = n_prm_stroop + 1;
                fprintf(['GLM',mod_inc.ROI_activity_GLM,...
                    ' ',mod_inc.ROI_activity_period,' period ',...
                    mod_inc.ROI_activity_ROI_nm,' activity ']);
            end
            %% number of pairs solved (effort #1) (incentive, to lose trials)
            if mod_inc.n_pairs_solved ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.n_pairs_solved
                    case 1
                        fprintf('Number of pairs solved per trial');
                end
                fprintf('/');
            end
            %% percentage of incongruent pairs (effort #2) (incentive, to lose trials)
            if mod_inc.n_incong ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.n_incong
                    case 1
                        fprintf('Number of incongruent pairs solved');
                    case 2
                        fprintf('Percentage of incongruent pairs solved per trial');
                end
                fprintf('/');
            end
            %% error numbers (incentive, to lose trials)
            if mod_inc.n_errors ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.n_errors
                    case 1
                        fprintf('Trial with at least one error (1) or no error at all (0)');
                    case 2
                        fprintf('Number of errors made per trial across pairs'); % includes errors repeated on the same pair
                    case 3
                        fprintf('Percentage of errors made per trial');
                end
                fprintf('/');
            end
            %% model variables
            if mod_inc.mdl_n ~= 0
                %% ressource (incentive, to lose trials)
                if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[3,4])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.E_pred
                        case 3
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                        case 4
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
                %% cost (incentive, to lose trials)
                if mod_inc.mdl_cost ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.mdl_cost
                        case 1
                            fprintf(['model ',num2str(mod_inc.mdl_n),' cost based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_inc.mdl_n),' cost based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_inc.mdl_n),' cost (before perf/X scaling)']);
                    end
                    fprintf('/');
                end
                %% benefit (incentive, to lose trials)
                if mod_inc.mdl_benefit ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.mdl_benefit
                        case 1
                            fprintf(['model ',num2str(mod_inc.mdl_n),' benefit based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_inc.mdl_n),' benefit based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_inc.mdl_n),' benefit (before perf/X scaling)']);
                        case 4
                            fprintf(['model ',num2str(mod_inc.mdl_n),' expected monetary payoff']);
                    end
                    fprintf('/');
                end
                %% EV (incentive, to lose trials)
                if mod_inc.mdl_EV ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.mdl_EV
                        case 1
                            fprintf(['model ',num2str(mod_inc.mdl_n),' EV based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_inc.mdl_n),' EV based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% trial number (middle regressor) (incentive, to lose trials)
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[5,6])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.trialN
                    case 5
                        fprintf('trialN');
                    case 6
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first pair (incentive, to lose trials)
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.RT_fp
                    case 1
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% mean RT all but first pair (incentive, to lose trials)
            if mod_inc.RT_mRT ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.RT_mRT
                    case 1
                        fprintf('mean(RT) (across pairs except 1st pair)');
                end
                fprintf('/');
            end
            %% cumulated total gain until now (incentive, to lose trials)
            if mod_inc.totalGain_prev ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.totalGain_prev
                    case 1
                        fprintf('total gain (until this trial)');
                end
                fprintf('/');
            end
            %% sum(performance previous trials) until now (incentive, to lose trials)
            if mod_inc.sumPerfPrev ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.sumPerfPrev
                    case 1
                        fprintf('sum(performance previous trials) (raw values)');
                    case 2
                        fprintf('sum(performance previous trials) (normalized values)');
                end
                fprintf('/');
            end
            %% luminance (last regressor) (incentive, to lose trials)
            if mod_inc.lum ~= 0 && ismember(mod_inc.lum,[1,3,5,7])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.lum
                    case 1
                        fprintf('luminance');
                    case 3
                        fprintf('luminance zPerRun');
                    case 5
                        fprintf('luminance zAcrossRuns');
                    case 7
                        fprintf('luminance zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (last regressor) (incentive, to lose trials)
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% incentive*trial number (incentive, to lose trials)
            if mod_inc.inc_x_trialN ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.inc_x_trialN
                    case 1
                        fprintf('inc_x_trialN');
                end
                fprintf('/');
            end
            %% ressource (incentive, to lose trials)
            if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[1,2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.E_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (incentive, to lose trials)
            if mod_inc.perf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.perf
                    case 1
                        fprintf('performance (1-100)');
                    case 2
                        fprintf('performance (0-1)');
                    case 3
                        fprintf('performance (0-1) corrected');
                    case 4
                        fprintf(['performance (0-1) predicted by model ',num2str(mod_inc.mdl_n)]);
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
        case 3 % split by no-error/error trial type
            %% incentive onset + duration (incentive, no error trials)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Incentive display (no error trials): ');
            switch dur_inc
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset to effort scale');
                case 2
                    fprintf('boxcar from onset to feedback');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% Gain/Loss condition (incentive, no error trials)
            if mod_inc.GL_cond ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.GL_cond
                    case 1
                        fprintf('Gain or Loss condition (1/0)');
                    case 2
                        fprintf('Gain or Loss condition (1/-1)');
                end
                fprintf('/');
            end
            %% Luminance (first regressor) (incentive, no error trials)
            if mod_inc.lum ~= 0 && ismember(mod_inc.lum,[2,4,6,8])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.lum
                    case 2
                        fprintf('luminance');
                    case 4
                        fprintf('luminance zPerRun');
                    case 6
                        fprintf('luminance zAcrossRuns');
                    case 8
                        fprintf('luminance zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (first regressor) (incentive, no error trials)
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first pair (incentive, no error trials)
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[1])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.RT_fp
                    case 2
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% model variables
            if mod_inc.mdl_n ~= 0
                %% ressource (incentive, no error trials)
                if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[5,6])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.E_pred
                        case 5
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                        case 6
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% piece/bill condition (incentive, no error trials)
            if mod_inc.R_type ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.R_type
                    case 1
                        fprintf('reward type: piece of money or bill (0/1)');
                    otherwise
                        error('not ready');
                end
                fprintf('/');
            end
            %% incentive (incentive, no error trials)
            if mod_inc.inc ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.inc
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivational value');
                    case 4
                        fprintf('incentive motivational value');
                end
                fprintf('/');
            end
            %% incentive bis (incentive, no error trials)
            if mod_inc.inc_bis ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.inc_bis
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivational value');
                    case 4
                        fprintf('incentive motivational value');
                end
                fprintf('/');
            end
            %% confidence proxy (incentive, no error trials)
            if mod_inc.conf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.conf
                    case 1
                        fprintf('confidence (incentive rank centered)');
                    case 2
                        fprintf('confidence (incentive money centered)');
                    otherwise
                        error(['mod_inc.conf = ',num2str(mod_inc.conf),' not ready yet']);
                end
                fprintf('/');
            end
            %% ROI activity as regressor (incentive, no error trials)
            if mod_inc.ROI_activity_yn == 1
                n_prm_stroop = n_prm_stroop + 1;
                fprintf(['GLM',mod_inc.ROI_activity_GLM,...
                    ' ',mod_inc.ROI_activity_period,' period ',...
                    mod_inc.ROI_activity_ROI_nm,' activity ']);
            end
            %% number of pairs solved (effort #1) (incentive, no error trials)
            if mod_inc.n_pairs_solved ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.n_pairs_solved
                    case 1
                        fprintf('Number of pairs solved per trial');
                end
                fprintf('/');
            end
            %% percentage of incongruent pairs (effort #2) (incentive, no error trials)
            if mod_inc.n_incong ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.n_incong
                    case 1
                        fprintf('Number of incongruent pairs solved');
                    case 2
                        fprintf('Percentage of incongruent pairs solved per trial');
                end
                fprintf('/');
            end
            %% error numbers (incentive, no error trials)
            % pointless for no-error case! be careful with contrasts
            %% model variables
            if mod_inc.mdl_n ~= 0
                %% ressource (incentive, no error trials)
                if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[3,4])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.E_pred
                        case 3
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                        case 4
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
                %% cost (incentive, no error trials)
                if mod_inc.mdl_cost ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.mdl_cost
                        case 1
                            fprintf(['model ',num2str(mod_inc.mdl_n),' cost based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_inc.mdl_n),' cost based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_inc.mdl_n),' cost (before perf/X scaling)']);
                    end
                    fprintf('/');
                end
                %% benefit (incentive, no error trials)
                if mod_inc.mdl_benefit ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.mdl_benefit
                        case 1
                            fprintf(['model ',num2str(mod_inc.mdl_n),' benefit based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_inc.mdl_n),' benefit based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_inc.mdl_n),' benefit (before perf/X scaling)']);
                        case 4
                            fprintf(['model ',num2str(mod_inc.mdl_n),' expected monetary payoff']);
                    end
                    fprintf('/');
                end
                %% EV (incentive, no error trials)
                if mod_inc.mdl_EV ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.mdl_EV
                        case 1
                            fprintf(['model ',num2str(mod_inc.mdl_n),' EV based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_inc.mdl_n),' EV based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% trial number (middle regressor) (incentive, no error trials)
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[5,6])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.trialN
                    case 5
                        fprintf('trialN');
                    case 6
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first pair (incentive, no error trials)
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.RT_fp
                    case 1
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% mean RT all but first pair (incentive, no error trials)
            if mod_inc.RT_mRT ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.RT_mRT
                    case 1
                        fprintf('mean(RT) (across pairs except 1st pair)');
                end
                fprintf('/');
            end
            %% luminance (last regressor) (incentive, no error trials)
            if mod_inc.lum ~= 0 && ismember(mod_inc.lum,[1,3,5,7])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.lum
                    case 1
                        fprintf('luminance');
                    case 3
                        fprintf('luminance zPerRun');
                    case 5
                        fprintf('luminance zAcrossRuns');
                    case 7
                        fprintf('luminance zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (last regressor) (incentive, no error trials)
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% incentive*trial number (incentive, no error trials)
            if mod_inc.inc_x_trialN ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.inc_x_trialN
                    case 1
                        fprintf('inc_x_trialN');
                end
                fprintf('/');
            end
            %% ressource (incentive, no error trials)
            if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[1,2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.E_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (incentive, no error trials)
            if mod_inc.perf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.perf
                    case 1
                        fprintf('performance (1-100)');
                    case 2
                        fprintf('performance (0-1)');
                    case 3
                        fprintf('performance (0-1) corrected');
                    case 4
                        fprintf(['performance (0-1) predicted by model ',num2str(mod_inc.mdl_n)]);
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            %% incentive onset + duration (incentive, error trials)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Incentive display (error trials): ');
            switch dur_inc
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset to effort scale');
                case 2
                    fprintf('boxcar from onset to feedback');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% Gain/Loss condition (incentive, error trials)
            if mod_inc.GL_cond ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.GL_cond
                    case 1
                        fprintf('Gain or Loss condition (1/0)');
                    case 2
                        fprintf('Gain or Loss condition (1/-1)');
                end
                fprintf('/');
            end
            %% Luminance (first regressor) (incentive, error trials)
            if mod_inc.lum ~= 0 && ismember(mod_inc.lum,[2,4,6,8])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.lum
                    case 2
                        fprintf('luminance');
                    case 4
                        fprintf('luminance zPerRun');
                    case 6
                        fprintf('luminance zAcrossRuns');
                    case 8
                        fprintf('luminance zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (first regressor) (incentive, error trials)
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first pair (incentive, error trials)
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[1])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.RT_fp
                    case 2
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% model variables
            if mod_inc.mdl_n ~= 0
                %% ressource (incentive, error trials)
                if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[5,6])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.E_pred
                        case 5
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                        case 6
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% piece/bill condition (incentive, error trials)
            if mod_inc.R_type ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.R_type
                    case 1
                        fprintf('reward type: piece of money or bill (0/1)');
                    otherwise
                        error('not ready');
                end
                fprintf('/');
            end
            %% incentive (incentive, error trials)
            if mod_inc.inc ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.inc
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivational value');
                    case 4
                        fprintf('incentive motivational value');
                end
                fprintf('/');
            end
            %% incentive bis (incentive, error trials)
            if mod_inc.inc_bis ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.inc_bis
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivational value');
                    case 4
                        fprintf('incentive motivational value');
                end
                fprintf('/');
            end
            %% confidence proxy (incentive, error trials)
            if mod_inc.conf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.conf
                    case 1
                        fprintf('confidence (incentive rank centered)');
                    case 2
                        fprintf('confidence (incentive money centered)');
                    otherwise
                        error(['mod_inc.conf = ',num2str(mod_inc.conf),' not ready yet']);
                end
                fprintf('/');
            end
            %% ROI activity as regressor (incentive, error trials)
            if mod_inc.ROI_activity_yn == 1
                n_prm_stroop = n_prm_stroop + 1;
                fprintf(['GLM',mod_inc.ROI_activity_GLM,...
                    ' ',mod_inc.ROI_activity_period,' period ',...
                    mod_inc.ROI_activity_ROI_nm,' activity ']);
            end
            %% number of pairs solved (incentive, error trials)
            if mod_inc.n_pairs_solved ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.n_pairs_solved
                    case 1
                        fprintf('Number of pairs solved per trial');
                end
                fprintf('/');
            end
            %% percentage of incongruent pairs (effort #2) (incentive, error trials)
            if mod_inc.n_incong ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.n_incong
                    case 1
                        fprintf('Number of incongruent pairs solved');
                    case 2
                        fprintf('Percentage of incongruent pairs solved per trial');
                end
                fprintf('/');
            end
            %% error numbers (incentive, error trials)
            if mod_inc.n_errors ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.n_errors
                    case 1
                        error('n_errors = 1 pointless if o_inc = 3');
                    case 2
                        fprintf('Number of errors made per trial across pairs'); % includes errors repeated on the same pair
                    case 3
                        fprintf('Percentage of errors made per trial');
                end
                fprintf('/');
            end
            %% model variables
            if mod_inc.mdl_n ~= 0
                %% ressource (incentive, error trials)
                if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[3,4])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.E_pred
                        case 3
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                        case 4
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
                %% cost (incentive, error trials)
                if mod_inc.mdl_cost ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.mdl_cost
                        case 1
                            fprintf(['model ',num2str(mod_inc.mdl_n),' cost based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_inc.mdl_n),' cost based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_inc.mdl_n),' cost (before perf/X scaling)']);
                    end
                    fprintf('/');
                end
                %% benefit (incentive, error trials)
                if mod_inc.mdl_benefit ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.mdl_benefit
                        case 1
                            fprintf(['model ',num2str(mod_inc.mdl_n),' benefit based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_inc.mdl_n),' benefit based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_inc.mdl_n),' benefit (before perf/X scaling)']);
                        case 4
                            fprintf(['model ',num2str(mod_inc.mdl_n),' expected monetary payoff']);
                    end
                    fprintf('/');
                end
                %% EV (incentive, error trials)
                if mod_inc.mdl_EV ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_inc.mdl_EV
                        case 1
                            fprintf(['model ',num2str(mod_inc.mdl_n),' EV based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_inc.mdl_n),' EV based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% trial number (middle regressor) (incentive, error trials)
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[5,6])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.trialN
                    case 5
                        fprintf('trialN');
                    case 6
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first pair (incentive, error trials)
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.RT_fp
                    case 1
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% mean RT all but first pair (incentive, error trials)
            if mod_inc.RT_mRT ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.RT_mRT
                    case 1
                        fprintf('mean(RT) (across pairs except 1st pair)');
                end
                fprintf('/');
            end
            %% luminance (last regressor) (incentive, error trials)
            if mod_inc.lum ~= 0 && ismember(mod_inc.lum,[1,3,5,7])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.lum
                    case 1
                        fprintf('luminance');
                    case 3
                        fprintf('luminance zPerRun');
                    case 5
                        fprintf('luminance zAcrossRuns');
                    case 7
                        fprintf('luminance zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (last regressor) (incentive, error trials)
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% incentive*trial number (incentive, error trials)
            if mod_inc.inc_x_trialN ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.inc_x_trialN
                    case 1
                        fprintf('inc_x_trialN');
                end
                fprintf('/');
            end
            %% ressource (incentive, error trials)
            if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[1,2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.E_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (incentive, error trials)
            if mod_inc.perf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_inc.perf
                    case 1
                        fprintf('performance (1-100)');
                    case 2
                        fprintf('performance (0-1)');
                    case 3
                        fprintf('performance (0-1) corrected');
                    case 4
                        fprintf(['performance (0-1) predicted by model ',num2str(mod_inc.mdl_n)]);
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
    end % o_inc
end % o_inc#0

%% effort scale display
o_dispE = stroop.o_dispE;
dur_dispE = stroop.dur_dispE;
mod_dispE = stroop.mod_dispE;
if o_dispE ~= 0
    switch o_dispE
        case 1 % trials pooled
            %% effort scale display onset + duration
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Effort Scale display (all): ');
            switch dur_dispE
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until feedback');
                case 2
                    fprintf('boxcar from onset until RT');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% Gain/Loss condition (effort scale, all trials)
            if mod_dispE.GL_cond ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.GL_cond
                    case 1
                        fprintf('Gain or Loss condition (1/0)');
                    case 2
                        fprintf('Gain or Loss condition (1/-1)');
                end
                fprintf('/');
            end
            %% luminance (first regressor) (effort scale, all trials)
            if mod_dispE.lum ~= 0 && ismember(mod_dispE.lum,[2,4,6,8])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.lum
                    case 2
                        fprintf('luminance (incentive without effort scale)');
                    case 4
                        fprintf('luminance (incentive without effort scale) zPerRun');
                    case 6
                        fprintf('luminance (incentive without effort scale) zAcrossRuns');
                    case 8
                        fprintf('luminance (incentive without effort scale) zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (first regressor) (all trials, effort scale)
            if mod_dispE.trialN ~= 0 && ismember(mod_dispE.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first pair (effort scale, all trials)
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[1])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.RT_fp
                    case 2
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% model variables
            if mod_dispE.mdl_n ~= 0
                %% ressource (effort scale, all trials)
                if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[5,6])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.E_pred
                        case 5
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                        case 6
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% piece/bill condition (effort scale, all trials)
            if mod_dispE.R_type ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.R_type
                    case 1
                        fprintf('reward type: piece of money or bill (0/1)');
                    otherwise
                        error('not ready');
                end
                fprintf('/');
            end
            %% incentive (effort scale, all trials)
            if mod_dispE.inc ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.inc
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivationnal value');
                    case 4
                        fprintf('incentive motivationnal value');
                end
                fprintf('/');
            end
            %% incentive bis (all trials, effort scale)
            if mod_dispE.inc_bis ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.inc_bis
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivational value');
                    case 4
                        fprintf('incentive motivational value');
                end
                fprintf('/');
            end
            %% confidence proxy (effort scale, all trials)
            if mod_dispE.conf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.conf
                    case 1
                        fprintf('confidence (incentive rank centered)');
                    case 2
                        fprintf('confidence (incentive money centered)');
                    otherwise
                        error(['mod_dispE.conf = ',num2str(mod_dispE.conf),' not ready yet']);
                end
                fprintf('/');
            end
            %% ROI activity as regressor (effort scale, all trials)
            if mod_dispE.ROI_activity_yn == 1
                n_prm_stroop = n_prm_stroop + 1;
                fprintf(['GLM',mod_dispE.ROI_activity_GLM,...
                    ' ',mod_dispE.ROI_activity_period,' period ',...
                    mod_dispE.ROI_activity_ROI_nm,' activity ']);
            end
            %% number of pairs solved (effort #1) (effort scale, all trials)
            if mod_dispE.n_pairs_solved ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_pairs_solved
                    case 1
                        fprintf('Number of pairs solved per trial');
                end
                fprintf('/');
            end
            %% percentage of incongruent pairs (effort #2) (effort scale, all trials)
            if mod_dispE.n_incong ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_incong
                    case 1
                        fprintf('Number of incongruent pairs solved');
                    case 2
                        fprintf('Percentage of incongruent pairs solved per trial');
                end
                fprintf('/');
            end
            %% number distance (effort #3) (effort scale, all trials)
            if mod_dispE.n_dist ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_dist
                    case 1
                        fprintf('average of the number distance of the pairs solved per trial');
                    case 2
                        fprintf('average of 1/(number distance of the pairs solved) per trial');
                end
                fprintf('/');
            end
            %% error numbers (effort scale, all trials)
            if mod_dispE.n_errors ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_errors
                    case 1
                        fprintf('Trial with at least one error (1) or no error at all (0)');
                    case 2
                        fprintf('Number of errors made per trial across pairs'); % includes errors repeated on the same pair
                    case 3
                        fprintf('Percentage of errors made per trial');
                end
                fprintf('/');
            end
            %% model variables
            if mod_dispE.mdl_n ~= 0
                %% ressource (effort scale, all trials)
                if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[3,4])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.E_pred
                        case 3
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                        case 4
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
                %% cost (all trials, effort scale)
                if mod_dispE.mdl_cost ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.mdl_cost
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' cost based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' cost based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' cost (before perf/X scaling)']);
                    end
                    fprintf('/');
                end
                %% benefit (all trials, effort scale)
                if mod_dispE.mdl_benefit ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.mdl_benefit
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' benefit based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' benefit based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' benefit (before perf/X scaling)']);
                        case 4
                            fprintf(['model ',num2str(mod_inc.mdl_n),' expected monetary payoff']);
                    end
                    fprintf('/');
                end
                %% EV (all trials, effort scale)
                if mod_dispE.mdl_EV ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.mdl_EV
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' EV based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' EV based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% RT first pair (effort scale, all trials)
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.RT_fp
                    case 1
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% mean RT across pairs (except first pair) (effort scale, all trials)
            if mod_dispE.RT_mRT ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.RT_mRT
                    case 1
                        fprintf('mean RT across pairs (except first pair)');
                    case 2
                        fprintf('mean RT across pairs (except first pair AND errors)');
                end
                fprintf('/');
            end
            %% luminance (last regressor) (effort scale, all trials)
            if mod_dispE.lum ~= 0 && ismember(mod_dispE.lum,[1,3,5,7])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.lum
                    case 1
                        fprintf('luminance (incentive without effort scale)');
                    case 3
                        fprintf('luminance (incentive without effort scale) zPerRun');
                    case 5
                        fprintf('luminance (incentive without effort scale) zAcrossRuns');
                    case 7
                        fprintf('luminance (incentive without effort scale) zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (last regressor) (all trials, effort scale)
            if mod_dispE.trialN ~= 0 && ismember(mod_dispE.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% incentive*trial number (effort scale, all trials)
            if mod_dispE.inc_x_trialN ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.inc_x_trialN
                    case 1
                        fprintf('inc_x_trialN');
                end
                fprintf('/');
            end
            %% ressource (effort scale, all trials)
            if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[1,2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.E_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (effort scale, all trials)
            if mod_dispE.perf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.perf
                    case 1
                        fprintf('performance (1-100)');
                    case 2
                        fprintf('performance (0-1)');
                    case 3
                        fprintf('performance (0-1) corrected');
                    case 4
                        fprintf(['performance (0-1) predicted by model ',num2str(mod_dispE.mdl_n)]);
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
        case 2 % split by gain/loss trials type
            %% Effort scale (to gain trials)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Effort Scale display (to gain trials): ');
            switch dur_dispE
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until feedback');
                case 2
                    fprintf('boxcar from onset until RT');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% luminance (first regressor) (effort scale, to gain trials)
            if mod_dispE.lum ~= 0 && ismember(mod_dispE.lum,[2,4,6,8])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.lum
                    case 2
                        fprintf('luminance (incentive without effort scale)');
                    case 4
                        fprintf('luminance (incentive without effort scale) zPerRun');
                    case 6
                        fprintf('luminance (incentive without effort scale) zAcrossRuns');
                    case 8
                        fprintf('luminance (incentive without effort scale) zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (first regressor) (to gain trials, effort scale)
            if mod_dispE.trialN ~= 0 && ismember(mod_dispE.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first pair (effort scale, to gain trials)
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[1])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.RT_fp
                    case 2
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% model variables
            if mod_dispE.mdl_n ~= 0
                %% ressource (effort scale, to gain trials)
                if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[5,6])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.E_pred
                        case 5
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                        case 6
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% piece/bill condition (effort scale, to gain trials)
            if mod_dispE.R_type ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.R_type
                    case 1
                        fprintf('reward type: piece of money or bill (0/1)');
                    otherwise
                        error('not ready');
                end
                fprintf('/');
            end
            %% incentive (effort scale, to gain trials)
            if mod_dispE.inc ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.inc
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivationnal value');
                    case 4
                        fprintf('incentive motivationnal value');
                end
                fprintf('/');
            end
            %% incentive bis (to gain trials, effort scale)
            if mod_dispE.inc_bis ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.inc_bis
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivationnal value');
                    case 4
                        fprintf('incentive motivationnal value');
                end
                fprintf('/');
            end
            %% confidence proxy (effort scale, to gain trials)
            if mod_dispE.conf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.conf
                    case 1
                        fprintf('confidence (incentive rank centered)');
                    case 2
                        fprintf('confidence (incentive money centered)');
                    otherwise
                        error(['mod_dispE.conf = ',num2str(mod_dispE.conf),' not ready yet']);
                end
                fprintf('/');
            end
            %% ROI activity as regressor (effort scale, to gain trials)
            if mod_dispE.ROI_activity_yn == 1
                n_prm_stroop = n_prm_stroop + 1;
                fprintf(['GLM',mod_dispE.ROI_activity_GLM,...
                    ' ',mod_dispE.ROI_activity_period,' period ',...
                    mod_dispE.ROI_activity_ROI_nm,' activity ']);
            end
            %% number of pairs solved (effort #1) (effort scale, to gain trials)
            if mod_dispE.n_pairs_solved ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_pairs_solved
                    case 1
                        fprintf('Number of pairs solved per trial');
                end
                fprintf('/');
            end
            %% percentage of incongruent pairs (effort #2) (effort scale, to gain trials)
            if mod_dispE.n_incong ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_incong
                    case 1
                        fprintf('Number of incongruent pairs solved');
                    case 2
                        fprintf('Percentage of incongruent pairs solved per trial');
                end
                fprintf('/');
            end
            %% number distance (effort #3) (effort scale, to gain trials)
            if mod_dispE.n_dist ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_dist
                    case 1
                        fprintf('average of the number distance of the pairs solved per trial');
                    case 2
                        fprintf('average of 1/(number distance of the pairs solved) per trial');
                end
                fprintf('/');
            end
            %% error numbers (effort scale, to gain trials)
            if mod_dispE.n_errors ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_errors
                    case 1
                        fprintf('Trial with at least one error (1) or no error at all (0)');
                    case 2
                        fprintf('Number of errors made per trial across pairs'); % includes errors repeated on the same pair
                    case 3
                        fprintf('Percentage of errors made per trial');
                end
                fprintf('/');
            end
            %% model variables
            if mod_dispE.mdl_n ~= 0
                %% ressource (effort scale, to gain trials)
                if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[3,4])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.E_pred
                        case 3
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                        case 4
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
                %% cost (to gain trials, effort scale)
                if mod_dispE.mdl_cost ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.mdl_cost
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' cost based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' cost based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' cost (before perf/X scaling)']);
                    end
                    fprintf('/');
                end
                %% benefit (to gain trials, effort scale)
                if mod_dispE.mdl_benefit ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.mdl_benefit
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' benefit based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' benefit based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' benefit (before perf/X scaling)']);
                        case 4
                            fprintf(['model ',num2str(mod_inc.mdl_n),' expected monetary payoff']);
                    end
                    fprintf('/');
                end
                %% EV (to gain trials, effort scale)
                if mod_dispE.mdl_EV ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.mdl_EV
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' EV based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' EV based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% RT first pair (effort scale, to gain trials)
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.RT_fp
                    case 1
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% mean RT across pairs (except first pair) (effort scale, to gain trials)
            if mod_dispE.RT_mRT ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.RT_mRT
                    case 1
                        fprintf('mean RT across pairs (except first pair)');
                    case 2
                        fprintf('mean RT across pairs (except first pair AND errors)');
                end
                fprintf('/');
            end
            %% luminance (last regressor) (effort scale, to gain trials)
            if mod_dispE.lum ~= 0 && ismember(mod_dispE.lum,[1,3,5,7])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.lum
                    case 1
                        fprintf('luminance (incentive without effort scale)');
                    case 3
                        fprintf('luminance (incentive without effort scale) zPerRun');
                    case 5
                        fprintf('luminance (incentive without effort scale) zAcrossRuns');
                    case 7
                        fprintf('luminance (incentive without effort scale) zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (last regressor) (to gain trials, effort scale)
            if mod_dispE.trialN ~= 0 && ismember(mod_dispE.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% incentive*trial number (effort scale, to gain trials)
            if mod_dispE.inc_x_trialN ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.inc_x_trialN
                    case 1
                        fprintf('inc_x_trialN');
                end
                fprintf('/');
            end
            %% ressource (effort scale, to gain trials)
            if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[1,2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.E_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (effort scale, to gain trials)
            if mod_dispE.perf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.perf
                    case 1
                        fprintf('performance (1-100)');
                    case 2
                        fprintf('performance (0-1)');
                    case 3
                        fprintf('performance (0-1) corrected');
                    case 4
                        fprintf(['performance (0-1) predicted by model ',num2str(mod_dispE.mdl_n)]);
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            %% Effort scale onset + duration (to lose trials)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Effort Scale display (to lose trials): ');
            switch dur_dispE
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until feedback');
                case 2
                    fprintf('boxcar from onset until RT');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% luminance (first regressor) (effort scale, to lose trials)
            if mod_dispE.lum ~= 0 && ismember(mod_dispE.lum,[2,4,6,8])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.lum
                    case 2
                        fprintf('luminance (incentive without effort scale)');
                    case 4
                        fprintf('luminance (incentive without effort scale) zPerRun');
                    case 6
                        fprintf('luminance (incentive without effort scale) zAcrossRuns');
                    case 8
                        fprintf('luminance (incentive without effort scale) zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (first regressor) (to lose trials, effort scale)
            if mod_dispE.trialN ~= 0 && ismember(mod_dispE.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first pair (effort scale, to lose trials)
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[1])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.RT_fp
                    case 2
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% model variables
            if mod_dispE.mdl_n ~= 0
                %% ressource (effort scale, to lose trials)
                if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[5,6])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.E_pred
                        case 5
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                        case 6
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% piece/bill condition (effort scale, to lose trials)
            if mod_dispE.R_type ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.R_type
                    case 1
                        fprintf('reward type: piece of money or bill (0/1)');
                    otherwise
                        error('not ready');
                end
                fprintf('/');
            end
            %% incentive (effort scale, to lose trials)
            if mod_dispE.inc ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.inc
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivationnal value');
                    case 4
                        fprintf('incentive motivationnal value');
                end
                fprintf('/');
            end
            %% incentive bis (to lose trials, effort scale)
            if mod_dispE.inc_bis ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.inc_bis
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivationnal value');
                    case 4
                        fprintf('incentive motivationnal value');
                end
                fprintf('/');
            end
            %% confidence proxy (effort scale, to lose trials)
            if mod_dispE.conf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.conf
                    case 1
                        fprintf('confidence (incentive rank centered)');
                    case 2
                        fprintf('confidence (incentive money centered)');
                    otherwise
                        error(['mod_dispE.conf = ',num2str(mod_dispE.conf),' not ready yet']);
                end
                fprintf('/');
            end
            %% ROI activity as regressor (effort scale, to lose trials)
            if mod_dispE.ROI_activity_yn == 1
                n_prm_stroop = n_prm_stroop + 1;
                fprintf(['GLM',mod_dispE.ROI_activity_GLM,...
                    ' ',mod_dispE.ROI_activity_period,' period ',...
                    mod_dispE.ROI_activity_ROI_nm,' activity ']);
            end
            %% number of pairs solved (effort #1) (effort scale, to lose trials)
            if mod_dispE.n_pairs_solved ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_pairs_solved
                    case 1
                        fprintf('Number of pairs solved per trial');
                end
                fprintf('/');
            end
            %% percentage of incongruent pairs (effort #2) (effort scale, to lose trials)
            if mod_dispE.n_incong ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_incong
                    case 1
                        fprintf('Number of incongruent pairs solved');
                    case 2
                        fprintf('Percentage of incongruent pairs solved per trial');
                end
                fprintf('/');
            end
            %% number distance (effort #3) (effort scale, to lose trials)
            if mod_dispE.n_dist ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_dist
                    case 1
                        fprintf('average of the number distance of the pairs solved per trial');
                    case 2
                        fprintf('average of 1/(number distance of the pairs solved) per trial');
                end
                fprintf('/');
            end
            %% error numbers (effort scale, to lose trials)
            if mod_dispE.n_errors ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_errors
                    case 1
                        fprintf('Trial with at least one error (1) or no error at all (0)');
                    case 2
                        fprintf('Number of errors made per trial across pairs'); % includes errors repeated on the same pair
                    case 3
                        fprintf('Percentage of errors made per trial');
                end
                fprintf('/');
            end
            %% model variables
            if mod_dispE.mdl_n ~= 0
                %% ressource (effort scale, to lose trials)
                if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[3,4])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.E_pred
                        case 3
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                        case 4
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
                %% cost (to lose trials, effort scale)
                if mod_dispE.mdl_cost ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.mdl_cost
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' cost based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' cost based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' cost (before perf/X scaling)']);
                    end
                    fprintf('/');
                end
                %% benefit (to lose trials, effort scale)
                if mod_dispE.mdl_benefit ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.mdl_benefit
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' benefit based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' benefit based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' benefit (before perf/X scaling)']);
                        case 4
                            fprintf(['model ',num2str(mod_inc.mdl_n),' expected monetary payoff']);
                    end
                    fprintf('/');
                end
                %% EV (to lose trials, effort scale)
                if mod_dispE.mdl_EV ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.mdl_EV
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' EV based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' EV based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% RT first pair (effort scale, to lose trials)
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.RT_fp
                    case 1
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% mean RT across pairs (except first pair) (effort scale, to lose trials)
            if mod_dispE.RT_mRT ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.RT_mRT
                    case 1
                        fprintf('mean RT across pairs (except first pair)');
                    case 2
                        fprintf('mean RT across pairs (except first pair AND errors)');
                end
                fprintf('/');
            end
            %% luminance (last regressor) (effort scale, to lose trials)
            if mod_dispE.lum ~= 0 && ismember(mod_dispE.lum,[1,3,5,7])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.lum
                    case 1
                        fprintf('luminance (incentive without effort scale)');
                    case 3
                        fprintf('luminance (incentive without effort scale) zPerRun');
                    case 5
                        fprintf('luminance (incentive without effort scale) zAcrossRuns');
                    case 7
                        fprintf('luminance (incentive without effort scale) zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (last regressor) (to lose trials, effort scale)
            if mod_dispE.trialN ~= 0 && ismember(mod_dispE.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% incentive*trial number (effort scale, to lose trials)
            if mod_dispE.inc_x_trialN ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.inc_x_trialN
                    case 1
                        fprintf('inc_x_trialN');
                end
                fprintf('/');
            end
            %% ressource (effort scale, to lose trials)
            if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[1,2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.E_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (effort scale, to lose trials)
            if mod_dispE.perf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.perf
                    case 1
                        fprintf('performance (1-100)');
                    case 2
                        fprintf('performance (0-1)');
                    case 3
                        fprintf('performance (0-1) corrected');
                    case 4
                        fprintf(['performance (0-1) predicted by model ',num2str(mod_dispE.mdl_n)]);
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
        case 3 % split by no-error/error trials type
            %% Effort scale (no error trials)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Effort Scale display (no error trials): ');
            switch dur_dispE
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until feedback');
                case 2
                    fprintf('boxcar from onset until RT');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% Gain/Loss condition (effort scale, no error trials)
            if mod_dispE.GL_cond ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.GL_cond
                    case 1
                        fprintf('Gain or Loss condition (1/0)');
                    case 2
                        fprintf('Gain or Loss condition (1/-1)');
                end
                fprintf('/');
            end
            %% luminance (first regressor) (effort scale, no error trials)
            if mod_dispE.lum ~= 0 && ismember(mod_dispE.lum,[2,4,6,8])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.lum
                    case 2
                        fprintf('luminance (incentive without effort scale)');
                    case 4
                        fprintf('luminance (incentive without effort scale) zPerRun');
                    case 6
                        fprintf('luminance (incentive without effort scale) zAcrossRuns');
                    case 8
                        fprintf('luminance (incentive without effort scale) zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (first regressor) (no error trials, effort scale)
            if mod_dispE.trialN ~= 0 && ismember(mod_dispE.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first pair (effort scale, no error trials)
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[1])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.RT_fp
                    case 2
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% model variables
            if mod_dispE.mdl_n ~= 0
                %% ressource (effort scale, no error trials)
                if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[5,6])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.E_pred
                        case 5
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                        case 6
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% piece/bill condition (effort scale, no error trials)
            if mod_dispE.R_type ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.R_type
                    case 1
                        fprintf('reward type: piece of money or bill (0/1)');
                    otherwise
                        error('not ready');
                end
                fprintf('/');
            end
            %% incentive (effort scale, no error trials)
            if mod_dispE.inc ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.inc
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivationnal value');
                    case 4
                        fprintf('incentive motivationnal value');
                end
                fprintf('/');
            end
            %% incentive bis (no error trials, effort scale)
            if mod_dispE.inc_bis ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.inc_bis
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivationnal value');
                    case 4
                        fprintf('incentive motivationnal value');
                end
                fprintf('/');
            end
            %% confidence proxy (effort scale, no error trials)
            if mod_dispE.conf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.conf
                    case 1
                        fprintf('confidence (incentive rank centered)');
                    case 2
                        fprintf('confidence (incentive money centered)');
                    otherwise
                        error(['mod_dispE.conf = ',num2str(mod_dispE.conf),' not ready yet']);
                end
                fprintf('/');
            end
            %% ROI activity as regressor (effort scale, no error trials)
            if mod_dispE.ROI_activity_yn == 1
                n_prm_stroop = n_prm_stroop + 1;
                fprintf(['GLM',mod_dispE.ROI_activity_GLM,...
                    ' ',mod_dispE.ROI_activity_period,' period ',...
                    mod_dispE.ROI_activity_ROI_nm,' activity ']);
            end
            %% number of pairs solved (effort #1) (effort scale, no error trials)
            if mod_dispE.n_pairs_solved ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_pairs_solved
                    case 1
                        fprintf('Number of pairs solved per trial');
                end
                fprintf('/');
            end
            %% percentage of incongruent pairs (effort #2) (effort scale, no error trials)
            if mod_dispE.n_incong ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_incong
                    case 1
                        fprintf('Number of incongruent pairs solved');
                    case 2
                        fprintf('Percentage of incongruent pairs solved per trial');
                end
                fprintf('/');
            end
            %% number distance (effort #3) (effort scale, no error trials)
            if mod_dispE.n_dist ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_dist
                    case 1
                        fprintf('average of the number distance of the pairs solved per trial');
                    case 2
                        fprintf('average of 1/(number distance of the pairs solved) per trial');
                end
                fprintf('/');
            end
            %% error numbers (effort scale, no error trials)
            % pointless for no-error case! be careful with contrasts
            %% model variables
            if mod_dispE.mdl_n ~= 0
                %% ressource (effort scale, no error trials)
                if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[3,4])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.E_pred
                        case 3
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                        case 4
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
                %% cost (no error trials, effort scale)
                if mod_dispE.mdl_cost ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.mdl_cost
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' cost based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' cost based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' cost (before perf/X scaling)']);
                    end
                    fprintf('/');
                end
                %% benefit (no error trials, effort scale)
                if mod_dispE.mdl_benefit ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.mdl_benefit
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' benefit based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' benefit based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' benefit (before perf/X scaling)']);
                        case 4
                            fprintf(['model ',num2str(mod_inc.mdl_n),' expected monetary payoff']);
                    end
                    fprintf('/');
                end
                %% EV (no error trials, effort scale)
                if mod_dispE.mdl_EV ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.mdl_EV
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' EV based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' EV based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% RT first pair (effort scale, no error trials)
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.RT_fp
                    case 1
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% mean RT across pairs (except first pair) (effort scale, no error trials)
            if mod_dispE.RT_mRT ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.RT_mRT
                    case 1
                        fprintf('mean RT across pairs (except first pair)');
                    case 2
                        fprintf('mean RT across pairs (except first pair AND errors)');
                end
                fprintf('/');
            end
            %% luminance (last regressor) (effort scale, no error trials)
            if mod_dispE.lum ~= 0 && ismember(mod_dispE.lum,[1,3,5,7])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.lum
                    case 1
                        fprintf('luminance (incentive without effort scale)');
                    case 3
                        fprintf('luminance (incentive without effort scale) zPerRun');
                    case 5
                        fprintf('luminance (incentive without effort scale) zAcrossRuns');
                    case 7
                        fprintf('luminance (incentive without effort scale) zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (last regressor) (no error trials, effort scale)
            if mod_dispE.trialN ~= 0 && ismember(mod_dispE.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% incentive*trial number (effort scale, no error trials)
            if mod_dispE.inc_x_trialN ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.inc_x_trialN
                    case 1
                        fprintf('inc_x_trialN');
                end
                fprintf('/');
            end
            %% ressource (effort scale, no error trials)
            if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[1,2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.E_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' EV based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' EV based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (effort scale, no error trials)
            if mod_dispE.perf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.perf
                    case 1
                        fprintf('performance (1-100)');
                    case 2
                        fprintf('performance (0-1)');
                    case 3
                        fprintf('performance (0-1) corrected');
                    case 4
                        fprintf(['performance (0-1) predicted by model ',num2str(mod_dispE.mdl_n)]);
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            %% Effort scale onset + duration (error trials)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Effort Scale display (error trials): ');
            switch dur_dispE
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until feedback');
                case 2
                    fprintf('boxcar from onset until RT');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% Gain/Loss condition (effort scale, error trials)
            if mod_dispE.GL_cond ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.GL_cond
                    case 1
                        fprintf('Gain or Loss condition (1/0)');
                    case 2
                        fprintf('Gain or Loss condition (1/-1)');
                end
                fprintf('/');
            end
            %% luminance (first regressor) (effort scale, error trials)
            if mod_dispE.lum ~= 0 && ismember(mod_dispE.lum,[2,4,6,8])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.lum
                    case 2
                        fprintf('luminance (incentive without effort scale)');
                    case 4
                        fprintf('luminance (incentive without effort scale) zPerRun');
                    case 6
                        fprintf('luminance (incentive without effort scale) zAcrossRuns');
                    case 8
                        fprintf('luminance (incentive without effort scale) zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (first regressor) (error trials, effort scale)
            if mod_dispE.trialN ~= 0 && ismember(mod_dispE.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first pair (effort scale, error trials)
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[1])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.RT_fp
                    case 2
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% model variables
            if mod_dispE.mdl_n ~= 0
                %% ressource (effort scale, error trials)
                if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[5,6])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.E_pred
                        case 5
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                        case 6
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% piece/bill condition (effort scale, error trials)
            if mod_dispE.R_type ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.R_type
                    case 1
                        fprintf('reward type: piece of money or bill (0/1)');
                    otherwise
                        error('not ready');
                end
                fprintf('/');
            end
            %% incentive (effort scale, error trials)
            if mod_dispE.inc ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.inc
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivationnal value');
                    case 4
                        fprintf('incentive motivationnal value');
                end
                fprintf('/');
            end
            %% incentive bis (error trials, effort scale)
            if mod_dispE.inc_bis ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.inc_bis
                    case 1
                        fprintf('incentive rank nominal value');
                    case 2
                        fprintf('incentive nominal value');
                    case 3
                        fprintf('incentive rank motivationnal value');
                    case 4
                        fprintf('incentive motivationnal value');
                end
                fprintf('/');
            end
            %% confidence proxy (effort scale, error trials)
            if mod_dispE.conf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.conf
                    case 1
                        fprintf('confidence (incentive rank centered)');
                    case 2
                        fprintf('confidence (incentive money centered)');
                    otherwise
                        error(['mod_dispE.conf = ',num2str(mod_dispE.conf),' not ready yet']);
                end
                fprintf('/');
            end
            %% ROI activity as regressor (effort scale, error trials)
            if mod_dispE.ROI_activity_yn == 1
                n_prm_stroop = n_prm_stroop + 1;
                fprintf(['GLM',mod_dispE.ROI_activity_GLM,...
                    ' ',mod_dispE.ROI_activity_period,' period ',...
                    mod_dispE.ROI_activity_ROI_nm,' activity ']);
            end
            %% number of pairs solved (effort #1) (effort scale, error trials)
            if mod_dispE.n_pairs_solved ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_pairs_solved
                    case 1
                        fprintf('Number of pairs solved per trial');
                end
                fprintf('/');
            end
            %% percentage of incongruent pairs (effort #2) (effort scale, error trials)
            if mod_dispE.n_incong ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_incong
                    case 1
                        fprintf('Number of incongruent pairs solved');
                    case 2
                        fprintf('Percentage of incongruent pairs solved per trial');
                end
                fprintf('/');
            end
            %% number distance (effort #3) (effort scale, error trials)
            if mod_dispE.n_dist ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_dist
                    case 1
                        fprintf('average of the number distance of the pairs solved per trial');
                    case 2
                        fprintf('average of 1/(number distance of the pairs solved) per trial');
                end
                fprintf('/');
            end
            %% error numbers (effort scale, error trials)
            if mod_dispE.n_errors ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.n_errors
                    case 1
                        error('n_errors = 1 incompatible with o_dispE = 3');
                    case 2
                        fprintf('Number of errors made per trial across pairs'); % includes errors repeated on the same pair
                    case 3
                        fprintf('Percentage of errors made per trial');
                end
                fprintf('/');
            end
            %% model variables
            if mod_dispE.mdl_n ~= 0
                %% ressource (effort scale, error trials)
                if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[3,4])
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.E_pred
                        case 3
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                        case 4
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
                %% cost (error trials, effort scale)
                if mod_dispE.mdl_cost ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.mdl_cost
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' cost based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' cost based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' cost (before perf/X scaling)']);
                    end
                    fprintf('/');
                end
                %% benefit (error trials, effort scale)
                if mod_dispE.mdl_benefit ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.mdl_benefit
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' benefit based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' benefit based on predicted perf']);
                        case 3
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' benefit (before perf/X scaling)']);
                        case 4
                            fprintf(['model ',num2str(mod_inc.mdl_n),' expected monetary payoff']);
                    end
                    fprintf('/');
                end
                %% EV (error trials, effort scale)
                if mod_dispE.mdl_EV ~= 0
                    n_prm_stroop = n_prm_stroop + 1;
                    switch mod_dispE.mdl_EV
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' EV based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' EV based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% RT first pair (effort scale, error trials)
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.RT_fp
                    case 1
                        fprintf('RT first pair');
                end
                fprintf('/');
            end
            %% mean RT across pairs (except first pair) (effort scale, error trials)
            if mod_dispE.RT_mRT ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.RT_mRT
                    case 1
                        fprintf('mean RT across pairs (except first pair)');
                    case 2
                        fprintf('mean RT across pairs (except first pair AND errors)');
                end
                fprintf('/');
            end
            %% luminance (last regressor) (effort scale, error trials)
            if mod_dispE.lum ~= 0 && ismember(mod_dispE.lum,[1,3,5,7])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.lum
                    case 1
                        fprintf('luminance (incentive without effort scale)');
                    case 3
                        fprintf('luminance (incentive without effort scale) zPerRun');
                    case 5
                        fprintf('luminance (incentive without effort scale) zAcrossRuns');
                    case 7
                        fprintf('luminance (incentive without effort scale) zAcrossRuns and zPerRun');
                end
                fprintf('/');
            end
            %% trial number (last regressor) (error trials, effort scale)
            if mod_dispE.trialN ~= 0 && ismember(mod_dispE.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% incentive*trial number (effort scale, error trials)
            if mod_dispE.inc_x_trialN ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.inc_x_trialN
                    case 1
                        fprintf('inc_x_trialN');
                end
                fprintf('/');
            end
            %% ressource (effort scale, error trials)
            if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[1,2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.E_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' EV based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' EV based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (effort scale, error trials)
            if mod_dispE.perf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_dispE.perf
                    case 1
                        fprintf('performance (1-100)');
                    case 2
                        fprintf('performance (0-1)');
                    case 3
                        fprintf('performance (0-1) corrected');
                    case 4
                        fprintf(['performance (0-1) predicted by model ',num2str(mod_dispE.mdl_n)]);
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
    end % o_dispE
end % o_dispE#0

%% effort performance onset
o_perfE = stroop.o_perfE;
dur_perfE = stroop.dur_perfE;
mod_perfE = stroop.mod_perfE;
if o_perfE ~= 0
    switch o_perfE
        case 1 % trials pooled
            %% Effort performance onset + duration (all trials)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Effort Performance onset (all trials pooled): ');
            switch dur_perfE
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% Gain/Loss condition (all trials, effort perf)
            if mod_perfE.GL_cond ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.GL_cond
                    case 1
                        fprintf('Gain or Loss condition (1/0)');
                    case 2
                        fprintf('Gain or Loss condition (1/-1)');
                end
                fprintf('/');
            end
            %% trial number (first regressor) (all trials, effort perf)
            if mod_perfE.trialN ~= 0 && ismember(mod_perfE.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            
            
            %% number of pairs solved (effort #1) (effort performance, all trials)
            if mod_perfE.n_pairs_solved ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_pairs_solved
                    case 1
                        fprintf('Number of pairs solved per trial');
                end
                fprintf('/');
            end
            
            %% percentage of incongruent pairs (effort #2) (effort performance, all trials)
            if mod_perfE.n_incong ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_incong
                    case 1
                        fprintf('Number of incongruent pairs solved');
                    case 2
                        fprintf('Percentage of incongruent pairs solved per trial');
                end
                fprintf('/');
            end
            
            %% number distance (effort #3) (effort performance, all trials)
            if mod_perfE.n_dist ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_dist
                    case 1
                        fprintf('average of the number distance of the pairs solved per trial');
                    case 2
                        fprintf('average of 1/(number distance of the pairs solved) per trial');
                end
                fprintf('/');
            end
            
            %% error numbers (effort performance, all trials)
            if mod_perfE.n_errors ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_errors
                    case 1
                        fprintf('Trial with at least one error (1) or no error at all (0)');
                    case 2
                        fprintf('Number of errors made per trial across pairs'); % includes errors repeated on the same pair
                    case 3
                        fprintf('Percentage of errors made per trial');
                end
                fprintf('/');
            end
            
            %% trial number (last regressor) (all trials, effort perf)
            if mod_perfE.trialN ~= 0 && ismember(mod_perfE.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            
            
        case 2 % split by to gain/to lose trial type
            %% Effort performance onset + duration (to gain trials)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Effort Performance onset (to gain trials): ');
            switch dur_perfE
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% trial number (first regressor) (to gain trials, effort perf)
            if mod_perfE.trialN ~= 0 && ismember(mod_perfE.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% number of pairs solved (effort #1) (effort performance, to gain trials)
            if mod_perfE.n_pairs_solved ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_pairs_solved
                    case 1
                        fprintf('Number of pairs solved per trial');
                end
                fprintf('/');
            end
            
            %% percentage of incongruent pairs (effort #2) (effort performance, to gain trials)
            if mod_perfE.n_incong ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_incong
                    case 1
                        fprintf('Number of incongruent pairs solved');
                    case 2
                        fprintf('Percentage of incongruent pairs solved per trial');
                end
                fprintf('/');
            end
            
            %% number distance (effort #3) (effort performance, to gain trials)
            if mod_perfE.n_dist ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_dist
                    case 1
                        fprintf('average of the number distance of the pairs solved per trial');
                    case 2
                        fprintf('average of 1/(number distance of the pairs solved) per trial');
                end
                fprintf('/');
            end
            
            %% error numbers (effort performance, to gain trials)
            if mod_perfE.n_errors ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_errors
                    case 1
                        fprintf('Trial with at least one error (1) or no error at all (0)');
                    case 2
                        fprintf('Number of errors made per trial across pairs'); % includes errors repeated on the same pair
                    case 3
                        fprintf('Percentage of errors made per trial');
                end
                fprintf('/');
            end
            %% trial number (last regressor) (to gain trials, effort perf)
            if mod_perfE.trialN ~= 0 && ismember(mod_perfE.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            %% Effort performance onset + duration (to lose trials)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Effort Performance onset (to lose trials): ');
            switch dur_perfE
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% trial number (first regressor) (to lose trials, effort perf)
            if mod_perfE.trialN ~= 0 && ismember(mod_perfE.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% number of pairs solved (effort #1) (effort performance, to lose trials)
            if mod_perfE.n_pairs_solved ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_pairs_solved
                    case 1
                        fprintf('Number of pairs solved per trial');
                end
                fprintf('/');
            end
            
            %% percentage of incongruent pairs (effort #2) (effort performance, to lose trials)
            if mod_perfE.n_incong ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_incong
                    case 1
                        fprintf('Number of incongruent pairs solved');
                    case 2
                        fprintf('Percentage of incongruent pairs solved per trial');
                end
                fprintf('/');
            end
            
            %% number distance (effort #3) (effort performance, to lose trials)
            if mod_perfE.n_dist ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_dist
                    case 1
                        fprintf('average of the number distance of the pairs solved per trial');
                    case 2
                        fprintf('average of 1/(number distance of the pairs solved) per trial');
                end
                fprintf('/');
            end
            
            %% error numbers (effort performance, to lose trials)
            if mod_perfE.n_errors ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_errors
                    case 1
                        fprintf('Trial with at least one error (1) or no error at all (0)');
                    case 2
                        fprintf('Number of errors made per trial across pairs'); % includes errors repeated on the same pair
                    case 3
                        fprintf('Percentage of errors made per trial');
                end
                fprintf('/');
            end
            %% trial number (last regressor) (to lose trials, effort perf)
            if mod_perfE.trialN ~= 0 && ismember(mod_perfE.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            
            %%
            fprintf('\n');
            
        case 3 % split by no-error/error trial type
            %% Effort performance onset + duration (no error trials)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Effort Performance onset (no error trials): ');
            switch dur_perfE
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% Gain/Loss condition (no error trials, effort perf)
            if mod_perfE.GL_cond ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.GL_cond
                    case 1
                        fprintf('Gain or Loss condition (1/0)');
                    case 2
                        fprintf('Gain or Loss condition (1/-1)');
                end
                fprintf('/');
            end
            %% trial number (first regressor) (no error trials, effort perf)
            if mod_perfE.trialN ~= 0 && ismember(mod_perfE.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            
            %% number of pairs solved (effort #1) (effort performance, no error trials)
            if mod_perfE.n_pairs_solved ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_pairs_solved
                    case 1
                        fprintf('Number of pairs solved per trial');
                end
                fprintf('/');
            end
            
            %% percentage of incongruent pairs (effort #2) (effort performance, no error trials)
            if mod_perfE.n_incong ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_incong
                    case 1
                        fprintf('Number of incongruent pairs solved');
                    case 2
                        fprintf('Percentage of incongruent pairs solved per trial');
                end
                fprintf('/');
            end
            
            %% number distance (effort #3) (effort performance, no error trials)
            if mod_perfE.n_dist ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_dist
                    case 1
                        fprintf('average of the number distance of the pairs solved per trial');
                    case 2
                        fprintf('average of 1/(number distance of the pairs solved) per trial');
                end
                fprintf('/');
            end
            
            %% error numbers (effort scale, no error trials)
            % pointless for no-error case! be careful with contrasts
            %% trial number (last regressor) (no error trials, effort perf)
            if mod_perfE.trialN ~= 0 && ismember(mod_perfE.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            
            %%
            fprintf('\n');
            
            %% Effort performance onset + duration (error trials)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Effort Performance onset (error trials): ');
            switch dur_perfE
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset until feedback');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% Gain/Loss condition (error trials, effort perf)
            if mod_perfE.GL_cond ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.GL_cond
                    case 1
                        fprintf('Gain or Loss condition (1/0)');
                    case 2
                        fprintf('Gain or Loss condition (1/-1)');
                end
                fprintf('/');
            end
            %% trial number (first regressor) (error trials, effort perf)
            if mod_perfE.trialN ~= 0 && ismember(mod_perfE.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            
            %% number of pairs solved (effort #1) (effort performance, error trials)
            if mod_perfE.n_pairs_solved ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_pairs_solved
                    case 1
                        fprintf('Number of pairs solved per trial');
                end
                fprintf('/');
            end
            
            %% percentage of incongruent pairs (effort #2) (effort performance, error trials)
            if mod_perfE.n_incong ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_incong
                    case 1
                        fprintf('Number of incongruent pairs solved');
                    case 2
                        fprintf('Percentage of incongruent pairs solved per trial');
                end
                fprintf('/');
            end
            
            %% number distance (effort #3) (effort performance, error trials)
            if mod_perfE.n_dist ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_dist
                    case 1
                        fprintf('average of the number distance of the pairs solved per trial');
                    case 2
                        fprintf('average of 1/(number distance of the pairs solved) per trial');
                end
                fprintf('/');
            end
            
            %% error numbers (effort performance, error trials)
            if mod_perfE.n_errors ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.n_errors
                    case 1
                        error('n_errors = 1 incompatible with o_perfE = 3');
                    case 2
                        fprintf('Number of errors made per trial across pairs'); % includes errors repeated on the same pair
                    case 3
                        fprintf('Percentage of errors made per trial');
                end
                fprintf('/');
            end
            %% trial number (last regressor) (error trials, effort perf)
            if mod_perfE.trialN ~= 0 && ismember(mod_perfE.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_perfE.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            
            %%
            fprintf('\n');
    end % o_perfE
end % o_perfE#0


%% feedback
o_fbk = stroop.o_fbk;
dur_fbk = stroop.dur_fbk;
mod_fbk = stroop.mod_fbk;

if o_fbk ~= 0
    switch o_fbk
        case 1 % all trials pooled
            %% feedback onset + duration (all trials)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Feedback (all trials): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% Gain/Loss condition (all trials, feedback)
            if mod_fbk.GL_cond ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.GL_cond
                    case 1
                        fprintf('Gain or Loss condition (1/0)');
                    case 2
                        fprintf('Gain or Loss condition (1/-1)');
                end
                fprintf('/');
            end
            %% luminance (first regressor) (all trials, feedback)
            if mod_fbk.lum ~= 0
                %                n_prm_stroop = n_prm_stroop + 1;
                error('luminance for feedback not ready yet.');
            end
            %% trial number (first regressor) (all trials, feedback)
            if mod_fbk.trialN ~= 0 && ismember(mod_fbk.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% ressource (feedback, all trials)
            if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[5,6])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.E_pred
                    case 5
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on actual perf']);
                    case 6
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% feedback (all trials, feedback)
            if mod_fbk.fbk ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.fbk
                    case 1
                        fprintf('feedback');
                end
                fprintf('/');
            end
            %% total gain (all trials, feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
            end
            %% ressource (feedback, all trials)
            if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[1,2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.E_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (feedback, all trials)
            if mod_fbk.perf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.perf
                    case 1
                        fprintf('performance (1-100)');
                    case 2
                        fprintf('performance (0-1)');
                    case 3
                        fprintf('performance (0-1) corrected');
                    case 4
                        fprintf(['performance (0-1) predicted by model ',num2str(mod_fbk.mdl_n)]);
                end
                fprintf('/');
            end
            %% ROI activity as regressor (feedback, all trials)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_stroop = n_prm_stroop + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %% luminance as last regressor
            if mod_fbk.lum ~= 0
                error('luminance for feedback not ready yet.');
            end
            %% trial number (last regressor) (all trials, feedback)
            if mod_fbk.trialN ~= 0 && ismember(mod_fbk.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
        case 2 % split by to gain/to lose trial type
            %% Feedback onset + duration (to gain trials, feedback)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Feedback (to gain trials): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% Gain/Loss condition (to gain trials, feedback)
            if mod_fbk.GL_cond ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.GL_cond
                    case 1
                        fprintf('Gain or Loss condition (1/0)');
                    case 2
                        fprintf('Gain or Loss condition (1/-1)');
                end
                fprintf('/');
            end
            %% luminance (first regressor) (to gain trials, feedback)
            if mod_fbk.lum ~= 0
                %                n_prm_stroop = n_prm_stroop + 1;
                error('luminance for feedback not ready yet.');
            end
            %% trial number (first regressor) (to gain trials, feedback)
            if mod_fbk.trialN ~= 0 && ismember(mod_fbk.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% ressource (feedback, to gain trials)
            if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[5,6])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.E_pred
                    case 5
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on actual perf']);
                    case 6
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% feedback (to gain trials, feedback)
            if mod_fbk.fbk ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.fbk
                    case 1
                        fprintf('feedback');
                end
                fprintf('/');
            end
            %% total gain (to gain trials, feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
            end
            %% ressource (feedback, to gain trials)
            if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[1,2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.E_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (feedback, to gain trials)
            if mod_fbk.perf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.perf
                    case 1
                        fprintf('performance (1-100)');
                    case 2
                        fprintf('performance (0-1)');
                    case 3
                        fprintf('performance (0-1) corrected');
                    case 4
                        fprintf(['performance (0-1) predicted by model ',num2str(mod_fbk.mdl_n)]);
                end
                fprintf('/');
            end
            %% ROI activity as regressor (feedback, to gain trials)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_stroop = n_prm_stroop + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %% luminance (last regressor) (to gain trials, feedback)
            if mod_fbk.lum ~= 0
                %                n_prm_stroop = n_prm_stroop + 1;
                error('luminance for feedback not ready yet.');
            end
            %% trial number (last regressor) (to gain trials, feedback)
            if mod_fbk.trialN ~= 0 && ismember(mod_fbk.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            %% feedback onset + duration (to lose trials, feedback)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Feedback (to lose trials): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% Gain/Loss condition (to lose trials, feedback)
            if mod_fbk.GL_cond ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.GL_cond
                    case 1
                        fprintf('Gain or Loss condition (1/0)');
                    case 2
                        fprintf('Gain or Loss condition (1/-1)');
                end
                fprintf('/');
            end
            %% luminance (first regressor) (to lose trials, feedback)
            if mod_fbk.lum ~= 0
                %                n_prm_stroop = n_prm_stroop + 1;
                error('luminance for feedback not ready yet.');
            end
            %% trial number (first regressor) (to lose trials, feedback)
            if mod_fbk.trialN ~= 0 && ismember(mod_fbk.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% ressource (feedback, to lose trials)
            if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[5,6])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.E_pred
                    case 5
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on actual perf']);
                    case 6
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% feedback (to lose trials, feedback)
            if mod_fbk.fbk ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.fbk
                    case 1
                        fprintf('feedback');
                end
                fprintf('/');
            end
            %% total gain (to lose trials, feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
            end
            %% ressource (feedback, to lose trials)
            if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[1,2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.E_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (feedback, to lose trials)
            if mod_fbk.perf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.perf
                    case 1
                        fprintf('performance (1-100)');
                    case 2
                        fprintf('performance (0-1)');
                    case 3
                        fprintf('performance (0-1) corrected');
                    case 4
                        fprintf(['performance (0-1) predicted by model ',num2str(mod_fbk.mdl_n)]);
                end
                fprintf('/');
            end
            %% ROI activity as regressor (feedback, to lose trials)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_stroop = n_prm_stroop + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %% luminance (last regressor) (to lose trials, feedback)
            if mod_fbk.lum ~= 0
                %                n_prm_stroop = n_prm_stroop + 1;
                error('luminance for feedback not ready yet.');
            end
            %% trial number (last regressor) (to lose trials, feedback)
            if mod_fbk.trialN ~= 0 && ismember(mod_fbk.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
        case 3 % split by no-error/error trial type
            %% Feedback onset + duration (no error trials)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Feedback (no error trials): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% Gain/Loss condition (no error trials, feedback)
            if mod_fbk.GL_cond ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.GL_cond
                    case 1
                        fprintf('Gain or Loss condition (1/0)');
                    case 2
                        fprintf('Gain or Loss condition (1/-1)');
                end
                fprintf('/');
            end
            %% luminance (first regressor) (no error trials, feedback)
            if mod_fbk.lum ~= 0
                %                n_prm_stroop = n_prm_stroop + 1;
                error('luminance for feedback not ready yet.');
            end
            %% trial number (first regressor) (no error trials, feedback)
            if mod_fbk.trialN ~= 0 && ismember(mod_fbk.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% ressource (feedback, no error trials)
            if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[5,6])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.E_pred
                    case 5
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on actual perf']);
                    case 6
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% feedback (no error trials, feedback)
            if mod_fbk.fbk ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.fbk
                    case 1
                        fprintf('feedback');
                end
                fprintf('/');
            end
            %% total gain (no error trials, feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
            end
            %% ressource (feedback, no error trials)
            if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[1,2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.E_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (feedback, no error trials)
            if mod_fbk.perf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.perf
                    case 1
                        fprintf('performance (1-100)');
                    case 2
                        fprintf('performance (0-1)');
                    case 3
                        fprintf('performance (0-1) corrected');
                    case 4
                        fprintf(['performance (0-1) predicted by model ',num2str(mod_fbk.mdl_n)]);
                end
                fprintf('/');
            end
            %% ROI activity as regressor (feedback, no error trials)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_stroop = n_prm_stroop + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %% luminance (last regressor) (no error trials, feedback)
            if mod_fbk.lum ~= 0
                %                n_prm_stroop = n_prm_stroop + 1;
                error('luminance for feedback not ready yet.');
            end
            %% trial number (last regressor) (no error trials, feedback)
            if mod_fbk.trialN ~= 0 && ismember(mod_fbk.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
            
            %% feedback onset + duration (error trials)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Feedback (error trials): ');
            switch dur_fbk
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar');
                otherwise
                    error('not possible');
            end
            %%
            fprintf('\n');
            
            %% Gain/Loss condition (error trials, feedback)
            if mod_fbk.GL_cond ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.GL_cond
                    case 1
                        fprintf('Gain or Loss condition (1/0)');
                    case 2
                        fprintf('Gain or Loss condition (1/-1)');
                end
                fprintf('/');
            end
            %% luminance (first regressor) (error trials, feedback)
            if mod_fbk.lum ~= 0
                %                n_prm_stroop = n_prm_stroop + 1;
                error('luminance for feedback not ready yet.');
            end
            %% trial number (first regressor) (error trials, feedback)
            if mod_fbk.trialN ~= 0 && ismember(mod_fbk.trialN,[2,4])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% ressource (feedback, error trials)
            if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[5,6])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.E_pred
                    case 5
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on actual perf']);
                    case 6
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% feedback (error trials, feedback)
            if mod_fbk.fbk ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.fbk
                    case 1
                        fprintf('feedback');
                end
                fprintf('/');
            end
            %% total gain (error trials, feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
            end
            %% ressource (feedback, error trials)
            if mod_fbk.E_pred ~= 0 && ismember(mod_fbk.E_pred,[1,2])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.E_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (feedback, error trials)
            if mod_fbk.perf ~= 0
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.perf
                    case 1
                        fprintf('performance (1-100)');
                    case 2
                        fprintf('performance (0-1)');
                    case 3
                        fprintf('performance (0-1) corrected');
                    case 4
                        fprintf(['performance (0-1) predicted by model ',num2str(mod_fbk.mdl_n)]);
                end
                fprintf('/');
            end
            %% ROI activity as regressor (feedback, error trials)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_stroop = n_prm_stroop + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %% luminance (last regressor) (error trials, feedback)
            if mod_fbk.lum ~= 0
                %                n_prm_stroop = n_prm_stroop + 1;
                error('luminance for feedback not ready yet.');
            end
            %% trial number (last regressor) (error trials, feedback)
            if mod_fbk.trialN ~= 0 && ismember(mod_fbk.trialN,[1,3])
                n_prm_stroop = n_prm_stroop + 1;
                switch mod_fbk.trialN
                    case 1
                        fprintf('trialN');
                    case 3
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %%
            fprintf('\n');
    end % o_fbk
end % o_fbk#0

%% cross
o_cross = stroop.o_cross;
dur_cross = stroop.dur_cross;
if o_cross ~= 0
    switch o_cross
        case 1
            %% cross onset + duration (all trials)
            n_prm_stroop = n_prm_stroop + 1;
            fprintf('**Cross: ');
            switch dur_cross
                case 0
                    fprintf('stick');
                case 1
                    fprintf('boxcar');
            end
            %%
            fprintf('\n');
    end
end

%% missed trials
o_missed_trials_dispE = stroop.o_missed_trials_dispE;
dur_missed_trials_dispE = stroop.dur_missed_trials_dispE;
if o_missed_trials_dispE ~= 0
    switch o_missed_trials_dispE
        case 1
            %% missed trials onset + duration (all missed trials)
            % number of parameters depends on subject and run
            fprintf('**Missed Trials: ');
            switch dur_missed_trials_dispE
                case 0
                    fprintf('stick function');
                case 1
                    fprintf('boxcar from onset of display effort scale to feedback');
            end
            %%
            fprintf('\n');
    end
end

end % function