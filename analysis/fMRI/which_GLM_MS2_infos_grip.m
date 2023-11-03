function [ n_prm_grip ] = which_GLM_MS2_infos_grip( GLMprm )
%which_GLM_MS2_infos_grip( GLMprm ) outputs the informations about the GLM based
% on the list of parameters designated inside GLMprm for the
% grip task.
%
% INPUTS
% GLMprm: structure with GLM parameters for each task extracted with
% which_GLM_MS2.
%
% OUTPUTS
% n_prm_grip: number of parameters in grip task (includes onsets and
% regressors) for a given run (without taking into account derivative)
%
% See also which_GLM_MS2

%% grip parameters info
grip    = GLMprm.grip;
disp('***Grip parameters');
n_prm_grip = 0;

%% incentive display
o_inc = grip.o_inc;
dur_inc = grip.dur_inc;
mod_inc = grip.mod_inc;
if o_inc ~= 0
    switch o_inc
        case 1 % all trials grouped
            %% incentive onset + duration
            n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first press (when starts squeezing the grip) (incentive, all trials)
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[1])
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.RT_fp
                    case 2
                        fprintf('RT first press');
                end
                fprintf('/');
            end
            %% model variables
            if mod_inc.mdl_n ~= 0
                %% ressource (incentive, all trials)
                if mod_inc.X_pred ~= 0 && ismember(mod_inc.X_pred,[5,6])
                    n_prm_grip = n_prm_grip + 1;
                    switch mod_inc.X_pred
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                fprintf(['GLM',mod_inc.ROI_activity_GLM,...
                    ' ',mod_inc.ROI_activity_period,' period ',...
                    mod_inc.ROI_activity_ROI_nm,' activity ']);
            end
            %% effort (incentive, all trials)
            if mod_inc.Effort ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.Effort
                    case 1
                        fprintf('peakForce');
                    case 2
                        fprintf('integral of Force');
                    case 3
                        fprintf('peakForce normalised by run final Fmax');
                    case 4
                        fprintf('integral of Force normalised by run final Fmax');
                end
                fprintf('/');
            end
            %% model variables
            if mod_inc.mdl_n ~= 0
                %% ressource (incentive, all trials)
                if mod_inc.X_pred ~= 0 && ismember(mod_inc.X_pred,[3,4])
                    n_prm_grip = n_prm_grip + 1;
                    switch mod_inc.X_pred
                        case 3
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                        case 4
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
                %% cost (incentive, all trials)
                if mod_inc.mdl_cost ~= 0
                    n_prm_grip = n_prm_grip + 1;
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
                    n_prm_grip = n_prm_grip + 1;
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
                    n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.trialN
                    case 5
                        fprintf('trialN');
                    case 6
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first press (when starts squeezing the grip) (incentive, all trials)
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[2])
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.RT_fp
                    case 1
                        fprintf('RT first press');
                end
                fprintf('/');
            end
            %% cumulated total gain until now (incentive, all trials)
            if mod_inc.totalGain_prev ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.totalGain_prev
                    case 1
                        fprintf('total gain (until this trial)');
                end
                fprintf('/');
            end
            %% sum(performance previous trials) until now (incentive, all trials)
            if mod_inc.sumPerfPrev ~= 0
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.inc_x_trialN
                    case 1
                        fprintf('inc_x_trialN');     
                end
                fprintf('/');
            end
            %% ressource (incentive, all trials)
            if mod_inc.X_pred ~= 0 && ismember(mod_inc.X_pred,[1,2])
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.X_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (incentive, all trials)
            if mod_inc.perf ~= 0
                n_prm_grip = n_prm_grip + 1;
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
            n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first press (when starts squeezing the grip) (incentive, to gain trials)
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[1])
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.RT_fp
                    case 2
                        fprintf('RT first press');
                end
                fprintf('/');
            end
            %% model variables
            if mod_inc.mdl_n ~= 0
                %% ressource (incentive, to gain trials)
                if mod_inc.X_pred ~= 0 && ismember(mod_inc.X_pred,[5,6])
                    n_prm_grip = n_prm_grip + 1;
                    switch mod_inc.X_pred
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                fprintf(['GLM',mod_inc.ROI_activity_GLM,...
                    ' ',mod_inc.ROI_activity_period,' period ',...
                    mod_inc.ROI_activity_ROI_nm,' activity ']);
            end
            %% Effort (incentive, to gain trials)
            if mod_inc.Effort ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.Effort
                    case 1
                        fprintf('peakForce');
                    case 2
                        fprintf('integral of Force');
                    case 3
                        fprintf('peakForce normalised by run final Fmax');
                    case 4
                        fprintf('integral of Force normalised by run final Fmax');
                end
                fprintf('/');
            end
            %% model variables
            if mod_inc.mdl_n ~= 0
                %% ressource (incentive, to gain trials)
                if mod_inc.X_pred ~= 0 && ismember(mod_inc.X_pred,[3,4])
                    n_prm_grip = n_prm_grip + 1;
                    switch mod_inc.X_pred
                        case 3
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                        case 4
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
                %% cost (incentive, to gain trials)
                if mod_inc.mdl_cost ~= 0
                    n_prm_grip = n_prm_grip + 1;
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
                    n_prm_grip = n_prm_grip + 1;
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
                    n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.trialN
                    case 5
                        fprintf('trialN');
                    case 6
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first press (when starts squeezing the grip) (incentive, to gain trials)
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[2])
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.RT_fp
                    case 1
                        fprintf('RT first press');
                end
                fprintf('/');
            end
            %% cumulated total gain until now (incentive, to gain trials)
            if mod_inc.totalGain_prev ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.totalGain_prev
                    case 1
                        fprintf('total gain (until this trial)');
                end
                fprintf('/');
            end
            %% sum(performance previous trials) until now (incentive, to gain trials)
            if mod_inc.sumPerfPrev ~= 0
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.inc_x_trialN
                    case 1
                        fprintf('inc_x_trialN');     
                end
                fprintf('/');
            end
            %% ressource (incentive, to gain trials)
            if mod_inc.X_pred ~= 0 && ismember(mod_inc.X_pred,[1,2])
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.X_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (incentive, to gain trials)
            if mod_inc.perf ~= 0
                n_prm_grip = n_prm_grip + 1;
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
            n_prm_grip = n_prm_grip + 1;
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
            fprintf('\n');
            
            %% luminance (first regressor) (incentive, to lose trials)
            if mod_inc.lum ~= 0 && ismember(mod_inc.lum,[2,4,6,8])
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first press (when starts squeezing the grip) (incentive, to lose trials)
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[1])
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.RT_fp
                    case 2
                        fprintf('RT first press');
                end
                fprintf('/');
            end
            %% model variables
            if mod_inc.mdl_n ~= 0
                %% ressource (incentive, to lose trials)
                if mod_inc.X_pred ~= 0 && ismember(mod_inc.X_pred,[5,6])
                    n_prm_grip = n_prm_grip + 1;
                    switch mod_inc.X_pred
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                fprintf(['GLM',mod_inc.ROI_activity_GLM,...
                    ' ',mod_inc.ROI_activity_period,' period ',...
                    mod_inc.ROI_activity_ROI_nm,' activity ']);
            end
            %% Effort (incentive, to lose trials)
            if mod_inc.Effort ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.Effort
                    case 1
                        fprintf('peakForce');
                    case 2
                        fprintf('integral of Force');
                    case 3
                        fprintf('peakForce normalised by run final Fmax');
                    case 4
                        fprintf('integral of Force normalised by run final Fmax');
                end
                fprintf('/');
            end
            %% model variables
            if mod_inc.mdl_n ~= 0
                %% ressource (incentive, to lose trials)
                if mod_inc.X_pred ~= 0 && ismember(mod_inc.X_pred,[3,4])
                    n_prm_grip = n_prm_grip + 1;
                    switch mod_inc.X_pred
                        case 3
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                        case 4
                            fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
                %% cost (incentive, to lose trials)
                if mod_inc.mdl_cost ~= 0
                    n_prm_grip = n_prm_grip + 1;
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
                    n_prm_grip = n_prm_grip + 1;
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
                    n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.trialN
                    case 5
                        fprintf('trialN');
                    case 6
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first press (when starts squeezing the grip) (incentive, to lose trials)
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[2])
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.RT_fp
                    case 1
                        fprintf('RT first press');
                end
                fprintf('/');
            end
            %% cumulated total gain until now (incentive, to lose trials)
            if mod_inc.totalGain_prev ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.totalGain_prev
                    case 1
                        fprintf('total gain (until this trial)');
                end
                fprintf('/');
            end
            %% sum(performance previous trials) until now (incentive, to lose trials)
            if mod_inc.sumPerfPrev ~= 0
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.inc_x_trialN
                    case 1
                        fprintf('inc_x_trialN');     
                end
                fprintf('/');
            end
            %% ressource (incentive, to lose trials)
            if mod_inc.X_pred ~= 0 && ismember(mod_inc.X_pred,[1,2])
                n_prm_grip = n_prm_grip + 1;
                switch mod_inc.X_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_inc.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (incentive, to lose trials)
            if mod_inc.perf ~= 0
                n_prm_grip = n_prm_grip + 1;
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
o_dispE = grip.o_dispE;
dur_dispE = grip.dur_dispE;
mod_dispE = grip.mod_dispE;
if o_dispE ~= 0
    switch o_dispE
        case 1 % trials pooled
            %% effort scale display onset + duration
            n_prm_grip = n_prm_grip + 1;
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
            
            %% Gain/Loss condition (all trials, effort scale)
            if mod_dispE.GL_cond ~= 0
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first press (when starts squeezing the grip) (all trials, effort scale)
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[1])
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.RT_fp
                    case 2
                        fprintf('RT first press');
                end
                fprintf('/');
            end
            %% model variables
            if mod_dispE.mdl_n ~= 0
                %% ressource (effort scale, all trials)
                if mod_dispE.X_pred ~= 0 && ismember(mod_dispE.X_pred,[5,6])
                    n_prm_grip = n_prm_grip + 1;
                    switch mod_dispE.X_pred
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.R_type
                    case 1
                        fprintf('reward type: piece of money or bill (0/1)');
                    otherwise
                        error('not ready');
                end
                fprintf('/');
            end
            %% incentive (all trials, effort scale)
            if mod_dispE.inc ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.inc
                    case 1
                        fprintf('incentive rank');
                    case 2
                        fprintf('incentive value');
                end
                fprintf('/');
            end
            %% incentive bis (all trials, effort scale)
            if mod_dispE.inc_bis ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.inc_bis
                    case 1
                        fprintf('incentive rank');
                    case 2
                        fprintf('incentive value');
                end
                fprintf('/');
            end
            %% confidence proxy (effort scale, all trials)
            if mod_dispE.conf ~= 0
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                fprintf(['GLM',mod_dispE.ROI_activity_GLM,...
                    ' ',mod_dispE.ROI_activity_period,' period ',...
                    mod_dispE.ROI_activity_ROI_nm,' activity ']);
            end
            %% effort intensity (all trials, effort scale)
            if mod_dispE.Effort ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.Effort
                    case 1
                        fprintf('peakForce');
                    case 2
                        fprintf('integral of Force');
                    case 3
                        fprintf('peakForce normalised by run final Fmax');
                    case 4
                        fprintf('integral of Force normalised by run final Fmax');
                end
                fprintf('/');
            end
            %% model variables
            if mod_dispE.mdl_n ~= 0
                %% ressource (effort scale, all trials)
                if mod_dispE.X_pred ~= 0 && ismember(mod_dispE.X_pred,[3,4])
                    n_prm_grip = n_prm_grip + 1;
                    switch mod_dispE.X_pred
                        case 3
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                        case 4
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
                %% cost (all trials, effort scale)
                if mod_dispE.mdl_cost ~= 0
                    n_prm_grip = n_prm_grip + 1;
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
                    n_prm_grip = n_prm_grip + 1;
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
                    n_prm_grip = n_prm_grip + 1;
                    switch mod_dispE.mdl_EV
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' EV based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' EV based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% RT first press (when starts squeezing the grip) (all trials, effort scale)
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[2])
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.RT_fp
                    case 1
                        fprintf('RT first press');
                end
                fprintf('/');
            end
            %% luminance (last regressor) (effort scale, all trials)
            if mod_dispE.lum ~= 0 && ismember(mod_dispE.lum,[1,3,5,7])
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.inc_x_trialN
                    case 1
                        fprintf('inc_x_trialN');     
                end
                fprintf('/');
            end
            %% ressource (effort scale, all trials)
            if mod_dispE.X_pred ~= 0 && ismember(mod_dispE.X_pred,[1,2])
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.X_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (effort scale, all trials)
            if mod_dispE.perf ~= 0
                n_prm_grip = n_prm_grip + 1;
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
            %% effort scale onset + duration (to gain trials)
            n_prm_grip = n_prm_grip + 1;
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
            
            %% luminance (first regressor) (to gain trials, effort scale)
            if mod_dispE.lum ~= 0 && ismember(mod_dispE.lum,[2,4,6,8])
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first press (when starts squeezing the grip) (to gain trials, effort scale)
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[1])
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.RT_fp
                    case 2
                        fprintf('RT first press');
                end
                fprintf('/');
            end
            %% model variables
            if mod_dispE.mdl_n ~= 0
                %% ressource (effort scale, to gain trials)
                if mod_dispE.X_pred ~= 0 && ismember(mod_dispE.X_pred,[5,6])
                    n_prm_grip = n_prm_grip + 1;
                    switch mod_dispE.X_pred
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.R_type
                    case 1
                        fprintf('reward type: piece of money or bill (0/1)');
                    otherwise
                        error('not ready');
                end
                fprintf('/');
            end
            %% incentive (to gain trials, effort scale)
            if mod_dispE.inc ~= 0
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                fprintf(['GLM',mod_dispE.ROI_activity_GLM,...
                    ' ',mod_dispE.ROI_activity_period,' period ',...
                    mod_dispE.ROI_activity_ROI_nm,' activity ']);
            end
            %% effort intensity (to gain trials, effort scale)
            if mod_dispE.Effort ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.Effort
                    case 1
                        fprintf('peakForce');
                    case 2
                        fprintf('integral of Force');
                    case 3
                        fprintf('peakForce normalised by run final Fmax');
                    case 4
                        fprintf('integral of Force normalised by run final Fmax');
                end
                fprintf('/');
            end
            %% model variables
            if mod_dispE.mdl_n ~= 0
                %% ressource (effort scale, to gain trials)
                if mod_dispE.X_pred ~= 0 && ismember(mod_dispE.X_pred,[3,4])
                    n_prm_grip = n_prm_grip + 1;
                    switch mod_dispE.X_pred
                        case 3
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                        case 4
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
                %% cost (to gain trials, effort scale)
                if mod_dispE.mdl_cost ~= 0
                    n_prm_grip = n_prm_grip + 1;
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
                    n_prm_grip = n_prm_grip + 1;
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
                    n_prm_grip = n_prm_grip + 1;
                    switch mod_dispE.mdl_EV
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' EV based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' EV based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% RT first press (when starts squeezing the grip) (to gain trials, effort scale)
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[2])
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.RT_fp
                    case 1
                        fprintf('RT first press');
                end
                fprintf('/');
            end
            %% luminance (last regressor) (to gain trials, effort scale)
            if mod_dispE.lum ~= 0 && ismember(mod_dispE.lum,[1,3,5,7])
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.inc_x_trialN
                    case 1
                        fprintf('inc_x_trialN');     
                end
                fprintf('/');
            end
            %% ressource (effort scale, to gain trials)
            if mod_dispE.X_pred ~= 0 && ismember(mod_dispE.X_pred,[1,2])
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.X_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (effort scale, to gain trials)
            if mod_dispE.perf ~= 0
                n_prm_grip = n_prm_grip + 1;
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
            
            %% effort scale onset + duration (to lose trials)
            n_prm_grip = n_prm_grip + 1;
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
            
            %% luminance (first regressor) (to lose trials, effort scale)
            if mod_dispE.lum ~= 0 && ismember(mod_dispE.lum,[2,4,6,8])
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% RT first press (when starts squeezing the grip) (to lose trials, effort scale)
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[1])
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.RT_fp
                    case 2
                        fprintf('RT first press');
                end
                fprintf('/');
            end
            %% model variables
            if mod_dispE.mdl_n ~= 0
                %% ressource (effort scale, to lose trials)
                if mod_dispE.X_pred ~= 0 && ismember(mod_dispE.X_pred,[5,6])
                    n_prm_grip = n_prm_grip + 1;
                    switch mod_dispE.X_pred
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.R_type
                    case 1
                        fprintf('reward type: piece of money or bill (0/1)');
                    otherwise
                        error('not ready');
                end
                fprintf('/');
            end
            %% incentive (to lose trials, effort scale)
            if mod_dispE.inc ~= 0
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                fprintf(['GLM',mod_dispE.ROI_activity_GLM,...
                    ' ',mod_dispE.ROI_activity_period,' period ',...
                    mod_dispE.ROI_activity_ROI_nm,' activity ']);
            end
            %% effort intensity (to lose trials, effort scale)
            if mod_dispE.Effort ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.Effort
                    case 1
                        fprintf('peakForce');
                    case 2
                        fprintf('integral of Force');
                    case 3
                        fprintf('peakForce normalised by run final Fmax');
                    case 4
                        fprintf('integral of Force normalised by run final Fmax');
                end
                fprintf('/');
            end
            %% model variables
            if mod_dispE.mdl_n ~= 0
                %% ressource (effort scale, to lose trials)
                if mod_dispE.X_pred ~= 0 && ismember(mod_dispE.X_pred,[3,4])
                    n_prm_grip = n_prm_grip + 1;
                    switch mod_dispE.X_pred
                        case 3
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                        case 4
                            fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                    end
                    fprintf('/');
                end
                %% cost (to lose trials, effort scale)
                if mod_dispE.mdl_cost ~= 0
                    n_prm_grip = n_prm_grip + 1;
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
                    n_prm_grip = n_prm_grip + 1;
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
                    n_prm_grip = n_prm_grip + 1;
                    switch mod_dispE.mdl_EV
                        case 1
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' EV based on actual perf']);
                        case 2
                            fprintf(['model ',num2str(mod_dispE.mdl_n),' EV based on predicted perf']);
                    end
                    fprintf('/');
                end
            end
            %% RT first press (when starts squeezing the grip) (to lose trials, effort scale)
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[2])
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.RT_fp
                    case 1
                        fprintf('RT first press');
                end
                fprintf('/');
            end
            %% luminance (last regressor) (to lose trials, effort scale)
            if mod_dispE.lum ~= 0 && ismember(mod_dispE.lum,[1,3,5,7])
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.inc_x_trialN
                    case 1
                        fprintf('inc_x_trialN');     
                end
                fprintf('/');
            end
            %% ressource (effort scale, to lose trials)
            if mod_dispE.X_pred ~= 0 && ismember(mod_dispE.X_pred,[1,2])
                n_prm_grip = n_prm_grip + 1;
                switch mod_dispE.X_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_dispE.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (effort scale, to lose trials)
            if mod_dispE.perf ~= 0
                n_prm_grip = n_prm_grip + 1;
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
    end
end

%% effort performance onset
o_perfE = grip.o_perfE;
dur_perfE = grip.dur_perfE;
mod_perfE = grip.mod_perfE;
if o_perfE ~= 0
    switch o_perfE
        case 1 % trials pooled
            %% effort performance onset + duration (all trials)
            n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_perfE.GL_cond
                    case 1
                        fprintf('Gain or Loss condition (1/0)');
                    case 2
                        fprintf('Gain or Loss condition (1/-1)');
                end
                fprintf('/');
            end
            %% piece/bill condition (effort scale, to lose trials)
            if mod_perfE.R_type ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_perfE.R_type
                    case 1
                        fprintf('reward type: piece of money or bill (0/1)');
                    otherwise
                        error('not ready');
                end
                fprintf('/');
            end
            %% trial number (first regressor) (all trials, effort perf)
            if mod_perfE.trialN ~= 0 && ismember(mod_perfE.trialN,[2,4])
                n_prm_grip = n_prm_grip + 1;
                switch mod_perfE.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% effort intensity (all trials)
            if mod_perfE.Effort ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_perfE.Effort
                    case 1
                        fprintf('peakForce');
                    case 2
                        fprintf('integral of Force');
                    case 3
                        fprintf('peakForce normalised by run final Fmax');
                    case 4
                        fprintf('integral of Force normalised by run final Fmax');
                end
                fprintf('/');
            end
            %% trial number (last regressor) (all trials, effort perf)
            if mod_perfE.trialN ~= 0 && ismember(mod_perfE.trialN,[1,3])
                n_prm_grip = n_prm_grip + 1;
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
            
        case 2 % split by to gain/to lose trial type
            %% effort performance onset + duration (to gain trials, effort perf)
            n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_perfE.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            
            %% effort intensity (to gain trials, effort perf)
            if mod_perfE.Effort ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_perfE.Effort
                    case 1
                        fprintf('peakForce');
                    case 2
                        fprintf('integral of Force');
                    case 3
                        fprintf('peakForce normalised by run final Fmax');
                    case 4
                        fprintf('integral of Force normalised by run final Fmax');
                end
                fprintf('/');
            end
            
            %% trial number (last regressor) (to gain trials, effort perf)
            if mod_perfE.trialN ~= 0 && ismember(mod_perfE.trialN,[1,3])
                n_prm_grip = n_prm_grip + 1;
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
            
            %% Effort performance (to lose trials, effort perf)
            n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                switch mod_perfE.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% effort intensity (to lose trials, effort perf)
            if mod_perfE.Effort ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_perfE.Effort
                    case 1
                        fprintf('peakForce');
                    case 2
                        fprintf('integral of Force');
                    case 3
                        fprintf('peakForce normalised by run final Fmax');
                    case 4
                        fprintf('integral of Force normalised by run final Fmax');
                end
                fprintf('/');
            end
            %% trial number (last regressor) (to lose trials, effort perf)
            if mod_perfE.trialN ~= 0 && ismember(mod_perfE.trialN,[1,3])
                n_prm_grip = n_prm_grip + 1;
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
    end
end

%% feedback
o_fbk = grip.o_fbk;
dur_fbk = grip.dur_fbk;
mod_fbk = grip.mod_fbk;

if o_fbk ~= 0
    switch o_fbk
        case 1 % all trials pooled
            %% feedback onset + duration (all trials)
            n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                %                n_prm_grip = n_prm_grip + 1;
                error('luminance for feedback not ready yet.');
            end
            %% trial number (first regressor) (all trials, feedback)
            if mod_fbk.trialN ~= 0 && ismember(mod_fbk.trialN,[2,4])
                n_prm_grip = n_prm_grip + 1;
                switch mod_fbk.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% ressource (feedback, all trials)
            if mod_fbk.X_pred ~= 0 && ismember(mod_fbk.X_pred,[5,6])
                n_prm_grip = n_prm_grip + 1;
                switch mod_fbk.X_pred
                    case 5
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on actual perf']);
                    case 6
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% feedback (all trials, feedback)
            if mod_fbk.fbk ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_fbk.fbk
                    case 1
                        fprintf('feedback');
                end
                fprintf('/');
            end
            %% total gain (all trials, feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
            end
            %% ressource (feedback, all trials)
            if mod_fbk.X_pred ~= 0 && ismember(mod_fbk.X_pred,[1,2])
                n_prm_grip = n_prm_grip + 1;
                switch mod_fbk.X_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (feedback, all trials)
            if mod_fbk.perf ~= 0
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
            %% Feedback onset + duration (to gain trials)
            n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                %                n_prm_grip = n_prm_grip + 1;
                error('luminance for feedback not ready yet.');
            end
            %% trial number (first regressor) (to gain trials, feedback)
            if mod_fbk.trialN ~= 0 && ismember(mod_fbk.trialN,[2,4])
                n_prm_grip = n_prm_grip + 1;
                switch mod_fbk.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% feedback (to gain trials, feedback)
            if mod_fbk.fbk ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_fbk.fbk
                    case 1
                        fprintf('feedback');
                end
                fprintf('/');
            end
            %% total gain (to gain trials, feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
            end
            %% ressource (feedback, to gain trials)
            if mod_fbk.X_pred ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_fbk.X_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (feedback, to gain trials)
            if mod_fbk.perf ~= 0
                n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %% luminance (last regressor) (to gain trials, feedback)
            if mod_fbk.lum ~= 0
                %                n_prm_grip = n_prm_grip + 1;
                error('luminance for feedback not ready yet.');
            end
            %% trial number (last regressor) (to gain trials, feedback)
            if mod_fbk.trialN ~= 0 && ismember(mod_fbk.trialN,[1,3])
                n_prm_grip = n_prm_grip + 1;
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
            
            %% feedback onset + duration (to lose trials)
            n_prm_grip = n_prm_grip + 1;
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
                n_prm_grip = n_prm_grip + 1;
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
                %                n_prm_grip = n_prm_grip + 1;
                error('luminance for feedback not ready yet.');
            end
            %% trial number (first regressor) (to lose trials, feedback)
            if mod_fbk.trialN ~= 0 && ismember(mod_fbk.trialN,[2,4])
                n_prm_grip = n_prm_grip + 1;
                switch mod_fbk.trialN
                    case 2
                        fprintf('trialN');
                    case 4
                        fprintf('trialN zPerRun');
                end
                fprintf('/');
            end
            %% feedback (to lose trials, feedback)
            if mod_fbk.fbk ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_fbk.fbk
                    case 1
                        fprintf('feedback');
                end
                fprintf('/');
            end
            %% total gain (to lose trials, feedback)
            if mod_fbk.totalGain ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_fbk.totalGain
                    case 1
                        fprintf('total gain');
                end
            end
            %% ressource (feedback, to lose trials)
            if mod_fbk.X_pred ~= 0
                n_prm_grip = n_prm_grip + 1;
                switch mod_fbk.X_pred
                    case 1
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on actual perf']);
                    case 2
                        fprintf(['ressource predicted by model ',num2str(mod_fbk.mdl_n),' based on predicted perf']);
                end
                fprintf('/');
            end
            %% performance (effort scale, to lose trials)
            if mod_fbk.perf ~= 0
                n_prm_grip = n_prm_grip + 1;
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
            %% ROI activity as regressor (effort scale, to lose trials)
            if mod_fbk.ROI_activity_yn == 1
                n_prm_grip = n_prm_grip + 1;
                fprintf(['GLM',mod_fbk.ROI_activity_GLM,...
                    ' ',mod_fbk.ROI_activity_period,' period ',...
                    mod_fbk.ROI_activity_ROI_nm,' activity ']);
            end
            %% luminance (last regressor) (to lose trials, feedback)
            if mod_fbk.lum ~= 0
                %                n_prm_grip = n_prm_grip + 1;
                error('luminance for feedback not ready yet.');
            end
            %% trial number (last regressor) (to lose trials, feedback)
            if mod_fbk.trialN ~= 0 && ismember(mod_fbk.trialN,[1,3])
                n_prm_grip = n_prm_grip + 1;
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
    end
end

%% cross
o_cross = grip.o_cross;
dur_cross = grip.dur_cross;
if o_cross ~= 0
    switch o_cross
        case 1
            %% cross onset + duration (all trials)
            n_prm_grip = n_prm_grip + 1;
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
o_missed_trials_dispE = grip.o_missed_trials_dispE;
dur_missed_trials_dispE = grip.dur_missed_trials_dispE;
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