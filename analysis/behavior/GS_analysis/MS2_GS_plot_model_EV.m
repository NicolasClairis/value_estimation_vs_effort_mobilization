function [ EV, benefit, cost,...
    F_pred, perf_pred, E_pred,...
    EV_pred, benefit_pred, cost_pred, E,...
    benefit_var, cost_var, R2,...
    expected_payoff, expected_payoff_pred] = MS2_GS_plot_model_EV(task_nm, model_nm, sub_nm, run_nm)
%[ EV, benefit, cost,...
%     F_pred, perf_pred, E_pred,...
%     EV_pred, benefit_pred, cost_pred, E,...
%     benefit_var, cost_var, R2,...
%     expected_payoff, expected_payoff_pred] = MS2_GS_plot_model_EV(task_nm, model_nm, sub_nm, run_nm)
% MS2_GS_plot_model_EV plots the expected value for each task, based on the
% parameters estimated by the model
%
% INPUTS
% task_nm: 'grip'/'stroop'
%
% model_nm: 'model_XX'
%
% sub_nm: 'sXX_DDMMYY'
%
% run_nm: '2'/'3'/'5'/'6'
%
% OUTPUTS
% EV: expected value as fitted by the model using the actual performance
%
% benefit: benefit term fitted by the model using the actual performance
%
% cost: cost term fitted by the model using the actual performance
%
% F_pred: Force level predicted by the model
%
% E_pred: level of ressource E* predicted by the model (when model predicts
% it)
%
% perf_pred: performance predicted by the model (based on the force
% predicted by the model F_pred and the way to normalize by Fmax)
%
% EV_pred: expected value as fitted by the model using the predicted performance
%
% benefit_pred: benefit term fitted by the model using the predicted performance
%
% cost_pred: cost term fitted by the model using the predicted performance
%
% E: level of ressource predicted E by the model (based on actual
% performance and not on predicted performance)
%
% benefit_var: benefit term before scaling with performance
%
% cost_var: cost term before scaling with performance
%
% R2: R² for the model
%
% expected_payoff: part multiplying the incentives alone
%
% expected_payoff_pred: expected payoff term fitted by the model using the predicted performance
%
% See also MS2_GS_perf_optimal_level.m

%% compute model and load parameters
[ prm_in_mdl ] = MS2_GS_Festimation_model_space();
mdl_prm = prm_in_mdl.(model_nm);
% compute model
[prm_quality, prm] = MS2_GS_perf_optimal_level( mdl_prm, sub_nm );
prm_kmax_fixed_or_free = mdl_prm.C_kmax_fixed_or_free;
prm_Fmax_fixed_or_free = mdl_prm.Fmax_fixed_or_free;
R2 = prm_quality.(task_nm).prm_qlty.R2;

prm_norm = mdl_prm.norm;
prm_incentive_type = mdl_prm.incentive_type;
prm_incentive_var = mdl_prm.incentive_var;
prm_incentive_type_bis = mdl_prm.incentive_type_bis;
prm_incentive_var_bis = mdl_prm.incentive_var_bis;
prm_fatigue_var = mdl_prm.fatigue_var;
prm_B_time_on_benef = mdl_prm.B_time_on_benef;
prm_C_fatigue = mdl_prm.C_fatigue;
prm_C_rest = mdl_prm.C_rest;
prm_B_perf_C_perf_or_force = mdl_prm.B_perf_C_perf_or_force;
switch task_nm
    case 'grip'
        Fmax_scale = 100/75;
end

%% subject identification
if strcmp(sub_nm(3),'_')
    subid   = sub_nm(2);
elseif ~strcmp(sub_nm(3),'_') && strcmp(sub_nm(4),'_')
    subid = sub_nm(2:3);
end

%% working directories
root = 'enter path here';
sub_folder             = fullfile(root, sub_nm);
sub_onsets_folder      = [fullfile(sub_folder,'fMRI_analysis'),filesep];

%% load model parameters
kCost       = prm.(task_nm).kCost;
kTcost      = prm.(task_nm).kTcost;
switch prm_incentive_var
    case {'inc','absInc'}
        kI          = prm.(task_nm).kI;
    case {'absInc_plus_nomInc','absInc_plus_cond',...
            'inc_perCond'}
        kI1          = prm.(task_nm).kI1;
        kI2          = prm.(task_nm).kI2;
end
kaI1 = 1;
kaI2 = 1;
switch prm_incentive_var_bis
    case {'inc','absInc' }
        kI_bis          = prm.(task_nm).kI;
    case {'absInc_plus_nomInc','absInc_plus_cond',...
            'inc_perCond'}
        kI_bis1          = prm.(task_nm).kI_bis1;
        kI_bis2          = prm.(task_nm).kI_bis2;
end
kRest       = prm.(task_nm).kRest;
kTreward    = prm.(task_nm).kTreward;
if strcmp(prm_kmax_fixed_or_free,'free')
    switch run_nm
        case {'2','3'}
            kmax_run = prm.(task_nm).kmax_estim(1);
        case {'5','6'}
            kmax_run = prm.(task_nm).kmax_estim(2);
    end
else
    kmax_run = 1;
end

switch prm_Fmax_fixed_or_free
    case {'fixed','free','free_bis'}
        switch run_nm
            case {'2','3'}
                Fmax_estim_run = prm.(task_nm).Fmax_estim(1);
            case {'5','6'}
                Fmax_estim_run = prm.(task_nm).Fmax_estim(2);
        end
end

if ismember(prm_B_perf_C_perf_or_force,...
        {'perf_f_X','perf_f_X_bis','perf_f_X_ter',...
        'perf_f_X_4','perf_f_X_5','perf_f_X_6','perf_f_X_7'})
    kX              = prm.(task_nm).kX;
    X_predicted     = prm.(task_nm).X_predicted;
end

%% load vars of interest for the current subject
max_nTrial = 60;
loadStruct = getfield(load([sub_onsets_folder,...
    'onsets_sub',subid,'_',task_nm,'_run',run_nm,'.mat']),task_nm);
switch prm_incentive_type
    case 'rank'
        absInc          = loadStruct.mod.all.absIncentiveRank;
        inc             = loadStruct.mod.all.incentiveRank;
    case 'value'
        absInc          = loadStruct.mod.all.absIncentive;
        inc             = loadStruct.mod.all.incentive;
end
switch prm_incentive_type_bis
    case 'rank'
        absInc_bis          = loadStruct.mod.all.absIncentiveRank;
        inc_bis             = loadStruct.mod.all.incentiveRank;
    case 'value'
        absInc_bis          = loadStruct.mod.all.absIncentive;
        inc_bis             = loadStruct.mod.all.incentive;
end
inc_type        = loadStruct.mod.all.inc_type;
trialValence    = loadStruct.mod.all.trialValence;
GLcond = trialValence;
GLcond(trialValence == -1) = 0;
[absIncGain, absIncLoss] = deal( absInc );
absIncGain(GLcond == 0) = 0;
absIncLoss(GLcond == 1) = 0;
if ~isempty(prm_incentive_type_bis) && ~isempty(prm_incentive_var_bis)
    [absIncGain_bis, absIncLoss_bis] = deal( absInc_bis );
    absIncGain_bis(GLcond == 0) = 0;
    absIncLoss_bis(GLcond == 1) = 0;
end
trialN          = loadStruct.mod.all.trialN;
totalGain_prev  = loadStruct.mod.all.totalGain_prev;
sumPerf_prev = loadStruct.mod.all.sumPerf_prev;
perf            = loadStruct.mod.all.perf./100; % set between 0 and 1
time_rest       = loadStruct.duration.all.ITI +...
    loadStruct.duration.all.incentive; % between [1.5; 4.95]seconds
switch prm_incentive_var
    case 'inc'
        inc_var = inc.*((GLcond == 1) - (GLcond == 0));
    case {'absInc'}
        inc_var = absInc;
    case 'absInc_plus_nomInc'
        inc_var1 = absInc;
        inc_var2 = inc.*((GLcond == 1) - (GLcond == 0));
    case {'absInc_plus_cond'}
        inc_var1 = absInc;
        inc_var2 = GLcond;
    case {'inc_perCond'}
        inc_var1 = absIncGain;
        inc_var2 = absIncLoss;
    otherwise
        error('inc_var name to fix');
end
switch prm_incentive_var_bis
    case ''
        inc_var_bis = [];
    case 'inc'
        inc_var_bis = inc_bis.*((GLcond == 1) - (GLcond == 0));
    case {'absInc'}
        inc_var_bis = absInc_bis;
    case 'absInc_plus_nomInc'
        inc_var_bis1 = absInc_bis;
        inc_var_bis2 = inc_bis.*((GLcond == 1) - (GLcond == 0));
    case {'absInc_plus_cond'}
        inc_var_bis1 = absInc_bis;
        inc_var_bis2 = GLcond;
    case {'inc_perCond'}
        inc_var_bis1 = absIncGain_bis;
        inc_var_bis2 = absIncLoss_bis;
    otherwise
        error('inc_var_bis name to fix');
end
if strcmp(task_nm, 'grip')
    Fmax = loadStruct.mod.all.Fmax;
end

%% Force and Performance predicted extraction

% extract force and performance level predicted by the model (across all
% trials)
perf_predicted  = prm.(task_nm).perf_predicted;

% filter only trials for the current run
run_nb = str2double(run_nm);
if ismember(run_nb, [2,3])
    trial_idx = 1:max_nTrial;
elseif ismember(run_nb,[5,6])
    trial_idx = (1:max_nTrial) + max_nTrial;
end
F_pred = [];
perf_pred   = perf_predicted(1, trial_idx);

if ismember(prm_B_perf_C_perf_or_force,...
        {'perf_f_X','perf_f_X_bis','perf_f_X_ter',...
        'perf_f_X_4','perf_f_X_5','perf_f_X_6','perf_f_X_7'})
    E_pred = X_predicted(1, trial_idx);
end

%% filter NaN trials
perf_pred   = perf_pred(trialN);

if ismember(prm_B_perf_C_perf_or_force,...
        {'perf_f_X','perf_f_X_bis','perf_f_X_ter',...
        'perf_f_X_4','perf_f_X_5','perf_f_X_6','perf_f_X_7'})
    E_pred = E_pred(trialN);
else
    E_pred = [];
end

%% normalize vars
if prm_norm == 1
    switch prm_incentive_type
        case 'value'
            inc_norm = 20; % max 20 €
        case 'rank'
            inc_norm = 6; % max rank = 6
    end
    switch prm_incentive_var
        case {'inc','absInc'}
            inc_var = inc_var./inc_norm;
        case {'absInc_plus_nomInc',...
                'inc_perCond'}
            inc_var1 = inc_var1./inc_norm;
            inc_var2 = inc_var2./inc_norm;
        case {'absInc_plus_cond'}
            % when 2nd variable = binary variable => should not be normalized
            inc_var1 = inc_var1./inc_norm;
    end


    switch prm_incentive_type_bis
        case 'value'
            inc_norm_bis = 20; % max 20 €
        case 'rank'
            inc_norm_bis = 6; % max rank = 6
    end
    switch prm_incentive_var_bis
        case {'inc','absInc'}
            inc_var_bis = inc_var_bis./inc_norm_bis;
        case {'absInc_plus_nomInc',...
                'inc_perCond'}
            inc_var_bis1 = inc_var_bis1./inc_norm_bis;
            inc_var_bis2 = inc_var_bis2./inc_norm_bis;
        case {'absInc_plus_cond'}
            % when 2nd variable = binary variable => should not be normalized
            inc_var_bis1 = inc_var_bis1./inc_norm_bis;
    end
    trialN = trialN./max_nTrial;
    time_rest = time_rest./4.95; % if you want to normalize it,
    %         maximum possible = 4.95 seconds of rest (0.5s cross + 3.95s
    %         incentive)
    maxTotalGain = 267.10;
    totalGain_prev = totalGain_prev./maxTotalGain;
    sumPerf_prev = sumPerf_prev./max_nTrial;
end

%% select which fatigue var to use
switch prm_fatigue_var
    case 'trialN'
        fatigue_var = trialN;
    case 'totalGainPrev'
        fatigue_var = totalGain_prev;
    case 'sumPerfPrev'
        fatigue_var = sumPerf_prev;
end
%% compute benefit and cost terms
%% cost term
switch prm_C_fatigue
    case 0
        C_fatigue_effect = 1;
    case 1
        C_fatigue_effect = 1 + kTcost.*fatigue_var;
    case 2
        C_fatigue_effect = 1 + exp(kTcost.*fatigue_var);
    case 3
        C_fatigue_effect = 1 + log(1 + kTcost.*fatigue_var);
    case 4
        C_fatigue_effect = 1 - exp(kTcost.*fatigue_var);
end
switch prm_C_rest
    case 0 % no rest effect
        cost_var    = kCost.*C_fatigue_effect;
    case 1 % subtract by rest
        cost_var    = kCost.*(C_fatigue_effect - kRest.*time_rest);
    case 2 % divide by rest
        cost_var    = kCost.*C_fatigue_effect./(1 + kRest.*time_rest);
    case 3 % multiply by rest
        cost_var    = kCost.*C_fatigue_effect.*(1 - kRest.*time_rest); %% multiplication = weird, should check range of values
end

%% benefit term
switch prm_incentive_var
    case {'inc','absInc'}
        switch prm_incentive_var_bis
            case ''
                switch prm_B_time_on_benef
                    case 0 % no rest effect
                        benefit_var = (1 +kI.*inc_var);
                    case 1 % subtract by fatigue
                        benefit_var = (1 +kI.*inc_var - kTreward.*fatigue_var);
                    case 2 % divide by fatigue
                        benefit_var = (1 +kI.*inc_var)./(1 + kTreward.*fatigue_var);
                    case 3 % multiply by fatigue
                        benefit_var = (1 +kI.*inc_var).*(1 - kTreward.*fatigue_var);
                end
                expected_payoff_var = kI.*inc_var;
            case {'inc','absInc'}
                switch prm_B_time_on_benef
                    case 0 % no rest effect
                        benefit_var = (1 +kI.*inc_var + kI_bis.*inc_var_bis);
                    case 1 % subtract by fatigue
                        benefit_var = (1 +kI.*inc_var + kI_bis.*inc_var_bis - kTreward.*fatigue_var);
                    case 2 % divide by fatigue
                        benefit_var = (1 +kI.*inc_var + kI_bis.*inc_var_bis)./(1 + kTreward.*fatigue_var);
                    case 3 % multiply by fatigue
                        benefit_var = (1 +kI.*inc_var + kI_bis.*inc_var_bis).*(1 - kTreward.*fatigue_var);
                end
                expected_payoff_var = kI.*inc_var + kI_bis.*inc_var_bis;
            case {'absInc_plus_nomInc','absInc_plus_cond',...
                    'inc_perCond'}
                switch prm_B_time_on_benef
                    case 0 % no rest effect
                        benefit_var = (1 +kI.*inc_var + kI_bis1.*inc_var_bis1 + kI_bis2.*inc_var_bis2);
                    case 1 % subtract by fatigue
                        benefit_var = (1 +kI.*inc_var + kI_bis1.*inc_var_bis1 + kI_bis2.*inc_var_bis2 - kTreward.*fatigue_var);
                    case 2 % divide by fatigue
                        benefit_var = (1 +kI.*inc_var + kI_bis1.*inc_var_bis1 + kI_bis2.*inc_var_bis2)./(1 + kTreward.*fatigue_var);
                    case 3 % multiply by fatigue
                        benefit_var = (1 +kI.*inc_var + kI_bis1.*inc_var_bis1 + kI_bis2.*inc_var_bis2).*(1 - kTreward.*fatigue_var);
                end
                expected_payoff_var = kI.*inc_var + kI_bis1.*inc_var_bis1 + kI_bis2.*inc_var_bis2;
        end % incentive bis

    case {'absInc_plus_nomInc','absInc_plus_cond',...
            'inc_perCond'}
        switch prm_incentive_var_bis
            case ''
                switch prm_B_time_on_benef
                    case 0 % no rest effect
                        benefit_var = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2));
                    case 1 % subtract by fatigue
                        benefit_var = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) - kTreward.*fatigue_var);
                    case 2 % divide by fatigue
                        benefit_var = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2))./(1 + kTreward.*fatigue_var);
                    case 3 % multiply by fatigue
                        benefit_var = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2)).*(1 - kTreward.*fatigue_var);
                end
                expected_payoff_var = kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2);

            case {'inc','absInc'}
                switch prm_B_time_on_benef
                    case 0 % no rest effect
                        benefit_var = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) + kI_bis.*inc_var_bis);
                    case 1 % subtract by fatigue
                        benefit_var = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) + kI_bis.*inc_var_bis - kTreward.*fatigue_var);
                    case 2 % divide by fatigue
                        benefit_var = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) + kI_bis.*inc_var_bis)./(1 + kTreward.*fatigue_var);
                    case 3 % multiply by fatigue
                        benefit_var = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) + kI_bis.*inc_var_bis).*(1 - kTreward.*fatigue_var);
                end
                expected_payoff_var= kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) + kI_bis.*inc_var_bis;
            case {'absInc_plus_nomInc','absInc_plus_cond',...
                    'inc_perCond'}
                switch prm_B_time_on_benef
                    case 0 % no rest effect
                        benefit_var = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) + kI_bis1.*inc_var_bis1 + kI_bis2.*inc_var_bis2);
                    case 1 % subtract by fatigue
                        benefit_var = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) + kI_bis1.*inc_var_bis1 + kI_bis2.*inc_var_bis2 - kTreward.*fatigue_var);
                    case 2 % divide by fatigue
                        benefit_var = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) + kI_bis1.*inc_var_bis1 + kI_bis2.*inc_var_bis2)./(1 + kTreward.*fatigue_var);
                    case 3 % multiply by fatigue
                        benefit_var = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) + kI_bis1.*inc_var_bis1 + kI_bis2.*inc_var_bis2).*(1 - kTreward.*fatigue_var);
                end
                expected_payoff_var= kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) + kI_bis1.*inc_var_bis1 + kI_bis2.*inc_var_bis2;
        end % incentive bis
end

%% X based on the data
if ismember(prm_B_perf_C_perf_or_force,...
        {'perf_f_X','perf_f_X_bis','perf_f_X_ter',...
        'perf_f_X_4','perf_f_X_5','perf_f_X_6','perf_f_X_7'})

    %% define PM value
    if ismember(prm_B_perf_C_perf_or_force,{'perf_f_X_ter','perf_f_X_7'}) % PM fixed at 1 for both tasks
        PM = 1;

    elseif strcmp(prm_B_perf_C_perf_or_force,'perf_f_X_5') % PM fixed at 0.75 for grip and 1 for stroop
        switch task_nm
            case 'grip'
                PM = 0.75;
            case 'stroop'
                PM = 1;
        end

    else
        switch task_nm
            case 'grip'
                switch prm_B_perf_C_perf_or_force
                    case {'perf_f_X','perf_f_X_bis','perf_f_X_6'}
                        PM = Fmax_estim_run./(Fmax*Fmax_scale);
                    case 'perf_f_X_4'
                        PM = Fmax_estim_run./Fmax;
                end
            case 'stroop'
                PM = Fmax_estim_run;
        end

        % check PM never higher than 1
        for iT = 1:length(PM)
            if PM(iT) > 1 % max perf can not be higher than 1
                PM(iT) = 1;
            end
        end

    end

    if ismember(prm_B_perf_C_perf_or_force,{'perf_f_X','perf_f_X_ter',...
            'perf_f_X_4','perf_f_X_5'})
        E = -log(1 - perf./PM)./kX; % can use X_pred instead as well

    elseif strcmp(prm_B_perf_C_perf_or_force,'perf_f_X_bis')
        error('not ready');

    elseif ismember(prm_B_perf_C_perf_or_force,{'perf_f_X_6','perf_f_X_7'})
        E = kX./( (PM./perf) - 1);
    end

else
    E = [];
end

%% expected value based on actual performance
perf_for_EV = perf; % perf_pred or perf you have to chose which one you want to use for EV computation
switch prm_B_perf_C_perf_or_force
    case 'perf'
        benefit = benefit_var.*((perf_for_EV.*(GLcond==1)) + (1-perf_for_EV).*(GLcond==0));
        cost = cost_var.*(perf_for_EV./(kmax_run - perf_for_EV));
        expected_payoff = expected_payoff_var.*((perf_for_EV.*(GLcond==1)) + (1-perf_for_EV).*(GLcond==0));
    case 'B_perf_C_force'
        benefit = benefit_var.*((perf_for_EV.*(GLcond==1)) + (1-perf_for_EV).*(GLcond==0));
        cost = cost_var.*(perf_for_EV./(kmax_run - perf_for_EV));
        expected_payoff = expected_payoff_var.*((perf_for_EV.*(GLcond==1)) + (1-perf_for_EV).*(GLcond==0));
    case {'perf_f_X','perf_f_X_ter','perf_f_X_4','perf_f_X_5',...
            'perf_f_X_6','perf_f_X_7'}
        benefit = benefit_var.*((perf_for_EV.*(GLcond==1)) + (1-perf_for_EV).*(GLcond==0));
        cost = cost_var.*E;
        expected_payoff = expected_payoff_var.*((perf_for_EV.*(GLcond==1)) + (1-perf_for_EV).*(GLcond==0));
    case 'perf_f_X_bis'
        benefit = benefit_var.*((perf_for_EV.*(GLcond==1)) + (1-perf_for_EV).*(GLcond==0));
        cost = cost_var.*(E.^2);
        expected_payoff = expected_payoff_var.*((perf_for_EV.*(GLcond==1)) + (1-perf_for_EV).*(GLcond==0));
end
EV = benefit - cost;

%% expected value based on predicted performance
switch prm_B_perf_C_perf_or_force
    case 'perf'
        benefit_pred = benefit_var.*((perf_pred.*(GLcond==1)) + (1-perf_pred).*(GLcond==0));
        cost_pred = cost_var.*(perf_pred./(kmax_run - perf_pred));
        expected_payoff_pred = expected_payoff_var.*((perf_pred.*(GLcond==1)) + (1-perf_pred).*(GLcond==0));
    case 'B_perf_C_force'
        benefit_pred = benefit_var.*((perf_pred.*(GLcond==1)) + (1-perf_pred).*(GLcond==0));
        cost_pred = cost_var.*(perf_pred./(kmax_run - perf_pred));
        expected_payoff_pred = expected_payoff_var.*((perf_pred.*(GLcond==1)) + (1-perf_pred).*(GLcond==0));
    case {'perf_f_X','perf_f_X_ter',...
            'perf_f_X_4','perf_f_X_5','perf_f_X_6','perf_f_X_7'}
        benefit_pred = benefit_var.*((perf_pred.*(GLcond==1)) + (1-perf_pred).*(GLcond==0));
        cost_pred = cost_var.*E_pred;
        expected_payoff_pred = expected_payoff_var.*((perf_pred.*(GLcond==1)) + (1-perf_pred).*(GLcond==0));
    case 'perf_f_X_bis'
        benefit_pred = perf_pred + benefit_var.*((perf_pred.*(GLcond==1)) + (1-perf_pred).*(GLcond==0));
        cost_pred = cost_var.*(E_pred.^2);
        expected_payoff_pred = expected_payoff_var.*((perf_pred.*(GLcond==1)) + (1-perf_pred).*(GLcond==0));
end
EV_pred = benefit_pred - cost_pred;

end % function