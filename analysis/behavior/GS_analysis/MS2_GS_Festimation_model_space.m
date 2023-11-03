function [ prm_in_mdl ] = MS2_GS_Festimation_model_space()
%[ prm_in_mdl ] = MS2_GS_Festimation_model_space()
%
% OUTPUTS
%prm_in_mdl: structure with all the parameters/model to know the
%characteristics of all the models included in the model comparison.
%
%   parameters:
%   .norm = normalize regressors or not to force them to be <= 1
%           (0) not-normalized
%           (1) normalized variables
%
%   .priors = decide whether you want to add constraints to the priors or
%   not
%       'free': values for fitted parameters (except when really necessary) can take all values from -Inf to +Inf
%       'pos': all the parameters are constrained to be positive (makes
%       sense if your goal is to find the best possible parameters and not
%       to see whether they are significant)
%
%   .incentive_type = incentive type for incentive_var (see below)
%           'value': actual value in euros (0.01€ to 20€)
%           'rank': rank of the incentive (1 to 6)
%
%   .incentive_var =
%           - 'absInc' motivational (absolute) incentive
%           - 'inc': nominal incentive (linear values)
%           - 'absInc_plus_nomInc': motivational incentive + nominal
%           incentive = linear value (|I|/I), 2 parameters
%           - 'absInc_plus_cond': motivational incentive + 1 binary variable
%           to know which condition it is (gain/loss)
%           - 'inc_perCond': 1 incentive parameter per condition (gain/loss)
%
%   .incentive_type_bis = incentive type for incentive_var_bis (see below)
%           'value': actual value in euros (0.01€ to 20€)
%           'rank': rank of the incentive (1 to 6)
%
%   .incentive_var_bis =
%           - 'absInc' motivational (absolute) incentive
%           - 'inc': nominal incentive (linear values)
%           - 'absInc_plus_nomInc': motivational incentive + nominal
%           incentive = linear value (|I|/I), 2 parameters
%           - 'absInc_plus_cond': motivational incentive + 1 binary variable
%           to know which condition it is (gain/loss)
%           - 'inc_perCond': 1 incentive parameter per condition (gain/loss)
%
%   .fatigue_var:
%       '': no fatigue term used
%       'trialN': trial number
%       'totalGainPrev': total cumulated gains
%       'sumPerfPrev': sum previous performance
%
%   .C_fatigue: fatigue effect on cost term
%       (0) no fatigue term inside cost term
%       (1) (1+kT*T) fatigue effect on cost term
%       (2) (1+exp(kT*T)) fatigue effect on cost term
%       (3) (1+log(1+kT*T) fatigue effect on cost term (requires kT>0)
%       (4) (1-exp(kT*T)) fatigue effect on cost term
%
%   .C_rest: rest effect on cost term
%       (0) no rest effect on cost term
%       (1) subtract fixation cross time (=rest time) to fatigue effect
%
%   .B_time_on_benef: benefit - time effect on benefit
%       (0) no time effect on benefit term
%       (1) subtract benefit term by fatigue term
%       (2) divide benefit term by fatigue term
%
%   .Fmax_fixed_or_free:
%       'fixed': use an (arbitrary) Fmax a little higher than the measured
%       Fmax and fix it
%       'free': estimate the Fmax of a given run for a given subject with
%       the model rather than using a predetermined value (=> this may be
%       penalized by model comparison since you add a free parameter)
%       'free_bis': Fmax = runFmax*k where k>1 so that you cannot have
%       a fitted Fmax lower than the runFmax
%
%   .C_kmax_fixed_or_free: determine value of kmax in Perf/(kmax-Perf) for
%   perf model
%       'fixed': fix kmax at 1
%       'free': estimate the kmax for each subject and each run (=> this may be
%       penalized by model comparison since you add a free parameter)
%
%   .B_perf_C_perf_or_force
%       'perf': use performance
%       'B_perf_C_force': model with benefit term depending solely on
%       performance and cost term depending solely on force exerted
%       'perf_f_X': model with performance term depending on a fictive
%       ressource X which investment has to be minimized to make the task
%       (=> EV derived based on X and not on P in that case)
%       'perf_f_X_bis': same but X^2 with cost instead of linear
%       'perf_f_X_ter': same as perf_f_X but Pmax erased from the model
%       (always equal to 1)
%       'perf_f_X_4': same as perf_f_X but Pmax replaced by runFmax/Fmax(t)
%       for grip (ignoring scaling factor) (no change for stroop)
%       'perf_f_X_5': same as perf_f_X but Pmax replaced by 0.75 for grip
%       (and 1 for Stroop which means it's like perf_f_X_ter for stroop)
%       'perf_f_X_6': different formula P(X) = PM*(X/(k+X)) instead of
%       P(X)=PM*(1-exp(-k*X)) like previous models. PM=runFmax/Fmax(t)
%       'perf_f_X_7': same as 'perf_f_X_6' but PM = 1
%       
%   .multisess
%       '': regroup all data together
%       'multisess_singlePrmSet': add a multisession factor (data of each
%       session might be considered with different levels of noise) but
%       will keep one single set of parameters across sessions for each
%       participant

%% model 1
mdl_nm = 'model_1';
prm_in_mdl.(mdl_nm).norm = 1;
prm_in_mdl.(mdl_nm).priors = 'pos';
prm_in_mdl.(mdl_nm).incentive_type = 'value';
prm_in_mdl.(mdl_nm).incentive_var = 'inc';
prm_in_mdl.(mdl_nm).incentive_type_bis = '';
prm_in_mdl.(mdl_nm).incentive_var_bis = '';
prm_in_mdl.(mdl_nm).fatigue_var = 'trialN';
prm_in_mdl.(mdl_nm).C_fatigue = 1;
prm_in_mdl.(mdl_nm).C_rest = 0; % no rest effect on cost term
prm_in_mdl.(mdl_nm).B_time_on_benef = 0; % no time effect on benefit term
prm_in_mdl.(mdl_nm).Fmax_fixed_or_free = '';
prm_in_mdl.(mdl_nm).C_kmax_fixed_or_free = 'fixed';
prm_in_mdl.(mdl_nm).B_perf_C_perf_or_force = 'perf_f_X_7';
prm_in_mdl.(mdl_nm).multisess = 'multisess_singlePrmSet';
end % function