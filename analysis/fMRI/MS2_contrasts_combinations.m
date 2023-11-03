function [con_names, con_vec] = MS2_contrasts_combinations(GLMprm, con_names, con_vec)
%[con_names, con_vec] = MS2_contrasts_combinations(GLMprm, prm_idx, con_names, con_vec);
% MS2_contrasts_combinations extracts the names and vectors constructions
% for the contrasts which need to combine multiple regressors (sum and/or
% comparison).
%
% INPUTS
% GLMprm: structure with GLM parameters
%
% con_names: list of contrast names in the order they will be applied at
% the first and second level.
%
% con_vec: contrast vectors to be used at the first level, their order
% follows con_names
%
% OUTPUTS
% con_names: list of contrast names in the order they will be applied at
% the first and second level after adding composite contrasts (based on
% combinations of regressors)
%
% con_vec: contrast vectors to be used at the first level, their order
% follows con_names after adding composite contrasts (based on
% combinations of regressors)
%
% prm_idx: structure with corresponding index for each regressor
%

posCon_nm = '_pos';

%% extract main relevant parameters for the script to work
% reinforcement-learning task
RLprm = GLMprm.RL;
% grip task
gripRPprm = GLMprm.grip;
% stroop task
stroopRPprm = GLMprm.stroop;

n_regs = size(con_vec,2);

%% contrast index
jCon = length(con_names); % start after the last contrast based on the single regressors

%% grip task
[ jCon, con_vec, con_names ] = MS2_contrasts_combinations_Grip( jCon, con_vec, con_names,...
    posCon_nm, n_regs, gripRPprm);

%% stroop task
[ jCon, con_vec, con_names ] = MS2_contrasts_combinations_Stroop( jCon, con_vec, con_names,...
    posCon_nm, n_regs, stroopRPprm);

%% grip and stroop pooled
[ jCon, con_vec, con_names ] = MS2_contrasts_combinations_GS( jCon, con_vec, con_names,...
    posCon_nm, n_regs, gripRPprm, stroopRPprm);

%% reinforcement-learning task
[ jCon, con_vec, con_names ] = MS2_contrasts_combinations_RL( jCon, con_vec, con_names,...
    posCon_nm, n_regs, RLprm);

%% grip + stroop + RL contrasts
[ ~, con_vec, con_names ] = MS2_contrasts_combinations_GSL( jCon, con_vec, con_names,...
    posCon_nm, n_regs, RLprm, gripRPprm, stroopRPprm);

end % function