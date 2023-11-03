function [ ] = which_GLM_MS2_infos_gal( GLMprm )
%which_GLM_MS2_infos_gal( GLMprm ) outputs the informations about the GLM based
% on the list of parameters designated inside GLMprm for the
% parameters global for all tasks
%
% INPUTS
% GLMprm: structure with GLM parameters for each task extracted with
% which_GLM_MS2.
%
% See also which_GLM_MS2

gal     = GLMprm.gal;

%% main parameters info
disp('***main parameters');

%% grey mask
switch gal.grey_mask
    case 1
        disp('Use of 1st level probability grey mask for each subject.');
    case 2
        disp('Use of 1st level probability grey mask based on the average T1 of all subjects.');
    case {3,4}
        disp('Use of >10% probability of being in grey matter mask based on SPM.');
    case 5
        disp('Use of exclusive ventricle mask based on SPM.');
end

%% derivative
switch gal.add_drv
    case 1
        disp('Use of temporal derivative => be aware that each regressor will be doubled in the matrix list.');
    case 2
        disp('Use of temporal and spatial derivative => be aware that each regressor will be tripled in the matrix list');
end

%% variables orthogonalized or not? 
switch gal.orth_vars
    case 0
        disp('variables not orthogonalized in the model.');
    case 1
        disp('all variables orthogonalized');
end

%% 1 regressor/trial
switch gal.onsets_only
    case 1
        disp('onsets_only: 1 regressor/trial');
end

%% FIR model
switch gal.FIR
    case 1
        disp('FIR model');
end

%% zscore
switch gal.zscore
    case 0
        disp('no zcore/raw values');
    case 1
        disp('all regressors zscored per run');
end

end % function