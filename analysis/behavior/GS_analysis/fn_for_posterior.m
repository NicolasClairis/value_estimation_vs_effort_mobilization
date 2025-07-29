function[muPosterior, varPosterior] = fn_for_posterior(muX, varX, prm_priors)
%[muPosterior, varPosterior] = fn_for_posterior(muX, varX, prm_priors)
% fn_for_posterior recovers the parameter of interest based on the mean
% (muX) and variance (varX) produced by the VBA for the parameter
% provided in the inputs and the presence (or not) of a positivity
% constraint to force the parameter to be within a certain range.
% If there was no positivity constraint, then the parameter is the same
%
% INPUTS
% muX: mean of the posterior obtained with the VBA toolbox
%
% varX: variance of the posterior obtained with the VBA toolbox
%
% prm_priors = parameter for the prior to know how to transform X
%   'free': no constraint on the prior => keeps param = X
%   'pos': prior is constrained to be positive using param K = exp(X)
%
% OUTPUTS
% muPosterior: mean of the posterior parameter, after taking into account the (potential)
% positivity constraint
% 
% varPosterior: variance of the posterior parameter, after taking into account
% the (potential) positivity constraint
%
% See also fn_for_prior.m and
% https://mbb-team.github.io/VBA-toolbox/wiki/param-transform/#positivity-constraint
% for details.
%
% Note: careful variance = (std)^2 so be careful when you use these values
% to not make an error on which you are using and reporting.

switch prm_priors
    case 'free' % no constraint => prm = X
        muPosterior = muX;
        varPosterior = varX;
    case 'pos' % positive constraint on prior: prm=exp(X)
        % formula detailled here: https://mbb-team.github.io/VBA-toolbox/wiki/param-transform/#positivity-constraint
        muPosterior = exp(muX + (varX/2));
        varPosterior = exp(2*muX + varX)*(exp(varX) - 1);
end

end % function