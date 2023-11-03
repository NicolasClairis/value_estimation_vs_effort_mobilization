function [gx] = MS2_GS_g_observation_perf_optimal(x_t, P, u_t, inG)
%[gx] = MS2_GS_g_observation_perf_optimal(x_t, P, u_t, inG)
% Predict performance for grip and stroop tasks.
%
% INPUT
% x_t : []
%
% P=parameters to fit
%
%
% u_t : [incentive; fatigue var (trial number/total gain)/time rest/run index]
%
% inG: structure which informs about which model to use (formula, variables to include,
% etc.)
%
% OUTPUT
% - gx : Scalar
%
% See also MS2_GS_perf_optimal_level

%% define parameters
prm_priors          = inG.priors;
C_fatigue           = inG.C_fatigue;
C_rest              = inG.C_rest;
B_time_on_benef     = inG.B_time_on_benef;
kmax_fixed_or_free  = inG.C_kmax_fixed_or_free;
incentive_var       = inG.incentive_var;
incentive_var_bis   = inG.incentive_var_bis;
fatigue_var_nm      = inG.fatigue_var;
task_nm             = inG.task_nm;
Fmax_fixed_or_free  = inG.Fmax_fixed_or_free;
B_perf_C_perf_or_force = inG.B_perf_C_perf_or_force;
run_Fmax            = inG.run_Fmax;

switch task_nm
    case 'grip'
        Fmax_scalingF = 100/75; % to compute performance, Fmax(t) is multiplied by 100/75 for grip
end

u_var_idx = 0;
% incentive
switch incentive_var
    case {'inc','absInc'}
        u_var_idx = u_var_idx + 1;
        inc_var         = u_t(u_var_idx);
    case {'absInc_plus_nomInc','absInc_plus_cond',...
            'inc_perCond'}
        u_var_idx = u_var_idx + 1;
        inc_var1        = u_t(u_var_idx);
        u_var_idx = u_var_idx + 1;
        inc_var2        = u_t(u_var_idx);
end
% incentive bis
switch incentive_var_bis
    case {'inc','absInc'}
        u_var_idx = u_var_idx + 1;
        inc_var_bis = u_t(u_var_idx);
    case {'absInc_plus_nomInc','absInc_plus_cond',...
            'inc_perCond'}
        u_var_idx = u_var_idx + 1;
        inc_var_bis1 = u_t(u_var_idx);
        u_var_idx = u_var_idx + 1;
        inc_var_bis2 = u_t(u_var_idx);
end

% fatigue
if ~strcmp(fatigue_var_nm,'')
    u_var_idx = u_var_idx + 1;
    fatigue_var = u_t(u_var_idx);
end

% time rest (always in u, even if not in the model)
u_var_idx = u_var_idx + 1;
timeRest_idx    = u_var_idx;
% include rest in the model
if C_rest > 0
    time_rest_var   = u_t(timeRest_idx);
end

% include rest in the model
u_var_idx = u_var_idx + 1;
runIdx_idx      = u_var_idx;
run_idx = u_t(runIdx_idx);

% compute Fmax rather than a fixed (always in u, even if not in the model)
u_var_idx = u_var_idx + 1;
Fmax_idx        = u_var_idx;
% Fmax trial by trial
Fmax_var = u_t(Fmax_idx);
run_Fmax_var = run_Fmax(run_idx);

% extract if gain or loss trial
u_var_idx = u_var_idx + 1;
GL_idx        = u_var_idx;
gain_or_loss = u_t(GL_idx);
% multiply by (+1) or (-1) depending if gain or loss (to account for the
% fact that monetary incentives vary as (1-P) in loss and as (P) in gains)
switch gain_or_loss
    case 0 % loss
        GL_term = -1;
    case 1 % gain
        GL_term = 1;
end
switch incentive_var
    case {'inc'}
        inc_var = inc_var.*GL_term;
    otherwise
        error('please fix the sign and do again later');
end
switch incentive_var_bis
    case {'inc'}
        inc_var_bis = inc_var_bis.*GL_term;
    otherwise
        error('please fix the sign and do again later');
end

%% parameters to fit (force them all to be superior to zero => adapt formula
% accordingly in gx)
% cost
iPrm = 1;
kCost       = fn_for_prior(P(iPrm), prm_priors);
% incentive
switch incentive_var
    case {'inc','absInc'}
        iPrm = iPrm + 1;
        kI          = fn_for_prior(P(iPrm), prm_priors);
    case {'absInc_plus_nomInc','absInc_plus_cond',...
            'inc_perCond'}
        iPrm = iPrm + 1;
        kI1          = fn_for_prior(P(iPrm), prm_priors);
        iPrm = iPrm + 1;
        kI2          = fn_for_prior(P(iPrm), prm_priors);
end


% demotivation effect on benefit term
if B_time_on_benef > 0
    iPrm = iPrm + 1;
    kTreward    = fn_for_prior(P(iPrm), prm_priors);
end

% incentive bis
switch incentive_var_bis
    case {'inc','absInc'}
        iPrm = iPrm + 1;
        kI_bis = fn_for_prior(P(iPrm), prm_priors);
    case {'absInc_plus_nomInc','absInc_plus_cond',...
            'inc_perCond'}
        iPrm = iPrm + 1;
        kI_bis1          = fn_for_prior(P(iPrm), prm_priors);
        iPrm = iPrm + 1;
        kI_bis2          = fn_for_prior(P(iPrm), prm_priors);
end

% fatigue effect on cost term
if C_fatigue > 0
    iPrm = iPrm + 1;
    if C_fatigue ~= 3
        kTcost      = fn_for_prior(P(iPrm), prm_priors);
    elseif C_fatigue == 3 % log(kT*T) => has to be positive
        kTcost      = fn_for_prior(P(iPrm), 'pos');
    end
end

% rest effect on cost term
if C_rest > 0
    iPrm = iPrm + 1;
    kRest       = fn_for_prior(P(iPrm), prm_priors);
end

% kmax/run
if strcmp(kmax_fixed_or_free,'free')
    switch run_idx
        case 1
            iPrm = iPrm + 1;
            kmax_estim_run = 1 + exp(P(iPrm)); % kmax>1
        case 2
            iPrm = iPrm + 2;
            kmax_estim_run = 1 + exp(P(iPrm)); % kmax>1
    end
end

% Fmax/run
switch Fmax_fixed_or_free
    case 'fixed'
        Fmax_estim_run = run_Fmax_var;
    case 'free'
        switch run_idx
            case 1
                iPrm = iPrm + 1;
                Fmax_estim_run = P(iPrm);
            case 2
                iPrm = iPrm + 2;
                Fmax_estim_run = P(iPrm);
        end
    case 'free_bis'
        switch run_idx
            case 1
                iPrm = iPrm + 1;
                Fmax_estim_run = (1 + sigmo(P(iPrm))).*run_Fmax_var; % Fmax=kF*runFmax with kF > 1
            case 2
                iPrm = iPrm + 2;
                Fmax_estim_run = (1 + sigmo(P(iPrm))).*run_Fmax_var; % Fmax=kF*runFmax with kF > 1
        end
end
if ismember(B_perf_C_perf_or_force,...
        {'perf_f_X','perf_f_X_bis','perf_f_X_ter',...
        'perf_f_X_4','perf_f_X_5','perf_f_X_6','perf_f_X_7'})
    iPrm = iPrm + 1;
    kX = fn_for_prior(P(iPrm), 'pos');
end

%% compute output
% cost term
switch C_fatigue
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
switch C_rest
    case 0 % no rest effect
        cost    = kCost.*C_fatigue_effect;
    case 1 % subtract by rest
        cost    = kCost.*(C_fatigue_effect - kRest.*time_rest_var);
    case 2 % divide by rest
        cost    = kCost.*C_fatigue_effect./(1 + kRest.*time_rest_var);
    case 3 % multiply by rest
        cost    = kCost.*C_fatigue_effect.*(1 - kRest.*time_rest_var); %% multiplication = weird, should check range of values
end

% benefit term
switch incentive_var
    case {'inc','absInc'}
        switch incentive_var_bis
            case ''
                switch B_time_on_benef
                    case 0 % no rest effect
                        benefit = (1 +kI.*inc_var);
                    case 1 % subtract by fatigue
                        benefit = (1 +kI.*inc_var - kTreward.*fatigue_var);
                    case 2 % divide by fatigue
                        benefit = (1 +kI.*inc_var)./(1 + kTreward.*fatigue_var);
                    case 3 % multiply by fatigue
                        benefit = (1 +kI.*inc_var).*(1 - kTreward.*fatigue_var);
                end
            case {'inc','absInc'}
                switch B_time_on_benef
                    case 0 % no rest effect
                        benefit = (1 +kI.*inc_var + kI_bis.*inc_var_bis);
                    case 1 % subtract by fatigue
                        benefit = (1 +kI.*inc_var + kI_bis.*inc_var_bis - kTreward.*fatigue_var);
                    case 2 % divide by fatigue
                        benefit = (1 +kI.*inc_var + kI_bis.*inc_var_bis)./(1 + kTreward.*fatigue_var);
                    case 3 % multiply by fatigue
                        benefit = (1 +kI.*inc_var + kI_bis.*inc_var_bis).*(1 - kTreward.*fatigue_var);
                end
            case {'absInc_plus_nomInc','absInc_plus_cond',...
                    'inc_perCond'}
                switch B_time_on_benef
                    case 0 % no rest effect
                        benefit = (1 + kI1.*inc_var_bis1 + kI2.*inc_var_bis2);
                    case 1 % subtract by fatigue
                        benefit = (1 + kI1.*inc_var_bis1 + kI2.*inc_var_bis2 - kTreward.*fatigue_var);
                    case 2 % divide by fatigue
                        benefit = (1 + kI1.*inc_var_bis1 + kI2.*inc_var_bis2)./(1 + kTreward.*fatigue_var);
                    case 3 % multiply by fatigue
                        benefit = (1 + kI1.*inc_var_bis1 + kI2.*inc_var_bis2).*(1 - kTreward.*fatigue_var);
                end
        end % incentive bis
        % expected benefit weighted by expected Fmax on the effort scale
        if strcmp(task_nm,'grip') &&...
                ismember(incentive_var,{'inc_GexpPerf','absInc_GexpPerf'})
            benefit = benefit./expPerfMax;
        end

    case {'absInc_plus_nomInc','absInc_plus_cond',...
            'inc_perCond'}
        switch incentive_var_bis
            case ''
                switch B_time_on_benef
                    case 0 % no rest effect
                        benefit = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2));
                    case 1 % subtract by fatigue
                        benefit = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) - kTreward.*fatigue_var);
                    case 2 % divide by fatigue
                        benefit = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2))./(1 + kTreward.*fatigue_var);
                    case 3 % multiply by fatigue
                        benefit = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2)).*(1 - kTreward.*fatigue_var);
                end
            case {'inc','absInc'}
                switch B_time_on_benef
                    case 0 % no rest effect
                        benefit = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) + kI_bis.*inc_var_bis);
                    case 1 % subtract by fatigue
                        benefit = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2 + kI_bis.*inc_var_bis) - kTreward.*fatigue_var);
                    case 2 % divide by fatigue
                        benefit = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) + kI_bis.*inc_var_bis)./(1 + kTreward.*fatigue_var);
                    case 3 % multiply by fatigue
                        benefit = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) + kI_bis.*inc_var_bis).*(1 - kTreward.*fatigue_var);
                end
            case {'absInc_plus_nomInc','absInc_plus_cond',...
                    'inc_perCond'}
                switch B_time_on_benef
                    case 0 % no rest effect
                        benefit = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) + kI_bis1.*inc_var_bis1 + kI_bis2.*inc_var_bis2);
                    case 1 % subtract by fatigue
                        benefit = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) + kI_bis1.*inc_var_bis1 + kI_bis2.*inc_var_bis2 - kTreward.*fatigue_var);
                    case 2 % divide by fatigue
                        benefit = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) + kI_bis1.*inc_var_bis1 + kI_bis2.*inc_var_bis2)./(1 + kTreward.*fatigue_var);
                    case 3 % multiply by fatigue
                        benefit = (1 +kI1.*(inc_var1.^kaI1) + kI2.*(inc_var2.^kaI2) + kI_bis1.*inc_var_bis1 + kI_bis2.*inc_var_bis2).*(1 - kTreward.*fatigue_var);
                end
        end % incentive bis
end % incentive

%% extract Fmax and runFmax data
if strcmp(B_perf_C_perf_or_force,'B_perf_C_force')
    run_maxForce_term = Fmax_estim_run;
    Force_term = run_maxForce_term./Fmax_var;
end

%% run perf max (based on fitted Fmax)
if ismember(B_perf_C_perf_or_force,...
        {'perf_f_X','perf_f_X_bis','perf_f_X_4','perf_f_X_6'})

    switch task_nm
        case 'grip'

            switch B_perf_C_perf_or_force
                case {'perf_f_X','perf_f_X_bis','perf_f_X_6'}

                    if Fmax_estim_run <= (Fmax_var*Fmax_scalingF) % perf max depends on the relation between force capacity and max of the scale (for each trial)
                        Pmax_estim_run = Fmax_estim_run./(Fmax_var*Fmax_scalingF);
                    else % perf max = 1 if capacity is higher than what the scale allows
                        Pmax_estim_run = 1;
                    end

                case {'perf_f_X_4'} % no scaling of Fmax

                    if Fmax_estim_run <= Fmax_var % perf max depends on the relation between force capacity and max of the scale (for each trial)
                        Pmax_estim_run = Fmax_estim_run./Fmax_var;
                    else % perf max = 1 if capacity is higher than what the scale allows
                        Pmax_estim_run = 1;
                    end

            end

        case 'stroop'
            if Fmax_estim_run <= 1
                Pmax_estim_run = Fmax_estim_run;
            else % perf max = 1 if capacity is higher than what the scale allows
                Pmax_estim_run = 1;
            end
    end

elseif ismember(B_perf_C_perf_or_force,{'perf_f_X_ter','perf_f_X_7'}) % Pmax fixed at 1 (=useless)
    Pmax_estim_run = 1;

elseif strcmp(B_perf_C_perf_or_force,'perf_f_X_5') % Pmax fixed at 1 for stroop and 0.75 for grip (=useless)
    switch task_nm
        case 'grip'
            Pmax_estim_run = 0.75;
        case 'stroop'
            Pmax_estim_run = 1;
    end
end

%% constraints on cost and benefit terms

switch prm_priors
    case 'pos' % when priors are constrained to be positive, cost and benefit will be positive by default => no need of additional constraint in principle

        % cost term has to be positive
        % for cases where rest effect is too strong compared to fatigue effect (and
        % subtraction)
%         if cost <= 0
%             cost = eps;
%         end

%         % benefit term has to be positive
%         % for cases where fatigue effect is too strong compared to reward effect (and
%         % subtraction)
%         if benefit <= 0
%             benefit = eps;
%         end

    case 'free' % no constraint on the priors to claim we found the good ones => only constraint on benefit to avoid dividing by zero

        % avoid dividing by zero => avoid point where benefits are equal
        % zero
        if benefit == 0
            benefit = eps;
        end

        if strcmp(B_perf_C_perf_or_force,'perf_f_X_bis') && cost == 0
            cost = eps;
        end

        % force cost and benefit to be positive
        if cost < 0 %|| benefit < 0
            cost = eps;
%             benefit = eps;
        end
end

%% make the fit for force
if strcmp(B_perf_C_perf_or_force,'B_perf_C_force') &&...
        strcmp(kmax_fixed_or_free,'fixed')
    gx = Force_term - sqrt( (cost./benefit).*Force_term );

elseif ismember(B_perf_C_perf_or_force,{'perf_f_X','perf_f_X_ter','perf_f_X_4','perf_f_X_5'})
    gx = Pmax_estim_run - cost./(benefit*kX);

elseif ismember(B_perf_C_perf_or_force,{'perf_f_X_6','perf_f_X_7'})
    gx = Pmax_estim_run.*(1 - sqrt( (kX*cost)./(benefit*Pmax_estim_run )) );

elseif strcmp(B_perf_C_perf_or_force,'perf_f_X_bis')
    gx = Pmax_estim_run.*(1 - exp( -lambertw( -(kX^2).*benefit.*Pmax_estim_run./(2*cost)) ) );
    %     error('weird values, with this model maybe needs fixing');

else
    switch kmax_fixed_or_free
        case 'fixed'
            gx = 1 - sqrt( cost./benefit ); % mathematically [1 + sqrt(-cost./benefit)] could also be a solution
            % but the perf has to be < 1 since can not be > 1 so only this solution
            % can be considered

        case 'free'
            gx = kmax_estim_run - sqrt( cost./benefit );
    end
end

%% check weird gx values
if isnan(gx) || isinf(gx) || ~isreal(gx)
    error(['gx = ',num2str(gx),', please fix it']);
end

end % function