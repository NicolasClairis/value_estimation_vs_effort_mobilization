function [m_prm, prm] = MS2_GS_perf_optimal_level( prm_in_mdl, sub_nm )
%[m_prm, prm_qlty] = MS2_GS_perf_optimal_level( prm_in_mdl, sub_nm )
% estimates parameters for grip and stroop tasks for each subject with VBA
% toolbox.
%
% INPUTS
% prm_in_mdl: structure with info about which model to use
%   see MS2_GS_Festimation_model_space.m for the details
%
% sub_nm: '' or 'sXX_DDMMYY' if left empty, computes for all subjects ('fMRI'), otherwise,
% will only fit the model on the designated subject
%
% OUTPUT
% m_prm: structure with model quality and posteriors of the different
% parameters averaged across subjects subject
% NOTE: if only one subject entered in the input, then m_prm will only
% reflect this subject's parameters and not the average across subjects..
%
% prm: structure with model quality and posteriors of the different
% parameters for each subject
%
% See also MS2_GS_g_observation_perf_optimal, MS2_GS_Festimation_model_space.m

%% working directories
root = 'enter path here';

%% check if VBA toolbox is in the path already or not
if exist('VBA_NLStateSpaceModel.m','file') == 0
    error(['Please install the VBA toolbox at ',...
        'https://mbb-team.github.io/VBA-toolbox/download/ to be able to ',...
        'apply the RFT on the pupil data.']);
end

%% subject list
if ~exist('sub_nm','var') || isempty(sub_nm)
    subject_id = {'enter list of subjects here'};
    NS = length(subject_id);
else
    subject_id = {sub_nm};
    NS = 1;
end

%% load parameters of the model
prm_norm = prm_in_mdl.norm;
prm_priors = prm_in_mdl.priors;
prm_incentive_type = prm_in_mdl.incentive_type;
prm_incentive_var = prm_in_mdl.incentive_var;
prm_incentive_type_bis = prm_in_mdl.incentive_type_bis;
prm_incentive_var_bis = prm_in_mdl.incentive_var_bis;
prm_kmax_fixed_or_free = prm_in_mdl.C_kmax_fixed_or_free;
prm_fatigue_var = prm_in_mdl.fatigue_var;
prm_B_time_on_benef = prm_in_mdl.B_time_on_benef;
prm_C_fatigue = prm_in_mdl.C_fatigue;
prm_C_rest = prm_in_mdl.C_rest;
prm_Fmax_fixed_or_free = prm_in_mdl.Fmax_fixed_or_free;
prm_B_perf_C_perf_or_force = prm_in_mdl.B_perf_C_perf_or_force;
prm_multisess = prm_in_mdl.multisess;

%% load variables of interest
nTrials_per_run = 60;
n_GS_runs = 2;
task_names = {'grip','stroop'};
nTasks = length(task_names);
n_Stroop_maxPairs = 10; % Fmax possible for stroop
grip_scalingFactor = 100/75;

for iGS = 1:nTasks
        task_nm = task_names{iGS};
        
        [kCost.(task_nm),...
            kTcost.(task_nm),...
            kRest.(task_nm),...
            kI.(task_nm),...
            kI1.(task_nm),...
            kI2.(task_nm),...
            kI3.(task_nm),...
            kI_bis.(task_nm),...
            kI_bis1.(task_nm),...
            kI_bis2.(task_nm),...
            kTreward.(task_nm),...
            kX.(task_nm),...
            prm_qlty.(task_nm).LL,...
            prm_qlty.(task_nm).F,...
            prm_qlty.(task_nm).R2] = deal( NaN(1, NS) );
        [kmax_estim.(task_nm),...
            Fmax_estim.(task_nm),...
            kFmax_estim.(task_nm)] = deal(NaN(2,NS)); % 1/run
    [perf_predicted, E_predicted] = deal(NaN(NS, nTrials_per_run*n_GS_runs));
    
    for iS = 1:NS
        
        % subject identification
        sub_nm = subject_id{iS};
        if strcmp(sub_nm(3),'_')
            subid   = sub_nm(2);
        elseif ~strcmp(sub_nm(3),'_') && strcmp(sub_nm(4),'_')
            subid = sub_nm(2:3);
        end
        
        onsets_folder = [root,filesep,sub_nm,filesep,'fMRI_analysis',filesep];
        
        [runs_idx] = MS2_task_runs_extraction(task_nm, sub_nm);
        
        [trialN, trialN_idx_aRuns, time_rest,...
            totalGain_prev, sumPerf_prev,...
            inc, absInc, R_type,...
            inc_bis, absInc_bis,...
            oneCent, cents, bills,...
            Fmax,...
            perf,...
            trialVal,...
            run_idx] = deal( NaN(1, nTrials_per_run*n_GS_runs) );
        run_Fmax = NaN(1,n_GS_runs);
        
        % loop through runs
        for iRun = 1:n_GS_runs
            jRun = runs_idx(iRun);
            run_nm = num2str(jRun);
            
            % load data
            loadStruct = getfield( load([onsets_folder,'onsets_sub',subid,'_',task_nm,'_run',run_nm,'.mat'],task_nm),task_nm);
            
            % extract relevant data
            trialN_tmp  = loadStruct.mod.all.trialN + nTrials_per_run*(iRun - 1);
            run_idx(trialN_tmp) = iRun;
            trialN(trialN_tmp)              = loadStruct.mod.all.trialN;
            trialN_idx_aRuns(trialN_tmp)    = trialN_tmp;
            switch prm_incentive_type
                case 'value'
                    inc(trialN_tmp)         = loadStruct.mod.all.incentive;
                    absInc(trialN_tmp)      = loadStruct.mod.all.absIncentive;
                case 'rank'
                    inc(trialN_tmp)         = loadStruct.mod.all.incentiveRank;
                    absInc(trialN_tmp)      = loadStruct.mod.all.absIncentiveRank;
            end
            switch prm_incentive_type_bis
                case 'value'
                    inc_bis(trialN_tmp)     = loadStruct.mod.all.incentive;
                    absInc_bis(trialN_tmp)  = loadStruct.mod.all.absIncentive;
                case 'rank'
                    inc_bis(trialN_tmp)     = loadStruct.mod.all.incentiveRank;
                    absInc_bis(trialN_tmp)  = loadStruct.mod.all.absIncentiveRank;
            end
            oneCent(trialN_tmp) = (loadStruct.mod.all.absIncentive == 0.01);
            cents(trialN_tmp) = (loadStruct.mod.all.absIncentive < 1);
            bills(trialN_tmp) = (loadStruct.mod.all.absIncentive >= 5);
            R_type(trialN_tmp)              = loadStruct.mod.all.inc_type;
            trialVal(trialN_tmp)            = loadStruct.mod.all.trialValence;
            trialVal(trialVal == -1)        = 0; % transform into binary variable 0/1
            totalGain_prev(trialN_tmp)      = loadStruct.mod.all.totalGain_prev;
            sumPerf_prev(trialN_tmp) = loadStruct.mod.all.sumPerf_prev;
            perf(trialN_tmp)                = loadStruct.mod.all.perf/100; % normalize between 0 and 1 (be careful otherwise you will have problems in the estimation for STROOP)
            time_rest(trialN_tmp)           = loadStruct.duration.all.ITI + loadStruct.duration.all.incentive; % between [1.5; 4.95]seconds
            
            % Fmax
            switch task_nm
                case 'grip'
                    Fmax_tmp = loadStruct.mod.all.Fmax;
                    Fmax(trialN_tmp) = Fmax_tmp;
                    run_Fmax(iRun) = max(Fmax_tmp,[],'omitnan');
                case 'stroop'
                    Fmax(trialN_tmp) = 1;
                    run_Fmax(iRun) = max(loadStruct.mod.all.n_pairs_solved,[],'omitnan')/n_Stroop_maxPairs; % max perf for this run
            end
        end % run loop
        
        % add incentive per condition
        [absIncGain, absIncLoss] = deal( absInc );
        absIncGain(trialVal == 0) = 0;
        absIncLoss(trialVal == 1) = 0;
        [absIncGain_bis, absIncLoss_bis] = deal( absInc_bis );
        absIncGain_bis(trialVal == 0) = 0;
        absIncLoss_bis(trialVal == 1) = 0;

        %% normalize vars
        if prm_norm == 1
            trialN  = trialN./nTrials_per_run;
            switch prm_incentive_type
                case 'value'
                    inc_norm = 20; % max value 20€
                case 'rank'
                    inc_norm = 6; % max rank = 6
            end
            inc     = inc./inc_norm;% linear value -6 => +6
            absInc  = absInc./inc_norm;% motiv value +1 => +6
            absIncGain = absIncGain./inc_norm;
            absIncLoss = absIncLoss./inc_norm;
            switch prm_incentive_type_bis
                case 'value'
                    inc_norm_bis = 20; % max value 20€
                case 'rank'
                    inc_norm_bis = 6; % max rank = 6
            end
            if ~isempty(prm_incentive_type_bis) && ~isempty(prm_incentive_var_bis)
                inc_bis     = inc_bis./inc_norm_bis;% linear value -6 => +6
                absInc_bis  = absInc_bis./inc_norm_bis;% motiv value +1 => +6
                absIncGain_bis = absIncGain_bis./inc_norm_bis;
                absIncLoss_bis = absIncLoss_bis./inc_norm_bis;
            end
            time_rest = time_rest./4.95; % if you want to normalize it,
            %         maximum possible = 4.95 seconds of rest (0.5s cross + 3.95s
            %         incentive)
            maxTotalGain = 267.10; % maximal amount possible to win if maximal perf in all trials for both runs
            totalGain_prev = totalGain_prev./maxTotalGain; % normalize by maximum possible amount to cumulate across runs
            sumPerf_prev = sumPerf_prev./nTrials_per_run;
        end
        
        %% fit the model
        % optimizing force levels
        y = perf;
        
        % inputs
        % incentive effect
        switch prm_incentive_var
            case {'inc'}
                inc_var = inc;
            case {'absInc'}
                inc_var = absInc;
            case 'absInc_plus_nomInc'
                inc_var = [absInc; inc];
            case {'absInc_plus_cond'}
                inc_var = [absInc; trialVal];
            case {'inc_perCond'}
                inc_var = [absIncGain; absIncLoss];
            otherwise
                error('inc_var name to fix');
        end
        % incentive bis effect
        switch prm_incentive_var_bis
            case ''
                inc_var_bis = [];
            case {'inc'}
                inc_var_bis = inc;
            case {'absInc'}
                inc_var_bis = absInc;
            case 'absInc_plus_nomInc'
                inc_var_bis = [absInc; inc];
            case {'absInc_plus_cond'}
                inc_var_bis = [absInc; trialVal];
            case {'inc_perCond'}
                inc_var_bis = [absIncGain; absIncLoss];
            otherwise
                error('inc_var_bis name to fix');
        end
        % fatigue effect
        switch prm_fatigue_var
            case 'trialN'
                fatigue_var = trialN;
            case 'totalGainPrev'
                fatigue_var = totalGain_prev;
            case 'sumPerfPrev'
                fatigue_var = sumPerf_prev;
            case ''
                fatigue_var = [];
            otherwise
                error('fatigue_var name to fix');
        end
        % time rest included or not
        u = [inc_var; inc_var_bis; fatigue_var; time_rest; run_idx; Fmax; trialVal];
        
        % store trial index of the trials included in the model
        trialN_kept = trialN_idx_aRuns;
        
        %% exlude NaNs (= bad trials excluded from onsets)
        bad_trials = isnan(trialN);
        y(bad_trials) = [];
        u(:,bad_trials) = [];
        trialN_kept(bad_trials) = []; % to keep track of the trials included
        trialN(bad_trials ) = []; % remove also from trialN for options.isYout to work properly
        
        %% options
        options = struct;
%         options.isYout = isnan(y);
        options.DisplayWin      = 0; % Do you want to display fig during inversion?
        options.verbose         = 0; % Do you want to display text during inversion ?
        options.inG             = prm_in_mdl;
        options.inG.task_nm     = task_nm; % store task name for case where difference between grip and stroop (expected perf mainly)
        options.inG.run_Fmax    = run_Fmax;
        
        % for case where Pmax based on Fmax(t-1) ignore the first trial for
        % the estimation of the parameters
        options.isYout = zeros(1,length(y)); % reset for each subject and each task
        if strcmp(task_nm,'grip') &&...
                ismember(prm_B_perf_C_perf_or_force,...
                {'perf_f_X','perf_f_X_bis','perf_f_X_4','perf_f_X_6'}) % ignore first trial for fitting (we don't know what prediction subjects make)
            options.isYout(trialN == 1) = 1;
        end
        
        % dimensions
        n_prm_to_fit = 1; % at least kCost
        theta_id_nm = {'kC'};
        switch prm_incentive_var
            case {'inc','absInc'}
                n_prm_to_fit = n_prm_to_fit + 1; % adding kI
                theta_id_nm{n_prm_to_fit} = 'kI';
            case {'absInc_plus_nomInc','absInc_plus_cond',...
                'inc_perCond'}
            n_prm_to_fit = n_prm_to_fit + 2; % adding kI1 AND kI2
            theta_id_nm{n_prm_to_fit-1} = 'kI1';
            theta_id_nm{n_prm_to_fit} = 'kI2';
        end
        % +demotivation effect on the benefit term
        if prm_B_time_on_benef > 0
            n_prm_to_fit = n_prm_to_fit + 1;
            theta_id_nm{n_prm_to_fit} = 'B_time';
        end
        % more incentive parameters
        switch prm_incentive_var_bis
            case {'absInc','inc'}
                n_prm_to_fit = n_prm_to_fit + 1; % kI_bis
                theta_id_nm{n_prm_to_fit} = 'kI_bis';
            case {'absInc_plus_nomInc','absInc_plus_cond',...
                'inc_perCond'}
            n_prm_to_fit = n_prm_to_fit + 2; % kI_bis1 AND kI_bis2
            theta_id_nm{n_prm_to_fit-1} = 'kI_bis1';
            theta_id_nm{n_prm_to_fit} = 'kI_bis2';
        end
        % +fatigue effect on the cost term
        if prm_C_fatigue > 0
            n_prm_to_fit = n_prm_to_fit + 1;
            theta_id_nm{n_prm_to_fit} = 'C_fatigue';
        end
        % +resting effect On the cost term
        if prm_C_rest > 0
            n_prm_to_fit = n_prm_to_fit + 1;
            theta_id_nm{n_prm_to_fit} = 'C_rest';
        end
        % kmax to estimate
        if strcmp(prm_kmax_fixed_or_free,'free')
            n_prm_to_fit = n_prm_to_fit + 2; % 1/run
            theta_id_nm{n_prm_to_fit-1} = 'kmax_r1';
            theta_id_nm{n_prm_to_fit} = 'kmax_r2';
        end
        % Fmax to estimate
        if ismember(prm_Fmax_fixed_or_free,{'free','free_bis'})
            n_prm_to_fit = n_prm_to_fit + 2; % 1/run
            theta_id_nm{n_prm_to_fit-1} = 'Fmax_r1';
            theta_id_nm{n_prm_to_fit} = 'Fmax_r2';
        end
        % perf=f(x) additional parameter to fit Perf=PM*(1+exp(-k*X))
        if ismember(prm_B_perf_C_perf_or_force,...
                {'perf_f_X','perf_f_X_bis','perf_f_X_ter',...
                'perf_f_X_4','perf_f_X_5','perf_f_X_6','perf_f_X_7'})
            n_prm_to_fit = n_prm_to_fit + 1; % 1/run
            theta_id_nm{n_prm_to_fit} = 'kX';
        end
        
        dim = struct('n',0,... % number of hidden states
            'n_theta',0 ,...    % number of evolution parameters
            'n_phi', n_prm_to_fit); % number of observation parameters
        
        % Define priors for observation parameters
        init_s_prior = 1; % one prior for all to make it simpler
        options.priors.muPhi = zeros(n_prm_to_fit, 1);
        options.priors.sigmaPhi = eye(n_prm_to_fit).*init_s_prior;
        if strcmp(prm_kmax_fixed_or_free,'free')
            options.priors.muPhi(ismember(theta_id_nm,{'kmax_r1','kmax_r2'})) = -10; % remember prior is on P in kmax =1+exp(P)
        end
        if ismember(prm_Fmax_fixed_or_free,{'free','free_bis'})
            Fmax_r1_idx = strcmp(theta_id_nm,'Fmax_r1');
            Fmax_r2_idx = strcmp(theta_id_nm,'Fmax_r2');
            switch prm_Fmax_fixed_or_free
                case 'free'
                    options.priors.muPhi(Fmax_r1_idx) = run_Fmax(1);
                    options.priors.muPhi(Fmax_r2_idx) = run_Fmax(2);
                    options.priors.sigmaPhi(Fmax_r1_idx,Fmax_r1_idx) = (10./100).*run_Fmax(1);
                    options.priors.sigmaPhi(Fmax_r2_idx,Fmax_r2_idx) = (10./100).*run_Fmax(2);
                case 'free_bis'
                    options.priors.muPhi(ismember(theta_id_nm,{'Fmax_r1','Fmax_r2'})) = -1; % remember prior is on P in Fmax =1+sigmoid(P)
            end
        end
        % multisession
        if strcmp(prm_multisess, 'multisess_singlePrmSet')
            options.multisession.split = [sum(run_idx == 1), sum(run_idx == 2)]; % info about which trials are in each session
            options.multisession.fixed.theta = []; % no theta (no f.function)
            options.multisession.fixed.phi = 1:n_prm_to_fit;
        end
        
        %% define f and g functions (evolution & observation)
        f_fname = [];
        g_fname = @MS2_GS_g_observation_perf_optimal;
        
        %% inversion
        [posterior,out] = VBA_NLStateSpaceModel(y, u, f_fname, g_fname, dim, options);
        
        %% retransform fitted paremeters
        iPrm = 1;
        kCost.(task_nm)(iS)       = fn_for_prior(posterior.muPhi(iPrm), prm_priors);
        switch prm_incentive_var
            case {'inc','absInc'}
                iPrm = iPrm + 1;
                kI.(task_nm)(iS)          = fn_for_prior(posterior.muPhi(iPrm), prm_priors);
            case {'absInc_plus_nomInc','absInc_plus_cond',...
                    'inc_perCond'}
                iPrm = iPrm + 1;
                kI1.(task_nm)(iS)          = fn_for_prior(posterior.muPhi(iPrm), prm_priors);
                iPrm = iPrm + 1;
                kI2.(task_nm)(iS)          = fn_for_prior(posterior.muPhi(iPrm), prm_priors);
        end
        % demotivation effect on the benefit term included
        if prm_B_time_on_benef > 0
            iPrm = iPrm + 1;
            kTreward.(task_nm)(iS)    = fn_for_prior(posterior.muPhi(iPrm), prm_priors);
        end

        % fatigue effect on the cost term included
        if prm_C_fatigue > 0
            iPrm = iPrm + 1;
            if prm_C_fatigue ~= 3
                kTcost.(task_nm)(iS)      = fn_for_prior(posterior.muPhi(iPrm), prm_priors);
            elseif prm_C_fatigue == 3
                kTcost.(task_nm)(iS)      = fn_for_prior(posterior.muPhi(iPrm), 'pos');
            end
        end
        % resting effect on the cost term included
        if prm_C_rest > 0
            iPrm = iPrm + 1;
            kRest.(task_nm)(iS)       = fn_for_prior(posterior.muPhi(iPrm), prm_priors);
        end
        % kmax
        if strcmp(prm_kmax_fixed_or_free,'free') % only variable where the constraint is always positive
            % run 1
            iPrm = iPrm + 1;
            kmax_estim.(task_nm)(1,iS)       = 1 + exp(posterior.muPhi(iPrm));
            % run 2
            iPrm = iPrm + 1;
            kmax_estim.(task_nm)(2,iS)       = 1 + exp(posterior.muPhi(iPrm));
        end
        % Fmax
        if ismember(prm_Fmax_fixed_or_free,{'fixed','free','free_bis'})
            switch prm_Fmax_fixed_or_free
                case 'fixed'
                    Fmax_estim.(task_nm)(1,iS)       = run_Fmax(1);
                    Fmax_estim.(task_nm)(2,iS)       = run_Fmax(2);
                case 'free'
                    % run 1
                    iPrm = iPrm + 1;
                    Fmax_estim.(task_nm)(1,iS)       = posterior.muPhi(iPrm);
                    % run 2
                    iPrm = iPrm + 1;
                    Fmax_estim.(task_nm)(2,iS)       = posterior.muPhi(iPrm);
                case 'free_bis'
                    % run 1
                    iPrm = iPrm + 1;
                    kFmax_estim.(task_nm)(1,iS)       = 1 + sigmo(posterior.muPhi(iPrm));
                    Fmax_estim.(task_nm)(1,iS)       = kFmax_estim.(task_nm)(1,iS).*run_Fmax(1);
                    % run 2
                    iPrm = iPrm + 1;
                    kFmax_estim.(task_nm)(2,iS)       = 1 + sigmo(posterior.muPhi(iPrm));
                    Fmax_estim.(task_nm)(2,iS)       = kFmax_estim.(task_nm)(2,iS).*run_Fmax(2);
            end
        end
        if ismember(prm_B_perf_C_perf_or_force,...
                {'perf_f_X','perf_f_X_bis','perf_f_X_ter',...
                'perf_f_X_4','perf_f_X_5','perf_f_X_6','perf_f_X_7'})
            iPrm = iPrm + 1;
            kX.(task_nm)(iS) = exp(posterior.muPhi(iPrm));
        end
        
        % predicted level of performance
        perf_predicted_tmp = out.suffStat.gx; % match with trials used in the modelling procedure
        perf_predicted(iS, trialN_kept) = perf_predicted_tmp;% match with trials used in the modelling procedure
        
        % if ressource fitted in the model extract variable
        if ismember(prm_B_perf_C_perf_or_force,...
                {'perf_f_X','perf_f_X_bis','perf_f_X_ter',...
                'perf_f_X_4','perf_f_X_5','perf_f_X_6','perf_f_X_7'})
            for iRun = 1:n_GS_runs
                run_trials_idx = ( trialN_kept >= 1 + nTrials_per_run*(iRun - 1) ).*( trialN_kept <= nTrials_per_run*iRun ) == 1;
                run_trials_tmp = trialN_kept(run_trials_idx);
                switch iRun
                    case 1
                        run_trials_pred_tmp = 1:length(run_trials_tmp);
                    case 2
                        first_trial_run2_idx = find(trialN_kept <= nTrials_per_run, 1, 'last');
                        run_trials_pred_tmp = (1:length(run_trials_tmp)) + first_trial_run2_idx;
                end
                % extract Perf max per run
                if ~ismember(prm_B_perf_C_perf_or_force,{'perf_f_X_ter','perf_f_X_5','perf_f_X_7'})
                    
                    switch task_nm
                        case 'grip'
                            switch prm_B_perf_C_perf_or_force
                                case {'perf_f_X','perf_f_X_bis','perf_f_X_6'}
                                    Pmax_estim_run = Fmax_estim.(task_nm)(iRun,iS)./(Fmax(run_trials_tmp).*grip_scalingFactor);
                                case 'perf_f_X_4'
                                    Pmax_estim_run = Fmax_estim.(task_nm)(iRun,iS)./Fmax(run_trials_tmp);
                            end
                        case 'stroop'
                            Pmax_estim_run = Fmax_estim.(task_nm)(iRun,iS); % no need to divide by 10 because already in percentage
                    end
                    
                elseif ismember(prm_B_perf_C_perf_or_force,{'perf_f_X_ter','perf_f_X_7'})
                    Pmax_estim_run = 1;
                    
                elseif strcmp(prm_B_perf_C_perf_or_force,'perf_f_X_5')
                    switch task_nm
                        case 'grip'
                            Pmax_estim_run = 0.75;
                        case 'stroop'
                            Pmax_estim_run = 1;
                    end
                end
                
                % Perf can not be higher than 1 => Pmax cannot go upper
                % than 1
                Pmax_estim_run(Pmax_estim_run > 1) = 1;
                % compute X predicted
                if ismember(prm_B_perf_C_perf_or_force,...
                        {'perf_f_X','perf_f_X_bis','perf_f_X_ter',...
                        'perf_f_X_4','perf_f_X_5'})
                    E_predicted(iS, run_trials_tmp) = (-1./(kX.(task_nm)(iS))).*log(1 - perf_predicted_tmp(run_trials_pred_tmp)./Pmax_estim_run);
                elseif ismember(prm_B_perf_C_perf_or_force,{'perf_f_X_6','perf_f_X_7'})
                    E_predicted(iS, run_trials_tmp) = -kX.(task_nm)(iS)./(1 - Pmax_estim_run./perf_predicted_tmp(run_trials_pred_tmp));
                else
                    error('not ready yet');
                end
            end % run loop
        end
        
        prm_qlty.(task_nm).LL(iS)   = out.fit.LL;
        prm_qlty.(task_nm).F(iS)    = out.F;
        prm_qlty.(task_nm).R2(iS)   = out.fit.R2;
    end % subject loop
    
    %% store var for each sub for each task
    prm.(task_nm).kCost           = kCost.(task_nm);
    prm.(task_nm).kTcost          = kTcost.(task_nm);
    prm.(task_nm).kRest           = kRest.(task_nm);
    switch prm_incentive_var
        case {'inc','absInc'}
                prm.(task_nm).kI              = kI.(task_nm);
        case {'absInc_plus_nomInc','absInc_plus_cond',...
                'inc_perCond'}
            prm.(task_nm).kI1             = kI1.(task_nm);
            prm.(task_nm).kI2             = kI2.(task_nm);
    end
    prm.(task_nm).kTreward        = kTreward.(task_nm);
    switch prm_incentive_var_bis
        case {'inc','absInc'}
            prm.(task_nm).kI_bis           = kI_bis.(task_nm);
        case {'absInc_plus_nomInc','absInc_plus_cond',...
                'inc_perCond'}
            prm.(task_nm).kI_bis1             = kI_bis1.(task_nm);
            prm.(task_nm).kI_bis2             = kI_bis2.(task_nm);
    end
    if strcmp(prm_kmax_fixed_or_free,'free')
        prm.(task_nm).kmax_estim      = kmax_estim.(task_nm);
    end
    if ismember(prm_Fmax_fixed_or_free,{'fixed','free','free_bis'})
        prm.(task_nm).Fmax_estim      = Fmax_estim.(task_nm);
        if strcmp(prm_Fmax_fixed_or_free,'free_bis')
            prm.(task_nm).kFmax_estim      = kFmax_estim.(task_nm);
        end
    end
    if ismember(prm_B_perf_C_perf_or_force,...
            {'perf_f_X','perf_f_X_bis','perf_f_X_ter',...
            'perf_f_X_4','perf_f_X_5','perf_f_X_6','perf_f_X_7'})
        prm.(task_nm).kX      = kX.(task_nm);
        prm.(task_nm).X_predicted  = E_predicted;
    end
    
    prm.(task_nm).F_predicted     = [];
    prm.(task_nm).perf_predicted  = perf_predicted;
    
    
    prm.(task_nm).prm_qlty.LL     = prm_qlty.(task_nm).LL;
    prm.(task_nm).prm_qlty.F      = prm_qlty.(task_nm).F;
    
    prm.(task_nm).prm_qlty.R2     = prm_qlty.(task_nm).R2;
    %% average across subs
    % Note: if only one subject fitted, of course this measure will be
    % redundant
    m_prm.(task_nm).kCost           = mean(kCost.(task_nm), 'omitnan');
    m_prm.(task_nm).kTcost          = mean(kTcost.(task_nm), 'omitnan');
    m_prm.(task_nm).kRest           = mean(kRest.(task_nm), 'omitnan');
    switch prm_incentive_var
        case {'inc','absInc'}
            m_prm.(task_nm).kI              = mean(kI.(task_nm), 'omitnan');
            
        case {'absInc_plus_nomInc','absInc_plus_cond',...
                'inc_perCond'}
            m_prm.(task_nm).kI1              = mean(kI1.(task_nm), 'omitnan');
            m_prm.(task_nm).kI2              = mean(kI2.(task_nm), 'omitnan');
    end
    m_prm.(task_nm).kTreward        = mean(kTreward.(task_nm), 'omitnan');
    switch prm_incentive_var_bis
        case {'inc','absInc'}
            m_prm.(task_nm).kI_bis      = mean(kI_bis.(task_nm), 'omitnan');
        case {'absInc_plus_nomInc','absInc_plus_cond',...
                'inc_perCond'}
            m_prm.(task_nm).kI_bis1     = mean(kI_bis1.(task_nm), 'omitnan');
            m_prm.(task_nm).kI_bis2     = mean(kI_bis2.(task_nm), 'omitnan');
    end
    if strcmp(prm_kmax_fixed_or_free,'free')
        m_prm.(task_nm).kmax_estim      = mean(kmax_estim.(task_nm),2, 'omitnan');
    end
    if ismember(prm_Fmax_fixed_or_free,{'fixed','free','free_bis'})
        m_prm.(task_nm).Fmax_estim      = mean(Fmax_estim.(task_nm),2, 'omitnan');
        if strcmp(prm_Fmax_fixed_or_free,'free_bis')
            m_prm.(task_nm).kFmax_estim      = mean(kFmax_estim.(task_nm),2, 'omitnan');
        end
    end
    if ismember(prm_B_perf_C_perf_or_force,...
            {'perf_f_X','perf_f_X_bis','perf_f_X_ter',...
            'perf_f_X_4','perf_f_X_5','perf_f_X_6','perf_f_X_7'})
        m_prm.(task_nm).kX      = mean(kX.(task_nm), 'omitnan');
    end
    
    m_prm.(task_nm).prm_qlty.LL     = mean(prm_qlty.(task_nm).LL, 'omitnan');
    m_prm.(task_nm).prm_qlty.F      = mean(prm_qlty.(task_nm).F, 'omitnan');
    m_prm.(task_nm).prm_qlty.R2     = mean(prm_qlty.(task_nm).R2, 'omitnan');
end % task


end % function