function[matlabbatch] = First_level_MS2_gripRP_prm(matlabbatch, GLMprm,...
    subid, sub_idx, iRun, subj_analysis_folder)
% [matlabbatch] = First_level_MS2_gripRP_prm(matlabbatch, GLMprm,...
%     subid, sub_idx, iRun, subj_analysis_folder)
% First_level_MS2_gripRP_prm laods the first level infos for
% First_level_MS2_megaconcatenation_NicoC_batch.m
%
% INPUTS
% matlabbatch: batch that needs to be filled with the info of the grip task
%
% GLMprm: GLM parameters extracted with which_GLM_MS2 to know which
% conditions (onsets/modulators) are to be used
%
% subid: string with subject number
%
% sub_idx: batch subject number for matlabbatch
%
% iRun: run to extract number
%
% subj_analysis_folder: subject analysis folder path
%
% OUTPUTS
% matlabbatch: matlabbatch for first level fMRI analysis filled with the
% infos of the current task
%
% See also First_level_MS2_megaconcatenation_NicoC_batch.m

path_split = strsplit(subj_analysis_folder,filesep);
sub_nm = path_split{end-2};
run_nm = num2str(iRun);
gripRPprm = GLMprm.grip;
gal = GLMprm.gal;
task_id = 'G_';
n_trials_perRun = 60;

%% load onsets and modulators for the GLM
grip   = getfield( load([subj_analysis_folder, filesep, 'onsets_sub',subid,'_grip_run',run_nm],'grip'), 'grip');
onset       = grip.onset; % onsets
duration    = grip.duration;
mod         = grip.mod; % modulators

%% see how many onsets need to be considered
n_cond = 0;
cond_nm = {};
modName = struct;
modulators = struct;

%% onset incentive
o_inc   = gripRPprm.o_inc;
dur_inc = gripRPprm.dur_inc;
mod_inc = gripRPprm.mod_inc;
switch o_inc
    case {1,2}
        if o_inc == 1 % all trials pooled
            curr_oInc_names = {'G_inc'};
            trialTypes_oInc = {'all'};
        elseif o_inc == 2 % to gain/to lose trials separated
            curr_oInc_names = {'G_inc_toGain','G_inc_toLose'};
            trialTypes_oInc = {'toGain','toLose'};
        end
        n_cond = n_cond + length(curr_oInc_names);
        cond_nm = [cond_nm, curr_oInc_names];
        
        for iCond_oInc = 1:length(curr_oInc_names)
            curr_oInc_name = curr_oInc_names{iCond_oInc};
            trialType_oInc = trialTypes_oInc{iCond_oInc};
            
            onsetCond.(curr_oInc_name) = onset.(trialType_oInc).incentive;
            
            % duration
            switch dur_inc
                case 0 % stick
                    durCond.(curr_oInc_name) = 0;
                case 1 % boxcar from incentive to effort scale display
                    durCond.(curr_oInc_name) = duration.(trialType_oInc).incentive;
                case 2 % boxcar from incentive to feedback
                    durCond.(curr_oInc_name) = duration.(trialType_oInc).incentive +...
                        duration.(trialType_oInc).effortScale;
            end
            
            % modulators
            Nmod.(curr_oInc_name) = 0;
            
            %% G/L condition
            if mod_inc.GL_cond ~= 0
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_GLcond = Nmod.(curr_oInc_name);
                switch mod_inc.GL_cond
                    case 1 % 0/1
                        modName.(curr_oInc_name)(idx_oInc_GLcond) = {'GL_cond'};
                        GL_cond = mod.(trialType_oInc).trialValence;
                        GL_cond(GL_cond == -1) = 0; % convert into 0/1 variable
                        modulators.(curr_oInc_name)(:,idx_oInc_GLcond) = GL_cond;
                    case 2 % -1/1
                        modName.(curr_oInc_name)(idx_oInc_GLcond) = {'GL_cond'};
                        modulators.(curr_oInc_name)(:,idx_oInc_GLcond) = mod.(trialType_oInc).trialValence;
                    otherwise
                        error('case not ready yet');
                end
            end
            
            %% luminance
            if mod_inc.lum ~= 0 && ismember(mod_inc.lum, [2,4,6,8])
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_lum = Nmod.(curr_oInc_name);
                switch mod_inc.lum
                    case 2
                        modName.(curr_oInc_name)(idx_oInc_lum) = {'luminance'};
                        modulators.(curr_oInc_name)(:,idx_oInc_lum) = mod.(trialType_oInc).lum_inc;
                    case 4
                        modName.(curr_oInc_name)(idx_oInc_lum) = {'luminance'};
                        modulators.(curr_oInc_name)(:,idx_oInc_lum) = nanzscore(mod.(trialType_oInc).lum_inc);
                    otherwise
                        error('case not ready yet');
                end
            end
            
            %% trial number
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[2,4])
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_trialN = Nmod.(curr_oInc_name);
                switch mod_inc.trialN
                    case 2
                        modName.(curr_oInc_name)(idx_oInc_trialN) = {'trialN'};
                        modulators.(curr_oInc_name)(:,idx_oInc_trialN) = mod.(trialType_oInc).trialN;
                    case 4
                        modName.(curr_oInc_name)(idx_oInc_trialN) = {'trialN'};
                        modulators.(curr_oInc_name)(:,idx_oInc_trialN) = nanzscore(mod.(trialType_oInc).trialN);
                end
            end

            %% RT
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[1])
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_RT_fp = Nmod.(curr_oInc_name);
                switch mod_inc.RT_fp
                    case 2
                        modName.(curr_oInc_name)(idx_oInc_RT_fp) = {'RT_fp'};
                        modulators.(curr_oInc_name)(:,idx_oInc_RT_fp) = mod.(trialType_oInc).RT_fp;
                    otherwise
                        error('not ready yet');
                end
            end

            %% model benefit/cost terms
            if mod_inc.mdl_n ~= 0
                loadStruct_mdl   = load([subj_analysis_folder, filesep, 'GS_model_sub',subid,'_grip_run',run_nm,'.mat']);
                
                %% predicted ressource level
                if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[5,6])
                    Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                    idx_oInc_Epred = Nmod.(curr_oInc_name);
                    modName.(curr_oInc_name)(idx_oInc_Epred) = {'E_pred'};
                    switch mod_inc.E_pred
                        case 5
                            E = loadStruct_mdl.E.(['model_',num2str(mod_inc.mdl_n)]);
                            modulators.(curr_oInc_name)(:,idx_oInc_Epred) = E;
                        case 6
                            E_pred = loadStruct_mdl.E_pred.(['model_',num2str(mod_inc.mdl_n)]);
                            modulators.(curr_oInc_name)(:,idx_oInc_Epred) = E_pred;
                    end
                end
            end
            
            %% reward type: piece of money or bill
            if mod_inc.R_type ~= 0
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_Rtype = Nmod.(curr_oInc_name);
                switch mod_inc.R_type
                    case 1
                        modName.(curr_oInc_name)(idx_oInc_Rtype) = {'reward piece of money or bill'};
                        modulators.(curr_oInc_name)(:,idx_oInc_Rtype) = mod.(trialType_oInc).inc_type;
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% incentive
            if mod_inc.inc ~= 0
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_inc = Nmod.(curr_oInc_name);
                switch mod_inc.inc
                    case 1
                        modName.(curr_oInc_name)(idx_oInc_inc) = {'incentive rank nominal value'};
                        modulators.(curr_oInc_name)(:,idx_oInc_inc) = mod.(trialType_oInc).incentiveRank;
                    case 2
                        modName.(curr_oInc_name)(idx_oInc_inc) = {'incentive nominal value'};
                        modulators.(curr_oInc_name)(:,idx_oInc_inc) = mod.(trialType_oInc).incentive;
                    case 3
                        modName.(curr_oInc_name)(idx_oInc_inc) = {'incentive rank motivational value'};
                        modulators.(curr_oInc_name)(:,idx_oInc_inc) = mod.(trialType_oInc).absIncentiveRank;
                    case 4
                        modName.(curr_oInc_name)(idx_oInc_inc) = {'incentive motivational value'};
                        modulators.(curr_oInc_name)(:,idx_oInc_inc) = mod.(trialType_oInc).absIncentive;
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% incentive bis
            if mod_inc.inc_bis ~= 0
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_inc_bis = Nmod.(curr_oInc_name);
                switch mod_inc.inc_bis
                    case 1
                        modName.(curr_oInc_name)(idx_oInc_inc_bis) = {'incentive rank nominal value'};
                        modulators.(curr_oInc_name)(:,idx_oInc_inc_bis) = mod.(trialType_oInc).incentiveRank;
                    case 2
                        modName.(curr_oInc_name)(idx_oInc_inc_bis) = {'incentive nominal value'};
                        modulators.(curr_oInc_name)(:,idx_oInc_inc_bis) = mod.(trialType_oInc).incentive;
                    case 3
                        modName.(curr_oInc_name)(idx_oInc_inc_bis) = {'incentive rank motivationnal value'};
                        modulators.(curr_oInc_name)(:,idx_oInc_inc_bis) = mod.(trialType_oInc).absIncentiveRank;
                    case 4
                        modName.(curr_oInc_name)(idx_oInc_inc_bis) = {'incentive motivational value'};
                        modulators.(curr_oInc_name)(:,idx_oInc_inc_bis) = mod.(trialType_oInc).absIncentive;
                    otherwise
                        error('not ready yet');
                end
            end

            %% confidence
            if mod_inc.conf ~= 0
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_conf = Nmod.(curr_oInc_name);
                switch mod_inc.conf
                    case 1
                        modName.(curr_oInc_name)(idx_oInc_conf) = {'confidence (rank)'};
                        modulators.(curr_oInc_name)(:,idx_oInc_conf) = (mod.(trialType_oInc).absIncentiveRank - 3.5).^2;
                    case 2
                        modName.(curr_oInc_name)(idx_oInc_conf) = {'confidence (money)'};
                        modulators.(curr_oInc_name)(:,idx_oInc_conf) = (mod.(trialType_oInc).absIncentive - 4.4517).^2;
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% ROI activity
            if mod_inc.ROI_activity_yn ~= 0
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_ROI = Nmod.(curr_oInc_name);
                modName.(curr_oInc_name)(idx_oInc_ROI) = {[]};
                % get the data
                trialN_idx = mod.(trialType_oInc).trialN; % trial index
                ROI_nm = mod_inc.ROI_activity_ROI_nm; % ROI name
                GLM_ROI_nm = mod_inc.ROI_activity_GLM;
                ROI_period_nm = mod_inc.ROI_activity_period; % trial period when ROI extracted
                jRun_nb = (iRun <= 3)*1 + (iRun > 3)*2;
                jRun_nm = ['run',num2str(jRun_nb)];
                ROI_data = getfield(load([subj_analysis_folder,'ROI',filesep,'GLM',GLM_ROI_nm,'_',ROI_nm,'.mat'],...
                    'sub_ROI_trial_b_trial'),'sub_ROI_trial_b_trial');
                modulators.(curr_oInc_name)(:,idx_oInc_ROI) = ROI_data.(ROI_nm).G.(ROI_period_nm).(jRun_nm).(sub_nm)(trialN_idx);
            end
            
            %% effort
            if mod_inc.Effort ~= 0
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_Eff = Nmod.(curr_oInc_name);
                switch mod_inc.Effort
                    case 1
                        modName.(curr_oInc_name)(idx_oInc_Eff) = {'F_peak'};
                        modulators.(curr_oInc_name)(:,idx_oInc_Eff) = mod.(trialType_oInc).peakForce;
                    case 2
                        modName.(curr_oInc_name)(idx_oInc_Eff) = {'F_integral'};
                        modulators.(curr_oInc_name)(:,idx_oInc_Eff) = mod.(trialType_oInc).intForce;
                    case 3
                        runFmax = mod.all.Fmax(end);
                        modName.(curr_oInc_name)(idx_oInc_Eff) = {'F_peak'};
                        modulators.(curr_oInc_name)(:,idx_oInc_Eff) = mod.(trialType_oInc).peakForce./runFmax;
                    case 4
                        runFmax = mod.all.Fmax(end);
                        modName.(curr_oInc_name)(idx_oInc_Eff) = {'F_integral'};
                        modulators.(curr_oInc_name)(:,idx_oInc_Eff) = mod.(trialType_oInc).intForce./runFmax;
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% model benefit/cost terms
            if mod_inc.mdl_n ~= 0
                loadStruct_mdl   = load([subj_analysis_folder, filesep, 'GS_model_sub',subid,'_grip_run',run_nm,'.mat']);
                EV              = loadStruct_mdl.EV.(['model_',num2str(mod_inc.mdl_n)]);
                cost            = loadStruct_mdl.cost.(['model_',num2str(mod_inc.mdl_n)]);
                benefit         = loadStruct_mdl.benefit.(['model_',num2str(mod_inc.mdl_n)]);
                EV_pred         = loadStruct_mdl.EV_pred.(['model_',num2str(mod_inc.mdl_n)]);
                cost_pred       = loadStruct_mdl.cost_pred.(['model_',num2str(mod_inc.mdl_n)]);
                benefit_pred    = loadStruct_mdl.benefit_pred.(['model_',num2str(mod_inc.mdl_n)]);
%                 cost_var        = loadStruct_mdl.cost_var.(['model_',num2str(mod_inc.mdl_n)]);
%                 benefit_var     = loadStruct_mdl.benefit_var.(['model_',num2str(mod_inc.mdl_n)]);
                expected_payoff = loadStruct_mdl.expected_payoff.(['model_',num2str(mod_inc.mdl_n)]);
                
                %% predicted ressource level
                if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[3,4])
                    Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                    idx_oInc_Epred = Nmod.(curr_oInc_name);
                    modName.(curr_oInc_name)(idx_oInc_Epred) = {'E_pred'};
                    switch mod_inc.E_pred
                        case 3
                            E = loadStruct_mdl.E.(['model_',num2str(mod_inc.mdl_n)]);
                            modulators.(curr_oInc_name)(:,idx_oInc_Epred) = E;
                        case 4
                            E_pred = loadStruct_mdl.E_pred.(['model_',num2str(mod_inc.mdl_n)]);
                            modulators.(curr_oInc_name)(:,idx_oInc_Epred) = E_pred;
                    end
                end
                
                %% cost term
                if mod_inc.mdl_cost ~= 0
                    Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                    idx_oInc_modelCost = Nmod.(curr_oInc_name);
                    modName.(curr_oInc_name)(idx_oInc_modelCost) = {'model_cost'};
                    switch mod_inc.mdl_cost
                        case 1
                            modulators.(curr_oInc_name)(:,idx_oInc_modelCost) = cost;
                        case 2
                            modulators.(curr_oInc_name)(:,idx_oInc_modelCost) = cost_pred;
                        case 3
                            modulators.(curr_oInc_name)(:,idx_oInc_modelCost) = cost_var;
                    end
                end
                
                %% benefit term
                if mod_inc.mdl_benefit ~= 0
                    Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                    idx_oInc_modelBenef = Nmod.(curr_oInc_name);
                    modName.(curr_oInc_name)(idx_oInc_modelBenef) = {'model_benefit'};
                    switch mod_inc.mdl_benefit
                        case 1
                            modulators.(curr_oInc_name)(:,idx_oInc_modelBenef) = benefit;
                        case 2
                            modulators.(curr_oInc_name)(:,idx_oInc_modelBenef) = benefit_pred;
                        case 3
                            modulators.(curr_oInc_name)(:,idx_oInc_modelBenef) = benefit_var;
                        case 4
                            modulators.(curr_oInc_name)(:,idx_oInc_modelBenef) = expected_payoff;
                    end
                end
                
                %% EV term
                if mod_inc.mdl_EV ~= 0
                    Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                    idx_oInc_modelEV = Nmod.(curr_oInc_name);
                    modName.(curr_oInc_name)(idx_oInc_modelEV) = {'model_EV'};
                    switch mod_inc.mdl_EV
                        case 1
                            modulators.(curr_oInc_name)(:,idx_oInc_modelEV) = EV;
                        case 2
                            modulators.(curr_oInc_name)(:,idx_oInc_modelEV) = EV_pred;
                    end
                end
                
            end
            
            %% trial number
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[5,6])
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_trialN = Nmod.(curr_oInc_name);
                switch mod_inc.trialN
                    case 5
                        modName.(curr_oInc_name)(idx_oInc_trialN) = {'trialN'};
                        modulators.(curr_oInc_name)(:,idx_oInc_trialN) = mod.(trialType_oInc).trialN;
                    case 6
                        modName.(curr_oInc_name)(idx_oInc_trialN) = {'trialN'};
                        modulators.(curr_oInc_name)(:,idx_oInc_trialN) = nanzscore(mod.(trialType_oInc).trialN);
                end
            end
            
            %% RT
            if mod_inc.RT_fp ~= 0 && ~ismember(mod_inc.RT_fp,[2])
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_RT_fp = Nmod.(curr_oInc_name);
                switch mod_inc.RT_fp
                    case 1
                        modName.(curr_oInc_name)(idx_oInc_RT_fp) = {'RT_fp'};
                        modulators.(curr_oInc_name)(:,idx_oInc_RT_fp) = mod.(trialType_oInc).RT_fp;
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% total gain until now
            if mod_inc.totalGain_prev ~= 0
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_totalGain_prev = Nmod.(curr_oInc_name);
                switch mod_inc.totalGain_prev
                    case 1
                        modName.(curr_oInc_name)(idx_oInc_totalGain_prev) = {'total gain until now'};
                        modulators.(curr_oInc_name)(:,idx_oInc_totalGain_prev) = mod.(trialType_oInc).totalGain_prev;
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% sum(performance previous trials)
            if mod_inc.sumPerfPrev ~= 0
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_sumPerfPrev = Nmod.(curr_oInc_name);
                switch mod_inc.sumPerfPrev
                    case 1
                        modName.(curr_oInc_name)(idx_oInc_sumPerfPrev) = {'sum(peformance previous trials)'};
                        modulators.(curr_oInc_name)(:,idx_oInc_sumPerfPrev) = mod.(trialType_oInc).sumPerf_prev;
                    case 2
                        modName.(curr_oInc_name)(idx_oInc_sumPerfPrev) = {'sum(peformance previous trials)'};
                        modulators.(curr_oInc_name)(:,idx_oInc_sumPerfPrev) = mod.(trialType_oInc).sumPerf_prev./n_trials_perRun; % normalize by total number of trials
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% luminance
            if mod_inc.lum ~= 0 && ismember(mod_inc.lum, [1,3,5,7])
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_lum = Nmod.(curr_oInc_name);
                switch mod_inc.lum
                    case 1
                        modName.(curr_oInc_name)(idx_oInc_lum) = {'luminance'};
                        modulators.(curr_oInc_name)(:,idx_oInc_lum) = mod.(trialType_oInc).lum_inc;
                    case 3
                        modName.(curr_oInc_name)(idx_oInc_lum) = {'luminance'};
                        modulators.(curr_oInc_name)(:,idx_oInc_lum) = nanzscore(mod.(trialType_oInc).lum_inc);
                    otherwise
                        error('case not ready yet');
                end
            end
            
            %% trial number
            if mod_inc.trialN ~= 0 && ismember(mod_inc.trialN,[1,3])
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_trialN = Nmod.(curr_oInc_name);
                switch mod_inc.trialN
                    case 1
                        modName.(curr_oInc_name)(idx_oInc_trialN) = {'trialN'};
                        modulators.(curr_oInc_name)(:,idx_oInc_trialN) = mod.(trialType_oInc).trialN;
                    case 3
                        modName.(curr_oInc_name)(idx_oInc_trialN) = {'trialN'};
                        modulators.(curr_oInc_name)(:,idx_oInc_trialN) = nanzscore(mod.(trialType_oInc).trialN);
                end
            end
            
            %% incentive*trial number
            if mod_inc.inc_x_trialN ~= 0
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_incXtrialN = Nmod.(curr_oInc_name);
                switch mod_inc.inc_x_trialN
                    case 1
                        modName.(curr_oInc_name)(idx_oInc_incXtrialN) = {'inc_x_trialN'};
                        modulators.(curr_oInc_name)(:,idx_oInc_incXtrialN) = (mod.(trialType_oInc).absIncentiveRank).*(mod.(trialType_oInc).trialN);
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% predicted ressource level
            if mod_inc.E_pred ~= 0 && ismember(mod_inc.E_pred,[1,2])
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_Epred = Nmod.(curr_oInc_name);
                loadStruct_mdl   = load([subj_analysis_folder, filesep, 'GS_model_sub',subid,'_grip_run',run_nm,'.mat']);
                modName.(curr_oInc_name)(idx_oInc_Epred) = {'E_pred'};
                switch mod_inc.E_pred
                    case 1
                        E = loadStruct_mdl.E.(['model_',num2str(mod_inc.mdl_n)]);
                        modulators.(curr_oInc_name)(:,idx_oInc_Epred) = E;
                    case 2
                        E_pred = loadStruct_mdl.E_pred.(['model_',num2str(mod_inc.mdl_n)]);
                        modulators.(curr_oInc_name)(:,idx_oInc_Epred) = E_pred;
                end
            end
            
            %% performance
            if mod_inc.perf ~= 0
                Nmod.(curr_oInc_name) = Nmod.(curr_oInc_name) + 1;
                idx_oInc_perf = Nmod.(curr_oInc_name);
                switch mod_inc.perf
                    case 1 % 1-100
                        modName.(curr_oInc_name)(idx_oInc_perf) = {'perf'};
                        modulators.(curr_oInc_name)(:,idx_oInc_perf) = mod.(trialType_oInc).perf;
                    case 2 % 0-1
                        modName.(curr_oInc_name)(idx_oInc_perf) = {'perf'};
                        modulators.(curr_oInc_name)(:,idx_oInc_perf) = mod.(trialType_oInc).perf./100;
                    case 3
                        modName.(curr_oInc_name)(idx_oInc_perf) = {'perf'};
                        modulators.(curr_oInc_name)(:,idx_oInc_perf) = mod.(trialType_oInc).perf_orth;
                    case 4
                        loadStruct_mdl   = load([subj_analysis_folder, filesep, 'GS_model_sub',subid,'_grip_run',run_nm,'.mat']);
                        perf_pred = loadStruct_mdl.perf_pred.(['model_',num2str(mod_inc.mdl_n)]);
                        modName.(curr_oInc_name)(idx_oInc_perf) = {'perf'};
                        modulators.(curr_oInc_name)(:,idx_oInc_perf) = perf_pred;
                end
            end
            
        end % modulators
    case 0 % nothing happens because not included in the GLM
    otherwise
        error(['case o_inc = ',num2str(o_inc),' not ready']);
end % onset incentive type

%% onset effort scale
o_dispE     = gripRPprm.o_dispE;
dur_dispE   = gripRPprm.dur_dispE;
mod_dispE   = gripRPprm.mod_dispE;
switch o_dispE
    case {1,2}
        if o_dispE == 1 % all trials pooled
            curr_oEScale_names = {'G_effortScale'};
            trialTypes_oEScale = {'all'};
        elseif o_dispE == 2 % to gain/to lose trials separated
            curr_oEScale_names = {'G_effortScale_toGain','G_effortScale_toLose'};
            trialTypes_oEScale = {'toGain','toLose'};
        end
        n_cond = n_cond + length(curr_oEScale_names);
        cond_nm = [cond_nm, curr_oEScale_names];
        
        for iCond_oDispE = 1:length(curr_oEScale_names)
            curr_oEScale_name = curr_oEScale_names{iCond_oDispE};
            trialType_oEScale = trialTypes_oEScale{iCond_oDispE};
            
            onsetCond.(curr_oEScale_name) = onset.(trialType_oEScale).effortScale;
            
            %% duration
            switch dur_dispE
                case 0 % stick
                    durCond.(curr_oEScale_name) = 0;
                case 1 % boxcar from effort scale to feedback
                    durCond.(curr_oEScale_name) = duration.(trialType_oEScale).effortScale;
                case 2 % boxcar from effort scale to RT
                    durCond.(curr_oEScale_name) = mod.(trialType_oEScale).RT_fp;
            end
            
            % modulators
            Nmod.(curr_oEScale_name) = 0;
            
            %% G/L condition
            if mod_dispE.GL_cond ~= 0
                Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                idx_oDispE_GLcond = Nmod.(curr_oEScale_name);
                switch mod_dispE.GL_cond
                    case 1 % 0/1
                        modName.(curr_oEScale_name)(idx_oDispE_GLcond) = {'GL_cond'};
                        GL_cond = mod.(trialType_oEScale).trialValence;
                        GL_cond(GL_cond == -1) = 0; % convert into 0/1 variable
                        modulators.(curr_oEScale_name)(:,idx_oDispE_GLcond) = GL_cond;
                    case 2 % -1/1
                        modName.(curr_oEScale_name)(idx_oDispE_GLcond) = {'GL_cond'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_GLcond) = mod.(trialType_oEScale).trialValence;
                    otherwise
                        error('case not ready yet');
                end
            end
            
            %% luminance
            if mod_dispE.lum ~= 0 && ismember(mod_dispE.lum, [2,4,6,8])
                Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                idx_oDispE_lum = Nmod.(curr_oEScale_name);
                switch mod_dispE.lum
                    case 2
                        modName.(curr_oEScale_name)(idx_oDispE_lum) = {'luminance'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_lum) = mod.(trialType_oEScale).lum_inc;
                    case 4
                        modName.(curr_oEScale_name)(idx_oDispE_lum) = {'luminance'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_lum) = nanzscore(mod.(trialType_oEScale).lum_inc);
                    otherwise
                        error('case not ready yet');
                end
            end
            
            %% trial number
            if mod_dispE.trialN ~= 0 && ismember(mod_dispE.trialN,[2,4])
                Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                idx_oDispE_trialN = Nmod.(curr_oEScale_name);
                switch mod_dispE.trialN
                    case 2
                        modName.(curr_oEScale_name)(idx_oDispE_trialN) = {'trialN'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_trialN) = mod.(trialType_oEScale).trialN;
                    case 4
                        modName.(curr_oEScale_name)(idx_oDispE_trialN) = {'trialN'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_trialN) = nanzscore(mod.(trialType_oEScale).trialN);
                end
            end

            %% RT (first regressor)
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[1])
                Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                idx_oDispE_RT_fp = Nmod.(curr_oEScale_name);
                switch mod_dispE.RT_fp
                    case 2
                        modName.(curr_oEScale_name)(idx_oDispE_RT_fp) = {'RT_fp'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_RT_fp) = mod.(trialType_oEScale).RT_fp;
                    otherwise
                        error('not ready yet');
                end
            end

            %% model benefit/cost terms
            if mod_dispE.mdl_n ~= 0
                loadStruct_mdl   = load([subj_analysis_folder, filesep, 'GS_model_sub',subid,'_grip_run',run_nm,'.mat']);

                %% predicted ressource level
                if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[5,6])
                    Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                    idx_oDispE_Epred = Nmod.(curr_oEScale_name);
                    modName.(curr_oEScale_name)(idx_oDispE_Epred) = {'E_pred'};
                    switch mod_dispE.E_pred
                        case 5
                            E = loadStruct_mdl.E.(['model_',num2str(mod_dispE.mdl_n)]);
                            modulators.(curr_oEScale_name)(:,idx_oDispE_Epred) = E;
                        case 6
                            E_pred = loadStruct_mdl.E_pred.(['model_',num2str(mod_dispE.mdl_n)]);
                            modulators.(curr_oEScale_name)(:,idx_oDispE_Epred) = E_pred;
                    end
                end
            end
            
            %% reward type: piece of money or bill
            if mod_dispE.R_type ~= 0
                Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                idx_oDispE_Rtype = Nmod.(curr_oEScale_name);
                switch mod_dispE.R_type
                    case 1
                        modName.(curr_oEScale_name)(idx_oDispE_Rtype) = {'reward piece of money or bill'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_Rtype) = mod.(trialType_oEScale).inc_type;
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% incentive
            if mod_dispE.inc ~= 0
                Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                idx_oDispE_inc = Nmod.(curr_oEScale_name);
                switch mod_dispE.inc
                    case 1
                        modName.(curr_oEScale_name)(idx_oDispE_inc) = {'incentive rank nominal value'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_inc) = mod.(trialType_oEScale).incentiveRank;
                    case 2
                        modName.(curr_oEScale_name)(idx_oDispE_inc) = {'incentive nominal value'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_inc) = mod.(trialType_oEScale).incentive;
                    case 3
                        modName.(curr_oEScale_name)(idx_oDispE_inc) = {'incentive rank motivationnal value'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_inc) = mod.(trialType_oEScale).absIncentiveRank;
                    case 4
                        modName.(curr_oEScale_name)(idx_oDispE_inc) = {'incentive motivational value'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_inc) = mod.(trialType_oEScale).absIncentive;
                        
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% incentive bis
            if mod_dispE.inc_bis ~= 0
                Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                idx_oDispE_inc_bis = Nmod.(curr_oEScale_name);
                switch mod_dispE.inc_bis
                    case 1
                        modName.(curr_oEScale_name)(idx_oDispE_inc_bis) = {'incentive rank nominal value'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_inc_bis) = mod.(trialType_oEScale).incentiveRank;
                    case 2
                        modName.(curr_oEScale_name)(idx_oDispE_inc_bis) = {'incentive nominal value'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_inc_bis) = mod.(trialType_oEScale).incentive;
                    case 3
                        modName.(curr_oEScale_name)(idx_oDispE_inc_bis) = {'incentive rank motivationnal value'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_inc_bis) = mod.(trialType_oEScale).absIncentiveRank;
                    case 4
                        modName.(curr_oEScale_name)(idx_oDispE_inc_bis) = {'incentive motivational value'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_inc_bis) = mod.(trialType_oEScale).absIncentive;
                        
                    otherwise
                        error('not ready yet');
                end
            end

            %% confidence
            if mod_dispE.conf ~= 0
                Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                idx_oDispE_conf = Nmod.(curr_oEScale_name);
                switch mod_dispE.conf
                    case 1
                        modName.(curr_oEScale_name)(idx_oDispE_conf) = {'confidence (rank)'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_conf) = (mod.(trialType_oEScale).absIncentiveRank - 3.5).^2;
                    case 2
                        modName.(curr_oEScale_name)(idx_oDispE_conf) = {'confidence (money)'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_conf) = (mod.(trialType_oEScale).absIncentive - 4.4517).^2;
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% ROI activity
            if mod_dispE.ROI_activity_yn ~= 0
                Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                idx_oDispE_ROI = Nmod.(curr_oEScale_name);
                modName.(curr_oEScale_name)(idx_oDispE_ROI) = {[]};
                % get the data
                trialN_idx = mod.(trialType_oEScale).trialN; % trial index
                ROI_nm          = mod_dispE.ROI_activity_ROI_nm; % ROI name
                GLM_ROI_nm      = mod_dispE.ROI_activity_GLM;
                ROI_period_nm   = mod_dispE.ROI_activity_period; % trial period when ROI extracted
                jRun_nb = (iRun <= 3)*1 + (iRun > 3)*2;
                jRun_nm = ['run',num2str(jRun_nb)];
                ROI_data = getfield(load([subj_analysis_folder,'ROI',filesep,'GLM',GLM_ROI_nm,'_',ROI_nm,'.mat'],...
                    'sub_ROI_trial_b_trial'),'sub_ROI_trial_b_trial');
                modulators.(curr_oEScale_name)(:,idx_oDispE_ROI) = ROI_data.(ROI_nm).G.(ROI_period_nm).(jRun_nm).(sub_nm)(trialN_idx);
            end
            
            %% effort
            if mod_dispE.Effort ~= 0
                Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                idx_oDispE_Effort = Nmod.(curr_oEScale_name);
                switch mod_dispE.Effort
                    case 1
                        modName.(curr_oEScale_name)(idx_oDispE_Effort) = {'F_peak'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_Effort) = mod.(trialType_oEScale).peakForce;
                    case 2
                        modName.(curr_oEScale_name)(idx_oDispE_Effort) = {'F_integral'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_Effort) = mod.(trialType_oEScale).intForce;
                    case 3
                        runFmax = mod.all.Fmax(end);
                        modName.(curr_oEScale_name)(idx_oDispE_Effort) = {'F_peak'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_Effort) = mod.(trialType_oEScale).peakForce./runFmax;
                    case 4
                        runFmax = mod.all.Fmax(end);
                        modName.(curr_oEScale_name)(idx_oDispE_Effort) = {'F_integral'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_Effort) = mod.(trialType_oEScale).intForce./runFmax;
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% model benefit/cost terms
            if mod_dispE.mdl_n ~= 0
                loadStruct_mdl   = load([subj_analysis_folder, filesep, 'GS_model_sub',subid,'_grip_run',run_nm,'.mat']);
                EV = loadStruct_mdl.EV.(['model_',num2str(mod_dispE.mdl_n)]);
                cost = loadStruct_mdl.cost.(['model_',num2str(mod_dispE.mdl_n)]);
                benefit = loadStruct_mdl.benefit.(['model_',num2str(mod_dispE.mdl_n)]);
                EV_pred = loadStruct_mdl.EV_pred.(['model_',num2str(mod_dispE.mdl_n)]);
                cost_pred = loadStruct_mdl.cost_pred.(['model_',num2str(mod_dispE.mdl_n)]);
                benefit_pred = loadStruct_mdl.benefit_pred.(['model_',num2str(mod_dispE.mdl_n)]);
%                 cost_var = loadStruct_mdl.cost_var.(['model_',num2str(mod_dispE.mdl_n)]);
%                 benefit_var = loadStruct_mdl.benefit_var.(['model_',num2str(mod_dispE.mdl_n)]);
                expected_payoff = loadStruct_mdl.expected_payoff.(['model_',num2str(mod_dispE.mdl_n)]);
                
                %% predicted ressource level
                if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[3,4])
                    Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                    idx_oDispE_Epred = Nmod.(curr_oEScale_name);
                    modName.(curr_oEScale_name)(idx_oDispE_Epred) = {'E_pred'};
                    switch mod_dispE.E_pred
                        case 3
                            E = loadStruct_mdl.E.(['model_',num2str(mod_dispE.mdl_n)]);
                            modulators.(curr_oEScale_name)(:,idx_oDispE_Epred) = E;
                        case 4
                            E_pred = loadStruct_mdl.E_pred.(['model_',num2str(mod_dispE.mdl_n)]);
                            modulators.(curr_oEScale_name)(:,idx_oDispE_Epred) = E_pred;
                    end
                end
                
                %% cost term
                if mod_dispE.mdl_cost ~= 0
                    Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                    idx_oDispE_modelCost = Nmod.(curr_oEScale_name);
                    modName.(curr_oEScale_name)(idx_oDispE_modelCost) = {'model_cost'};
                    switch mod_dispE.mdl_cost
                        case 1
                            modulators.(curr_oEScale_name)(:,idx_oDispE_modelCost) = cost;
                        case 2
                            modulators.(curr_oEScale_name)(:,idx_oDispE_modelCost) = cost_pred;
                        case 3
                            modulators.(curr_oEScale_name)(:,idx_oDispE_modelCost) = cost_var;
                    end
                end
                
                %% benefit term
                if mod_dispE.mdl_benefit ~= 0
                    Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                    idx_oDispE_modelBenef = Nmod.(curr_oEScale_name);
                    modName.(curr_oEScale_name)(idx_oDispE_modelBenef) = {'model_benefit'};
                    switch mod_dispE.mdl_benefit
                        case 1
                            modulators.(curr_oEScale_name)(:,idx_oDispE_modelBenef) = benefit;
                        case 2
                            modulators.(curr_oEScale_name)(:,idx_oDispE_modelBenef) = benefit_pred;
                        case 3
                            modulators.(curr_oEScale_name)(:,idx_oDispE_modelBenef) = benefit_var;
                        case 4
                            modulators.(curr_oEScale_name)(:,idx_oDispE_modelBenef) = expected_payoff;
                    end
                end
                
                %% EV term
                if mod_dispE.mdl_EV ~= 0
                    Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                    idx_oDispE_modelEV = Nmod.(curr_oEScale_name);
                    modName.(curr_oEScale_name)(idx_oDispE_modelEV) = {'model_EV'};
                    switch mod_dispE.mdl_EV
                        case 1     
                            modulators.(curr_oEScale_name)(:,idx_oDispE_modelEV) = EV;
                        case 2
                            modulators.(curr_oEScale_name)(:,idx_oDispE_modelEV) = EV_pred;
                    end
                end
                
            end
            
            %% RT
            if mod_dispE.RT_fp ~= 0 && ~ismember(mod_dispE.RT_fp,[2])
                Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                idx_oDispE_RT_fp = Nmod.(curr_oEScale_name);
                switch mod_dispE.RT_fp
                    case 1
                        modName.(curr_oEScale_name)(idx_oDispE_RT_fp) = {'RT_fp'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_RT_fp) = mod.(trialType_oEScale).RT_fp;
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% luminance
            if mod_dispE.lum ~= 0 && ismember(mod_dispE.lum, [1,3,5,7])
                Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                idx_oDispE_lum = Nmod.(curr_oEScale_name);
                switch mod_dispE.lum
                    case 1
                        modName.(curr_oEScale_name)(idx_oDispE_lum) = {'luminance'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_lum) = mod.(trialType_oEScale).lum_inc;
                    case 3
                        modName.(curr_oEScale_name)(idx_oDispE_lum) = {'luminance'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_lum) = nanzscore(mod.(trialType_oEScale).lum_inc);
                    otherwise
                        error('case not ready yet');
                end
            end
            
            %% trial number
            if mod_dispE.trialN ~= 0 && ismember(mod_dispE.trialN,[1,3])
                Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                idx_oDispE_trialN = Nmod.(curr_oEScale_name);
                switch mod_dispE.trialN
                    case 1
                        modName.(curr_oEScale_name)(idx_oDispE_trialN) = {'trialN'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_trialN) = mod.(trialType_oEScale).trialN;
                    case 3
                        modName.(curr_oEScale_name)(idx_oDispE_trialN) = {'trialN'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_trialN) = nanzscore(mod.(trialType_oEScale).trialN);
                end
            end
            
            %% incentive*trial number
            if mod_dispE.inc_x_trialN ~= 0
                Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                idx_oDispE_incXtrialN = Nmod.(curr_oEScale_name);
                switch mod_dispE.inc_x_trialN
                    case 1
                        modName.(curr_oEScale_name)(idx_oDispE_incXtrialN) = {'inc_x_trialN'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_incXtrialN) = (mod.(trialType_oEScale).absIncentiveRank).*(mod.(trialType_oEScale).trialN);
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% predicted ressource level
            if mod_dispE.E_pred ~= 0 && ismember(mod_dispE.E_pred,[1,2])
                Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                idx_oDispE_Epred = Nmod.(curr_oEScale_name);
                loadStruct_mdl   = load([subj_analysis_folder, filesep, 'GS_model_sub',subid,'_grip_run',run_nm,'.mat']);
                modName.(curr_oEScale_name)(idx_oDispE_Epred) = {'E_pred'};
                switch mod_dispE.E_pred
                    case 1
                        E = loadStruct_mdl.E.(['model_',num2str(mod_dispE.mdl_n)]);
                        modulators.(curr_oEScale_name)(:,idx_oDispE_Epred) = E;
                    case 2
                        E_pred = loadStruct_mdl.E_pred.(['model_',num2str(mod_dispE.mdl_n)]);
                        modulators.(curr_oEScale_name)(:,idx_oDispE_Epred) = E_pred;
                end
            end
            
            %% performance
            if mod_dispE.perf ~= 0
                Nmod.(curr_oEScale_name) = Nmod.(curr_oEScale_name) + 1;
                idx_oDispE_perf = Nmod.(curr_oEScale_name);
                switch mod_dispE.perf
                    case 1 % 1-100
                        modName.(curr_oEScale_name)(idx_oDispE_perf) = {'perf'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_perf) = mod.(trialType_oEScale).perf;
                    case 2 % 0-1
                        modName.(curr_oEScale_name)(idx_oDispE_perf) = {'perf'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_perf) = mod.(trialType_oEScale).perf./100;
                    case 3
                        modName.(curr_oEScale_name)(idx_oDispE_perf) = {'perf'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_perf) = mod.(trialType_oEScale).perf_orth;
                    case 4
                        loadStruct_mdl   = load([subj_analysis_folder, filesep, 'GS_model_sub',subid,'_grip_run',run_nm,'.mat']);
                        perf_pred = loadStruct_mdl.perf_pred.(['model_',num2str(mod_dispE.mdl_n)]);
                        modName.(curr_oEScale_name)(idx_oDispE_perf) = {'perf'};
                        modulators.(curr_oEScale_name)(:,idx_oDispE_perf) = perf_pred;
                end
            end
            
        end % condition
    case 0 % nothing happens because not included in the GLM
    otherwise
        error(['case o_dispE = ',num2str(o_dispE),' not ready']);
end % o_dispE

%% onset effort performance
o_perfE       = gripRPprm.o_perfE;
% dur_perfE     = gripRPprm.dur_perfE;
% mod_perfE     = gripRPprm.mod_perfE;
if o_perfE ~= 0
    error('modelling effort performance not ready yet');
end

%% onset feedback
o_fbk       = gripRPprm.o_fbk;
dur_fbk     = gripRPprm.dur_fbk;
mod_fbk     = gripRPprm.mod_fbk;
switch o_fbk
    case {1,2}
        if o_fbk == 1 % all trials pooled
            curr_oFbk_names = {'G_feedback'};
            trialTypes_oFbk = {'all'};
        elseif o_fbk == 2 % to gain/to lose trials separated
            curr_oFbk_names = {'G_feedback_toGain','G_feedback_toLose'};
            trialTypes_oFbk = {'toGain','toLose'};
        end
        n_cond = n_cond + length(curr_oFbk_names);
        cond_nm = [cond_nm, curr_oFbk_names];
        
        for iCond_oFbk = 1:length(curr_oFbk_names)
            curr_oFbk_name = curr_oFbk_names{iCond_oFbk};
            trialType_oFbk = trialTypes_oFbk{iCond_oFbk};
            
            onsetCond.(curr_oFbk_name) = onset.(trialType_oFbk).feedback;
            
            %% duration
            switch dur_fbk
                case 0 % stick
                    durCond.(curr_oFbk_name) = 0;
                case 1 % boxcar from effort scale to feedback
                    durCond.(curr_oFbk_name) = duration.(trialType_oFbk).feedback;
            end
            
            %% modulators
            Nmod.(curr_oFbk_name) = 0;
            
            %% G/L condition
            if mod_fbk.GL_cond ~= 0
                Nmod.(curr_oFbk_name) = Nmod.(curr_oFbk_name) + 1;
                idx_oFbk_GLcond = Nmod.(curr_oFbk_name);
                switch mod_fbk.GL_cond
                    case 1 % 0/1
                        modName.(curr_oFbk_name)(idx_oFbk_GLcond) = {'GL_cond'};
                        GL_cond = mod.(trialType_oFbk).trialValence;
                        GL_cond(GL_cond == -1) = 0; % convert into 0/1 variable
                        modulators.(curr_oFbk_name)(:,idx_oFbk_GLcond) = GL_cond;
                    case 2 % -1/1
                        modName.(curr_oFbk_name)(idx_oFbk_GLcond) = {'GL_cond'};
                        modulators.(curr_oFbk_name)(:,idx_oFbk_GLcond) = mod.(trialType_oFbk).trialValence;
                    otherwise
                        error('case not ready yet');
                end
            end
            
            %% luminance
            
            %% trial number
            if mod_fbk.trialN ~= 0 && ismember(mod_fbk.trialN,[2,4])
                Nmod.(curr_oFbk_name) = Nmod.(curr_oFbk_name) + 1;
                idx_oFbk_trialN = Nmod.(curr_oFbk_name);
                switch mod_fbk.trialN
                    case 2
                        modName.(curr_oFbk_name)(idx_oFbk_trialN) = {'trialN'};
                        modulators.(curr_oFbk_name)(:,idx_oFbk_trialN) = mod.(trialType_oFbk).trialN;
                    case 4
                        modName.(curr_oFbk_name)(idx_oFbk_trialN) = {'trialN'};
                        modulators.(curr_oFbk_name)(:,idx_oFbk_trialN) = nanzscore(mod.(trialType_oFbk).trialN);
                end
            end

            %% predicted ressource level
            if mod_fbk.E_pred ~= 0
                Nmod.(curr_oFbk_name) = Nmod.(curr_oFbk_name) + 1;
                idx_oFbk_Epred = Nmod.(curr_oFbk_name);
                loadStruct_mdl   = load([subj_analysis_folder, filesep, 'GS_model_sub',subid,'_grip_run',run_nm,'.mat']);
                modName.(curr_oFbk_name)(idx_oFbk_Epred) = {'E_pred'};
                switch mod_fbk.E_pred
                    case 5
                        E = loadStruct_mdl.E.(['model_',num2str(mod_fbk.mdl_n)]);
                        modulators.(curr_oFbk_name)(:,idx_oFbk_Epred) = E;
                    case 6
                        E_pred = loadStruct_mdl.E_pred.(['model_',num2str(mod_fbk.mdl_n)]);
                        modulators.(curr_oFbk_name)(:,idx_oFbk_Epred) = E_pred;
                end
            end
            
            %% feedback
            if mod_fbk.fbk ~= 0
                switch mod_fbk.fbk
                    case 1
                        Nmod.(curr_oFbk_name) = Nmod.(curr_oFbk_name) + 1;
                        idx_oFbk_Fbk = Nmod.(curr_oFbk_name);
                        modName.(curr_oFbk_name)(idx_oFbk_Fbk) = {'trial feedback'};
                        modulators.(curr_oFbk_name)(:,idx_oFbk_Fbk) = mod.(trialType_oFbk).gain;
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% total feedback
            if mod_fbk.totalGain ~= 0
                switch mod_fbk.totalGain
                    case 1
                        Nmod.(curr_oFbk_name) = Nmod.(curr_oFbk_name) + 1;
                        idx_oFbk_Fbk = Nmod.(curr_oFbk_name);
                        modName.(curr_oFbk_name)(idx_oFbk_Fbk) = {'total gain'};
                        modulators.(curr_oFbk_name)(:,idx_oFbk_Fbk) = mod.(trialType_oFbk).totalFeedback;
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% predicted ressource level
            if mod_fbk.E_pred ~= 0
                Nmod.(curr_oFbk_name) = Nmod.(curr_oFbk_name) + 1;
                idx_oFbk_Epred = Nmod.(curr_oFbk_name);
                loadStruct_mdl   = load([subj_analysis_folder, filesep, 'GS_model_sub',subid,'_grip_run',run_nm,'.mat']);
                modName.(curr_oFbk_name)(idx_oFbk_Epred) = {'E_pred'};
                switch mod_fbk.E_pred
                    case 1
                        E = loadStruct_mdl.E.(['model_',num2str(mod_fbk.mdl_n)]);
                        modulators.(curr_oFbk_name)(:,idx_oFbk_Epred) = E;
                    case 2
                        E_pred = loadStruct_mdl.E_pred.(['model_',num2str(mod_fbk.mdl_n)]);
                        modulators.(curr_oFbk_name)(:,idx_oFbk_Epred) = E_pred;
                end
            end
            
            %% performance
            if mod_fbk.perf ~= 0
                Nmod.(curr_oFbk_name) = Nmod.(curr_oFbk_name) + 1;
                idx_oFbk_perf = Nmod.(curr_oFbk_name);
                switch mod_fbk.perf
                    case 1 % 1-100
                        modName.(curr_oFbk_name)(idx_oFbk_perf) = {'perf'};
                        modulators.(curr_oFbk_name)(:,idx_oFbk_perf) = mod.(trialType_oFbk).perf;
                    case 2 % 0-1
                        modName.(curr_oFbk_name)(idx_oFbk_perf) = {'perf'};
                        modulators.(curr_oFbk_name)(:,idx_oFbk_perf) = mod.(trialType_oFbk).perf./100;
                    case 3
                        modName.(curr_oFbk_name)(idx_oFbk_perf) = {'perf'};
                        modulators.(curr_oFbk_name)(:,idx_oFbk_perf) = mod.(trialType_oFbk).perf_orth;
                    case 4
                        loadStruct_mdl   = load([subj_analysis_folder, filesep, 'GS_model_sub',subid,'_grip_run',run_nm,'.mat']);
                        perf_pred = loadStruct_mdl.perf_pred.(['model_',num2str(mod_fbk.mdl_n)]);
                        modName.(curr_oFbk_name)(idx_oFbk_perf) = {'perf'};
                        modulators.(curr_oFbk_name)(:,idx_oFbk_perf) = perf_pred;
                end
            end
            
            %% ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                Nmod.(curr_oFbk_name) = Nmod.(curr_oFbk_name) + 1;
                idx_oFbk_ROI = Nmod.(curr_oFbk_name);
                modName.(curr_oFbk_name)(idx_oFbk_ROI) = {[]};
                % get the data
                trialN_idx = mod.(trialType_oFbk).trialN; % trial index
                ROI_nm          = mod_fbk.ROI_activity_ROI_nm; % ROI name
                GLM_ROI_nm      = mod_fbk.ROI_activity_GLM;
                ROI_period_nm   = mod_fbk.ROI_activity_period; % trial period when ROI extracted
                jRun_nb = (iRun <= 3)*1 + (iRun > 3)*2;
                jRun_nm = ['run',num2str(jRun_nb)];
                ROI_data = getfield(load([subj_analysis_folder,'ROI',filesep,'GLM',GLM_ROI_nm,'_',ROI_nm,'.mat'],...
                    'sub_ROI_trial_b_trial'),'sub_ROI_trial_b_trial');
                modulators.(curr_oFbk_name)(:,idx_oFbk_ROI) = ROI_data.(ROI_nm).G.(ROI_period_nm).(jRun_nm).(sub_nm)(trialN_idx);
            end
            
            %% trial number
            if mod_fbk.trialN ~= 0 && ismember(mod_fbk.trialN,[1,3])
                Nmod.(curr_oFbk_name) = Nmod.(curr_oFbk_name) + 1;
                idx_oFbk_trialN = Nmod.(curr_oFbk_name);
                switch mod_fbk.trialN
                    case 1
                        modName.(curr_oFbk_name)(idx_oFbk_trialN) = {'trialN'};
                        modulators.(curr_oFbk_name)(:,idx_oFbk_trialN) = mod.(trialType_oFbk).trialN;
                    case 3
                        modName.(curr_oFbk_name)(idx_oFbk_trialN) = {'trialN'};
                        modulators.(curr_oFbk_name)(:,idx_oFbk_trialN) = nanzscore(mod.(trialType_oFbk).trialN);
                end
            end
            
        end
end

%% onset cross
o_cross = gripRPprm.o_cross;
dur_cross = gripRPprm.dur_cross;
switch o_cross
    case 1 % onset for fixation cross across all trials and conditions
        n_cond = n_cond + 1;
        cond_nm = [cond_nm, 'G_cross'];
        onsetCond.G_cross = onset.ITI;
        
        % duration
        switch dur_cross
            case 0 % stick
                durCond.(curr_oFbk_name) = 0;
            case 1 % boxcar from effort scale to feedback
                durCond.(curr_oFbk_name) = duration.ITI;
        end
        
        Nmod.G_cross = 0;
end

%% preparing matlab batch
[ matlabbatch ] = First_level_MS2_matlabbatch_onsets_modulators( matlabbatch,...
    cond_nm, onsetCond, durCond,...
    sub_idx, iRun, Nmod,...
    task_id, modName, modulators,...
    gal);


end % function