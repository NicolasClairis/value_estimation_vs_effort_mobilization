function [ jCon, con_vec, con_names ] = MS2_contrasts_combinations_Stroop( jCon, con_vec, con_names,...
    posCon_nm, n_regs, stroopRPprm)
% [ jCon, con_vec, con_names ] = MS2_contrasts_combinations_Stroop( jCon, con_vec, con_names,...
%     posCon_nm, n_regs, stroopRPprm) combine similar contrasts in the
%     Stroop task.

%% incentive
S_o_inc = stroopRPprm.o_inc;
S_mod_inc = stroopRPprm.mod_inc;
switch S_o_inc
    case 1
        % one contrast for run 1 and one other for run 2
        S_oInc_idx = strcmp(con_names,['S_o_inc',posCon_nm]);
        S_oInc_vec = con_vec(S_oInc_idx,:);
        S_oInc_r1_idx = find(S_oInc_vec ~= 0, 1, 'first');
        S_oInc_r2_idx = find(S_oInc_vec ~= 0, 1, 'last');
        
        % run 1
        curr_vec_nm = 'S_o_inc_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'S_o_inc_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run 1
        curr_vec_nm = 'S_o_inc_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_r2_idx) = 1;
        con_vec(jCon, S_oInc_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % model variables compute EV = benefit - cost
        if S_mod_inc.mdl_n ~= 0 && S_mod_inc.mdl_cost ~= 0 && S_mod_inc.mdl_benefit ~= 0
            curr_vec_nm = ['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_Benefit_min_Cost'];
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_mdlBenef_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit',posCon_nm]);
            S_inc_mdlCost_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_mdlBenef_idx,:) - con_vec(S_inc_mdlCost_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
    case 2 % trials split according to trial type => contrast across trial types (pool to gain + to lose)
        
        % to gain only
        % one contrast for run 1 and one other for run 2
        S_oInc_toGain_idx = strcmp(con_names,['S_o_inc_toGain',posCon_nm]);
        S_oInc_toGain_vec = con_vec(S_oInc_toGain_idx,:);
        S_oInc_toGain_r1_idx = find(S_oInc_toGain_vec ~= 0, 1, 'first');
        S_oInc_toGain_r2_idx = find(S_oInc_toGain_vec ~= 0, 1, 'last');
        % run 1
        curr_vec_nm = 'S_o_inc_toGain_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_toGain_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'S_o_inc_toGain_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_toGain_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run 1 - to Gain
        curr_vec_nm = 'S_o_inc_toGain_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_toGain_r2_idx) = 1;
        con_vec(jCon, S_oInc_toGain_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % to lose only
        % one contrast for run 1 and one other for run 2
        S_oInc_toLose_idx = strcmp(con_names,['S_o_inc_toLose',posCon_nm]);
        S_oInc_toLose_vec = con_vec(S_oInc_toLose_idx,:);
        S_oInc_toLose_r1_idx = find(S_oInc_toLose_vec ~= 0, 1, 'first');
        S_oInc_toLose_r2_idx = find(S_oInc_toLose_vec ~= 0, 1, 'last');
        % run 1
        curr_vec_nm = 'S_o_inc_toLose_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_toLose_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'S_o_inc_toLose_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_toLose_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run 1 - to Lose
        curr_vec_nm = 'S_o_inc_toLose_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_toLose_r2_idx) = 1;
        con_vec(jCon, S_oInc_toLose_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % pool both trial types
        % run 1
        curr_vec_nm = 'S_o_inc_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_toGain_r1_idx) = 1/2;
        con_vec(jCon, S_oInc_toLose_r1_idx) = 1/2;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'S_o_inc_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_toGain_r2_idx) = 1/2;
        con_vec(jCon, S_oInc_toLose_r2_idx) = 1/2;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run 1
        curr_vec_nm = 'S_o_inc_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_toGain_r2_idx) = 1;
        con_vec(jCon, S_oInc_toLose_r2_idx) = 1;
        con_vec(jCon, S_oInc_toGain_r1_idx) = -1;
        con_vec(jCon, S_oInc_toLose_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % compare incentive Gain vs incentive Loss period
        % pool gain and loss condition for incentive modulation
        curr_vec_nm = 'S_o_inc_Gain_min_Loss';
        
        % positive contrast
        jCon = jCon + 1;
        S_Oinc_toGain_idx = strcmp(con_names,['S_o_inc_toGain',posCon_nm]);
        S_Oinc_toLose_idx = strcmp(con_names,['S_o_inc_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(S_Oinc_toGain_idx,:) - con_vec(S_Oinc_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % reward type
        if S_mod_inc.R_type ~= 0
            curr_vec_nm = 'S_mod_inc_R_type';
            % positive contrast
            jCon = jCon + 1;
            S_R_type_toGain_idx = strcmp(con_names,['S_mod_inc_R_type_toGain',posCon_nm]);
            S_R_type_toLose_idx = strcmp(con_names,['S_mod_inc_R_type_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_R_type_toGain_idx,:) + con_vec(S_R_type_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % comparison gain and loss condition for incentive modulation
            curr_vec_nm = 'S_mod_inc_R_type_Gain_min_Loss';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(S_R_type_toGain_idx,:) - con_vec(S_R_type_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive
        if S_mod_inc.inc ~= 0
            curr_vec_nm = 'S_mod_inc_inc';
            % positive contrast
            jCon = jCon + 1;
            S_inc_toGain_idx = strcmp(con_names,['S_mod_inc_inc_toGain',posCon_nm]);
            S_inc_toLose_idx = strcmp(con_names,['S_mod_inc_inc_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_toGain_idx,:) + con_vec(S_inc_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % comparison gain and loss condition for incentive modulation
            curr_vec_nm = 'S_mod_inc_inc_Gain_min_Loss';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_toGain_idx,:) - con_vec(S_inc_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive bis
        if S_mod_inc.inc_bis ~= 0
            curr_vec_nm = 'S_mod_inc_inc_bis';
            % positive contrast
            jCon = jCon + 1;
            S_inc_bis_toGain_idx = strcmp(con_names,['S_mod_inc_inc_bis_toGain',posCon_nm]);
            S_inc_bis_toLose_idx = strcmp(con_names,['S_mod_inc_inc_bis_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_bis_toGain_idx,:) + con_vec(S_inc_bis_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % comparison gain and loss condition for incentive modulation
            curr_vec_nm = 'S_mod_inc_inc_bis_Gain_min_Loss';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_bis_toGain_idx,:) - con_vec(S_inc_bis_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ROI activity
        if S_mod_inc.ROI_activity_yn ~= 0
            ROI_full_nm = ['S_mod_inc_GLM',S_mod_inc.ROI_activity_GLM,...
                    '_',S_mod_inc.ROI_activity_period,'_period_',S_mod_inc.ROI_activity_ROI_nm,'_activity'];
            curr_vec_nm = ROI_full_nm;
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_ROIactivity_toGain_idx = strcmp(con_names,[ROI_full_nm,'_toGain',posCon_nm]);
            S_inc_ROIactivity_toLose_idx = strcmp(con_names,[ROI_full_nm,'_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_ROIactivity_toGain_idx,:) + con_vec(S_inc_ROIactivity_toLose_idx,:); % pool effort to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % number of pairs solved
        if S_mod_inc.n_pairs_solved ~= 0
            curr_vec_nm = 'S_mod_inc_n_pairs_solved';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_nPairs_solved_toGain_idx = strcmp(con_names,['S_mod_inc_n_pairs_solved_toGain',posCon_nm]);
            S_inc_nPairs_solved_toLose_idx = strcmp(con_names,['S_mod_inc_n_pairs_solved_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_nPairs_solved_toGain_idx,:) + con_vec(S_inc_nPairs_solved_toLose_idx,:); % pool nb pairs solved to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % errors
        if S_mod_inc.n_errors ~= 0
            curr_vec_nm = 'S_mod_inc_n_errors';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_n_errors_toGain_idx = strcmp(con_names,['S_mod_inc_n_errors_toGain',posCon_nm]);
            S_inc_n_errors_toLose_idx = strcmp(con_names,['S_mod_inc_n_errors_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_n_errors_toGain_idx,:) + con_vec(S_inc_n_errors_toLose_idx,:); % pool errors nber to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % model variables
        if S_mod_inc.mdl_n ~= 0
            
            % pool effort to gain and to lose trials
            % cost
            if S_mod_inc.mdl_cost ~= 0
                curr_vec_nm = ['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost'];
                
                % positive contrast
                jCon = jCon + 1;
                S_inc_cost_toGain_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost_toGain',posCon_nm]);
                S_inc_cost_toLose_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost_toLose',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_inc_cost_toGain_idx,:) + con_vec(S_inc_cost_toLose_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % benefit
            if S_mod_inc.mdl_benefit ~= 0
                curr_vec_nm = ['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit'];
                
                % positive contrast
                jCon = jCon + 1;
                S_inc_benefit_toGain_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit_toGain',posCon_nm]);
                S_inc_benefit_toLose_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit_toLose',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_inc_benefit_toGain_idx,:) + con_vec(S_inc_benefit_toLose_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % EV
            if S_mod_inc.mdl_EV ~= 0
                curr_vec_nm = ['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_EV'];
                
                % positive contrast
                jCon = jCon + 1;
                S_inc_EV_toGain_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_EV_toGain',posCon_nm]);
                S_inc_EV_toLose_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_EV_toLose',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_inc_EV_toGain_idx,:) + con_vec(S_inc_EV_toLose_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % model variables compute EV = benefit - cost
            if S_mod_inc.mdl_cost ~= 0 && S_mod_inc.mdl_benefit ~= 0
                
                % to gain
                curr_vec_nm = ['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_Benefit_min_Cost_toGain'];
                % positive contrast
                jCon = jCon + 1;
                S_inc_toGain_mdlBenef_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit_toGain',posCon_nm]);
                S_inc_toGain_mdlCost_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost_toGain',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_inc_toGain_mdlBenef_idx,:) - con_vec(S_inc_toGain_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
                
                % to lose
                curr_vec_nm = ['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_Benefit_min_Cost_toLose'];
                % positive contrast
                jCon = jCon + 1;
                S_inc_toLose_mdlBenef_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit_toLose',posCon_nm]);
                S_inc_toLose_mdlCost_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost_toLose',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_inc_toLose_mdlBenef_idx,:) - con_vec(S_inc_toLose_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
                
                % pool gain and losses
                curr_vec_nm = ['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_Benefit_min_Cost'];
                % positive contrast
                jCon = jCon + 1;
                S_inc_mdlBenef_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit',posCon_nm]);
                S_inc_mdlCost_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_inc_mdlBenef_idx,:) - con_vec(S_inc_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
        end
        
        % RT first pair
        if S_mod_inc.RT_fp ~= 0
            curr_vec_nm = 'S_mod_inc_RT_fp';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_RTfp_toGain_idx = strcmp(con_names,['S_mod_inc_RT_fp_toGain',posCon_nm]);
            S_inc_RTfp_toLose_idx = strcmp(con_names,['S_mod_inc_RT_fp_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_RTfp_toGain_idx,:) + con_vec(S_inc_RTfp_toLose_idx,:); % pool RTfp to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % mean(RT) except first pair
        if S_mod_inc.RT_mRT ~= 0
            curr_vec_nm = 'S_mod_inc_RT_mRT';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_mRT_toGain_idx = strcmp(con_names,['S_mod_inc_RT_mRT_toGain',posCon_nm]);
            S_inc_mRT_toLose_idx = strcmp(con_names,['S_mod_inc_RT_mRT_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_mRT_toGain_idx,:) + con_vec(S_inc_mRT_toLose_idx,:); % pool mRT to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % total gain until now
        if S_mod_inc.totalGain_prev ~= 0
            curr_vec_nm = 'S_mod_inc_totalGain_prev';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_totalGainPrev_toGain_idx = strcmp(con_names,['S_mod_inc_totalGain_prev_toGain',posCon_nm]);
            S_inc_totalGainPrev_toLose_idx = strcmp(con_names,['S_mod_inc_totalGain_prev_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_totalGainPrev_toGain_idx,:) + con_vec(S_inc_totalGainPrev_toLose_idx,:); % pool to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % sum(performance previous trials)
         if S_mod_inc.sumPerfPrev ~= 0
            curr_vec_nm = 'S_mod_inc_sumPerfPrev';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_sumPerfPrev_toGain_idx = strcmp(con_names,['S_mod_inc_sumPerfPrev_toGain',posCon_nm]);
            S_inc_sumPerfPrev_toLose_idx = strcmp(con_names,['S_mod_inc_sumPerfPrev_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_sumPerfPrev_toGain_idx,:) + con_vec(S_inc_sumPerfPrev_toLose_idx,:); % pool to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % trial number
        if S_mod_inc.trialN ~= 0
            curr_vec_nm = 'S_mod_inc_trialN';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_trialN_toGain_idx = strcmp(con_names,['S_mod_inc_trialN_toGain',posCon_nm]);
            S_inc_trialN_toLose_idx = strcmp(con_names,['S_mod_inc_trialN_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_trialN_toGain_idx,:) + con_vec(S_inc_trialN_toLose_idx,:); % pool RTfp to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ressources
        if S_mod_inc.X_pred ~= 0
            curr_vec_nm = 'S_mod_inc_X_pred';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_X_pred_toGain_idx = strcmp(con_names,['S_mod_inc_X_pred_toGain',posCon_nm]);
            S_inc_X_pred_toLose_idx = strcmp(con_names,['S_mod_inc_X_pred_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_X_pred_toGain_idx,:) + con_vec(S_inc_X_pred_toLose_idx,:); % pool trialN to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % performance
        if S_mod_inc.perf ~= 0
            curr_vec_nm = 'S_mod_inc_perf';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_perf_toGain_idx = strcmp(con_names,['S_mod_inc_perf_toGain',posCon_nm]);
            S_inc_perf_toLose_idx = strcmp(con_names,['S_mod_inc_perf_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_perf_toGain_idx,:) + con_vec(S_inc_perf_toLose_idx,:); % pool trialN to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
    case 3 % trials split according to trial type => contrast across trial types (pool error+no-error trials)
        
        % no-error trials only
        % one contrast for run 1 and one other for run 2
        S_oInc_noErrorTrials_idx = strcmp(con_names,['S_o_inc_noErrorTrials',posCon_nm]);
        S_oInc_noErrorTrials_vec = con_vec(S_oInc_noErrorTrials_idx,:);
        S_oInc_noErrorTrials_r1_idx = find(S_oInc_noErrorTrials_vec ~= 0, 1, 'first');
        S_oInc_noErrorTrials_r2_idx = find(S_oInc_noErrorTrials_vec ~= 0, 1, 'last');
        % run 1
        curr_vec_nm = 'S_o_inc_noErrorTrials_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_noErrorTrials_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'S_o_inc_noErrorTrials_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_noErrorTrials_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run 1 - no Error trials
        curr_vec_nm = 'S_o_inc_noErrorTrials_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_noErrorTrials_r2_idx) = 1;
        con_vec(jCon, S_oInc_noErrorTrials_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % error trials only
        % one contrast for run 1 and one other for run 2
        S_oInc_ErrorTrials_idx = strcmp(con_names,['S_o_inc_ErrorTrials',posCon_nm]);
        S_oInc_ErrorTrials_vec = con_vec(S_oInc_ErrorTrials_idx,:);
        S_oInc_ErrorTrials_r1_idx = find(S_oInc_ErrorTrials_vec ~= 0, 1, 'first');
        S_oInc_ErrorTrials_r2_idx = find(S_oInc_ErrorTrials_vec ~= 0, 1, 'last');
        % run 1
        curr_vec_nm = 'S_o_inc_ErrorTrials_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_ErrorTrials_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'S_o_inc_ErrorTrials_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_ErrorTrials_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run 1 - error trials
        curr_vec_nm = 'S_o_inc_ErrorTrials_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_ErrorTrials_r2_idx) = 1;
        con_vec(jCon, S_oInc_ErrorTrials_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % pool both trial types
        % run 1
        curr_vec_nm = 'S_o_inc_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_noErrorTrials_r1_idx) = 1/2;
        con_vec(jCon, S_oInc_ErrorTrials_r1_idx) = 1/2;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'S_o_inc_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_noErrorTrials_r2_idx) = 1/2;
        con_vec(jCon, S_oInc_ErrorTrials_r2_idx) = 1/2;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run 1 - both trial types
        curr_vec_nm = 'S_o_inc_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oInc_noErrorTrials_r2_idx)  = 1;
        con_vec(jCon, S_oInc_ErrorTrials_r2_idx)    = 1;
        con_vec(jCon, S_oInc_noErrorTrials_r1_idx)  = -1;
        con_vec(jCon, S_oInc_ErrorTrials_r1_idx)    = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % GL condition
        if S_mod_inc.GL_cond ~= 0
            curr_vec_nm = 'S_mod_inc_GL_cond_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_GL_cond_noErrorTrials_idx = strcmp(con_names,['S_mod_inc_GL_cond_noErrorTrials',posCon_nm]);
            S_inc_GL_cond_ErrorTrials_idx = strcmp(con_names,['S_mod_inc_GL_cond_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_GL_cond_noErrorTrials_idx,:) + con_vec(S_inc_GL_cond_ErrorTrials_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % reward type
        if S_mod_inc.R_type ~= 0
            curr_vec_nm = 'S_mod_inc_R_type_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_R_type_noErrorTrials_idx = strcmp(con_names,['S_mod_inc_R_type_noErrorTrials',posCon_nm]);
            S_R_type_ErrorTrials_idx = strcmp(con_names,['S_mod_inc_R_type_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_R_type_noErrorTrials_idx,:) + con_vec(S_R_type_ErrorTrials_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive
        if S_mod_inc.inc ~= 0
            curr_vec_nm = 'S_mod_inc_inc_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_noErrorTrials_idx = strcmp(con_names,['S_mod_inc_inc_noErrorTrials',posCon_nm]);
            S_inc_ErrorTrials_idx = strcmp(con_names,['S_mod_inc_inc_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_noErrorTrials_idx,:) + con_vec(S_inc_ErrorTrials_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive bis
        if S_mod_inc.inc_bis ~= 0
            curr_vec_nm = 'S_mod_inc_inc_bis_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_bis_noErrorTrials_idx = strcmp(con_names,['S_mod_inc_inc_bis_noErrorTrials',posCon_nm]);
            S_inc_bis_ErrorTrials_idx = strcmp(con_names,['S_mod_inc_inc_bis_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_bis_noErrorTrials_idx,:) + con_vec(S_inc_bis_ErrorTrials_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ROI activity
        if S_mod_inc.ROI_activity_yn ~= 0
            ROI_full_nm = ['S_mod_inc_GLM',S_mod_inc.ROI_activity_GLM,...
                    '_',S_mod_inc.ROI_activity_period,'_period_',S_mod_inc.ROI_activity_ROI_nm,'_activity'];
            curr_vec_nm = [ROI_full_nm,'_ETnoETpool'];
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_ROIactivity_noErrorTrials_idx = strcmp(con_names,[ROI_full_nm,'_noErrorTrials',posCon_nm]);
            S_inc_ROIactivity_ErrorTrials_idx = strcmp(con_names,[ROI_full_nm,'_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_ROIactivity_noErrorTrials_idx,:) + con_vec(S_inc_ROIactivity_ErrorTrials_idx,:); % pool effort to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % number of pairs solved
        if S_mod_inc.n_pairs_solved ~= 0
            curr_vec_nm = 'S_mod_inc_n_pairs_solved_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_nPairs_solved_noErrorTrials_idx = strcmp(con_names,['S_mod_inc_n_pairs_solved_noErrorTrials',posCon_nm]);
            S_inc_nPairs_solved_ErrorTrials_idx = strcmp(con_names,['S_mod_inc_n_pairs_solved_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_nPairs_solved_noErrorTrials_idx,:) + con_vec(S_inc_nPairs_solved_ErrorTrials_idx,:); % pool nb pairs solved to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % model variables
        if S_mod_inc.mdl_n ~= 0
            
            % pool effort to gain and to lose trials
            % cost
            if S_mod_inc.mdl_cost ~= 0
                curr_vec_nm = ['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost'];
                
                % positive contrast
                jCon = jCon + 1;
                S_inc_cost_noErrorTrials_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost_noErrorTrials',posCon_nm]);
                S_inc_cost_ErrorTrials_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost_ErrorTrials',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_inc_cost_noErrorTrials_idx,:) + con_vec(S_inc_cost_ErrorTrials_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % benefit
            if S_mod_inc.mdl_benefit ~= 0
                curr_vec_nm = ['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit'];
                
                % positive contrast
                jCon = jCon + 1;
                S_inc_benefit_noErrorTrials_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit_noErrorTrials',posCon_nm]);
                S_inc_benefit_ErrorTrials_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit_ErrorTrials',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_inc_benefit_noErrorTrials_idx,:) + con_vec(S_inc_benefit_ErrorTrials_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % EV
            if S_mod_inc.mdl_EV ~= 0
                curr_vec_nm = ['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_EV'];
                
                % positive contrast
                jCon = jCon + 1;
                S_inc_EV_noErrorTrials_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_EV_noErrorTrials',posCon_nm]);
                S_inc_EV_ErrorTrials_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_EV_ErrorTrials',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_inc_EV_noErrorTrials_idx,:) + con_vec(S_inc_EV_ErrorTrials_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % model variables compute EV = benefit - cost
            if S_mod_inc.mdl_cost ~= 0 && S_mod_inc.mdl_benefit ~= 0
                
                % no error trials
                curr_vec_nm = ['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_Benefit_min_Cost_noErrorTrials'];
                % positive contrast
                jCon = jCon + 1;
                S_inc_noErrorTrials_mdlBenef_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit_noErrorTrials',posCon_nm]);
                S_inc_noErrorTrials_mdlCost_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost_noErrorTrials',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_inc_noErrorTrials_mdlBenef_idx,:) - con_vec(S_inc_noErrorTrials_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
                
                % error trials
                curr_vec_nm = ['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_Benefit_min_Cost_ErrorTrials'];
                % positive contrast
                jCon = jCon + 1;
                S_inc_ErrorTrials_mdlBenef_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit_ErrorTrials',posCon_nm]);
                S_inc_ErrorTrials_mdlCost_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost_ErrorTrials',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_inc_ErrorTrials_mdlBenef_idx,:) - con_vec(S_inc_ErrorTrials_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
                
                % pool error and no-error
                curr_vec_nm = ['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_Benefit_min_Cost_ETnoETpool'];
                % positive contrast
                jCon = jCon + 1;
                S_inc_mdlBenef_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit',posCon_nm]);
                S_inc_mdlCost_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_inc_mdlBenef_idx,:) - con_vec(S_inc_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
        end
        
        % RT first pair
        if S_mod_inc.RT_fp ~= 0
            curr_vec_nm = 'S_mod_inc_RT_fp_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_RTfp_noErrorTrials_idx = strcmp(con_names,['S_mod_inc_RT_fp_noErrorTrials',posCon_nm]);
            S_inc_RTfp_ErrorTrials_idx = strcmp(con_names,['S_mod_inc_RT_fp_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_RTfp_noErrorTrials_idx,:) + con_vec(S_inc_RTfp_ErrorTrials_idx,:); % pool nb RTfp to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % mean(RT) except first pair
        if S_mod_inc.RT_mRT ~= 0
            curr_vec_nm = 'S_mod_inc_RT_mRT_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_mRT_noErrorTrials_idx = strcmp(con_names,['S_mod_inc_RT_mRT_noErrorTrials',posCon_nm]);
            S_inc_mRT_ErrorTrials_idx = strcmp(con_names,['S_mod_inc_RT_mRT_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_mRT_noErrorTrials_idx,:) + con_vec(S_inc_mRT_ErrorTrials_idx,:); % pool mRT to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % total gain until now
        if S_mod_inc.totalGain_prev ~= 0
            curr_vec_nm = 'S_mod_inc_totalGain_prev_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_totalGainPrev_noErrorTrials_idx = strcmp(con_names,['S_mod_inc_totalGain_prev_noErrorTrials',posCon_nm]);
            S_inc_totalGainPrev_ErrorTrials_idx = strcmp(con_names,['S_mod_inc_totalGain_prev_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_totalGainPrev_noErrorTrials_idx,:) + con_vec(S_inc_totalGainPrev_ErrorTrials_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % sum(performance previous trials)
         if S_mod_inc.sumPerfPrev ~= 0
            curr_vec_nm = 'S_mod_inc_sumPerfPrev_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_sumPerfPrev_noErrorTrials_idx = strcmp(con_names,['S_mod_inc_sumPerfPrev_noErrorTrials',posCon_nm]);
            S_inc_sumPerfPrev_ErrorTrials_idx = strcmp(con_names,['S_mod_inc_sumPerfPrev_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_sumPerfPrev_noErrorTrials_idx,:) + con_vec(S_inc_sumPerfPrev_ErrorTrials_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        
        % trial number
        if S_mod_inc.trialN ~= 0
            curr_vec_nm = 'S_mod_inc_trialN_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_trialN_noErrorTrials_idx = strcmp(con_names,['S_mod_inc_trialN_noErrorTrials',posCon_nm]);
            S_inc_trialN_ErrorTrials_idx = strcmp(con_names,['S_mod_inc_trialN_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_trialN_noErrorTrials_idx,:) + con_vec(S_inc_trialN_ErrorTrials_idx,:); % pool nb RTfp to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ressource
        if S_mod_inc.X_pred ~= 0
            curr_vec_nm = 'S_mod_inc_X_pred_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_X_pred_noErrorTrials_idx = strcmp(con_names,['S_mod_inc_X_pred_noErrorTrials',posCon_nm]);
            S_inc_X_pred_ErrorTrials_idx = strcmp(con_names,['S_mod_inc_X_pred_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_X_pred_noErrorTrials_idx,:) + con_vec(S_inc_X_pred_ErrorTrials_idx,:); % pool trialN to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % performance
        if S_mod_inc.perf ~= 0
            curr_vec_nm = 'S_mod_inc_perf_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_inc_perf_noErrorTrials_idx = strcmp(con_names,['S_mod_inc_perf_noErrorTrials',posCon_nm]);
            S_inc_perf_ErrorTrials_idx = strcmp(con_names,['S_mod_inc_perf_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_inc_perf_noErrorTrials_idx,:) + con_vec(S_inc_perf_ErrorTrials_idx,:); % pool trialN to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
end % o_inc


%% effort scale
S_o_dispE = stroopRPprm.o_dispE;
S_mod_dispE = stroopRPprm.mod_dispE;
switch S_o_dispE
    case 1
        % one contrast for run 1 and one other for run 2
        S_oDispE_idx = strcmp(con_names,['S_o_dispE',posCon_nm]);
        S_oDispE_vec = con_vec(S_oDispE_idx,:);
        S_oDispE_r1_idx = find(S_oDispE_vec ~= 0, 1, 'first');
        S_oDispE_r2_idx = find(S_oDispE_vec ~= 0, 1, 'last');
        
        % run 1
        curr_vec_nm = 'S_o_dispE_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'S_o_dispE_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run 1
        curr_vec_nm = 'S_o_dispE_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_r2_idx) = 1;
        con_vec(jCon, S_oDispE_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % model variables compute EV = benefit - cost
        if S_mod_dispE.mdl_n ~= 0 && S_mod_dispE.mdl_cost ~= 0 && S_mod_dispE.mdl_benefit ~= 0
            curr_vec_nm = ['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_Benefit_min_Cost'];
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_mdlBenef_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit',posCon_nm]);
            S_dispE_mdlCost_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_mdlBenef_idx,:) - con_vec(S_dispE_mdlCost_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
    case 2 % trials split according to trial type => contrast across trial types (pool to gain + to lose)
        % to gain only
        % one contrast for run 1 and one other for run 2
        S_oDispE_toGain_idx = strcmp(con_names,['S_o_dispE_toGain',posCon_nm]);
        S_oDispE_toGain_vec = con_vec(S_oDispE_toGain_idx,:);
        S_oDispE_toGain_r1_idx = find(S_oDispE_toGain_vec ~= 0, 1, 'first');
        S_oDispE_toGain_r2_idx = find(S_oDispE_toGain_vec ~= 0, 1, 'last');
        % run 1
        curr_vec_nm = 'S_o_dispE_toGain_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_toGain_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'S_o_dispE_toGain_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_toGain_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run1 - toGain
        curr_vec_nm = 'S_o_dispE_toGain_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_toGain_r2_idx) = 1;
        con_vec(jCon, S_oDispE_toGain_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % to lose only
        % one contrast for run 1 and one other for run 2
        S_oDispE_toLose_idx = strcmp(con_names,['S_o_dispE_toLose',posCon_nm]);
        S_oDispE_toLose_vec = con_vec(S_oDispE_toLose_idx,:);
        S_oDispE_toLose_r1_idx = find(S_oDispE_toLose_vec ~= 0, 1, 'first');
        S_oDispE_toLose_r2_idx = find(S_oDispE_toLose_vec ~= 0, 1, 'last');
        
        % run 1
        curr_vec_nm = 'S_o_dispE_toLose_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_toLose_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'S_o_dispE_toLose_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_toLose_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run 1 - toLose
        curr_vec_nm = 'S_o_dispE_toLose_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_toLose_r2_idx) = 1;
        con_vec(jCon, S_oDispE_toLose_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % pool both trial types
        % run 1
        curr_vec_nm = 'S_o_dispE_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_toGain_r1_idx) = 1;
        con_vec(jCon, S_oDispE_toLose_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'S_o_dispE_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_toGain_r2_idx) = 1/2;
        con_vec(jCon, S_oDispE_toLose_r2_idx) = 1/2;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 min run 1
        curr_vec_nm = 'S_o_dispE_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_toGain_r2_idx) = 1;
        con_vec(jCon, S_oDispE_toLose_r2_idx) = 1;
        con_vec(jCon, S_oDispE_toGain_r1_idx) = -1;
        con_vec(jCon, S_oDispE_toLose_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % compare incentive Gain vs incentive Loss period
        % pool gain and loss condition for incentive modulation
        curr_vec_nm = 'S_o_dispE_Gain_min_Loss';
        
        % positive contrast
        jCon = jCon + 1;
        S_OdispE_toGain_idx = strcmp(con_names,['S_o_dispE_toGain',posCon_nm]);
        S_OdispE_toLose_idx = strcmp(con_names,['S_o_dispE_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(S_OdispE_toGain_idx,:) - con_vec(S_OdispE_toLose_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % reward type
        if S_mod_dispE.R_type ~= 0
            curr_vec_nm = 'S_mod_dispE_R_type';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_R_type_toGain_idx = strcmp(con_names,['S_mod_dispE_R_type_toGain',posCon_nm]);
            S_dispE_R_type_toLose_idx = strcmp(con_names,['S_mod_dispE_R_type_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_R_type_toGain_idx,:) + con_vec(S_dispE_R_type_toLose_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive
        if S_mod_dispE.inc ~= 0
            curr_vec_nm = 'S_mod_dispE_inc';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_inc_toGain_idx = strcmp(con_names,['S_mod_dispE_inc_toGain',posCon_nm]);
            S_dispE_inc_toLose_idx = strcmp(con_names,['S_mod_dispE_inc_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_inc_toGain_idx,:) + con_vec(S_dispE_inc_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive bis
        if S_mod_dispE.inc_bis ~= 0
            curr_vec_nm = 'S_mod_dispE_inc_bis';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_inc_bis_toGain_idx = strcmp(con_names,['S_mod_dispE_inc_bis_toGain',posCon_nm]);
            S_dispE_inc_bis_toLose_idx = strcmp(con_names,['S_mod_dispE_inc_bis_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_inc_bis_toGain_idx,:) + con_vec(S_dispE_inc_bis_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ROI activity
        if S_mod_dispE.ROI_activity_yn ~= 0
            ROI_full_nm = ['S_mod_dispE_GLM',S_mod_dispE.ROI_activity_GLM,...
                    '_',S_mod_dispE.ROI_activity_period,'_period_',S_mod_dispE.ROI_activity_ROI_nm,'_activity'];
            curr_vec_nm = ROI_full_nm;
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_ROIactivity_toGain_idx = strcmp(con_names,[ROI_full_nm,'_toGain',posCon_nm]);
            S_dispE_ROIactivity_toLose_idx = strcmp(con_names,[ROI_full_nm,'_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_ROIactivity_toGain_idx,:) + con_vec(S_dispE_ROIactivity_toLose_idx,:); % pool effort to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % number of pairs solved
        if S_mod_dispE.n_pairs_solved ~= 0
            curr_vec_nm = 'S_mod_dispE_n_pairs';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_nPairs_toGain_idx = strcmp(con_names,['S_mod_dispE_n_pairs_solved_toGain',posCon_nm]);
            S_dispE_nPairs_toLose_idx = strcmp(con_names,['S_mod_dispE_n_pairs_solved_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_nPairs_toGain_idx,:) + con_vec(S_dispE_nPairs_toLose_idx,:); % pool nb pairs solved to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incongruent pairs solved
        if S_mod_dispE.n_incong ~= 0
            curr_vec_nm = 'S_mod_dispE_n_incong';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_nIncong_toGain_idx = strcmp(con_names,['S_mod_dispE_n_incong_toGain',posCon_nm]);
            S_dispE_nIncong_toLose_idx = strcmp(con_names,['S_mod_dispE_n_incong_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_nIncong_toGain_idx,:) + con_vec(S_dispE_nIncong_toLose_idx,:); % pool incongruent pairs solved to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % distance between numbers of pairs solved
        if S_mod_dispE.n_dist ~= 0
            curr_vec_nm = 'S_mod_dispE_n_dist';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_nDist_toGain_idx = strcmp(con_names,['S_mod_dispE_n_dist_toGain',posCon_nm]);
            S_dispE_nDist_toLose_idx = strcmp(con_names,['S_mod_dispE_n_dist_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_nDist_toGain_idx,:) + con_vec(S_dispE_nDist_toLose_idx,:); % pool dNbers to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % errors
        if S_mod_dispE.n_errors ~= 0
            curr_vec_nm = 'S_mod_dispE_n_errors';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_nErrors_toGain_idx = strcmp(con_names,['S_mod_dispE_n_errors_toGain',posCon_nm]);
            S_dispE_nErrors_toLose_idx = strcmp(con_names,['S_mod_dispE_n_errors_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_nErrors_toGain_idx,:) + con_vec(S_dispE_nErrors_toLose_idx,:); % pool nErrors to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % model variables compute EV = benefit - cost
        % model variables
        if S_mod_dispE.mdl_n ~= 0
            
            % pool effort to gain and to lose trials
            % cost
            if S_mod_dispE.mdl_cost ~= 0
                curr_vec_nm = ['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost'];
                
                % positive contrast
                jCon = jCon + 1;
                S_dispE_cost_toGain_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost_toGain',posCon_nm]);
                S_dispE_cost_toLose_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost_toLose',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_dispE_cost_toGain_idx,:) + con_vec(S_dispE_cost_toLose_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % benefit
            if S_mod_dispE.mdl_benefit ~= 0
                curr_vec_nm = ['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit'];
                
                % positive contrast
                jCon = jCon + 1;
                S_dispE_benefit_toGain_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit_toGain',posCon_nm]);
                S_dispE_benefit_toLose_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit_toLose',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_dispE_benefit_toGain_idx,:) + con_vec(S_dispE_benefit_toLose_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % EV
            if S_mod_dispE.mdl_EV ~= 0
                curr_vec_nm = ['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_EV'];
                
                % positive contrast
                jCon = jCon + 1;
                S_dispE_EV_toGain_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_EV_toGain',posCon_nm]);
                S_dispE_EV_toLose_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_EV_toLose',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_dispE_EV_toGain_idx,:) + con_vec(S_dispE_EV_toLose_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % model variables compute EV = benefit - cost
            if S_mod_dispE.mdl_cost ~= 0 && S_mod_dispE.mdl_benefit ~= 0
                
                % to gain
                curr_vec_nm = ['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_Benefit_min_Cost_toGain'];
                % positive contrast
                jCon = jCon + 1;
                S_dispE_toGain_mdlBenef_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit_toGain',posCon_nm]);
                S_dispE_toGain_mdlCost_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost_toGain',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_dispE_toGain_mdlBenef_idx,:) - con_vec(S_dispE_toGain_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
                
                % to lose
                curr_vec_nm = ['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_Benefit_min_Cost_toLose'];
                % positive contrast
                jCon = jCon + 1;
                S_dispE_toLose_mdlBenef_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit_toLose',posCon_nm]);
                S_dispE_toLose_mdlCost_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost_toLose',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_dispE_toLose_mdlBenef_idx,:) - con_vec(S_dispE_toLose_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
                
                % pool gain and loss
                curr_vec_nm = ['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_Benefit_min_Cost'];
                % positive contrast
                jCon = jCon + 1;
                S_dispE_mdlBenef_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit',posCon_nm]);
                S_dispE_mdlCost_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_dispE_mdlBenef_idx,:) - con_vec(S_dispE_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
        end
        
        % RT first pair
        if S_mod_dispE.RT_fp ~= 0
            curr_vec_nm = 'S_mod_dispE_RT_fp';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_RTfp_toGain_idx = strcmp(con_names,['S_mod_dispE_RT_fp_toGain',posCon_nm]);
            S_dispE_RTfp_toLose_idx = strcmp(con_names,['S_mod_dispE_RT_fp_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_RTfp_toGain_idx,:) + con_vec(S_dispE_RTfp_toLose_idx,:); % pool RTfp to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % mean RT across pairs (excluding first pair)
        if S_mod_dispE.RT_mRT ~= 0
            curr_vec_nm = 'S_mod_dispE_RT_mRT';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_mRT_toGain_idx = strcmp(con_names,['S_mod_dispE_RT_mRT_toGain',posCon_nm]);
            S_dispE_mRT_toLose_idx = strcmp(con_names,['S_mod_dispE_RT_mRT_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_mRT_toGain_idx,:) + con_vec(S_dispE_mRT_toLose_idx,:); % pool mRT to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % trial number
        if S_mod_dispE.trialN ~= 0
            curr_vec_nm = 'S_mod_dispE_trialN';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_trialN_toGain_idx = strcmp(con_names,['S_mod_dispE_trialN_toGain',posCon_nm]);
            S_dispE_trialN_toLose_idx = strcmp(con_names,['S_mod_dispE_trialN_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_trialN_toGain_idx,:) + con_vec(S_dispE_trialN_toLose_idx,:); % pool RTfp to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ressource
        if S_mod_dispE.X_pred ~= 0
            curr_vec_nm = 'S_mod_dispE_X_pred';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_X_pred_toGain_idx = strcmp(con_names,['S_mod_dispE_X_pred_toGain',posCon_nm]);
            S_dispE_X_pred_toLose_idx = strcmp(con_names,['S_mod_dispE_X_pred_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_X_pred_toGain_idx,:) + con_vec(S_dispE_X_pred_toLose_idx,:); % pool trialN to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % performance
        if S_mod_dispE.perf ~= 0
            curr_vec_nm = 'S_mod_dispE_perf';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_perf_toGain_idx = strcmp(con_names,['S_mod_dispE_perf_toGain',posCon_nm]);
            S_dispE_perf_toLose_idx = strcmp(con_names,['S_mod_dispE_perf_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_perf_toGain_idx,:) + con_vec(S_dispE_perf_toLose_idx,:); % pool trialN to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
    case 3 % trials split according to trial type => contrast across trial types (pool no-error and error trials)
        % no-error trials only
        % one contrast for run 1 and one other for run 2
        S_oDispE_noErrorTrials_idx = strcmp(con_names,['S_o_dispE_noErrorTrials',posCon_nm]);
        S_oDispE_noErrorTrials_vec = con_vec(S_oDispE_noErrorTrials_idx,:);
        S_oDispE_noErrorTrials_r1_idx = find(S_oDispE_noErrorTrials_vec ~= 0, 1, 'first');
        S_oDispE_noErrorTrials_r2_idx = find(S_oDispE_noErrorTrials_vec ~= 0, 1, 'last');
        % run 1
        curr_vec_nm = 'S_o_dispE_noErrorTrials_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_noErrorTrials_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'S_o_dispE_noErrorTrials_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_noErrorTrials_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run 1 - no Error trials
        curr_vec_nm = 'S_o_dispE_noErrorTrials_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_noErrorTrials_r2_idx) = 1;
        con_vec(jCon, S_oDispE_noErrorTrials_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % error trials only
        % one contrast for run 1 and one other for run 2
        S_oDispE_ErrorTrials_idx = strcmp(con_names,['S_o_dispE_ErrorTrials',posCon_nm]);
        S_oDispE_ErrorTrials_vec = con_vec(S_oDispE_ErrorTrials_idx,:);
        S_oDispE_ErrorTrials_r1_idx = find(S_oDispE_ErrorTrials_vec ~= 0, 1, 'first');
        S_oDispE_ErrorTrials_r2_idx = find(S_oDispE_ErrorTrials_vec ~= 0, 1, 'last');
        % run 1
        curr_vec_nm = 'S_o_dispE_ErrorTrials_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_ErrorTrials_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'S_o_dispE_ErrorTrials_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_ErrorTrials_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run 1 - error trials
        curr_vec_nm = 'S_o_dispE_ErrorTrials_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_ErrorTrials_r2_idx) = 1;
        con_vec(jCon, S_oDispE_ErrorTrials_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % pool both trial types
        % run 1
        curr_vec_nm = 'S_o_dispE_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_noErrorTrials_r1_idx) = 1/2;
        con_vec(jCon, S_oDispE_ErrorTrials_r1_idx) = 1/2;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'S_o_dispE_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_noErrorTrials_r2_idx) = 1/2;
        con_vec(jCon, S_oDispE_ErrorTrials_r2_idx) = 1/2;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run 1 - both trial types
        curr_vec_nm = 'S_o_dispE_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, S_oDispE_noErrorTrials_r2_idx)  = 1;
        con_vec(jCon, S_oDispE_ErrorTrials_r2_idx)    = 1;
        con_vec(jCon, S_oDispE_noErrorTrials_r1_idx)  = -1;
        con_vec(jCon, S_oDispE_ErrorTrials_r1_idx)    = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % GL condition
        if S_mod_dispE.inc ~= 0
            curr_vec_nm = 'S_mod_dispE_GL_cond_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_GL_cond_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_GL_cond_noErrorTrials',posCon_nm]);
            S_dispE_GL_cond_ErrorTrials_idx = strcmp(con_names,['S_mod_dispE_GL_cond_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_GL_cond_noErrorTrials_idx,:) + con_vec(S_dispE_GL_cond_ErrorTrials_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % reward type
        if S_mod_dispE.R_type ~= 0
            curr_vec_nm = 'S_mod_dispE_R_type_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_R_type_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_R_type_noErrorTrials',posCon_nm]);
            S_dispE_R_type_ErrorTrials_idx = strcmp(con_names,['S_mod_dispE_R_type_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_R_type_noErrorTrials_idx,:) + con_vec(S_dispE_R_type_ErrorTrials_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive
        if S_mod_dispE.inc ~= 0
            curr_vec_nm = 'S_mod_dispE_inc_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_inc_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_inc_noErrorTrials',posCon_nm]);
            S_dispE_inc_ErrorTrials_idx = strcmp(con_names,['S_mod_dispE_inc_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_inc_noErrorTrials_idx,:) + con_vec(S_dispE_inc_ErrorTrials_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive bis
        if S_mod_dispE.inc_bis ~= 0
            curr_vec_nm = 'S_mod_dispE_inc_bis_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_inc_bis_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_inc_bis_noErrorTrials',posCon_nm]);
            S_dispE_inc_bis_ErrorTrials_idx = strcmp(con_names,['S_mod_dispE_inc_bis_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_inc_bis_noErrorTrials_idx,:) + con_vec(S_dispE_inc_bis_ErrorTrials_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ROI activity
        if S_mod_dispE.ROI_activity_yn ~= 0
            ROI_full_nm = ['S_mod_dispE_GLM',S_mod_dispE.ROI_activity_GLM,...
                    '_',S_mod_dispE.ROI_activity_period,'_period_',S_mod_dispE.ROI_activity_ROI_nm,'_activity'];
            curr_vec_nm = [ROI_full_nm,'_ETnoETpool'];
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_ROIactivity_noErrorTrials_idx = strcmp(con_names,[ROI_full_nm,'_noErrorTrials',posCon_nm]);
            S_dispE_ROIactivity_ErrorTrials_idx = strcmp(con_names,[ROI_full_nm,'_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_ROIactivity_noErrorTrials_idx,:) + con_vec(S_dispE_ROIactivity_ErrorTrials_idx,:); % pool effort to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % number of pairs solved
        if S_mod_dispE.n_pairs_solved ~= 0
            curr_vec_nm = 'S_mod_dispE_n_pairs';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_nPairs_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_n_pairs_solved_noErrorTrials',posCon_nm]);
            S_dispE_nPairs_ErrorTrials_idx = strcmp(con_names,['S_mod_dispE_n_pairs_solved_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_nPairs_noErrorTrials_idx,:) + con_vec(S_dispE_nPairs_ErrorTrials_idx,:); % pool nb pairs to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incongruent pairs solved
        if S_mod_dispE.n_incong ~= 0
            curr_vec_nm = 'S_mod_dispE_n_incong';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_nIncong_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_n_incong_noErrorTrials',posCon_nm]);
            S_dispE_nIncong_ErrorTrials_idx = strcmp(con_names,['S_mod_dispE_n_incong_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_nIncong_noErrorTrials_idx,:) + con_vec(S_dispE_nIncong_ErrorTrials_idx,:); % pool incongruent pairs to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % distance between numbers of pairs solved
        if S_mod_dispE.n_dist ~= 0
            curr_vec_nm = 'S_mod_dispE_n_dist';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_nDist_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_n_dist_noErrorTrials',posCon_nm]);
            S_dispE_nDist_ErrorTrials_idx = strcmp(con_names,['S_mod_dispE_n_dist_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_nDist_noErrorTrials_idx,:) + con_vec(S_dispE_nDist_ErrorTrials_idx,:); % pool dNbers to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % model variables compute EV = benefit - cost
        % model variables
        if S_mod_dispE.mdl_n ~= 0
            
            % pool effort to gain and to lose trials
            % cost
            if S_mod_dispE.mdl_cost ~= 0
                curr_vec_nm = ['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost_ETnoETpool'];
                
                % positive contrast
                jCon = jCon + 1;
                S_dispE_cost_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost__noErrorTrials',posCon_nm]);
                S_dispE_cost_ErrorTrials_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost_ErrorTrials',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_dispE_cost_noErrorTrials_idx,:) + con_vec(S_dispE_cost_ErrorTrials_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % benefit
            if S_mod_dispE.mdl_benefit ~= 0
                curr_vec_nm = ['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit_ETnoETpool'];
                
                % positive contrast
                jCon = jCon + 1;
                S_dispE_benefit_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit__noErrorTrials',posCon_nm]);
                S_dispE_benefit_ErrorTrials_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit_ErrorTrials',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_dispE_benefit_noErrorTrials_idx,:) + con_vec(S_dispE_benefit_ErrorTrials_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % EV
            if S_mod_dispE.mdl_EV ~= 0
                curr_vec_nm = ['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_EV_ETnoETpool'];
                
                % positive contrast
                jCon = jCon + 1;
                S_dispE_EV_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_EV__noErrorTrials',posCon_nm]);
                S_dispE_EV_ErrorTrials_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_EV_ErrorTrials',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_dispE_EV_noErrorTrials_idx,:) + con_vec(S_dispE_EV_ErrorTrials_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % model variables compute EV = benefit - cost
            if S_mod_dispE.mdl_cost ~= 0 && S_mod_dispE.mdl_benefit ~= 0
                
                % no error trials
                curr_vec_nm = ['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_Benefit_min_Cost_noErrorTrials'];
                % positive contrast
                jCon = jCon + 1;
                S_dispE_noErrorTrials_mdlBenef_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit_noErrorTrials',posCon_nm]);
                S_dispE_noErrorTrials_mdlCost_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost_noErrorTrials',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_dispE_noErrorTrials_mdlBenef_idx,:) - con_vec(S_dispE_noErrorTrials_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
                
                % error trials
                curr_vec_nm = ['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_Benefit_min_Cost_ErrorTrials'];
                % positive contrast
                jCon = jCon + 1;
                S_dispE_ErrorTrials_mdlBenef_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit_ErrorTrials',posCon_nm]);
                S_dispE_ErrorTrials_mdlCost_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost_ErrorTrials',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_dispE_ErrorTrials_mdlBenef_idx,:) - con_vec(S_dispE_ErrorTrials_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
                
                % pool gains and losses
                curr_vec_nm = ['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_Benefit_min_Cost_ETnoETpool'];
                % positive contrast
                jCon = jCon + 1;
                S_dispE_mdlBenef_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit',posCon_nm]);
                S_dispE_mdlCost_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(S_dispE_mdlBenef_idx,:) - con_vec(S_dispE_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
        end
        
        % RT first pair
        if S_mod_dispE.RT_fp ~= 0
            curr_vec_nm = 'S_mod_dispE_RT_fp_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_RTfp_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_RT_fp_noErrorTrials',posCon_nm]);
            S_dispE_RTfp_ErrorTrials_idx = strcmp(con_names,['S_mod_dispE_RT_fp_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_RTfp_noErrorTrials_idx,:) + con_vec(S_dispE_RTfp_ErrorTrials_idx,:); % pool RTfp to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % mean RT across pairs (excluding first pair)
        if S_mod_dispE.RT_mRT ~= 0
            curr_vec_nm = 'S_mod_dispE_RT_mRT';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_mRT_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_RT_mRT_noErrorTrials',posCon_nm]);
            S_dispE_mRT_ErrorTrials_idx = strcmp(con_names,['S_mod_dispE_RT_mRT_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_mRT_noErrorTrials_idx,:) + con_vec(S_dispE_mRT_ErrorTrials_idx,:); % pool mRT to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % trial number
        if S_mod_dispE.trialN ~= 0
            curr_vec_nm = 'S_mod_dispE_trialN_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_RTfp_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_RT_fp_noErrorTrials',posCon_nm]);
            S_dispE_RTfp_ErrorTrials_idx = strcmp(con_names,['S_mod_dispE_RT_fp_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_RTfp_noErrorTrials_idx,:) + con_vec(S_dispE_RTfp_ErrorTrials_idx,:); % pool RTfp to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ressource
        if S_mod_dispE.X_pred ~= 0
            curr_vec_nm = 'S_mod_dispE_X_pred_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_X_pred_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_X_pred_noErrorTrials',posCon_nm]);
            S_dispE_X_pred_ErrorTrials_idx = strcmp(con_names,['S_mod_dispE_X_pred_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_X_pred_noErrorTrials_idx,:) + con_vec(S_dispE_X_pred_ErrorTrials_idx,:); % pool trialN to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % performance
        if S_mod_dispE.perf ~= 0
            curr_vec_nm = 'S_mod_dispE_perf_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_dispE_perf_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_perf_noErrorTrials',posCon_nm]);
            S_dispE_perf_ErrorTrials_idx = strcmp(con_names,['S_mod_dispE_perf_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_dispE_perf_noErrorTrials_idx,:) + con_vec(S_dispE_perf_ErrorTrials_idx,:); % pool trialN to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
end


%% effort performance
S_o_perfE = stroopRPprm.o_perfE;
S_mod_perfE = stroopRPprm.mod_perfE;
switch S_o_perfE
    case 2 % trials split according to trial type => contrast across trial types (pool to gain + to lose)
        
        % reward type
        if S_mod_perfE.R_type ~= 0
            curr_vec_nm = 'S_mod_perfE_R_type';
            
            % positive contrast
            jCon = jCon + 1;
            S_perfE_R_type_toGain_idx = strcmp(con_names,['S_mod_perfE_R_type_toGain',posCon_nm]);
            S_perfE_R_type_toLose_idx = strcmp(con_names,['S_mod_perfE_R_type_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_perfE_R_type_toGain_idx,:) + con_vec(S_perfE_R_type_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive
        if S_mod_perfE.inc ~= 0
            curr_vec_nm = 'S_mod_perfE_inc';
            
            % positive contrast
            jCon = jCon + 1;
            S_perfE_inc_toGain_idx = strcmp(con_names,['S_mod_perfE_inc_toGain',posCon_nm]);
            S_perfE_inc_toLose_idx = strcmp(con_names,['S_mod_perfE_inc_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_perfE_inc_toGain_idx,:) + con_vec(S_perfE_inc_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive bis
        if S_mod_perfE.inc_bis ~= 0
            curr_vec_nm = 'S_mod_perfE_inc_bis';
            
            % positive contrast
            jCon = jCon + 1;
            S_perfE_inc_bis_toGain_idx = strcmp(con_names,['S_mod_perfE_inc_bis_toGain',posCon_nm]);
            S_perfE_inc_bis_toLose_idx = strcmp(con_names,['S_mod_perfE_inc_bis_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_perfE_inc_bis_toGain_idx,:) + con_vec(S_perfE_inc_bis_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % number of pairs solved
        if S_mod_perfE.n_pairs_solved ~= 0
            curr_vec_nm = 'S_mod_perfE_n_pairs';
            
            % positive contrast
            jCon = jCon + 1;
            S_perfE_nPairs_toGain_idx = strcmp(con_names,['S_mod_perfE_n_pairs_solved_toGain',posCon_nm]);
            S_perfE_nPairs_toLose_idx = strcmp(con_names,['S_mod_perfE_n_pairs_solved_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_perfE_nPairs_toGain_idx,:) + con_vec(S_perfE_nPairs_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incongruent pairs solved
        if S_mod_perfE.n_incong ~= 0
            curr_vec_nm = 'S_mod_perfE_n_incong';
            
            % positive contrast
            jCon = jCon + 1;
            S_perfE_nIncong_toGain_idx = strcmp(con_names,['S_mod_perfE_n_incong_toGain',posCon_nm]);
            S_perfE_nIncong_toLose_idx = strcmp(con_names,['S_mod_perfE_n_incong_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_perfE_nIncong_toGain_idx,:) + con_vec(S_perfE_nIncong_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % distance between numbers of pairs solved
        if S_mod_perfE.n_dist ~= 0
            curr_vec_nm = 'S_mod_perfE_n_dist';
            
            % positive contrast
            jCon = jCon + 1;
            S_perfE_nDist_toGain_idx = strcmp(con_names,['S_mod_perfE_n_dist_toGain',posCon_nm]);
            S_perfE_nDist_toLose_idx = strcmp(con_names,['S_mod_perfE_n_dist_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_perfE_nDist_toGain_idx,:) + con_vec(S_perfE_nDist_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % errors
        if S_mod_perfE.n_errors ~= 0
            curr_vec_nm = 'S_mod_perfE_n_errors';
            
            % positive contrast
            jCon = jCon + 1;
            S_perfE_nErrors_toGain_idx = strcmp(con_names,['S_mod_perfE_n_errors_toGain',posCon_nm]);
            S_perfE_nErrors_toLose_idx = strcmp(con_names,['S_mod_perfE_n_errors_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_perfE_nErrors_toGain_idx,:) + con_vec(S_perfE_nErrors_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % trial number
        if S_mod_perfE.trialN ~= 0
            curr_vec_nm = 'S_mod_perfE_trialN';
            
            % positive contrast
            jCon = jCon + 1;
            S_perfE_trialN_toGain_idx = strcmp(con_names,['S_mod_perfE_trialN_toGain',posCon_nm]);
            S_perfE_trialN_toLose_idx = strcmp(con_names,['S_mod_perfE_trialN_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_perfE_trialN_toGain_idx,:) + con_vec(S_perfE_trialN_toLose_idx,:); % pool trialN to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
    case 3 % trials split according to trial type => contrast across trial types (pool no-error and error trials)
        % GL condition
        if S_mod_perfE.GL_cond ~= 0
            curr_vec_nm = 'S_mod_perfE_GL_cond';
            
            % positive contrast
            jCon = jCon + 1;
            S_perfE_GL_cond_noErrorTrials_idx = strcmp(con_names,['S_mod_perfE_GL_cond_noErrorTrials',posCon_nm]);
            S_perfE_GL_cond_ErrorTrials_idx = strcmp(con_names,['S_mod_perfE_GL_cond_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_perfE_GL_cond_noErrorTrials_idx,:) + con_vec(S_perfE_GL_cond_ErrorTrials_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % reward type
        if S_mod_perfE.R_type ~= 0
            curr_vec_nm = 'S_mod_perfE_R_type';
            
            % positive contrast
            jCon = jCon + 1;
            S_perfE_R_type_noErrorTrials_idx = strcmp(con_names,['S_mod_perfE_R_type_noErrorTrials',posCon_nm]);
            S_perfE_R_type_ErrorTrials_idx = strcmp(con_names,['S_mod_perfE_R_type_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_perfE_R_type_noErrorTrials_idx,:) + con_vec(S_perfE_R_type_ErrorTrials_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive
        if S_mod_perfE.inc ~= 0
            curr_vec_nm = 'S_mod_perfE_inc';
            
            % positive contrast
            jCon = jCon + 1;
            S_perfE_inc_noErrorTrials_idx = strcmp(con_names,['S_mod_perfE_inc_noErrorTrials',posCon_nm]);
            S_perfE_inc_ErrorTrials_idx = strcmp(con_names,['S_mod_perfE_inc_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_perfE_inc_noErrorTrials_idx,:) + con_vec(S_perfE_inc_ErrorTrials_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive bis
        if S_mod_perfE.inc_bis ~= 0
            curr_vec_nm = 'S_mod_perfE_inc_bis';
            
            % positive contrast
            jCon = jCon + 1;
            S_perfE_inc_bis_noErrorTrials_idx = strcmp(con_names,['S_mod_perfE_inc_bis_noErrorTrials',posCon_nm]);
            S_perfE_inc_bis_ErrorTrials_idx = strcmp(con_names,['S_mod_perfE_inc_bis_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_perfE_inc_bis_noErrorTrials_idx,:) + con_vec(S_perfE_inc_bis_ErrorTrials_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % number of pairs solved
        if S_mod_perfE.n_pairs_solved ~= 0
            curr_vec_nm = 'S_mod_perfE_n_pairs';
            
            % positive contrast
            jCon = jCon + 1;
            S_perfE_nPairs_noErrorTrials_idx = strcmp(con_names,['S_mod_perfE_n_pairs_solved_noErrorTrials',posCon_nm]);
            S_perfE_nPairs_ErrorTrials_idx = strcmp(con_names,['S_mod_perfE_n_pairs_solved_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_perfE_nPairs_noErrorTrials_idx,:) + con_vec(S_perfE_nPairs_ErrorTrials_idx,:); % pool nb pairs to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incongruent pairs solved
        if S_mod_perfE.n_incong ~= 0
            curr_vec_nm = 'S_mod_perfE_n_incong';
            
            % positive contrast
            jCon = jCon + 1;
            S_perfE_nIncong_noErrorTrials_idx = strcmp(con_names,['S_mod_perfE_n_incong_noErrorTrials',posCon_nm]);
            S_perfE_nIncong_ErrorTrials_idx = strcmp(con_names,['S_mod_perfE_n_incong_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_perfE_nIncong_noErrorTrials_idx,:) + con_vec(S_perfE_nIncong_ErrorTrials_idx,:); % pool incongruent pairs to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % distance between numbers of pairs solved
        if S_mod_perfE.n_dist ~= 0
            curr_vec_nm = 'S_mod_perfE_n_dist';
            
            % positive contrast
            jCon = jCon + 1;
            S_perfE_nDist_noErrorTrials_idx = strcmp(con_names,['S_mod_perfE_n_dist_noErrorTrials',posCon_nm]);
            S_perfE_nDist_ErrorTrials_idx = strcmp(con_names,['S_mod_perfE_n_dist_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_perfE_nDist_noErrorTrials_idx,:) + con_vec(S_perfE_nDist_ErrorTrials_idx,:); % pool dNbers to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % trial number
        if S_mod_perfE.trialN ~= 0
            curr_vec_nm = 'S_mod_perfE_trialN';
            
            % positive contrast
            jCon = jCon + 1;
            S_perfE_trialN_noErrorTrials_idx = strcmp(con_names,['S_mod_perfE_trialN_noErrorTrials',posCon_nm]);
            S_perfE_trialN_ErrorTrials_idx = strcmp(con_names,['S_mod_perfE_trialN_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_perfE_trialN_noErrorTrials_idx,:) + con_vec(S_perfE_trialN_ErrorTrials_idx,:); % pool dNbers to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
end


%% feedback
S_o_fbk = stroopRPprm.o_fbk;
S_mod_fbk = stroopRPprm.mod_fbk;
switch S_o_fbk
    case 2 % to gain/to lose trials
        
        % feedback
        if S_mod_fbk.fbk ~= 0
            curr_vec_nm = 'S_mod_fbk_fbk';
            
            % positive contrast
            jCon = jCon + 1;
            S_fbk_fbk_toGain_idx = strcmp(con_names,['S_mod_fbk_fbk_toGain',posCon_nm]);
            S_fbk_fbk_toLose_idx = strcmp(con_names,['S_mod_fbk_fbk_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_fbk_fbk_toGain_idx,:) + con_vec(S_fbk_fbk_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % total gain
        if S_mod_fbk.totalGain ~= 0
            curr_vec_nm = 'S_mod_fbk_totalGain';
            
            % positive contrast
            jCon = jCon + 1;
            S_fbk_totalGain_toGain_idx = strcmp(con_names,['S_mod_fbk_totalGain_toGain',posCon_nm]);
            S_fbk_totalGain_toLose_idx = strcmp(con_names,['S_mod_fbk_totalGain_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_fbk_totalGain_toGain_idx,:) + con_vec(S_fbk_totalGain_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ressource
        if S_mod_fbk.X_pred ~= 0
            curr_vec_nm = 'S_mod_fbk_X_pred';
            
            % positive contrast
            jCon = jCon + 1;
            S_fbk_X_pred_toGain_idx = strcmp(con_names,['S_mod_fbk_X_pred_toGain',posCon_nm]);
            S_fbk_X_pred_toLose_idx = strcmp(con_names,['S_mod_fbk_X_pred_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_fbk_X_pred_toGain_idx,:) + con_vec(S_fbk_X_pred_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % performance
        if S_mod_fbk.perf ~= 0
            curr_vec_nm = 'S_mod_fbk_perf';
            
            % positive contrast
            jCon = jCon + 1;
            S_fbk_perf_toGain_idx = strcmp(con_names,['S_mod_fbk_perf_toGain',posCon_nm]);
            S_fbk_perf_toLose_idx = strcmp(con_names,['S_mod_fbk_perf_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_fbk_perf_toGain_idx,:) + con_vec(S_fbk_perf_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ROI activity
        if S_mod_fbk.ROI_activity_yn ~= 0
            ROI_full_nm = ['S_mod_fbk_GLM',S_mod_fbk.ROI_activity_GLM,...
                    '_',S_mod_fbk.ROI_activity_period,'_period_',S_mod_fbk.ROI_activity_ROI_nm,'_activity'];
            curr_vec_nm = [ROI_full_nm,''];
            
            % positive contrast
            jCon = jCon + 1;
            S_fbk_ROIactivity_toGain_idx = strcmp(con_names,[ROI_full_nm,'_toGain',posCon_nm]);
            S_fbk_ROIactivity_toLose_idx = strcmp(con_names,[ROI_full_nm,'_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_fbk_ROIactivity_toGain_idx,:) + con_vec(S_fbk_ROIactivity_toLose_idx,:); % pool effort to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % trial number
        if S_mod_fbk.trialN ~= 0
            curr_vec_nm = 'S_mod_fbk_trialN';
            
            % positive contrast
            jCon = jCon + 1;
            S_fbk_trialN_toGain_idx = strcmp(con_names,['S_mod_fbk_trialN_toGain',posCon_nm]);
            S_fbk_trialN_toLose_idx = strcmp(con_names,['S_mod_fbk_trialN_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_fbk_trialN_toGain_idx,:) + con_vec(S_fbk_trialN_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
    case 3 % no-error/error trials
        % GL condition
        if S_mod_fbk.GL_cond ~= 0
            curr_vec_nm = 'S_mod_fbk_GL_cond_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_fbk_GL_cond_noErrorTrials_idx = strcmp(con_names,['S_mod_fbk_GL_cond_noErrorTrials',posCon_nm]);
            S_fbk_GL_cond_ErrorTrials_idx = strcmp(con_names,['S_mod_fbk_GL_cond_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_fbk_GL_cond_noErrorTrials_idx,:) + con_vec(S_fbk_GL_cond_ErrorTrials_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % feedback
        if S_mod_fbk.fbk ~= 0
            curr_vec_nm = 'S_mod_fbk_fbk_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_fbk_fbk_noErrorTrials_idx = strcmp(con_names,['S_mod_fbk_fbk_noErrorTrials',posCon_nm]);
            S_fbk_fbk_ErrorTrials_idx = strcmp(con_names,['S_mod_fbk_fbk_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_fbk_fbk_noErrorTrials_idx,:) + con_vec(S_fbk_fbk_ErrorTrials_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % total gain
        if S_mod_fbk.totalGain ~= 0
            curr_vec_nm = 'S_mod_fbk_totalGain_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_fbk_totalGain_noErrorTrials_idx = strcmp(con_names,['S_mod_fbk_totalGain_noErrorTrials',posCon_nm]);
            S_fbk_totalGain_ErrorTrials_idx = strcmp(con_names,['S_mod_fbk_totalGain_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_fbk_totalGain_noErrorTrials_idx,:) + con_vec(S_fbk_totalGain_ErrorTrials_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ressource
        if S_mod_fbk.X_pred ~= 0
            curr_vec_nm = 'S_mod_fbk_X_pred_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_fbk_X_pred_noErrorTrials_idx = strcmp(con_names,['S_mod_fbk_X_pred_noErrorTrials',posCon_nm]);
            S_fbk_X_pred_ErrorTrials_idx = strcmp(con_names,['S_mod_fbk_X_pred_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_fbk_X_pred_noErrorTrials_idx,:) + con_vec(S_fbk_X_pred_ErrorTrials_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % performance
        if S_mod_fbk.perf ~= 0
            curr_vec_nm = 'S_mod_fbk_perf_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_fbk_perf_noErrorTrials_idx = strcmp(con_names,['S_mod_fbk_perf_noErrorTrials',posCon_nm]);
            S_fbk_perf_ErrorTrials_idx = strcmp(con_names,['S_mod_fbk_perf_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_fbk_perf_noErrorTrials_idx,:) + con_vec(S_fbk_perf_ErrorTrials_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ROI activity
        if S_mod_fbk.ROI_activity_yn ~= 0
            ROI_full_nm = ['S_mod_fbk_GLM',S_mod_fbk.ROI_activity_GLM,...
                    '_',S_mod_fbk.ROI_activity_period,'_period_',S_mod_fbk.ROI_activity_ROI_nm,'_activity'];
            curr_vec_nm = [ROI_full_nm,'_ETnoETpool'];
            
            % positive contrast
            jCon = jCon + 1;
            S_fbk_ROIactivity_noErrorTrials_idx = strcmp(con_names,[ROI_full_nm,'_noErrorTrials',posCon_nm]);
            S_fbk_ROIactivity_ErrorTrials_idx = strcmp(con_names,[ROI_full_nm,'_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_fbk_ROIactivity_noErrorTrials_idx,:) + con_vec(S_fbk_ROIactivity_ErrorTrials_idx,:); % pool effort to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % trial number
        if S_mod_fbk.trialN ~= 0
            curr_vec_nm = 'S_mod_fbk_trialN_ETnoETpool';
            
            % positive contrast
            jCon = jCon + 1;
            S_fbk_trialN_noErrorTrials_idx = strcmp(con_names,['S_mod_fbk_trialN_noErrorTrials',posCon_nm]);
            S_fbk_trialN_ErrorTrials_idx = strcmp(con_names,['S_mod_fbk_trialN_ErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(S_fbk_trialN_noErrorTrials_idx,:) + con_vec(S_fbk_trialN_ErrorTrials_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
end


end % function