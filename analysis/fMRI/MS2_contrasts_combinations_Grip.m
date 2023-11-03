function [ jCon, con_vec, con_names ] = MS2_contrasts_combinations_Grip( jCon, con_vec, con_names,...
    posCon_nm, n_regs, gripRPprm)
%[ jCon, con_vec, con_names ] = MS2_contrasts_combinations_Grip( jCon, con_vec, con_names,...
%     posCon_nm, n_regs, gripRPprm) combine similar contrasts in the
%     grip task.


%% incentive
G_o_inc = gripRPprm.o_inc;
G_mod_inc = gripRPprm.mod_inc;
switch G_o_inc
    case 1
        % one contrast for run 1 and one other for run 2
        G_oInc_idx = strcmp(con_names,['G_o_inc',posCon_nm]);
        G_oInc_vec = con_vec(G_oInc_idx,:);
        G_oInc_r1_idx = find(G_oInc_vec ~= 0, 1, 'first');
        G_oInc_r2_idx = find(G_oInc_vec ~= 0, 1, 'last');
        
        % run 1
        curr_vec_nm = 'G_o_inc_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oInc_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'G_o_inc_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oInc_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run 1
        curr_vec_nm = 'G_o_inc_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oInc_r2_idx) = 1;
        con_vec(jCon, G_oInc_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % model variables compute EV = benefit - cost
        if G_mod_inc.mdl_n ~= 0 && G_mod_inc.mdl_cost ~= 0 && G_mod_inc.mdl_benefit ~= 0
            curr_vec_nm = ['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_Benefit_min_Cost'];
            
            % positive contrast
            jCon = jCon + 1;
            G_inc_mdlBenef_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_benefit',posCon_nm]);
            G_inc_mdlCost_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_cost',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_mdlBenef_idx,:) - con_vec(G_inc_mdlCost_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
    case 2 % trials split according to trial type => contrast across trial types (pool to gain + to lose)
        
        % to gain only
        % one contrast for run 1 and one other for run 2
        G_oInc_toGain_idx = strcmp(con_names,['G_o_inc_toGain',posCon_nm]);
        G_oInc_toGain_vec = con_vec(G_oInc_toGain_idx,:);
        G_oInc_toGain_r1_idx = find(G_oInc_toGain_vec ~= 0, 1, 'first');
        G_oInc_toGain_r2_idx = find(G_oInc_toGain_vec ~= 0, 1, 'last');
        % run 1
        curr_vec_nm = 'G_o_inc_toGain_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oInc_toGain_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'G_o_inc_toGain_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oInc_toGain_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run1 - toGain
        curr_vec_nm = 'G_o_inc_toGain_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oInc_toGain_r2_idx) = 1;
        con_vec(jCon, G_oInc_toGain_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % to lose only
        % one contrast for run 1 and one other for run 2
        G_oInc_toLose_idx = strcmp(con_names,['G_o_inc_toLose',posCon_nm]);
        G_oInc_toLose_vec = con_vec(G_oInc_toLose_idx,:);
        G_oInc_toLose_r1_idx = find(G_oInc_toLose_vec ~= 0, 1, 'first');
        G_oInc_toLose_r2_idx = find(G_oInc_toLose_vec ~= 0, 1, 'last');
        
        % run 1
        curr_vec_nm = 'G_o_inc_toLose_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oInc_toLose_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'G_o_inc_toLose_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oInc_toLose_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run 1 - toLose
        curr_vec_nm = 'G_o_inc_toLose_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oInc_toLose_r2_idx) = 1;
        con_vec(jCon, G_oInc_toLose_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % pool both trial types
        % run 1
        curr_vec_nm = 'G_o_inc_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oInc_toGain_r1_idx) = 1;
        con_vec(jCon, G_oInc_toLose_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'G_o_inc_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oInc_toGain_r2_idx) = 1/2;
        con_vec(jCon, G_oInc_toLose_r2_idx) = 1/2;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 min run 1
        curr_vec_nm = 'G_o_inc_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oInc_toGain_r2_idx) = 1;
        con_vec(jCon, G_oInc_toLose_r2_idx) = 1;
        con_vec(jCon, G_oInc_toGain_r1_idx) = -1;
        con_vec(jCon, G_oInc_toLose_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % compare incentive Gain vs incentive Loss period
        % pool gain and loss condition for incentive modulation
        curr_vec_nm = 'G_o_inc_Gain_min_Loss';
        
        % positive contrast
        jCon = jCon + 1;
        G_Oinc_toGain_idx = strcmp(con_names,['G_o_inc_toGain',posCon_nm]);
        G_Oinc_toLose_idx = strcmp(con_names,['G_o_inc_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_Oinc_toGain_idx,:) - con_vec(G_Oinc_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % reward type
        if G_mod_inc.R_type ~= 0
            % pool gain and loss condition for incentive modulation
            curr_vec_nm = 'G_mod_inc_R_type';
            
            % positive contrast
            jCon = jCon + 1;
            G_R_type_toGain_idx = strcmp(con_names,['G_mod_inc_R_type_toGain',posCon_nm]);
            G_R_type_toLose_idx = strcmp(con_names,['G_mod_inc_R_type_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_R_type_toGain_idx,:) + con_vec(G_R_type_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % comparison gain and loss condition for incentive modulation
            curr_vec_nm = 'G_mod_inc_R_type_Gain_min_Loss';
            
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(G_R_type_toGain_idx,:) - con_vec(G_R_type_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive
        if G_mod_inc.inc ~= 0
            % pool gain and loss condition for incentive modulation
            curr_vec_nm = 'G_mod_inc_inc';
            
            % positive contrast
            jCon = jCon + 1;
            G_inc_toGain_idx = strcmp(con_names,['G_mod_inc_inc_toGain',posCon_nm]);
            G_inc_toLose_idx = strcmp(con_names,['G_mod_inc_inc_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_toGain_idx,:) + con_vec(G_inc_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % comparison gain and loss condition for incentive modulation
            curr_vec_nm = 'G_mod_inc_inc_Gain_min_Loss';
            
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_toGain_idx,:) - con_vec(G_inc_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive bis
        if G_mod_inc.inc_bis ~= 0
            curr_vec_nm = 'G_mod_inc_inc_bis';
            
            % positive contrast
            jCon = jCon + 1;
            G_inc_bis_toGain_idx = strcmp(con_names,['G_mod_inc_inc_bis_toGain',posCon_nm]);
            G_inc_bis_toLose_idx = strcmp(con_names,['G_mod_inc_inc_bis_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_bis_toGain_idx,:) + con_vec(G_inc_bis_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % comparison gain and loss condition for incentive modulation
            curr_vec_nm = 'G_mod_inc_inc_bis_Gain_min_Loss';
            
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_bis_toGain_idx,:) - con_vec(G_inc_bis_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ROI activity
        if G_mod_inc.ROI_activity_yn ~= 0
            ROI_full_nm = ['G_mod_inc_GLM',G_mod_inc.ROI_activity_GLM,...
                    '_',G_mod_inc.ROI_activity_period,'_period_',G_mod_inc.ROI_activity_ROI_nm,'_activity'];
            curr_vec_nm = ROI_full_nm;
            
            % positive contrast
            jCon = jCon + 1;
            G_inc_ROIactivity_toGain_idx = strcmp(con_names,[ROI_full_nm,'_toGain',posCon_nm]);
            G_inc_ROIactivity_toLose_idx = strcmp(con_names,[ROI_full_nm,'_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_ROIactivity_toGain_idx,:) + con_vec(G_inc_ROIactivity_toLose_idx,:); % pool effort to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % effort
        if G_mod_inc.Effort ~= 0
            curr_vec_nm = 'G_mod_inc_effort';
            
            % positive contrast
            jCon = jCon + 1;
            G_inc_eff_toGain_idx = strcmp(con_names,['G_mod_inc_effort_toGain',posCon_nm]);
            G_inc_eff_toLose_idx = strcmp(con_names,['G_mod_inc_effort_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_eff_toGain_idx,:) + con_vec(G_inc_eff_toLose_idx,:); % pool effort to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % model variables
        if G_mod_inc.mdl_n ~= 0
            
            % pool effort to gain and to lose trials
            % cost
            if G_mod_inc.mdl_cost ~= 0
                curr_vec_nm = ['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_cost'];
                
                % positive contrast
                jCon = jCon + 1;
                G_inc_cost_toGain_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_cost_toGain',posCon_nm]);
                G_inc_cost_toLose_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_cost_toLose',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(G_inc_cost_toGain_idx,:) + con_vec(G_inc_cost_toLose_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % benefit
            if G_mod_inc.mdl_benefit ~= 0
                curr_vec_nm = ['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_benefit'];
                
                % positive contrast
                jCon = jCon + 1;
                G_inc_benefit_toGain_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_benefit_toGain',posCon_nm]);
                G_inc_benefit_toLose_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_benefit_toLose',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(G_inc_benefit_toGain_idx,:) + con_vec(G_inc_benefit_toLose_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % EV
            if G_mod_inc.mdl_EV ~= 0
                curr_vec_nm = ['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_EV'];
                
                % positive contrast
                jCon = jCon + 1;
                G_inc_EV_toGain_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_EV_toGain',posCon_nm]);
                G_inc_EV_toLose_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_EV_toLose',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(G_inc_EV_toGain_idx,:) + con_vec(G_inc_EV_toLose_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % model variables compute EV = benefit - cost
            if G_mod_inc.mdl_cost ~= 0 && G_mod_inc.mdl_benefit ~= 0
                
                % to gain
                curr_vec_nm = ['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_Benefit_min_Cost_toGain'];
                % positive contrast
                jCon = jCon + 1;
                G_inc_toGain_mdlBenef_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_benefit_toGain',posCon_nm]);
                G_inc_toGain_mdlCost_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_cost_toGain',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(G_inc_toGain_mdlBenef_idx,:) - con_vec(G_inc_toGain_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
                
                
                % to Lose
                curr_vec_nm = ['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_Benefit_min_Cost_toLose'];
                % positive contrast
                jCon = jCon + 1;
                G_inc_toLose_mdlBenef_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_benefit_toLose',posCon_nm]);
                G_inc_toLose_mdlCost_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_cost_toLose',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(G_inc_toLose_mdlBenef_idx,:) - con_vec(G_inc_toLose_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
                
                
                % pool both
                curr_vec_nm = ['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_Benefit_min_Cost'];
                % positive contrast
                jCon = jCon + 1;
                G_inc_mdlBenef_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_benefit',posCon_nm]);
                G_inc_mdlCost_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_cost',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(G_inc_mdlBenef_idx,:) - con_vec(G_inc_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
        end
        
        % RT first squeeze
        if G_mod_inc.RT_fp ~= 0
            curr_vec_nm = 'G_mod_inc_RT_fp';
            
            % positive contrast
            jCon = jCon + 1;
            G_inc_RTfp_toGain_idx = strcmp(con_names,['G_mod_inc_RT_fp_toGain',posCon_nm]);
            G_inc_RTfp_toLose_idx = strcmp(con_names,['G_mod_inc_RT_fp_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_RTfp_toGain_idx,:) + con_vec(G_inc_RTfp_toLose_idx,:); % pool RT to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % total gain until now
        if G_mod_inc.totalGain_prev ~= 0
            curr_vec_nm = 'G_mod_inc_totalGain_prev';
            
            % positive contrast
            jCon = jCon + 1;
            G_inc_totalGainPrev_toGain_idx = strcmp(con_names,['G_mod_inc_totalGain_prev_toGain',posCon_nm]);
            G_inc_totalGainPrev_toLose_idx = strcmp(con_names,['G_mod_inc_totalGain_prev_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_totalGainPrev_toGain_idx,:) + con_vec(G_inc_totalGainPrev_toLose_idx,:); % pool to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % sum(performance previous trials)
         if G_mod_inc.sumPerfPrev ~= 0
            curr_vec_nm = 'G_mod_inc_sumPerfPrev';
            
            % positive contrast
            jCon = jCon + 1;
            G_inc_sumPerfPrev_toGain_idx = strcmp(con_names,['G_mod_inc_sumPerfPrev_toGain',posCon_nm]);
            G_inc_sumPerfPrev_toLose_idx = strcmp(con_names,['G_mod_inc_sumPerfPrev_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_sumPerfPrev_toGain_idx,:) + con_vec(G_inc_sumPerfPrev_toLose_idx,:); % pool to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % trial number
        if G_mod_inc.trialN ~= 0
            curr_vec_nm = 'G_mod_inc_trialN';
            
            % positive contrast
            jCon = jCon + 1;
            G_inc_trialN_toGain_idx = strcmp(con_names,['G_mod_inc_trialN_toGain',posCon_nm]);
            G_inc_trialN_toLose_idx = strcmp(con_names,['G_mod_inc_trialN_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_trialN_toGain_idx,:) + con_vec(G_inc_trialN_toLose_idx,:); % pool trialN to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ressource
        if G_mod_inc.X_pred ~= 0
            curr_vec_nm = 'G_mod_inc_X_pred';
            
            % positive contrast
            jCon = jCon + 1;
            G_inc_X_pred_toGain_idx = strcmp(con_names,['G_mod_inc_X_pred_toGain',posCon_nm]);
            G_inc_X_pred_toLose_idx = strcmp(con_names,['G_mod_inc_X_pred_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_X_pred_toGain_idx,:) + con_vec(G_inc_X_pred_toLose_idx,:); % pool to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % performance
        if G_mod_inc.perf ~= 0
            curr_vec_nm = 'G_mod_inc_perf';
            
            % positive contrast
            jCon = jCon + 1;
            G_inc_perf_toGain_idx = strcmp(con_names,['G_mod_inc_perf_toGain',posCon_nm]);
            G_inc_perf_toLose_idx = strcmp(con_names,['G_mod_inc_perf_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_perf_toGain_idx,:) + con_vec(G_inc_perf_toLose_idx,:); % pool trialN to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
end


%% effort scale
G_o_dispE = gripRPprm.o_dispE;
G_mod_dispE = gripRPprm.mod_dispE;
switch G_o_dispE
    case 1
        % one contrast for run 1 and one other for run 2
        G_oDispE_idx = strcmp(con_names,['G_o_dispE',posCon_nm]);
        G_oDispE_vec = con_vec(G_oDispE_idx,:);
        G_oDispE_r1_idx = find(G_oDispE_vec ~= 0, 1, 'first');
        G_oDispE_r2_idx = find(G_oDispE_vec ~= 0, 1, 'last');
        
        % run 1
        curr_vec_nm = 'G_o_dispE_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oDispE_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'G_o_dispE_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oDispE_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run 1
        curr_vec_nm = 'G_o_dispE_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oDispE_r2_idx) = 1;
        con_vec(jCon, G_oDispE_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % model variables compute EV = benefit - cost
        if G_mod_dispE.mdl_n ~= 0 && G_mod_dispE.mdl_cost ~= 0 && G_mod_dispE.mdl_benefit ~= 0
            curr_vec_nm = ['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_Benefit_min_Cost'];
            
            % positive contrast
            jCon = jCon + 1;
            G_dispE_mdlBenef_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_benefit',posCon_nm]);
            G_dispE_mdlCost_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_cost',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_mdlBenef_idx,:) - con_vec(G_dispE_mdlCost_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
    case 2 % trials split according to trial type => contrast across trial types (pool to gain + to lose)
        
        % to gain only
        % one contrast for run 1 and one other for run 2
        G_oDispE_toGain_idx = strcmp(con_names,['G_o_dispE_toGain',posCon_nm]);
        G_oDispE_toGain_vec = con_vec(G_oDispE_toGain_idx,:);
        G_oDispE_toGain_r1_idx = find(G_oDispE_toGain_vec ~= 0, 1, 'first');
        G_oDispE_toGain_r2_idx = find(G_oDispE_toGain_vec ~= 0, 1, 'last');
        % run 1
        curr_vec_nm = 'G_o_dispE_toGain_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oDispE_toGain_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'G_o_dispE_toGain_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oDispE_toGain_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run1 - toGain
        curr_vec_nm = 'G_o_dispE_toGain_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oDispE_toGain_r2_idx) = 1;
        con_vec(jCon, G_oDispE_toGain_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % to lose only
        % one contrast for run 1 and one other for run 2
        G_oDispE_toLose_idx = strcmp(con_names,['G_o_dispE_toLose',posCon_nm]);
        G_oDispE_toLose_vec = con_vec(G_oDispE_toLose_idx,:);
        G_oDispE_toLose_r1_idx = find(G_oDispE_toLose_vec ~= 0, 1, 'first');
        G_oDispE_toLose_r2_idx = find(G_oDispE_toLose_vec ~= 0, 1, 'last');
        
        % run 1
        curr_vec_nm = 'G_o_dispE_toLose_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oDispE_toLose_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'G_o_dispE_toLose_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oDispE_toLose_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 - run 1 - toLose
        curr_vec_nm = 'G_o_dispE_toLose_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oDispE_toLose_r2_idx) = 1;
        con_vec(jCon, G_oDispE_toLose_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % pool both trial types
        % run 1
        curr_vec_nm = 'G_o_dispE_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oDispE_toGain_r1_idx) = 1;
        con_vec(jCon, G_oDispE_toLose_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'G_o_dispE_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oDispE_toGain_r2_idx) = 1/2;
        con_vec(jCon, G_oDispE_toLose_r2_idx) = 1/2;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2 min run 1
        curr_vec_nm = 'G_o_dispE_run2_min_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, G_oDispE_toGain_r2_idx) = 1;
        con_vec(jCon, G_oDispE_toLose_r2_idx) = 1;
        con_vec(jCon, G_oDispE_toGain_r1_idx) = -1;
        con_vec(jCon, G_oDispE_toLose_r1_idx) = -1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % compare incentive Gain vs incentive Loss period
        % pool gain and loss condition for incentive modulation
        curr_vec_nm = 'G_o_dispE_Gain_min_Loss';
        
        % positive contrast
        jCon = jCon + 1;
        G_OdispE_toGain_idx = strcmp(con_names,['G_o_dispE_toGain',posCon_nm]);
        G_OdispE_toLose_idx = strcmp(con_names,['G_o_dispE_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_OdispE_toGain_idx,:) - con_vec(G_OdispE_toLose_idx,:); % dispE to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % reward type
        if G_mod_dispE.R_type ~= 0
            curr_vec_nm = 'G_mod_dispE_R_type';
            
            % positive contrast
            jCon = jCon + 1;
            G_dispE_R_type_toGain_idx = strcmp(con_names,['G_mod_dispE_R_type_toGain',posCon_nm]);
            G_dispE_R_type_toLose_idx = strcmp(con_names,['G_mod_dispE_R_type_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_R_type_toGain_idx,:) + con_vec(G_dispE_R_type_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive
        if G_mod_dispE.inc ~= 0
            curr_vec_nm = 'G_mod_dispE_inc';
            
            % positive contrast
            jCon = jCon + 1;
            G_dispE_inc_toGain_idx = strcmp(con_names,['G_mod_dispE_inc_toGain',posCon_nm]);
            G_dispE_inc_toLose_idx = strcmp(con_names,['G_mod_dispE_inc_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_inc_toGain_idx,:) + con_vec(G_dispE_inc_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive bis
        if G_mod_dispE.inc_bis ~= 0
            curr_vec_nm = 'G_mod_dispE_inc_bis';
            
            % positive contrast
            jCon = jCon + 1;
            G_dispE_inc_bis_toGain_idx = strcmp(con_names,['G_mod_dispE_inc_bis_toGain',posCon_nm]);
            G_dispE_inc_bis_toLose_idx = strcmp(con_names,['G_mod_dispE_inc_bis_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_inc_bis_toGain_idx,:) + con_vec(G_dispE_inc_bis_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ROI activity
        if G_mod_dispE.ROI_activity_yn ~= 0
            ROI_full_nm = ['G_mod_dispE_GLM',G_mod_dispE.ROI_activity_GLM,...
                    '_',G_mod_dispE.ROI_activity_period,'_period_',G_mod_dispE.ROI_activity_ROI_nm,'_activity'];
            curr_vec_nm = ROI_full_nm;
            
            % positive contrast
            jCon = jCon + 1;
            G_dispE_ROIactivity_toGain_idx = strcmp(con_names,[ROI_full_nm,'_toGain',posCon_nm]);
            G_dispE_ROIactivity_toLose_idx = strcmp(con_names,[ROI_full_nm,'_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_ROIactivity_toGain_idx,:) + con_vec(G_dispE_ROIactivity_toLose_idx,:); % pool effort to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % effort exertion
        if G_mod_dispE.Effort ~= 0
            curr_vec_nm = 'G_mod_dispE_Effort';
            
            % positive contrast
            jCon = jCon + 1;
            G_dispE_E_toGain_idx = strcmp(con_names,['G_mod_dispE_Effort_toGain',posCon_nm]);
            G_dispE_E_toLose_idx = strcmp(con_names,['G_mod_dispE_Effort_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_E_toGain_idx,:) + con_vec(G_dispE_E_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % model variables compute EV = benefit - cost
        % model variables
        if G_mod_dispE.mdl_n ~= 0
            
            % pool effort to gain and to lose trials
            % cost
            if G_mod_dispE.mdl_cost ~= 0
                curr_vec_nm = ['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_cost'];
                
                % positive contrast
                jCon = jCon + 1;
                G_dispE_cost_toGain_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_cost_toGain',posCon_nm]);
                G_dispE_cost_toLose_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_cost_toLose',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(G_dispE_cost_toGain_idx,:) + con_vec(G_dispE_cost_toLose_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % benefit
            if G_mod_dispE.mdl_benefit ~= 0
                curr_vec_nm = ['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_benefit'];
                
                % positive contrast
                jCon = jCon + 1;
                G_dispE_benefit_toGain_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_benefit_toGain',posCon_nm]);
                G_dispE_benefit_toLose_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_benefit_toLose',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(G_dispE_benefit_toGain_idx,:) + con_vec(G_dispE_benefit_toLose_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % EV
            if G_mod_dispE.mdl_EV ~= 0
                curr_vec_nm = ['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_EV'];
                
                % positive contrast
                jCon = jCon + 1;
                G_dispE_EV_toGain_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_EV_toGain',posCon_nm]);
                G_dispE_EV_toLose_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_EV_toLose',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(G_dispE_EV_toGain_idx,:) + con_vec(G_dispE_EV_toLose_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
            
            % model variables compute EV = benefit - cost
            if G_mod_dispE.mdl_cost ~= 0 && G_mod_dispE.mdl_benefit ~= 0
                % to gain
                curr_vec_nm = ['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_Benefit_min_Cost_toGain'];
                % positive contrast
                jCon = jCon + 1;
                G_dispE_toGain_mdlBenef_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_benefit_toGain',posCon_nm]);
                G_dispE_toGain_mdlCost_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_cost_toGain',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(G_dispE_toGain_mdlBenef_idx,:) - con_vec(G_dispE_toGain_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
                
                
                % to Lose
                curr_vec_nm = ['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_Benefit_min_Cost_toLose'];
                % positive contrast
                jCon = jCon + 1;
                G_dispE_toLose_mdlBenef_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_benefit_toLose',posCon_nm]);
                G_dispE_toLose_mdlCost_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_cost_toLose',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(G_dispE_toLose_mdlBenef_idx,:) - con_vec(G_dispE_toLose_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
                
                % pool gain and loss
                curr_vec_nm = ['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_Benefit_min_Cost'];
                % positive contrast
                jCon = jCon + 1;
                G_dispE_mdlBenef_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_benefit',posCon_nm]);
                G_dispE_mdlCost_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_cost',posCon_nm]);
                con_vec(jCon, 1:n_regs) = con_vec(G_dispE_mdlBenef_idx,:) - con_vec(G_dispE_mdlCost_idx,:);
                con_names{jCon} = [curr_vec_nm,posCon_nm];
                % negative contrast
                [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            end
        end
        
        % RT first press
        if G_mod_dispE.RT_fp ~= 0
            curr_vec_nm = 'G_mod_dispE_RT_fp';
            
            % positive contrast
            jCon = jCon + 1;
            G_dispE_RTfp_toGain_idx = strcmp(con_names,['G_mod_dispE_RT_fp_toGain',posCon_nm]);
            G_dispE_RTfp_toLose_idx = strcmp(con_names,['G_mod_dispE_RT_fp_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_RTfp_toGain_idx,:) + con_vec(G_dispE_RTfp_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % trial number
        if G_mod_dispE.trialN ~= 0
            curr_vec_nm = 'G_mod_dispE_trialN';
            
            % positive contrast
            jCon = jCon + 1;
            G_dispE_trialN_toGain_idx = strcmp(con_names,['G_mod_dispE_trialN_toGain',posCon_nm]);
            G_dispE_trialN_toLose_idx = strcmp(con_names,['G_mod_dispE_trialN_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_trialN_toGain_idx,:) + con_vec(G_dispE_trialN_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ressource
        if G_mod_dispE.X_pred ~= 0
            curr_vec_nm = 'G_mod_dispE_X_pred';
            
            % positive contrast
            jCon = jCon + 1;
            G_mod_dispE_X_pred_toGain_idx = strcmp(con_names,['G_mod_dispE_X_pred_toGain',posCon_nm]);
            G_mod_dispE_X_pred_toLose_idx = strcmp(con_names,['G_mod_dispE_X_pred_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_mod_dispE_X_pred_toGain_idx,:) + con_vec(G_mod_dispE_X_pred_toLose_idx,:); % pool to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % performance
        if G_mod_dispE.perf ~= 0
            curr_vec_nm = 'G_mod_dispE_perf';
            
            % positive contrast
            jCon = jCon + 1;
            G_dispE_perf_toGain_idx = strcmp(con_names,['G_mod_dispE_perf_toGain',posCon_nm]);
            G_dispE_perf_toLose_idx = strcmp(con_names,['G_mod_dispE_perf_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_perf_toGain_idx,:) + con_vec(G_dispE_perf_toLose_idx,:); % pool trialN to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
end


%% effort performance
G_o_perfE = gripRPprm.o_perfE;
G_mod_perfE = gripRPprm.mod_perfE;
switch G_o_perfE
    case 2 % trials split according to trial type => contrast across trial types (pool to gain + to lose)
        
        % reward type
        if G_mod_perfE.R_type ~= 0
            curr_vec_nm = 'G_mod_perfE_R_type';
            
            % positive contrast
            jCon = jCon + 1;
            G_perfE_R_type_toGain_idx = strcmp(con_names,['G_mod_perfE_R_type_toGain',posCon_nm]);
            G_perfE_R_type_toLose_idx = strcmp(con_names,['G_mod_perfE_R_type_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_perfE_R_type_toGain_idx,:) + con_vec(G_perfE_R_type_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive
        if G_mod_perfE.inc ~= 0
            curr_vec_nm = 'G_mod_perfE_inc';
            
            % positive contrast
            jCon = jCon + 1;
            G_perfE_inc_toGain_idx = strcmp(con_names,['G_mod_perfE_inc_toGain',posCon_nm]);
            G_perfE_inc_toLose_idx = strcmp(con_names,['G_mod_perfE_inc_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_perfE_inc_toGain_idx,:) + con_vec(G_perfE_inc_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % incentive bis
        if G_mod_perfE.inc_bis ~= 0
            curr_vec_nm = 'G_mod_perfE_inc_bis';
            
            % positive contrast
            jCon = jCon + 1;
            G_perfE_inc_bis_toGain_idx = strcmp(con_names,['G_mod_perfE_inc_bis_toGain',posCon_nm]);
            G_perfE_inc_bis_toLose_idx = strcmp(con_names,['G_mod_perfE_inc_bis_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_perfE_inc_bis_toGain_idx,:) + con_vec(G_perfE_inc_bis_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % effort exertion
        if G_mod_perfE.Effort ~= 0
            curr_vec_nm = 'G_mod_perfE_Effort';
            
            % positive contrast
            jCon = jCon + 1;
            G_perfE_E_toGain_idx = strcmp(con_names,['G_mod_perfE_Effort_toGain',posCon_nm]);
            G_perfE_E_toLose_idx = strcmp(con_names,['G_mod_perfE_Effort_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_perfE_E_toGain_idx,:) + con_vec(G_perfE_E_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % trial number
        if G_mod_perfE.trialN ~= 0
            curr_vec_nm = 'G_mod_perfE_trialN';
            
            % positive contrast
            jCon = jCon + 1;
            G_perfE_trialN_toGain_idx = strcmp(con_names,['G_mod_perfE_trialN_toGain',posCon_nm]);
            G_perfE_trialN_toLose_idx = strcmp(con_names,['G_mod_perfE_trialN_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_perfE_trialN_toGain_idx,:) + con_vec(G_perfE_trialN_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
end


%% feedback
G_o_fbk = gripRPprm.o_fbk;
G_mod_fbk = gripRPprm.mod_fbk;
switch G_o_fbk
    case 2
        % feedback
        if G_mod_fbk.fbk ~= 0
            curr_vec_nm = 'G_mod_fbk_fbk';
            
            % positive contrast
            jCon = jCon + 1;
            G_fbk_fbk_toGain_idx = strcmp(con_names,['G_mod_fbk_fbk_toGain',posCon_nm]);
            G_fbk_fbk_toLose_idx = strcmp(con_names,['G_mod_fbk_fbk_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_fbk_fbk_toGain_idx,:) + con_vec(G_fbk_fbk_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % total gain
        if G_mod_fbk.totalGain ~= 0
            curr_vec_nm = 'G_mod_fbk_totalGain';
            
            % positive contrast
            jCon = jCon + 1;
            G_fbk_totalGain_toGain_idx = strcmp(con_names,['G_mod_fbk_totalGain_toGain',posCon_nm]);
            G_fbk_totalGain_toLose_idx = strcmp(con_names,['G_mod_fbk_totalGain_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_fbk_totalGain_toGain_idx,:) + con_vec(G_fbk_totalGain_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ressource
        if G_mod_fbk.X_pred ~= 0
            curr_vec_nm = 'G_mod_fbk_X_pred';
            
            % positive contrast
            jCon = jCon + 1;
            G_fbk_X_pred_toGain_idx = strcmp(con_names,['G_mod_fbk_X_pred_toGain',posCon_nm]);
            G_fbk_X_pred_toLose_idx = strcmp(con_names,['G_mod_fbk_X_pred_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_fbk_X_pred_toGain_idx,:) + con_vec(G_fbk_X_pred_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % performance
        if G_mod_fbk.perf ~= 0
            curr_vec_nm = 'G_mod_fbk_perf';
            
            % positive contrast
            jCon = jCon + 1;
            G_fbk_perf_toGain_idx = strcmp(con_names,['G_mod_fbk_perf_toGain',posCon_nm]);
            G_fbk_perf_toLose_idx = strcmp(con_names,['G_mod_fbk_perf_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_fbk_perf_toGain_idx,:) + con_vec(G_fbk_perf_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ROI activity
        if G_mod_fbk.ROI_activity_yn ~= 0
            ROI_full_nm = ['G_mod_fbk_GLM',G_mod_fbk.ROI_activity_GLM,...
                    '_',G_mod_fbk.ROI_activity_period,'_period_',G_mod_fbk.ROI_activity_ROI_nm,'_activity'];
            curr_vec_nm = ROI_full_nm;
            
            % positive contrast
            jCon = jCon + 1;
            G_fbk_ROIactivity_toGain_idx = strcmp(con_names,[ROI_full_nm,'_toGain',posCon_nm]);
            G_fbk_ROIactivity_toLose_idx = strcmp(con_names,[ROI_full_nm,'_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_fbk_ROIactivity_toGain_idx,:) + con_vec(G_fbk_ROIactivity_toLose_idx,:); % pool effort to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % trial number
        if G_mod_fbk.trialN ~= 0
            curr_vec_nm = 'G_mod_fbk_trialN';
            
            % positive contrast
            jCon = jCon + 1;
            G_fbk_trialN_toGain_idx = strcmp(con_names,['G_mod_fbk_trialN_toGain',posCon_nm]);
            G_fbk_trialN_toLose_idx = strcmp(con_names,['G_mod_fbk_trialN_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_fbk_trialN_toGain_idx,:) + con_vec(G_fbk_trialN_toLose_idx,:); % pool incentive to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
end


end % function