function [ jCon, con_vec, con_names ] = MS2_contrasts_combinations_GS( jCon, con_vec, con_names,...
    posCon_nm, n_regs, gripRPprm, stroopRPprm)
% [ jCon, con_vec, con_names ] = MS2_contrasts_combinations_GS( jCon, con_vec, con_names,...
%     posCon_nm, n_regs, gripRPprm, stroopRPprm) combine similar contrasts in the
%     grip+Stroop task.


%% incentive
G_o_inc = gripRPprm.o_inc;
G_mod_inc = gripRPprm.mod_inc;
S_o_inc = stroopRPprm.o_inc;
S_mod_inc = stroopRPprm.mod_inc;
%% incentive
if G_o_inc == 1 && S_o_inc == 1 % all trials pooled => pool tasks
    
    % run2 - run1
    curr_vec_nm = 'GSpool_o_inc_run2_min_run1';
    % positive contrast
    jCon = jCon + 1;
    G_oInc_r2_min_r1_idx = strcmp(con_names,['G_o_inc_run2_min_run1',posCon_nm]);
    S_oInc_r2_min_r1_idx = strcmp(con_names,['S_o_inc_run2_min_run1',posCon_nm]);
    con_vec(jCon, 1:n_regs) = (con_vec(G_oInc_r2_min_r1_idx,:) + con_vec(S_oInc_r2_min_r1_idx,:))/2; % pool incentive to gain and to lose trials
    con_names{jCon} = [curr_vec_nm,posCon_nm];
    % negative contrast
    [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    
    
    % compare Grip and Stroop
    curr_vec_nm = 'G_min_S_o_inc';
    
    % positive contrast
    jCon = jCon + 1;
    G_OInc_idx = strcmp(con_names,['G_o_inc',posCon_nm]);
    S_OInc_idx = strcmp(con_names,['S_o_inc',posCon_nm]);
    con_vec(jCon, 1:n_regs) = (con_vec(G_OInc_idx,:) - con_vec(S_OInc_idx,:))/2; % pool incentive to gain and to lose trials
    con_names{jCon} = [curr_vec_nm,posCon_nm];
    
    % negative contrast
    [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    
    
    % GL condition
    if G_mod_inc.GL_cond ~= 0 && S_mod_inc.GL_cond ~= 0
        curr_vec_nm = 'GSpool_mod_inc_GL_cond';
        
        % positive contrast
        jCon = jCon + 1;
        G_inc_GLcond_idx = strcmp(con_names,['G_mod_inc_GL_cond',posCon_nm]);
        S_inc_GLcond_idx = strcmp(con_names,['S_mod_inc_GL_cond',posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_GLcond_idx,:) + con_vec(S_inc_GLcond_idx,:))/2; % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % reward type
    if G_mod_inc.R_type ~= 0 && S_mod_inc.R_type ~= 0
        
        curr_vec_nm = 'GSpool_mod_inc_R_type';
        
        % positive contrast
        jCon = jCon + 1;
        G_Rtype_idx = strcmp(con_names,['G_mod_inc_R_type',posCon_nm]);
        S_Rtype_idx = strcmp(con_names,['S_mod_inc_R_type',posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_Rtype_idx,:) + con_vec(S_Rtype_idx,:))/2; % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % incentive
    if G_mod_inc.inc ~= 0 && S_mod_inc.inc ~= 0
        
        curr_vec_nm = 'GSpool_mod_inc_inc';
        
        % positive contrast
        jCon = jCon + 1;
        G_inc_idx = strcmp(con_names,['G_mod_inc_inc',posCon_nm]);
        S_inc_idx = strcmp(con_names,['S_mod_inc_inc',posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_idx,:) + con_vec(S_inc_idx,:))/2; % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % incentive bis
    if G_mod_inc.inc_bis ~= 0 && S_mod_inc.inc_bis ~= 0
        
        curr_vec_nm = 'GSpool_mod_inc_inc_bis';
        
        % positive contrast
        jCon = jCon + 1;
        G_inc_bis_idx = strcmp(con_names,['G_mod_inc_inc_bis',posCon_nm]);
        S_inc_bis_idx = strcmp(con_names,['S_mod_inc_inc_bis',posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_bis_idx,:) + con_vec(S_inc_bis_idx,:))/2; % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end

    %% confidence
    if G_mod_inc.conf ~= 0 && S_mod_inc.conf ~= 0
        
        curr_vec_nm = 'GSpool_mod_inc_conf';
        
        % positive contrast
        jCon = jCon + 1;
        G_conf_idx = strcmp(con_names,['G_mod_inc_conf',posCon_nm]);
        S_conf_idx = strcmp(con_names,['S_mod_inc_conf',posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_conf_idx,:) + con_vec(S_conf_idx,:))/2; % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end

    %% ROI activity
    if G_mod_inc.ROI_activity_yn ~= 0 && S_mod_inc.ROI_activity_yn ~= 0
        G_ROI_full_nm = ['mod_inc_GLM',G_mod_inc.ROI_activity_GLM,...
            '_',G_mod_inc.ROI_activity_period,'_period_',G_mod_inc.ROI_activity_ROI_nm,'_activity'];
        S_ROI_full_nm = ['mod_inc_GLM',S_mod_inc.ROI_activity_GLM,...
            '_',S_mod_inc.ROI_activity_period,'_period_',S_mod_inc.ROI_activity_ROI_nm,'_activity'];
        if strcmp(G_ROI_full_nm,S_ROI_full_nm)
            ROI_full_nm = G_ROI_full_nm;
        else
            ROI_full_nm = ['G_',G_ROI_full_nm,'_and_S_',S_ROI_full_nm];
        end
        
        curr_vec_nm = ['GSpool_',ROI_full_nm];
        
        % positive contrast
        jCon = jCon + 1;
        G_inc_ROIactivity_idx = strcmp(con_names,['G_',G_ROI_full_nm,posCon_nm]);
        S_inc_ROIactivity_idx = strcmp(con_names,['S_',S_ROI_full_nm,posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_ROIactivity_idx,:) + con_vec(S_inc_ROIactivity_idx,:))/2; % pool to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% effort: force+number of pairs solved
    if G_mod_inc.Effort ~= 0 && S_mod_inc.n_pairs_solved ~= 0
        curr_vec_nm = 'GSpool_mod_inc_perf_G_Eff_S_nPairs';
        
        % positive contrast
        jCon = jCon + 1;
        G_inc_Effort_idx = strcmp(con_names,['G_mod_inc_effort',posCon_nm]);
        S_inc_nPairs_idx = strcmp(con_names,['S_mod_inc_n_pairs_solved',posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_Effort_idx,:) + con_vec(S_inc_nPairs_idx,:))/2; % pool to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% effort bis: force*RT first pair
    if G_mod_inc.Effort ~= 0 && S_mod_inc.RT_fp ~= 0
        curr_vec_nm = 'GSpool_mod_inc_perf_G_Eff_S_RTfp';
        
        % positive contrast
        jCon = jCon + 1;
        G_inc_Effort_idx = strcmp(con_names,['G_mod_inc_effort',posCon_nm]);
        S_inc_RTfp_idx = strcmp(con_names,['S_mod_inc_RT_fp',posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_Effort_idx,:) + con_vec(S_inc_RTfp_idx,:))/2; % pool to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% effort ter: force*mean RT all but first pair
    if G_mod_inc.Effort ~= 0 && S_mod_inc.RT_mRT ~= 0
        curr_vec_nm = 'GSpool_mod_inc_perf_G_Eff_S_mRT';
        
        % positive contrast
        jCon = jCon + 1;
        G_inc_Effort_idx = strcmp(con_names,['G_mod_inc_effort',posCon_nm]);
        S_inc_mRT_idx = strcmp(con_names,['S_mod_inc_RT_mRT',posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_Effort_idx,:) + con_vec(S_inc_mRT_idx,:))/2; % pool to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% model variables compute EV = benefit - cost
    if G_mod_inc.mdl_n ~= 0 && G_mod_inc.mdl_cost ~= 0 && G_mod_inc.mdl_benefit ~= 0 &&...
            S_mod_inc.mdl_n ~= 0 && S_mod_inc.mdl_cost ~= 0 && S_mod_inc.mdl_benefit ~= 0
        
        % benefit
        curr_vec_nm = ['GSpool_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_benefit'];
        % positive contrast
        jCon = jCon + 1;
        G_inc_mdlBenef_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_benefit',posCon_nm]);
        S_inc_mdlBenef_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit',posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_mdlBenef_idx,:) + con_vec(S_inc_mdlBenef_idx,:))/2;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % cost
        curr_vec_nm = ['GSpool_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_cost'];
        % positive contrast
        jCon = jCon + 1;
        G_inc_mdlCost_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_cost',posCon_nm]);
        S_inc_mdlCost_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost',posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_mdlCost_idx,:) + con_vec(S_inc_mdlCost_idx,:))/2;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % benefit - cost
        curr_vec_nm = ['GSpool_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_Benefit_min_Cost'];
        % positive contrast
        jCon = jCon + 1;
        G_inc_mdlBenef_min_Cost_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_Benefit_min_Cost',posCon_nm]);
        S_inc_mdlBenef_min_Cost_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_Benefit_min_Cost',posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_mdlBenef_min_Cost_idx,:) + con_vec(S_inc_mdlBenef_min_Cost_idx,:))/2;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    
    %% RT first squeeze/first pair
    if G_mod_inc.RT_fp ~= 0 && S_mod_inc.RT_fp ~= 0
        
        curr_vec_nm = 'GSpool_mod_inc_RT_fp';
        
        % positive contrast
        jCon = jCon + 1;
        G_inc_RT_fp_idx = strcmp(con_names,['G_mod_inc_RT_fp',posCon_nm]);
        S_inc_RT_fp_idx = strcmp(con_names,['S_mod_inc_RT_fp',posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_RT_fp_idx,:) + con_vec(S_inc_RT_fp_idx,:))/2; % pool to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% total gain until now
    if G_mod_inc.totalGain_prev ~= 0 && S_mod_inc.totalGain_prev ~= 0
        
        curr_vec_nm = 'GSpool_mod_inc_totalGain_prev';
        
        % positive contrast
        jCon = jCon + 1;
        G_inc_totalGainPrev_idx = strcmp(con_names,['G_mod_inc_totalGain_prev',posCon_nm]);
        S_inc_totalGainPrev_idx = strcmp(con_names,['S_mod_inc_totalGain_prev',posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_totalGainPrev_idx,:) + con_vec(S_inc_totalGainPrev_idx,:))/2; % pool to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% sum(performance previous trials)
    if G_mod_inc.sumPerfPrev ~= 0 && S_mod_inc.sumPerfPrev ~= 0
        
        curr_vec_nm = 'GSpool_mod_inc_sumPerfPrev';
        
        % positive contrast
        jCon = jCon + 1;
        G_inc_sumPerfPrev_idx = strcmp(con_names,['G_mod_inc_sumPerfPrev',posCon_nm]);
        S_inc_sumPerfPrev_idx = strcmp(con_names,['S_mod_inc_sumPerfPrev',posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_sumPerfPrev_idx,:) + con_vec(S_inc_sumPerfPrev_idx,:))/2; % pool to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% trial number
    if G_mod_inc.trialN ~= 0 && S_mod_inc.trialN ~= 0
        
        curr_vec_nm = 'GSpool_mod_inc_trialN';
        
        % positive contrast
        jCon = jCon + 1;
        G_inc_trialN_idx = strcmp(con_names,['G_mod_inc_trialN',posCon_nm]);
        S_inc_trialN_idx = strcmp(con_names,['S_mod_inc_trialN',posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_trialN_idx,:) + con_vec(S_inc_trialN_idx,:))/2; % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% ressource
    if G_mod_inc.X_pred ~= 0 && S_mod_inc.X_pred ~= 0
        %% pool grip+stroop
        curr_vec_nm = 'GSpool_mod_inc_X_pred';
        
        % positive contrast
        jCon = jCon + 1;
        G_inc_X_pred_idx = strcmp(con_names,['G_mod_inc_X_pred',posCon_nm]);
        S_inc_X_pred_idx = strcmp(con_names,['S_mod_inc_X_pred',posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_X_pred_idx,:) + con_vec(S_inc_X_pred_idx,:))/2; % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        %% compare Grip>Stroop
        curr_vec_nm = 'G_min_S_mod_inc_X_pred';
        
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_X_pred_idx,:) - con_vec(S_inc_X_pred_idx,:))/2; % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% performance
    if G_mod_inc.perf ~= 0 && S_mod_inc.perf ~= 0
        %% Grip+Stroop
        curr_vec_nm = 'GSpool_mod_inc_perf';
        
        % positive contrast
        jCon = jCon + 1;
        G_inc_perf_idx = strcmp(con_names,['G_mod_inc_perf',posCon_nm]);
        S_inc_perf_idx = strcmp(con_names,['S_mod_inc_perf',posCon_nm]);
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_perf_idx,:) + con_vec(S_inc_perf_idx,:))/2; % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% Grip>Stroop
        curr_vec_nm = 'G_min_S_mod_inc_perf';
        
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_perf_idx,:) - con_vec(S_inc_perf_idx,:))/2; % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
elseif G_o_inc == 2 && S_o_inc == 2 % to gain/to loss => pool tasks and conditions
    
    % run2 - run1 - toGain
    curr_vec_nm = 'GSpool_o_inc_toGain_run2_min_run1';
    % positive contrast
    jCon = jCon + 1;
    G_oInc_toGain_r2_min_r1_idx = strcmp(con_names,['G_o_inc_toGain_run2_min_run1',posCon_nm]);
    S_oInc_toGain_r2_min_r1_idx = strcmp(con_names,['S_o_inc_toGain_run2_min_run1',posCon_nm]);
    con_vec(jCon, 1:n_regs) = con_vec(G_oInc_toGain_r2_min_r1_idx,:) + con_vec(S_oInc_toGain_r2_min_r1_idx,:); % pool incentive to gain and to lose trials
    con_names{jCon} = [curr_vec_nm,posCon_nm];
    % negative contrast
    [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    
    % run2 - run1 - toLose
    curr_vec_nm = 'GSpool_o_inc_toLose_run2_min_run1';
    % positive contrast
    jCon = jCon + 1;
    G_oInc_toLose_r2_min_r1_idx = strcmp(con_names,['G_o_inc_toLose_run2_min_run1',posCon_nm]);
    S_oInc_toLose_r2_min_r1_idx = strcmp(con_names,['S_o_inc_toLose_run2_min_run1',posCon_nm]);
    con_vec(jCon, 1:n_regs) = con_vec(G_oInc_toLose_r2_min_r1_idx,:) + con_vec(S_oInc_toLose_r2_min_r1_idx,:); % pool incentive to gain and to lose trials
    con_names{jCon} = [curr_vec_nm,posCon_nm];
    % negative contrast
    [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    
    % run2 - run1
    curr_vec_nm = 'GSpool_o_inc_run2_min_run1';
    % positive contrast
    jCon = jCon + 1;
    G_oInc_r2_min_r1_idx = strcmp(con_names,['G_o_inc_run2_min_run1',posCon_nm]);
    S_oInc_r2_min_r1_idx = strcmp(con_names,['S_o_inc_run2_min_run1',posCon_nm]);
    con_vec(jCon, 1:n_regs) = con_vec(G_oInc_r2_min_r1_idx,:) + con_vec(S_oInc_r2_min_r1_idx,:); % pool incentive to gain and to lose trials
    con_names{jCon} = [curr_vec_nm,posCon_nm];
    % negative contrast
    [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    
    % compare incentive Gain vs incentive Loss period
    % pool gain and loss condition for incentive modulation
    curr_vec_nm = 'GSpool_o_inc_Gain_min_Loss';
    % positive contrast
    jCon = jCon + 1;
    G_Oinc_toGain_idx = strcmp(con_names,['G_o_inc_toGain',posCon_nm]);
    G_Oinc_toLose_idx = strcmp(con_names,['G_o_inc_toLose',posCon_nm]);
    S_Oinc_toGain_idx = strcmp(con_names,['S_o_inc_toGain',posCon_nm]);
    S_Oinc_toLose_idx = strcmp(con_names,['S_o_inc_toLose',posCon_nm]);
    con_vec(jCon, 1:n_regs) = con_vec(G_Oinc_toGain_idx,:) + con_vec(S_Oinc_toGain_idx,:)+...
        -con_vec(G_Oinc_toLose_idx,:) -con_vec(S_Oinc_toLose_idx,:); % pool incentive to gain and to lose trials
    con_names{jCon} = [curr_vec_nm,posCon_nm];
    % negative contrast
    [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    
    % incentive
    if G_mod_inc.R_type ~= 0 && S_mod_inc.R_type ~= 0
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_inc_R_type';
        % positive contrast
        jCon = jCon + 1;
        G_R_type_toGain_idx = strcmp(con_names,['G_mod_inc_R_type_toGain',posCon_nm]);
        G_R_type_toLose_idx = strcmp(con_names,['G_mod_inc_R_type_toLose',posCon_nm]);
        S_R_type_toGain_idx = strcmp(con_names,['S_mod_inc_R_type_toGain',posCon_nm]);
        S_R_type_toLose_idx = strcmp(con_names,['S_mod_inc_R_type_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_R_type_toGain_idx,:) + con_vec(S_R_type_toGain_idx,:) +...
            con_vec(G_R_type_toLose_idx,:) + con_vec(S_R_type_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % pool tasks together but keep valence (to gain/to lose trials)
        curr_vec_nm = 'GSpool_mod_inc_R_type_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_R_type_toGain_idx,:) + con_vec(S_R_type_toGain_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        curr_vec_nm = 'GSpool_mod_inc_R_type_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_R_type_toLose_idx,:) + con_vec(S_R_type_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % comparison gain and loss condition for incentive modulation
        curr_vec_nm = 'GSpool_mod_inc_R_type_Gain_min_Loss';
        
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_R_type_toGain_idx,:) + con_vec(S_R_type_toGain_idx,:)+...
            -con_vec(G_R_type_toLose_idx,:) -con_vec(S_R_type_toLose_idx,:); % compare incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end % incentive
    
    %% incentive
    if G_mod_inc.inc ~= 0 && S_mod_inc.inc ~= 0
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_inc_inc';
        % positive contrast
        jCon = jCon + 1;
        G_inc_toGain_idx = strcmp(con_names,['G_mod_inc_inc_toGain',posCon_nm]);
        G_inc_toLose_idx = strcmp(con_names,['G_mod_inc_inc_toLose',posCon_nm]);
        S_inc_toGain_idx = strcmp(con_names,['S_mod_inc_inc_toGain',posCon_nm]);
        S_inc_toLose_idx = strcmp(con_names,['S_mod_inc_inc_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_toGain_idx,:) + con_vec(S_inc_toGain_idx,:) +...
            con_vec(G_inc_toLose_idx,:) + con_vec(S_inc_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % pool tasks together but keep valence (to gain/to lose trials)
        curr_vec_nm = 'GSpool_mod_inc_inc_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_toGain_idx,:) + con_vec(S_inc_toGain_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        curr_vec_nm = 'GSpool_mod_inc_inc_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_toLose_idx,:) + con_vec(S_inc_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % comparison gain and loss condition for incentive modulation
        curr_vec_nm = 'GSpool_mod_inc_inc_Gain_min_Loss';
        
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_inc_toGain_idx,:) + con_vec(S_inc_toGain_idx,:)+...
            -con_vec(G_inc_toLose_idx,:) -con_vec(S_inc_toLose_idx,:); % compare incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end % incentive
    
    %% incentive bis
    if G_mod_inc.inc_bis ~= 0 && S_mod_inc.inc_bis ~= 0
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_inc_inc_bis';
        % positive contrast
        jCon = jCon + 1;
        G_inc_bis_toGain_idx = strcmp(con_names,['G_mod_inc_inc_bis_toGain',posCon_nm]);
        G_inc_bis_toLose_idx = strcmp(con_names,['G_mod_inc_inc_bis_toLose',posCon_nm]);
        S_inc_bis_toGain_idx = strcmp(con_names,['S_mod_inc_inc_bis_toGain',posCon_nm]);
        S_inc_bis_toLose_idx = strcmp(con_names,['S_mod_inc_inc_bis_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_bis_toGain_idx,:) + con_vec(S_inc_bis_toGain_idx,:) +...
            con_vec(G_inc_bis_toLose_idx,:) + con_vec(S_inc_bis_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % pool tasks together but keep valence (to gain/to lose trials)
        curr_vec_nm = 'GSpool_mod_inc_inc_bis_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_bis_toGain_idx,:) + con_vec(S_inc_bis_toGain_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        curr_vec_nm = 'GSpool_mod_inc_inc_bis_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_bis_toLose_idx,:) + con_vec(S_inc_bis_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % comparison gain and loss condition for incentive modulation
        curr_vec_nm = 'GSpool_mod_inc_inc_bis_Gain_min_Loss';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_inc_bis_toGain_idx,:) + con_vec(S_inc_bis_toGain_idx,:)+...
            -con_vec(G_inc_bis_toLose_idx,:) -con_vec(S_inc_bis_toLose_idx,:); % compare incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end % incentive

    %% confidence   
    if G_mod_inc.conf ~= 0 && S_mod_inc.conf ~= 0
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_inc_conf';
        % positive contrast
        jCon = jCon + 1;
        G_conf_toGain_idx = strcmp(con_names,['G_mod_inc_conf_toGain',posCon_nm]);
        G_conf_toLose_idx = strcmp(con_names,['G_mod_inc_conf_toLose',posCon_nm]);
        S_conf_toGain_idx = strcmp(con_names,['S_mod_inc_conf_toGain',posCon_nm]);
        S_conf_toLose_idx = strcmp(con_names,['S_mod_inc_conf_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_conf_toGain_idx,:) + con_vec(S_conf_toGain_idx,:) +...
            con_vec(G_conf_toLose_idx,:) + con_vec(S_conf_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % pool tasks together but keep valence (to gain/to lose trials)
        curr_vec_nm = 'GSpool_mod_inc_conf_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_conf_toGain_idx,:) + con_vec(S_conf_toGain_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        curr_vec_nm = 'GSpool_mod_inc_conf_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_conf_toLose_idx,:) + con_vec(S_conf_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % comparison gain and loss condition for incentive modulation
        curr_vec_nm = 'GSpool_mod_inc_conf_Gain_min_Loss';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_conf_toGain_idx,:) + con_vec(S_conf_toGain_idx,:)+...
            -con_vec(G_conf_toLose_idx,:) -con_vec(S_conf_toLose_idx,:); % compare incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end % incentive
    
    %% ROI activity
    if G_mod_inc.ROI_activity_yn ~= 0 && S_mod_inc.ROI_activity_yn ~= 0
        G_ROI_full_nm = ['mod_inc_GLM',G_mod_inc.ROI_activity_GLM,...
            '_',G_mod_inc.ROI_activity_period,'_period_',G_mod_inc.ROI_activity_ROI_nm,'_activity'];
        S_ROI_full_nm = ['mod_inc_GLM',S_mod_inc.ROI_activity_GLM,...
            '_',S_mod_inc.ROI_activity_period,'_period_',S_mod_inc.ROI_activity_ROI_nm,'_activity'];
        if strcmp(G_ROI_full_nm,S_ROI_full_nm)
            ROI_full_nm = G_ROI_full_nm;
        else
            ROI_full_nm = ['G_',G_ROI_full_nm,'_and_S_',S_ROI_full_nm];
        end
        
        G_inc_ROIactivity_toGain_idx = strcmp(con_names,['G_',G_ROI_full_nm,posCon_nm]);
        S_inc_ROIactivity_toGain_idx = strcmp(con_names,['S_',S_ROI_full_nm,posCon_nm]);
        G_inc_ROIactivity_toLose_idx = strcmp(con_names,['G_',G_ROI_full_nm,posCon_nm]);
        S_inc_ROIactivity_toLose_idx = strcmp(con_names,['S_',S_ROI_full_nm,posCon_nm]);
        
        % Gain+Loss
        curr_vec_nm = ['GSpool_',ROI_full_nm];
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_inc_ROIactivity_toGain_idx,:) + con_vec(S_inc_ROIactivity_toGain_idx,:) +...
            con_vec(G_inc_ROIactivity_toLose_idx,:) + con_vec(S_inc_ROIactivity_toLose_idx,:); % pool effort to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % Gain
        curr_vec_nm = ['GSpool_',ROI_full_nm,'_toGain'];
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_inc_ROIactivity_toGain_idx,:) + con_vec(S_inc_ROIactivity_toGain_idx,:); % pool effort to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % Loss
        curr_vec_nm = ['GSpool_',ROI_full_nm,'_toLose'];
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_inc_ROIactivity_toLose_idx,:) + con_vec(S_inc_ROIactivity_toLose_idx,:); % pool effort to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % Gain - Loss
        curr_vec_nm = ['GSpool_',ROI_full_nm,'_Gain_min_Loss'];
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_inc_ROIactivity_toLose_idx,:) + con_vec(S_inc_ROIactivity_toLose_idx,:); % pool effort to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% effort: force*number of pairs solved
    if G_mod_inc.Effort ~= 0 && S_mod_inc.n_pairs_solved ~= 0
        %%
        G_inc_eff_toGain_idx = strcmp(con_names,['G_mod_inc_effort_toGain',posCon_nm]);
        G_inc_eff_toLose_idx = strcmp(con_names,['G_mod_inc_effort_toLose',posCon_nm]);
        S_inc_nPairs_solved_toGain_idx = strcmp(con_names,['S_mod_inc_n_pairs_solved_toGain',posCon_nm]);
        S_inc_nPairs_solved_toLose_idx = strcmp(con_names,['S_mod_inc_n_pairs_solved_toLose',posCon_nm]);
        
        curr_vec_nm = 'GSpool_mod_inc_perf_G_eff_S_nPairs';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_eff_toGain_idx,:) + con_vec(S_inc_nPairs_solved_toGain_idx,:) +...
            con_vec(G_inc_eff_toLose_idx,:) + con_vec(S_inc_nPairs_solved_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% pool tasks together but keep valence (to gain/to lose trials)
        %% to gain
        curr_vec_nm = 'GSpool_mod_inc_perf_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_eff_toGain_idx,:) + con_vec(S_inc_nPairs_solved_toGain_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% to lose
        curr_vec_nm = 'GSpool_mod_inc_perf_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_eff_toLose_idx,:) + con_vec(S_inc_nPairs_solved_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% effort bis: force*RT first pair
    if G_mod_inc.Effort ~= 0 && S_mod_inc.RT_fp ~= 0
        %%
        curr_vec_nm = 'GSpool_mod_inc_perf_G_Eff_S_RTfp';
        
        % positive contrast
        jCon = jCon + 1;
        G_inc_eff_toGain_idx = strcmp(con_names,['G_mod_inc_effort_toGain',posCon_nm]);
        G_inc_eff_toLose_idx = strcmp(con_names,['G_mod_inc_effort_toLose',posCon_nm]);
        S_inc_RTfp_toGain_idx = strcmp(con_names,['S_mod_inc_RT_fp_toGain',posCon_nm]);
        S_inc_RTfp_toLose_idx = strcmp(con_names,['S_mod_inc_RT_fp_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_inc_eff_toGain_idx,:) + con_vec(S_inc_RTfp_toGain_idx,:) +...
            con_vec(G_inc_eff_toLose_idx,:) + con_vec(S_inc_RTfp_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% to gain
        curr_vec_nm = 'GSpool_mod_inc_perf_G_eff_S_RTfp_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_eff_toGain_idx,:) + con_vec(S_inc_RTfp_toGain_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% to lose
        curr_vec_nm = 'GSpool_mod_inc_perf_G_eff_S_RTfp_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_eff_toLose_idx,:) + con_vec(S_inc_RTfp_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% effort ter: force*mean RT (except first pair)
    if G_mod_inc.Effort ~= 0 && S_mod_inc.RT_mRT ~= 0
        %% pool
        curr_vec_nm = 'GSpool_mod_inc_perf_G_eff_S_mRT';
        
        % positive contrast
        jCon = jCon + 1;
        G_inc_eff_toGain_idx = strcmp(con_names,['G_mod_inc_effort_toGain',posCon_nm]);
        G_inc_eff_toLose_idx = strcmp(con_names,['G_mod_inc_effort_toLose',posCon_nm]);
        S_inc_mRT_toGain_idx = strcmp(con_names,['S_mod_inc_RT_mRT_toGain',posCon_nm]);
        S_inc_mRT_toLose_idx = strcmp(con_names,['S_mod_inc_RT_mRT_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_inc_eff_toGain_idx,:) + con_vec(S_inc_mRT_toGain_idx,:) +...
            con_vec(G_inc_eff_toLose_idx,:) + con_vec(S_inc_mRT_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% to gain
        curr_vec_nm = 'GSpool_mod_inc_perf_G_eff_S_mRT_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_eff_toGain_idx,:) + con_vec(S_inc_mRT_toGain_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% to lose
        curr_vec_nm = 'GSpool_mod_inc_perf_G_eff_S_mRT_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_eff_toLose_idx,:) + con_vec(S_inc_mRT_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% model variables
    if G_mod_inc.mdl_n ~= 0 && S_mod_inc.mdl_n ~= 0
        
        % pool effort to gain and to lose trials
        % cost
        if G_mod_inc.mdl_cost ~= 0 && S_mod_inc.mdl_cost ~= 0
            %% to gain
            curr_vec_nm = ['GSpool_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_cost_toGain'];
            % positive contrast
            jCon = jCon + 1;
            G_inc_cost_toGain_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_cost_toGain',posCon_nm]);
            S_inc_cost_toGain_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost_toGain',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_cost_toGain_idx,:) + con_vec(S_inc_cost_toGain_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            %% to lose
            curr_vec_nm = ['GSpool_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_cost_toLose'];
            % positive contrast
            jCon = jCon + 1;
            G_inc_cost_toLose_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_cost_toLose',posCon_nm]);
            S_inc_cost_toLose_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_cost_toLose_idx,:) + con_vec(S_inc_cost_toLose_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            %% gain + loss
            curr_vec_nm = ['GSpool_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_cost'];
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_cost_toGain_idx,:) + con_vec(G_inc_cost_toLose_idx,:) +...
                con_vec(S_inc_cost_toGain_idx,:) + con_vec(S_inc_cost_toLose_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        %% benefit
        if G_mod_inc.mdl_benefit ~= 0 && S_mod_inc.mdl_benefit ~= 0
            % to gain
            curr_vec_nm = ['GSpool_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_benefit_toGain'];
            % positive contrast
            jCon = jCon + 1;
            G_inc_benefit_toGain_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_benefit_toGain',posCon_nm]);
            S_inc_benefit_toGain_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit_toGain',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_benefit_toGain_idx,:) + con_vec(S_inc_benefit_toGain_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % to lose
            curr_vec_nm = ['GSpool_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_benefit_toLose'];
            % positive contrast
            jCon = jCon + 1;
            G_inc_benefit_toLose_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_benefit_toLose',posCon_nm]);
            S_inc_benefit_toLose_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_benefit_toLose_idx,:) + con_vec(S_inc_benefit_toLose_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % gain + loss
            curr_vec_nm = ['GSpool_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_benefit'];
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_benefit_toGain_idx,:) + con_vec(G_inc_benefit_toLose_idx,:) +...
                con_vec(S_inc_benefit_toGain_idx,:) + con_vec(S_inc_benefit_toLose_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        %% EV
        if G_mod_inc.mdl_EV ~= 0 && S_mod_inc.mdl_EV ~= 0
            %% to gain
            curr_vec_nm = ['GSpool_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_EV_toGain'];
            % positive contrast
            jCon = jCon + 1;
            G_inc_EV_toGain_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_EV_toGain',posCon_nm]);
            S_inc_EV_toGain_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_EV_toGain',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_EV_toGain_idx,:) + con_vec(S_inc_EV_toGain_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            %% to lose
            curr_vec_nm = ['GSpool_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_EV_toLose'];
            % positive contrast
            jCon = jCon + 1;
            G_inc_EV_toLose_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_EV_toLose',posCon_nm]);
            S_inc_EV_toLose_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_EV_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_EV_toLose_idx,:) + con_vec(S_inc_EV_toLose_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            %% gain + loss
            curr_vec_nm = ['GSpool_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_EV'];
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_EV_toGain_idx,:) + con_vec(G_inc_EV_toLose_idx,:) +...
                con_vec(S_inc_EV_toGain_idx,:) + con_vec(S_inc_EV_toLose_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            %% Grip>Stroop gain + loss
            curr_vec_nm = ['G_min_S_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_EV'];
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = (con_vec(G_inc_EV_toGain_idx,:) + con_vec(G_inc_EV_toLose_idx,:)) +...
                -(con_vec(S_inc_EV_toGain_idx,:) + con_vec(S_inc_EV_toLose_idx,:));
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        %% model variables compute EV = benefit - cost
        if G_mod_inc.mdl_cost ~= 0 && G_mod_inc.mdl_benefit ~= 0
            
            % to gain
            curr_vec_nm = ['GSpool_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_Benefit_min_Cost_toGain'];
            % positive contrast
            jCon = jCon + 1;
            G_inc_toGain_mdlBenef_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_benefit_toGain',posCon_nm]);
            G_inc_toGain_mdlCost_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_cost_toGain',posCon_nm]);
            S_inc_toGain_mdlBenef_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit_toGain',posCon_nm]);
            S_inc_toGain_mdlCost_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost_toGain',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_toGain_mdlBenef_idx,:) - con_vec(G_inc_toGain_mdlCost_idx,:) +...
                con_vec(S_inc_toGain_mdlBenef_idx,:) - con_vec(S_inc_toGain_mdlCost_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            
            % to Lose
            curr_vec_nm = ['GSpool_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_Benefit_min_Cost_toLose'];
            % positive contrast
            jCon = jCon + 1;
            G_inc_toLose_mdlBenef_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_benefit_toLose',posCon_nm]);
            G_inc_toLose_mdlCost_idx = strcmp(con_names,['G_mod_inc_model',num2str(G_mod_inc.mdl_n),'_cost_toLose',posCon_nm]);
            S_inc_toLose_mdlBenef_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_benefit_toLose',posCon_nm]);
            S_inc_toLose_mdlCost_idx = strcmp(con_names,['S_mod_inc_model',num2str(S_mod_inc.mdl_n),'_cost_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_inc_toLose_mdlBenef_idx,:) - con_vec(G_inc_toLose_mdlCost_idx,:) +...
                con_vec(S_inc_toLose_mdlBenef_idx,:) - con_vec(S_inc_toLose_mdlCost_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            
            % pool both
            curr_vec_nm = ['GSpool_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_Benefit_min_Cost'];
            % positive contrast
            jCon = jCon + 1;
            GS_inc_mdlBenef_min_Cost_toGain_idx = strcmp(con_names,...
                ['GSpool_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_Benefit_min_Cost_toGain']);
            GS_inc_mdlBenef_min_Cost_toLose_idx = strcmp(con_names,...
                ['GSpool_mod_inc_Gmodel',num2str(G_mod_inc.mdl_n),'_Smodel',num2str(S_mod_inc.mdl_n),'_Benefit_min_Cost_toLose']);
            con_vec(jCon, 1:n_regs) = con_vec(GS_inc_mdlBenef_min_Cost_toGain_idx,:) + con_vec(GS_inc_mdlBenef_min_Cost_toLose_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
    end
    
    
    %% RT first squeeze/first pair
    if G_mod_inc.RT_fp ~= 0 && S_mod_inc.RT_fp ~= 0
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_inc_RT_fp';
        % positive contrast
        jCon = jCon + 1;
        G_inc_RTfp_toGain_idx = strcmp(con_names,['G_mod_inc_RT_fp_toGain',posCon_nm]);
        G_inc_RTfp_toLose_idx = strcmp(con_names,['G_mod_inc_RT_fp_toLose',posCon_nm]);
        S_inc_RTfp_toGain_idx = strcmp(con_names,['S_mod_inc_RT_fp_toGain',posCon_nm]);
        S_inc_RTfp_toLose_idx = strcmp(con_names,['S_mod_inc_RT_fp_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_RTfp_toGain_idx,:) + con_vec(S_inc_RTfp_toGain_idx,:) +...
            con_vec(G_inc_RTfp_toLose_idx,:) + con_vec(S_inc_RTfp_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % pool tasks together but keep valence (to gain/to lose trials)
        curr_vec_nm = 'GSpool_mod_inc_RT_fp_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_RTfp_toGain_idx,:) + con_vec(S_inc_RTfp_toGain_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        curr_vec_nm = 'GSpool_mod_inc_RT_fp_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_RTfp_toLose_idx,:) + con_vec(S_inc_RTfp_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end % RT fp
    
    %% total gain until now
    if G_mod_inc.totalGain_prev ~= 0 && S_mod_inc.totalGain_prev ~= 0
        
        G_inc_totalGainPrev_toGain_idx = strcmp(con_names,['G_mod_inc_totalGain_prev_toGain',posCon_nm]);
        G_inc_totalGainPrev_toLose_idx = strcmp(con_names,['G_mod_inc_totalGain_prev_toLose',posCon_nm]);
        S_inc_totalGainPrev_toGain_idx = strcmp(con_names,['S_mod_inc_totalGain_prev_toGain',posCon_nm]);
        S_inc_totalGainPrev_toLose_idx = strcmp(con_names,['S_mod_inc_totalGain_prev_toLose',posCon_nm]);
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_inc_totalGain_prev';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_totalGainPrev_toGain_idx,:) + con_vec(S_inc_totalGainPrev_toGain_idx,:) +...
            con_vec(G_inc_totalGainPrev_toLose_idx,:) + con_vec(S_inc_totalGainPrev_toLose_idx,:); % pool to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % pool tasks together but keep valence (to gain/to lose trials)
        curr_vec_nm = 'GSpool_mod_inc_totalGain_prev_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_totalGainPrev_toGain_idx,:) + con_vec(S_inc_totalGainPrev_toGain_idx,:); % pool to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        curr_vec_nm = 'GSpool_mod_inc_totalGain_prev_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_totalGainPrev_toLose_idx,:) + con_vec(S_inc_totalGainPrev_toLose_idx,:); % pool to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end % total gain until now
    
    %% sum(performance previous trials)
    if G_mod_inc.sumPerfPrev ~= 0 && S_mod_inc.sumPerfPrev ~= 0
        
        G_inc_sumPerfPrev_toGain_idx = strcmp(con_names,['G_mod_inc_sumPerfPrev_toGain',posCon_nm]);
        G_inc_sumPerfPrev_toLose_idx = strcmp(con_names,['G_mod_inc_sumPerfPrev_toLose',posCon_nm]);
        S_inc_sumPerfPrev_toGain_idx = strcmp(con_names,['S_mod_inc_sumPerfPrev_toGain',posCon_nm]);
        S_inc_sumPerfPrev_toLose_idx = strcmp(con_names,['S_mod_inc_sumPerfPrev_toLose',posCon_nm]);
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_inc_sumPerfPrev';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_sumPerfPrev_toGain_idx,:) + con_vec(S_inc_sumPerfPrev_toGain_idx,:) +...
            con_vec(G_inc_sumPerfPrev_toLose_idx,:) + con_vec(S_inc_sumPerfPrev_toLose_idx,:); % pool to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % pool tasks together but keep valence (to gain/to lose trials)
        curr_vec_nm = 'GSpool_mod_inc_sumPerfPrev_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_sumPerfPrev_toGain_idx,:) + con_vec(S_inc_sumPerfPrev_toGain_idx,:); % pool to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        curr_vec_nm = 'GSpool_mod_inc_sumPerfPrev_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_sumPerfPrev_toLose_idx,:) + con_vec(S_inc_sumPerfPrev_toLose_idx,:); % pool to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end % total gain until now
    
    %% trial number
    if G_mod_inc.trialN ~= 0 && S_mod_inc.trialN ~= 0
        %%
        G_inc_trialN_toGain_idx = strcmp(con_names,['G_mod_inc_trialN_toGain',posCon_nm]);
        G_inc_trialN_toLose_idx = strcmp(con_names,['G_mod_inc_trialN_toLose',posCon_nm]);
        S_inc_trialN_toGain_idx = strcmp(con_names,['S_mod_inc_trialN_toGain',posCon_nm]);
        S_inc_trialN_toLose_idx = strcmp(con_names,['S_mod_inc_trialN_toLose',posCon_nm]);
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_inc_trialN';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_trialN_toGain_idx,:) + con_vec(S_inc_trialN_toGain_idx,:) +...
            con_vec(G_inc_trialN_toLose_idx,:) + con_vec(S_inc_trialN_toLose_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% pool tasks together but keep valence (to gain/to lose trials)
        curr_vec_nm = 'GSpool_mod_inc_trialN_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_trialN_toGain_idx,:) + con_vec(S_inc_trialN_toGain_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %%
        curr_vec_nm = 'GSpool_mod_inc_trialN_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_inc_trialN_toLose_idx,:) + con_vec(S_inc_trialN_toLose_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end % trial number
    
    %% ressource
    if G_mod_inc.X_pred ~= 0 && S_mod_inc.X_pred ~= 0
        G_inc_X_pred_toGain_idx = strcmp(con_names,['G_mod_inc_X_pred_toGain',posCon_nm]);
        G_inc_X_pred_toLose_idx = strcmp(con_names,['G_mod_inc_X_pred_toLose',posCon_nm]);
        S_inc_X_pred_toGain_idx = strcmp(con_names,['S_mod_inc_X_pred_toGain',posCon_nm]);
        S_inc_X_pred_toLose_idx = strcmp(con_names,['S_mod_inc_X_pred_toLose',posCon_nm]);
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_inc_X_pred';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_inc_X_pred_toGain_idx,:) + con_vec(G_inc_X_pred_toLose_idx,:) +...
            con_vec(S_inc_X_pred_toGain_idx,:) + con_vec(S_inc_X_pred_toLose_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% to gain
        curr_vec_nm = 'GSpool_mod_inc_X_pred_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_inc_X_pred_toGain_idx,:) + con_vec(S_inc_X_pred_toGain_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% to lose
        curr_vec_nm = 'GSpool_mod_inc_X_pred_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_inc_X_pred_toLose_idx,:) + con_vec(S_inc_X_pred_toLose_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% gain - loss
        curr_vec_nm = 'GSpool_mod_inc_X_pred_Gain_min_Loss';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_inc_X_pred_toGain_idx,:) - con_vec(G_inc_X_pred_toLose_idx,:) +...
            con_vec(S_inc_X_pred_toGain_idx,:) - con_vec(S_inc_X_pred_toLose_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% Grip>Stroop gain - loss
        curr_vec_nm = 'Gmin_S_mod_inc_X_pred_Gain_min_Loss';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_X_pred_toGain_idx,:) - con_vec(G_inc_X_pred_toLose_idx,:)) +...
            -(con_vec(S_inc_X_pred_toGain_idx,:) - con_vec(S_inc_X_pred_toLose_idx,:)); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% performance
    if G_mod_inc.perf ~= 0 && S_mod_inc.perf ~= 0
        G_inc_perf_toGain_idx = strcmp(con_names,['G_mod_inc_perf_toGain',posCon_nm]);
        G_inc_perf_toLose_idx = strcmp(con_names,['G_mod_inc_perf_toLose',posCon_nm]);
        S_inc_perf_toGain_idx = strcmp(con_names,['S_mod_inc_perf_toGain',posCon_nm]);
        S_inc_perf_toLose_idx = strcmp(con_names,['S_mod_inc_perf_toLose',posCon_nm]);
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_inc_perf';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_inc_perf_toGain_idx,:) + con_vec(G_inc_perf_toLose_idx,:) +...
            con_vec(S_inc_perf_toGain_idx,:) + con_vec(S_inc_perf_toLose_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% to gain
        curr_vec_nm = 'GSpool_mod_inc_perf_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_inc_perf_toGain_idx,:) + con_vec(S_inc_perf_toGain_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% to lose
        curr_vec_nm = 'GSpool_mod_inc_perf_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_inc_perf_toLose_idx,:) + con_vec(S_inc_perf_toLose_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% gain - loss
        curr_vec_nm = 'GSpool_mod_inc_perf_Gain_min_Loss';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_inc_perf_toGain_idx,:) - con_vec(G_inc_perf_toLose_idx,:) +...
            con_vec(S_inc_perf_toGain_idx,:) - con_vec(S_inc_perf_toLose_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% Grip>Stroop gain - loss
        curr_vec_nm = 'G_min_S_mod_inc_perf_Gain_min_Loss';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = (con_vec(G_inc_perf_toGain_idx,:) - con_vec(G_inc_perf_toLose_idx,:)) +...
            -(con_vec(S_inc_perf_toGain_idx,:) - con_vec(S_inc_perf_toLose_idx,:)); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
end


%% effort scale
G_o_dispE = gripRPprm.o_dispE;
G_mod_dispE = gripRPprm.mod_dispE;
S_o_dispE = stroopRPprm.o_dispE;
S_mod_dispE = stroopRPprm.mod_dispE;
if G_o_dispE == 1 && S_o_dispE == 1 % all trials pooled => pool tasks
    
    % run2 - run1
    curr_vec_nm = 'GSpool_o_dispE_run2_min_run1';
    % positive contrast
    jCon = jCon + 1;
    G_oDispE_r2_min_r1_idx = strcmp(con_names,['G_o_dispE_run2_min_run1',posCon_nm]);
    S_oDispE_r2_min_r1_idx = strcmp(con_names,['S_o_dispE_run2_min_run1',posCon_nm]);
    con_vec(jCon, 1:n_regs) = con_vec(G_oDispE_r2_min_r1_idx,:) + con_vec(S_oDispE_r2_min_r1_idx,:); % pool incentive to gain and to lose trials
    con_names{jCon} = [curr_vec_nm,posCon_nm];
    % negative contrast
    [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    
    
    % compare Grip and Stroop
    curr_vec_nm = 'G_min_S_o_dispE';
    
    % positive contrast
    jCon = jCon + 1;
    G_OdispE_idx = strcmp(con_names,['G_o_dispE',posCon_nm]);
    S_OdispE_idx = strcmp(con_names,['S_o_dispE',posCon_nm]);
    con_vec(jCon, 1:n_regs) = con_vec(G_OdispE_idx,:) - con_vec(S_OdispE_idx,:); % pool incentive to gain and to lose trials
    con_names{jCon} = [curr_vec_nm,posCon_nm];
    
    % negative contrast
    [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    
    %% GL condition
    if G_mod_dispE.GL_cond ~= 0 && S_mod_dispE.GL_cond ~= 0
        %%
        curr_vec_nm = 'GSpool_mod_dispE_GL_cond';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_GLcond_idx = strcmp(con_names,['G_mod_dispE_GL_cond',posCon_nm]);
        S_dispE_GLcond_idx = strcmp(con_names,['S_mod_dispE_GL_cond',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_GLcond_idx,:) + con_vec(S_dispE_GLcond_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% reward type (piece money/bill)
    if G_mod_dispE.R_type ~= 0 && S_mod_dispE.R_type ~= 0
        
        curr_vec_nm = 'GSpool_mod_dispE_R_type';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_R_type_idx = strcmp(con_names,['G_mod_dispE_R_type',posCon_nm]);
        S_dispE_R_type_idx = strcmp(con_names,['S_mod_dispE_R_type',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_R_type_idx,:) + con_vec(S_dispE_R_type_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% incentive
    if G_mod_dispE.inc ~= 0 && S_mod_dispE.inc ~= 0
        
        curr_vec_nm = 'GSpool_mod_dispE_inc';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_inc_idx = strcmp(con_names,['G_mod_dispE_inc',posCon_nm]);
        S_dispE_inc_idx = strcmp(con_names,['S_mod_dispE_inc',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_inc_idx,:) + con_vec(S_dispE_inc_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% incentive bis
    if G_mod_dispE.inc_bis ~= 0 && S_mod_dispE.inc_bis ~= 0
        
        curr_vec_nm = 'GSpool_mod_dispE_inc_bis';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_inc_bis_idx = strcmp(con_names,['G_mod_dispE_inc_bis',posCon_nm]);
        S_dispE_inc_bis_idx = strcmp(con_names,['S_mod_dispE_inc_bis',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_inc_bis_idx,:) + con_vec(S_dispE_inc_bis_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end

    %% confidence
    if G_mod_dispE.conf ~= 0 && S_mod_dispE.conf ~= 0
        
        curr_vec_nm = 'GSpool_mod_dispE_conf';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_conf_idx = strcmp(con_names,['G_mod_dispE_conf',posCon_nm]);
        S_dispE_conf_idx = strcmp(con_names,['S_mod_dispE_conf',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_conf_idx,:) + con_vec(S_dispE_conf_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% ROI activity
    if G_mod_dispE.ROI_activity_yn ~= 0 && S_mod_dispE.ROI_activity_yn ~= 0
        G_ROI_full_nm = ['mod_dispE_GLM',G_mod_dispE.ROI_activity_GLM,...
            '_',G_mod_dispE.ROI_activity_period,'_period_',G_mod_dispE.ROI_activity_ROI_nm,'_activity'];
        S_ROI_full_nm = ['mod_dispE_GLM',S_mod_dispE.ROI_activity_GLM,...
            '_',S_mod_dispE.ROI_activity_period,'_period_',S_mod_dispE.ROI_activity_ROI_nm,'_activity'];
        if strcmp(G_ROI_full_nm,S_ROI_full_nm)
            ROI_full_nm = G_ROI_full_nm;
        else
            ROI_full_nm = ['G_',G_ROI_full_nm,'_and_S_',S_ROI_full_nm];
        end
        
        curr_vec_nm = ['GSpool_',ROI_full_nm];
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_ROIactivity_idx = strcmp(con_names,['G_',G_ROI_full_nm,posCon_nm]);
        S_dispE_ROIactivity_idx = strcmp(con_names,['S_',S_ROI_full_nm,posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_ROIactivity_idx,:) + con_vec(S_dispE_ROIactivity_idx,:); % pool effort to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% model variables compute EV = benefit - cost
    if G_mod_dispE.mdl_n ~= 0 && G_mod_dispE.mdl_cost ~= 0 && G_mod_dispE.mdl_benefit ~= 0 &&...
            S_mod_dispE.mdl_n ~= 0 && S_mod_dispE.mdl_cost ~= 0 && S_mod_dispE.mdl_benefit ~= 0
        
        curr_vec_nm = ['GSpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_Benefit_min_Cost'];
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_mdlBenef_min_Cost_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_Benefit_min_Cost',posCon_nm]);
        S_dispE_mdlBenef_min_Cost_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_Benefit_min_Cost',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_mdlBenef_min_Cost_idx,:) + con_vec(S_dispE_mdlBenef_min_Cost_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% RT first press
    if G_mod_dispE.RT_fp ~= 0 && S_mod_dispE.RT_fp ~= 0
        curr_vec_nm = 'GSpool_mod_dispE_RT_fp';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_RTfp_idx = strcmp(con_names,['G_mod_dispE_RT_fp',posCon_nm]);
        S_dispE_RTfp_idx = strcmp(con_names,['S_mod_dispE_RT_fp',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_RTfp_idx,:) + con_vec(S_dispE_RTfp_idx,:); % pool RT first press grip + stroop trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% trial number
    if G_mod_dispE.trialN ~= 0 && S_mod_dispE.trialN ~= 0
        curr_vec_nm = 'GSpool_mod_dispE_trialN';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_trialN_idx = strcmp(con_names,['G_mod_dispE_trialN',posCon_nm]);
        S_dispE_trialN_idx = strcmp(con_names,['S_mod_dispE_trialN',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_trialN_idx,:) + con_vec(S_dispE_trialN_idx,:); % pool trial number grip + stroop trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end

    if G_mod_dispE.X_pred ~= 0 && S_mod_dispE.X_pred ~= 0
        %% pool grip+stroop
        curr_vec_nm = 'GSpool_mod_dispE_X_pred';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_X_pred_idx = strcmp(con_names,['G_mod_dispE_X_pred',posCon_nm]);
        S_dispE_X_pred_idx = strcmp(con_names,['S_mod_dispE_X_pred',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_X_pred_idx,:) + con_vec(S_dispE_X_pred_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        %% compare Grip>Stroop
        curr_vec_nm = 'G_min_S_mod_dispE_X_pred';
        
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_X_pred_idx,:) - con_vec(S_dispE_X_pred_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% performance
    if G_mod_dispE.perf ~= 0 && S_mod_dispE.perf ~= 0
        curr_vec_nm = 'GSpool_mod_dispE_perf';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_perf_idx = strcmp(con_names,['G_mod_dispE_perf',posCon_nm]);
        S_dispE_perf_idx = strcmp(con_names,['S_mod_dispE_perf',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_perf_idx,:) + con_vec(S_dispE_perf_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
elseif G_o_dispE == 2 && S_o_dispE == 2 % to gain/to lose => pool tasks and conditions
    
    % run2 - run1 - toGain
    curr_vec_nm = 'GSpool_o_dispE_toGain_run2_min_run1';
    % positive contrast
    jCon = jCon + 1;
    G_oDispE_toGain_r2_min_r1_idx = strcmp(con_names,['G_o_dispE_toGain_run2_min_run1',posCon_nm]);
    S_oDispE_toGain_r2_min_r1_idx = strcmp(con_names,['S_o_dispE_toGain_run2_min_run1',posCon_nm]);
    con_vec(jCon, 1:n_regs) = con_vec(G_oDispE_toGain_r2_min_r1_idx,:) + con_vec(S_oDispE_toGain_r2_min_r1_idx,:); % pool incentive to gain and to lose trials
    con_names{jCon} = [curr_vec_nm,posCon_nm];
    % negative contrast
    [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    
    % run2 - run1 - toLose
    curr_vec_nm = 'GSpool_o_dispE_toLose_run2_min_run1';
    % positive contrast
    jCon = jCon + 1;
    G_oDispE_toLose_r2_min_r1_idx = strcmp(con_names,['G_o_dispE_toLose_run2_min_run1',posCon_nm]);
    S_oDispE_toLose_r2_min_r1_idx = strcmp(con_names,['S_o_dispE_toLose_run2_min_run1',posCon_nm]);
    con_vec(jCon, 1:n_regs) = con_vec(G_oDispE_toLose_r2_min_r1_idx,:) + con_vec(S_oDispE_toLose_r2_min_r1_idx,:); % pool incentive to gain and to lose trials
    con_names{jCon} = [curr_vec_nm,posCon_nm];
    % negative contrast
    [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    
    % run2 - run1
    curr_vec_nm = 'GSpool_o_dispE_run2_min_run1';
    % positive contrast
    jCon = jCon + 1;
    G_oDispE_r2_min_r1_idx = strcmp(con_names,['G_o_dispE_run2_min_run1',posCon_nm]);
    S_oDispE_r2_min_r1_idx = strcmp(con_names,['S_o_dispE_run2_min_run1',posCon_nm]);
    con_vec(jCon, 1:n_regs) = con_vec(G_oDispE_r2_min_r1_idx,:) + con_vec(S_oDispE_r2_min_r1_idx,:); % pool incentive to gain and to lose trials
    con_names{jCon} = [curr_vec_nm,posCon_nm];
    % negative contrast
    [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    
    
    % compare incentive Gain vs incentive Loss period
    % pool gain and loss condition for incentive modulation
    curr_vec_nm = 'GSpool_o_dispE_Gain_min_Loss';
    
    % positive contrast
    jCon = jCon + 1;
    G_OdispE_toGain_idx = strcmp(con_names,['G_o_dispE_toGain',posCon_nm]);
    G_OdispE_toLose_idx = strcmp(con_names,['G_o_dispE_toLose',posCon_nm]);
    S_OdispE_toGain_idx = strcmp(con_names,['S_o_dispE_toGain',posCon_nm]);
    S_OdispE_toLose_idx = strcmp(con_names,['S_o_dispE_toLose',posCon_nm]);
    con_vec(jCon, 1:n_regs) = con_vec(G_OdispE_toGain_idx,:) + con_vec(S_OdispE_toGain_idx,:)+...
        -con_vec(G_OdispE_toLose_idx,:) -con_vec(S_OdispE_toLose_idx,:);
    con_names{jCon} = [curr_vec_nm,posCon_nm];
    
    % negative contrast
    [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    
    % compare Grip and Stroop
    curr_vec_nm = 'G_min_S_o_dispE';
    % positive contrast
    jCon = jCon + 1;
    con_vec(jCon, 1:n_regs) = con_vec(G_OdispE_toGain_idx,:) + con_vec(G_OdispE_toLose_idx,:) +...
        -( con_vec(S_OdispE_toGain_idx,:) + con_vec(S_OdispE_toLose_idx,:) ); % pool incentive to gain and to lose trials
    con_names{jCon} = [curr_vec_nm,posCon_nm];
    % negative contrast
    [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    
    %% reward type (piece/bill)
    if G_mod_dispE.R_type ~= 0 && S_mod_dispE.R_type ~= 0
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_dispE_R_type';
        % positive contrast
        jCon = jCon + 1;
        G_dispE_R_type_toGain_idx = strcmp(con_names,['G_mod_dispE_R_type_toGain',posCon_nm]);
        G_dispE_R_type_toLose_idx = strcmp(con_names,['G_mod_dispE_R_type_toLose',posCon_nm]);
        S_dispE_R_type_toGain_idx = strcmp(con_names,['S_mod_dispE_R_type_toGain',posCon_nm]);
        S_dispE_R_type_toLose_idx = strcmp(con_names,['S_mod_dispE_R_type_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_R_type_toGain_idx,:) + con_vec(S_dispE_R_type_toGain_idx,:) +...
            con_vec(G_dispE_R_type_toLose_idx,:) + con_vec(S_dispE_R_type_toLose_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % pool tasks together but keep valence (to gain/to lose trials)
        curr_vec_nm = 'GSpool_mod_dispE_R_type_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_R_type_toGain_idx,:) + con_vec(S_dispE_R_type_toGain_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        curr_vec_nm = 'GSpool_mod_dispE_R_type_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_R_type_toLose_idx,:) + con_vec(S_dispE_R_type_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% incentive
    if G_mod_dispE.inc ~= 0 && S_mod_dispE.inc ~= 0
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_dispE_inc';
        % positive contrast
        jCon = jCon + 1;
        G_dispE_inc_toGain_idx = strcmp(con_names,['G_mod_dispE_inc_toGain',posCon_nm]);
        G_dispE_inc_toLose_idx = strcmp(con_names,['G_mod_dispE_inc_toLose',posCon_nm]);
        S_dispE_inc_toGain_idx = strcmp(con_names,['S_mod_dispE_inc_toGain',posCon_nm]);
        S_dispE_inc_toLose_idx = strcmp(con_names,['S_mod_dispE_inc_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_inc_toGain_idx,:) + con_vec(S_dispE_inc_toGain_idx,:) +...
            con_vec(G_dispE_inc_toLose_idx,:) + con_vec(S_dispE_inc_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % pool tasks together but keep valence (to gain/to lose trials)
        curr_vec_nm = 'GSpool_mod_dispE_inc_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_inc_toGain_idx,:) + con_vec(S_dispE_inc_toGain_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        curr_vec_nm = 'GSpool_mod_dispE_inc_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_inc_toLose_idx,:) + con_vec(S_dispE_inc_toLose_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% incentive bis
    if G_mod_dispE.inc_bis ~= 0 && S_mod_dispE.inc_bis ~= 0
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_dispE_inc_bis';
        % positive contrast
        jCon = jCon + 1;
        G_dispE_inc_bis_toGain_idx = strcmp(con_names,['G_mod_dispE_inc_bis_toGain',posCon_nm]);
        G_dispE_inc_bis_toLose_idx = strcmp(con_names,['G_mod_dispE_inc_bis_toLose',posCon_nm]);
        S_dispE_inc_bis_toGain_idx = strcmp(con_names,['S_mod_dispE_inc_bis_toGain',posCon_nm]);
        S_dispE_inc_bis_toLose_idx = strcmp(con_names,['S_mod_dispE_inc_bis_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_inc_bis_toGain_idx,:) + con_vec(S_dispE_inc_bis_toGain_idx,:) +...
            con_vec(G_dispE_inc_bis_toLose_idx,:) + con_vec(S_dispE_inc_bis_toLose_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % pool tasks together but keep valence (to gain/to lose trials)
        curr_vec_nm = 'GSpool_mod_dispE_inc_bis_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_inc_bis_toGain_idx,:) + con_vec(S_dispE_inc_bis_toGain_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        curr_vec_nm = 'GSpool_mod_dispE_inc_bis_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_inc_bis_toLose_idx,:) + con_vec(S_dispE_inc_bis_toLose_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% confidence
    if G_mod_dispE.conf ~= 0 && S_mod_dispE.conf ~= 0
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_dispE_conf';
        % positive contrast
        jCon = jCon + 1;
        G_dispE_conf_toGain_idx = strcmp(con_names,['G_mod_dispE_conf_toGain',posCon_nm]);
        G_dispE_conf_toLose_idx = strcmp(con_names,['G_mod_dispE_conf_toLose',posCon_nm]);
        S_dispE_conf_toGain_idx = strcmp(con_names,['S_mod_dispE_conf_toGain',posCon_nm]);
        S_dispE_conf_toLose_idx = strcmp(con_names,['S_mod_dispE_conf_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_conf_toGain_idx,:) + con_vec(S_dispE_conf_toGain_idx,:) +...
            con_vec(G_dispE_conf_toLose_idx,:) + con_vec(S_dispE_conf_toLose_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % pool tasks together but keep valence (to gain/to lose trials)
        curr_vec_nm = 'GSpool_mod_dispE_conf_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_conf_toGain_idx,:) + con_vec(S_dispE_conf_toGain_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        curr_vec_nm = 'GSpool_mod_dispE_conf_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_conf_toLose_idx,:) + con_vec(S_dispE_conf_toLose_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    %% ROI activity
    if G_mod_dispE.ROI_activity_yn ~= 0 && S_mod_dispE.ROI_activity_yn ~= 0
        G_ROI_full_nm = ['mod_dispE_GLM',G_mod_dispE.ROI_activity_GLM,...
            '_',G_mod_dispE.ROI_activity_period,'_period_',G_mod_dispE.ROI_activity_ROI_nm,'_activity'];
        S_ROI_full_nm = ['mod_dispE_GLM',S_mod_dispE.ROI_activity_GLM,...
            '_',S_mod_dispE.ROI_activity_period,'_period_',S_mod_dispE.ROI_activity_ROI_nm,'_activity'];
        if strcmp(G_ROI_full_nm,S_ROI_full_nm)
            ROI_full_nm = G_ROI_full_nm;
        else
            ROI_full_nm = ['G_',G_ROI_full_nm,'_and_S_',S_ROI_full_nm];
        end
        
        G_dispE_ROIactivity_toGain_idx = strcmp(con_names,['G_',G_ROI_full_nm,posCon_nm]);
        S_dispE_ROIactivity_toGain_idx = strcmp(con_names,['S_',S_ROI_full_nm,posCon_nm]);
        G_dispE_ROIactivity_toLose_idx = strcmp(con_names,['G_',G_ROI_full_nm,posCon_nm]);
        S_dispE_ROIactivity_toLose_idx = strcmp(con_names,['S_',S_ROI_full_nm,posCon_nm]);
        
        curr_vec_nm = ['GSpool_',ROI_full_nm];
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_ROIactivity_toGain_idx,:) + con_vec(S_dispE_ROIactivity_toGain_idx,:) +...
            con_vec(G_dispE_ROIactivity_toLose_idx,:) + con_vec(S_dispE_ROIactivity_toLose_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        curr_vec_nm = ['GSpool_',ROI_full_nm,'_toGain'];
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_ROIactivity_toGain_idx,:) + con_vec(S_dispE_ROIactivity_toGain_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        curr_vec_nm = ['GSpool_',ROI_full_nm,'_toLose'];
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_ROIactivity_toLose_idx,:) + con_vec(S_dispE_ROIactivity_toLose_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % Gain - Loss
        curr_vec_nm = ['GSpool_',ROI_full_nm,'_Gain_min_Loss'];
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_ROIactivity_toGain_idx,:) + con_vec(S_dispE_ROIactivity_toGain_idx,:) +...
            -con_vec(G_dispE_ROIactivity_toLose_idx,:) -con_vec(S_dispE_ROIactivity_toLose_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% model variables
    if G_mod_dispE.mdl_n ~= 0 && S_mod_dispE.mdl_n ~= 0
        
        % pool effort to gain and to lose trials
        % cost
        if G_mod_dispE.mdl_cost ~= 0 && S_mod_dispE.mdl_cost ~= 0
            % to gain
            curr_vec_nm = ['GSpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_cost_toGain'];
            % positive contrast
            jCon = jCon + 1;
            G_dispE_cost_toGain_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_cost_toGain',posCon_nm]);
            S_dispE_cost_toGain_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost_toGain',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_cost_toGain_idx,:) + con_vec(S_dispE_cost_toGain_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % to lose
            curr_vec_nm = ['GSpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_cost_toLose'];
            % positive contrast
            jCon = jCon + 1;
            G_dispE_cost_toLose_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_cost_toLose',posCon_nm]);
            S_dispE_cost_toLose_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_cost_toLose_idx,:) + con_vec(S_dispE_cost_toLose_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % gain + loss
            curr_vec_nm = ['GSpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_cost'];
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_cost_toGain_idx,:) + con_vec(G_dispE_cost_toLose_idx,:) +...
                con_vec(S_dispE_cost_toGain_idx,:) + con_vec(S_dispE_cost_toLose_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        %% benefit
        if G_mod_dispE.mdl_benefit ~= 0 && S_mod_dispE.mdl_benefit ~= 0
            % to gain
            curr_vec_nm = ['GSpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_benefit_toGain'];
            % positive contrast
            jCon = jCon + 1;
            G_dispE_benefit_toGain_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_benefit_toGain',posCon_nm]);
            S_dispE_benefit_toGain_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit_toGain',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_benefit_toGain_idx,:) + con_vec(S_dispE_benefit_toGain_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % to lose
            curr_vec_nm = ['GSpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_benefit_toLose'];
            % positive contrast
            jCon = jCon + 1;
            G_dispE_benefit_toLose_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_benefit_toLose',posCon_nm]);
            S_dispE_benefit_toLose_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_benefit_toLose_idx,:) + con_vec(S_dispE_benefit_toLose_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % gain + loss
            curr_vec_nm = ['GSpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_benefit'];
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_benefit_toGain_idx,:) + con_vec(G_dispE_benefit_toLose_idx,:) +...
                con_vec(S_dispE_benefit_toGain_idx,:) + con_vec(S_dispE_benefit_toLose_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        %% EV
        if G_mod_dispE.mdl_EV ~= 0 && S_mod_dispE.mdl_EV ~= 0
            % to gain
            curr_vec_nm = ['GSpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_EV_toGain'];
            % positive contrast
            jCon = jCon + 1;
            G_dispE_EV_toGain_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_EV_toGain',posCon_nm]);
            S_dispE_EV_toGain_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_EV_toGain',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_EV_toGain_idx,:) + con_vec(S_dispE_EV_toGain_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % to lose
            curr_vec_nm = ['GSpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_EV_toLose'];
            % positive contrast
            jCon = jCon + 1;
            G_dispE_EV_toLose_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_EV_toLose',posCon_nm]);
            S_dispE_EV_toLose_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_EV_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_EV_toLose_idx,:) + con_vec(S_dispE_EV_toLose_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % gain + loss
            curr_vec_nm = ['GSpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_EV'];
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_EV_toGain_idx,:) + con_vec(G_dispE_EV_toLose_idx,:) +...
                con_vec(S_dispE_EV_toGain_idx,:) + con_vec(S_dispE_EV_toLose_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        %% model variables compute EV = benefit - cost
        if G_mod_dispE.mdl_cost ~= 0 && G_mod_dispE.mdl_benefit ~= 0
            
            % to gain
            curr_vec_nm = ['GSpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_Benefit_min_Cost_toGain'];
            % positive contrast
            jCon = jCon + 1;
            G_dispE_toGain_mdlBenef_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_benefit_toGain',posCon_nm]);
            G_dispE_toGain_mdlCost_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_cost_toGain',posCon_nm]);
            S_dispE_toGain_mdlBenef_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit_toGain',posCon_nm]);
            S_dispE_toGain_mdlCost_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost_toGain',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_toGain_mdlBenef_idx,:) - con_vec(G_dispE_toGain_mdlCost_idx,:) +...
                con_vec(S_dispE_toGain_mdlBenef_idx,:) - con_vec(S_dispE_toGain_mdlCost_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            
            % to Lose
            curr_vec_nm = ['GSpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_Benefit_min_Cost_toLose'];
            % positive contrast
            jCon = jCon + 1;
            G_dispE_toLose_mdlBenef_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_benefit_toLose',posCon_nm]);
            G_dispE_toLose_mdlCost_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_cost_toLose',posCon_nm]);
            S_dispE_toLose_mdlBenef_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit_toLose',posCon_nm]);
            S_dispE_toLose_mdlCost_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost_toLose',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_toLose_mdlBenef_idx,:) - con_vec(G_dispE_toLose_mdlCost_idx,:) +...
                con_vec(S_dispE_toLose_mdlBenef_idx,:) - con_vec(S_dispE_toLose_mdlCost_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            
            % pool both
            curr_vec_nm = ['GSpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_Benefit_min_Cost'];
            % positive contrast
            jCon = jCon + 1;
            GS_dispE_mdlBenef_min_Cost_toGain_idx = strcmp(con_names,...
                ['GSpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_Benefit_min_Cost_toGain']);
            GS_dispE_mdlBenef_min_Cost_toLose_idx = strcmp(con_names,...
                ['GSpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_Benefit_min_Cost_toLose']);
            con_vec(jCon, 1:n_regs) = con_vec(GS_dispE_mdlBenef_min_Cost_toGain_idx,:) + con_vec(GS_dispE_mdlBenef_min_Cost_toLose_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
    end
    
    %% RT first press
    if G_mod_dispE.RT_fp ~= 0 && S_mod_dispE.RT_fp ~= 0
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_dispE_RT_fp';
        % positive contrast
        jCon = jCon + 1;
        G_dispE_RTfp_toGain_idx = strcmp(con_names,['G_mod_dispE_RT_fp_toGain',posCon_nm]);
        G_dispE_RTfp_toLose_idx = strcmp(con_names,['G_mod_dispE_RT_fp_toLose',posCon_nm]);
        S_dispE_RTfp_toGain_idx = strcmp(con_names,['S_mod_dispE_RT_fp_toGain',posCon_nm]);
        S_dispE_RTfp_toLose_idx = strcmp(con_names,['S_mod_dispE_RT_fp_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_RTfp_toGain_idx,:) + con_vec(S_dispE_RTfp_toGain_idx,:) +...
            con_vec(G_dispE_RTfp_toLose_idx,:) + con_vec(S_dispE_RTfp_toLose_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % pool tasks together but keep valence (to gain/to lose trials)
        curr_vec_nm = 'GSpool_mod_dispE_RT_fp_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_RTfp_toGain_idx,:) + con_vec(S_dispE_RTfp_toGain_idx,:); % pool RT to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        curr_vec_nm = 'GSpool_mod_dispE_RT_fp_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_RTfp_toLose_idx,:) + con_vec(S_dispE_RTfp_toLose_idx,:); % pool RT to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% trial number
    if G_mod_dispE.trialN ~= 0 && S_mod_dispE.trialN ~= 0
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_dispE_trialN';
        % positive contrast
        jCon = jCon + 1;
        G_dispE_trialN_toGain_idx = strcmp(con_names,['G_mod_dispE_trialN_toGain',posCon_nm]);
        G_dispE_trialN_toLose_idx = strcmp(con_names,['G_mod_dispE_trialN_toLose',posCon_nm]);
        S_dispE_trialN_toGain_idx = strcmp(con_names,['S_mod_dispE_trialN_toGain',posCon_nm]);
        S_dispE_trialN_toLose_idx = strcmp(con_names,['S_mod_dispE_trialN_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_trialN_toGain_idx,:) + con_vec(S_dispE_trialN_toGain_idx,:) +...
            con_vec(G_dispE_trialN_toLose_idx,:) + con_vec(S_dispE_trialN_toLose_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % pool tasks together but keep valence (to gain/to lose trials)
        curr_vec_nm = 'GSpool_mod_dispE_trialN_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_trialN_toGain_idx,:) + con_vec(S_dispE_trialN_toGain_idx,:); % pool RT to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        curr_vec_nm = 'GSpool_mod_dispE_trialN_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_dispE_trialN_toLose_idx,:) + con_vec(S_dispE_trialN_toLose_idx,:); % pool RT to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end

    %% ressource
    if G_mod_dispE.X_pred ~= 0 && S_mod_dispE.X_pred ~= 0
        G_dispE_X_pred_toGain_idx = strcmp(con_names,['G_mod_dispE_X_pred_toGain',posCon_nm]);
        G_dispE_X_pred_toLose_idx = strcmp(con_names,['G_mod_dispE_X_pred_toLose',posCon_nm]);
        S_dispE_X_pred_toGain_idx = strcmp(con_names,['S_mod_dispE_X_pred_toGain',posCon_nm]);
        S_dispE_X_pred_toLose_idx = strcmp(con_names,['S_mod_dispE_X_pred_toLose',posCon_nm]);
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_dispE_X_pred';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_X_pred_toGain_idx,:) + con_vec(G_dispE_X_pred_toLose_idx,:) +...
            con_vec(S_dispE_X_pred_toGain_idx,:) + con_vec(S_dispE_X_pred_toLose_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% to gain
        curr_vec_nm = 'GSpool_mod_dispE_X_pred_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_X_pred_toGain_idx,:) + con_vec(S_dispE_X_pred_toGain_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% to lose
        curr_vec_nm = 'GSpool_mod_dispE_X_pred_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_X_pred_toLose_idx,:) + con_vec(S_dispE_X_pred_toLose_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% gain - loss
        curr_vec_nm = 'GSpool_mod_dispE_X_pred_Gain_min_Loss';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_X_pred_toGain_idx,:) - con_vec(G_dispE_X_pred_toLose_idx,:) +...
            con_vec(S_dispE_X_pred_toGain_idx,:) - con_vec(S_dispE_X_pred_toLose_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        %% Grip>Stroop gain - loss
        curr_vec_nm = 'Gmin_S_mod_dispE_X_pred_Gain_min_Loss';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = (con_vec(G_dispE_X_pred_toGain_idx,:) - con_vec(G_dispE_X_pred_toLose_idx,:)) +...
            -(con_vec(S_dispE_X_pred_toGain_idx,:) - con_vec(S_dispE_X_pred_toLose_idx,:)); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% performance
    if G_mod_dispE.perf ~= 0 && S_mod_dispE.perf ~= 0
        G_dispE_perf_toGain_idx = strcmp(con_names,['G_mod_dispE_perf_toGain',posCon_nm]);
        G_dispE_perf_toLose_idx = strcmp(con_names,['G_mod_dispE_perf_toLose',posCon_nm]);
        S_dispE_perf_toGain_idx = strcmp(con_names,['S_mod_dispE_perf_toGain',posCon_nm]);
        S_dispE_perf_toLose_idx = strcmp(con_names,['S_mod_dispE_perf_toLose',posCon_nm]);
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_dispE_perf';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_perf_toGain_idx,:) + con_vec(G_dispE_perf_toLose_idx,:) +...
            con_vec(S_dispE_perf_toGain_idx,:) + con_vec(S_dispE_perf_toLose_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % to gain
        curr_vec_nm = 'GSpool_mod_dispE_perf_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_perf_toGain_idx,:) + con_vec(S_dispE_perf_toGain_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % to lose
        curr_vec_nm = 'GSpool_mod_dispE_perf_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_perf_toLose_idx,:) + con_vec(S_dispE_perf_toLose_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
elseif G_o_dispE == 1 && S_o_dispE == 3 % pool grip trials with stroop no-error trials
    
    %% pool grip and stroop error AND no error trials
    % run2 - run1
    curr_vec_nm = 'GSpool_o_dispE_run2_min_run1';
    % positive contrast
    jCon = jCon + 1;
    G_oDispE_r2_min_r1_idx = strcmp(con_names,['G_o_dispE_run2_min_run1',posCon_nm]);
    S_oDispE_r2_min_r1_idx = strcmp(con_names,['S_o_dispE_run2_min_run1_ETnoETpool',posCon_nm]);
    con_vec(jCon, 1:n_regs) = con_vec(G_oDispE_r2_min_r1_idx,:) + con_vec(S_oDispE_r2_min_r1_idx,:); % pool incentive to gain and to lose trials
    con_names{jCon} = [curr_vec_nm,posCon_nm];
    % negative contrast
    [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    
    %% GL condition
    if G_mod_dispE.GL_cond ~= 0 && S_mod_dispE.GL_cond ~= 0
        
        curr_vec_nm = 'GSpool_mod_dispE_GL_cond';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_GLcond_idx = strcmp(con_names,['G_mod_dispE_GL_cond',posCon_nm]);
        S_dispE_GLcond_idx = strcmp(con_names,['S_mod_dispE_GL_cond_ETnoETpool',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_GLcond_idx,:) + con_vec(S_dispE_GLcond_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% reward type
    if G_mod_dispE.R_type ~= 0 && S_mod_dispE.R_type ~= 0
        
        curr_vec_nm = 'GSpool_mod_dispE_R_type';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_R_type_idx = strcmp(con_names,['G_mod_dispE_R_type',posCon_nm]);
        S_dispE_R_type_idx = strcmp(con_names,['S_mod_dispE_R_type_ETnoETpool',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_R_type_idx,:) + con_vec(S_dispE_R_type_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% incentive
    if G_mod_dispE.inc ~= 0 && S_mod_dispE.inc ~= 0
        
        curr_vec_nm = 'GSpool_mod_dispE_inc';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_inc_idx = strcmp(con_names,['G_mod_dispE_inc',posCon_nm]);
        S_dispE_inc_idx = strcmp(con_names,['S_mod_dispE_inc_ETnoETpool',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_inc_idx,:) + con_vec(S_dispE_inc_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% incentive bis
    if G_mod_dispE.inc_bis ~= 0 && S_mod_dispE.inc_bis ~= 0
        
        curr_vec_nm = 'GSpool_mod_dispE_inc_bis';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_inc_bis_idx = strcmp(con_names,['G_mod_dispE_inc_bis',posCon_nm]);
        S_dispE_inc_bis_idx = strcmp(con_names,['S_mod_dispE_inc_bis_ETnoETpool',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_inc_bis_idx,:) + con_vec(S_dispE_inc_bis_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % ROI activity
    if G_mod_dispE.ROI_activity_yn ~= 0 && S_mod_dispE.ROI_activity_yn ~= 0
        G_ROI_full_nm = ['mod_dispE_GLM',G_mod_dispE.ROI_activity_GLM,...
            '_',G_mod_dispE.ROI_activity_period,'_period_',G_mod_dispE.ROI_activity_ROI_nm,'_activity'];
        S_ROI_full_nm = ['mod_dispE_GLM',S_mod_dispE.ROI_activity_GLM,...
            '_',S_mod_dispE.ROI_activity_period,'_period_',S_mod_dispE.ROI_activity_ROI_nm,'_activity'];
        if strcmp(G_ROI_full_nm,S_ROI_full_nm)
            ROI_full_nm = G_ROI_full_nm;
        else
            ROI_full_nm = ['G_',G_ROI_full_nm,'_and_S_',S_ROI_full_nm,'_ETnoETpool'];
        end
        
        curr_vec_nm = ['GSpool_',ROI_full_nm];
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_ROIactivity_idx = strcmp(con_names,['G_',G_ROI_full_nm,posCon_nm]);
        S_dispE_ROIactivity_idx = strcmp(con_names,['S_',S_ROI_full_nm,'_ETnoETpool',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_ROIactivity_idx,:) + con_vec(S_dispE_ROIactivity_idx,:); % pool effort to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % RT first press
    if G_mod_dispE.RT_fp ~= 0 && S_mod_dispE.RT_fp ~= 0
        curr_vec_nm = 'GSpool_mod_dispE_RT_fp';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_RTfp_idx = strcmp(con_names,['G_mod_dispE_RT_fp',posCon_nm]);
        S_dispE_RTfp_idx = strcmp(con_names,['S_mod_dispE_RT_fp_ETnoETpool',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_RTfp_idx,:) + con_vec(S_dispE_RTfp_idx,:); % pool RT first press grip + stroop trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % trial number
    if G_mod_dispE.trialN ~= 0 && S_mod_dispE.trialN ~= 0
        curr_vec_nm = 'GSpool_mod_dispE_trialN';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_trialN_idx = strcmp(con_names,['G_mod_dispE_trialN',posCon_nm]);
        S_dispE_trialN_idx = strcmp(con_names,['S_mod_dispE_trialN_ETnoETpool',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_trialN_idx,:) + con_vec(S_dispE_trialN_idx,:); % pool trialN grip + stroop trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % performance
    if G_mod_dispE.perf ~= 0 && S_mod_dispE.perf ~= 0
        curr_vec_nm = 'GSpool_mod_dispE_perf';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_perf_idx = strcmp(con_names,['G_mod_dispE_perf',posCon_nm]);
        S_dispE_perf_idx = strcmp(con_names,['S_mod_dispE_perf_ETnoETpool',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_perf_idx,:) + con_vec(S_dispE_perf_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% pool grip and stroop no error trials only
    % incentive
    if G_mod_dispE.R_type ~= 0 && S_mod_dispE.R_type ~= 0
        
        curr_vec_nm = 'GSnoEpool_mod_dispE_R_type';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_R_type_idx     = strcmp(con_names,['G_mod_dispE_R_type',posCon_nm]);
        S_dispE_R_type_noE_idx = strcmp(con_names,['S_mod_dispE_R_type_noErrorTrials',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_R_type_idx,:) + con_vec(S_dispE_R_type_noE_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % incentive
    if G_mod_dispE.inc ~= 0 && S_mod_dispE.inc ~= 0
        
        curr_vec_nm = 'GSnoEpool_mod_dispE_inc';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_inc_idx     = strcmp(con_names,['G_mod_dispE_inc',posCon_nm]);
        S_dispE_inc_noE_idx = strcmp(con_names,['S_mod_dispE_inc_noErrorTrials',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_inc_idx,:) + con_vec(S_dispE_inc_noE_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % incentive bis
    if G_mod_dispE.inc_bis ~= 0 && S_mod_dispE.inc_bis ~= 0
        
        curr_vec_nm = 'GSnoEpool_mod_dispE_inc_bis';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_inc_bis_idx     = strcmp(con_names,['G_mod_dispE_inc_bis',posCon_nm]);
        S_dispE_inc_bis_noE_idx = strcmp(con_names,['S_mod_dispE_inc_bis_noErrorTrials',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_inc_bis_idx,:) + con_vec(S_dispE_inc_bis_noE_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % ROI activity
    if G_mod_dispE.ROI_activity_yn ~= 0 && S_mod_dispE.ROI_activity_yn ~= 0
        G_ROI_full_nm = ['mod_dispE_GLM',G_mod_dispE.ROI_activity_GLM,...
            '_',G_mod_dispE.ROI_activity_period,'_period_',G_mod_dispE.ROI_activity_ROI_nm,'_activity'];
        S_ROI_full_nm = ['mod_dispE_GLM',S_mod_dispE.ROI_activity_GLM,...
            '_',S_mod_dispE.ROI_activity_period,'_period_',S_mod_dispE.ROI_activity_ROI_nm,'_activity'];
        if strcmp(G_ROI_full_nm,S_ROI_full_nm)
            ROI_full_nm = G_ROI_full_nm;
        else
            ROI_full_nm = ['G_',G_ROI_full_nm,'_and_S_',S_ROI_full_nm,'_noErrorTrials'];
        end
        
        curr_vec_nm = ['GSnoEpool_',ROI_full_nm];
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_ROIactivity_idx = strcmp(con_names,['G_',G_ROI_full_nm,posCon_nm]);
        S_dispE_ROIactivity_idx = strcmp(con_names,['S_',S_ROI_full_nm,'_noErrorTrials',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_ROIactivity_idx,:) + con_vec(S_dispE_ROIactivity_idx,:); % pool effort to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % model variables
    if G_mod_dispE.mdl_n ~= 0 && S_mod_dispE.mdl_n ~= 0
        
        % pool effort to gain and to lose trials
        % cost
        if G_mod_dispE.mdl_cost ~= 0 && S_mod_dispE.mdl_cost ~= 0
            % to gain
            curr_vec_nm = ['GSnoEpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_cost'];
            % positive contrast
            jCon = jCon + 1;
            G_dispE_cost_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_cost',posCon_nm]);
            S_dispE_cost_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost_noErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_cost_idx,:) + con_vec(S_dispE_cost_noErrorTrials_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % benefit
        if G_mod_dispE.mdl_benefit ~= 0 && S_mod_dispE.mdl_benefit ~= 0
            % to gain
            curr_vec_nm = ['GSnoEpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_benefit'];
            % positive contrast
            jCon = jCon + 1;
            G_dispE_benefit_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_benefit',posCon_nm]);
            S_dispE_benefit_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit_noErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_benefit_idx,:) + con_vec(S_dispE_benefit_noErrorTrials_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % EV
        if G_mod_dispE.mdl_EV ~= 0 && S_mod_dispE.mdl_EV ~= 0
            curr_vec_nm = ['GSnoEpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_EV'];
            % positive contrast
            jCon = jCon + 1;
            G_dispE_EV_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_EV',posCon_nm]);
            S_dispE_EV_noErrorTrials_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_EV_noErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_EV_idx,:) + con_vec(S_dispE_EV_noErrorTrials_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % model variables compute EV = benefit - cost
        if G_mod_dispE.mdl_cost ~= 0 && G_mod_dispE.mdl_benefit ~= 0
            % to gain
            curr_vec_nm = ['GSpool_mod_dispE_Gmodel',num2str(G_mod_dispE.mdl_n),'_Smodel',num2str(S_mod_dispE.mdl_n),'_Benefit_min_Cost'];
            % positive contrast
            jCon = jCon + 1;
            G_dispE_mdlBenef_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_benefit',posCon_nm]);
            G_dispE_mdlCost_idx = strcmp(con_names,['G_mod_dispE_model',num2str(G_mod_dispE.mdl_n),'_cost',posCon_nm]);
            S_dispE_noErrorTrials_mdlBenef_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_benefit_noErrorTrials',posCon_nm]);
            S_dispE_noErrorTrials_mdlCost_idx = strcmp(con_names,['S_mod_dispE_model',num2str(S_mod_dispE.mdl_n),'_cost_noErrorTrials',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(G_dispE_mdlBenef_idx,:) - con_vec(G_dispE_mdlCost_idx,:) +...
                con_vec(S_dispE_noErrorTrials_mdlBenef_idx,:) - con_vec(S_dispE_noErrorTrials_mdlCost_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
    end
    
    % RT first press
    if G_mod_dispE.RT_fp ~= 0 && S_mod_dispE.RT_fp ~= 0
        curr_vec_nm = 'GSnoEpool_mod_dispE_RT_fp';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_RTfp_idx        = strcmp(con_names,['G_mod_dispE_RT_fp',posCon_nm]);
        S_dispE_RTfp_noE_idx    = strcmp(con_names,['S_mod_dispE_RT_fp_noErrorTrials',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_RTfp_idx,:) + con_vec(S_dispE_RTfp_noE_idx,:); % pool RT first press grip + stroop trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % trial number
    if G_mod_dispE.trialN ~= 0 && S_mod_dispE.trialN ~= 0
        curr_vec_nm = 'GSnoEpool_mod_dispE_trialN';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_trialN_idx        = strcmp(con_names,['G_mod_dispE_trialN',posCon_nm]);
        S_dispE_trialN_noE_idx    = strcmp(con_names,['S_mod_dispE_trialN_noErrorTrials',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_trialN_idx,:) + con_vec(S_dispE_trialN_noE_idx,:); % pool trialN grip + stroop trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % performance
    if G_mod_dispE.perf ~= 0 && S_mod_dispE.perf ~= 0
        curr_vec_nm = 'GSnoEpool_mod_dispE_perf';
        
        % positive contrast
        jCon = jCon + 1;
        G_dispE_perf_idx = strcmp(con_names,['G_mod_dispE_perf',posCon_nm]);
        S_dispE_perf_noE_idx = strcmp(con_names,['S_mod_dispE_perf_noErrorTrials',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_dispE_perf_idx,:) + con_vec(S_dispE_perf_noE_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
end


%% feedback
G_o_fbk = gripRPprm.o_fbk;
G_mod_fbk = gripRPprm.mod_fbk;
S_o_fbk = stroopRPprm.o_fbk;
S_mod_fbk = stroopRPprm.mod_fbk;
if G_o_fbk == 1 && S_o_fbk == 1 % all trials pooled => pool tasks
    % GL condition
    if G_mod_fbk.GL_cond ~= 0 && S_mod_fbk.GL_cond ~= 0
        
        curr_vec_nm = 'GSpool_mod_fbk_GL_cond';
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_GLcond_idx = strcmp(con_names,['G_mod_fbk_GL_cond',posCon_nm]);
        S_fbk_GLcond_idx = strcmp(con_names,['S_mod_fbk_GL_cond',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_GLcond_idx,:) + con_vec(S_fbk_GLcond_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % feedback
    if G_mod_fbk.fbk ~= 0 && S_mod_fbk.fbk ~= 0
        
        curr_vec_nm = 'GSpool_mod_fbk_fbk';
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_fbk_idx = strcmp(con_names,['G_mod_fbk_fbk',posCon_nm]);
        S_fbk_fbk_idx = strcmp(con_names,['S_mod_fbk_fbk',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_fbk_idx,:) + con_vec(S_fbk_fbk_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % total gain
    if G_mod_fbk.totalGain ~= 0 && S_mod_fbk.totalGain ~= 0
        
        curr_vec_nm = 'GSpool_mod_fbk_totalGain';
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_totalGain_idx = strcmp(con_names,['G_mod_fbk_totalGain',posCon_nm]);
        S_fbk_totalGain_idx = strcmp(con_names,['S_mod_fbk_totalGain',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_totalGain_idx,:) + con_vec(S_fbk_totalGain_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % performance
    if G_mod_fbk.perf ~= 0 && S_mod_fbk.perf ~= 0
        curr_vec_nm = 'GSpool_mod_fbk_perf';
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_perf_idx = strcmp(con_names,['G_mod_fbk_perf',posCon_nm]);
        S_fbk_perf_idx = strcmp(con_names,['S_mod_fbk_perf',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_perf_idx,:) + con_vec(S_fbk_perf_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % ROI activity
    if G_mod_fbk.ROI_activity_yn ~= 0 && S_mod_fbk.ROI_activity_yn ~= 0
        G_ROI_full_nm = ['mod_fbk_GLM',G_mod_fbk.ROI_activity_GLM,...
            '_',G_mod_fbk.ROI_activity_period,'_period_',G_mod_fbk.ROI_activity_ROI_nm,'_activity'];
        S_ROI_full_nm = ['mod_fbk_GLM',S_mod_fbk.ROI_activity_GLM,...
            '_',S_mod_fbk.ROI_activity_period,'_period_',S_mod_fbk.ROI_activity_ROI_nm,'_activity'];
        if strcmp(G_ROI_full_nm,S_ROI_full_nm)
            ROI_full_nm = G_ROI_full_nm;
        else
            ROI_full_nm = ['G_',G_ROI_full_nm,'_and_S_',S_ROI_full_nm];
        end
        curr_vec_nm = ['GSpool_',ROI_full_nm];
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_ROIactivity_idx = strcmp(con_names,['G_',ROI_full_nm,posCon_nm]);
        S_fbk_ROIactivity_idx = strcmp(con_names,['S_',ROI_full_nm,posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_ROIactivity_idx,:) + con_vec(S_fbk_ROIactivity_idx,:); % pool effort to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % trialN
    if G_mod_fbk.trialN ~= 0 && S_mod_fbk.trialN ~= 0
        
        curr_vec_nm = 'GSpool_mod_fbk_trialN';
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_trialN_idx = strcmp(con_names,['G_mod_fbk_trialN',posCon_nm]);
        S_fbk_trialN_idx = strcmp(con_names,['S_mod_fbk_trialN',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_trialN_idx,:) + con_vec(S_fbk_trialN_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
elseif G_o_fbk == 2 && S_o_fbk == 2 % to gain/to loss => pool tasks and conditions
    
    % feedback
    if G_mod_fbk.fbk ~= 0 && S_mod_fbk.fbk ~= 0
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_fbk_fbk';
        % positive contrast
        jCon = jCon + 1;
        G_fbk_fbk_toGain_idx = strcmp(con_names,['G_mod_fbk_fbk_toGain',posCon_nm]);
        G_fbk_fbk_toLose_idx = strcmp(con_names,['G_mod_fbk_fbk_toLose',posCon_nm]);
        S_fbk_fbk_toGain_idx = strcmp(con_names,['S_mod_fbk_fbk_toGain',posCon_nm]);
        S_fbk_fbk_toLose_idx = strcmp(con_names,['S_mod_fbk_fbk_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_fbk_fbk_toGain_idx,:) + con_vec(S_fbk_fbk_toGain_idx,:) +...
            con_vec(G_fbk_fbk_toLose_idx,:) + con_vec(S_fbk_fbk_toLose_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % pool tasks together but keep valence (to gain/to lose trials)
        % to gain
        curr_vec_nm = 'GSpool_mod_fbk_fbk_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_fbk_fbk_toGain_idx,:) + con_vec(S_fbk_fbk_toGain_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % to lose
        curr_vec_nm = 'GSpool_mod_fbk_fbk_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_fbk_fbk_toLose_idx,:) + con_vec(S_fbk_fbk_toLose_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % total gain
    if G_mod_fbk.totalGain ~= 0 && S_mod_fbk.totalGain ~= 0
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_fbk_totalGain';
        % positive contrast
        jCon = jCon + 1;
        G_fbk_totalGain_toGain_idx = strcmp(con_names,['G_mod_fbk_totalGain_toGain',posCon_nm]);
        G_fbk_totalGain_toLose_idx = strcmp(con_names,['G_mod_fbk_totalGain_toLose',posCon_nm]);
        S_fbk_totalGain_toGain_idx = strcmp(con_names,['S_mod_fbk_totalGain_toGain',posCon_nm]);
        S_fbk_totalGain_toLose_idx = strcmp(con_names,['S_mod_fbk_totalGain_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_fbk_totalGain_toGain_idx,:) + con_vec(S_fbk_totalGain_toGain_idx,:) +...
            con_vec(G_fbk_totalGain_toLose_idx,:) + con_vec(S_fbk_totalGain_toLose_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % pool tasks together but keep valence (to gain/to lose trials)
        % to gain
        curr_vec_nm = 'GSpool_mod_fbk_totalGain_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_fbk_totalGain_toGain_idx,:) + con_vec(S_fbk_totalGain_toGain_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % to lose
        curr_vec_nm = 'GSpool_mod_fbk_totalGain_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_fbk_totalGain_toLose_idx,:) + con_vec(S_fbk_totalGain_toLose_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % performance
    if G_mod_fbk.perf ~= 0 && S_mod_fbk.perf ~= 0
        G_fbk_perf_toGain_idx = strcmp(con_names,['G_mod_fbk_perf_toGain',posCon_nm]);
        G_fbk_perf_toLose_idx = strcmp(con_names,['G_mod_fbk_perf_toLose',posCon_nm]);
        S_fbk_perf_toGain_idx = strcmp(con_names,['S_mod_fbk_perf_toGain',posCon_nm]);
        S_fbk_perf_toLose_idx = strcmp(con_names,['S_mod_fbk_perf_toLose',posCon_nm]);
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_fbk_perf';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_perf_toGain_idx,:) + con_vec(G_fbk_perf_toLose_idx,:) +...
            con_vec(S_fbk_perf_toGain_idx,:) + con_vec(S_fbk_perf_toLose_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % to gain
        curr_vec_nm = 'GSpool_mod_fbk_perf_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_perf_toGain_idx,:) + con_vec(S_fbk_perf_toGain_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % to lose
        curr_vec_nm = 'GSpool_mod_fbk_perf_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_perf_toLose_idx,:) + con_vec(S_fbk_perf_toLose_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % ROI activity
    if G_mod_fbk.ROI_activity_yn ~= 0 && S_mod_fbk.ROI_activity_yn ~= 0
        G_ROI_full_nm = ['mod_fbk_GLM',G_mod_fbk.ROI_activity_GLM,...
            '_',G_mod_fbk.ROI_activity_period,'_period_',G_mod_fbk.ROI_activity_ROI_nm,'_activity'];
        S_ROI_full_nm = ['mod_fbk_GLM',S_mod_fbk.ROI_activity_GLM,...
            '_',S_mod_fbk.ROI_activity_period,'_period_',S_mod_fbk.ROI_activity_ROI_nm,'_activity'];
        if strcmp(G_ROI_full_nm,S_ROI_full_nm)
            ROI_full_nm = G_ROI_full_nm;
        else
            ROI_full_nm = ['G_',G_ROI_full_nm,'_and_S_',S_ROI_full_nm];
        end
        
        G_fbk_ROIactivity_toGain_idx = strcmp(con_names,['G_',G_ROI_full_nm,posCon_nm]);
        S_fbk_ROIactivity_toGain_idx = strcmp(con_names,['S_',S_ROI_full_nm,posCon_nm]);
        G_fbk_ROIactivity_toLose_idx = strcmp(con_names,['G_',G_ROI_full_nm,posCon_nm]);
        S_fbk_ROIactivity_toLose_idx = strcmp(con_names,['S_',S_ROI_full_nm,posCon_nm]);
        
        curr_vec_nm = ['GSpool_',ROI_full_nm];
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_ROIactivity_toGain_idx,:) + con_vec(S_fbk_ROIactivity_toGain_idx,:) +...
            con_vec(G_fbk_ROIactivity_toLose_idx,:) + con_vec(S_fbk_ROIactivity_toLose_idx,:); % pool effort to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        curr_vec_nm = ['GSpool_',ROI_full_nm,'_toGain'];
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_ROIactivity_toGain_idx,:) + con_vec(S_fbk_ROIactivity_toGain_idx,:); % pool effort to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        curr_vec_nm = ['GSpool_',ROI_full_nm,'_toLose'];
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_ROIactivity_toLose_idx,:) + con_vec(S_fbk_ROIactivity_toLose_idx,:); % pool effort to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % trial number
    if G_mod_fbk.trialN ~= 0 && S_mod_fbk.trialN ~= 0
        
        % pool all conditions and tasks together
        curr_vec_nm = 'GSpool_mod_fbk_trialN';
        % positive contrast
        jCon = jCon + 1;
        G_fbk_trialN_toGain_idx = strcmp(con_names,['G_mod_fbk_trialN_toGain',posCon_nm]);
        G_fbk_trialN_toLose_idx = strcmp(con_names,['G_mod_fbk_trialN_toLose',posCon_nm]);
        S_fbk_trialN_toGain_idx = strcmp(con_names,['S_mod_fbk_trialN_toGain',posCon_nm]);
        S_fbk_trialN_toLose_idx = strcmp(con_names,['S_mod_fbk_trialN_toLose',posCon_nm]);
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_fbk_trialN_toGain_idx,:) + con_vec(S_fbk_trialN_toGain_idx,:) +...
            con_vec(G_fbk_trialN_toLose_idx,:) + con_vec(S_fbk_trialN_toLose_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % pool tasks together but keep valence (to gain/to lose trials)
        % to gain
        curr_vec_nm = 'GSpool_mod_fbk_trialN_toGain';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_fbk_trialN_toGain_idx,:) + con_vec(S_fbk_trialN_toGain_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % to lose
        curr_vec_nm = 'GSpool_mod_fbk_trialN_toLose';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) =...
            con_vec(G_fbk_trialN_toLose_idx,:) + con_vec(S_fbk_trialN_toLose_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
elseif G_o_fbk == 1 && S_o_fbk == 3 % pool grip trials with stroop no-error trials
    % GL condition
    if G_mod_fbk.GL_cond ~= 0 && S_mod_fbk.GL_cond ~= 0
        
        curr_vec_nm = 'GSpool_mod_fbk_GL_cond';
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_GLcond_idx = strcmp(con_names,['G_mod_fbk_GL_cond',posCon_nm]);
        S_fbk_GLcond_idx = strcmp(con_names,['S_mod_fbk_GL_cond_ETnoETpool',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_GLcond_idx,:) + con_vec(S_fbk_GLcond_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % feedback
    if G_mod_fbk.fbk ~= 0 && S_mod_fbk.fbk ~= 0
        
        curr_vec_nm = 'GSpool_mod_fbk_fbk';
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_fbk_idx = strcmp(con_names,['G_mod_fbk_fbk',posCon_nm]);
        S_fbk_fbk_idx = strcmp(con_names,['S_mod_fbk_fbk_ETnoETpool',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_fbk_idx,:) + con_vec(S_fbk_fbk_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % total gain
    if G_mod_fbk.totalGain ~= 0 && S_mod_fbk.totalGain ~= 0
        
        curr_vec_nm = 'GSpool_mod_fbk_totalGain';
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_totalGain_idx = strcmp(con_names,['G_mod_fbk_totalGain',posCon_nm]);
        S_fbk_totalGain_idx = strcmp(con_names,['S_mod_fbk_totalGain_ETnoETpool',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_totalGain_idx,:) + con_vec(S_fbk_totalGain_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % performance
    if G_mod_fbk.perf ~= 0 && S_mod_fbk.perf ~= 0
        curr_vec_nm = 'GSpool_mod_fbk_perf';
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_perf_idx = strcmp(con_names,['G_mod_fbk_perf',posCon_nm]);
        S_fbk_perf_idx = strcmp(con_names,['S_mod_fbk_perf_ETnoETpool',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_perf_idx,:) + con_vec(S_fbk_perf_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % ROI activity
    if G_mod_fbk.ROI_activity_yn ~= 0 && S_mod_fbk.ROI_activity_yn ~= 0
        G_ROI_full_nm = ['mod_fbk_GLM',G_mod_fbk.ROI_activity_GLM,...
            '_',G_mod_fbk.ROI_activity_period,'_period_',G_mod_fbk.ROI_activity_ROI_nm,'_activity'];
        S_ROI_full_nm = ['mod_fbk_GLM',S_mod_fbk.ROI_activity_GLM,...
            '_',S_mod_fbk.ROI_activity_period,'_period_',S_mod_fbk.ROI_activity_ROI_nm,'_activity'];
        if strcmp(G_ROI_full_nm,S_ROI_full_nm)
            ROI_full_nm = G_ROI_full_nm;
        else
            ROI_full_nm = ['G_',G_ROI_full_nm,'_and_S_',S_ROI_full_nm,'_ETnoETpool'];
        end
        
        curr_vec_nm = ['GSpool_',ROI_full_nm];
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_ROIactivity_idx = strcmp(con_names,['G_',G_ROI_full_nm,posCon_nm]);
        S_fbk_ROIactivity_idx = strcmp(con_names,['S_',S_ROI_full_nm,'_ETnoETpool',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_ROIactivity_idx,:) + con_vec(S_fbk_ROIactivity_idx,:); % pool effort to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % trial number
    if G_mod_fbk.trialN ~= 0 && S_mod_fbk.trialN ~= 0
        
        curr_vec_nm = 'GSpool_mod_fbk_trialN';
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_trialN_idx = strcmp(con_names,['G_mod_fbk_trialN',posCon_nm]);
        S_fbk_trialN_idx = strcmp(con_names,['S_mod_fbk_trialN_ETnoETpool',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_trialN_idx,:) + con_vec(S_fbk_trialN_idx,:); % pool trialN to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    %% no error trials only
    % GL condition
    if G_mod_fbk.GL_cond ~= 0 && S_mod_fbk.GL_cond ~= 0
        
        curr_vec_nm = 'GSpool_mod_fbk_GL_cond';
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_GLcond_idx = strcmp(con_names,['G_mod_fbk_GL_cond',posCon_nm]);
        S_fbk_GLcond_noE_idx = strcmp(con_names,['S_mod_fbk_GL_cond_noErrorTrials',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_GLcond_idx,:) + con_vec(S_fbk_GLcond_noE_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % feedback
    if G_mod_fbk.fbk ~= 0 && S_mod_fbk.fbk ~= 0
        
        curr_vec_nm = 'GSnoEpool_mod_fbk_fbk';
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_fbk_idx = strcmp(con_names,['G_mod_fbk_fbk',posCon_nm]);
        S_fbk_fbk_noE_idx = strcmp(con_names,['S_mod_fbk_fbk_noErrorTrials',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_fbk_idx,:) + con_vec(S_fbk_fbk_noE_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % total gain
    if G_mod_fbk.totalGain ~= 0 && S_mod_fbk.totalGain ~= 0
        
        curr_vec_nm = 'GSnoEpool_mod_fbk_totalGain';
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_totalGain_idx     = strcmp(con_names,['G_mod_fbk_totalGain',posCon_nm]);
        S_fbk_noE_totalGain_idx = strcmp(con_names,['S_mod_fbk_totalGain_noErrorTrials',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_totalGain_idx,:) + con_vec(S_fbk_noE_totalGain_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % performance
    if G_mod_fbk.perf ~= 0 && S_mod_fbk.perf ~= 0
        curr_vec_nm = 'GSpool_mod_fbk_perf';
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_perf_idx = strcmp(con_names,['G_mod_fbk_perf',posCon_nm]);
        S_fbk_perf_noE_idx = strcmp(con_names,['S_mod_fbk_perf_noErrorTrials',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_perf_idx,:) + con_vec(S_fbk_perf_noE_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % ROI activity
    if G_mod_fbk.ROI_activity_yn ~= 0 && S_mod_fbk.ROI_activity_yn ~= 0
        G_ROI_full_nm = ['mod_fbk_GLM',G_mod_fbk.ROI_activity_GLM,...
            '_',G_mod_fbk.ROI_activity_period,'_period_',G_mod_fbk.ROI_activity_ROI_nm,'_activity'];
        S_ROI_full_nm = ['mod_fbk_GLM',S_mod_fbk.ROI_activity_GLM,...
            '_',S_mod_fbk.ROI_activity_period,'_period_',S_mod_fbk.ROI_activity_ROI_nm,'_activity'];
        if strcmp(G_ROI_full_nm,S_ROI_full_nm)
            ROI_full_nm = G_ROI_full_nm;
        else
            ROI_full_nm = ['G_',G_ROI_full_nm,'_and_S_',S_ROI_full_nm,'_noErrorTrials'];
        end
        
        curr_vec_nm = ['GSpool_',ROI_full_nm];
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_ROIactivity_idx = strcmp(con_names,['G_',G_ROI_full_nm,posCon_nm]);
        S_fbk_ROIactivity_noE_idx = strcmp(con_names,['S_',S_ROI_full_nm,'_noErrorTrials',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_ROIactivity_idx,:) + con_vec(S_fbk_ROIactivity_noE_idx,:); % pool effort to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % trial number
    if G_mod_fbk.trialN ~= 0 && S_mod_fbk.trialN ~= 0
        
        curr_vec_nm = 'GSnoEpool_mod_fbk_trialN';
        
        % positive contrast
        jCon = jCon + 1;
        G_fbk_trialN_idx     = strcmp(con_names,['G_mod_fbk_trialN',posCon_nm]);
        S_fbk_noE_trialN_idx = strcmp(con_names,['S_mod_fbk_trialN_noErrorTrials',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(G_fbk_trialN_idx,:) + con_vec(S_fbk_noE_trialN_idx,:); % pool incentive to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
end


end % function