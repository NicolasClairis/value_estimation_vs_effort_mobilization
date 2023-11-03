function [ jCon, con_vec, con_names ] = MS2_contrasts_combinations_GSL( jCon, con_vec, con_names,...
    posCon_nm, n_regs, RLprm, gripRPprm, stroopRPprm)
% [ jCon, con_vec, con_names ] = MS2_contrasts_combinations_GSL( jCon, con_vec, con_names,...
%     posCon_nm, n_regs, RLprm, gripRPprm, stroopRPprm) combine similar contrasts in the
%     grip+Stroop+learning task.


%% trial number
% at incentive display
G_o_inc = gripRPprm.o_inc;
G_mod_inc = gripRPprm.mod_inc;
S_o_inc = stroopRPprm.o_inc;
S_mod_inc = stroopRPprm.mod_inc;
L_o_stim = RLprm.o_stim;
L_mod_stim = RLprm.mod_stim;
if G_o_inc ~= 0 && S_o_inc ~= 0 && L_o_stim ~= 0
    if G_mod_inc.trialN ~= 0 && S_mod_inc.trialN ~= 0 && L_mod_stim.trialN ~=0
        
        % main regressor
        curr_vec_nm = 'GSLpool_inc_mod_trialN';
        % positive contrast
        jCon = jCon + 1;
        GSpool_inc_trialN_idx = strcmp(con_names,['GSpool_mod_inc_trialN',posCon_nm]);
        L_stim_trialN_idx = strcmp(con_names,['L_mod_stim_trialN',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(GSpool_inc_trialN_idx,:) +...
            con_vec(L_stim_trialN_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
end

% at display of effort scale
G_o_dispE = gripRPprm.o_dispE;
G_mod_dispE = gripRPprm.mod_dispE;
S_o_dispE = stroopRPprm.o_dispE;
S_mod_dispE = stroopRPprm.mod_dispE;
if G_o_dispE ~= 0 && S_o_dispE ~= 0 && L_o_stim ~= 0
    if G_mod_dispE.trialN ~= 0 && S_mod_dispE.trialN ~= 0 && L_mod_stim.trialN ~= 0
        
        % main regressor
        curr_vec_nm = 'GSLpool_dispE_mod_trialN';
        % positive contrast
        jCon = jCon + 1;
        GSpool_dispE_trialN_idx = strcmp(con_names,['GSpool_mod_dispE_trialN',posCon_nm]);
        L_stim_trialN_idx = strcmp(con_names,['L_mod_stim_trialN',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(GSpool_dispE_trialN_idx,:) +...
            con_vec(L_stim_trialN_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
end

% at performance onset
G_o_perfE = gripRPprm.o_perfE;
G_mod_perfE = gripRPprm.mod_perfE;
S_o_perfE = stroopRPprm.o_perfE;
S_mod_perfE = stroopRPprm.mod_perfE;
L_o_answer = RLprm.o_answer;
% L_mod_answer = RLprm.mod_answer;
if G_o_perfE ~= 0 && S_o_perfE ~= 0 && L_o_answer ~= 0
    if G_mod_perfE.trialN ~= 0 && S_mod_perfE.trialN ~= 0 && L_mod_answer.trialN ~= 0
        
        % main regressor
        curr_vec_nm = 'GSLpool_perfE_mod_trialN';
        % positive contrast
        jCon = jCon + 1;
        GSpool_perfE_trialN_idx = strcmp(con_names,['GSpool_mod_perfE_trialN',posCon_nm]);
        L_answer_trialN_idx = strcmp(con_names,['L_mod_answer_trialN',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(GSpool_perfE_trialN_idx,:) +...
            con_vec(L_answer_trialN_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
end

% at feedback
G_o_fbk = gripRPprm.o_fbk;
G_mod_fbk = gripRPprm.mod_fbk;
S_o_fbk = stroopRPprm.o_fbk;
S_mod_fbk = stroopRPprm.mod_fbk;
L_o_fbk = RLprm.o_fbk;
L_mod_fbk = RLprm.mod_fbk;
if G_o_fbk ~= 0 && S_o_fbk ~= 0 && L_o_fbk ~= 0
    if G_mod_fbk.trialN ~= 0 && S_mod_fbk.trialN ~= 0 && L_mod_fbk.trialN ~= 0
        % main regressor
        curr_vec_nm = 'GSLpool_fbk_mod_trialN';
        % positive contrast
        jCon = jCon + 1;
        GSpool_fbk_trialN_idx = strcmp(con_names,['GSpool_mod_fbk_trialN',posCon_nm]);
        L_fbk_trialN_idx = strcmp(con_names,['L_mod_fbk_trialN',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(GSpool_fbk_trialN_idx,:) +...
            con_vec(L_fbk_trialN_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
end

%% RT
if G_o_dispE ~= 0 && S_o_dispE ~= 0 && L_o_stim ~= 0
    if G_mod_dispE.RT_fp ~= 0 && S_mod_dispE.RT_fp ~= 0 && L_mod_stim.RT ~= 0
        
        % main regressor
        curr_vec_nm = 'GSLpool_mod_RT_fp';
        % positive contrast
        jCon = jCon + 1;
        GSpool_dispE_RTfp_idx   = strcmp(con_names,['GSpool_mod_dispE_RT_fp',posCon_nm]);
        L_stim_RTfp_idx         = strcmp(con_names,['L_mod_stim_RT',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(GSpool_dispE_RTfp_idx,:) +...
            con_vec(L_stim_RTfp_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
end

%% feedback
if L_o_fbk ~= 0 && G_o_fbk ~= 0 && S_o_fbk ~= 0
    
    if G_mod_fbk.fbk ~= 0 && S_mod_fbk.fbk ~= 0 && L_mod_fbk.fbk ~= 0
        
        % main regressor
        curr_vec_nm = 'GSLpool_mod_fbk_fbk';
        % positive contrast
        jCon = jCon + 1;
        GSpool_fbk_fbk_idx = strcmp(con_names,['GSpool_mod_fbk_fbk',posCon_nm]);
        L_fbk_fbk_idx = strcmp(con_names,['L_mod_fbk_fbk',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(GSpool_fbk_fbk_idx,:) +...
            con_vec(L_fbk_fbk_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
end

%% total gain
if L_o_fbk ~= 0 && G_o_fbk ~= 0 && S_o_fbk ~= 0
    if G_mod_fbk.totalGain ~= 0 && S_mod_fbk.totalGain ~= 0 && L_mod_fbk.totalGain ~= 0
        
        % main regressor
        curr_vec_nm = 'GSLpool_mod_fbk_totalGain';
        % positive contrast
        jCon = jCon + 1;
        GSpool_fbk_totalGain_idx = strcmp(con_names,['GSpool_mod_fbk_totalGain',posCon_nm]);
        L_fbk_totalGain_idx = strcmp(con_names,['L_mod_fbk_totalGain',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(GSpool_fbk_totalGain_idx,:) +...
            con_vec(L_fbk_totalGain_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
end

%% ROI activity

% at RL-stim/GS-incentive display
if L_mod_stim.ROI_activity_yn ~= 0 &&...
        G_mod_inc.ROI_activity_yn ~= 0 &&...
        S_mod_inc.ROI_activity_yn ~= 0
    
    L_ROI_full_nm = ['GLM',L_mod_stim.ROI_activity_GLM,...
        '_',L_mod_stim.ROI_activity_period,'_period_',L_mod_stim.ROI_activity_ROI_nm,'_activity'];
    G_ROI_full_nm = ['GLM',G_mod_inc.ROI_activity_GLM,...
        '_',G_mod_inc.ROI_activity_period,'_period_',G_mod_inc.ROI_activity_ROI_nm,'_activity'];
    S_ROI_full_nm = ['GLM',S_mod_inc.ROI_activity_GLM,...
        '_',S_mod_inc.ROI_activity_period,'_period_',S_mod_inc.ROI_activity_ROI_nm,'_activity'];
    if strcmp(G_ROI_full_nm,S_ROI_full_nm) && strcmp(L_ROI_full_nm, G_ROI_full_nm) && strcmp(L_ROI_full_nm, S_ROI_full_nm)
        ROI_full_nm = G_ROI_full_nm;
    else
        ROI_full_nm = ['G_',G_ROI_full_nm,'_and_S_',S_ROI_full_nm,'_and_L_',L_ROI_full_nm];
    end
    
    curr_vec_nm = ['GSLpool_GSinc_Lstim',ROI_full_nm];
    
    % positive contrast
    jCon = jCon + 1;
    G_inc_ROIactivity_idx = strcmp(con_names,['G_mod_inc_',G_ROI_full_nm,posCon_nm]);
    S_inc_ROIactivity_idx = strcmp(con_names,['S_mod_inc_',S_ROI_full_nm,posCon_nm]);
    L_stim_ROIactivity_idx = strcmp(con_names,['L_mod_stim_',L_ROI_full_nm,posCon_nm]);
    con_vec(jCon, 1:n_regs) = con_vec(G_inc_ROIactivity_idx,:) +...
        con_vec(S_inc_ROIactivity_idx,:) +...
        con_vec(L_stim_ROIactivity_idx,:); % pool effort to gain and to lose trials
    con_names{jCon} = [curr_vec_nm,posCon_nm];
    
    % negative contrast
    [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
end

% at RL-stim/GS-dispE display
if L_mod_stim.ROI_activity_yn ~= 0 &&...
        G_mod_dispE.ROI_activity_yn ~= 0 &&...
        S_mod_dispE.ROI_activity_yn ~= 0
    
    L_ROI_full_nm = ['GLM',L_mod_stim.ROI_activity_GLM,...
        '_',L_mod_stim.ROI_activity_period,'_period_',L_mod_stim.ROI_activity_ROI_nm,'_activity'];
    G_ROI_full_nm = ['GLM',G_mod_dispE.ROI_activity_GLM,...
        '_',G_mod_dispE.ROI_activity_period,'_period_',G_mod_dispE.ROI_activity_ROI_nm,'_activity'];
    S_ROI_full_nm = ['GLM',S_mod_dispE.ROI_activity_GLM,...
        '_',S_mod_dispE.ROI_activity_period,'_period_',S_mod_dispE.ROI_activity_ROI_nm,'_activity'];
    if strcmp(G_ROI_full_nm,S_ROI_full_nm) && strcmp(L_ROI_full_nm, G_ROI_full_nm) && strcmp(L_ROI_full_nm, S_ROI_full_nm)
        ROI_full_nm = G_ROI_full_nm;
    else
        ROI_full_nm = ['G_',G_ROI_full_nm,'_and_S_',S_ROI_full_nm,'_and_L_',L_ROI_full_nm];
    end
    
    curr_vec_nm = ['GSLpool_GSdispE_Lstim',ROI_full_nm];
    
    % positive contrast
    jCon = jCon + 1;
    G_dispE_ROIactivity_idx = strcmp(con_names,['G_mod_dispE_',G_ROI_full_nm,posCon_nm]);
    S_dispE_ROIactivity_idx = strcmp(con_names,['S_mod_dispE_',S_ROI_full_nm,posCon_nm]);
    L_stim_ROIactivity_idx = strcmp(con_names,['L_mod_stim_',L_ROI_full_nm,posCon_nm]);
    con_vec(jCon, 1:n_regs) = con_vec(G_dispE_ROIactivity_idx,:) +...
        con_vec(S_dispE_ROIactivity_idx,:) +...
        con_vec(L_stim_ROIactivity_idx,:); % pool effort to gain and to lose trials
    con_names{jCon} = [curr_vec_nm,posCon_nm];
    
    % negative contrast
    [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
end

% at feedback
if L_mod_fbk.ROI_activity_yn ~= 0 &&...
        G_mod_fbk.ROI_activity_yn ~= 0 &&...
        S_mod_fbk.ROI_activity_yn ~= 0
    
    L_ROI_full_nm = ['mod_fbk_GLM',L_mod_stim.ROI_activity_GLM,...
        '_',L_mod_stim.ROI_activity_period,'_period_',L_mod_stim.ROI_activity_ROI_nm,'_activity'];
    G_ROI_full_nm = ['mod_fbk_GLM',G_mod_fbk.ROI_activity_GLM,...
        '_',G_mod_fbk.ROI_activity_period,'_period_',G_mod_fbk.ROI_activity_ROI_nm,'_activity'];
    S_ROI_full_nm = ['mod_fbk_GLM',S_mod_fbk.ROI_activity_GLM,...
        '_',S_mod_fbk.ROI_activity_period,'_period_',S_mod_fbk.ROI_activity_ROI_nm,'_activity'];
    if strcmp(G_ROI_full_nm,S_ROI_full_nm) &&...
            strcmp(L_ROI_full_nm, G_ROI_full_nm) &&...
            strcmp(L_ROI_full_nm, S_ROI_full_nm)
        ROI_full_nm = G_ROI_full_nm;
    else
        ROI_full_nm = ['G_',G_ROI_full_nm,'_and_S_',S_ROI_full_nm,'_and_L_',L_ROI_full_nm];
    end
    
    curr_vec_nm = ['GSLpool_fbk_',ROI_full_nm];
    
    % positive contrast
    jCon = jCon + 1;
    G_fbk_ROIactivity_idx = strcmp(con_names,['G_',G_ROI_full_nm,posCon_nm]);
    S_fbk_ROIactivity_idx = strcmp(con_names,['S_',S_ROI_full_nm,posCon_nm]);
    L_fbk_ROIactivity_idx = strcmp(con_names,['L_',L_ROI_full_nm,posCon_nm]);
    con_vec(jCon, 1:n_regs) = con_vec(G_fbk_ROIactivity_idx,:) +...
        con_vec(S_fbk_ROIactivity_idx,:) +...
        con_vec(L_fbk_ROIactivity_idx,:); % pool effort to gain and to lose trials
    con_names{jCon} = [curr_vec_nm,posCon_nm];
    
    % negative contrast
    [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
end

end % function

