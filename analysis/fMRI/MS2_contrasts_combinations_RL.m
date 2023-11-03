function [ jCon, con_vec, con_names ] = MS2_contrasts_combinations_RL( jCon, con_vec, con_names,...
    posCon_nm, n_regs, RLprm)
% [ jCon, con_vec, con_names ] = MS2_contrasts_combinations_RL( jCon, con_vec, con_names,...
%     posCon_nm, n_regs, RLprm) combine similar contrasts in the
%     RL task.


% stimuli
L_o_stim = RLprm.o_stim;
L_mod_stim = RLprm.mod_stim;
switch L_o_stim
    case 1
        % one contrast per run
        L_oStim_idx = strcmp(con_names,['L_o_stim',posCon_nm]);
        L_oStim_vec = con_vec(L_oStim_idx,:);
        L_oStim_runs_idx = find(L_oStim_vec ~= 0);
        L_oStim_r1_idx = L_oStim_runs_idx(1);
        L_oStim_r2_idx = L_oStim_runs_idx(2);
        L_oStim_r3_idx = L_oStim_runs_idx(3);
        
        % run 1
        curr_vec_nm = 'L_o_stim_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, L_oStim_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'L_o_stim_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, L_oStim_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 3
        curr_vec_nm = 'L_o_stim_run3';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, L_oStim_r3_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
    case 2 % trials split by pair type
        % gain pair only
        % one contrast per run
        L_oStim_gPair_idx = strcmp(con_names,['L_o_stim_gainPair',posCon_nm]);
        L_oStim_gPair_vec = con_vec(L_oStim_gPair_idx,:);
        L_oStim_gPair_runs_idx = find(L_oStim_gPair_vec ~= 0);
        L_oStim_gPair_r1_idx = L_oStim_gPair_runs_idx(1);
        L_oStim_gPair_r2_idx = L_oStim_gPair_runs_idx(2);
        L_oStim_gPair_r3_idx = L_oStim_gPair_runs_idx(3);
        
        % run 1
        curr_vec_nm = 'L_o_stim_gainPair_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, L_oStim_gPair_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'L_o_stim_gainPair_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, L_oStim_gPair_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 3
        curr_vec_nm = 'L_o_stim_gainPair_run3';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, L_oStim_gPair_r3_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % neutral pair only
        % one contrast per run
        L_oStim_nPair_idx = strcmp(con_names,['L_o_stim_ntalPair',posCon_nm]);
        L_oStim_nPair_vec = con_vec(L_oStim_nPair_idx,:);
        L_oStim_nPair_runs_idx = find(L_oStim_nPair_vec ~= 0);
        L_oStim_nPair_r1_idx = L_oStim_nPair_runs_idx(1);
        L_oStim_nPair_r2_idx = L_oStim_nPair_runs_idx(2);
        L_oStim_nPair_r3_idx = L_oStim_nPair_runs_idx(3);
        
        % run 1
        curr_vec_nm = 'L_o_stim_ntalPair_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, L_oStim_nPair_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'L_o_stim_ntalPair_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, L_oStim_nPair_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 3
        curr_vec_nm = 'L_o_stim_ntalPair_run3';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, L_oStim_nPair_r3_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % loss pair only
        % one contrast per run
        L_oStim_lPair_idx = strcmp(con_names,['L_o_stim_lossPair',posCon_nm]);
        L_oStim_lPair_vec = con_vec(L_oStim_lPair_idx,:);
        L_oStim_lPair_runs_idx = find(L_oStim_lPair_vec ~= 0);
        L_oStim_lPair_r1_idx = L_oStim_lPair_runs_idx(1);
        L_oStim_lPair_r2_idx = L_oStim_lPair_runs_idx(2);
        L_oStim_lPair_r3_idx = L_oStim_lPair_runs_idx(3);
        
        % run 1
        curr_vec_nm = 'L_o_stim_lossPair_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, L_oStim_lPair_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'L_o_stim_lossPair_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, L_oStim_lPair_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 3
        curr_vec_nm = 'L_o_stim_lossPair_run3';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, L_oStim_lPair_r3_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % pool all pairs together per run
        % run 1
        curr_vec_nm = 'L_o_stim_run1';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, L_oStim_gPair_r1_idx) = 1;
        con_vec(jCon, L_oStim_nPair_r1_idx) = 1;
        con_vec(jCon, L_oStim_lPair_r1_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 2
        curr_vec_nm = 'L_o_stim_run2';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, L_oStim_gPair_r2_idx) = 1;
        con_vec(jCon, L_oStim_nPair_r2_idx) = 1;
        con_vec(jCon, L_oStim_lPair_r2_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % run 3
        curr_vec_nm = 'L_o_stim_run3';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = zeros(1,n_regs);
        con_vec(jCon, L_oStim_gPair_r3_idx) = 1;
        con_vec(jCon, L_oStim_nPair_r3_idx) = 1;
        con_vec(jCon, L_oStim_lPair_r3_idx) = 1;
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % trial number
        if L_mod_stim.trialN ~= 0
            curr_vec_nm = 'L_mod_stim_trialN';
            % positive contrast
            jCon = jCon + 1;
            L_stim_trialN_gainPair_idx = strcmp(con_names,['L_mod_stim_trialN_gainPair',posCon_nm]);
            L_stim_trialN_lossPair_idx = strcmp(con_names,['L_mod_stim_trialN_lossPair',posCon_nm]);
            L_stim_trialN_ntalPair_idx = strcmp(con_names,['L_mod_stim_trialN_ntalPair',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_trialN_gainPair_idx,:) +...
                con_vec(L_stim_trialN_lossPair_idx,:) + con_vec(L_stim_trialN_ntalPair_idx,:); % pool gain, neutral and loss pair trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % SV
        if L_mod_stim.SV ~= 0
            % SV gain + loss
            curr_vec_nm = 'L_mod_stim_SV';
            % positive contrast
            jCon = jCon + 1;
            L_stim_SV_gainPair_idx = strcmp(con_names,['L_mod_stim_SV_gainPair',posCon_nm]);
            L_stim_SV_lossPair_idx = strcmp(con_names,['L_mod_stim_SV_lossPair',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_SV_gainPair_idx,:) +...
                con_vec(L_stim_SV_lossPair_idx,:); % pool gain, neutral and loss pair trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % SV gain - loss
            curr_vec_nm = 'L_mod_stim_SV_gain_min_loss';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_SV_gainPair_idx,:) +...
                -con_vec(L_stim_SV_lossPair_idx,:); % pool gain, neutral and loss pair trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % dQ
        if L_mod_stim.dQ ~= 0
            curr_vec_nm = 'L_mod_stim_dQ';
            % positive contrast
            jCon = jCon + 1;
            L_stim_dQ_gainPair_idx = strcmp(con_names,['L_mod_stim_dQ_gainPair',posCon_nm]);
            L_stim_dQ_lossPair_idx = strcmp(con_names,['L_mod_stim_dQ_lossPair',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_dQ_gainPair_idx,:) +...
                con_vec(L_stim_dQ_lossPair_idx,:); % pool gain, neutral and loss pair trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % p(choice = best option)
        if L_mod_stim.pBest ~= 0
            curr_vec_nm = 'L_mod_stim_pBest';
            % positive contrast
            jCon = jCon + 1;
            L_stim_pBest_gainPair_idx = strcmp(con_names,['L_mod_stim_pBest_gainPair',posCon_nm]);
            L_stim_pBest_lossPair_idx = strcmp(con_names,['L_mod_stim_pBest_lossPair',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_pBest_gainPair_idx,:) +...
                con_vec(L_stim_pBest_lossPair_idx,:); % pool gain, neutral and loss pair trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ROI activity
        if L_mod_stim.ROI_activity_yn ~= 0
            ROI_full_nm = ['L_mod_stim_GLM',L_mod_stim.ROI_activity_GLM,...
                    '_',L_mod_stim.ROI_activity_period,'_period_',L_mod_stim.ROI_activity_ROI_nm,'_activity'];
            curr_vec_nm = ROI_full_nm;
            
            % positive contrast
            jCon = jCon + 1;
            L_stim_ROIactivity_gainPair_idx = strcmp(con_names,[ROI_full_nm,'_gainPair',posCon_nm]);
            L_stim_ROIactivity_lossPair_idx = strcmp(con_names,[ROI_full_nm,'_lossPair',posCon_nm]);
            L_stim_ROIactivity_ntalPair_idx = strcmp(con_names,[ROI_full_nm,'_ntalPair',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_ROIactivity_gainPair_idx,:) +...
                con_vec(L_stim_ROIactivity_lossPair_idx,:) +...
                con_vec(L_stim_ROIactivity_ntalPair_idx,:); % pool effort to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % RT
        if L_mod_stim.RT ~= 0
            curr_vec_nm = 'L_mod_stim_RT';
            % positive contrast
            jCon = jCon + 1;
            L_stim_RT_gainPair_idx = strcmp(con_names,['L_mod_stim_RT_gainPair',posCon_nm]);
            L_stim_RT_lossPair_idx = strcmp(con_names,['L_mod_stim_RT_lossPair',posCon_nm]);
            L_stim_RT_ntalPair_idx = strcmp(con_names,['L_mod_stim_RT_ntalPair',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_RT_gainPair_idx,:) +...
                con_vec(L_stim_RT_lossPair_idx,:) +...
                con_vec(L_stim_RT_ntalPair_idx,:); % pool gain, neutral and loss pair trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
    case 3 % trials split by pair type Gain and Loss pairs pooled but Neutral Pair apart
        
        % trial number
        if L_mod_stim.trialN ~= 0
            curr_vec_nm = 'L_mod_stim_trialN';
            % positive contrast
            jCon = jCon + 1;
            L_stim_trialN_GL_Pairs_idx = strcmp(con_names,['L_mod_stim_trialN_GL_Pairs',posCon_nm]);
            L_stim_trialN_ntalPair_idx = strcmp(con_names,['L_mod_stim_trialN_ntalPair',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_trialN_GL_Pairs_idx,:) +...
                con_vec(L_stim_trialN_ntalPair_idx,:); % pool gain, neutral and loss pair trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ROI activity
        if L_mod_stim.ROI_activity_yn ~= 0
            ROI_full_nm = ['L_mod_stim_GLM',L_mod_stim.ROI_activity_GLM,...
                '_',L_mod_stim.ROI_activity_period,'_period_',L_mod_stim.ROI_activity_ROI_nm,'_activity'];
            curr_vec_nm = ROI_full_nm;
            
            % positive contrast
            jCon = jCon + 1;
            L_stim_ROIactivity_GL_Pairs_idx = strcmp(con_names,[ROI_full_nm,'_GL_Pairs',posCon_nm]);
            L_stim_ROIactivity_ntalPair_idx = strcmp(con_names,[ROI_full_nm,'_ntalPair',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_ROIactivity_GL_Pairs_idx,:) +...
                con_vec(L_stim_ROIactivity_ntalPair_idx,:); % pool effort to gain and to lose trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % RT
        if L_mod_stim.RT ~= 0
            curr_vec_nm = 'L_mod_stim_RT';
            % positive contrast
            jCon = jCon + 1;
            L_stim_RT_GL_Pairs_idx = strcmp(con_names,['L_mod_stim_RT_GL_Pairs',posCon_nm]);
            L_stim_RT_ntalPair_idx = strcmp(con_names,['L_mod_stim_RT_ntalPair',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_RT_GL_Pairs_idx,:) +...
                con_vec(L_stim_RT_ntalPair_idx,:); % pool gain, neutral and loss pair trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
    case 5 % split by trials type and trial number
        % pool first and second half of trials for all contrasts
        % and also pool pairs
        
        % indexes for onsets
        L_oStim_gainPair_first_idx = strcmp(con_names,['L_o_stim_gainPair_first',posCon_nm]);
        L_oStim_gainPair_last_idx  = strcmp(con_names,['L_o_stim_gainPair_last',posCon_nm]);
        L_oStim_lossPair_first_idx = strcmp(con_names,['L_o_stim_lossPair_first',posCon_nm]);
        L_oStim_lossPair_last_idx  = strcmp(con_names,['L_o_stim_lossPair_last',posCon_nm]);
        L_oStim_ntalPair_first_idx = strcmp(con_names,['L_o_stim_ntalPair_first',posCon_nm]);
        L_oStim_ntalPair_last_idx  = strcmp(con_names,['L_o_stim_ntalPair_last',posCon_nm]);
        
        % Gain Pair > Neutral Pair (first)
        curr_vec_nm = 'L_o_stim_gainPair_min_ntalPair_first';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(L_oStim_gainPair_first_idx,:) +...
            -con_vec(L_oStim_ntalPair_first_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % Gain Pair > Neutral Pair (last)
        curr_vec_nm = 'L_o_stim_gainPair_min_ntalPair_last';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(L_oStim_gainPair_last_idx,:) +...
            -con_vec(L_oStim_ntalPair_last_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        
        % Loss Pair > Neutral Pair (first)
        curr_vec_nm = 'L_o_stim_lossPair_min_ntalPair_first';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(L_oStim_lossPair_first_idx,:) +...
            -con_vec(L_oStim_ntalPair_first_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % Loss Pair > Neutral Pair (last)
        curr_vec_nm = 'L_o_stim_lossPair_min_ntalPair_last';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(L_oStim_lossPair_last_idx,:) +...
            -con_vec(L_oStim_ntalPair_last_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % Gain Pair > Loss Pair (first)
        curr_vec_nm = 'L_o_stim_gainPair_min_lossPair_first';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(L_oStim_gainPair_first_idx,:) +...
            -con_vec(L_oStim_lossPair_first_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % Gain Pair > Loss Pair (last)
        curr_vec_nm = 'L_o_stim_gainPair_min_lossPair_last';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(L_oStim_gainPair_last_idx,:) +...
            -con_vec(L_oStim_lossPair_last_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % trial number
        if L_mod_stim.trialN ~= 0
            L_stim_trialN_gainPair_first_idx = strcmp(con_names,['L_mod_stim_trialN_gainPair_first',posCon_nm]);
            L_stim_trialN_gainPair_last_idx  = strcmp(con_names,['L_mod_stim_trialN_gainPair_last',posCon_nm]);
            L_stim_trialN_lossPair_first_idx = strcmp(con_names,['L_mod_stim_trialN_lossPair_first',posCon_nm]);
            L_stim_trialN_lossPair_last_idx  = strcmp(con_names,['L_mod_stim_trialN_lossPair_last',posCon_nm]);
            L_stim_trialN_ntalPair_first_idx = strcmp(con_names,['L_mod_stim_trialN_ntalPair_first',posCon_nm]);
            L_stim_trialN_ntalPair_last_idx  = strcmp(con_names,['L_mod_stim_trialN_ntalPair_last',posCon_nm]);
            
            % pool first and second half
            % gain
            curr_vec_nm = 'L_mod_stim_trialN_gainPair';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_trialN_gainPair_first_idx,:) +...
                con_vec(L_stim_trialN_gainPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % neutral
            curr_vec_nm = 'L_mod_stim_trialN_neutralPair';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_trialN_ntalPair_first_idx,:) +...
                con_vec(L_stim_trialN_ntalPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % loss
            curr_vec_nm = 'L_mod_stim_trialN_lossPair';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_trialN_lossPair_first_idx,:) +...
                con_vec(L_stim_trialN_lossPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % all pairs
            L_stim_trialN_gainPair_idx = strcmp(con_names,['L_mod_stim_trialN_gainPair',posCon_nm]);
            L_stim_trialN_lossPair_idx = strcmp(con_names,['L_mod_stim_trialN_lossPair',posCon_nm]);
            L_stim_trialN_ntalPair_idx = strcmp(con_names,['L_mod_stim_trialN_ntalPair',posCon_nm]);
            %
            curr_vec_nm = 'L_mod_stim_trialN';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_trialN_gainPair_idx,:) +...
                con_vec(L_stim_trialN_lossPair_idx,:) +...
                con_vec(L_stim_trialN_ntalPair_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % pool pairs first
            curr_vec_nm = 'L_mod_stim_trialN_first';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_trialN_gainPair_first_idx,:) +...
                con_vec(L_stim_trialN_lossPair_first_idx,:) +...
                con_vec(L_stim_trialN_ntalPair_first_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % pool pairs second
            curr_vec_nm = 'L_mod_stim_trialN_last';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_trialN_gainPair_last_idx,:) +...
                con_vec(L_stim_trialN_lossPair_last_idx,:) +...
                con_vec(L_stim_trialN_ntalPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % SV
        if L_mod_stim.SV ~= 0
            L_stim_SV_gainPair_first_idx = strcmp(con_names,['L_mod_stim_SV_gainPair_first',posCon_nm]);
            L_stim_SV_gainPair_last_idx  = strcmp(con_names,['L_mod_stim_SV_gainPair_last',posCon_nm]);
            L_stim_SV_lossPair_first_idx = strcmp(con_names,['L_mod_stim_SV_lossPair_first',posCon_nm]);
            L_stim_SV_lossPair_last_idx  = strcmp(con_names,['L_mod_stim_SV_lossPair_last',posCon_nm]);
            L_stim_SV_ntalPair_first_idx = strcmp(con_names,['L_mod_stim_SV_ntalPair_first',posCon_nm]);
            L_stim_SV_ntalPair_last_idx  = strcmp(con_names,['L_mod_stim_SV_ntalPair_last',posCon_nm]);
            
            % pool first and second half
            % gain
            curr_vec_nm = 'L_mod_stim_SV_gainPair';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_SV_gainPair_first_idx,:) +...
                con_vec(L_stim_SV_gainPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % neutral
            curr_vec_nm = 'L_mod_stim_SV_neutralPair';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_SV_ntalPair_first_idx,:) +...
                con_vec(L_stim_SV_ntalPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % loss
            curr_vec_nm = 'L_mod_stim_SV_lossPair';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_SV_lossPair_first_idx,:) +...
                con_vec(L_stim_SV_lossPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % all pairs
            L_stim_SV_gainPair_idx = strcmp(con_names,['L_mod_stim_SV_gainPair',posCon_nm]);
            L_stim_SV_lossPair_idx = strcmp(con_names,['L_mod_stim_SV_lossPair',posCon_nm]);
            L_stim_SV_ntalPair_idx = strcmp(con_names,['L_mod_stim_SV_ntalPair',posCon_nm]);
            %
            curr_vec_nm = 'L_mod_stim_SV';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_SV_gainPair_idx,:) +...
                con_vec(L_stim_SV_lossPair_idx,:) +...
                con_vec(L_stim_SV_ntalPair_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % pool pairs first
            curr_vec_nm = 'L_mod_stim_SV_first';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_SV_gainPair_first_idx,:) +...
                con_vec(L_stim_SV_lossPair_first_idx,:) +...
                con_vec(L_stim_SV_ntalPair_first_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % pool pairs second
            curr_vec_nm = 'L_mod_stim_SV_last';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_SV_gainPair_last_idx,:) +...
                con_vec(L_stim_SV_lossPair_last_idx,:) +...
                con_vec(L_stim_SV_ntalPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % dQ
        if L_mod_stim.dQ ~= 0
            L_stim_dQ_gainPair_first_idx = strcmp(con_names,['L_mod_stim_dQ_gainPair_first',posCon_nm]);
            L_stim_dQ_gainPair_last_idx  = strcmp(con_names,['L_mod_stim_dQ_gainPair_last',posCon_nm]);
            L_stim_dQ_lossPair_first_idx = strcmp(con_names,['L_mod_stim_dQ_lossPair_first',posCon_nm]);
            L_stim_dQ_lossPair_last_idx  = strcmp(con_names,['L_mod_stim_dQ_lossPair_last',posCon_nm]);
            L_stim_dQ_ntalPair_first_idx = strcmp(con_names,['L_mod_stim_dQ_ntalPair_first',posCon_nm]);
            L_stim_dQ_ntalPair_last_idx  = strcmp(con_names,['L_mod_stim_dQ_ntalPair_last',posCon_nm]);
            
            % pool first and second half
            % gain
            curr_vec_nm = 'L_mod_stim_dQ_gainPair';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_dQ_gainPair_first_idx,:) +...
                con_vec(L_stim_dQ_gainPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % neutral
            curr_vec_nm = 'L_mod_stim_dQ_neutralPair';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_dQ_ntalPair_first_idx,:) +...
                con_vec(L_stim_dQ_ntalPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % loss
            curr_vec_nm = 'L_mod_stim_dQ_lossPair';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_dQ_lossPair_first_idx,:) +...
                con_vec(L_stim_dQ_lossPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % all pairs
            L_stim_dQ_gainPair_idx = strcmp(con_names,['L_mod_stim_dQ_gainPair',posCon_nm]);
            L_stim_dQ_lossPair_idx = strcmp(con_names,['L_mod_stim_dQ_lossPair',posCon_nm]);
            L_stim_dQ_ntalPair_idx = strcmp(con_names,['L_mod_stim_dQ_ntalPair',posCon_nm]);
            %
            curr_vec_nm = 'L_mod_stim_dQ';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_dQ_gainPair_idx,:) +...
                con_vec(L_stim_dQ_lossPair_idx,:) +...
                con_vec(L_stim_dQ_ntalPair_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % pool pairs first
            curr_vec_nm = 'L_mod_stim_dQ_first';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_dQ_gainPair_first_idx,:) +...
                con_vec(L_stim_dQ_lossPair_first_idx,:) +...
                con_vec(L_stim_dQ_ntalPair_first_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % pool pairs second
            curr_vec_nm = 'L_mod_stim_dQ_last';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_dQ_gainPair_last_idx,:) +...
                con_vec(L_stim_dQ_lossPair_last_idx,:) +...
                con_vec(L_stim_dQ_ntalPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % p(choice = best option)
        if L_mod_stim.pBest ~= 0
            L_stim_pBest_gainPair_first_idx = strcmp(con_names,['L_mod_stim_pBest_gainPair_first',posCon_nm]);
            L_stim_pBest_gainPair_last_idx  = strcmp(con_names,['L_mod_stim_pBest_gainPair_last',posCon_nm]);
            L_stim_pBest_lossPair_first_idx = strcmp(con_names,['L_mod_stim_pBest_lossPair_first',posCon_nm]);
            L_stim_pBest_lossPair_last_idx  = strcmp(con_names,['L_mod_stim_pBest_lossPair_last',posCon_nm]);
            L_stim_pBest_ntalPair_first_idx = strcmp(con_names,['L_mod_stim_pBest_ntalPair_first',posCon_nm]);
            L_stim_pBest_ntalPair_last_idx  = strcmp(con_names,['L_mod_stim_pBest_ntalPair_last',posCon_nm]);
            
            % pool first and second half
            % gain
            curr_vec_nm = 'L_mod_stim_pBest_gainPair';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_pBest_gainPair_first_idx,:) +...
                con_vec(L_stim_pBest_gainPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % neutral
            curr_vec_nm = 'L_mod_stim_pBest_neutralPair';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_pBest_ntalPair_first_idx,:) +...
                con_vec(L_stim_pBest_ntalPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % loss
            curr_vec_nm = 'L_mod_stim_pBest_lossPair';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_pBest_lossPair_first_idx,:) +...
                con_vec(L_stim_pBest_lossPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % all pairs
            L_stim_pBest_gainPair_idx = strcmp(con_names,['L_mod_stim_pBest_gainPair',posCon_nm]);
            L_stim_pBest_lossPair_idx = strcmp(con_names,['L_mod_stim_pBest_lossPair',posCon_nm]);
            L_stim_pBest_ntalPair_idx = strcmp(con_names,['L_mod_stim_pBest_ntalPair',posCon_nm]);
            %
            curr_vec_nm = 'L_mod_stim_pBest';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_pBest_gainPair_idx,:) +...
                con_vec(L_stim_pBest_lossPair_idx,:) +...
                con_vec(L_stim_pBest_ntalPair_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % pool pairs first
            curr_vec_nm = 'L_mod_stim_pBest_first';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_pBest_gainPair_first_idx,:) +...
                con_vec(L_stim_pBest_lossPair_first_idx,:) +...
                con_vec(L_stim_pBest_ntalPair_first_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % pool pairs second
            curr_vec_nm = 'L_mod_stim_pBest_last';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_pBest_gainPair_last_idx,:) +...
                con_vec(L_stim_pBest_lossPair_last_idx,:) +...
                con_vec(L_stim_pBest_ntalPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % ROI activity
        if L_mod_stim.ROI_activity_yn ~= 0
            ROI_full_nm = ['L_mod_stim_GLM',L_mod_stim.ROI_activity_GLM,...
                    '_',L_mod_stim.ROI_activity_period,'_period_',L_mod_stim.ROI_activity_ROI_nm,'_activity'];
            
            L_stim_ROIactivity_gainPair_first_idx   = strcmp(con_names,[ROI_full_nm,'_gainPair_first',posCon_nm]);
            L_stim_ROIactivity_gainPair_last_idx    = strcmp(con_names,[ROI_full_nm,'_gainPair_last',posCon_nm]);
            L_stim_ROIactivity_lossPair_first_idx   = strcmp(con_names,[ROI_full_nm,'_lossPair_first',posCon_nm]);
            L_stim_ROIactivity_lossPair_last_idx    = strcmp(con_names,[ROI_full_nm,'_lossPair_last',posCon_nm]);
            L_stim_ROIactivity_ntalPair_first_idx   = strcmp(con_names,[ROI_full_nm,'_ntalPair_first',posCon_nm]);
            L_stim_ROIactivity_ntalPair_last_idx    = strcmp(con_names,[ROI_full_nm,'_ntalPair_last',posCon_nm]);
            
            % pool first and second half
            % gain
            curr_vec_nm = [ROI_full_nm,'_gainPair'];
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_ROIactivity_gainPair_first_idx,:) +...
                con_vec(L_stim_ROIactivity_gainPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % neutral
            curr_vec_nm = [ROI_full_nm,'_neutralPair'];
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_ROIactivity_ntalPair_first_idx,:) +...
                con_vec(L_stim_ROIactivity_ntalPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % loss
            curr_vec_nm = [ROI_full_nm,'_lossPair'];
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_ROIactivity_lossPair_first_idx,:) +...
                con_vec(L_stim_ROIactivity_lossPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % all pairs
            L_stim_ROIactivity_gainPair_idx = strcmp(con_names,[ROI_full_nm,'_gainPair',posCon_nm]);
            L_stim_ROIactivity_lossPair_idx = strcmp(con_names,[ROI_full_nm,'_lossPair',posCon_nm]);
            L_stim_ROIactivity_ntalPair_idx = strcmp(con_names,[ROI_full_nm,'_ntalPair',posCon_nm]);
            %
            curr_vec_nm = ROI_full_nm;
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_ROIactivity_gainPair_idx,:) +...
                con_vec(L_stim_ROIactivity_lossPair_idx,:) +...
                con_vec(L_stim_ROIactivity_ntalPair_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % pool pairs first
            curr_vec_nm = [ROI_full_nm,'_first'];
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_ROIactivity_gainPair_first_idx,:) +...
                con_vec(L_stim_ROIactivity_lossPair_first_idx,:) +...
                con_vec(L_stim_ROIactivity_ntalPair_first_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % pool pairs second
            curr_vec_nm = [ROI_full_nm,'_last'];
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_ROIactivity_gainPair_last_idx,:) +...
                con_vec(L_stim_ROIactivity_lossPair_last_idx,:) +...
                con_vec(L_stim_ROIactivity_ntalPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
        % RT
        if L_mod_stim.RT ~= 0
            L_stim_RT_gainPair_first_idx = strcmp(con_names,['L_mod_stim_RT_gainPair_first',posCon_nm]);
            L_stim_RT_gainPair_last_idx  = strcmp(con_names,['L_mod_stim_RT_gainPair_last',posCon_nm]);
            L_stim_RT_lossPair_first_idx = strcmp(con_names,['L_mod_stim_RT_lossPair_first',posCon_nm]);
            L_stim_RT_lossPair_last_idx  = strcmp(con_names,['L_mod_stim_RT_lossPair_last',posCon_nm]);
            L_stim_RT_ntalPair_first_idx = strcmp(con_names,['L_mod_stim_RT_ntalPair_first',posCon_nm]);
            L_stim_RT_ntalPair_last_idx  = strcmp(con_names,['L_mod_stim_RT_ntalPair_last',posCon_nm]);
            
            % pool first and second half
            % gain
            curr_vec_nm = 'L_mod_stim_RT_gainPair';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_RT_gainPair_first_idx,:) +...
                con_vec(L_stim_RT_gainPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % neutral
            curr_vec_nm = 'L_mod_stim_RT_neutralPair';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_RT_ntalPair_first_idx,:) +...
                con_vec(L_stim_RT_ntalPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % loss
            curr_vec_nm = 'L_mod_stim_RT_lossPair';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_RT_lossPair_first_idx,:) +...
                con_vec(L_stim_RT_lossPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            % all pairs
            L_stim_RT_gainPair_idx = strcmp(con_names,['L_mod_stim_RT_gainPair',posCon_nm]);
            L_stim_RT_lossPair_idx = strcmp(con_names,['L_mod_stim_RT_lossPair',posCon_nm]);
            L_stim_RT_ntalPair_idx = strcmp(con_names,['L_mod_stim_RT_ntalPair',posCon_nm]);
            %
            curr_vec_nm = 'L_mod_stim_RT';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_RT_gainPair_idx,:) +...
                con_vec(L_stim_RT_lossPair_idx,:) +...
                con_vec(L_stim_RT_ntalPair_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % pool pairs first
            curr_vec_nm = 'L_mod_stim_RT_first';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_RT_gainPair_first_idx,:) +...
                con_vec(L_stim_RT_lossPair_first_idx,:) +...
                con_vec(L_stim_RT_ntalPair_first_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % pool pairs second
            curr_vec_nm = 'L_mod_stim_RT_last';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_RT_gainPair_last_idx,:) +...
                con_vec(L_stim_RT_lossPair_last_idx,:) +...
                con_vec(L_stim_RT_ntalPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end
        
    case 6 % pool gain+loss/neutral split first/second half trials
        
        % indexes for onsets
        L_oStim_GL_Pairs_first_idx = strcmp(con_names,['L_o_stim_GL_Pairs_first',posCon_nm]);
        L_oStim_GL_Pairs_last_idx  = strcmp(con_names,['L_o_stim_GL_Pairs_last',posCon_nm]);
        L_oStim_ntalPair_first_idx = strcmp(con_names,['L_o_stim_ntalPair_first',posCon_nm]);
        L_oStim_ntalPair_last_idx  = strcmp(con_names,['L_o_stim_ntalPair_last',posCon_nm]);
        
        % Gain+Loss Pair > Neutral Pair (first)
        curr_vec_nm = 'L_o_stim_GL_Pairs_min_ntalPair_first';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(L_oStim_GL_Pairs_first_idx,:) +...
            -con_vec(L_oStim_ntalPair_first_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % Gain+Loss Pair > Neutral Pair (last)
        curr_vec_nm = 'L_o_stim_GL_Pairs_min_ntalPair_last';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(L_oStim_GL_Pairs_last_idx,:) +...
            -con_vec(L_oStim_ntalPair_last_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % Gain+Loss Pair (first) > Gain+Loss Pair (last)
        curr_vec_nm = 'L_o_stim_GL_Pairs_first_min_last';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(L_oStim_GL_Pairs_first_idx,:) +...
            -con_vec(L_oStim_GL_Pairs_last_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % Neutral Pair (first) > Neutral Pair (last)
        curr_vec_nm = 'L_o_stim_ntalPair_first_min_last';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(L_oStim_ntalPair_first_idx,:) +...
            -con_vec(L_oStim_ntalPair_last_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % all pairs (first) > all pairs (last)
        curr_vec_nm = 'L_o_stim_first_min_last';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(L_oStim_GL_Pairs_first_idx,:) +...
            con_vec(L_oStim_ntalPair_first_idx,:) +...
            -con_vec(L_oStim_GL_Pairs_last_idx,:) +...
            -con_vec(L_oStim_ntalPair_last_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        if L_mod_stim.pBest ~= 0
            L_stim_pBest_GL_Pairs_first_idx = strcmp(con_names,['L_mod_stim_pBest_GL_Pairs_first',posCon_nm]);
            L_stim_pBest_GL_Pairs_last_idx  = strcmp(con_names,['L_mod_stim_pBest_GL_Pairs_last',posCon_nm]);
            
            % pool first and second half
            % gain + loss
            curr_vec_nm = 'L_mod_stim_pBest_GL_Pairs';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_pBest_GL_Pairs_first_idx,:) +...
                con_vec(L_stim_pBest_GL_Pairs_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % compare confidence first vs second half of trials (confidence
            % should be stronger at the end)
            % gain + loss
            curr_vec_nm = 'L_mod_stim_pBest_GL_Pairs_first_min_last';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_pBest_GL_Pairs_first_idx,:) +...
                -con_vec(L_stim_pBest_GL_Pairs_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
        end % confidence
        
        if L_mod_stim.RT ~= 0
            L_stim_RT_GL_Pairs_first_idx = strcmp(con_names,['L_mod_stim_RT_GL_Pairs_first',posCon_nm]);
            L_stim_RT_GL_Pairs_last_idx  = strcmp(con_names,['L_mod_stim_RT_GL_Pairs_last',posCon_nm]);
            L_stim_RT_ntalPair_first_idx = strcmp(con_names,['L_mod_stim_RT_ntalPair_first',posCon_nm]);
            L_stim_RT_ntalPair_last_idx  = strcmp(con_names,['L_mod_stim_RT_ntalPair_last',posCon_nm]);
            
            % pool first and second half
            % gain + loss
            curr_vec_nm = 'L_mod_stim_RT_GL_Pairs';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_RT_GL_Pairs_first_idx,:) +...
                con_vec(L_stim_RT_GL_Pairs_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % neutral
            curr_vec_nm = 'L_mod_stim_RT_ntalPair';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_RT_ntalPair_first_idx,:) +...
                con_vec(L_stim_RT_ntalPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % all pairs
            L_stim_RT_GL_Pairs_idx = strcmp(con_names,['L_mod_stim_RT_GL_Pairs',posCon_nm]);
            L_stim_RT_ntalPair_idx = strcmp(con_names,['L_mod_stim_RT_ntalPair',posCon_nm]);
            %
            curr_vec_nm = 'L_mod_stim_RT';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_RT_GL_Pairs_idx,:) +...
                con_vec(L_stim_RT_ntalPair_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % pool pairs first
            curr_vec_nm = 'L_mod_stim_RT_first';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_RT_GL_Pairs_first_idx,:) +...
                con_vec(L_stim_RT_ntalPair_first_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % pool pairs second
            curr_vec_nm = 'L_mod_stim_RT_last';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_RT_GL_Pairs_last_idx,:) +...
                con_vec(L_stim_RT_ntalPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % compare RT first vs second half of trials (effort beginning
            % learning, not at the end?)
            % gain + loss
            curr_vec_nm = 'L_mod_stim_RT_GL_Pairs_first_min_last';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_RT_GL_Pairs_first_idx,:) +...
                -con_vec(L_stim_RT_GL_Pairs_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
            
            % neutral
            curr_vec_nm = 'L_mod_stim_RT_ntalPair_first_min_last';
            % positive contrast
            jCon = jCon + 1;
            con_vec(jCon, 1:n_regs) = con_vec(L_stim_RT_ntalPair_first_idx,:) +...
                -con_vec(L_stim_RT_ntalPair_last_idx,:);
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        end % RT
end

% choice


% feedback
L_o_fbk = RLprm.o_fbk;
L_mod_fbk = RLprm.mod_fbk;
if L_o_fbk == 2 % split by pair type
    % trial number
    if L_mod_fbk.trialN ~= 0
            curr_vec_nm = 'L_mod_fbk_trialN';
            % positive contrast
            jCon = jCon + 1;
            L_fbk_trialN_gainPair_idx = strcmp(con_names,['L_mod_fbk_trialN_gainPair',posCon_nm]);
            L_fbk_trialN_lossPair_idx = strcmp(con_names,['L_mod_fbk_trialN_lossPair',posCon_nm]);
            L_fbk_trialN_ntalPair_idx = strcmp(con_names,['L_mod_fbk_trialN_ntalPair',posCon_nm]);
            con_vec(jCon, 1:n_regs) = con_vec(L_fbk_trialN_gainPair_idx,:) +...
                con_vec(L_fbk_trialN_lossPair_idx,:) + con_vec(L_fbk_trialN_ntalPair_idx,:); % pool gain, neutral and loss pair trials
            con_names{jCon} = [curr_vec_nm,posCon_nm];
            % negative contrast
            [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
        
    % feedback
    if L_mod_fbk.fbk ~= 0
        curr_vec_nm = 'L_mod_fbk_fbk';
        % positive contrast
        jCon = jCon + 1;
        L_fbk_fbk_gainPair_idx = strcmp(con_names,['L_mod_fbk_fbk_gainPair',posCon_nm]);
        L_fbk_fbk_lossPair_idx = strcmp(con_names,['L_mod_fbk_fbk_lossPair',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(L_fbk_fbk_gainPair_idx,:) +...
            con_vec(L_fbk_fbk_lossPair_idx,:); % pool gain, neutral and loss pair trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % PE
    if L_mod_fbk.PE ~= 0
        curr_vec_nm = 'L_mod_fbk_PE';
        % positive contrast
        jCon = jCon + 1;
        L_fbk_PE_gainPair_idx = strcmp(con_names,['L_mod_fbk_PE_gainPair',posCon_nm]);
        L_fbk_PE_lossPair_idx = strcmp(con_names,['L_mod_fbk_PE_lossPair',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(L_fbk_PE_gainPair_idx,:) +...
            con_vec(L_fbk_PE_lossPair_idx,:); % pool gain, neutral and loss pair trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % PE bis
    if L_mod_fbk.PE_bis ~= 0
        curr_vec_nm = 'L_mod_fbk_PE_bis';
        % positive contrast
        jCon = jCon + 1;
        L_fbk_PE_bis_gainPair_idx = strcmp(con_names,['L_mod_fbk_PE_bis_gainPair',posCon_nm]);
        L_fbk_PE_bis_lossPair_idx = strcmp(con_names,['L_mod_fbk_PE_bis_lossPair',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(L_fbk_PE_bis_gainPair_idx,:) +...
            con_vec(L_fbk_PE_bis_lossPair_idx,:); % pool gain, neutral and loss pair trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % total gain
    if L_mod_fbk.totalGain ~= 0
        curr_vec_nm = 'L_mod_fbk_totalGain';
        % positive contrast
        jCon = jCon + 1;
        L_fbk_totalGain_gainPair_idx = strcmp(con_names,['L_mod_fbk_totalGain_gainPair',posCon_nm]);
        L_fbk_totalGain_lossPair_idx = strcmp(con_names,['L_mod_fbk_totalGain_lossPair',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(L_fbk_totalGain_gainPair_idx,:) +...
            con_vec(L_fbk_totalGain_lossPair_idx,:); % pool gain, neutral and loss pair trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % ROI activity
    if L_mod_fbk.ROI_activity_yn ~= 0
        ROI_full_nm = ['L_mod_fbk_GLM',L_mod_fbk.ROI_activity_GLM,...
            '_',L_mod_fbk.ROI_activity_period,'_period_',L_mod_fbk.ROI_activity_ROI_nm,'_activity'];
        curr_vec_nm = ROI_full_nm;
        
        % positive contrast
        jCon = jCon + 1;
        L_fbk_ROIactivity_gainPair_idx = strcmp(con_names,[ROI_full_nm,'_gainPair',posCon_nm]);
        L_fbk_ROIactivity_lossPair_idx = strcmp(con_names,[ROI_full_nm,'_lossPair',posCon_nm]);
        L_fbk_ROIactivity_ntalPair_idx = strcmp(con_names,[ROI_full_nm,'_ntalPair',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(L_fbk_ROIactivity_gainPair_idx,:) +...
            con_vec(L_fbk_ROIactivity_lossPair_idx,:) +...
            con_vec(L_fbk_ROIactivity_ntalPair_idx,:); % pool effort to gain and to lose trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
elseif L_o_fbk == 6 % regroup & compare first and second half trials for gain+loss
    
    % PE
    if L_mod_fbk.PE ~= 0
        % regroup PE bis first and second half
        curr_vec_nm = 'L_mod_fbk_PE';
        % positive contrast
        jCon = jCon + 1;
        L_fbk_PE_GLPair_first_idx = strcmp(con_names,['L_mod_fbk_PE_GL_Pairs_first',posCon_nm]);
        L_fbk_PE_GLPair_last_idx = strcmp(con_names,['L_mod_fbk_PE_GL_Pairs_last',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(L_fbk_PE_GLPair_first_idx,:) +...
            con_vec(L_fbk_PE_GLPair_last_idx,:); % pool gain, neutral and loss pair trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % compare PE bis first and second half
        curr_vec_nm = 'L_mod_fbk_PE_first_min_last';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(L_fbk_PE_GLPair_first_idx,:) +...
            -con_vec(L_fbk_PE_GLPair_last_idx,:); % pool gain, neutral and loss pair trials
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
    
    % PE bis
    if L_mod_fbk.PE_bis ~= 0
        % regroup PE bis first and second half
        curr_vec_nm = 'L_mod_fbk_PE_bis';
        % positive contrast
        jCon = jCon + 1;
        L_fbk_PE_bis_GLPairs_first_idx = strcmp(con_names,['L_mod_fbk_PE_bis_GL_Pairs_first',posCon_nm]);
        L_fbk_PE_bis_GLPairs_last_idx = strcmp(con_names,['L_mod_fbk_PE_bis_GL_Pairs_last',posCon_nm]);
        con_vec(jCon, 1:n_regs) = con_vec(L_fbk_PE_bis_GLPairs_first_idx,:) +...
            con_vec(L_fbk_PE_bis_GLPairs_last_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
        
        % compare PE bis first and second half
        curr_vec_nm = 'L_mod_fbk_PE_bis_first_min_last';
        % positive contrast
        jCon = jCon + 1;
        con_vec(jCon, 1:n_regs) = con_vec(L_fbk_PE_bis_GLPairs_first_idx,:) +...
            -con_vec(L_fbk_PE_bis_GLPairs_last_idx,:);
        con_names{jCon} = [curr_vec_nm,posCon_nm];
        % negative contrast
        [jCon, con_vec, con_names] = MS2_write_neg_con(jCon, con_vec, con_names, curr_vec_nm, n_regs);
    end
        
end

end % function