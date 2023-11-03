function [  ] = MS2_GS_perf_fig( model_EV_n, n_T_bins )
%[  ] = MS2_GS_perf_fig( model_EV_n, n_T_bins )
%MS2_GS_perf_fig will display figure of performance according to
% incentive value for gains and losses and according to trial number
%
% INPUTS
% model_EV_n: number of the expected value model you want to use for the fit
% (if set to 0 or empty, won't be used at all)
%
% n_T_bins: number bins for trial number (5 by default since trials divided
% into 5 blocks with all levels of incentive in each block)
%
% OUTPUTS

%% working directories
root = 'enter path here';
save_dir = [root, filesep, 'behavior_summary' filesep];

%% subject identification
subject_id = {'enter list of subjects here'};
NS = length(subject_id);

%% load variables of interest
is_EV_mdl_used = ~isempty(model_EV_n) && model_EV_n > 0;
if is_EV_mdl_used
    model_EV_nm = ['model_',num2str(model_EV_n)];
else
    model_EV_nm = '';
end
n_GS_runs = 2;
task_names = {'grip','stroop'};
nTasks = length(task_names);
inc_RorV = 'rank'; % 'rank'/'value'
switch inc_RorV
    case 'rank'
        mInc_levels = 1:6;
    case 'value'
        mInc_levels = [0.01 0.2 0.5 1 5 20];
end
n_mInc = length(mInc_levels);
nTrials_per_run = 60;

if ~exist('n_T_bins','var') || isempty(n_T_bins)
    n_T_bins = 5;
end

lossColor = [251 106 74]./255;
gainColor = [165 15 21]./255;
pSize = 60;
lSize = 5;
mSize = 30;
dSpace = 0.05;
overwrite_allowed = 1; % remove previous versions of the figure when launched

%% analyze
for iGS = 1:nTasks
    task_nm = task_names{iGS};
    
    [perf_pred_EV.(task_nm).all,...
        perf_pred_GLM.(task_nm).all,...
        perf.(task_nm).all,...
        absInc.(task_nm).all,...
        GL.(task_nm).all,...
        trialN.(task_nm).all,...
        perf_pred_EV.(task_nm).gainTrials,...
        perf_pred_GLM.(task_nm).gainTrials,...
        perf.(task_nm).gainTrials,...
        absInc.(task_nm).gainTrials,...
        trialN.(task_nm).gainTrials,...
        perf_pred_EV.(task_nm).lossTrials,...
        perf_pred_GLM.(task_nm).lossTrials,...
        perf.(task_nm).lossTrials,...
        absInc.(task_nm).lossTrials,...
        trialN.(task_nm).lossTrials] = deal( NaN(nTrials_per_run*n_GS_runs, NS) );
    [perf_pred_EV_bin.(task_nm).f_inc.all,...
        perf_pred_EV_bin.(task_nm).f_inc.gainTrials,...
        perf_pred_EV_bin.(task_nm).f_inc.lossTrials,...
        perf_pred_GLM_bin.(task_nm).f_inc.all,...
        perf_pred_GLM_bin.(task_nm).f_inc.gainTrials,...
        perf_pred_GLM_bin.(task_nm).f_inc.lossTrials,...
        perf_bin.(task_nm).f_inc.all,...
        perf_bin.(task_nm).f_inc.gainTrials,...
        perf_bin.(task_nm).f_inc.lossTrials] = deal( NaN(n_mInc, NS) );
    [perf_pred_EV_bin.(task_nm).f_trialN.all,...
        perf_pred_EV_bin.(task_nm).f_trialN.gainTrials,...
        perf_pred_EV_bin.(task_nm).f_trialN.lossTrials,...
        perf_pred_GLM_bin.(task_nm).f_trialN.all,...
        perf_pred_GLM_bin.(task_nm).f_trialN.gainTrials,...
        perf_pred_GLM_bin.(task_nm).f_trialN.lossTrials,...
        perf_bin.(task_nm).f_trialN.all,...
        perf_bin.(task_nm).f_trialN.gainTrials,...
        perf_bin.(task_nm).f_trialN.lossTrials,...
        trialN_bin.(task_nm).f_trialN.all,...
        trialN_bin.(task_nm).f_trialN.gainTrials,...
        trialN_bin.(task_nm).f_trialN.lossTrials] = deal( NaN(n_T_bins, NS) );
    
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
        
        % loop through runs
        for iRun = 1:n_GS_runs
            jRun = runs_idx(iRun);
            run_nm = num2str(jRun);
            
            % load data
            loadStruct = getfield( load([onsets_folder,'onsets_sub',subid,'_',task_nm,'_run',run_nm,'.mat'],task_nm),task_nm);
            trialN_tmp = loadStruct.mod.all.trialN;
            trialN_tmp_bis = trialN_tmp + nTrials_per_run*(iRun - 1);
            % split gain/loss
            trialN_gain_tmp = loadStruct.mod.toGain.trialN; % trial number corresponding to the gain trials
            trialN_loss_tmp = loadStruct.mod.toLose.trialN;
            trialN_gain_tmp_bis = trialN_gain_tmp + nTrials_per_run*(iRun - 1);
            trialN_loss_tmp_bis = trialN_loss_tmp + nTrials_per_run*(iRun - 1);
            
            % trial number
            trialN.(task_nm).all(trialN_tmp_bis, iS) = trialN_tmp;
            trialN.(task_nm).gainTrials(trialN_gain_tmp_bis, iS) = trialN_gain_tmp;
            trialN.(task_nm).lossTrials(trialN_loss_tmp_bis, iS) = trialN_loss_tmp;
            
            % gain/loss trials indexes
            gain_trials = loadStruct.mod.all.trialValence == 1;
            GL.(task_nm).all(trialN_tmp_bis, iS) = gain_trials;
            
            % extract incentive
            switch inc_RorV
                case 'rank'
                    absInc_tmp = loadStruct.mod.all.absIncentiveRank;
                    absInc_G_tmp = loadStruct.mod.toGain.absIncentiveRank;
                    absInc_L_tmp = loadStruct.mod.toLose.absIncentiveRank;
                case 'value'
                    absInc_tmp = loadStruct.mod.all.absIncentive;
                    absInc_G_tmp = loadStruct.mod.toGain.absIncentive;
                    absInc_L_tmp = loadStruct.mod.toLose.absIncentive;
            end
            absInc.(task_nm).all(trialN_tmp_bis, iS) = absInc_tmp;
            absInc.(task_nm).gainTrials(trialN_gain_tmp_bis, iS) = absInc_G_tmp;
            absInc.(task_nm).lossTrials(trialN_loss_tmp_bis, iS) = absInc_L_tmp;
            
            % extract perf
            perf_tmp = loadStruct.mod.all.perf./100;
            perf_G_tmp = loadStruct.mod.toGain.perf./100;
            perf_L_tmp = loadStruct.mod.toLose.perf./100;
            perf.(task_nm).all(trialN_tmp_bis, iS) = perf_tmp;
            perf.(task_nm).gainTrials(trialN_gain_tmp_bis, iS) = perf_G_tmp;
            perf.(task_nm).lossTrials(trialN_loss_tmp_bis, iS) = perf_L_tmp;
            
            %% load fitted performance EV model
            if is_EV_mdl_used
                loadFit = load([onsets_folder,'GS_model_sub',subid,'_',task_nm,'_run',run_nm,'.mat']);
                perf_pred_EV_tmp = loadFit.perf_pred.(model_EV_nm);
                perf_pred_EV.(task_nm).all(trialN_tmp_bis, iS) = perf_pred_EV_tmp;
                perf_pred_EV.(task_nm).gainTrials(trialN_gain_tmp_bis, iS) = perf_pred_EV_tmp(gain_trials);
                perf_pred_EV.(task_nm).lossTrials(trialN_loss_tmp_bis, iS) = perf_pred_EV_tmp(gain_trials == 0);
            end
            
        end % run loop
        
        %% perform the bins
        for iInc = 1:n_mInc
            curr_I_trials = absInc.(task_nm).all(:,iS) == mInc_levels(iInc);
            if is_EV_mdl_used
                perf_pred_EV_bin.(task_nm).f_inc.all(iInc,iS) = mean(perf_pred_EV.(task_nm).all(curr_I_trials,iS),'omitnan');
            end
            perf_bin.(task_nm).f_inc.all(iInc,iS) = mean(perf.(task_nm).all(curr_I_trials,iS),'omitnan');
            % split gain/loss
            curr_I_gainTrials = absInc.(task_nm).gainTrials(:,iS) == mInc_levels(iInc);
            curr_I_lossTrials = absInc.(task_nm).lossTrials(:,iS) == mInc_levels(iInc);
            if is_EV_mdl_used
                perf_pred_EV_bin.(task_nm).f_inc.gainTrials(iInc,iS)   = mean(perf_pred_EV.(task_nm).gainTrials(curr_I_gainTrials,iS),'omitnan');
                perf_pred_EV_bin.(task_nm).f_inc.lossTrials(iInc,iS)   = mean(perf_pred_EV.(task_nm).lossTrials(curr_I_lossTrials,iS),'omitnan');
            end
            perf_bin.(task_nm).f_inc.gainTrials(iInc,iS)        = mean(perf.(task_nm).gainTrials(curr_I_gainTrials,iS),'omitnan');
            perf_bin.(task_nm).f_inc.lossTrials(iInc,iS)        = mean(perf.(task_nm).lossTrials(curr_I_lossTrials,iS),'omitnan');
        end
        
        if n_T_bins == 5 % split into 5 blocks per run since each run is subdivided into 5 blocks per run
            for iTrialBlock = 1:n_T_bins
                trial_idx = (1:12) + 12*(iTrialBlock - 1);
                trial_idx_bis = ismember(trialN.(task_nm).all(:,iS), trial_idx);
                if is_EV_mdl_used
                    perf_pred_EV_bin.(task_nm).f_trialN.all(iTrialBlock,iS) = mean(perf_pred_EV.(task_nm).all(trial_idx_bis,iS),'omitnan');
                end
                perf_bin.(task_nm).f_trialN.all(iTrialBlock,iS)         = mean(perf.(task_nm).all(trial_idx_bis,iS),'omitnan');
                trialN_bin.(task_nm).f_trialN.all(iTrialBlock,iS)       = mean(trialN.(task_nm).all(trial_idx_bis,iS),'omitnan');
                
                % gain/loss split
                trial_G_idx = ismember(trialN.(task_nm).gainTrials(:,iS), trial_idx);
                trial_L_idx = ismember(trialN.(task_nm).lossTrials(:,iS), trial_idx);
                %gain
                if is_EV_mdl_used
                    perf_pred_EV_bin.(task_nm).f_trialN.gainTrials(iTrialBlock,iS) = mean(perf_pred_EV.(task_nm).gainTrials(trial_G_idx,iS),'omitnan');
                end
                perf_bin.(task_nm).f_trialN.gainTrials(iTrialBlock,iS)      = mean(perf.(task_nm).gainTrials(trial_G_idx,iS),'omitnan');
                trialN_bin.(task_nm).f_trialN.gainTrials(iTrialBlock,iS)    = mean(trialN.(task_nm).gainTrials(trial_G_idx,iS),'omitnan');
                %loss
                if is_EV_mdl_used
                    perf_pred_EV_bin.(task_nm).f_trialN.lossTrials(iTrialBlock,iS) = mean(perf_pred_EV.(task_nm).lossTrials(trial_L_idx,iS),'omitnan');
                end
                perf_bin.(task_nm).f_trialN.lossTrials(iTrialBlock,iS)      = mean(perf.(task_nm).lossTrials(trial_L_idx,iS),'omitnan');
                trialN_bin.(task_nm).f_trialN.lossTrials(iTrialBlock,iS)    = mean(trialN.(task_nm).lossTrials(trial_L_idx,iS),'omitnan');
            end
        else % if you want another number of bins, ignore the block structure of the runs
            if is_EV_mdl_used
                [perf_pred_EV_bin.(task_nm).f_trialN.all(:,iS),...
                    trialN_bin.(task_nm).f_trialN.all(:,iS)] = do_bin2(perf_pred_EV.(task_nm).all(:,iS),...
                    trialN.(task_nm).all(:,iS),...
                    n_T_bins, 0);
            end
            perf_bin.(task_nm).f_trialN.all(:,iS) = do_bin2(perf.(task_nm).all(:,iS),...
                trialN.(task_nm).all(:,iS),...
                n_T_bins, 0);
            % split gain/loss
            if is_EV_mdl_used
                [perf_pred_EV_bin.(task_nm).f_trialN.gainTrials(:,iS),...
                    trialN_bin.(task_nm).f_trialN.gainTrials(:,iS)] = do_bin2(perf_pred_EV.(task_nm).gainTrials(:,iS),...
                    trialN.(task_nm).gainTrials(:,iS),...
                    n_T_bins, 0);
                [perf_pred_EV_bin.(task_nm).f_trialN.lossTrials(:,iS),...
                    trialN_bin.(task_nm).f_trialN.lossTrials(:,iS)] = do_bin2(perf_pred_EV.(task_nm).lossTrials(:,iS),...
                    trialN.(task_nm).lossTrials(:,iS),...
                    n_T_bins, 0);
            end
            perf_bin.(task_nm).f_trialN.gainTrials(:,iS) = do_bin2(perf.(task_nm).gainTrials(:,iS),...
                trialN.(task_nm).gainTrials(:,iS),...
                n_T_bins, 0);
            perf_bin.(task_nm).f_trialN.lossTrials(:,iS) = do_bin2(perf.(task_nm).lossTrials(:,iS),...
                trialN.(task_nm).lossTrials(:,iS),...
                n_T_bins, 0);
        end
    end % subject loop
    
    %% average bins across subjects
    % mean
    if is_EV_mdl_used
        m_perf_pred_EV_bin.(task_nm).f_inc.all             = mean(perf_pred_EV_bin.(task_nm).f_inc.all, 2,'omitnan');
        m_perf_pred_EV_bin.(task_nm).f_inc.gainTrials      = mean(perf_pred_EV_bin.(task_nm).f_inc.gainTrials, 2,'omitnan');
        m_perf_pred_EV_bin.(task_nm).f_inc.lossTrials      = mean(perf_pred_EV_bin.(task_nm).f_inc.lossTrials, 2,'omitnan');
    end
    m_perf_bin.(task_nm).f_inc.all                  = mean(perf_bin.(task_nm).f_inc.all, 2,'omitnan');
    m_perf_bin.(task_nm).f_inc.gainTrials           = mean(perf_bin.(task_nm).f_inc.gainTrials, 2,'omitnan');
    m_perf_bin.(task_nm).f_inc.lossTrials           = mean(perf_bin.(task_nm).f_inc.lossTrials, 2,'omitnan');
    if is_EV_mdl_used
        m_perf_pred_EV_bin.(task_nm).f_trialN.all          = mean(perf_pred_EV_bin.(task_nm).f_trialN.all, 2,'omitnan');
        m_perf_pred_EV_bin.(task_nm).f_trialN.gainTrials   = mean(perf_pred_EV_bin.(task_nm).f_trialN.gainTrials, 2,'omitnan');
        m_perf_pred_EV_bin.(task_nm).f_trialN.lossTrials   = mean(perf_pred_EV_bin.(task_nm).f_trialN.lossTrials, 2,'omitnan');
    end
    m_perf_bin.(task_nm).f_trialN.all               = mean(perf_bin.(task_nm).f_trialN.all, 2,'omitnan');
    m_perf_bin.(task_nm).f_trialN.gainTrials        = mean(perf_bin.(task_nm).f_trialN.gainTrials, 2,'omitnan');
    m_perf_bin.(task_nm).f_trialN.lossTrials        = mean(perf_bin.(task_nm).f_trialN.lossTrials, 2,'omitnan');
    m_trialN_bin.(task_nm).f_trialN.all             = mean(trialN_bin.(task_nm).f_trialN.all, 2,'omitnan');
    m_trialN_bin.(task_nm).f_trialN.gainTrials      = mean(trialN_bin.(task_nm).f_trialN.gainTrials, 2,'omitnan');
    m_trialN_bin.(task_nm).f_trialN.lossTrials      = mean(trialN_bin.(task_nm).f_trialN.lossTrials, 2,'omitnan');
    
    % SEM
    if is_EV_mdl_used
        sem_perf_pred_EV_bin.(task_nm).f_inc.all           = sem(perf_pred_EV_bin.(task_nm).f_inc.all, 2);
        sem_perf_pred_EV_bin.(task_nm).f_inc.gainTrials    = sem(perf_pred_EV_bin.(task_nm).f_inc.gainTrials, 2);
        sem_perf_pred_EV_bin.(task_nm).f_inc.lossTrials    = sem(perf_pred_EV_bin.(task_nm).f_inc.lossTrials, 2);
    end
    sem_perf_bin.(task_nm).f_inc.all                = sem(perf_bin.(task_nm).f_inc.all, 2);
    sem_perf_bin.(task_nm).f_inc.gainTrials         = sem(perf_bin.(task_nm).f_inc.gainTrials, 2);
    sem_perf_bin.(task_nm).f_inc.lossTrials         = sem(perf_bin.(task_nm).f_inc.lossTrials, 2);
    if is_EV_mdl_used
        sem_perf_pred_EV_bin.(task_nm).f_trialN.all        = sem(perf_pred_EV_bin.(task_nm).f_trialN.all, 2);
        sem_perf_pred_EV_bin.(task_nm).f_trialN.gainTrials = sem(perf_pred_EV_bin.(task_nm).f_trialN.gainTrials, 2);
        sem_perf_pred_EV_bin.(task_nm).f_trialN.lossTrials = sem(perf_pred_EV_bin.(task_nm).f_trialN.lossTrials, 2);
    end
    sem_perf_bin.(task_nm).f_trialN.all             = sem(perf_bin.(task_nm).f_trialN.all, 2);
    sem_perf_bin.(task_nm).f_trialN.gainTrials      = sem(perf_bin.(task_nm).f_trialN.gainTrials, 2);
    sem_perf_bin.(task_nm).f_trialN.lossTrials      = sem(perf_bin.(task_nm).f_trialN.lossTrials, 2);
    sem_trialN_bin.(task_nm).f_trialN.all           = sem(trialN_bin.(task_nm).f_trialN.all, 2);
    sem_trialN_bin.(task_nm).f_trialN.gainTrials    = sem(trialN_bin.(task_nm).f_trialN.gainTrials, 2);
    sem_trialN_bin.(task_nm).f_trialN.lossTrials    = sem(trialN_bin.(task_nm).f_trialN.lossTrials, 2);
    %% graphs
    dirPath = [save_dir, task_nm];
    perf_col            = lossColor;
    perf_pred_EV_col       = [1 0.6 0];
    gain_perf_col       = gainColor;
    loss_perf_col        = lossColor;
    gain_perf_pred_EV_col  = gainColor;
    loss_perf_pred_EV_col  = lossColor;
    %% performance = f(|incentive|)
    fig_perf_hdl = fig();
    perf_hdl = errorbar(mInc_levels,...
        100.*m_perf_bin.(task_nm).f_inc.all,...
        100.*sem_perf_bin.(task_nm).f_inc.all,...
        'r *',...
        'LineWidth',5);
    perf_hdl.Color = perf_col;
    perf_hdl.MarkerSize = mSize;
    hold on;
    if is_EV_mdl_used
        perf_pred_EV_hdl = plot(mInc_levels,...
            100.*m_perf_pred_EV_bin.(task_nm).f_inc.all,...
            '--',...
            'LineWidth',lSize);
        perf_pred_EV_hdl.Color = perf_pred_EV_col;
    end
    xlabel('Incentive (€)');
    ylabel('Performance (%)');
    switch inc_RorV
        case 'rank'
            xlim([0.8 6.2]);
        case 'value'
            xlim([0 20.5]);
    end
    if is_EV_mdl_used
        legend([perf_hdl, perf_pred_EV_hdl,],...
            'Performance','Predicted performance');
    else
        legend(perf_hdl,...
            'Performance');
    end
    legend('Location','northwest');
    legend('boxoff');
    legend_size(pSize);
    img_nm = [task_nm,'_perf_and_',model_EV_nm,'_',model_GLM_nm,'_perf_f_absInc_',N_subs,'subs.png'];
    save_fig(fig_perf_hdl, dirPath, img_nm, overwrite_allowed);
    
    %% performance = f(|incentive| split gain/loss)
    fig_perf_GL_hdl = fig();
    perf_G_hdl = errorbar(mInc_levels-dSpace,...
        100.*m_perf_bin.(task_nm).f_inc.gainTrials,...
        100.*sem_perf_bin.(task_nm).f_inc.gainTrials,...
        'g *',...
        'LineWidth',lSize);
    perf_G_hdl.Color = gain_perf_col;
    perf_G_hdl.Marker = 'o';
    perf_G_hdl.MarkerSize = mSize;
    hold on;
    perf_L_hdl = errorbar(mInc_levels+dSpace,...
        100.*m_perf_bin.(task_nm).f_inc.lossTrials,...
        100.*sem_perf_bin.(task_nm).f_inc.lossTrials,...
        'k *',...
        'LineWidth',lSize);
    perf_L_hdl.Color = loss_perf_col;
    perf_L_hdl.Marker = 's';
    perf_L_hdl.MarkerSize = mSize;
    if is_EV_mdl_used
        perf_pred_EV_G_hdl = plot(mInc_levels,...
            100.*m_perf_pred_EV_bin.(task_nm).f_inc.gainTrials,...
            '--',...
            'LineWidth',lSize);
        perf_pred_EV_G_hdl.Color = gain_perf_pred_EV_col;
        perf_pred_EV_L_hdl = plot(mInc_levels,...
            100.*m_perf_pred_EV_bin.(task_nm).f_inc.lossTrials,...
            '--',...
            'LineWidth',lSize);
        perf_pred_EV_L_hdl.Color = loss_perf_pred_EV_col;
    end
    xticks(1:6);
    xticklabels({'0.01','0.20','0.50','1.00','5.00','20.00'});
    xlabel('Incentive (€)');
    ylabel('Performance (%)');
    switch inc_RorV
        case 'rank'
            xlim([0.8 6.2]);
        case 'value'
            xlim([0 20.5]);
    end
    legend_size(pSize);
    img_nm = [task_nm,'_perf_and_',model_EV_nm,'_',model_GLM_nm,'_perf_f_absInc_GL_split',N_subs,'subs.png'];
    save_fig(fig_perf_GL_hdl, dirPath, img_nm, overwrite_allowed);
    
    %% performance = f(trial number)
    fig_trialN_hdl = fig();
    perf_T_hdl = errorbar(m_trialN_bin.(task_nm).f_trialN.all,...
        100.*m_perf_bin.(task_nm).f_trialN.all,...
        100.*sem_perf_bin.(task_nm).f_trialN.all,...
        'b *',...
        'LineWidth',lSize);
    perf_T_hdl.Color = perf_col;
    perf_T_hdl.MarkerSize = mSize;
    hold on;
    if is_EV_mdl_used
        perf_pred_EV_T_hdl = plot(m_trialN_bin.(task_nm).f_trialN.all,...
            100.*m_perf_pred_EV_bin.(task_nm).f_trialN.all,...
            '--',...
            'LineWidth',lSize);
        perf_pred_EV_T_hdl.Color = perf_pred_EV_col;
    end
    xlabel('Trial');
    ylabel('Performance (%)');
    if is_EV_mdl_used
        legend([perf_T_hdl, perf_pred_EV_T_hdl,],...
            'Performance','Predicted performance');
    else
        legend(perf_T_hdl,...
            'Performance');
    end
    legend('boxoff');
    legend_size(pSize);
    xlim([5 55]);
    img_nm = [task_nm,'_perf_and_',model_EV_nm,'_',model_GLM_nm,'_perf_f_trialN_',N_subs,'subs.png'];
    save_fig(fig_trialN_hdl, dirPath, img_nm, overwrite_allowed);
    
    %% performance = f(trial number split gain/loss)
    fig_trialN_GL_hdl_bis = fig();
    perf_T_G_hdl = errorbar(m_trialN_bin.(task_nm).f_trialN.gainTrials-dSpace,...
        100.*m_perf_bin.(task_nm).f_trialN.gainTrials,...
        100.*sem_perf_bin.(task_nm).f_trialN.gainTrials,...
        'g *',...
        'LineWidth',lSize);
    perf_T_G_hdl.Color = gain_perf_col;
    perf_T_G_hdl.MarkerSize = mSize;
    hold on;
    perf_T_L_hdl = errorbar(m_trialN_bin.(task_nm).f_trialN.lossTrials+dSpace,...
        100.*m_perf_bin.(task_nm).f_trialN.lossTrials,...
        100.*sem_perf_bin.(task_nm).f_trialN.lossTrials,...
        'k *',...
        'LineWidth',lSize);
    perf_T_L_hdl.Color = loss_perf_col;
    perf_T_L_hdl.MarkerSize = mSize;
    if is_EV_mdl_used
        perf_pred_EV_T_G_hdl = plot(m_trialN_bin.(task_nm).f_trialN.gainTrials,...
            100.*m_perf_pred_EV_bin.(task_nm).f_trialN.gainTrials,...
            '--',...
            'LineWidth',lSize);
        perf_pred_EV_T_G_hdl.Color = gain_perf_pred_EV_col;
        perf_pred_EV_T_L_hdl = plot(m_trialN_bin.(task_nm).f_trialN.lossTrials,...
            100.*m_perf_pred_EV_bin.(task_nm).f_trialN.lossTrials,...
            '--',...
            'LineWidth',lSize);
        perf_pred_EV_T_L_hdl.Color = loss_perf_pred_EV_col;
    end
    xlabel('Trial');
    ylabel('Performance (%)');
    legend_size(pSize);
    xlim([5 55]);
    img_nm = [task_nm,'_perf_and_',model_EV_nm,'_',model_GLM_nm,'_perf_f_trialN_GL_split',N_subs,'subs.png'];
    save_fig(fig_trialN_GL_hdl_bis, dirPath, img_nm, overwrite_allowed);
    
end % task loop

end % function