function [ betas, pval ] = MS2_GS_RT(  )
%[ betas, pval ] = MS2_GS_RT(  )
% computes correlation between RT and variables of interest of grip and
% stroop task (|incentive|, incentive and trial number mainly)
%
% OUTPUT
% betas: structure with the betas across participants
%
% pval: structure with corresponding p.values of how significant the betas
% are
%

%% working directories
root = 'enter path here';

%% subject selection
subject_id = {'enter list of subjects here'};
NS = length(subject_id);

%% load variables of interest
n_GS_runs = 2;
task_names = {'grip','stroop'};
nTasks = length(task_names);
mInc_levels = 1:6;
n_mInc = length(mInc_levels);
nInc_levels = [-6:-1, 1:6];
n_nInc = length(nInc_levels);
nTrials_per_run = 60;
trialN_bins = 5;
% norm = 1; % normalize the variables?
pSize = 18;

figRT = fig();
for iGS = 1:nTasks
    task_nm = task_names{iGS};
    
    [b_r1.(task_nm),...
        b_r2.(task_nm),...
        b_mInc.(task_nm),...
        b_nInc.(task_nm),...
        b_GLcond.(task_nm),...
        b_trialN.(task_nm)] = deal(NaN(1,NS));
    [RT_f_mInc.(task_nm), RT_fit_f_mInc.(task_nm)] = deal(NaN(n_mInc, NS));
    [RT_f_nInc.(task_nm), RT_fit_f_nInc.(task_nm)] = deal(NaN(n_nInc, NS));
    RT_f_GLcond.(task_nm) = NaN(2, NS);
    [RT_f_trialN.(task_nm), trialN_f_trialN.(task_nm)] = deal( NaN(trialN_bins, NS) );
    mRT.(task_nm) = NaN(1,NS);
    
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
        
        [trialN, trialN_idx_aRuns, time_rest, totalGain_prev,...
            inc, absInc,...
            incRank, absIncRank,...
            perf,...
            GL_cond,...
            RT_fp,...
            run1_cstt, run2_cstt] = deal( NaN(nTrials_per_run*n_GS_runs,1) );
        
        % loop through runs
        for iRun = 1:n_GS_runs
            jRun = runs_idx(iRun);
            run_nm = num2str(jRun);
            
            % load data
            loadStruct = getfield( load([onsets_folder,'onsets_sub',subid,'_',task_nm,'_run',run_nm,'.mat'],task_nm),task_nm);
            % extract relevant data
            trialN_tmp  = loadStruct.mod.all.trialN + nTrials_per_run*(iRun - 1);
            trialN(trialN_tmp)              = loadStruct.mod.all.trialN;
            trialN_idx_aRuns(trialN_tmp)    = trialN_tmp;
            incRank(trialN_tmp)                 = loadStruct.mod.all.incentiveRank;
            absIncRank(trialN_tmp)              = loadStruct.mod.all.absIncentiveRank;
            inc(trialN_tmp)                 = loadStruct.mod.all.incentive;
            absInc(trialN_tmp)              = loadStruct.mod.all.absIncentive;
            GL_cond(trialN_tmp)             = loadStruct.mod.all.trialValence;
            GL_cond(GL_cond == -1)          = 0; % transform into binary variable 0/1
            totalGain_prev(trialN_tmp)      = loadStruct.mod.all.totalGain_prev;
            perf(trialN_tmp)                = loadStruct.mod.all.perf;
            time_rest(trialN_tmp)           = loadStruct.duration.all.ITI + loadStruct.duration.all.incentive; % between [1.5; 4.95]seconds
            %             switch task_nm
            %                 case 'grip'
            RT_fp(trialN_tmp)  = loadStruct.mod.all.RT_fp;
            %                 case 'stroop'
            %                     RT_fp(trialN_tmp)  = loadStruct.mod.all.RT_fp;
            %             end
            if jRun < 4
                run1_cstt(trialN_tmp) = 1;
                run2_cstt(trialN_tmp) = 0;
            else
                run1_cstt(trialN_tmp) = 0;
                run2_cstt(trialN_tmp) = 1;
            end
        end % run loop
        
        % add incentive per condition
        [absIncGain, absIncLoss] = deal( absInc );
        [absIncRankGain, absIncRankLoss] = deal( absIncRank );
        absIncGain(GL_cond == 0) = 0;
        absIncRankGain(GL_cond == 0) = 0;
        absIncLoss(GL_cond == 1) = 0;
        absIncRankLoss(GL_cond == 1) = 0;
        
%         %% normalize vars
%         if norm == 1
%             trialN  = trialN./nTrials_per_run;
%             inc     = inc./6;% linear value -6 => +61
%             absInc  = absInc./6;% motiv value +1 => +6
%             absIncGain = absIncGain./6;
%             absIncLoss = absIncLoss./6;
%             time_rest = time_rest./4.95; % if you want to normalize it,
%             %         maximum possible = 4.95 seconds of rest (0.5s cross + 3.95s
%             %         incentive)
%             maxTotalGain = 267.10; % maximal amount possible to win if maximal perf in all trials for both runs
%             totalGain_prev = totalGain_prev./maxTotalGain; % normalize by maximum possible amount to cumulate across runs
%         end
        
        %% perform the glm
        x_reg = [run1_cstt, run2_cstt, absInc, trialN];
        % orthogonalize the regressors (especially for incentive linear to
        % condition)
        %
        % remove NaN for orthogonalization
        ok_trials = ~isnan(absInc);
        x_reg = x_reg(ok_trials,:);
        % orthogonalize the regressors
        x_reg = mtrx_orthog(x_reg);
        % extract the betas
        betas_tmp = glmfit(x_reg, RT_fp(ok_trials), 'normal','constant','off');
        RT_fp_fit = glmval(betas_tmp,x_reg,'identity','constant','off');
        
        b_r1.(task_nm)(iS)       = betas_tmp(1);
        b_r2.(task_nm)(iS)       = betas_tmp(2);
        iB = 3;
        b_mInc.(task_nm)(iS)    = betas_tmp(iB);
        iB = iB + 1;
        b_trialN.(task_nm)(iS)  = betas_tmp(iB);
        
        %% bins
        RT_f_mInc.(task_nm)(:,iS)   = do_bin2(RT_fp, absInc, n_mInc, 0);
        RT_fit_f_mInc.(task_nm)(:,iS) = do_bin2(RT_fp_fit, absInc(ok_trials), n_mInc, 0);
        RT_f_nInc.(task_nm)(:,iS)   = do_bin2(RT_fp, inc, n_nInc, 0);
        RT_fit_f_nInc.(task_nm)(:,iS)   = do_bin2(RT_fp_fit, inc(ok_trials), n_nInc, 0);
        RT_f_trialN.(task_nm)(:,iS) = do_bin2(RT_fp, trialN, trialN_bins, 0);
        trialN_f_trialN.(task_nm)(:,iS) = do_bin2(trialN, trialN, trialN_bins, 0);
        RT_f_GLcond.(task_nm)(1,iS)   = mean(RT_fp(GL_cond == 0),1,'omitnan');
        RT_f_GLcond.(task_nm)(2,iS)   = mean(RT_fp(GL_cond == 1),1,'omitnan');
        
        %% mean across trials
        mRT.(task_nm)(iS) = mean(RT_fp,1,'omitnan');
    end % subject loop
    
    %% average across subjects
    % betas
    [mB_r1.(task_nm),semB_r1.(task_nm)]      = mean_sem_sd(b_r1.(task_nm),2);
    [mB_r2.(task_nm),semB_r2.(task_nm)]      = mean_sem_sd(b_r2.(task_nm),2);
    [mB_mInc.(task_nm),semB_mInc.(task_nm)]   = mean_sem_sd(b_mInc.(task_nm),2);
    [mB_trialN.(task_nm),semB_trialN.(task_nm)] = mean_sem_sd(b_trialN.(task_nm),2);
    [~, pval.b_r1.(task_nm)]         = ttest( b_r1.(task_nm) );
    [~, pval.b_r2.(task_nm)]         = ttest( b_r2.(task_nm) );
    [~, pval.b_mInc.(task_nm)]      = ttest( b_mInc.(task_nm) );
    [~, pval.b_trialN.(task_nm)]    = ttest( b_trialN.(task_nm) );
    % bins
    [mRT_f_mInc.(task_nm),semRT_f_mInc.(task_nm)]        = mean_sem_sd(RT_f_mInc.(task_nm),2);
    mRT_fit_f_mInc.(task_nm)        = mean_sem_sd(RT_fit_f_mInc.(task_nm),2);
    [mRT_f_GLcond.(task_nm),semRT_f_GLcond.(task_nm)]      = mean_sem_sd(RT_f_GLcond.(task_nm),2);
    [mRT_f_nInc.(task_nm),semRT_f_nInc.(task_nm)]        = mean_sem_sd(RT_f_nInc.(task_nm),2);
    mRT_fit_f_nInc.(task_nm)        = mean_sem_sd(RT_fit_f_nInc.(task_nm),2);
    [mRT_f_trialN.(task_nm),semRT_f_trialN.(task_nm)]      = mean_sem_sd(RT_f_trialN.(task_nm),2);
    mTrialN_f_trialN.(task_nm)  = mean_sem_sd(trialN_f_trialN.(task_nm),2);

    
    %% store
    betas.b_r1.(task_nm)         = b_r1.(task_nm);
    betas.b_r2.(task_nm)         = b_r2.(task_nm);
    betas.b_mInc.(task_nm)      = b_mInc.(task_nm);
    betas.b_trialN.(task_nm)    = b_trialN.(task_nm);
    betas.mB_r1.(task_nm)        = mB_r1.(task_nm);
    betas.mB_mInc.(task_nm)     = mB_mInc.(task_nm);
    betas.mB_trialN.(task_nm)   = mB_trialN.(task_nm);
    
    % store also SEM
    betas.semB_mInc.(task_nm)     = semB_mInc.(task_nm);
    betas.semB_trialN.(task_nm)   = semB_trialN.(task_nm);
        
    %% graph
    figure(figRT);
    
    % RT = f(motiv inc)
    subplot(2, 4, 1 + 4*(iGS-1));
    errorbar(mInc_levels, mRT_f_mInc.(task_nm), semRT_f_mInc.(task_nm),...
        'b* ', 'LineWidth',3);
    hold on;
    xlim([0.8 6.2]);
    ylabel([task_nm,' RT (s)']);
    xlabel('|Incentive|');
    set(gca,'fontsize',pSize,'FontWeight','bold');
    set(findall(gcf,'type','text'),'FontSize',pSize,'fontWeight','bold');
    
    % RT = f(inc linear)
    subplot(2, 4, 3 + 4*(iGS-1));
    errorbar(nInc_levels, mRT_f_nInc.(task_nm), semRT_f_nInc.(task_nm),...
        'b* ', 'LineWidth',3);
    hold on;
%     plot(nInc_levels, mb_r1.(task_nm) + mB_nInc.(task_nm).*nInc_levels,...
%         'r','LineWidth',3);
    xlim([-6.2 6.2]);
    ylabel([task_nm,' RT (s)']);
    xlabel('Incentive');
    set(gca,'fontsize',pSize,'FontWeight','bold');
    set(findall(gcf,'type','text'),'FontSize',pSize,'fontWeight','bold');
    
    % RT = f(trialN)
    subplot(2, 4, 4 + 4*(iGS-1));
    errorbar(mTrialN_f_trialN.(task_nm), mRT_f_trialN.(task_nm), semRT_f_trialN.(task_nm),...
        'b* ', 'LineWidth',3);
    hold on;
    xlim([1 60]);
    ylabel([task_nm,' RT (s)']);
    xlabel('Trial number');
    legend_size(pSize);
    
    %% figure apart with RT=f(trialN) alone
    figRT_bis = fig();
    % RT = f(trialN)
    errorbar(mTrialN_f_trialN.(task_nm), mRT_f_trialN.(task_nm), semRT_f_trialN.(task_nm),...
        'b* ', 'LineWidth',3);
    hold on;
    xlim([1 60]);
    ylabel([task_nm,' RT (s)']);
    xlabel('Trial number');
    legend_size(pSize);
    % save
    img_name = [root,filesep, 'behavior_summary', filesep task_nm filesep 'RT_f_trialN_' num2str(NS),'subs.png'];
    set(gcf,'PaperPosition',[0 0 1 1]);
    set(gcf,'PaperPositionMode','auto');
    saveas(figRT_bis,img_name);
    
end % task loop

% save
img_name = [root,filesep, 'behavior_summary', filesep 'RT_f_absInc_Inc_trialN_' num2str(NS),'subs.png'];
set(gcf,'PaperPosition',[0 0 1 1]);
set(gcf,'PaperPositionMode','auto');
saveas(figRT,img_name);

%% plot correlation of the betas between tasks

% beta0
% check correlation
[betas_r1_tmp,~,stats_r1] = glmfit(betas.b_r1.grip, betas.b_r1.stroop, 'normal');
% store results in function output
betas.b_r1_G_vs_S = betas_r1_tmp;
pval.b_r1_G_vs_S = stats_r1.p;
% check correlation
[betas_r2_tmp,~,stats_r2] = glmfit(betas.b_r2.grip, betas.b_r2.stroop, 'normal');
% fit_b0 = glmval(betas0_tmp, betas.b_r1.grip, 'identity');
% store results in function output
betas.b_r2_G_vs_S = betas_r2_tmp;
pval.b_r2_G_vs_S = stats_r2.p;


fig_betas = fig();

% beta |incentive|
% check correlation
[betas_mInc_tmp,~,stats_mInc] = glmfit(betas.b_mInc.grip, betas.b_mInc.stroop, 'normal');
fit_mInc = glmval(betas_mInc_tmp, betas.b_mInc.grip, 'identity');
% store results in function output
betas.b_mInc_G_vs_S = betas_mInc_tmp;
pval.b_mInc_G_vs_S = stats_mInc.p;

subplot(1,4,1);
scatter(betas.b_mInc.grip, betas.b_mInc.stroop,...
    'LineWidth',3);
hold on;
plot(betas.b_mInc.grip, fit_mInc, 'r--','LineWidth',3);
xlabel('Grip |incentive| betas on RT');
ylabel('Stroop |incentive| betas on RT');
legend_size(pSize);

% beta trial number
% check correlation
[betas_trialN_tmp,~,stats_trialN] = glmfit(betas.b_trialN.grip, betas.b_trialN.stroop, 'normal');
fit_trialN = glmval(betas_trialN_tmp, betas.b_trialN.grip, 'identity');
% store results in function output
betas.b_trialN_G_vs_S = betas_trialN_tmp;
pval.b_trialN_G_vs_S = stats_trialN.p;

subplot(1,4,4);
scatter(betas.b_trialN.grip, betas.b_trialN.stroop,...
    'LineWidth',3);
hold on;
plot(betas.b_trialN.grip, fit_trialN, 'r--','LineWidth',3);
xlabel('Grip trial number betas on RT');
ylabel('Stroop trial number betas on RT');
legend_size(pSize);

% save
img_name = [root,filesep, 'behavior_summary', filesep 'RT_betas_correl_btw_G_and_S_a' num2str(NS),'Subs.png'];
set(gcf,'PaperPosition',[0 0 1 1]);
set(gcf,'PaperPositionMode','auto');
saveas(fig_betas,img_name);

%% correl between RT
[betas_tmp,~,stats_tmp] = glmfit(mRT.grip, mRT.stroop,'normal');
betas.mRT = betas_tmp;
pval.mRT = stats_tmp.p;
fit_RT = glmval(betas_tmp, mRT.grip, 'identity');
fig_RT = fig();
scatter(mRT.grip, mRT.stroop,'b','LineWidth',3);
hold on;
plot(mRT.grip, fit_RT, 'r--','LineWidth',3);
xlabel('mean RT in grip');
ylabel('mean RT in stroop');
legend_size(pSize);

% save
img_name = [root,filesep, 'behavior_summary', filesep 'RT_correl_btw_G_and_S_a' num2str(NS),'Subs.png'];
set(gcf,'PaperPosition',[0 0 1 1]);
set(gcf,'PaperPositionMode','auto');
saveas(fig_RT,img_name);

%% RT=f(I)
xIncentives = [-6:-1, 1:6];
grey = [143 143 143]./255;
fig;

% grip
subplot(1,2,1);
jbfill(xIncentives,...
    (mRT_f_nInc.grip + semRT_f_nInc.grip)',...
    (mRT_f_nInc.grip - semRT_f_nInc.grip)',...
    (mRT_f_nInc.grip)',...
    grey,grey,1,0.5);
hold on;
plot(xIncentives, mRT_fit_f_nInc.grip,...
    'Color','k','LineWidth',3,'LineStyle','--');
xlim([-6.2 6.2]);

ylabel('RT (s)');
xlabel('Incentive (€)');
xticks([-6:-1, 1:6]);
xticklabels({'-20.00','-5.00','-1.00','-0.50','-0.20','-0.01','0.01','0.20','0.50','1.00','5.00','20.00'})
ylabel('RT (s)');
legend_size(30);

% stroop
subplot(1,2,2);
jbfill(xIncentives,...
    (mRT_f_nInc.stroop + semRT_f_nInc.stroop)',...
    (mRT_f_nInc.stroop - semRT_f_nInc.stroop)',...
    (mRT_f_nInc.stroop)',...
    grey,grey,1,0.5);
hold on;
plot(xIncentives, mRT_fit_f_nInc.stroop,...
    'Color','k','LineWidth',3,'LineStyle','--');
ylabel('RT (s)');
xlabel('Incentive (€)');
xticks(xIncentives);
xticklabels({'-20.00','-5.00','-1.00','-0.50','-0.20','-0.01','0.01','0.20','0.50','1.00','5.00','20.00'})
ylabel('RT (s)');
legend_size(30);

%% show RT = f(|incentive|) separately for gains and losses
close all;
xIncentives = 1:6;
gain_col = [165 15 21]./255;
loss_col = [251 106 74]./255;
gain_idx = 7:12;
loss_idx = 6:(-1):1; % start with smaller amount -0.01 and go to -20€
lSize = 5;
mSize = 30;
dSpace = 0.05;
lWidth = 3;

for iTask = 1:nTasks
    task_nm = task_names{iTask};
    fig;
    perf_G_hdl = errorbar(xIncentives-dSpace,...
        mRT_f_nInc.(task_nm)(gain_idx), semRT_f_nInc.(task_nm)(gain_idx),...
        'g *','LineWidth',lSize,'LineStyle','none');
    perf_G_hdl.Color = gain_col;
    perf_G_hdl.Marker = 'o';
    perf_G_hdl.MarkerSize = mSize;
    hold on;
    perf_G_fit_hdl = plot(xIncentives-dSpace,...
        mRT_fit_f_nInc.(task_nm)(gain_idx));
    perf_G_fit_hdl.Color = gain_col;
    perf_G_fit_hdl.LineStyle = '--';
    perf_G_fit_hdl.LineWidth = lWidth;
    perf_L_hdl = errorbar(xIncentives+dSpace,...
        mRT_f_nInc.(task_nm)(loss_idx), semRT_f_nInc.(task_nm)(loss_idx),...
        'g *','LineWidth',lSize,'LineStyle','none');
    perf_L_hdl.Color = loss_col;
    perf_L_hdl.Marker = 's';
    perf_L_hdl.MarkerSize = mSize;
    perf_L_fit_hdl = plot(xIncentives+dSpace,...
        mRT_fit_f_nInc.(task_nm)(loss_idx));
    perf_L_fit_hdl.Color = loss_col;
    perf_L_fit_hdl.LineStyle = '--';
    perf_L_fit_hdl.LineWidth = lWidth;
    xlim([0.8 6.2]);

    ylabel('RT (s)');
    xlabel('Incentive (€)');
    xticks(xIncentives);
    xticklabels({'0.01','0.20','0.50','1.00','5.00','20.00'})
    ylabel('RT (s)');
    legend_size(60);
end % task loop

end % function