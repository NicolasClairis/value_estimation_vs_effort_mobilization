%% script to load ROI data of GLM 156, 157 and 158 for Motiscan 2


%% define working directories
mainResultsFolder = 'enter path here';
%% load data
GLMS = {'GLM2','GLM3','GLM4'};
nGLMs = length(GLMS);
for iGLM = 1:nGLMs
    GLM_nm = GLMS{iGLM};
    GS_dataStruct.(GLM_nm) = load([mainResultsFolder, filesep,...
        GLM_nm,'_vmPFC_dmPFC_GS_22subs.mat']);
end % GLM loop
%% main parameters
red = [1 0 0];
green = [0 1 0];
blue = [0 0 1];
black = [0 0 0];
lSize = 3;
pSize = 25;

ROIs = {'vmPFC','dmPFC'};
nROIs = length(ROIs);
regs_GS = {'E*','RT'};
nRegs_GS = length(regs_GS);

%% extract GS data
for iGLM = 1:nGLMs
    GLM_nm = GLMS{iGLM};
    switch GLM_nm
        case {'GLM2','GLM3'}
            period_nm = 'mod_inc';
        case 'GLM4'
            period_nm = 'mod_dispE';
    end
    % split per task
    % resource E
    E_S_idx.(GLM_nm) = strcmp(GS_dataStruct.(GLM_nm).con_names,['S_',period_nm,'_E_pred_pos']);
    E_G_idx.(GLM_nm) = strcmp(GS_dataStruct.(GLM_nm).con_names,['G_',period_nm,'_E_pred_pos']);
    % reaction time
    RT_S_idx.(GLM_nm) = strcmp(GS_dataStruct.(GLM_nm).con_names,['S_',period_nm,'_RT_fp_pos']);
    RT_G_idx.(GLM_nm) = strcmp(GS_dataStruct.(GLM_nm).con_names,['G_',period_nm,'_RT_fp_pos']);
end
% ROI indexes
ROI_idx.vmPFC   = 1;
ROI_idx.dmPFC   = 2;
% number of subjects
NS_GS = length(GS_dataStruct.(GLM_nm).subject_id);

%% GS X/RT for the GLM with no orthogonalization spliting each task
fig1 = fig;
for iGLM = 1:nGLMs
    GLM_nm = GLMS{iGLM};
    
    %% X plot
    xlim_vals = [0.5 2.5];
    % loop through regressors and GLM
    for iReg = 1:nRegs_GS
        GS_reg_nm = regs_GS{iReg};
        subplot(nGLMs, 3, (nRegs_GS+1)*(iGLM-1) + iReg); hold on;
        % determine color and index
        switch GS_reg_nm
            case 'E*'
                con_idx_S_tmp = E_S_idx.(GLM_nm);
                con_idx_G_tmp = E_G_idx.(GLM_nm);
                col = green;
            case 'RT'
                con_idx_S_tmp = RT_S_idx.(GLM_nm);
                con_idx_G_tmp = RT_G_idx.(GLM_nm);
                col = black;
        end
        
        % convert 4D data to 2D vector
        [dmPFC_data_G_tmp,...
            dmPFC_data_S_tmp] = deal(NaN(NS_GS,1));
        dmPFC_data_G_tmp(:,1) = GS_dataStruct.(GLM_nm).con_vec_all(con_idx_G_tmp,1,:,ROI_idx.dmPFC);
        dmPFC_data_S_tmp(:,1) = GS_dataStruct.(GLM_nm).con_vec_all(con_idx_S_tmp,1,:,ROI_idx.dmPFC);
        violinplot([dmPFC_data_G_tmp, dmPFC_data_S_tmp],...
            {'G','S'},...
            'ViolinColor',[col;col]);
        hold on;
        line(xlim_vals,[0 0],'LineWidth',1,'Color','k');
        xticks(1:4);
        ylim([-5 15]);
        xlim(xlim_vals);
        ylabel('Regression estimate');
        xticklabels({'G','S'});
        legend_size(pSize);
    end % regressor
end % GLM loop

%% extract mean and SEM beta and p.value
for iGLM = 1:nGLMs
    GLM_nm = GLMS{iGLM};
    for iROI  = 1:2
        ROI_nm = ROIs{iROI};
        % extract mean and SEM for beta
        m_beta.(GLM_nm).E_S.(ROI_nm) = GS_dataStruct.(GLM_nm).con_avg(E_S_idx.(GLM_nm),1,1,ROI_idx.(ROI_nm));
        sem_beta.(GLM_nm).E_S.(ROI_nm) = GS_dataStruct.(GLM_nm).con_sem(E_S_idx.(GLM_nm),ROI_idx.(ROI_nm));
        m_beta.(GLM_nm).E_G.(ROI_nm) = GS_dataStruct.(GLM_nm).con_avg(E_G_idx.(GLM_nm),1,1,ROI_idx.(ROI_nm));
        sem_beta.(GLM_nm).E_G.(ROI_nm) = GS_dataStruct.(GLM_nm).con_sem(E_G_idx.(GLM_nm),ROI_idx.(ROI_nm));
        m_beta.(GLM_nm).RT_S.(ROI_nm) = GS_dataStruct.(GLM_nm).con_avg(RT_S_idx.(GLM_nm),1,1,ROI_idx.(ROI_nm));
        sem_beta.(GLM_nm).RT_S.(ROI_nm) = GS_dataStruct.(GLM_nm).con_sem(RT_S_idx.(GLM_nm),ROI_idx.(ROI_nm));
        m_beta.(GLM_nm).RT_G.(ROI_nm) = GS_dataStruct.(GLM_nm).con_avg(RT_G_idx.(GLM_nm),1,1,ROI_idx.(ROI_nm));
        sem_beta.(GLM_nm).RT_G.(ROI_nm) = GS_dataStruct.(GLM_nm).con_sem(RT_G_idx.(GLM_nm),ROI_idx.(ROI_nm));
        % extract p.value
        pval.(GLM_nm).E_S.(ROI_nm) = GS_dataStruct.(GLM_nm).ttest_pval(E_S_idx.(GLM_nm), ROI_idx.(ROI_nm));
        pval.(GLM_nm).E_G.(ROI_nm) = GS_dataStruct.(GLM_nm).ttest_pval(E_G_idx.(GLM_nm), ROI_idx.(ROI_nm));
        pval.(GLM_nm).RT_S.(ROI_nm) = GS_dataStruct.(GLM_nm).ttest_pval(RT_S_idx.(GLM_nm), ROI_idx.(ROI_nm));
        pval.(GLM_nm).RT_G.(ROI_nm) = GS_dataStruct.(GLM_nm).ttest_pval(RT_G_idx.(GLM_nm), ROI_idx.(ROI_nm));
    end % ROI loop
end % GLM loop