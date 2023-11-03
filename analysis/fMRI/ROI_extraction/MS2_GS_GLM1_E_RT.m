

GLM1_folder = 'enter path here';
GS_dataStruct.GLM1 = load([GLM1_folder,...
    'GLM1_vmPFC_dmPFC_GS_22subs.mat']);

%% main parameters
red = [1 0 0];
green = [0 1 0];
blue = [0 0 1];
black = [0 0 0];
lSize = 3;
pSize = 30;

GLMs = {'GLM1'};
nGLM = length(GLMs);
ROIs = {'vmPFC','dmPFC'};
nROIs = length(ROIs);
regs_GS = {'E','RT'};
nRegs_GS = length(regs_GS);

%% extract GS data
for iGLM = 1:nGLM
    GLM_nm = GLMs{iGLM};
    % resource E
    E_idx.(GLM_nm) = strcmp(GS_dataStruct.(GLM_nm).con_names,'GSpool_mod_inc_E_pred_pos');
    % reaction time
    RT_idx.(GLM_nm) = strcmp(GS_dataStruct.(GLM_nm).con_names,'GSpool_mod_inc_RT_fp_pos');

    % split per task
    % resource E
    E_S_idx.(GLM_nm) = strcmp(GS_dataStruct.(GLM_nm).con_names,'S_mod_inc_E_pred_pos');
    E_G_idx.(GLM_nm) = strcmp(GS_dataStruct.(GLM_nm).con_names,'G_mod_inc_E_pred_pos');
    % reaction time
    RT_S_idx.(GLM_nm) = strcmp(GS_dataStruct.(GLM_nm).con_names,'S_mod_inc_RT_fp_pos');
    RT_G_idx.(GLM_nm) = strcmp(GS_dataStruct.(GLM_nm).con_names,'G_mod_inc_RT_fp_pos');
end
% ROI indexes
ROI_idx.vmPFC   = 1;
ROI_idx.dmPFC   = 2;
% number of subjects
NS_GS = length(GS_dataStruct.(GLM_nm).subject_id);

%% GS X/RT for the GLM with no orthogonalization spliting each task
fig1 = fig;
jIdx = 0;
xlim_vals = [0.5 4.5];
% loop through regressors and GLM
for iReg = 1:nRegs_GS
    GS_reg_nm = regs_GS{iReg};
    GLM_nm = 'GLM1';
    % determine color and index
    switch GS_reg_nm
        case 'E'
            con_idx_S_tmp = E_S_idx.(GLM_nm);
            con_idx_G_tmp = E_G_idx.(GLM_nm);
            col = green;
        case 'RT'
            con_idx_S_tmp = RT_S_idx.(GLM_nm);
            con_idx_G_tmp = RT_G_idx.(GLM_nm);
            col = black;
    end

    % reload figure
    figure(fig1);
    % figure index
    jIdx = jIdx + 1;
    subplot(1,2,jIdx);
    % convert 4D data to 2D vector
    [vmPFC_data_G_tmp,...
        vmPFC_data_S_tmp,...
        dmPFC_data_G_tmp,...
        dmPFC_data_S_tmp,...
        VS_data_tmp] = deal(NaN(NS_GS,1));
    vmPFC_data_G_tmp(:,1) = GS_dataStruct.(GLM_nm).con_vec_all(con_idx_G_tmp,1,:,ROI_idx.vmPFC);
    vmPFC_data_S_tmp(:,1) = GS_dataStruct.(GLM_nm).con_vec_all(con_idx_S_tmp,1,:,ROI_idx.vmPFC);
    dmPFC_data_G_tmp(:,1) = GS_dataStruct.(GLM_nm).con_vec_all(con_idx_G_tmp,1,:,ROI_idx.dmPFC);
    dmPFC_data_S_tmp(:,1) = GS_dataStruct.(GLM_nm).con_vec_all(con_idx_S_tmp,1,:,ROI_idx.dmPFC);
    %         VS_data_tmp(:,1)    = GS_VS_dataStruct.(GLM_nm).con_vec_all(con_idx_tmp,1,:,1);
    % display value/Conf/DT
    %     violinplot([vmPFC_data_S_tmp, vmPFC_data_G_tmp,...
    %         mPFC_data_S_tmp, mPFC_data_G_tmp,...
    %         dmPFC_data_S_tmp, dmPFC_data_G_tmp],...
    %         {'S','G','S','G','S','G'},...
    %         'ViolinColor',[col;col;col;col;col;col]);
    violinplot([vmPFC_data_G_tmp, vmPFC_data_S_tmp,...
        dmPFC_data_G_tmp, dmPFC_data_S_tmp],...
        {'G','S','G','S'},...
        'ViolinColor',[col;col;col;col]);
    %         violinplot([vmPFC_data_tmp, mPFC_data_tmp, dmPFC_data_tmp, VS_data_tmp],...
    %             {'vmPFC','mPFC','dmPFC','VS'},...
    %             'ViolinColor',[col;col;col;col]);
    hold on;
    line(xlim_vals,[0 0],'LineWidth',1,'Color','k');
    xticks(1:6);
    ylim([-4 15]);
    xlim(xlim_vals);
    ylabel('Regression estimate');
    %         xticklabels({'vmPFC','mPFC','dmPFC','VS'});
    xticklabels({'G','S','G','S'});
    legend_size(pSize);
end % regressor

%% redo the same without the vmPFC
fig2 = fig;
jIdx = 0;
xlim_vals = [0.5 2.5];
% loop through regressors and GLM
for iReg = 1:nRegs_GS
    GS_reg_nm = regs_GS{iReg};
    GLM_nm = 'GLM1';
    % determine color and index
    switch GS_reg_nm
        case 'E'
            con_idx_S_tmp = E_S_idx.(GLM_nm);
            con_idx_G_tmp = E_G_idx.(GLM_nm);
            col = green;
        case 'RT'
            con_idx_S_tmp = RT_S_idx.(GLM_nm);
            con_idx_G_tmp = RT_G_idx.(GLM_nm);
            col = black;
    end

    % reload figure
    figure(fig2);
    % figure index
    jIdx = jIdx + 1;
    subplot(1,2,jIdx);
    % convert 4D data to 2D vector
    [dmPFC_data_G_tmp,...
        dmPFC_data_S_tmp,...
        VS_data_tmp] = deal(NaN(NS_GS,1));
    dmPFC_data_G_tmp(:,1) = GS_dataStruct.(GLM_nm).con_vec_all(con_idx_G_tmp,1,:,ROI_idx.dmPFC);
    dmPFC_data_S_tmp(:,1) = GS_dataStruct.(GLM_nm).con_vec_all(con_idx_S_tmp,1,:,ROI_idx.dmPFC);
    violinplot([dmPFC_data_G_tmp, dmPFC_data_S_tmp],...
        {'G','S'},...
        'ViolinColor',[col;col]);
    hold on;
    line(xlim_vals,[0 0],'LineWidth',1,'Color','k');
    xticks(1:2);
    ylim([-3 10]);
    xlim(xlim_vals);
    ylabel('Regression estimate');
    xticklabels({'G','S'});
    legend_size(pSize);
end % regressor

%% extract mean and SEM beta and p.value
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
end