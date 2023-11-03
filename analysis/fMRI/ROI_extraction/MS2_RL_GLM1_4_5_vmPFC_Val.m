%% script to load ROI data of GLM 156, 1 and 158 for Motiscan 2


%% define working directories
mainResultsFolder = 'enter path here';

%% load data
RL_dataStruct.GLM1 = load([mainResultsFolder,'ROI_data_GLM1_22subs.mat']);
RL_dataStruct.GLM4 = load([mainResultsFolder,'ROI_data_GLM4_22subs.mat']);
RL_dataStruct.GLM5 = load([mainResultsFolder,'ROI_data_GLM5_22subs.mat']);
close all;

%% main parameters
red = [1 0 0];
green = [0 1 0];
blue = [0 0 1];
black = [0 0 0];
lSize = 3;
pSize = 30;

GLMs = {'GLM1','GLM4','GLM5'};
nGLM = length(GLMs);
ROIs = {'vmPFC','mmPFC','dmPFC'};
nROIs = length(ROIs);
regs_RL = {'Val','Conf','DT'};
nRegs_RL = length(regs_RL);

%% extract RL data
for iGLM = 1:nGLM
    GLM_nm = GLMs{iGLM};
    % value
    Val_idx.(GLM_nm) = strcmp(RL_dataStruct.(GLM_nm).con_names,'L_mod_stim_SV_GL_Pairs_pos');
    % confidence
    Conf_idx.(GLM_nm) = strcmp(RL_dataStruct.(GLM_nm).con_names,'L_mod_stim_pBest_GL_Pairs_pos');
    % reaction time
    DT_idx.(GLM_nm) = strcmp(RL_dataStruct.(GLM_nm).con_names,'L_mod_stim_RT_GL_Pairs_pos');
end

% ROI indexes
ROI_idx.vmPFC   = 1;
ROI_idx.mmPFC    = 2;
ROI_idx.dmPFC   = 3;
% number of subjects
NS_RL = length(RL_dataStruct.(GLM_nm).subject_id);

%% RL Val/Conf/DT for GLM with no orthogonalization for each GLM
fig1 = fig;
jIdx = 0;
% loop through regressors and GLM
for iReg = 1:nRegs_RL
    reg_nm = regs_RL{iReg};
    for iGLM = 1:nGLM
        GLM_nm = GLMs{iGLM};
        % determine color and index
        switch reg_nm
            case 'Val'
                con_idx_tmp = Val_idx.(GLM_nm);
                col = red;
            case 'Conf'
                con_idx_tmp = Conf_idx.(GLM_nm);
                col = blue;
            case 'DT'
                con_idx_tmp = DT_idx.(GLM_nm);
                col = green;
        end

        % reload violin plot figure grouped by variable instead of ROI
        figure(fig1);
        % figure index
        jIdx = jIdx + 1;
        subplot(nGLM,3,jIdx);
        % convert 4D data to 2D vector
        [vmPFC_data_tmp, mmPFC_data_tmp, dmPFC_data_tmp] = deal(NaN(NS_RL,1));
        vmPFC_data_tmp(:,1) = RL_dataStruct.(GLM_nm).con_vec_all(con_idx_tmp,1,:,ROI_idx.vmPFC);
        mmPFC_data_tmp(:,1) = RL_dataStruct.(GLM_nm).con_vec_all(con_idx_tmp,1,:,ROI_idx.mmPFC);
        dmPFC_data_tmp(:,1) = RL_dataStruct.(GLM_nm).con_vec_all(con_idx_tmp,1,:,ROI_idx.dmPFC);
        % display value/Conf/DT
        violinplot([vmPFC_data_tmp, mmPFC_data_tmp, dmPFC_data_tmp],...
            {'vmPFC','mmPFC','dmPFC'},...
            'ViolinColor',[col;col;col]);
        hold on;
        line([0.7 3.3],[0 0],'LineWidth',1,'Color','k');
        xticks(1:3);
        ylim([-2 2]);
        ylabel('Regression estimate');
        xticklabels({'vmPFC','mmPFC','dmPFC'});
        legend_size(pSize);
    end % GLM
end % regressor
%% RL Val for each GLM to compare
fig;
subplot(1,3,1);hold on;
reg_nm = 'Val';
col = red;
% extract val for each GLM and then pool them together in one plot
for iGLM = 1:nGLM
    GLM_nm = GLMs{iGLM};
    vmPFC_data.(GLM_nm) = NaN(NS_RL,1);
    vmPFC_data.(GLM_nm)(:) = RL_dataStruct.(GLM_nm).con_vec_all(Val_idx.(GLM_nm),1,:,ROI_idx.vmPFC);
    % extract p.values
    pval.vmPFC_Val.(GLM_nm) = RL_dataStruct.(GLM_nm).ttest_pval(Val_idx.(GLM_nm), ROI_idx.vmPFC);
end
% violinplot([vmPFC_data.GLM1, vmPFC_data.GLM4, vmPFC_data.GLM5],...
%             {'V_c_h + V_u_n_c_h','V_c_h','V_c_h - V_u_n_c_h'},...
%             'ViolinColor',[col;col;col]);
% hold on;
% line([0.5 3.5],[0 0],'LineWidth',1,'Color','k');
% xlim([0.5 3.5]);

violinplot([vmPFC_data.GLM1, vmPFC_data.GLM4, vmPFC_data.GLM5],...
            {'pCHQCH + pUCQUC','QCH + QUC','QCH - QUC'},...
            'ViolinColor',[col;col;col]);
line([0.5 3.5],[0 0],'LineWidth',1,'Color','k');
xlim([0.5 3.5]);

ylim([-1 1.5]);
ylabel('Regression estimate');
legend_size(pSize);

% compare GLMs with paired t.test
[~,pval.vmPFC_Val.GLM1_vs_4] = ttest(vmPFC_data.GLM1, vmPFC_data.GLM4);
[~,pval.vmPFC_Val.GLM1_vs_5] = ttest(vmPFC_data.GLM1, vmPFC_data.GLM5);
[~,pval.vmPFC_Val.GLM4_vs_5] = ttest(vmPFC_data.GLM4, vmPFC_data.GLM5);
