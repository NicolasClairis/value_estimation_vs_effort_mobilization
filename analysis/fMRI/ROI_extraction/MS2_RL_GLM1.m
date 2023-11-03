% script to display main results RL: Val/Conf/DT in GLM1


mainResultsFolder = 'enter path here';
RL_GLM1folder = [mainResultsFolder,filesep];
%% load data
RL_dataStruct.GLM1 = load([RL_GLM1folder,'MS2_GLM1_vmPFC_mPFC_dmPFC_ROI_22subs.mat']);

%% main parameters
x = 1:3;
delta = 0.3;
lSize = 3;
pSize = 25;
red = [1 0 0];
green = [0 1 0];
blue = [0 0 1];
black = [0 0 0];
ROI_idx.vmPFC   = 1;
ROI_idx.mPFC    = 2;
ROI_idx.dmPFC   = 3;
ROIs = {'vmPFC','mPFC','dmPFC'};
nROIs = length(ROIs);
regs_RL = {'Val','Conf','DT'};
nRegs_RL = length(regs_RL);

%% extract RL data
GLM_nm = 'GLM1';
% value
Val_idx.(GLM_nm) = strcmp(RL_dataStruct.(GLM_nm).con_names,'L_mod_stim_SV_GL_Pairs_pos');
% confidence
Conf_idx.(GLM_nm) = strcmp(RL_dataStruct.(GLM_nm).con_names,'L_mod_stim_pBest_GL_Pairs_pos');
% reaction time
DT_idx.(GLM_nm) = strcmp(RL_dataStruct.(GLM_nm).con_names,'L_mod_stim_RT_GL_Pairs_pos');
NS_RL = length(RL_dataStruct.(GLM_nm).subject_id);

%% RL Val/Conf/DT for GLM with no orthogonalization
fig1 = fig;
jIdx = 0;
% loop through regressors and GLM
for iReg = 1:nRegs_RL
    reg_nm = regs_RL{iReg};
    GLM_nm = 'GLM1';
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
    subplot(1,3,jIdx);
    % convert 4D data to 2D vector
    [vmPFC_data_tmp, mPFC_data_tmp, dmPFC_data_tmp] = deal(NaN(NS_RL,1));
    vmPFC_data_tmp(:,1) = RL_dataStruct.(GLM_nm).con_vec_all(con_idx_tmp,1,:,ROI_idx.vmPFC);
    mPFC_data_tmp(:,1) = RL_dataStruct.(GLM_nm).con_vec_all(con_idx_tmp,1,:,ROI_idx.mPFC);
    dmPFC_data_tmp(:,1) = RL_dataStruct.(GLM_nm).con_vec_all(con_idx_tmp,1,:,ROI_idx.dmPFC);
    % display value/Conf/DT
    violinplot([vmPFC_data_tmp, mPFC_data_tmp, dmPFC_data_tmp],...
        {'vmPFC','mPFC','dmPFC'},...
        'ViolinColor',[col;col;col]);
    hold on;
    line([0.7 3.3],[0 0],'LineWidth',1,'Color','k');
    xticks(1:3);
    ylim([-1 1.5]);
    ylabel('Regression estimate');
    xticklabels({'vmPFC','mPFC','dmPFC'});
    xlim([x(1)-delta x(end)+delta]);
    legend_size(pSize);
end % regressor

%% check p.values
% Val
disp('avg beta for Val (vmPFC/mPFC/dmPFC)');
[mVal, semVal]=deal(NaN(1,nROIs));
mVal(1,:)=con_avg(Val_idx.GLM1,1,1,:)
semVal(1,:)=con_sem(Val_idx.GLM1,:)
disp('p.values for Val (vmPFC/mPFC/dmPFC)');
ttest_pval(Val_idx.GLM1,:)

% conf
disp('avg beta for Conf (vmPFC/mPFC/dmPFC)');
[mConf, semConf]=deal(NaN(1,nROIs));
mConf(1,:)=con_avg(Conf_idx.GLM1,1,1,:)
semConf(1,:)=con_sem(Conf_idx.GLM1,:)
disp('p.values for Conf (vmPFC/mPFC/dmPFC)');
ttest_pval(Conf_idx.GLM1,:)

% DT
disp('avg beta for DT (vmPFC/mPFC/dmPFC)');
[mDT, semDT]=deal(NaN(1,nROIs));
mDT(1,:)=con_avg(DT_idx.GLM1,1,1,:)
semDT(1,:)=con_sem(DT_idx.GLM1,:)
disp('p.values for DT (vmPFC/mPFC/dmPFC)');
ttest_pval(DT_idx.GLM1,:)