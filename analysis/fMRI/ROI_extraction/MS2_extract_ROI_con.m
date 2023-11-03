function [ con_vec_all, con_names, ttest_pval ] = MS2_extract_ROI_con(GLM )
%[ con_vec_all, con_names, ttest_pval ] = MS2_extract_ROI_con( GLM )
%
% INPUTS
% GLM: GLM number
%
% OUTPUTS
% con_vec_all (number of contrasts)*1*(number subjects)*(number of ROI)
% matrix containing the contrast value for each contrast, each subject and
% each ROI asked
%
% con_names: name of each contrast
%
% ttest_pval: matrix with p.value for each contrast
%

%% add paths
mainRoot = 'enter path here';

%% directories
pc_cluster = 'pc';
if ~exist('pc_cluster','var') || isempty(pc_cluster)
    pc_cluster = input(' ''pc'' or ''cluster''?','s');
end
disp(['Extracting GLM in the ',pc_cluster,'.']);
root = 'enter path here';
ROI_path = [root,filesep,'behavior_summary',filesep,'ROI',filesep];
ROI_masks_path = 'enter path here';
%% select the ROI you want to use
[ ROI_xyz, ROI_sphere_or_mask, ROI_nm, n_ROI ] = ROI_selection(ROI_masks_path);

%% which GLM
if ~exist('GLM','var') || isempty(GLM)
    GLM = spm_input('GLM number?',1,'e');
    GLMstr = num2str(GLM);
    GLMprm = which_GLM_MS2(GLM);
end

%% load subject names
subject_id = {'enter list of subjects here'};
NS = length(subject_id);

%% smoothed or unsmoothed contrasts?

%% beta or t value?
beta_or_t_values = {'beta_value','t_value'};
b_or_t = spm_input('beta value or t value?',1,'b','beta value| t.value',[1 2]);
beta_or_t_value  = beta_or_t_values{b_or_t};

%% how many figures do you want to plot
figNumber = spm_input('How many figures?',1,'e');

% prepare contrast names for question
[con_names, ~] = MS2_load_con(GLMprm, subject_id{1});
nb_max_con = length(con_names);
list_con_question = con_names{1};
for iCon = 2:nb_max_con
    list_con_question = [list_con_question,' | ',con_names{iCon}];
end
list_con_question = [list_con_question,' | If all wished contrasts for this figure selected click here'];
% select contrasts to display at the end
which_con = zeros(figNumber,nb_max_con); % column with 0/1 to see which contrasts to display at the end

%% name + contrasts to display for each figure
conName = cell(figNumber,1);
for iFig = 1:figNumber
    %% what name for each figure?
    conName{iFig} = spm_input(['Name for figure ',num2str(iFig),' please?'],1,'s');
    
    %% select which contrasts you want to display for each figure
    stop_con_loop = 0;
    while stop_con_loop == 0
        selectedContrast = spm_input(['What contrast for fig.',num2str(iFig),':',conName{iFig},' ?'],1,'m',...
            list_con_question, ...
            1:(nb_max_con+1), 0);
        if selectedContrast < nb_max_con + 1
            which_con(iFig,selectedContrast) = 1; % 1 for selected contrasts
        elseif selectedContrast == nb_max_con + 1 % stop contrast selection loop
            stop_con_loop = 1;
        end
    end % contrast loop
end % figure loop

%% create big matrix to store ROI values
con_vec_all = NaN(nb_max_con, 1, NS, n_ROI);

%% ROI loop through areas of interest
for iROI = 1:n_ROI
    
    fprintf('Region %d\n',iROI);
    
    sxyz_ROI = ROI_xyz.(['ROI_',num2str(iROI)]);
    
    %% loop through subjects for each study
    for iS = 1:NS
        
        sub_nm = subject_id{iS};
        
        subj_data_folder = fullfile(root, sub_nm, 'fMRI_analysis', 'functional', ['GLM',GLMstr]);
        
        % extract contrast for the current subject in the current ROI
        % for each contrast
        
        for iCon = 1:nb_max_con
            [ con_value ] = ROI_extraction(  iCon, subj_data_folder, sxyz_ROI, beta_or_t_value );
            con_vec_all(iCon, 1, iS, iROI) = con_value; % save mean beta for the selected ROI inside the big resulting matrix
        end % contrast
        
        disp(['Subject ',num2str(iS),'/',num2str(NS),' extracted']);
    end % subject
end % ROI

%% average across subjects
con_avg = mean(con_vec_all, 3,'omitnan');
con_sem = NaN(nb_max_con, n_ROI, 1);
for iCon = 1:nb_max_con
    for iROI = 1:n_ROI
        con_sem(iCon, iROI, 1) = sem(con_vec_all(iCon, 1, :, iROI), 3);
    end
end

%% test significativity of each contrast
ttest_pval = NaN(nb_max_con, n_ROI);
for iROI = 1:n_ROI
    for iCon = 1:nb_max_con
        curr_con = NaN(1,NS);
        curr_con(1,:) = con_vec_all(iCon, 1, :, iROI);
        [~, ttest_pval(iCon, iROI)] = ttest( curr_con );
    end % contrast loop
end % ROI loop

%% save all the data
filename = [ROI_path,beta_or_t_value,'_GLM',GLMstr,'_',num2str(n_ROI),'ROI_',num2str(NS),'subs'];
if ~exist([filename,'.mat'],'file')
        save([filename,'.mat'])
    else
        while exist([filename,'.mat'],'file')
            filename = [filename,'_bis'];
        end
        save([filename,'.mat'])
end

% %% loop through figures as defined when script launched at the
% % beginning
% for iFig = 1:figNumber
%     conFig = which_con(iFig,:);
%     selectedCon = find(conFig == 1); % extract index of contrasts selected
%     % display figure
% %     [roi_fig] = roi_graph(selectedCon, n_ROI,...
% %         con_vec_all, con_avg, con_sem, con_names, ttest_pval, ROI_nm, ROI_sphere_or_mask);
%     % save image
%     cd(ROI_path);
%     img_name = ['GLM',num2str(GLM),'_',beta_or_t_value,'_',conName{iFig},'_',num2str(n_ROI),'ROIs_',num2str(NS),'_subs.png'];
%     if exist(img_name,'file') ~= 2
%         set(gcf,'PaperPosition',[0 0 1 1]);
%         set(gcf,'PaperPositionMode','auto');
%         saveas(roi_fig,img_name);
%     else
%         warning(['There is already a file with the name ', img_name,' so the figure has not been saved yet.',...
%             ' Please do it manually or change the script to be able to do it automatically in the future.']);
%     end
% end % figure loop

end % function