function [ mR2 ] = MS2_GS_group_model_EV( )
%[ mR2 ] = MS2_GS_group_model_EV( )
%MS2_GS_group_model_EV launches MS2_GS_plot_model_EV.m to create one file
%per subject per task per run with the expected value, benefit and cost for
%each subject, task and run and each model.
%
% OUTPUTS
% mR2: structure with R² computed for each model, each run, each subject.
% Allows you to quick check how changing the model or the constraints
% improves the model.

%% working directories
root = 'enter path here';

%% subject selection
subject_id = {'enter list of subjects here'};
NS = length(subject_id);

%% tasks
task_names = {'grip','stroop'};
nTasks = length(task_names);

%% extract number of models to check
[ prm_in_mdl ] = MS2_GS_Festimation_model_space();
all_models = fieldnames(prm_in_mdl);
mdl_to_include_idx = listdlg('ListString',all_models);
models = all_models(mdl_to_include_idx);
nModels = length(models);

%% 
for iTask = 1:nTasks
    task_nm = task_names{iTask};
    for iMdl = 1:nModels
        mdl_nm = models{iMdl};
        mR2.(task_nm).(mdl_nm) = NaN(2,NS);
    end
    
    for iS = 1:NS
        sub_nm = subject_id{iS};
        if strcmp(sub_nm(3),'_')
            subid = sub_nm(2);
        elseif ~strcmp(sub_nm(3),'_') && strcmp(sub_nm(4),'_')
            subid = sub_nm(2:3);
        end
        
        run_idx = MS2_task_runs_extraction(task_nm, sub_nm);
        
        for iRun = run_idx
            run_nm = num2str(iRun);
            switch iRun
                case {2,3}
                    jRun = 1;
                case {5,6}
                    jRun = 2;
            end
            %% run the script
            for iMdl = 1:nModels
                mdl_nm = models{iMdl};
                [ EV_tmp, benefit_tmp, cost_tmp,...
                    F_pred_tmp, perf_pred_tmp, E_pred_tmp,...
                    EV_pred_tmp, benefit_pred_tmp, cost_pred_tmp, E_tmp,...
                    benefit_var_tmp, cost_var_tmp, R2_tmp,...
                    ~, expected_payoff_pred_tmp] =...
                    MS2_GS_plot_model_EV(task_nm, mdl_nm, sub_nm, run_nm);
                EV.(mdl_nm)             = EV_tmp;
                benefit.(mdl_nm)        = benefit_tmp;
                cost.(mdl_nm)           = cost_tmp;
                F_pred.(mdl_nm)         = F_pred_tmp;
                perf_pred.(mdl_nm)      = perf_pred_tmp;
                EV_pred.(mdl_nm)        = EV_pred_tmp;
                benefit_pred.(mdl_nm)   = benefit_pred_tmp;
                cost_pred.(mdl_nm)      = cost_pred_tmp;
                benefit_var.(mdl_nm)    = benefit_var_tmp;
                cost_var.(mdl_nm)       = cost_var_tmp;
                R2.(mdl_nm)             = R2_tmp;
                mR2.(task_nm).(mdl_nm)(jRun,iS) = R2_tmp;
                expected_payoff.(mdl_nm) = expected_payoff_pred_tmp;
                E_pred.(mdl_nm) = E_pred_tmp;
                E.(mdl_nm)     	= E_tmp;
            end % model
            
            %% save the data
            save_folder = [root, filesep sub_nm, filesep 'fMRI_analysis' filesep];
            file_nm = [save_folder, 'GS_model_sub',subid,'_',task_nm,'_run',run_nm,'.mat'];
                        save(file_nm,'EV','benefit','cost',...
                            'F_pred','perf_pred','E_pred','EV_pred','benefit_pred','cost_pred',...
                            'E','benefit_var','cost_var','R2','expected_payoff');
            
            %% clear after each run to avoid confusions
            clear('EV','cost','benefit',...
                'F_pred','perf_pred','E_pred','EV_pred','benefit_pred','cost_pred',...
                'E','benefit_var','cost_var','expected_payoff');
        end % run
        
    end % subject
end % task

for iMdl = 1:nModels
    mdl_nm = models{iMdl};
    avgGrip = mean(mR2.grip.(mdl_nm),1,'omitnan');
    avgStroop = mean(mR2.stroop.(mdl_nm),1,'omitnan');
    mR2.mean.grip.(mdl_nm) = mean(avgGrip,2);
    mR2.mean.stroop.(mdl_nm) = mean(avgStroop,2);
end

end % function