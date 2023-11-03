function [ mR2, semR2, sdR2 ] = MS2_RL_choice_model_quality( model_n )
%[ mR2, semR2, sdR2 ] = MS2_RL_choice_model_quality( model_n )
%MS2_RL_choice_model_quality extracts average R² for the model and the
%subjects selected in the input.
%
% INPUTS
% model_n: model number
%
% OUTPUTS
% mR2: mean R2 of the model across subjects
%
% semR2: sem of the R2 across subjects
%
% sdR2: standard deviation of the R2 across subjects
%

%% working directories
root = 'define path here';

%% subject identification
subject_id = {'enter list of subjects here'};
NS = length(subject_id);

%% load data
RL_data_folder = [fullfile(root, 'behavior_summary','RL_model'), filesep];
RL_mdl_quality = getfield(load([RL_data_folder,'RL_model_bis_',num2str(model_n),'_',num2str(NS),'subs.mat'],...
    'model_quality'),'model_quality');

R2 = NaN(1,NS);
for iS = 1:NS
    R2(iS) = RL_mdl_quality{1,iS}.fit.R2;
end
mR2     = mean(R2,2,'omitnan');
semR2 = sem(R2, 2);
sdR2 = std(R2);

end % function