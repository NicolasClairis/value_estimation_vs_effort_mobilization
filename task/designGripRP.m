function [outcomeValues, valence, training, task_inc, trials, blocks] = designGripRP(nrun)
%% designGripRP
%
% INPUT
% nrun: check whether training (nrun = 0) or actual task run (nrun > 0)
% 
% OUTPUT
% outcomeValues: values of possible rewards
% valence: no idea
% training: no idea
% task_inc: no idea
% trials: no idea
% blocks: no idea
%
% function developed by Nicolas Borderies, improved by Nicolas Clairis
% to be more flexible and adapted for fMRI specifically

%% Incentives
outcomeValues = [0.01, 0.2, 0.5, 1, 5, 20];
nOutcomeValuesNumber = length(outcomeValues);

%% Trial vectors
blocks = 1:10;
trials = 1:12;
trialNumber = [1:numel(blocks)*numel(trials)];
task_inc = [];
for nblock = blocks
    cond = [];
    for i = 1:length(trials)/6
        cond = [cond randperm(6)];
    end
    task_inc = [task_inc mod(cond-1,6) + 1]; % 1cents = 1, 20cents = 2, 50cents = 3, 1euro = 4, 5euros = 5, 20euros = 6
end

valence = mod(ceil(trialNumber/12),2)*2-1 ;

end