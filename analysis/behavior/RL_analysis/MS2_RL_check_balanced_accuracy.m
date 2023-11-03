root = 'define path';
dataPath = [fullfile(root,'behavior_summary','RL_model'),filesep];
dataSummary = load([dataPath,'RL_model_bis_6_22subs.mat']);
NS = length(model_quality);
bAcc_allSubs = NaN(1,NS);
for iS = 1:NS
    bAcc_allSubs(iS) = dataSummary.model_quality{1,iS}.fit.bacc;
end

[m_bAcc, sem_bAcc] = mean_sem_sd(bAcc_allSubs,2);
disp(['Balanced accuracy is of ',num2str(round(m_bAcc,3)),' +/- ',num2str(round(sem_bAcc,3)),' (mean+/-sem)']);