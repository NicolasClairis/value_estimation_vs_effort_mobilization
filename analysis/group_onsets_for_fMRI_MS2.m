function [] = group_onsets_for_fMRI_MS2()
% Launch onsets_for_fMRI_MS2.m for all subjects of MotiScan2 study.
% Then you will have one file per run, to be used for the 1st level GLM.
%
% Developed by Nicolas Clairis
% See also onsets_for_fMRI_MS2.m

subject_id = {'list of subjects'};
NS = length(subject_id);

% loop through subjects
for iSubject = 1:NS
    onsets_for_fMRI_MS2(subject_id{iSubject})
    disp(['Subject ',subject_id{iSubject},' done']);
end

end

