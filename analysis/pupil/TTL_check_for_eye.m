function [TTL] = TTL_check_for_eye(fileName, TTL, subjectName, runName)
% This script cleans run/run and subject/subject TTL data received on the
% eye-tracker.
%
% Checks TTL and filter first (or last) wrong TTLs in corresponding files
% if no TTL or problem with TTL too big to be fixed, the run/subject should
% be abandoned and thus there is no use in spending time with this
% run/subject.
%
% This bug happened in particular when the eye-tracker was launched while
% the fMRI was still running on a previous run aimed to correct for AP-PA
% deformations or when it was stopped while the functional run ended and a
% new run for this purpose was launched thus the eye-tracker received some
% triggers that are independent of the functional run and which should in
% consequence be excluded from the analysis to avoid mis-aligning the data.
%
% INPUTS
% fileName: string containing the name of the run and the subject currently
% analyzed
%
% TTL: vector containing all the triggers received by the eye-tracker
% during the scan
%
% subjectName: string with the full name of the subject currently analyzed
%
% runName: string containing the number of the run currently analyzed
%


%% indicate fMRI TR in ms
TR = 1100;

%% tolerance range for TR variation
diff_TR = 20; % (ms)
TRmin = TR - diff_TR;
TRmax = TR + diff_TR;

%% check if TTL ok compared to TR
TTLok = diff(TTL);

%% correct TTL list for each run of each subject
if (sum(TTLok < TRmin) > 0 || sum(TTLok > TRmax) > 0)
    disp('problem with TTL/TR : please check the data');
    disp(['sujet ',subjectName,', run ',runName]);
    return;
else
    disp(['TTL ok in ',fileName]);
end

end

