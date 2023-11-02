function [] = onsets_for_fMRI_MS2(subject_id)
%onsets_for_fMRI_MS2 function extracts all the usefull onsets and durations
% and behavioral components that may have to be used for the fMRI GLM at 
% the first level.
% Use group_onsets_for_fMRI_MS2 if you want to launch all subjects in a
% raw.
%
% INPUT
% subject_id: "sX_XXyyZZ" subject number + date of experiment
%
% Developed by Nicolas Clairis - 2017
%
% See also group_onsets_for_fMRI_MS2, onsets_for_fMRI_MS2_RL,
% onsets_for_fMRI_MS2_grip, onsets_for_fMRI_MS2_stroop

close all;

%% subject identification
if isempty(subject_id)
    subject_id = input('Identifiant du sujet? sX_yyZZnn','s'); % subject id: X: number, yy: day, ZZ month, nn: year
end

if strcmp(subject_id(3),'_')
%     subj_id = str2double(subject_id(2));
    subid = subject_id(2);
elseif strcmp(subject_id(3),'_') == 0 && strcmp(subject_id(4),'_')
%     subj_id = str2double(subject_id(2:3));
    subid = subject_id(2:3);
else
    warning('Weird bug with your subject identification name. Please correct.')
    return;
end

%% run numbers
% theoretical total number of runs (2 grip, 2 stroop, 3 RL)
nb_runs = 7; % 7 runs/subject

%% working directories
root = 'root path';
subDir = [root subject_id filesep];
behaviorDir = [subDir 'behavior' filesep];
fMRI_Dir = [subDir 'fMRI_analysis' filesep];

%% check correspondance btw number of theoretical runs and particular
% subject

% 2 stars required in old matlab versions but now 1 star should be enough
% fMRI_runnames = ls([behaviorDir, filesep, 'MBB_battery_*','*_onsets_sub',subid,'_sess*','*.mat']);
fMRI_runnames = ls([behaviorDir, 'MBB_battery_*_onsets_sub',subid,'_sess*.mat']);
check_nb_runs = size(fMRI_runnames,1);
if nb_runs == check_nb_runs
    disp('nb of runs ok');
else
    warning(['Mismatch: ',num2str(nb_runs), ' runs in theory but ',num2str(check_nb_runs),' runs in practical']);
    warning(['Please check the behavior directory for the subject ',subject_id])
    return;
end

%% learning
onsets_for_fMRI_MS2_RL(root, behaviorDir, fMRI_Dir, subid, subject_id);
disp(['Subject ',subject_id,' - learning runs done']);


%% Grip
onsets_for_fMRI_MS2_grip(behaviorDir, fMRI_Dir, subid);
disp(['Subject ',subject_id,' - grip runs done']);


%% Stroop
onsets_for_fMRI_MS2_stroop(behaviorDir, fMRI_Dir, subid);
disp(['Subject ',subject_id,' - stroop runs done']);


end % function end