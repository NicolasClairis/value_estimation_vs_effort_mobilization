function [] = preproc_eyeGazePupil_MotiScan(which_study)
%preproc_eyeGazePupil function to make the preprocessing of the eye gaze
%and pupil data before using analytical scripts to check any correlation
%with the variables of interest.
%
% See also read_eyelink_1_asc_kk.m (This is a modified version of
% read_eye_link_asc.m which is a Fieldtrip function. Feel free to contact
% me if you want that I share it with you)


%% select subjects to preprocess
subject_id = {'enter here list of subjects to check'};
NS_study = length(subject_id);

%% define main acquisition parameters (eye-tracker sampling rate, screen
% size, task names)
if which_study <= 3
    % sampling rate in Hz for CENIR fMRI eye-tracker
    samplingRate = 1000;
    
    % define screen size (in order to remove data when gaze is out)
    x_screen = 1024;
    y_screen = 768;

    taskName = {'RL','grip','stroop'};
    nTasks = length(taskName);
end

disp_gr = 1; % display resume of pupil graph for each run of each subject (1) or not (0)?

%% path definition
root = 'enter path where data is stored here';
list_skipped_files = {};
%% loop across subjects and tasks
for iSubject = 1:NS_study
    subName = subject_id{iSubject};
    % take particular subject number
    if strcmp(subName(3),'_')
        subid = subName(2);
    elseif strcmp(subName(4),'_')
        subid = subName(2:3);
    end
    
    % store data into one folder for each subject (and create it if it
    % does not exist already)
    eyeAnalysis_folder = [root, subName, filesep, 'eye_analysis', filesep];
    if exist(eyeAnalysis_folder,'dir') == 0
        mkdir(eyeAnalysis_folder);
    end
    % same for preprocessing folder
    eyePreproc_folder = [eyeAnalysis_folder, 'preprocessed_files', filesep];
    if exist(eyePreproc_folder,'dir') == 0
        mkdir(eyePreproc_folder);
    end
    
    %% loop through tasks
    for iTask = 1:nTasks
        %% you need a function indicating which run numbers correspond to each task here
        %(taskName input should be 'grip'/'stroop'/'RL'
        % and output will be a (1*n) vector of numbers corresponding to the runs of this task
        runs = task_runs_extraction(taskName{iTask}, subName);
        
        %% loop through runs
        for iRun = runs
            runName = num2str(iRun);
            
            % extract raw eye signal
            eyeDataFolder = [root, subName, filesep, 'eye_data', filesep];
            
            % extract file name
            fileName = ['MS2s',subid,'r',runName,'.asc'];
            full_fileName = [eyeDataFolder, fileName];
            
            if ~ismember(fileName, list_skipped_files)
                
                % extract raw signal of eyes
                asc = read_eyelink_1_asc_kk(full_fileName);
                
                %% check TTL and filter first wrong TTLs in corresponding files
                % if no TTL or problem with TTL too big to be fixed, the run/subject should
                % be abandoned and thus there is no use in spending time with this run/subject
                
                % extract TTL from asc file
                TTL = cell2mat(struct2cell(asc.button)');
                TTL = TTL(TTL(:,3) == 1,1);
                
                % fix TTL for each run and each subject
                [TTL] = TTL_check_for_eye(fileName, TTL, iStudy, subName, runName);
                
                
                %% use preprocessing script
                [eyeTime, X_eyeCoord, Y_eyeCoord, pupil_diam, z_pupil_diam, nanEpisodes] = preproc_eyeGazePupil(asc, TTL(1), x_screen, y_screen, samplingRate, disp_gr);
                %% save the resulting files
                % save resulting files (if data usable)
                if ~isempty(eyeTime) && ~isempty(X_eyeCoord) && ~isempty(Y_eyeCoord) && ~isempty(pupil_diam) && ~isempty(z_pupil_diam)
                    saveNameFiles = [eyePreproc_folder, 'eye_preproc_', subName,'_',taskName{iTask},'_run',runName,'.mat'];
                    save(saveNameFiles, 'eyeTime','X_eyeCoord', 'Y_eyeCoord', 'pupil_diam', 'z_pupil_diam', 'nanEpisodes');
                    disp([subName,' ',taskName{iTask},' run',runName,' done.']);
                else
                    warning([subName,' ',taskName{iTask},' run',runName,' skipped because data not usable.']);
                    list_skipped_files = [list_skipped_files, fileName];
                end
                % save resulting resume figure
                if disp_gr == 1
                    saveNameFig = [eyePreproc_folder, 'eye_preproc_', subName,'_',taskName{iTask},'_run',runName,'.png'];
                    saveas(gcf,saveNameFig)
                    close all;
                end
                
            end % run filter
        end % run loop
    end % task loop
end % subject loop

%% save files ignored because signal was too bad
save([root, filesep, 'behavior_summary', filesep 'eye' filesep 'list_skipped_runs.mat'], 'list_skipped_files');

end % function end

