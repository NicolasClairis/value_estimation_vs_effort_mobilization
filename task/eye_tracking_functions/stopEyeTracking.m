function[] = stopEyeTracking(subdir,eyeFileName,taskName,subid,runname)
% function based on Benoit Beranger functions to stop the fMRI eye-tracker
% created by Nicolas Clairis - february 2017
% to be used in taskGripRP.m, taskLearning75.m and taskMentalRP.m

% send a final trigger to the eye-tracker to know when the task stopped
    % and estimate global timing
    OpenParPort; % open parallel port
    WriteParPort(255) % send a trigger
    WaitSecs(0.003) % display triggers during 3ms at least (otherwise not detected since sampling rate of eye-tracker = 1ms) (use minium samplingRate x2, usually x3)
    WriteParPort(0) % stop trigger
    CloseParPort; % close parallel port (this line is not mandatory given what Benoit says)
    
    % Stop the record, record it in the results folder and change the name
    % to make a bigger and clearer one
    % check if eye folder exists for this subject, if not create it
    eyeFolder = [subdir,'\eye_data'];
    if exist(eyeFolder,'dir') == 0
        mkdir(eyeFolder);
    end
    Eyelink.StopRecording(eyeFileName,eyeFolder) % stop recording, close file, download the file into the directory
%     movefile([eyeFolder, filesep, eyeFileName,'.edf'],[eyeFolder, filesep, 'MS2_',taskName,'_s',subid,'r',runname,'.edf']) % rename / move file (rename for having a title longer than just 8 characters)
    
end