function[eyeFileName] = startEyeTracking(subid,runname)
% function to launch eye-tracker, based on Benoit Beranger functions
% made by Nicolas Clairis - february 2017

try
    Eyelink.Initialize; % Start connexion of stim computer with eye-tracking computer
    % if initialize doesn't work, use Eyelink.ForceShutDown and re-do
    % Eyelink.Initialize;
catch err1
    try % if EyeLink.Initialize, tries to check whether connection is working properly with EyeLink.IsConnected
        Eyelink.IsConnected
    catch err2 % if this doesn't work either, catch error
    end
end

% load eye-tracking parameters (number of points for calibration,
% left/right eye, etc.)
Eyelink.LoadParameters; % Load default parameters (record TTL, right eye, etc.) for each run
% (because the default parameters of the eye-tracker don't take all
% these vars into account and they reset everytime after a (re-)calibration)

% define name of the file where you will store eye data
eyeFileName = ['MS2s',subid,'r',runname];

% start recording in this file
Eyelink.StartRecording(eyeFileName), % open file, start recording (BE CAREFUL no more than 8 characters for eye-tracker)

end