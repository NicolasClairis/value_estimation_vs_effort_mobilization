% gripDisplay
% Benoit Beranger GUI to check grip values displayed on screen as if it was
% an oscilloscope

close all
clear all
clc

%% addpath to grip functions folder in case not there
addpath([pwd filesep 'grip_functions' filesep 'Basic_grip_functions' filesep 'GripCompatFunc']);
% main parameters
DeviceName = 'SerialMBB';
Channels = input('channel 1 or 2 ?');
COM = input('COM 1 or 3');

Handle = InitializeGrip(DeviceName, Channels, COM);
RefreshPeriod = 1/60;

Window              = 10; % secondes
Scope       = nan( round(Window/RefreshPeriod) , 1 );
Time        = RefreshPeriod : RefreshPeriod : Window;
L           = length(Scope);
Freq        = (0:(L/2))/L/RefreshPeriod;



figure;

while ~KbCheck
   
    Scope  = circshift(Scope ,1);
   
    % Read ADC
    [Values, Times] = ReadGripValue(DeviceName, Handle, Channels);
    Scope(1) = Values ;
   
    plot(Time,Scope)
    ylim([0 1024])
   
    drawnow
   
end

CloseGripDevice(DeviceName, Handle)
