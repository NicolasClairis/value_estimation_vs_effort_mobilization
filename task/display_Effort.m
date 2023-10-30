function [] = display_Effort( WindowPtr, xcenter, ycenter, Levels, trialValence, Physical_Or_Mental)
% Adaptation of the Emotional grip force task (Schmidt et al.,
% JNeurosci2009; Clery-Melin et al., Plos One, 2011) fot Matlab and
% PsychToolbox.
%
% Display_Effort draws the vertical effort scale, and an orange cursor indicating the
% force level.
%
%
% INPUTS
% WindowPtr: window where to display graph
% xcenter, ycenter: coordinates of the center of the screen and of the graph to display
% Levels: structure with data (actual value, max value obtained for this
% trial and incentives, all optional)
%       - Levels.actualLevel : actual level of performance (force for ex.)
%       - Levels.actualFmax : actual max performance (for this trial for
%       ex.)
%       - Levels.incentive: max gain/loss value (for this trial) =>
%       displayed on y axis
% trialValence: is it a Gain or a Loss trial? depending on which scale will
% or won't be reverted
% Physical_Or_Mental: structure with information about each task:
%       - Physical_Or_Mental.Phys_Or_Ment: physical (1) or mental (2) effort scale?
%           - Physical_Or_Mental.Mental.list: list of numbers for the
%           stroop
%           - Physical_Or_Mental.Mental.fontSize: font size for each item
%           of the stroop
%       - Physical_Or_Mental.Mental.error: error trial number in case it
%       happened, ignored if empty
%
%
%
% written by Raphael Le Bouc - May 2013.
% Adapted by Nicolas Clairis - February 2017 for fMRI + display of negative
% values for loss trials + scale which adapts to screen size instead of raw
% pixel values

% color parameters
red = [255 0 0];
weakRed = [100 0 0];
white = [255 255 255];
orange = [255 153 0];
% screen coordinates for effort scale
leftScaleLimit      = 3*xcenter/4; % left limit of the scale
rightScaleLimit     = 5*xcenter/4; % right limit of the scale
bottomScaleLimit    = 3*ycenter/2; % bottom limit of the scale
topScaleLimit       = ycenter/2; % upper limit of the scale
graphSize = bottomScaleLimit - topScaleLimit;
% screen coordinates for bar at the center
leftBarLimit = xcenter - xcenter/12;
rightBarLimit = xcenter + xcenter/12;
% size and coordinates of half of the effort scale
yMetrics = ycenter/2;
% distance between graduations
bigGrad = yMetrics/5;
smallGrad = bigGrad/4;

%% draw the scale (vertical bar)
Screen('DrawLine', WindowPtr, white, leftScaleLimit, topScaleLimit, leftScaleLimit, bottomScaleLimit, 3);

%% draw incentive ticks on the y-axis
if trialValence == 1 % gain trials
    tickValence = '';
elseif trialValence == -1 % loss trials
    tickValence = '-';
end
if isfield(Levels,'incentive') == 1
    incentive = Levels.incentive;
    Screen('TextSize', WindowPtr, 16);
    if trialValence == 1
        yaxis = (-yMetrics):bigGrad:yMetrics;
    elseif trialValence == -1
        yaxis = yMetrics:(-bigGrad):(-yMetrics);
    end
    for ntick = 1:11
        tick = [tickValence, num2str(incentive*(ntick-1)/10)];
        [w,h] = RectSize(Screen('TextBounds', WindowPtr, tick));
        Screen('DrawText', WindowPtr, tick, leftScaleLimit-10-w, ycenter-yaxis(ntick)-h/2, white);
    end
end

%% for mental effort, also add the pairs of symbols on the right of the y axis
if Physical_Or_Mental.Phys_Or_Ment == 2
    list = Physical_Or_Mental.Mental.list;
    fontSize = Physical_Or_Mental.Mental.fontSize;
    % center position of stroop numbers
    xStroop = 6*xcenter/4;
    % space between 2 items of the same pair (btw columns)
    nStroopPairDist = xcenter/12;
    
    % extract error pair in case an error was made
    if isfield(Physical_Or_Mental.Mental,'error') == 1
        nnErrorTrial = Physical_Or_Mental.Mental.error;
    else
        nnErrorTrial = [];
    end
    
    % draw stroop numbers
    Screen('TextSize', WindowPtr, 16);
    yaxis = (-yMetrics):bigGrad:yMetrics;
    for nPairNumber = 2:11
        for iItemOfPair = 1:2
            nn = num2str(list(iItemOfPair,nPairNumber-1));
            fontS = fontSize(iItemOfPair,nPairNumber-1);
            Screen('TextSize', WindowPtr, fontS);
            [w,h] = RectSize(Screen('TextBounds',WindowPtr,nn));
            if exist('nnErrorTrial','var') == 0 || isempty(nnErrorTrial)    % no error made in general
                Screen('DrawText', WindowPtr, nn, xStroop - w/2 + ((iItemOfPair-1)*2-1)*nStroopPairDist, ycenter - yaxis(nPairNumber) - h/2, white);
            elseif exist('nnErrorTrial','var') == 1 && isempty(nnErrorTrial) == 0 % at least one error was made
                if (nPairNumber - 1) ~= nnErrorTrial     % error made but not for this pair => display in white
                    Screen('DrawText', WindowPtr, nn, xStroop - w/2 + ((iItemOfPair-1)*2-1)*nStroopPairDist, ycenter - yaxis(nPairNumber) - h/2, white);
                elseif (nPairNumber - 1) == nnErrorTrial % error made for this pair => display in red
                    Screen('DrawText', WindowPtr, nn, xStroop - w/2 + ((iItemOfPair-1)*2-1)*nStroopPairDist, ycenter - yaxis(nPairNumber) - h/2, red);
                end
            end
        end
    end
end

%% draw the scale (horizontal bars)
for yaxis = -yMetrics:smallGrad:yMetrics
    Screen('DrawLine', WindowPtr, weakRed, leftScaleLimit, (ycenter+yaxis), rightScaleLimit, (ycenter+yaxis), 1);
end

%% update display of values depending on trialValence
switch trialValence
    %% gain trials
    case 1
        for yaxis = (-4*yMetrics/5):(yMetrics/5):yMetrics
            Screen('DrawLine', WindowPtr, white, leftScaleLimit, (ycenter+yaxis), rightScaleLimit, (ycenter+yaxis), 3);
        end
        %% loss trials
    case -1
        for yaxis = (-yMetrics):(yMetrics/5):(4*yMetrics/5)
            Screen('DrawLine', WindowPtr, white, leftScaleLimit, (ycenter+yaxis), rightScaleLimit, (ycenter+yaxis), 3);
        end
end

%% draw the Fmax for this trial until now
if isfield(Levels,'actualFmax') == 1
    FmaxTrialByNow = Levels.actualFmax;
    yFmaxUntilNow = bottomScaleLimit - (FmaxTrialByNow/100)*graphSize;
    Screen('DrawLine', WindowPtr, red, leftScaleLimit, yFmaxUntilNow, rightScaleLimit, yFmaxUntilNow,3);
end

%% draw an orange bar with the actual level of force
if isfield(Levels,'actualLevel') == 1
    actualLevel = Levels.actualLevel;
    yActualLevelBottom = bottomScaleLimit + 10;
    yActualLevelTop = bottomScaleLimit - (actualLevel/100)*graphSize;
    Screen('FillRect', WindowPtr, orange, [leftBarLimit, yActualLevelTop, rightBarLimit, yActualLevelBottom]);
end
end % function