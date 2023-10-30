% function [] = taskGripRP(subid,nsession)
% MBB battery
% Grip force task (Reward and punishments)
%
% Raphael Le Bouc - april 2014
% Re-adapted for fMRI - February 2017 - Nicolas Clairis

clear all; close all; clc;

testing_script = 0; % parameters of screen to check on my own computer

%% add eye-tracking functions folder
addpath(genpath('eye_tracking_functions'));
%% addpath to grip functions folder in case not there
addpath([pwd filesep 'grip_functions' filesep 'Basic_grip_functions' filesep 'GripCompatFunc']);

%% Identification
% identification_batmotiv; % i never liked this NicoB function
subid = input('subject identification number? ','s');
subject = str2double(subid);
nbTotalRuns = 8; % total number of runs (all tasks considered)
nrun = input(['Run number (training:0;task:1-',num2str(nbTotalRuns),') ?']);
runname = num2str(nrun);
if nrun == 0 % training outside fMRI
    IRM = 0;
else % task inside fMRI
    IRM = 1;
end
nbTotalSessions = 8; % 8 runs (4 learning75, 2 gripRP, 2 mentalRP
taskName = 'gripRP';

% identification_batmotiv: define name of the main results file = resultname
[resultname] = identification_batmotiv(taskName, subid, runname);

% setDir: create directory with results for all subjects (resultdir) if did
% not exist + folder for this subject (subdir) + folder where behavioral
% results stored for this subject (behaviordir)
% + define "root" folder
[root, resultdir, subdir, behaviordir, fMRIScansDir] = setDir(subid, IRM);

%% check if run already made or not (for the task only, not training)
[MadeOrNot] = checkIfRunMade(behaviordir, subid, nrun, nbTotalSessions);
if MadeOrNot == 1
    return; % stop task if run number already made
else
    clear('MadeOrNot')
end

%% Eye-tracking
% open and calibFmaxrate eye-tracker
if IRM == 1 % no eye-tracking outside fMRI
    addEyeTrack = input('Eye-tracking? 1 = yes / 0 = no');
    if addEyeTrack == 1
        [eyeFileName] = startEyeTracking(subid,runname);
    end
else
    addEyeTrack = 0;
end


%% Grip Configuration
% % IOPort Reset
% IOPort('CloseAll');
% initialize grip
GripPort = 'SerialMBB';
GripChannel = input('Do you want to use channel (1) or channel (2)? Please press corresponding number.'); % define whether you use channel 1 or channel 2 from the transducer
% GripChannel = 1;
if IRM == 0 || testing_script == 1 % outside fMRI scanner use default COM of the computer (COM1 for windows)
    COM = 1;
elseif IRM == 1 % inside fMRI scanner use COM3 because COM1 already used by button response device
    COM = 3;
end
[Handle] = InitializeGrip(GripPort, GripChannel, COM);

%% Screen configuration
[x, y, window, baselineTextSize] = ScreenConfiguration(IRM, testing_script);

%% design configuraton
% common parameters between gripRP and mentalRP
nCalibrationTrials = 3; % number of trials used for initial calibration of force
% Incentives
outcomeValues = [0.01, 0.2, 0.5, 1, 5, 20]; % positive = gain, same as negative = loss
nConditions = length(outcomeValues);
if nrun == 0 % training
    nBlocks = 1; % number of blocks (1 block = all conditions (rewards and losses) seen once)
elseif nrun > 0
    nBlocks = 5; % number of blocks (1 block = all conditions (rewards and losses) seen once)
end
nTrialsPerBlock = nConditions*2; % number of trials per block (2 per condition since that each condition may appear as a gain or a loss)
blocks = 1:nBlocks; % vector of blocks
totalTrial = nBlocks*nTrialsPerBlock; % total number of trials in this run
trialNumber = 1:totalTrial; % vector with all trials

% attribute which incentive goes for which trial
task_inc = nan(1,totalTrial); % vector containing the incentive corresponding to each trial
trialValence = nan(1,totalTrial); % valence is a vector indicating if it is a Rewarding or a Punishment trial
for iblock = blocks % loop through blocks
    blockDistribution = randperm(nTrialsPerBlock);
    blockIndex = (1 + (iblock - 1)*nTrialsPerBlock):(iblock*nTrialsPerBlock);
    trialValence(blockIndex) = (blockDistribution > nConditions).*(-1) +  (blockDistribution <= nConditions); % valence is negative for values higher than nConditions, positive for lower values
    task_incBlockDistribution = blockDistribution;
    task_incBlockDistribution(task_incBlockDistribution > nConditions) = task_incBlockDistribution(task_incBlockDistribution > nConditions) - nConditions; % convert blockDistribution to values that may be used without changing the original form of the script
    task_inc(blockIndex) = task_incBlockDistribution; % only values ranking from 1 to nConditions
end
% total of virtual money earned by the subject in the task
if nrun >= 3
    totalEarned = 0;
elseif nrun >= 4 && nrun <= 6
    if exist([behaviordir,'\global_sub_',subid,'_session_2_',taskName,'.mat'])
        load([behaviordir,'\global_sub_',subid,'_session_2_',taskName,'.mat'],'totalGain')
        totalEarned = totalGain(2);
    elseif exist([behaviordir,'\global_sub_',subid,'_session_3_',taskName,'.mat'])
        load([behaviordir,'\global_sub_',subid,'_session_3_',taskName,'.mat'],'totalGain')
        totalEarned = totalGain(3);
    end
end

% variables grip-specific
PercFmax = 100/75; % calibrate task so that real Fmax is at (1/PercFmax)% of actual Fmax displayed on screen
Physical_Or_Mental.Phys_Or_Ment = 1; % physical effort (for Display_Effort function)
Fmax = zeros(1,totalTrial); % Fmax used for each trial (calibrated during calibration, but may change during the task)
FmaxDisp = zeros(1,totalTrial); % Fmax used for each trial which is actually higher than the real one
peakForce = zeros(1,totalTrial); % store Fmax (peak force) exerted at each trial
sumforce = zeros(1,totalTrial);
perf = zeros(1,totalTrial);
gain = zeros(1,totalTrial);
feedback = zeros(1,totalTrial);
% calibration grip data
baseline_Calibration = zeros(1,nCalibrationTrials); % store baseline used for each calibration calibration trial in this variable
baselineGripData_Calibration = {}; % store baseline grip values and timings (in all moments when the grip should not be used but may show some variations)
baselineGripData_Calibration.instructions = []; % grip values for each calibration trial during instructions
gripData_Calibration = {}; % store grip values during the action phase of the calibration trial
gripData_Calibration.grip = []; % grip (raw) values for each calibration trial during action phase
gripData_Calibration.time = []; % grip timings for each calibration trial during action phase
gripData_Calibration.gripB = []; % grip values for each calibration trial with baseline correction during action phase
% task grip data
baseline = nan(1,totalTrial); % store baseline used for each trial in this variable
baselineGripData = {}; % store baseline grip values and timings (in all moments when the grip should not be used but may show some variations)
baselineGripData.fixationCross = []; % grip values for each trial during fixation cross
baselineGripData.incentive = []; % grip values for each trial during incentive display
baselineGripData.feedback = []; % grip values for each trial during feedback
gripData = {}; % store grip values during the action phase of the trial
gripData.grip = []; % grip (raw) values for each trial during action phase
gripData.time = []; % grip timings for each trial during action phase
gripData.gripB = []; % grip values for each trial with baseline correction during action phase
gripData.level = []; % level of force as it is displayed on screen at any moment of the trial
gripData.FmaxUntilNow = []; % level of force max for each trial as it is displayed on screen

%% Loading images
[ picCross, rectCross, pic_inc, rect_inc] = load_incentive_images(window,x,y, root, nConditions);

%% Generator reset
rand('state',sum(100*clock));

%% Time variables
instruction_T = 3;
calibStartMessage_T = 1; % "C'est parti" message duration
instructionCalibration_T = 2;
cross_block_T = 0.5;
block_T = 5;
cross_ITI_T = 0.5;
% final cross
if nrun == 0 % shorter for training
    ending_T = 5;
else
    ending_T = 10;
end
release_Wait = 0.5; % duration of cross after "please release"
PressTime_Wait = 5; % max time for waiting a response
% incentive jitter display between 1 to 3.95 seconds
jitters = 1:0.05:3.95;
varIncPerms = randperm(60);
varInc_Wait = jitters(varIncPerms);
% feedback jitter display between 1 to 3.95 seconds
varFdbkPerms = randperm(60);
varFdbk_Wait = jitters(varFdbkPerms);
% display of total gain at the end of the task
totalGain_Wait = 2;

%% events onsets and durations definition
% timings present/possibly present in all trials
listB = {'cross_ITI',... % fixation cross for ITI
    'displayIncentive',... % display options on screen
    'displayEffortScale',... % display effort scale on screen
    'gainFeedback',... % feedback win/lose/neutral
    'lossFeedback',...
    'neutralFeedback',...
    'feedback',...
    'pleaseRelease',... % for cases where subject presses some button during the ITI
    'crossRelease',... % fixation cross after a "Please Release" message
    'tooSlowTrial_displayIncentive',... % too slow trials
    'tooSlowTrial_displayEffortScale',...
    'tooSlow',...
    'tooSlowTrial_feedback'};

for l = 1 : length(listB)
    onset.(listB{l}) = NaN(totalTrial,1);
    duration.(listB{l}) = NaN(totalTrial,1);
end
onset.TTL = []; % record TTL in that variable also to have all timings here
onset.calibStart = []; % record when calibration "C'est parti" message starts
duration.calibStart = [];
onset.instruction = [];
duration.instruction = [];
onset.cross = []; % fixation cross when no ITI (start of task, end of task)
duration.cross = [];
onset.totalGain = []; % display of total gain at the end of the run
duration.totalGain = [];

%% keyboard keys configuration
if IRM == 0
    % key configuration
    KbName('UnifyKeyNames');
    key.escape = KbName('escape');
elseif IRM == 1
    % fMRI key configuration
    KbName('UnifyKeyNames');
    key.escape = KbName('escape');
    
    % TTL trigger
    dummy_scan = 4; % 4 trash volumes at fMRI start
    trigger = 53; % TTL encoded as the "5([" key in the keyboard
end

%==========================================================================
% CALIBRATION
%==========================================================================

% give/remind instructions
[timenow1] = disp_instructions(window, x,y, IRM, taskName);
onset.instruction = timenow1;
if IRM == 0
    KbWait;
elseif IRM == 1
    WaitSecs(instruction_T);
end
timenow2 = GetSecs;
duration.instruction = timenow2 - timenow1;

% outside fMRI wait for a button press before starting calib
if IRM == 0
    wrapat = 30; % nb of characters bf breaking line
    vspacing = 2; % space between lines
    DrawFormattedText(window,'D''abord quelques essais pour s''echauffer','center','center',[255 255 255],wrapat,[],[],vspacing)
    Screen('TextSize', window, 20);
    textPleasePress = 'appuyer sur une touche pour continuer';
   [width,hight] = RectSize(Screen('TextBounds',window,textPleasePress));
    Screen('DrawText',window,textPleasePress,x-width/2,9*y/5,[100 100 100]); 
    % onset
    [~,timenow1,~,~,~] = Screen(window,'Flip');
    onset.instruction = [onset.instruction; timenow1];
    KbWait;
    timenow2 = GetSecs;
    duration.instruction = [duration.instruction; timenow2 - timenow1];
    Screen('TextSize', window, baselineTextSize);
end

% calibFmax_made_or_not = input('Calibration already made (1) or not (0)?');
% if calibFmax_made_or_not == 0

% calibration Fmax initial value
calibFmax = 0;

% signal that it starts
Screen('TextSize', window, 60);
[width,hight] = RectSize(Screen('TextBounds',window,'C''est parti!'));
Screen('DrawText',window,'C''est parti!',x-width/2,y-hight/2-100,[255 0 0]);
Screen('TextSize', window, 20);
[~,timenow1,~,~,~] = Screen(window,'Flip');
onset.calibStart = timenow1;
WaitSecs(calibStartMessage_T);
timenow2 = GetSecs;
duration.calibStart = timenow2 - timenow1;

% calibration
for calibTrial = 1:nCalibrationTrials
    
    %% Display instruction
    instructime = GetSecs;
    iBaselineGrip = 0;
    onset.instruction_Calibration = timenow1;
    while GetSecs <= instructime + instructionCalibration_T
        iBaselineGrip = iBaselineGrip + 1;
        Screen('TextSize', window, baselineTextSize);
        [width,hight] = RectSize(Screen('TextBounds',window,[num2str(calibTrial) 'e essai']));
        Screen('DrawText',window,[num2str(calibTrial) 'e essai'],x-width/2,y-hight/2-100,[255 0 0]);
        Screen('TextSize', window, 30);
        [width,hight] = RectSize(Screen('TextBounds',window,'Vous devez serrer le plus fort possible'));
        Screen('DrawText',window,'Vous devez serrer le plus fort possible',x-width/2,y-hight/2,[255 0 0]);
        Screen('TextSize', window, baselineTextSize);
        Screen(window,'Flip');
        % extract grip values during fixation cross to extract baseline
        [grip,Tgrip] = ReadGripValue(GripPort,Handle,GripChannel);
        baselineGripData_Calibration.instructions.grip{calibTrial}(iBaselineGrip) = grip;
        baselineGripData_Calibration.instructions.time{calibTrial}(iBaselineGrip) = Tgrip;
    end
    % duration starting calibration
    timenow2 = GetSecs;
    duration.instruction_Calibration =  timenow2 - timenow1;
    % extract mean baseline for this calibration trial
    baseline_Calibration(calibTrial) = mean(baselineGripData_Calibration.instructions.grip{calibTrial});
    
    %% Display Effort score
    calibFmaxrationtime = GetSecs;
    keep = 0;
    igripCalib = 0;
    while GetSecs <= (calibFmaxrationtime + PressTime_Wait)
        igripCalib = igripCalib + 1;
        [grip,Tgrip] = ReadGripValue(GripPort,Handle,GripChannel);
        gripData_Calibration.grip{calibTrial}(igripCalib) = grip;
        gripData_Calibration.time{calibTrial}(igripCalib) = Tgrip;
        % correct for baseline
        gripBCalib = grip - baseline_Calibration(calibTrial);
        gripData_Calibration.gripB{calibTrial}(igripCalib) = gripBCalib;
        if gripBCalib < 0 % if for some strange reason grip is lower than baseline => threshold at 0 (weird to have negative values..)
            gripBCalib = 0;
        end
        keep = max([keep gripBCalib]);
        
        Screen('TextSize', window, 60);
        [width,hight] = RectSize(Screen('TextBounds',window,'Go!'));
        Screen('DrawText',window,'Go!',x-width/2,y-hight/2-100,[255 0 0]);
        [width,hight] = RectSize(Screen('TextBounds',window,num2str(round(gripBCalib))));
        Screen('DrawText',window,num2str(round(gripBCalib)),x-width/2,y-hight/2,[255 255 255]);
        
        Screen('TextSize', window, 40);
        [width,hight] = RectSize(Screen('TextBounds',window,'Record'));
        Screen('DrawText',window,'Record',x-width/2+300,y-hight/2-100,[175 175 175]);
        Screen('TextSize', window, 60);
        [width,hight] = RectSize(Screen('TextBounds',window,num2str(round(max([calibFmax keep])))));
        Screen('DrawText',window,num2str(round(max([calibFmax keep]))),x-width/2+300,y-hight/2,[175 175 175]);
        
        Screen(window,'Flip');
        WaitSecs(0.02); % useless WaitSecs?
    end
    
    % Record the maximal force
    calibFmax = max([calibFmax keep]);
end
    
% elseif calibFmax_made_or_not == 1
%     disp(['Loading Calibration parameters for subject ',subid],'calibFmax');
%     %     load([,subid,]);
% end

Fmax(1) = calibFmax;
FmaxDisp(1) = calibFmax*PercFmax;

%==========================================================================
% TASK
%==========================================================================

%% start recording all TTL
% KbQueue records every time an fMRI trigger is received.
% Note that KbCheck or KbWait commands will not be impacted by the fact
% that KbQueue only checks for TTL

%% waiting and recording first TTL for fMRI
if IRM == 1
    % add a cross at the beginning (better for eye tracking)
    Screen('DrawTexture',window,picCross,[],rectCross);
    % onset
    [~,timenow1,~,~,~] = Screen(window,'Flip');
    onset.cross = timenow1;
    
    %% TRIGGER RAPHAEL
    if nrun > 0 % useless for training
        next = 0;
        TTL = []; % TTL TIMES
        % wait dummy_scan number of volumes before starting the task
        while next < dummy_scan
            [keyisdown, T0IRM, keycode] = KbCheck;
            
            if keyisdown == 1 && keycode(trigger) == 1
                if next == 0
                    T0 = T0IRM;
                end
                next = next+1;
                disp(next);
                TTL = [TTL; T0IRM];
                while keycode(trigger) == 1
                    [keyisdown, T, keycode] = KbCheck;
                end
            end
        end
        
        % record all subsequent TTL in the whole task
        keysOfInterest = zeros(1,256);
        keysOfInterest(trigger) = 1;
        KbQueueCreate(0,keysOfInterest); % checks TTL only
        KbQueueStart; % starts checking
    end
    
    % duration starting cross
    timenow2 = GetSecs;
    duration.cross = timenow2 - timenow1;
elseif IRM == 0
    % add an instruction between calib and task
    wrapat = 30; % nb of characters bf breaking line
    vspacing = 2; % space between lines
    Screen('TextSize', window, 40);
    DrawFormattedText(window,'On continue des que vous etes pret(e)','center','center',[255 255 255],wrapat,[],[],vspacing)
    Screen('TextSize', window, 20);
    [width,hight] = RectSize(Screen('TextBounds',window,textPleasePress));
    Screen('DrawText',window,textPleasePress,x-width/2,9*y/5,[100 100 100]);
    % onset
    [~,timenow1,~,~,~] = Screen(window,'Flip');
    onset.instruction = timenow1;
    KbWait;
    timenow2 = GetSecs;
    duration.instruction = timenow2 - timenow1;
    Screen('TextSize', window, baselineTextSize);
end

%% TASK
block = 0; % block index
for ntrial = 1:totalTrial
    
    %% fixation cross
    Screen('DrawTexture',window,picCross,[],rectCross);
    [~,timenow1,~,~,~] = Screen(window,'Flip');
    onset.cross_ITI(ntrial) = timenow1;
    dur = 0;
    iBaselineGrip = 0;
    while dur <= cross_ITI_T
        iBaselineGrip = iBaselineGrip + 1;
        % keep displaying cross
        Screen('DrawTexture',window,picCross,[],rectCross);
        [~,timenow2,~,~,~] = Screen(window,'Flip');
        % extract grip values during fixation cross to extract baseline
        [grip,Tgrip] = ReadGripValue(GripPort,Handle,GripChannel);
        baselineGripData.fixationCross.grip{ntrial}(iBaselineGrip) = grip;
        baselineGripData.fixationCross.time{ntrial}(iBaselineGrip) = Tgrip;
        % check duration
        dur = timenow2 - timenow1;
    end
%     WaitSecs(cross_ITI_T);
    % duration
%     timenow2 = GetSecs;
    duration.cross_ITI(ntrial) = dur;
    
    %% Display incentive
    Screen('DrawTexture',window,pic_inc{task_inc(ntrial), 3-(trialValence(ntrial)+3)/2},[],rect_inc{task_inc(ntrial)});
    incentiveLeft = rect_inc{task_inc(ntrial)}(1);
    incentiveBottom = rect_inc{task_inc(ntrial)}(4);
    Screen('TextSize', window, 32);
    if trialValence(ntrial) == 1
        Screen('DrawText', window, 'A GAGNER', incentiveLeft, incentiveBottom, [255 255 255], [], 0);
    elseif trialValence(ntrial) == -1
        Screen('DrawText', window, 'A PERDRE', incentiveLeft, incentiveBottom, [255 255 255], [], 0);
    end
    [~,timenow1,~,~,~] = Screen(window,'Flip');
    onset.displayIncentive(ntrial) = timenow1;
    dur = 0;
    iBaselineGrip = 0;
    while dur <= varInc_Wait(ntrial)
        iBaselineGrip = iBaselineGrip + 1;
        % keep displaying incentive
        Screen('DrawTexture',window,pic_inc{task_inc(ntrial), 3-(trialValence(ntrial)+3)/2},[],rect_inc{task_inc(ntrial)});
        Screen('TextSize', window, 32);
        if trialValence(ntrial) == 1
            Screen('DrawText', window, 'A GAGNER', incentiveLeft, incentiveBottom, [255 255 255], [], 0);
        elseif trialValence(ntrial) == -1
            Screen('DrawText', window, 'A PERDRE', incentiveLeft, incentiveBottom, [255 255 255], [], 0);
        end
        [~,timenow2,~,~,~] = Screen(window,'Flip');
        % extract grip values during fixation cross to extract baseline
        [grip,Tgrip] = ReadGripValue(GripPort,Handle,GripChannel);
        baselineGripData.incentive.grip{ntrial}(iBaselineGrip) = grip;
        baselineGripData.incentive.time{ntrial}(iBaselineGrip) = Tgrip;
        % check duration
        dur = timenow2 - timenow1;
    end
    %     WaitSecs(varInc_Wait);
    % duration
%     timenow2 = GetSecs;
    duration.displayIncentive(ntrial) = dur;
    
    
    %% Display Physical Effort scale
    % extract baseline
    if ntrial == 1 % for the first trial: only take fixation cross grip baseline
        baseline(ntrial) = mean(baselineGripData.fixationCross.grip{ntrial});
    else % afterwards: extract baseline during fixation cross + during the display of the feedback
        baseline(ntrial) = mean([baselineGripData.fixationCross.grip{ntrial}, baselineGripData.feedback.grip{ntrial - 1}]);
    end
    % display scale
    igrip = 0;
    onset.displayEffortScale(ntrial) = GetSecs;
    trialTime = 0;
    while trialTime <= PressTime_Wait
        % grip extract current value
        igrip = igrip + 1;
        [grip,Tgrip] = ReadGripValue(GripPort,Handle,GripChannel);
        gripData.grip{ntrial}(igrip) = grip;
        gripData.time{ntrial}(igrip) = Tgrip;
        % correct for baseline
        gripB = grip - baseline(ntrial);
        gripData.gripB{ntrial}(igrip) = gripB;
        % extract maximum value (corrected for baseline) until now for this
        % trial in order to display that value on screen
        gripBUntilNowMax = max(gripData.gripB{ntrial});
        
        % extract ACTUAL level of force relative to calibrated Fmax
        if gripB > 0 && gripB <= FmaxDisp(ntrial)
            gripData.level{ntrial}(igrip) = (gripB/FmaxDisp(ntrial))*100; % convert F in percentage of Fmax
        elseif gripB > FmaxDisp(ntrial) % threshold high level to 100%
            gripData.level{ntrial}(igrip) = 100;
        elseif gripB <= 0 % threshold low level to 0%
            gripData.level{ntrial}(igrip) = 0;
        end
        % extract Fmax for this trial until now relative to calibrated Fmax
        if gripBUntilNowMax > 0 && gripBUntilNowMax <= FmaxDisp(ntrial)
            gripData.FmaxUntilNow{ntrial}(igrip) = (gripBUntilNowMax/FmaxDisp(ntrial))*100; % convert F in percentage of Fmax
        elseif gripBUntilNowMax > FmaxDisp(ntrial) % threshold high level to 100%
            gripData.FmaxUntilNow{ntrial}(igrip) = 100;
        elseif gripBUntilNowMax <= 0 % threshold low level to 0%
            gripData.FmaxUntilNow{ntrial}(igrip) = 0;
        end
        % display feedback to the subject about the force exerted
        Screen('DrawTexture', window, pic_inc{task_inc(ntrial), 3-(trialValence(ntrial)+3)/2}, [], rect_inc{task_inc(ntrial)});
        Screen('TextSize', window, 32);
        if trialValence(ntrial) == 1
            Screen('DrawText', window, 'A GAGNER', incentiveLeft, incentiveBottom, [255 255 255], [], 0);
        elseif trialValence(ntrial) == -1
            Screen('DrawText', window, 'A PERDRE', incentiveLeft, incentiveBottom, [255 255 255], [], 0);
        end
        % parameters for display_Effort
        Levels.actualLevel = gripData.level{ntrial}(igrip);
        Levels.incentive = outcomeValues(task_inc(ntrial));
        Levels.actualFmax = gripData.FmaxUntilNow{ntrial}(igrip);
        display_Effort(window, x, y, Levels, trialValence(ntrial), Physical_Or_Mental);
        [~,timenow3,~,~,~] = Screen(window,'Flip');
        trialTime = timenow3 - timenow2;
    end
    duration.displayEffortScale(ntrial) = trialTime;
    
    %% Data to save
    peakForce(ntrial) = max(gripData.gripB{ntrial});
    sumforce(ntrial) = sum(gripData.gripB{ntrial});
    % extract performance for this trial but thresholding at min: 0%/max: 100%
    if peakForce(ntrial) > 0 && peakForce(ntrial) <= FmaxDisp(ntrial)
        perf(ntrial) = (peakForce(ntrial)/FmaxDisp(ntrial))*100; % convert F in percentage of FmaxDisp
    elseif peakForce(ntrial) > FmaxDisp(ntrial) % threshold high level to 100%
        perf(ntrial) = 100;
    elseif peakForce(ntrial) <= 0 % threshold low level to 0%
        perf(ntrial) = 0;
    end
    % redefine Fmax if Fmax exerted in this trial is superior to the
    % previously recorded value 
    if ntrial < totalTrial % useless for the last trial, since no other trial afterwards
        if peakForce(ntrial) > Fmax(ntrial) % if peak force > previous Fmax => update it for the next trial
            Fmax(ntrial + 1) = peakForce(ntrial);
            FmaxDisp(ntrial + 1) = peakForce(ntrial)*PercFmax;
        elseif peakForce(ntrial) <= Fmax(ntrial) % if peak force <= previous Fmax => keep previous Fmax for the next trial
            Fmax(ntrial + 1) = Fmax(ntrial);
            FmaxDisp(ntrial + 1) = FmaxDisp(ntrial);
        end
    end
    % calculate gain depending on trialValence
    switch trialValence(ntrial)
        case 1
            gain(ntrial) = perf(ntrial)/100*outcomeValues(task_inc(ntrial));
        case -1
            gain(ntrial) = (100 - perf(ntrial))/100*-outcomeValues(task_inc(ntrial));
            gain(ntrial) = gain(ntrial)*(gain(ntrial)<0);   % in the loss condition, the best outcome is zero.
    end
    previoustotal = totalEarned;
    totalEarned = totalEarned + gain(ntrial); % add winnings/losses of this trial to totalEarned
    feedback(ntrial) = round(1000*totalEarned)/1000;
    
    %% Display Feedback about money won/lost
    timenow1 = GetSecs;
    onset.feedback(ntrial) = timenow1;
    iBaselineGrip = 0;
    display_FeedbackRP( window, trialValence(ntrial), previoustotal, totalEarned, x, y, IRM )
    [~,timenow1,~,~,~] = Screen(window,'Flip');
    onset.feedback(ntrial) = timenow1;
    while GetSecs <= onset.feedback(ntrial) + varFdbk_Wait(ntrial)
        display_FeedbackRP( window, trialValence(ntrial), previoustotal, totalEarned, x, y, IRM )
        Screen(window,'Flip');
        iBaselineGrip = iBaselineGrip + 1;
        % extract grip values during feedback to extract baseline
        [grip,Tgrip] = ReadGripValue(GripPort,Handle,GripChannel);
        baselineGripData.feedback.grip{ntrial}(iBaselineGrip) = grip;
        baselineGripData.feedback.time{ntrial}(iBaselineGrip) = Tgrip;
    end
    % duration
    timenow2 = GetSecs;
    duration.feedback(ntrial) = timenow2 - timenow1;
    
    %% display number of the trial so that the experimenter can keep track
    % on how much there is left at the end of each trial
    disp(['Trial ',num2str(ntrial),'/',num2str(totalTrial)]);
    %% save all stuff at the end of each trial
    save([behaviordir,'\global_sub_',subid,'_session_',runname,'_',taskName,'.mat'])
end

%% fixation cross (or black screen) at the end of the script
Screen('DrawTexture',window,picCross,[],rectCross);
% onset
[~,timenow1,~,~,~] = Screen(window,'Flip');
onset.cross = [onset.cross; timenow1];
WaitSecs(ending_T);
% duration
timenow2 = GetSecs;
dur = timenow2 - timenow1;
duration.cross = [duration.cross; dur];

%% announce total gain after this run here
%% here or at the beginning you should load previous totalGain if this is not the first run
if nrun > 0 % display gain only for actual fMRI runs, useless for training(or not?)
    if nrun == 1 % for the first run, define the totalgain variable
        totalGain = nan(1,nbTotalRuns);
    elseif nrun > 1 % otherwise load totalGain and add gain of this run to it
        load([behaviordir,'\totalGain_sub_',subid,'.mat'],'totalGain')
    end
    totalGain(nrun) = sum(gain);
    save([behaviordir,'\totalGain_sub_',subid,'.mat'],'totalGain')
    timenow1 = totalGainDisp(window, totalGain(nrun));
    onset.totalGain = timenow1;
    WaitSecs(totalGain_Wait);
    % duration
    timenow2 = GetSecs;
    duration.totalGain = timenow2 - timenow1;
end

%% save data
trialtime = onset.displayIncentive';
efforttime = onset.displayEffortScale';
feedbacktime = onset.feedback';

data = [1:totalTrial; task_inc.*trialValence; peakForce; sumforce; perf; gain; feedback; trialtime; efforttime; feedbacktime]';

cd(behaviordir);
save(resultname,...
    'subid','nrun','Fmax','data','baselineGripData','gripData','peakForce',...
    'trialValence','nBlocks','nTrialsPerBlock');

%% save all TTL in a results file
% extract all the TTLs recorded inside KbQueue
if IRM == 1 && nrun > 0
    TTL = [];
    while KbEventAvail
        [event, n] = KbEventGet;
        TTL = [TTL; event.Time];
    end
    KbQueueStop;
    KbQueueRelease;
    save([behaviordir,'\TTL_sub',subid,'_sess',runname,'.mat'],'T0','TTL');
    onset.TTL = TTL;
    save([behaviordir,'\MBB_battery_',taskName,'_onsets_sub',subid,'_sess',runname,'.mat'],'onset','duration');
end

%% stop and record eye-tracking
if addEyeTrack == 1
    stopEyeTracking(subdir,eyeFileName,taskName,subid,runname);
end

%% stop grip and record all values
CloseGripDevice(GripPort,Handle);
gripFullData.grip = []; % variable where all grip data pooled to have the continuous values
gripFullData.time = []; % variable where all grip data timings pooled to have all corresponding timings
for ntrial = 1:totalTrial
    gripFullData.grip = [gripFullData.grip,...
        baselineGripData.fixationCross.grip{ntrial},...
        baselineGripData.incentive.grip{ntrial},...
        gripData.grip{ntrial},...
        baselineGripData.feedback.grip{ntrial}];
gripFullData.time = [gripFullData.time,...
    baselineGripData.fixationCross.time{ntrial},...
    baselineGripData.incentive.time{ntrial},...
    gripData.time{ntrial},...
    baselineGripData.feedback.time{ntrial}];
end
%% save all stuff
save([behaviordir,'\global_sub_',subid,'_session_',runname,'_',taskName,'.mat']) % save all variables here
cd(root); % go back to scripts folder
Screen('CloseAll'); % close PTB

