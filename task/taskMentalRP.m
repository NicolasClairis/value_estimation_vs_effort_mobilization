% function [] = taskMentalRP(varargin)
%taskMentalRP - execute the motivational numeric stroop task (from motiscan battery)
% Launch the task with this function for testing one subject.
% The subject has to compare sequences of number pairs, each pair contains
% numbers with different values and font size, the rule is to select the
% number with the highest numerical value. Each trial is rettributed
% depending on the performance of the subject (% of correct comparisons).
%
%   task specifications: 
%       - structure: instructions -> calibration -> instructions -> training -> testing
%       - experimental conditions: between-trials incentive(0.01,0.20,0.50,1,5,20 euros),
%         valence(gain/loss), within-trials feature congruency(1/0),numerical values (from 0 to 10)
%       - randomization: within-block random incentive, switching valence blocks , 
%         within-trial balanced random permutation of congruency(50%),
%         within-trial balanced random permutation of numerical distances(1,2,3,4,5)
%       - ntrial = 120
%       - trial structure: 10 number pairs to compare within a time-limit,
%                          time-penalty for incorrect responses
%       - calibration: estimate best time performance for 10 correct
%         comparisons (5 trials) -> adapt subject time-limit to 70% of best time performance
%       - training: exposure to all between-trials condition levels (12 trials)
%
%
% Syntax:  taskMentalRP(subid,nsession,'optionName',optionValue)
%
% Inputs:
%    subid - subject identification number (double)
%    nsession - session number (double)
%    options:
%       calibtime - predefined subject's time-limit, no default value (double)
% Outputs:
%
% Example: 
%   taskMentalRP(10,2,'calibtime',5.2)
%
% Requirements: 
%   Subfunctions:   display_FeedbackRP.m , display_mentalEffort.m ,
%                   identification_batmotiv.m , setDir.m , 
%   MAT-files: study.mat
%   MATLAB products: MATLAB, Statistics and Machine Learning Toolbox,
%                    Psychtoolbox
%
% See also:
%
% Author: Raphael Le Bouc (first version), Nicolas Borderies
% (arrangements), Nicolas Clairis (re-adapted for CENIR fMRI)
% email address: nico.borderies@gmail.com 
% April 2014; Last revision: February 2017 by Nicolas Borderies and Nicolas
% Clairis

clear all; close all; clc;

testing_script = 0; % parameters of screen to check on my own computer (1) or (0) when at CENIR

%% add eye-tracking functions folder
addpath(genpath('eye_tracking_functions'));

%% Identification
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
taskName = 'mentalRP';

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
% open and calibrate eye-tracker
if IRM == 1 % no eye-tracking outside fMRI
    addEyeTrack = input('Eye-tracking? 1 = yes / 0 = no');
    if addEyeTrack == 1
        [eyeFileName] = startEyeTracking(subid,runname);
    end
else
    addEyeTrack = 0;
end

%% Screen configuration
[x, y, window, baselineTextSize] = ScreenConfiguration(IRM, testing_script);
yGoPosition = y/3;

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

taskPerf = nan(1,totalTrial);
gain = nan(1,totalTrial);
feedback = nan(1,totalTrial);

% variables mental_stroop-specific
Physical_Or_Mental.Phys_Or_Ment = 2; % mental stroop and not physical effort
nNumber = 10; % Number of pairs displayed on screen
distance = sort(repmat(1:5,1, nNumber/5)); % numerical distance within number pairs
fontMin = 24;
fontMax = 40;
keepMax = 0; % interim calibration value
% calibration variables
totalAnswers_Calibration = zeros(1,nCalibrationTrials); % variable to determine what is the maximum of answers possible to do in PressTime_Wait seconds
goodAnswers_Calibration = zeros(1,nCalibrationTrials); % variable to determine what is the maximum of correct answers possible to do in PressTime_Wait seconds
% task variables
totalAnswers = zeros(1,totalTrial); % variable to determine what is the maximum of answers possible to do in PressTime_Wait seconds
goodAnswers = zeros(1,totalTrial); % variable to determine what is the maximum of correct answers possible to do in PressTime_Wait seconds

%% Loads images and creates positions
[ picCross, rectCross, pic_inc, rect_inc] = load_incentive_images(window,x,y, root, nConditions);

%% Generator reset
rand('state',sum(100*clock));

%% Time variables
% calibration
calibStartMessage_T = 1; % "C'est parti" message duration
instruction_T = 3;
instructionCalibration_T = 2;
feedbackCalibration_T = 1;
penaltyCalibration_T = 0.5; % calibration time penalty
penaltyFactor = 0.1; % scaling factor of time penalty with respect to time-limit
limitFactor = 0.7; % scaling factor of time-limit with respect to calibration time
% task
cross_ITI_T = 0.5;
% final cross
if nrun == 0 % shorter for training
    ending_T = 5;
else
    ending_T = 10;
end
release_Wait = 0.5; % duration of cross after "please release"
PressTime_Wait = 3; % max time for waiting a response
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
% timings specific to calibration
listA = {'instruction_Calibration','displayEffortScale_Calibration','feedback_Calibration'}; % fixation cross before announcing the switch of the block
for l = 1 : length(listA)
    onset.(listA{l}) = NaN(nCalibrationTrials,1);
    duration.(listA{l}) = NaN(nCalibrationTrials,1);
end

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
onset.errorFeedback_Calibration = []; % onset error made during calibration
duration.errorFeedback_Calibration = []; % onset error made during calibration
onset.errorFeedback = []; % onset error made during the task
duration.errorFeedback = []; % onset error made during the task
onset.instruction = [];
duration.instruction = [];
onset.cross = []; % fixation cross when no ITI (start of task, end of task)
duration.cross = [];
onset.totalGain = []; % display of total gain at the end of the run
duration.totalGain = [];

%% keyboard keys configuration + waiting and recording first TTL for fMRI
if IRM == 0
    %% key configuration
    KbName('UnifyKeyNames');
    key.left = KbName('LeftArrow');
    key.right = KbName('RightArrow');
    key.escape= KbName('escape');
elseif IRM == 1
    %% fMRI key configuration
    KbName('UnifyKeyNames');
    key.left = 49;        % 49 % DROITE  bleu, left press
    key.right = 50;       %50   %% GAUCHE JAUNE
    key.escape = KbName('escape');
end

%-----------------------------------------------
%% Calibration
%-----------------------------------------------

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
    Screen('TextSize', window, 40);
    DrawFormattedText(window,'D''abord quelques essais pour s''echauffer','center','center',[255 255 255],wrapat,[],[],vspacing)
    textPleasePress = 'appuyer sur une touche pour continuer';
    Screen('TextSize', window, 20);
    [width,hight] = RectSize(Screen('TextBounds',window,textPleasePress));
    Screen('DrawText',window, textPleasePress, x-width/2, 9*y/5, [100 100 100]);
    % onset
    [~,timenow1,~,~,~] = Screen(window,'Flip');
    onset.instruction = [onset.instruction; timenow1];
    KbWait;
    timenow2 = GetSecs;
    duration.instruction = [duration.instruction; timenow2 - timenow1];
    Screen('TextSize', window, baselineTextSize);
end

% signal that it starts
Screen('TextSize', window, 60);
[width,hight] = RectSize(Screen('TextBounds',window,'C''est parti!'));
Screen('DrawText',window,'C''est parti!',x-width/2,y-hight/2-100,[255 0 0]);
Screen('TextSize', window, 20);
[~,timenow1,~,~,~] = Screen(window,'Flip');
onset.instruction = [onset.instruction; timenow1];
WaitSecs(calibStartMessage_T);
timenow2 = GetSecs;
duration.instruction = [duration.instruction; timenow2 - timenow1];

% calibration
for calibTrial = 1:nCalibrationTrials
    
    %% Display instruction
    Screen('TextSize', window, baselineTextSize);
    [width,hight] = RectSize(Screen('TextBounds',window,[num2str(calibTrial) 'e essai']));
    Screen('DrawText',window,[num2str(calibTrial) 'e essai'],x-width/2,y-hight/2-100,[255 0 0]);
    Screen('TextSize', window, 20);
    [width,hight] = RectSize(Screen('TextBounds',window,'Appuyez le plus vite possible du cote du nombre'));
    Screen('DrawText',window,'Appuyez le plus vite possible du cote du nombre',x-width/2,y-hight/2,[255 0 0]);
    [width,hight] = RectSize(Screen('TextBounds',window,'le plus grand numeriquement'));
    Screen('DrawText',window,'le plus grand numeriquement',x-width/2,y-hight/2+50,[255 0 0]);
    Screen('TextSize', window, baselineTextSize);
    [~,timenow1,~,~,~] = Screen(window,'Flip');
    onset.instruction_Calibration = timenow1;
    WaitSecs(instructionCalibration_T);
    % duration starting calibration
    timenow2 = GetSecs;
    duration.instruction_Calibration =  timenow2 - timenow1;
    
    %% Display mental stroop
    % Trial structure
    % generate within-trial conditions
    %%% numbers
    order = randperm(nNumber);
    distance = distance(order);
    rep = ceil(nNumber/10);
    aga = [];
    for iLoop = 1:rep
        aga = [aga (0:4) + distance(1:5) (5:9) - distance(6:10)];
    end
    prelist = [repmat(0:9,1,rep); aga];
    prelist = prelist(:,1:nNumber);
    
    %%% randomization
    order = randperm(nNumber); % randomize order
    prelist = prelist(:,order);
    side = randperm(nNumber); % randomize side
    side(side<(nNumber/2+1)) = 1;
    side(side>(nNumber/2)) = 2;
    list = prelist;
    list(1,side == 1) = max(prelist(:,side == 1));
    list(2,side == 1) = min(prelist(:,side == 1));
    list(1,side == 2) = min(prelist(:,side == 2));
    list(2,side == 2) = max(prelist(:,side == 2));
    
    %%% font size
    fontSize = zeros(2,nNumber) + fontMin;
    fontSize(1,side == 1) = fontMax;
    fontSize(2,side == 2) = fontMax;
    font = randperm(nNumber); % randomize congruency
    font = (font>(nNumber/2));
    for iPair = 1:nNumber
        if font(iPair) == 1
            temp = fontSize(1,iPair);
            fontSize(1,iPair) = fontSize(2,iPair);
            fontSize(2,iPair) = temp;
        end
    end
    %%% conditions storage
    numbers{calibTrial}.list = list;
    numbers{calibTrial}.fontSize = fontSize;
    
% start display
    calibrationtime = GetSecs;
    goodAnswers_Calibration(calibTrial) = 0;
    totalAnswers_Calibration(calibTrial) = 0;
    choice = 0;
    isWrongAnswer = 0;
    Physical_Or_Mental.Mental.list = list;
    Physical_Or_Mental.Mental.fontSize = fontSize;
    onset.displayEffortScale_Calibration(calibTrial) = GetSecs;
    % Monitor responses
    while goodAnswers_Calibration(calibTrial) < nNumber
        
        % display conditions & performance
        currentPair = goodAnswers_Calibration(calibTrial) + 1;
        nn = list(:, currentPair);
        fs = fontSize(:, currentPair);
        level = (goodAnswers_Calibration(calibTrial)/nNumber)*100;
        Screen('TextSize', window, 60);
        [width,hight] = RectSize(Screen('TextBounds',window,'Go!'));
        Screen('DrawText',window,'Go!',x-width/2,yGoPosition-hight/2,[255 0 0]);
        %         display_mentalEffort(window,x,y,level, list, fontSize);
        Levels.actualLevel = level;
        Levels.actualFmax = level;
        display_Effort(window, x, y, Levels, 1, Physical_Or_Mental);
        
        % error monitoring
        if isWrongAnswer == 1 % if error made => keep on same pair and display it in red
            Physical_Or_Mental.Mental.error = currentPair; % extract pair where error made
            display_Effort(window, x, y, Levels, 1, Physical_Or_Mental);
%             yaxis = [-300:60:300];
%             for iPairItem = 1:2
%                 Screen('TextSize', window, fs(iPairItem));
%                 [w,h] = RectSize(Screen('TextBounds',window,num2str(nn(iPairItem))));
%                 Screen('DrawText', window, num2str(nn(iPairItem)), x+150-w/2 + ((iPairItem-1)*2-1)*25, y-yaxis(goodAnswers_Calibration(calibTrial)+2)-h/2, [255 0 0]);
%             end
            [~,timenow1,~,~,~] = Screen(window,'Flip');
            onset.errorFeedback_Calibration = [onset.errorFeedback_Calibration; timenow1];
            WaitSecs(penaltyCalibration_T);
            timenow2 = GetSecs;
            dur = timenow2 - timenow1;
            duration.errorFeedback_Calibration = [duration.errorFeedback_Calibration; dur];
            Physical_Or_Mental.Mental = rmfield(Physical_Or_Mental.Mental,'error'); % delete it from the structure after using it for the display
            isWrongAnswer = 0;
        else
            Screen(window,'Flip');
        end
        
        % performance monitoring
        [keyisdown, secs, keycode] = KbCheck;
        if keyisdown == 1 && keycode(key.right) == 1
            timedown = GetSecs;
            choice = 2;
        elseif keyisdown == 1 && keycode(key.left) == 1
            timedown = GetSecs;
            choice = 1;
        end
        KbReleaseWait();
        if choice == side(goodAnswers_Calibration(calibTrial) + 1)
            totalAnswers_Calibration(calibTrial)    = totalAnswers_Calibration(calibTrial) + 1;
            goodAnswers_Calibration(calibTrial)     = goodAnswers_Calibration(calibTrial)  + 1;
            choice = 0;
            data_Calibration.answer{calibTrial}(totalAnswers_Calibration(calibTrial))       = 1;
            data_Calibration.answertime{calibTrial}(totalAnswers_Calibration(calibTrial))	= secs;
        elseif choice == 3 - side(goodAnswers_Calibration(calibTrial) + 1)
            totalAnswers_Calibration(calibTrial)	= totalAnswers_Calibration(calibTrial) + 1;
            choice = 0;
            data_Calibration.answer{calibTrial}(totalAnswers_Calibration(calibTrial))       = 0;
            data_Calibration.answertime{calibTrial}(totalAnswers_Calibration(calibTrial))	= secs;
            isWrongAnswer = 1; % answer is wrong this time
        end
    end
    duration.displayEffortScale_Calibration(calibTrial) = GetSecs - onset.displayEffortScale_Calibration(calibTrial);
    
    % update performance
    keepMax = max(keepMax, goodAnswers_Calibration(calibTrial));
    % update calibration
    keepmin = nanmin(duration.displayEffortScale_Calibration);
    
    %% Show feedback
    onset.feedback_Calibration(calibTrial) = GetSecs;
    dur = 0;
    while dur <= feedbackCalibration_T
        %         display_mentalEffort(window,x,y,level, list, fontSize);
        % display full bar at the end of the trial
        Levels.actualLevel = 100;
        Levels.actualFmax = 100;
        display_Effort(window, x, y, Levels, 1, Physical_Or_Mental)
        Screen('TextSize', window, 60);
        [width,hight] = RectSize(Screen('TextBounds',window,'Go!'));
        Screen('DrawText',window,'Go!',x-width/2,yGoPosition-hight/2,[255 0 0]);
        Screen('TextSize', window, 40);
        [width,hight] = RectSize(Screen('TextBounds',window,'Record'));
        Screen('DrawText',window,'Record',2*x/5-width/2,y-hight/2-100,[175 175 175]);
        Screen('TextSize', window, 60);
        [width,hight] = RectSize(Screen('TextBounds', window, [num2str(keepmin) 's']));
        Screen('DrawText',window,[num2str(keepmin) 's'], 2*x/5-width/2, y-hight/2, [175 175 175]);
        [~,timenow2,~,~,~] = Screen(window,'Flip');
        dur = timenow2 - onset.feedback_Calibration(calibTrial);
    end
    % duration starting cross
    duration.feedback_Calibration = dur;
end


%-----------------------------------------------
%% Time-limit adaptation
%-----------------------------------------------
calibTime = limitFactor*keepmin; % keep faster performance as time for 
timePenalty = penaltyFactor*calibTime;

%-----------------------------------------------
%% actual task
%-----------------------------------------------
%% waiting and recording first TTL for fMRI
if IRM == 1
    % TTL trigger
    dummy_scan = 4;
    trigger = 53;
    
    % add a cross at the beginning (better for eye tracking)
    Screen('DrawTexture',window,picCross,[],rectCross);
    % onset
    [~,timenow1,~,~,~] = Screen(window,'Flip');
    onset.cross = [onset.cross; timenow1];
    
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
        keysOfInterest(trigger) = 1; % check TTL
        keysOfInterest(key.left) = 1; % check all left key press
        keysOfInterest(key.right) = 1; % check all right key press
        KbQueueCreate(0,keysOfInterest); % checks TTL and keys of pad
        KbQueueStart; % starts checking
    end
    
    % duration starting cross
    timenow2 = GetSecs;
    dur = timenow2 - timenow1;
    duration.cross = [duration.cross; dur];
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
    onset.instruction = [onset.instruction; timenow1];
    KbWait;
    timenow2 = GetSecs;
    duration.instruction = [duration.instruction; timenow2 - timenow1];
    Screen('TextSize', window, baselineTextSize);
end

%% TASK
block = 0; % block index
for ntrial = 1:totalTrial
    
% Trial structure
    %% generate within-trial conditions
    %% numbers
    order = randperm(nNumber);
    distance = distance(order);
    rep = ceil(nNumber/10);
    aga = [];
    for iLoop = 1:rep
        aga = [aga (0:4) + distance(1:5) (5:9) - distance(6:10)];
    end
    prelist = [repmat(0:9,1,rep); aga];
    prelist = prelist(:,1:nNumber);
    
    %% randomization
    order = randperm(nNumber); % randomize order
    prelist = prelist(:,order);
    side = randperm(nNumber); % randomize side
    side(side<(nNumber/2+1)) = 1;
    side(side>(nNumber/2)) = 2;
    list = prelist;
    list(1,side == 1) = max(prelist(:,side == 1));
    list(2,side == 1) = min(prelist(:,side == 1));
    list(1,side == 2) = min(prelist(:,side == 2));
    list(2,side == 2) = max(prelist(:,side == 2));
    
    %% font size
    fontSize = zeros(2,nNumber) + fontMin;
    fontSize(1,side == 1) = fontMax;
    fontSize(2,side == 2) = fontMax;
    font = randperm(nNumber); % randomize congruency
    font = (font>(nNumber/2));
    for iPair = 1:nNumber
        if font(iPair) == 1
            temp = fontSize(1,iPair);
            fontSize(1,iPair) = fontSize(2,iPair);
            fontSize(2,iPair) = temp;
        end
    end
    %% conditions storage
    taskNumbers{ntrial}.list = list;
    taskNumbers{ntrial}.fontSize = fontSize;
    
    %% fixation cross (while to be homogeneous with grip but useless, could be replaced by a WaitSecs)
    Screen('DrawTexture',window,picCross,[],rectCross);
    [~,timenow1,~,~,~] = Screen(window,'Flip');
    onset.cross_ITI(ntrial) = timenow1;
    dur = 0;
    while dur <= cross_ITI_T
        % keep displaying cross
        Screen('DrawTexture',window,picCross,[],rectCross);
        [~,timenow2,~,~,~] = Screen(window,'Flip');
        % check duration
        dur = timenow2 - timenow1;
    end
%     WaitSecs(cross_ITI_T);
    % duration
%     timenow2 = GetSecs;
    duration.cross_ITI(ntrial) = dur;
    
    %% display incentive (while to be homogeneous with grip but useless, could be replaced by a WaitSecs)
    Screen('DrawTexture',window,pic_inc{task_inc(ntrial), 3-(trialValence(ntrial)+3)/2},[],rect_inc{task_inc(ntrial)});
    Screen('TextSize', window, 32);
    incentiveLeft = rect_inc{task_inc(ntrial)}(1);
    incentiveBottom = rect_inc{task_inc(ntrial)}(4);
    if trialValence(ntrial) == 1
        Screen('DrawText', window, 'A GAGNER', incentiveLeft, incentiveBottom, [255 255 255], [], 0);
    elseif trialValence(ntrial) == -1
        Screen('DrawText', window, 'A PERDRE', incentiveLeft, incentiveBottom, [255 255 255], [], 0);
    end
    [~,timenow1,~,~,~] = Screen(window,'Flip');
    onset.displayIncentive(ntrial) = timenow1;
    dur = 0;
    while dur <= varInc_Wait(ntrial)
        % keep displaying incentive
        Screen('DrawTexture',window,pic_inc{task_inc(ntrial), 3-(trialValence(ntrial)+3)/2},[],rect_inc{task_inc(ntrial)});
        Screen('TextSize', window, 32);
        if trialValence(ntrial) == 1
            Screen('DrawText', window, 'A GAGNER', incentiveLeft, incentiveBottom, [255 255 255], [], 0);
        elseif trialValence(ntrial) == -1
            Screen('DrawText', window, 'A PERDRE', incentiveLeft, incentiveBottom, [255 255 255], [], 0);
        end
        [~,timenow2,~,~,~] = Screen(window,'Flip');
        % check duration
        dur = timenow2 - timenow1;
    end
    duration.displayIncentive(ntrial) = dur;

    %% Display Mental Effort Scale
    % start display
    choice = 0;
    goodAnswers(ntrial) = 0;
    totalAnswers(ntrial) = 0;
    isWrongAnswer = 0;
    onset.displayEffortScale(ntrial) = GetSecs;

    % Monitor responses
    while GetSecs <= (onset.displayEffortScale(ntrial) + calibTime)
        if goodAnswers(ntrial) < nNumber
            currentPair = goodAnswers(ntrial) + 1;
            % display conditions & performance
            nn = list(:, currentPair);
            fs = fontSize(:, currentPair);
            level = goodAnswers(ntrial)/nNumber*100;
            Screen('DrawTexture',window,pic_inc{task_inc(ntrial), 3-(trialValence(ntrial)+3)/2},[],rect_inc{task_inc(ntrial)});
            Screen('TextSize', window, 32);
            if trialValence(ntrial) == 1
                Screen('DrawText', window, 'A GAGNER', incentiveLeft, incentiveBottom, [255 255 255], [], 0);
            elseif trialValence(ntrial) == -1
                Screen('DrawText', window, 'A PERDRE', incentiveLeft, incentiveBottom, [255 255 255], [], 0);
            end
            %             display_mentalEffort(window, x, y, level, list, fontSize, outcomeValues(task_inc(ntrial)));
            % load display_Effort parameters
            Levels.actualLevel = level;
            Levels.incentive = outcomeValues(task_inc(ntrial));
            Levels.actualFmax = level;
            Physical_Or_Mental.Mental.list = list;
            Physical_Or_Mental.Mental.fontSize = fontSize;
            display_Effort(window, x, y, Levels, trialValence(ntrial), Physical_Or_Mental);

            % error monitoring
            if isWrongAnswer == 1 % if error made => keep on same pair and display it in red
                Physical_Or_Mental.Mental.error = currentPair; % extract pair where error made
                display_Effort(window, x, y, Levels, trialValence(ntrial), Physical_Or_Mental);
                %             yaxis = [-300:60:300];
                %             for iPairItem = 1:2
                %                 Screen('TextSize', window, fs(iPairItem));
                %                 [w,h] = RectSize(Screen('TextBounds',window,num2str(nn(iPairItem))));
                %                 Screen('DrawText', window, num2str(nn(iPairItem)), x+150-w/2 + ((iPairItem-1)*2-1)*25, y-yaxis(goodAnswers_Calibration(calibTrial)+2)-h/2, [255 0 0]);
                %             end
                [~,timenow1,~,~,~] = Screen(window,'Flip');
                onset.errorFeedback = [onset.errorFeedback; timenow1];
                WaitSecs(timePenalty);
                timenow2 = GetSecs;
                dur = timenow2 - timenow1;
                duration.errorFeedback = [duration.errorFeedback; dur];
                Physical_Or_Mental.Mental = rmfield(Physical_Or_Mental.Mental,'error'); % delete it from the structure after using it for the display
                isWrongAnswer = 0;
            else
                Screen(window,'Flip');
            end
            
            % performance monitoring
            [keyisdown, secs, keycode] = KbCheck;
            if keyisdown == 1 && keycode(key.right) == 1
                timedown = GetSecs;
                choice = 2;
            elseif keyisdown == 1 && keycode(key.left) == 1
                timedown = GetSecs;
                choice = 1;
            end
            KbReleaseWait();
            if choice == side(goodAnswers(ntrial) + 1)
                totalAnswers(ntrial) = totalAnswers(ntrial) + 1;
                goodAnswers(ntrial) = goodAnswers(ntrial) + 1;
                choice = 0;
                data_Task.answer{ntrial}(totalAnswers(ntrial)) = 1;
                data_Task.answertime{ntrial}(totalAnswers(ntrial)) = secs;
            elseif choice == 3 - side(goodAnswers(ntrial) + 1)
                totalAnswers(ntrial) = totalAnswers(ntrial) + 1;
                choice = 0;
                data_Task.answer{ntrial}(totalAnswers(ntrial)) = 0;
                data_Task.answertime{ntrial}(totalAnswers(ntrial)) = secs;
                isWrongAnswer = 1;
            end
        else
            % display conditions & performance
            level = goodAnswers(ntrial)/nNumber*100;
            Screen('TextSize', window, 60);
            [width,hight] = RectSize(Screen('TextBounds',window,'Go!'));
            Screen('DrawText',window,'Go!',x-width/2,yGoPosition-hight/2,[255 0 0]);
%             display_mentalEffort(window,x,y,level, list, fontSize);
            display_Effort(window, x, y, Levels, trialValence(ntrial), Physical_Or_Mental);
            Screen(window,'Flip');
        end
    end
    duration.displayEffortScale(ntrial) = GetSecs - onset.displayEffortScale(ntrial);
    
    % Data to save
    taskPerf(ntrial) = goodAnswers(ntrial)/nNumber*100;
    
    % update monetary outcome
    switch trialValence(ntrial)
        case 1
            gain(ntrial) = taskPerf(ntrial)/100*outcomeValues(task_inc(ntrial));
        case -1
            gain(ntrial) = (100-taskPerf(ntrial))/100*-outcomeValues(task_inc(ntrial));
            gain(ntrial) = gain(ntrial)*(gain(ntrial)<0);                    % in the loss condition, the best outcom is zero.
    end
    previoustotal = totalEarned;
    totalEarned = totalEarned + gain(ntrial);
    feedback(ntrial) = round(1000*totalEarned)/1000;
    
    %% Display Feedback
    display_FeedbackRP( window, trialValence(ntrial), previoustotal, totalEarned, x, y, IRM )
    [~,timenow1,~,~,~] = Screen(window,'Flip');
    onset.feedback(ntrial) = timenow1;
    WaitSecs(varFdbk_Wait(ntrial));
    % duration
    timenow2 = GetSecs;
    duration.feedback(ntrial) = timenow2 - timenow1;
    
    % display number of the trial so that the experimenter can keep track
    % on how much there is left at the end of each trial
    disp(['Trial ',num2str(ntrial),'/',num2str(totalTrial)]);
    % save all stuff at the end of each trial in the subject's results directory
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
        totalGain = nan(1, nbTotalRuns);
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

%% Data saving
%-----------------------------------------------

trialtime = onset.displayIncentive';
efforttime = onset.displayEffortScale';
feedbacktime = onset.feedback';

% data aggregation
data = [1:totalTrial; task_inc.*trialValence; goodAnswers; totalAnswers; taskPerf; gain; feedback; trialtime; efforttime; feedbacktime]';
% calibration variables stored in one structure
calibration.trials = 1:nCalibrationTrials;
calibration.goodAnswers = goodAnswers_Calibration; % number of good answers per calibration trial (useless in that frame, but useful if time is limited)
calibration.totalAnswers = totalAnswers_Calibration; % number of total answers per calibration trial (useless in that frame, but useful if time is limited)
calibration.answer = data_Calibration.answer; % extract all subject's answers during calibration
calibration.answertime = data_Calibration.answertime; % extract all subject's answers timing during calibration
calibration.calibTime = calibTime; % time calibrated for nNumber pairs

% data directory
cd(subdir);
save([behaviordir, filesep, resultname],'data','subid','nrun',...
                'calibration','data_Task',...
                'numbers','taskNumbers');

            
%% save all TTL and onsets in a results file
if IRM == 1 && nrun > 0
    KbQueueStop; % stop recording key presses
    keyLeft.Start = []; % time when starts pressing left key
    keyLeft.Release = []; % time when releases right key
    keyRight.Start = []; % time when starts pressing right key
    keyRight.Release = []; % time when releases right key
    while KbEventAvail
        [event, n] = KbEventGet;
        if event.Keycode == trigger
            TTL = [TTL; event.Time];
        elseif event.Keycode == key.left % if left key pressed
            if event.Pressed == 1 % record start of press
                keyLeft.Start = [keyLeft.Start; event.Time];
            elseif event.Pressed == 0 % record time when release
                keyLeft.Release = [keyLeft.Release; event.Time];
            end
        elseif event.Keycode == key.right % if right key pressed
            if event.Pressed == 1 % record start of press
                keyRight.Start = [keyRight.Start; event.Time];
            elseif event.Pressed == 0 % record time when release
                keyRight.Release = [keyRight.Release; event.Time];
            end
        end
    end
    KbQueueRelease;
    save([behaviordir,'\TTL_sub',subid,'_sess',runname,'.mat'],'T0','TTL');
    onset.TTL = TTL;
    save([behaviordir,'\MBB_battery_',taskName,'_onsets_sub',subid,'_sess',runname,'.mat'],'onset','duration');
end
           
%% stop and record eye-tracking
if addEyeTrack == 1
    stopEyeTracking(subdir, eyeFileName, taskName, subid, runname);
end

%% save all stuff
save([behaviordir,'\global_sub_',subid,'_session_',runname,'_',taskName,'.mat']) % save all variables here
cd(root); % go back to scripts folder
Screen('CloseAll'); % close PTB