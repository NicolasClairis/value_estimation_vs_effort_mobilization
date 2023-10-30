function [ picCross, rectCross, pic_inc, rect_inc] = load_incentive_images(window,x,y, root, nConditions)
%load_incentive_images loads images used for incentive in taskMentalRP and
%in taskGripRP
%
% INPUTS
% window: window where will be displayed
% root: folder where scripts are
% nConditions: number of stimuli
%
% developed by Nicolas Clairis - february 2017 for fMRI task

%% Load fixation cross image
picCross = Screen('MakeTexture',window,imread('Cross.bmp'));
rectCross = CenterRectOnPoint(Screen('Rect',picCross),x,y);

%% extract coordinates of the biggest of the incentive images (by now 20euros, 6th image)
BigImage = 6; % the 6th picture = 20euros is the biggest one
pic_incBigImage = Screen('MakeTexture',window,imread([root, filesep, 'gripim', filesep, 'pic_inc_' num2str(BigImage) '.bmp']));
[wrectBigIncentiveImage,hrectBigIncentiveImage] = RectSize(Screen('Rect',pic_incBigImage));

%% load incentive images
% coordinates for the center of the images
xImageCenter = x/3;
yImageCenter = y/3;
for iIncentiveImage = 1:nConditions
    pic_inc{iIncentiveImage,1} = Screen('MakeTexture',window,imread([root, filesep, 'gripim', filesep, 'pic_inc_' num2str(iIncentiveImage) '.bmp']));
    [wrect{iIncentiveImage},hrect{iIncentiveImage}] = RectSize(Screen('Rect',pic_inc{iIncentiveImage}));
    % re-scaling factor to adapt image size to any screen type (adapt based
    % on the biggest of the images to keep the difference of size between
    % the images)
    reScalingFactor = 3*x/(6*wrectBigIncentiveImage);
    % extract rectangle coordinates where image should be displayed
    % old and ugly way to do it:
    % rect_inc{iIncentiveImage} = CenterRectOnPoint(Screen('Rect',pic_inc{iIncentiveImage}),x-wrect{iIncentiveImage}/2-300,y+hrect{iIncentiveImage}/2-300);
    % actual way:
    rect_inc{iIncentiveImage} = CenterRectOnPoint(Screen('Rect',pic_inc{iIncentiveImage})*reScalingFactor, xImageCenter, yImageCenter);
    % also draw loss incentive
    pic_inc{iIncentiveImage,2} = Screen('MakeTexture',window,imread([root, filesep, 'gripim', filesep, 'pic_inc_' num2str(iIncentiveImage) 'neg.bmp']));
end

end