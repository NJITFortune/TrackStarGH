function [eodata, struct, im] = swimmy
% Usage: [eodata, struct, im] = swimmy
% eodata is the Spike2 recordings of the EOD and the timing pulses
% struct is the tracking data from the video images
% im are the video frames
% swimmy will ask the necessary questions... you need at least one Matlab
% file from Spike2 and one (or two) video files.  Good luck.

oneortwo = input('How many videos (1 or 2 are the only acceptable answers!): ');

%if oneortwo==1;
%aorb = input('Which tank is the video filmed in A or B (a or b are the
%only acceptable answers!): ');
%if aorb=a;
%% User picks the video and Spike2 files for a trial.
fprintf('Pick the Tank A video file');
    [vidfile, vidpath] = uigetfile('*.avi', 'Pick the video file');
    f = fullfile(vidpath, vidfile);
    
    vidA = VideoReader(f);
%if aorb = b;
%fprintf('Pick the Tank B video file');
    %[vidfile, vidpath] = uigetfile('*.avi', 'Pick the video file');
   % f = fullfile(vidpath, vidfile);
    
   % vidB = VideoReader(f);
    
if oneortwo == 2;
fprintf('Pick the Tank B video file');
    [vidfile, vidpath] = uigetfile('*.avi', 'Pick the video file');
    f = fullfile(vidpath, vidfile);
    
    vidB = VideoReader(f);
end;

fprintf('Pick the Spike2 Matlab file');    
    [spkfile, spkpath] = uigetfile('*.mat', 'Pick the Spike2 Matlab file');
    f = fullfile(spkpath, spkfile);
    load(f);
    
%% Simplify the variable names from Spike2
    curvars = who; 

    % Find any variables with 'Ch' in them and simplify to the title
    for x=1:1:length(curvars)
        if (strfind(curvars{x},'Ch'))
                s = eval([curvars{x},'.title']);
                changeit = ['eodata.',s , ' = ', curvars{x}];
                eval(changeit);
         end
    end
     
%% Track the fish in the video

% Show a random frame for VidA
    rr = randi([1 vidA.NumberOfFrames], 1);
    figure(1); clf; imshow(rgb2gray(read(vidA, rr)));
    
    numFish = input('How many fish are in tank A? ');
    
    [struct.A, im.A] = traxer(vidA, numFish, 100);

    pause(3);

    if oneortwo == 2;
% Show a random frame for VidB
    rr = randi([1 vidB.NumberOfFrames], 1);
    figure(1); clf; imshow(rgb2gray(read(vidB, rr)));
    
    numFish = input('How many fish are in tank B? ');
    
    [struct.B, im.B] = traxer(vidB, numFish, 100, 9999);
    
    pause(3);
    end    
    %% At the end, close all of the figures and plot the data
    
%   close all;
   
%    figure(1); clf; imshow(im.A(1).gray);
%         hold on; plot(struct.A.x, struct.A.y, 'm-*');
%    if oneortwo == 2;
%         hold on; plot(struct.B.x, struct.B.y, 'c-*');
%    end

   pause(2);
   
  fprintf('Remember to save the data to prevent heartache and shame.\n');
   