function [eodata, struct, im] = BrachyFUN
% Usage: [eodata, struct, im] = FUN
% eodata is the Spike2 recordings of the EOD and the timing pulses
% struct is the tracking data from the video images
% im are the video frames
% swimmy will ask the necessary questions... you need at least one Matlab
% file from Spike2 and one (or two) video files.  Good luck.


%% User picks the video and Spike2 files for a trial.
fprintf('Pick the video file');
    [vidfile, vidpath] = uigetfile('*.avi', 'Pick the video file');
    f = fullfile(vidpath, vidfile);
    
    vid = VideoReader(f);

fprintf('Pick the Spike2 Matlab file');    
    [spkfile, spkpath] = uigetfile('*.mat', 'Pick the Spike2 Matlab file');
    f = fullfile(spkpath, spkfile);
    load(f);

% Show a random frame for the video 
    rr = randi([1 vid.NumberOfFrames], 1);
    figure(1); clf; imshow(rgb2gray(read(vid, rr)));

    
numFish = input('How many Brachyhypopomus in the video (1 or 2): ');

close(1); % Close the video frame

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
        
    [struct, im] = bestofriends(vid, numFish, 100);

    pause(3);

%% At the end, close all of the figures and plot the data
    
%   close all;
   
%    figure(1); clf; imshow(im.A(1).gray);
%         hold on; plot(struct.A.x, struct.A.y, 'm-*');
%    if oneortwo == 2;
%         hold on; plot(struct.B.x, struct.B.y, 'c-*');
%    end

   pause(2);
   
  fprintf('Remember to save the data to prevent heartache and shame.\n');
   