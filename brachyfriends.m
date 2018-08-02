function [fish, imdata] = brachyfriends(vdata, fishNum, dil, rango)
% Usage: [fish imdata] = traxer(vid, fishNum, dilation-num, rango);
% This is a simple electric fish tracker based on the Matlab object-tracking script example. Nothing fancy here.
% vid is the output from VideoReader.
% fishNum is the number of fish.  Currently "supports" up to 4
% dilation-num is the number of frames for computing the background - 50 is a good number
% rango (optional) is the frame range and can be used with the threshold 
%   level e.g. [startframe endframe level] 
%
% fish is a structure of output data, imdata is the images structure(grayscale and thresholded)
%

% figure(1); clf; figure(3); clf;

%% Preparations

finalframe = vdata.NumberOfFrames;
dst = [];  
jumpthresh = 15; % Number of pixels a fish may move between frames before we think something went wrong

if fishNum > 6; fprintf('The lazy troll who wrote this code has limited you to 6 fish. \n'); end

boxy = [10 10];

% A list of colors for plotting. This is why we have a limit of 6 fish at 
% the moment. Obviously we can add additional colors if necessary
colr(1,:)='r*'; colr(2,:)='b*'; colr(3,:)='m*'; colr(4,:)='g*'; colr(5,:)='c*'; colr(6,:)='k*';
colrl(1,:)='r-'; colrl(2,:)='b-'; colrl(3,:)='m-'; colrl(4,:)='g-'; colrl(5,:)='c-'; colrl(6,:)='k-';
colo(1,:)='ro'; colo(2,:)='bo'; colo(3,:)='mo'; colo(4,:)='go'; colo(5,:)='co'; colo(6,:)='ko';

fprintf('Number of frames in video: %i \n', finalframe);

% If the user gave a frame range, then we use that
if nargin > 3
    if rango(1) == 9999
            flipper = 1;
            startframe = 1;
            endframe = vdata.NumberOfFrames;
    else
            flipper = 0;
        startframe = rango(1);
        endframe = rango(2);
        if length(rango) == 3
            userpick = rango(3); % This is the threshold level
        end
    end
else
% Otherwise we as or do the whole video.
    startframe = 1;
    endframe = vdata.NumberOfFrames;
    rango = []; % input('Enter range (e.g. [1000 1500]) or hit return to do entire video: ');
    flipper = 0;
end

% Convert to grayscale and reduce to 1/2 size
% grae is the processed grayscale video - this takes time
fprintf('Importing video and reducing size.\n');
if flipper == 0
    for i = endframe:-1:startframe
        tmp = rgb2gray(read(vdata, i));
        grae(:,:,1+(i-startframe)) = tmp(1:2:end,1:2:end);
    end
end
if flipper == 1 % ROTATE 180 degrees
    for i = endframe:-1:startframe
        tmp = rgb2gray(read(vdata, i));
        grae(:,:,1+(i-startframe)) = tmp(end:-2:1,end:-2:1);
    end
end
    clear tmp;

% Initialize variables for first fish
fprintf('Allocating memory.\n');

    fish(fishNum).x = zeros(endframe-startframe,1); % centroid X
    fish(fishNum).y = zeros(endframe-startframe,1); % centroid Y
    fish(fishNum).orient = zeros(endframe-startframe,1); % Orientation of major axis
    fish(fishNum).majorLength = zeros(endframe-startframe,1); % Length of major axis
    fish(fishNum).minorLength = zeros(endframe-startframe,1); % Length of minor axis
    fish(fishNum).majorXs = zeros(endframe-startframe,2); % X pairs for major axis line
    fish(fishNum).majorYs = zeros(endframe-startframe,2); % Y pairs for major axis line
    fish(fishNum).minorXs = zeros(endframe-startframe,2); % X pairs for minor axis line
    fish(fishNum).minorYs = zeros(endframe-startframe,2); % Y pairs for minor axis line
    fish(fishNum).frameno = zeros(endframe-startframe,1);
    
% Copy initialization of fish for any additional fish
     if (fishNum > 1) 
         for j=fishNum:-1:2; fish(j) = fish(end); end 
     end

%% Prepare the data for analysis

% Compute the background
fprintf('Computing background. \n');
    bg = imdilate(grae, ones(1, 1, dil));

% Get the differences between each frame and the background
fprintf('Taking image differences. \n');
    df = imabsdiff(grae, bg);

    threshval = [-0.9 -0.7 -0.5 -0.3 -0.1 0.1];

    userpick = 99;

    while userpick == 99
    
% Take different thresholds of the background level
fprintf('Applying thresholds. \n');
    % Thresholds for closed door only IR
    thresh(1) = graythresh(df) + graythresh(df)* threshval(1);
    thresh(2) = graythresh(df) + graythresh(df)* threshval(2);
    thresh(3) = graythresh(df) + graythresh(df)* threshval(3);
    thresh(4) = graythresh(df) + graythresh(df)* threshval(4);
    thresh(5) = graythresh(df) + graythresh(df)* threshval(5);
    thresh(6) = graythresh(df) + graythresh(df)* threshval(6);

    % Thresholds for open door with IR
%     thresh(1) = graythresh(df) + graythresh(df)* 0.9;
%     thresh(2) = graythresh(df) + graythresh(df)* 0.7;
%     thresh(3) = graythresh(df) + graythresh(df)* 0.5;
%     thresh(4) = graythresh(df) + graythresh(df)* 0.3;
%     thresh(5) = graythresh(df) + graythresh(df)* 0.1;
%     thresh(6) = graythresh(df) + graythresh(df)* 0.0;

    lt = length(thresh); % There are presently 6 thresholds defined above
    N = sort(fix((endframe-startframe)*rand(1,lt))); % We'll display 'lt' random frames from our movie

if length(rango) < 3 % User did not provide which gray level they want
    % Show the User what we've got
    for pp = 1:lt % 1 2 3 4 5 6 
        figure(1);
        nrw = lt * (pp - 1);
        % Let's do 1 row
        for rw = 1+nrw:(lt + nrw)
            wb = (df >= thresh(pp) * 255);
            subplot(lt,lt,rw); imshow(wb(:,:,N(rw-nrw))); 
        end
    end
    
    userpick = input('Which threshold? If all are poopy, use 99: ');
    close(1);
    
end
    if userpick == 99 
        fprintf('The old numbers were: ');
        threshval
        threshval = input('Gimme the new numbers! e.g. [1 2 3 4 5 6]: ');
    end

    end
    
% Make the image binary - white for above thresh and 0 black below
    bw = (df >= thresh(userpick) * 255);
    
%% Get user picks for first fish positions

fprintf('Click initial fish position(s) in Figure 1. Click largest fish to smallest fish. \n');
    
    % Show the "end" image
       figure(1); imshow(grae(:,:,end)); 
    % Now click the fish
    
        fishi = round(ginput(fishNum)); % These are the x and y integers from the clicks on the video frame

        hormax = length(bw(:,1,1)); % hormax and vermax are the size of the video frame
        vermax = length(bw(1,:,1));
% Import initial clicks.        
for kk =1:fishNum
     fish(kk).x(endframe) = fishi(kk,1); 
     fish(kk).y(endframe) = fishi(kk,2); 
end
        

%% Track the fish        
        
stationaryfish = 0;

for jj = endframe-startframe-1:-1:1       % Cycle through every frame of our sample
    
    for kk = 1:fishNum % Cycle tracking for each fish

% We cheat - we only track in a box around where the fish was previously located.        
L = logical(bw(max([1 fishi(kk,2)-boxy(1)]):min([fishi(kk,2)+boxy(1) hormax]), max([1 fishi(kk,1)-boxy(2)]):min([fishi(kk,1)+boxy(2) vermax]), jj));
        
        s = regionprops(L, 'Area', 'Centroid', 'Orientation', 'MajorAxisLength', 'MinorAxisLength'); % Get the moving regions
        area_vector = [s.Area];
        % Take the largest blob
        [~, blobidx] = max(area_vector); 
        
            if isempty(blobidx) == 1 % We have no blob so we need to do something...
                
                if kk == stationaryfish % User already told us that this fish is stationary, so use previous location
                    
                    % Take the previous known location and other data of the fish
                    fish(kk).x(jj) = fish(kk).x(jj+1); 
                    fish(kk).y(jj) = fish(kk).y(jj+1);
                    fish(kk).majorLength(jj) = fish(kk).majorLength(jj+1);
                    fish(kk).minorLength(jj) = fish(kk).minorLength(jj+1);
                    fish(kk).majorXs(jj,:) = fish(kk).majorXs(jj+1,:);
                    fish(kk).majorYs(jj,:) = fish(kk).majorYs(jj+1,:);
                    fish(kk).minorXs(jj,:) = fish(kk).minorXs(jj+1,:);
                    fish(kk).minorYs(jj,:) = fish(kk).minorYs(jj+1,:);   
                    
                end
                
                if kk ~= stationaryfish % Lost fish is not a known stationary fish
                
                % Show the fish
                
                    figure(1); clf; imshow(grae(:,:,jj)); hold on;
                    plot(fish(kk).x(jj:min([endframe, jj+20])), fish(kk).y(jj:min([endframe, jj+20])), 'g*');
                    plot(fish(kk).x(jj+1), fish(kk).y(jj+1), 'go', 'MarkerSize', 4); % Last known location of lost fish
                    
                    yn = input('Is this fish stationary? 1 for yes. \n');
                    
                    if ~isempty(yn)
                        if yn == 1
                            
                    stationaryfish = kk; 
                    
                    fish(kk).x(jj) = fish(kk).x(jj+1); 
                    fish(kk).y(jj) = fish(kk).y(jj+1);
                    fish(kk).majorLength(jj) = fish(kk).majorLength(jj+1);
                    fish(kk).minorLength(jj) = fish(kk).minorLength(jj+1);
                    fish(kk).majorXs(jj,:) = fish(kk).majorXs(jj+1,:);
                    fish(kk).majorYs(jj,:) = fish(kk).majorYs(jj+1,:);
                    fish(kk).minorXs(jj,:) = fish(kk).minorXs(jj+1,:);
                    fish(kk).minorYs(jj,:) = fish(kk).minorYs(jj+1,:);   
                        
                        end
                    end
                    
                    if isempty(yn) == 0                    
                
                [xTMP, yTMP] = ginput(1); % Get our click
                fish(kk).x(jj) = round(xTMP); fish(kk).y(jj) = round(yTMP); % We need integers for indices

                    % We don't have these data because no blob was found.  
                    fish(kk).majorLength(jj) = [];
                    fish(kk).minorLength(jj) = [];
                    fish(kk).majorXs(jj,:) = [];
                    fish(kk).majorYs(jj,:) = [];
                    fish(kk).minorXs(jj,:) = [];
                    fish(kk).minorYs(jj,:) = [];   
                    
                    end
                    
                end % Lost fish is not stationary

            end % Lost the fish
        
        if isempty(blobidx) == 0 % Yay!  We have a blob.
        % Save the x (centroid) position for the largest blob
            fish(kk).x(jj) = round(s(blobidx).Centroid(1)) + fishi(kk,1)-boxy(1);
            fish(kk).y(jj) = round(s(blobidx).Centroid(2)) + fishi(kk,2)-boxy(2);
            
            % Save the rest of the data
                fish(kk).orient(jj) = s(blobidx).Orientation;
                fish(kk).majorLength(jj) = s(blobidx).MajorAxisLength;
                fish(kk).minorLength(jj) = s(blobidx).MinorAxisLength;
                    XX = (s(blobidx).MajorAxisLength * cosd(s(blobidx).Orientation))/2;
                    YY = (s(blobidx).MajorAxisLength * sind(s(blobidx).Orientation))/2;
                fish(kk).majorXs(jj,:) = [(s(blobidx).Centroid(1)-XX) (s(blobidx).Centroid(1)+XX)];
                fish(kk).majorYs(jj,:) = [(s(blobidx).Centroid(2)+YY) (s(blobidx).Centroid(2)-YY)];
                    XX = (s(blobidx).MinorAxisLength * cosd(s(blobidx).Orientation-90))/2;
                    YY = (s(blobidx).MinorAxisLength * sind(s(blobidx).Orientation-90))/2;
                fish(kk).minorXs(jj,:) = [(s(blobidx).Centroid(1)-XX) (s(blobidx).Centroid(1)+XX)];
                fish(kk).minorYs(jj,:) = [(s(blobidx).Centroid(2)+YY) (s(blobidx).Centroid(2)-YY)];
                fish(kk).frameno(jj) = jj+startframe-1;
            
        end
        
        fishi(kk,:) = [fish(kk).x(jj) fish(kk).y(jj)]; % fishi is the current location - set it for the next iteration
        
    end
             if rem(jj,10) == 0 % Update the track every 10 frames
                figure(3); clf; imshow(grae(:,:,jj)); hold on; 
                plot(fish(1).x(jj:min([endframe, jj+50])), fish(1).y(jj:min([endframe, jj+50])), 'r*');
                if fishNum == 2
                    plot(fish(2).x(jj:min([endframe, jj+50])), fish(2).y(jj:min([endframe, jj+50])), 'b*');
                end
                pause(0.1); 
             end

    % If fish are in the same location, we have a problem that needs to be
    % clicked
        if fishNum == 2
        if fish(1).x(jj) == fish(2).x(jj) && fish(1).y(jj) == fish(2).y(jj)
            
            figure(1); clf; imshow(grae(:,:,jj)); hold on; 
                plot(fish(1).x(jj+1), fish(1).y(jj+1), 'r*'); 
                plot(fish(2).x(jj+1), fish(2).y(jj+1), 'b*'); 
                fprintf('Use separate clicks, large fish first. \n');
            fishi = round(ginput(fishNum));
                         
        end
        end
    
    
end        


%% When all is said and done, transfer the images to our output variable    
fprintf('Taking final images. \n');
for i = (endframe-startframe):-1:1 
    imdata(i).bw = bw(:,:,i); 
    imdata(i).gray = grae(:,:,i);
end


%% Award Winning Code
%for j=16:2:length(asdf(1).x); plot(asdf(1).x(j-15:j), asdf(1).y(j-15:j), 'g*'); axis([0,350,0,350]); hold on; plot(asdf(2).x(j-15:j), asdf(2).y(j-15:j), 'm*'); pause(0.01); hold off; end;
%
% colr(1,:)='r*'; colr(2,:)='b*'; colr(3,:)='m*'; colr(4,:)='g*'; colr(5,:)='c*'; colr(6,:)='k*'; figure; set(gcf, 'Position', [838 639 759 541]);
% for i=1:length(fish(1).x); clf; imshow(im(i).gray); hold on; for j=1:length(fish); plot(fish(j).x(i),fish(j).y(i),colr(j,:)); end; pause(0.01); end;
%
% colrl(1,:)='r-'; colrl(2,:)='b-'; colrl(3,:)='m-'; colrl(4,:)='g-'; colrl(5,:)='c-'; colrl(6,:)='k-';
% colr(1,:)='r*'; colr(2,:)='b*'; colr(3,:)='m*'; colr(4,:)='g*'; colr(5,:)='c*'; colr(6,:)='k*'; figure;
% for i=21:length(fish(1).x); clf; imshow(im(i).gray); hold on; for j=1:length(fish); plot(fish(j).x(i-20:i),fish(j).y(i-20:i),colrl(j,:)); plot(fish(j).x(i),fish(j).y(i),colr(j,:)); end; pause(0.01); end;
