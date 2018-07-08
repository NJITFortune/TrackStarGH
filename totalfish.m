function [out, im, data] = totalfish(old, img)
% totalfish [out newim] = totalfish(data, im) where data and im are the outputs of traxer
% This script checks, cleans, and processes the data. It also takes clips around the centroid to make
% fish-o-centric movies.
% 

data = old;

%% Preparations
% How many fish do we have?
fishNum = length(data);

% This is a stupid way to get the x and y dimensions of the images
    %pp = img(1).gray;
    xlen = length(img(1).gray(1,:));
    ylen = length(img(1).gray(:,1));
    
% We make the area that we are going to snip be 1.8x larger than the
% longest length of the fish, and we'll pad the image to be 2.5x larger.
    maxlen = round(1.6 * max([data(1).majorLength])); 
    padlen = round(2.5 * max([data(1).majorLength])); 

% A big old white image.
    background = uint8(zeros((2*padlen)+ylen,(2*padlen)+xlen));

%% Fix errors in the data

for j = fishNum:-1:1; % This loop could be in the main analysis loop...
    
% Error type 1 - find corner points, if any, and fix as for type 1.

    zCs = find(data(j).x == 0 & data(j).y == 0);
    
    if ~isempty(zCs); 
        fprintf('Found %i zeros for fish %i. \n', length(zCs), j);
        
        % For now we will assume that the first and last values are not 0.
        % Too lazy to write the code that solves those special cases.
        
        for kk = 1:length(zCs); % for each zero
            % Get previous value
            preX = data(j).x(zCs(kk)-1); 
            preY = data(j).y(zCs(kk)-1);
            % Get next real actual value, skipping zeros
            zippX = 0; zippY = 0; adm = 1;
            while zippX == 0;
                zippY = data(j).y(zCs(kk)+adm); 
                zippX = data(j).x(zCs(kk)+adm); 
                adm = adm+1;
            end;
            data(j).x(zCs(kk)) = round((preX+zippX)/2);
            data(j).y(zCs(kk)) = round((preY+zippY)/2);
                    % These values are definitely incorrect, so we erase them.
                     data(j).majorLength(zCs(kk)) = 0;
                     data(j).minorLength(zCs(kk)) = 0;
                     data(j).majorXs(zCs(kk),:) = [0 0];
                     data(j).majorYs(zCs(kk),:) = [0 0];
                     data(j).minorXs(zCs(kk),:) = [0 0];
                     data(j).minorYs(zCs(kk),:) = [0 0];                    

        end;
    end;

% Error type 2 - find clicked points and add guestimates of centroid.

        zCs = find(data(j).majorLength == 0);
        
    if ~isempty(zCs); 
        fprintf('Found %i missing data for fish %i. \n', length(zCs), j);

        % For now we will assume that the first and last values are not 0.
        % Too lazy to write the code that solves those special cases.
        
        for kk = 1:length(zCs); % for each zero
            % Get previous values for each piece of data
            preMajL = data(j).x(zCs(kk)-1); preMinL = data(j).y(zCs(kk)-1);
            preMajX = data(j).majorXs(zCs(kk)-1,:);
            preMinX = data(j).minorXs(zCs(kk)-1,:);
            preMajY = data(j).majorYs(zCs(kk)-1,:);
            preMinY = data(j).minorYs(zCs(kk)-1,:);
            
            % Get next real actual value, skipping zeros
            nMajL = 0; adm = 1;
            while nMajL == 0;
                nMajX = data(j).majorXs(zCs(kk)+adm,:);
                nMinX = data(j).minorXs(zCs(kk)+adm,:);
                nMajY = data(j).majorYs(zCs(kk)+adm,:);
                nMinY = data(j).minorYs(zCs(kk)+adm,:);
                nMinL = data(j).x(zCs(kk)+adm);
                nMajL = data(j).x(zCs(kk)+adm);
                adm = adm+1;
            end;
            
data(j).majorXs(zCs(kk),:) = [round((preMajX(1)+nMajX(1))/2) round((preMajX(2)+nMajX(2))/2)];
data(j).minorXs(zCs(kk),:) = [round((preMinX(1)+nMinX(1))/2) round((preMinX(2)+nMinX(2))/2)];
data(j).majorYs(zCs(kk),:) = [round((preMajY(1)+nMajY(1))/2) round((preMajY(2)+nMajY(2))/2)];
data(j).minorYs(zCs(kk),:) = [round((preMinY(1)+nMinY(1))/2) round((preMinY(2)+nMinY(2))/2)];

        end;
        
    end;


% Error type 3 - filter the track and check for misguided tracks.

% medfilt1();

end;


%% Rotate the fish and get head direction

for j = fishNum:-1:1;

    % This is just shifting the orientation measure so that "0" points
    % north/up. traxer gives a range of -90 to 90.
    rot = data(j).orient(:) - 90; 
    
    % The problem with orient is that you get the 0/360 jump - this
    % "unwinds" that. This algorithm is dumber than that, however,
    % as it simply looks for jumps of 120 degrees or more. 
    for pp = 1:length(rot)-1;
        if ((rot(pp) - rot(pp+1)) > 120); % Really? 
            rot(pp+1:end) = rot(pp+1:end) + 180; 
        end;
        if ((rot(pp) - rot(pp+1)) < -120); 
            rot(pp+1:end) = rot(pp+1:end) - 180; 
        end;
    end;

    % Put these data into the output structure.
    out(j).rot = rot;
    
    %     % Now we build the fish-o-centric video data
        xmin = round(data(j).x(end) - maxlen) + padlen;
        xmax = round(data(j).x(end) + maxlen) + padlen;
        ymin = round(data(j).y(end) - maxlen) + padlen;
        ymax = round(data(j).y(end) + maxlen) + padlen;
 
     imagepadded = uint8(background);
%     %bwpadded = background;
     
     imagepadded(padlen:padlen+ylen-1,padlen:padlen+xlen-1) = img(end).gray;
%     %bwpadded(padlen:padlen+ylen-1,padlen:padlen+xlen-1) = img(i).bw;
% 
     tmp_fish(:,:) = imagepadded(ymin:ymax,xmin:xmax);
     tmp_rot(:,:) = imrotate(tmp_fish(:,:), -rot(end), 'crop');
     figure(2); imshow(tmp_rot); hold on; plot(length(tmp_rot(:,1))/2, length(tmp_rot(1,:))/2, 'yo');

     headsup = input('1 for fish head up and 2 for down: ');
    
     if headsup == 1;
            [hd(2), tmpidx] = max([data(j).majorYs(end)]);
            hd(1) = data(j).majorXs(tmpidx);
     end;
     if headsup == 2; 
         rot = rot + 180;
         [hd(2), tmpidx] = min([data(j).majorYs(end)]);
         hd(1) = data(j).majorXs(tmpidx);
     end;
    
% Ask the user for a click on the head so that we know where the head might
% be. 
%    figure(1); imshow(img(end).gray); 
%        hold on; 
%        plot(data(j).majorXs(end,:), data(j).majorYs(end,:),'r*');
%        hd = ginput(1);

%   if hd(1) < data(j).x(end); rot = rot + 180; end; % was y
    
%% Perform analysis and generate rotated video    

for i = length(img):-1:1;

% Find which of the two major axis vertices are closer to the last marked head.
% htl is 1 and 2 for each of the major axis ends.
for htl = 1:2; 
    dst(htl)=pdist([hd(1),hd(2);data(j).majorXs(i,htl),data(j).majorYs(i,htl)]);
end;
        
    % find the shorter distance of the two - this should work most of the time    
    [~, closeIdx] = min(dst); 

    out(j).head(i,:) = [data(j).majorXs(i,closeIdx) data(j).majorYs(i,closeIdx)];
    hd = out(j).head(i,:); % keep track of the last spot for next iteration

if i == length(img); % First image in the analysis    
    
end;

    % If we can, we calculate the velocity by comparing to the distances of
    % the prior frame centroid and the current centroid to the current frame head.

    if i > 1;
        curlen = pdist([out(j).head(i,:); data(j).x(i) data(j).y(i)]);
        prelen = pdist([out(j).head(i,:); data(j).x(i-1) data(j).y(i-1)]);
        if prelen > curlen; fr = 1; end;
        if prelen <= curlen; fr = -1; end;
        out(j).vel(i-1) = fr * pdist([data(j).x(i) data(j).y(i); data(j).x(i-1) data(j).y(i-1)]);
    end;

    % Now we build the fish-o-centric video data
        xmin = round(data(j).x(i) - maxlen) + padlen;
        xmax = round(data(j).x(i) + maxlen) + padlen;
        ymin = round(data(j).y(i) - maxlen) + padlen;
        ymax = round(data(j).y(i) + maxlen) + padlen;

    imagepadded = uint8(background);
    %bwpadded = background;
    
    imagepadded(padlen:padlen+ylen-1,padlen:padlen+xlen-1) = img(i).gray;
    %bwpadded(padlen:padlen+ylen-1,padlen:padlen+xlen-1) = img(i).bw;

    % This is error checking to correct the rare rounding error
    if i < length(img);
        if ymax-ymin > maxx; ymax = ymax-1; end;
        if xmax-xmin > maxy; xmax = xmax-1; end;
    end;

    im(j).fish(:,:,i) = imagepadded(ymin:ymax,xmin:xmax);
    im(j).rotFish(:,:,i) = imrotate(im(j).fish(:,:,i), -rot(i), 'crop');

    figure(2); clf; imshow(im(j).rotFish(:,:,i)); drawnow;
    % Get the size of the rotated image for error checking (above)
    if i == length(img);
        maxx = length(im(j).fish(:,1,1))-1;
        maxy = length(im(j).fish(1,:,1))-1;
    end;
end;

end;
im = 0;
%%%% USE THIS TO SHOW THE ROTATED DATA for both fish at the same time
%
% for i=1:length(img); 
%    subplot(1,2,1); imshow(img(1).rotFish(:,:,i)); 
%    subplot(1,2,2); imshow(img(2).rotFish(:,:,i));     
% pause(0.05); end

% figure(1); for i=1:length(im); 
% subplot(1,2,1); imshow(fish(1).rotFish(:,:,i)); 
% subplot(1,2,2); imshow(fish(2).rotFish(:,:,i)); pause(0.05); end;



% Award Winning Code
%v = [-1 -12];
%theta = 5;
%R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
%vR = v*R;



