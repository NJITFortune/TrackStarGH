function out = brachyJAABA(in, im, eodata)
% Usage out = proJAABA(spdata, im, eodata)
% Where spdata is the 'spadata' from fun or swimmy
% Should be either spadata.A or spadata.B
% im is 'im' from fun or swimmy
% Should be either im.A or im.B

%% Preset values

xmid = 375; ymid = 300; scale2mm = 2.344;

%% Transfer the easy stuff - not angle data

% Our video output is 320 x 256 pixels, and the dimensions of the tank
% are 750 x 600 mm.  We need to subtract 375 x 300 and multiply by 2.344.

for fishnum = 1:length(in)
    
    out(fishnum).sex = 'u';
    out(fishnum).fps = 30;  % CHECK TO SEE!!!
    out(fishnum).dt = 1/out(fishnum).fps * ones(1, length(in(fishnum).x)-1); % Time between frames, one less than total number of frames

    out(fishnum).firstframe = 1; % Hmm?
    out(fishnum).endframe = length(in(fishnum).x);
    out(fishnum).nframes = length(in(fishnum).x);

    out(fishnum).id = fishnum-1; % They start at zero, so should we.
    out(fishnum).off = 1 - out(fishnum).firstframe; % This is by definition, but have no idea how this is used.

% Copy the X and Y pixel data without modification

% out(j).x = in(j).x - xmid; % 0,0 in middle of tank
% out(j).y = in(j).y - ymid;  
    out(fishnum).x = in(fishnum).x; 
    out(fishnum).y = in(fishnum).y; 
    out(fishnum).a = in(fishnum).majorLength * 0.25;
    out(fishnum).b = in(fishnum).minorLength * 0.25;

% out(j).x_mm = (in(j).x - xmid) * scale2mm; % 0,0 in middle of tank
% out(j).y_mm = (in(j).y - ymid) * scale2mm;
    out(fishnum).x_mm = (in(fishnum).x) * scale2mm;
    out(fishnum).y_mm = (in(fishnum).y) * scale2mm;
    out(fishnum).a_mm = in(fishnum).majorLength * scale2mm * 0.25;
    out(fishnum).b_mm = in(fishnum).minorLength * scale2mm * 0.25;

% out(j).timestamps = []; % Placeholder - should be a value for each position

% Copy the orientation data into a temporary structure
    tmpfish(fishnum).orient = in(fishnum).orient;
    tmpfish(fishnum).majorXs = in(fishnum).majorXs;
    tmpfish(fishnum).majorYs = in(fishnum).majorYs;
    tmpfish(fishnum).majorLength = in(fishnum).majorLength;
    tmpfish(fishnum).minorXs = in(fishnum).minorXs;
    tmpfish(fishnum).minorYs = in(fishnum).minorYs;
    tmpfish(fishnum).minorLength = in(fishnum).minorLength;

end


%% Fish segregation check

if length(in) == 2 % we have two fish   
    
    figure(2); clf; imshow(im(1).gray, 'InitialMagnification', 300);
    hold on; plot(in(1).x(1), in(1).y(1), 'g*'); plot(in(2).x(1), in(2).y(1), 'm*');
    lrg = input('Which fish is larger (1 for green and 2 for magenta)? ');
    
        % vidplay(im, out, rango)
        scrollerdoller(in, im, [1 length(im)]);
        
        xs = input('Were there one or more fish switches? 0-no, numswaps: ');
        
        if xs > 0 
          
%             clsups = input('Give the centers of whatever closeups you need. ');
%             for zz = 1:length(clsups)
%               scrollerdoller(in, im, max([1, clsups-50]), min([clsups+50, length(im)]));
%             end
          
          for ii = 1:xs
 
            swpidx = input('Range to swap? : ');
            
            out(1).x(swpidx(1):swpidx(2)) = in(2).x(swpidx(1):swpidx(2));
            out(1).y(swpidx(1):swpidx(2)) = in(2).y(swpidx(1):swpidx(2));
            tmpfish(1).orient(swpidx(1):swpidx(2)) = in(2).orient(swpidx(1):swpidx(2));
            tmpfish(1).majorXs(swpidx(1):swpidx(2)) = in(2).majorXs(swpidx(1):swpidx(2));
            tmpfish(1).majorYs(swpidx(1):swpidx(2)) = in(2).majorYs(swpidx(1):swpidx(2));
            tmpfish(1).majorLength(swpidx(1):swpidx(2)) = in(2).majorLength(swpidx(1):swpidx(2));
            tmpfish(1).minorXs(swpidx(1):swpidx(2)) = in(2).minorXs(swpidx(1):swpidx(2));
            tmpfish(1).minorYs(swpidx(1):swpidx(2)) = in(2).minorYs(swpidx(1):swpidx(2));
            tmpfish(1).minorLength(swpidx(1):swpidx(2)) = in(2).minorLength(swpidx(1):swpidx(2));

            out(2).x(swpidx(1):swpidx(2)) = in(1).x(swpidx(1):swpidx(2));
            out(2).y(swpidx(1):swpidx(2)) = in(1).y(swpidx(1):swpidx(2));
            tmpfish(2).orient(swpidx(1):swpidx(2)) = in(1).orient(swpidx(1):swpidx(2));
            tmpfish(2).majorXs(swpidx(1):swpidx(2)) = in(1).majorXs(swpidx(1):swpidx(2));
            tmpfish(2).majorYs(swpidx(1):swpidx(2)) = in(1).majorYs(swpidx(1):swpidx(2));
            tmpfish(2).majorLength(swpidx(1):swpidx(2)) = in(1).majorLength(swpidx(1):swpidx(2));
            tmpfish(2).minorXs(swpidx(1):swpidx(2)) = in(1).minorXs(swpidx(1):swpidx(2));
            tmpfish(2).minorYs(swpidx(1):swpidx(2)) = in(1).minorYs(swpidx(1):swpidx(2));
            tmpfish(2).minorLength(swpidx(1):swpidx(2)) = in(1).minorLength(swpidx(1):swpidx(2));
          end
                                
        end
            
end

%% Rotate the animal

boxlen = 25; % This is the number of pixels to show around the animal

% For each fish
for fishnum = 1:length(in)

    % Plot the initial position to determine head
    figure; clf; 

    currXs = [max([1 out(fishnum).x(1)-boxlen]), min([320 out(fishnum).x(1)+boxlen])];
    currYs = [max([1 out(fishnum).y(1)-boxlen]), min([256 out(fishnum).y(1)+boxlen])];

    imshow(im(fishnum).gray(round(currYs(1)):round(currYs(2)),round(currXs(1)):round(currXs(2)))); 

    hold on;
  
    plot(out(fishnum).x(1)-currXs(1), out(fishnum).y(1)-currYs(1), 'ob', 'LineWidth', 2);
    plot(tmpfish(fishnum).majorXs(1,:)-currXs(1), tmpfish(fishnum).majorYs(1,:)-currYs(1), 'g*');
    plot(tmpfish(fishnum).majorXs(1,:)-currXs(1), tmpfish(fishnum).majorYs(1,:)-currYs(1), 'r-');
    plot([out(fishnum).x(1)-currXs(1) out(fishnum).x(1)-currXs(1)], [1 boxlen*2], 'm-');
    plot([1 boxlen*2], [out(fishnum).y(1)-currYs(1) out(fishnum).y(1)-currYs(1)], 'm-');
    
    truesize(1,[200 200])

fprintf('Please click in appropriate quadrent for the head. \n');

    [cx, cy] = ginput(1);
 
%     ddist(1) = pdist([cx, cy; in(j).majorXs(1,1)-currXs(1), in(j).majorYs(1,1)-currYs(1)]);
%     ddist(2) = pdist([cx, cy; in(j).majorXs(1,2)-currXs(1), in(j).majorYs(1,2)-currYs(1)]);
    
% end

% Calculate original orientation 

% X to the right of midpoint
if cx > boxlen
    % Y in bottom
    if cy < boxlen
        orient = 90 - tmpfish(fishnum).orient;
    end
    % Y in top 
    if cy > in(fishnum).y(1) - currYs(1)
        orient = 90 - tmpfish(fishnum).orient;
    end
end
% X to the left of midpoint
if cx < boxlen
    % Y in bottom
    if cy < boxlen
        orient = 270 - tmpfish(fishnum).orient;
    end
    % Y in top
    if cy > in(fishnum).y(1) - currYs(1)
        orient = 270 - tmpfish(fishnum).orient;
    end
end

% Correct subsequent problems
% 1) Unwind jumps 
% 2) Make new jumps at 360 (value goes above) and 0 (value goes below)

% UNWIND JUMPS
dd = diff(tmpfish(fishnum).orient);
jumpIDX = find(abs(dd) > 150); % esfesf

if isempty(jumpIDX) == 0
    for jumps = 1:length(jumpIDX)
        if dd(jumpIDX(jumps)) > 0
           orient(jumpIDX(jumps)+1:end) = orient(jumpIDX(jumps)+1:end)+180;
        end
        if dd(jumpIDX(jumps)) < 0
           orient(jumpIDX(jumps)+1:end) = orient(jumpIDX(jumps)+1:end)-180;
        end
    end
end

% ROTATE THROUGH 0/360

for alldat = 1:length(orient)
   if orient(alldat) > 360; orient(alldat:end) = orient(alldat:end) - 360; end
   if orient(alldat) < 0; orient(alldat:end) = orient(alldat:end) + 360; end

% Video
% imshow(im(alldat).gray); hold on; text(30, 30, num2str(orient(alldat)), 'Color', 'green');
% pause(0.01); clf;

end

    out(fishnum).theta = orient;
    out(fishnum).theta_mm = orient;
    
end

figure(2); clf;
        ax(1) = subplot(212); plot(out(1).theta,'.', 'MarkerSize', 3); % fish 1, unwound
    if length(in) == 2
        ax(2) = subplot(211); plot(out(2).theta,'.', 'MarkerSize', 3); % fish 2, unwound
    end
    linkaxes(ax, 'x');

    pause(0.5);
    
%% Get the arena corners

figure(27); clf; imshow(im(1).gray, 'InitialMagnification', 300);
hold on;
    fprintf('Click top right corner of tank \n');
        out(1).corners(1,:) = ginput(1);
        plot(out(1).corners(1,1), out(1).corners(1,2), 'g*');
    fprintf('Click bottom right corner of tank \n');
        out(1).corners(2,:) = ginput(1);
        plot(out(1).corners(:,1), out(1).corners(:,2), 'g-*');
    fprintf('Click bottom left corner of tank \n');
        out(1).corners(3,:) = ginput(1);
        plot(out(1).corners(:,1), out(1).corners(:,2), 'g-*');
    fprintf('Click top left corner of tank \n');
        out(1).corners(4,:) = ginput(1);
        plot(out(1).corners(:,1), out(1).corners(:,2), 'g-*');
        plot([out(1).corners(4,1), out(1).corners(1,1)], [out(1).corners(4,2), out(1).corners(1,2)], 'g-');    
    
        pause(0.2); % Should use force draw, but this works too.
    
%% Track the pulses of the fish

[b,a] = butter(3, 600/((1/eodata.interval)/2), 'high');
fpulses = filtfilt(b,a, eodata.values);


% If a single fish
if length(in) == 1
    tmp = singletrackEOD(eodata.values, length(im), 1/eodata.interval);
    out.eodfreq = tmp.freq;
    out.eodtim = tmp.tim;
    out.eodamp = tmp.amp;
    clear tmp;
end

% If two fish, we need to figure out who is who...

if length(in) == 2
    
    tmp = dualtrackEOD(eodata.values, length(im), 1/eodata.interval, species);

    if lrg ~= 0 % If the fish are different size...
        
     % Hopefully the two fish will have different amplitude EODs   
     pleaseplease = ttest(tmp(1).amp, tmp(2).amp);

    if pleaseplease == 1
        
        if mean(tmp(1).amp) > mean(tmp(2).amp)
            if lrg == 1
                out(1).eodfreq = tmp(1).freq; out(1).eodtim = tmp(1).tim; out(1).eodamp = tmp(1).amp;
                out(2).eodfreq = tmp(2).freq; out(2).eodtim = tmp(2).tim; out(2).eodamp = tmp(2).amp;
            end
            if lrg == 2
                out(2).eodfreq = tmp(1).freq; out(2).eodtim = tmp(1).tim; out(2).eodamp = tmp(1).amp;
                out(1).eodfreq = tmp(2).freq; out(1).eodtim = tmp(2).tim; out(1).eodamp = tmp(2).amp;
            end
        end
        
        if mean(tmp(1).amp) < mean(tmp(2).amp)
            if lrg == 2
                out(1).eodfreq = tmp(1).freq; out(1).eodtim = tmp(1).tim; out(1).eodamp = tmp(1).amp;
                out(2).eodfreq = tmp(2).freq; out(2).eodtim = tmp(2).tim; out(2).eodamp = tmp(2).amp;
            end
            if lrg == 1
                out(2).eodfreq = tmp(1).freq; out(2).eodtim = tmp(1).tim; out(2).eodamp = tmp(1).amp;
                out(1).eodfreq = tmp(2).freq; out(1).eodtim = tmp(2).tim; out(1).eodamp = tmp(2).amp;
            end
        end
                
    end
     
    if pleaseplease == 0
        
        fprintf('We have to re-assign the frequencies later.\n');
        
        out(1).eodfreq = tmp(1).freq; out(1).eodtim = tmp(1).tim; out(1).eodamp = tmp(1).amp;
        out(2).eodfreq = tmp(2).freq; out(2).eodtim = tmp(2).tim; out(2).eodamp = tmp(2).amp;
    
    end
        
    end
    
        if lrg == 0 % If the fish are different size...
            
        fprintf('We have to re-assign the frequencies later.\n');
        
        out(1).eodfreq = tmp(1).freq; out(1).eodtim = tmp(1).tim; out(1).eodamp = tmp(1).amp;
        out(2).eodfreq = tmp(2).freq; out(2).eodtim = tmp(2).tim; out(2).eodamp = tmp(2).amp;

        end
end
    
    % close(27); % Close the previous figure after we've tracked the frequency of the fish.

    figure(3); clf; subplot(211);
    specgram(eodata.values, 4096, 1/eodata.interval, [], 4000); 
    ylim(species); 
    
    caxis([-10 30]); colormap('HOT');
    hold on;
    plot(out(1).eodtim, out(1).eodfreq, 'g-', 'LineWidth', 1);
    subplot(212); hold on;
    plot(out(1).eodtim, out(1).eodamp, 'g-', 'LineWidth', 2); xlim([0 max(out(1).eodtim)]);

if length(in) == 2
    figure(3); subplot(211);
    plot(out(2).eodtim, out(2).eodfreq, 'c-', 'LineWidth', 1);
    figure(3); subplot(212); hold on;
    plot(out(2).eodtim, out(2).eodamp, 'c-', 'LineWidth', 2); xlim([0 max(out(2).eodtim)]);
    
end

pause(1); % An extra second to enjoy our handiwork.

%% Fix potential problems using sidewinder
% for j=1:length(out)
%     
%    out(j) = sidewinder(out(j), im); 
%     
% end


%% Replay spatial data
% figure(4); 
% for j=1:4:length(im) 
%     clf;
%     imshow(im(j).gray, 'InitialMagnification', 300); hold on;
%     
% % Fish 1
%     if out(1).theta(j) < 90
%         plot(out(1).x(j), out(1).y(j), 'g*', 'MarkerSize', 10);
%     end
%     if out(1).theta(j) >= 90 && out(1).theta(j) < 180
%         plot(out(1).x(j), out(1).y(j), 'c*', 'MarkerSize', 10);
%     end
%     if out(1).theta(j) >= 270
%         plot(out(1).x(j), out(1).y(j), 'r*', 'MarkerSize', 10);
%     end
%     if out(1).theta(j) >= 180 && out(1).theta(j) < 270
%         plot(out(1).x(j), out(1).y(j), 'm*', 'MarkerSize', 10);
%     end
%     plot([out(1).x(j) out(1).x(j)+(20*sin((out(1).theta(j)/360)*2*pi))], ...
%         [out(1).y(j) out(1).y(j)-(20*cos((out(1).theta(j)/360)*2*pi))], '-g');
% % Fish 2
% if length(in) == 2
%     if out(2).theta(j) < 90
%         plot(out(2).x(j), out(2).y(j), 'go', 'MarkerSize', 10);
%     end
%     if out(2).theta(j) >= 90 && out(2).theta(j) < 180
%         plot(out(2).x(j), out(2).y(j), 'co', 'MarkerSize', 10);
%     end
%     if out(2).theta(j) >= 270
%         plot(out(2).x(j), out(2).y(j), 'ro', 'MarkerSize', 10);
%     end
%     if out(2).theta(j) >= 180 && out(2).theta(j) < 270
%         plot(out(2).x(j), out(2).y(j), 'mo', 'MarkerSize', 10);
%     end
%     plot([out(2).x(j) out(2).x(j)+(20*sin((out(2).theta(j)/360)*2*pi))], ...
%         [out(2).y(j) out(2).y(j)-(20*cos((out(2).theta(j)/360)*2*pi))], '-m');
% end
% 
%     text(30,30, num2str(out(1).theta(j)), 'Color', 'w');
%     text(100,30, num2str(j), 'Color', 'g');
%     
%     text(30,40, '*', 'Color', 'r', 'FontSize', 18);
%     text(40,40, '*', 'Color', 'g', 'FontSize', 18);
%     text(30,50, '*', 'Color', 'm', 'FontSize', 18);
%     text(40,50, '*', 'Color', 'c', 'FontSize', 18);
% 
%     pause(0.05);  
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    imshow(im(j).gray, 'InitialMagnification', 'fit'); hold on;
%    plot(in(1).x(j), in(1).y(j), 'w*'); text(30,30, num2str(out(1).theta(j)), 'Color', 'w');


% A video replay
% for j=1:4:length(im.A) 
%     imshow(im.A(j).gray, 'InitialMagnification', 'fit'); hold on;
%     plot(spadata.A.x(j), spadata.A.y(j), 'g*'); text(30,30, num2str(qwer.theta(j)), 'Color', 'w');
%     pause(0.05);  
% end
%     plot(spadata.A.x(j), spadata.A.y(j), 'g*'); text(30,30, num2str(qwer.theta(j)), 'Color', 'w');

% A line plot replay of single fish in two simultaneous tanks
% figure(4); clf
% hold on;
% plot(jdat.A.corners(:,1), jdat.A.corners(:,2), 'b*');
% plot(jdat.B.corners(:,1), jdat.B.corners(:,2), 'm*');
% for j=1:10:length(spadata.A.x); plot(spadata.A.x(1:j), spadata.A.y(1:j), 'b', spadata.B.x(1:j), spadata.B.y(1:j), 'm'); pause(0.1); end;



%% JAABA

% Structure Field	Description
%- x	x-coordinate of the animal in pixels (1 x nframes).
%- y	y-coordinate of the animal in pixels (1 x nframes).
%- theta	Orientation of the animal (head) (1 x nframes).
% a	1/4 of the major-axis length in pixels (1 x nframes).
% b	1/4 of the minor-axis length in pixels (1 x nframes).
%- nframes	Number of frames in the trajectory of the current animal (scalar).
%- firstframe	First frame of the animal's trajectory (scalar).
%- endframe	Last frame of the animal's trajectory (scalar).
%- off	Offset for computing index into x, y, etc. Always equal to 1 - firstframe (scalar).
%- id	Identity number of the trajectory (scalar).
%- x_mm	x-coordinate of the animal in mm (1 x nframes).
%- y_mm	y-coordinate of the animal in mm (1 x nframes).
%- theta_mm	Orientation of the animal in real coordinates. This is often the same as theta, if no transformation other than translation and scaling is performed between pixels and real coordinates (1 x nframes).
% a_mm	1/4 of the major-axis length in mm (1 x nframes).
% b_mm	1/4 of the major-axis length in mm (1 x nframes).
%- sex	Sex of the animal. Can be just one value ('M' or 'F' or '?') or a cell array of 'M' and 'F' giving the sex for each frame. The size of the cell array should be nframes.
%- dt	Difference in timestamps of the current frame and next frame, in seconds (1 x nframes-1).
%- fps	Average frames-per-second (scalar).
% timestamps	Timestamp of each frame (optional), in days (1 x nframes). 

%% JAABA Larvae (we don't use)
% area	Larvae only. Area of the larva in pixels (1 x nframes).
% xcontour	Larvae only. x-coordinates of the contour of the larva in pixels (1 x nframes cell array).
% ycontour	Larvae only. y-coordinates of the contour of the larva in pixels (1 x nframes cell array).
% xspine	Larvae only. x-coordinates of the spine fit to the larva in pixels (nspinepts=11 x nframes).
% yspine	Larvae only. y-coordinates of the spine fit to the larva in pixels(nspinepts=11 x nframes).
% area_mm	Larvae only. Area of the larva in mm2 (1 x nframes).
% xspine_mm	Larvae only. x-coordinates of the spine fit to the larva in mm (nspinepts=11 x nframes).
% yspine_mm	Larvae only. y-coordinates of the spine fit to the larva in mm(nspinepts=11 x nframes).
