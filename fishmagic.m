function [im,x,y] = fishmagic(filename, frameno, outfile)
% [im,x,y] = fishmagic(filename, frameno, outfile)
% where im is the stabilized fish
% x and y are the coordinates in the tube
% filename is the input tiff file from CamWare64
% frameno is the number of frames that you want to process
% outfile is the outpuf filename (be sure to use single quotes)
vidObj = VideoWriter(outfile);
open(vidObj);

for i=1:frameno;
    stk(:,:,i) = imadjust(imread(filename,i),[0 0.075],[]);
end;    

thresh = 0.67; % This is the grey level for black/white
yspace = 90;
xspace = [150 950];

% Threshold
    bw = (stk >= thresh * 2^16);

    imshow(bw(:,:,1));
    [x(1), y(1)] = ginput(1);
    
    im(:,:,1) = stk(y(1)-yspace:y(1)+yspace,x(1)-xspace(1):x(1)+xspace(2),1);

for i = 2:frameno    
    coop = bw(y(i-1)-20:y(i-1)+20,x(i-1)-60:x(i-1)+60,i);
    L = logical(coop);
    s = regionprops(L, 'Area', 'Centroid', 'Orientation', 'MajorAxisLength', 'MinorAxisLength');
    % This is our list of blob sizes. 
        area_vector = [s.Area];
    % Sort the list so that we can take the top N+2 largest spots
    %    [~,sizeIndex] = sort(area_vector(:),'descend');  
    % Get the largest blob
        [~, idx] = max(area_vector); 
    % Save the x and y (centroid) position for the largest
        x(i) = s(idx).Centroid(1) + x(i-1)-60;
        y(i) = s(idx).Centroid(2) + y(i-1)-20;
        
end;        
ff = 5;
x = medfilt1(x,ff); y = medfilt1(y,ff);

close all;
        
% Compute the background
%    bg = imdilate(stk, ones(1, 1, dil));

% Get the differences between each frame and the background
%    df = imabsdiff(stk, bg);

for i=1:frameno;

    im(:,:,i) = stk(y(i)-yspace:y(i)+yspace, x(i)-xspace(1):x(i)+xspace(2), i);

    imshow(im(:,:,i));
    writeVideo(vidObj,getframe(gca));
%    
%    tmp = imread(filename, i);
%    stk(:,:,i) = imadjust(tmp,[0.0 0.05],[]);
%    subplot(211); imshow(stk(:,:,i));
%    writeVideo(vidObj,getframe(gca));
    
%    df = imabsdiff(stk(:,:,i-1), stk(:,:,i));
%    df = imadjust(df,[0.02 0.2]);
%    figure(1); imshow(df);

%    writeVideo(vidObj,getframe(gca));
   
end;




close(vidObj);
%close(dObj);
