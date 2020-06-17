function out = ritaclick(in)
% out = ritaclick(in);
% run in = VideoReader('/Volumes/Cgate/Data2019/4_12_2019/individual_trials/trial1414_000hz/BammBamm_1414_000Hz.mov');
% before using this script.
% save: save filename.mat out


numframes = in.FrameRate * in.Duration; % Calculate total frames in video

fprintf('There are %i frames./n', numframes);

figure(1);

for j=2:2:numframes % for the entire video
% for j=2:2:20 % For testing
   
   currframe = read(in, j); % The reads the j frame of the video file

   clf; imshow(currframe); hold on; % Clear the figure, show the frame

   [out.shuttlex(j/2), out.shuttley(j/2)] = ginput(1); % User click the shuttle
        plot(out.shuttlex(j/2), out.shuttley(j/2), '.m', 'MarkerSize', 12); drawnow;

        [out.fishx(j/2), out.fishy(j/2)] = ginput(1); % User click the fish
        plot(out.fishx(j/2), out.fishy(j/2), '.g', 'MarkerSize', 12); drawnow;
        
   out.tim(j/2) = (1/in.FrameRate) * j; % Add the time of the current frame

   pause(0.2); % Allows the user to see the second click.  Could be 0.1 for faster work.

    if rem(j / 20) == 0

       fprintf('Current Frame is %i./n', j);
       save temp.mat out

    end
   
end

% Plot the click data when done for sanity check.
figure(2); clf; hold on;
    plot(out.tim, out.shuttlex, '*-m');
    plot(out.tim, out.fishx, '*-g');