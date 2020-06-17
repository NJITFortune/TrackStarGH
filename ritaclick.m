function out = ritaclick(in)

% run in = VideoReader('/Volumes/Cgate/Data2019/4_12_2019/individual_trials/trial1414_000hz/BammBamm_1414_000Hz.mov');

numframes = in.FrameRate * in.Duration;

figure(1);

for j=2:2:numframes % for the entire video
% for j=2:2:20 % For testing
   
   currframe = read(in, j);
   clf; imshow(currframe); hold on;
   [out.shuttlex(j/2), out.shuttley(j/2)] = ginput(1);
        plot(out.shuttlex(j/2), out.shuttley(j/2), '.m', 'MarkerSize', 8); drawnow;
   [out.fishx(j/2), out.fishy(j/2)] = ginput(1);
        plot(out.fishx(j/2), out.fishy(j/2), '.g', 'MarkerSize', 8); drawnow;
   out.tim(j/2) = (1/in.FrameRate) * j;
   pause(0.2);
    
end

figure(2); clf; hold on;
    plot(out.tim, out.shuttlex, '*-m');
    plot(out.tim, out.fishx, '*-g');