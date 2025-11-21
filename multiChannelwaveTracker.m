function out = multiChannelwaveTracker(dataStruct, vtims)
% Usage out = multiChannelwaveTracker(data, videotimes, Fs)
% This tracks the EOD frequency of a single fish
% providing an 'instantaneous' value of the amplitude and frequency
% for each video frame.

numChans = length(dataStruct);
Fs = round(1/dataStruct(1).interval);

% Filter the low frequency information out of the data
cutoffreq = 320; % Default was 200 Hz

[b,a] = butter(3,cutoffreq/(Fs/2),'high');
dataStruct = filtfilt(b,a,dataStruct);
wid=1; % Width (in seconds) of the FFT
de = 50; % Width in Hz that we can drift from original frequency

% Time sequence for electrical data
etim = 1/Fs:1/Fs:dataStruct(1).length/Fs;

%% Plot all channels


%% Get clicks

numFreqs = input('Number of frequencies?');

for kk = numFreqs:-1:1

    % Get target frequency based on user click
    figure(27); clf; 
    freqres = 2; % multiple of 2 between 2 and 32
    specgram(dataStruct, 1024*freqres, Fs, [], round(1024*freqres*0.85));
    ylim([200 1100]); colormap('HOT'); clim([-20 30]);
    fprintf('Click on the EOD frequency.\n');
    [~, precisefreq] = ginput(1);
    rango = [precisefreq-de, precisefreq+de];

end


%% Trace each frequency across all channels
% Our strategy is to find the electrode with the highest power at that
% frequency.  We need to add a threshold for when the fish is not in the
% grid.

for jj = numChans:-1:1

    
end

% Parallel processor way
parfor j = 1:length(videotimes)
   
    % tt = find(etim < vtims(j) & etim > vtims(j)-wid)
    tmp = fftmachine(dataStruct(etim < vtims(j) & etim > vtims(j)-wid), Fs);
    listOfreqs = tmp.fftfreq(tmp.fftfreq > rango(1) & tmp.fftfreq < rango(2));
        [Pamp(j), idx] = max(tmp.fftdata(tmp.fftfreq > rango(1) & tmp.fftfreq < rango(2)));
        Pfreq(j) = listOfreqs(idx);
        Ptim(j) = vtims(j);
end

    out.amp = Pamp;
    out.freq = Pfreq;
    out.tim = Ptim;





%% Plot the data to make the user happy - comment this out if you don't need happiness
figure(1); clf; 
    ax(1) = subplot(211);
    specgram(dataStruct, 4096, Fs, [], 4000); ylim([100 1000]); clim([-10 30]); colormap('HOT');
    hold on;
    plot(out.tim, out.freq, 'g.', 'MarkerSize', 2);
    ax(2) = subplot(212);
    plot(out.tim, out.amp, 'g.', 'MarkerSize', 1);
    linkaxes(ax, 'x'); xlim([0, out.tim(end)]);
    
    

