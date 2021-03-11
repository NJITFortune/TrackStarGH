function out = singletrackEOD(data, videotimes, Fs)
% Usage out = singletrackEOD(data, videotimes, Fs)
% This tracks the EOD frequency of a single fish
% providing an 'instantaneous' value of the amplitude and frequency
% for each video frame.

% Filter the low frequency information out of the data
[b,a] = butter(3,200/(Fs/2),'high');
data = filtfilt(b,a,data);
wid=1; % Width (in seconds) of the FFT
de = 50; % Width in Hz that we can drift from original frequency

% Time sequence for electrical data
etim = 1/Fs:1/Fs:length(data)/Fs;

%     vidrate = max(etim)/videotimes;
%     vtims = vidrate:vidrate:videotimes*vidrate;

    vtims = videotimes; % times, not length !!!

    % Get target frequency based on max amplitude
%     tmp = fftmachine(data, Fs); 
%     [~, idx] = max(tmp.fftdata);
%     precisefreq = tmp.fftfreq(idx);
%     rango = [precisefreq-de, precisefreq+de];

    % Get target frequency based on user click
    figure(27); clf; 
    freqres = 2; % multiple of 2 between 2 and 32
    specgram(data, 1024*freqres, Fs, [], round(1024*freqres*0.85));
    ylim([200 1100]); colormap('HOT'); caxis([-20 30]);
    fprintf('Click on the EOD frequency.\n');
    [~, precisefreq] = ginput(1);
    rango = [precisefreq-de, precisefreq+de];
    
% Get frequencies

% % Single processor way
% for j = totalframes:-1:1
%    
%     % tt = find(etim < vtims(j) & etim > vtims(j)-wid)
%     tmp = fftmachine(data(etim < vtims(j) & etim > vtims(j)-wid), Fs);
%     listOfreqs = tmp.fftfreq(tmp.fftfreq > rango(1) & tmp.fftfreq < rango(2));
%         [out.amp(j), idx] = max(tmp.fftdata(tmp.fftfreq > rango(1) & tmp.fftfreq < rango(2)));
%         out.freq(j) = listOfreqs(idx);
%         out.tim(j) = vtims(j);
% end

% Parallel processor way
parfor j = 1:length(videotimes)
   
    % tt = find(etim < vtims(j) & etim > vtims(j)-wid)
    tmp = fftmachine(data(etim < vtims(j) & etim > vtims(j)-wid), Fs);
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
    specgram(data, 4096, Fs, [], 4000); ylim([100 1000]); caxis([-10 30]); colormap('HOT');
    hold on;
    plot(out.tim, out.freq, 'g.', 'MarkerSize', 2);
    ax(2) = subplot(212);
    plot(out.tim, out.amp, 'g.', 'MarkerSize', 1);
    linkaxes(ax, 'x'); xlim([0, out.tim(end)]);
    
    

