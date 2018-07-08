function out = singletrackEOD(data, totalframes, Fs)
% Usage out = singletrackEOD(chan1, Fs)
% This tracks the EOD frequency of a single fish
% providing an 'instantaneous' value for each video
% frame.

[b,a] = butter(3,400/(Fs/2),'high');

data = filtfilt(b,a,data);
wid=1;
de = 50;

etim = 1/Fs:1/Fs:length(data)/Fs;

    vidrate = max(etim)/totalframes;
    
vtims = vidrate:vidrate:totalframes*vidrate;

tmp = fftmachine(data, Fs); 
[~, idx] = max(tmp.fftdata);
precisefreq = tmp.fftfreq(idx);
rango = [precisefreq-de, precisefreq+de];

% Get frequencies

for j = totalframes:-1:1
   
    % tt = find(etim < vtims(j) & etim > vtims(j)-wid)
    tmp = fftmachine(data(etim < vtims(j) & etim > vtims(j)-wid), Fs);
    listOfreqs = tmp.fftfreq(tmp.fftfreq > rango(1) & tmp.fftfreq < rango(2));
        [out.amp(j), idx] = max(tmp.fftdata(tmp.fftfreq > rango(1) & tmp.fftfreq < rango(2)));
        out.freq(j) = listOfreqs(idx);
        out.tim(j) = vtims(j);
end

    
% figure(1); clf; 
% subplot(211);
% specgram(data, 4096, Fs, [], 4000); ylim([200 1000]); caxis([-10 30]); colormap('HOT');
% hold on;
% plot(out.tim, out.freq, 'g-', 'LineWidth', 1);
% subplot(212);
% plot(out.tim, out.amp, 'g-', 'LineWidth', 2);


