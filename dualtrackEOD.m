function out = dualtrackEOD(data, totalframes, Fs, species)
% Usage out = singletrackEOD(data, totalframes, Fs, species)
% This tracks the EOD frequencies of two fish
% providing 'instantaneous' values for each video frame.
% data is the samples, e.g. V170720_030_Ch1.values
% totalframes is the number of video frames, which we can get from the
% trigger channel, e.g. V170720_030_Ch4.length
% Fs can be extracted from the data, e.g. 1/V170720_030_Ch1.interval
% species is the Fourier transform (specgram) range for that species, e.g.
% species = [200 800]; % Eigenmannia
% species = [500 1100]; % Apteronotus

% Filter the low frequency information out of the data
[b,a] = butter(3,200/(Fs/2),'high');
data = filtfilt(b,a,data);
wid=1; % Width (in seconds) of the FFT
% de = 50;

% Time sequence for electrical data
etim = 1/Fs:1/Fs:length(data)/Fs;

    %vidrate = max(etim)/totalframes;
    %vtims = vidrate:vidrate:totalframes*vidrate;

    vtims = totalframes; % times, not length !!!
    
% tmp = fftmachine(data, Fs); 
% 
% semilogy(tmp.fftfreq, tmp.fftdata); 
% pause(10);
% 
% [~, idx] = max(tmp.fftdata);
% precisefreq = tmp.fftfreq(idx);
% rango = [precisefreq-de, precisefreq+de];

% Ask the users for a range

figure(27); clf; 
specgram(data, 1024*32, Fs, [], round(1024*32*0.85));
ylim(species); colormap('HOT'); caxis([-20 30]);

[~, splitterfreq] = ginput(1);

% Get frequencies

% % Single core mechanism
% for j = totalframes:-1:1
%    
%     % tt = find(etim < vtims(j) & etim > vtims(j)-wid)
%     tmp = fftmachine(data(etim < vtims(j) & etim > vtims(j)-wid), Fs);
%    
%     % Lower fish
%     
%     listOfreqs = tmp.fftfreq(tmp.fftfreq < splitterfreq & tmp.fftfreq > species(1));
%         [out(1).amp(j), idx] = max(tmp.fftdata(tmp.fftfreq < splitterfreq & tmp.fftfreq > species(1)));
%         out(1).freq(j) = listOfreqs(idx);
%         out(1).tim(j) = vtims(j);
% 
%     % Higher fish
%     
%     listOfreqs = tmp.fftfreq(tmp.fftfreq > splitterfreq & tmp.fftfreq < species(2));
%         [out(2).amp(j), idx] = max(tmp.fftdata(tmp.fftfreq > splitterfreq & tmp.fftfreq < species(2)));
%         out(2).freq(j) = listOfreqs(idx);
%         out(2).tim(j) = vtims(j);
%         
% end

% Parallel processing mechanism
parfor j = 1:totalframes
   
    % tt = find(etim < vtims(j) & etim > vtims(j)-wid)
    tmp = fftmachine(data(etim < vtims(j) & etim > vtims(j)-wid), Fs);
   
    % Lower fish
    
    listOfreqs = tmp.fftfreq(tmp.fftfreq < splitterfreq & tmp.fftfreq > species(1));
        [Pamp1(j), idx] = max(tmp.fftdata(tmp.fftfreq < splitterfreq & tmp.fftfreq > species(1)));
        Pfreq1(j) = listOfreqs(idx);
        Ptim1(j) = vtims(j);

    % Higher fish
    
    listOfreqs = tmp.fftfreq(tmp.fftfreq > splitterfreq & tmp.fftfreq < species(2));
        [Pamp2(j), idx] = max(tmp.fftdata(tmp.fftfreq > splitterfreq & tmp.fftfreq < species(2)));
        Pfreq2(j) = listOfreqs(idx);
        Ptim2(j) = vtims(j);
        
end

    out(1).amp = Pamp1;
    out(1).freq = Pfreq1;
    out(1).tim = Ptim1;
    out(2).amp = Pamp2;
    out(2).freq = Pfreq2;
    out(2).tim = Ptim2;



 close(27);

% figure(1); clf; 
% subplot(211);
% specgram(data, 4096, Fs, [], 4000); ylim([200 1000]); caxis([-10 30]); colormap('HOT');
% hold on;
% plot(out.tim, out.freq, 'g-', 'LineWidth', 1);
% subplot(212);
% plot(out.tim, out.amp, 'g-', 'LineWidth', 2);


