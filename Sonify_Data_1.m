%% clear workspace
clc
clear all
close all
%% Data Import
% loading time series (can be anything: seismic, potential fields, climate record)
load('C:\Users\alexamin\Dropbox (UiO)\GoNorth\Processing\data\mat\geus\s101ch1')
t=D.t;% time vector
dt = t(2)-t(1);% time interval
Fs = 1/dt; % sampling frequency (in)
X = D.d(:,100); % data vector
%% Data Processing
% Bandpass filtering to remove noise
freq_min = 4; % low frequency for filtering 
freq_max = 40; % high frequency for filtering 
f_Nyquist = Fs/2; % Nyquist (max) frequency in the signal
[B,A] = butter(6,[freq_min/f_Nyquist freq_max/f_Nyquist]);% setup Butterworth filter
X = filtfilt(B,A,X-mean(X(~isnan(X))));% apply filter (make sure that its zero phase- no shifts of signal)
X=X/max(abs(X));% normalize signal to 1
% Compute the spectrogram using 
window = 128;  % Window length
noverlap = 120; % Overlap between windows
nfft = 256;    % Number of DFT points
% Compute the spectrogram using the short-time Fourier transform 
% signal is divided into segments "window" with "noverlap" points overlap, 
% and each segment is windowed with a Hamming window. 
% "nfft" is the number of frequency points used to calculate the discrete Fourier transform
% ff - output frequency vector
% tt - oupput time vector
[S, ff, tt] = spectrogram(X, window, noverlap, nfft, Fs);
S = medfilt2(abs(S),[7 7]); % get amplitude spectrum and smooth spectrogram using median filter
% Find the peak frequency in each time slice
[~, maxIndex] = max(abs(S), [], 1);
peakFrequencies = ff(maxIndex);% find peak frequency
FreqCoeff = 20; % scaling frequency coefficient to get audible signal  
%% Convert frequencies to MIDI note numbers
midiNotes = ceil(58 + 12 * log2(FreqCoeff*peakFrequencies / 440));
% initialize matrix:
N = length(midiNotes);  % number of notes
%loudness ('velocity') of signal normalized in range 0 to 127 (2^7)
loudness = 127*max(abs(S)/max(abs(S(:))));
%% Create MIDI structure and write into disk
M = zeros(N,6);
M(:,1) = 1;         % track 1
M(:,2) = 1;         % channel 1
M(:,3) = midiNotes(:);      % note numbers: one ocatave starting at middle C (60)
M(:,4) = loudness; 
M(:,5) = 0.1*(1:N);  % note on:  notes start every .3 seconds
M(:,6) = M(:,5) + .15;   % note off: each note has duration .3 seconds
% Optional - use other channels to play chords 
% M = [M; M; M];
% M(1:N,2)=2;
% M(1:N,3)=midiNotes(:)+4;
% M(N+1:2*N,2)=3;
% M(N+1:2*N,3)=midiNotes(:)+7;
midi_new = matrix2midi(M); % convert matrix to MIDI structure
writemidi(midi_new, 'AudioSeismogram-1.mid'); % write MIDI file onn disk
%% Plot and play sonified seismic trace
figure('WindowState','maximized')% create and maximize figure
y = midi2audio(midi_new,Fs*FreqCoeff,'sine');  % convert MIDI structure to audio to play in MATLAB
player = audioplayer(y, Fs*FreqCoeff); % create audio object
play(player); %start playing 
for it = 1:N
    subplot(311) % 3-rows and 1-column panel figure
    plot(t,X)% plot signal in the upper plane 
    hold on, 
    plot([tt(it), tt(it)],[min(X) max(X)],'k','LineWidth',2) % plot pointer for time 
    hold off
    ylabel('Normalize trace amplitude (n.d.)')
    subplot(312)
    imagesc(tt,ff,abs(S)), colormap(flipud(bone))% instaneous spectrum plot
    hold on,
    plot([tt(it), tt(it)],[min(ff) max(ff)],'k','LineWidth',2)
    hold off
    ylim([0 50])
    ylabel('Frequency (Hz)')
    subplot(313)
    plot(tt,loudness(:)), hold on,
    plot([tt(it), tt(it)],[0 150],'k','LineWidth',2)
    hold off
    xlabel('Trace time (s)'), ylabel('Sound amplitude (pnt)')
    drawnow
    pause(0.02)
end
%
stop(player)