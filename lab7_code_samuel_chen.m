%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MAE 107 LAB 7
% Samuel Chen
% 704-453-763
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% NOTE: This code does take a little bit of time to run (about 30 seconds
%%% before the song starts playing

%%%%% NOTE:
% I made my own song -- has seven instruments which are: voice, violin1,
% violin2, choir, taiko, horn, and piano
% inspired by sinusoidal waves

% opening ritual
close all; clc; clear all;

% import audio files
voice_mtx = audioread('voice.wav');
piano_mtx = audioread('piano.wav');
violin1_mtx = audioread('violin1.wav');
violin2_mtx = audioread('violin2.wav');
taiko_mtx = audioread('taiko.wav');
choir_mtx = audioread('choir.wav');
horn_mtx = audioread('horn.wav');

% define variables we will use later
sample_rate = 44100;  % typical WAV sampling
sample_time = 1/sample_rate;

% get the length of all of the matrices
voice_length = length(voice_mtx);
piano_length = length(piano_mtx);
violin1_length = length(violin1_mtx);
violin2_length = length(violin2_mtx);
taiko_length = length(taiko_mtx);
choir_length = length(choir_mtx);
horn_length = length(horn_mtx);
% we need the lengths to convert to time and also later to truncate
% to the nyquist frequency

% Get the time for all of the samples
voice_time = sample_time.*(1:voice_length);
piano_time = sample_time.*(1:piano_length);
violin1_time = sample_time.*(1:violin1_length);
violin2_time = sample_time.*(1:violin2_length);
taiko_time = sample_time.*(1:taiko_length);
choir_time = sample_time.*(1:choir_length);
horn_time = sample_time.*(1:horn_length);


%%% INPUT PLOTS %%%

figure(1)
plot(voice_time, voice_mtx);
xlim([0 5]);
grid on;
xlabel('Time (seconds)'); ylabel('Amplitude (Input)');
title('Voice Amplitude versus Time');

figure(2)
plot(piano_time, piano_mtx);
xlim([0 5]);
grid on;
xlabel('Time (seconds)'); ylabel('Amplitude (Input)');
title('Piano Amplitude versus Time');

figure(3)
plot(violin1_time, violin1_mtx);
xlim([0 5]);
grid on;
xlabel('Time (seconds)'); ylabel('Amplitude (Input)');
title('Violin 1 Amplitude versus Time');

figure(4)
plot(violin2_time, violin2_mtx);
xlim([0 5]);
grid on;
xlabel('Time (seconds)'); ylabel('Amplitude (Input)');
title('Violin 2 Amplitude versus Time');

figure(5)
plot(taiko_time, taiko_mtx);
xlim([15 20]);   % note here that in the song, the taiko doesn't start
                 % until 15 seconds in so the first 15 seconds starts here
grid on;
xlabel('Time (seconds)'); ylabel('Amplitude (Input)');
title('Taiko Amplitude versus Time');

figure(6)
plot(choir_time, choir_mtx);
xlim([15 20]);  % same thing with choir, doesn't start until 15 seconds
grid on;
xlabel('Time (seconds)'); ylabel('Amplitude (Input)');
title('Choir Amplitude versus Time');

figure(7)
plot(horn_time, horn_mtx);
xlim([17 22]);  % horn starts later at 17 seconds (about half the song);
grid on;
xlabel('Time (seconds)'); ylabel('Amplitude (Input)');
title('French Horn Amplitude versus Time');


%%%% SECOND PART - FOURIER TRANSFORM

nyquist = 22050;  % needed later

% take the step size
voice_step = sample_rate./voice_length;
piano_step = sample_rate./piano_length;
violin1_step = sample_rate./violin1_length;
violin2_step = sample_rate./violin2_length;
taiko_step = sample_rate./taiko_length;
choir_step = sample_rate./choir_length;
horn_step = sample_rate./horn_length;
% these are the step sizes after convering to frequency
% note through dimensional analysis (samples/sec) divided by (samples)
% which equals 1/sec which equals Hz

% take the fourier transform
voice_fft = (fft(voice_mtx)).*sample_time;
piano_fft = (fft(piano_mtx)).*sample_time;
violin1_fft = (fft(violin1_mtx)).*sample_time;
violin2_fft = (fft(violin2_mtx)).*sample_time;
taiko_fft = (fft(taiko_mtx)).*sample_time;
choir_fft = (fft(choir_mtx)).*sample_time;
horn_fft = (fft(horn_mtx)).*sample_time;
% fourier transform times the sample time


% truncate to the nyquist frequency
voice_length_new = 0.5*voice_length;
piano_length_new = 0.5*piano_length;
violin1_length_new = 0.5*violin1_length;
violin2_length_new = 0.5*violin2_length;
taiko_length_new = 0.5*taiko_length;
choir_length_new = 0.5*choir_length;
horn_length_new = 0.5*horn_length;
% corresponds to 22050
% this means rate of 44100 in the half, so through dimensional analysis
% and logic we can just take half of our samples and that should truncate
% right or very close to the Nyquist frequency


% truncate the fourier transform to the nyquist frequency
voice_fft = voice_fft(1:voice_length_new);
piano_fft = piano_fft(1:piano_length_new);
violin1_fft = violin1_fft(1:floor(violin1_length_new));  % floor for fractions
violin2_fft = violin2_fft(1:floor(violin2_length_new));
taiko_fft = taiko_fft(1:floor(taiko_length_new));
choir_fft = choir_fft(1:choir_length_new);
horn_fft = horn_fft(1:horn_length_new);

% truncate omega to the nyquist frequency and establish it
% note that the matrix dimensions have to agree
voice_om = voice_step.*(1:voice_length_new);
piano_om = piano_step.*(1:piano_length_new);
violin1_om = violin1_step.*(1:floor(violin1_length_new));  % add floor or else
                                                           % we have
                                                           % fractions
violin2_om = violin2_step.*(1:floor(violin2_length_new));
taiko_om = taiko_step.*(1:floor(taiko_length_new));
choir_om = choir_step.*(1:choir_length_new);
horn_om = horn_step.*(1:horn_length_new);


%%% FOURIER TRANSFORM PLOTS
figure(8)
loglog(voice_om, abs(voice_fft), '.');
grid on;
xlabel('Frequency (Hz)'); ylabel ('Spectrum Magnitude');
title('Voice Magnitude Spectra vs Frequency (Hz)');

figure(9)
loglog(piano_om, abs(piano_fft), '.');
grid on;
xlabel('Frequency (Hz)'); ylabel ('Spectrum Magnitude');
title('Piano Magnitude Spectra vs Frequency (Hz)');

figure(10)
loglog(violin1_om, abs(violin1_fft), '.');
grid on;
xlabel('Frequency (Hz)'); ylabel ('Spectrum Magnitude');
title('Violin 1 Magnitude Spectra vs Frequency (Hz)');

figure(11)
loglog(violin2_om, abs(violin2_fft), '.');
grid on;
xlabel('Frequency (Hz)'); ylabel ('Spectrum Magnitude');
title('Violin 2 Magnitude Spectra vs Frequency (Hz)');

figure(12)
loglog(taiko_om, abs(taiko_fft), '.');
grid on;
xlabel('Frequency (Hz)'); ylabel ('Spectrum Magnitude');
title('Taiko Magnitude Spectra vs Frequency (Hz)');

figure(13)
loglog(choir_om, abs(choir_fft), '.');
grid on;
xlabel('Frequency (Hz)'); ylabel ('Spectrum Magnitude');
title('Choir Magnitude Spectra vs Frequency (Hz)');

figure(14)
loglog(horn_om, abs(horn_fft), '.');
grid on;
xlabel('Frequency (Hz)'); ylabel ('Spectrum Magnitude');
title('French Horn Magnitude Spectra vs Frequency (Hz)');



%%% FILTERS
% start the filters based on the spectrum magnitude
[voice_num, voice_den] = butter(4, [80 10000]./nyquist, 'bandpass');
voice_filtered = filter(voice_num, voice_den, voice_mtx);

[piano_num, piano_den] = butter(4, [100 10000]./nyquist, 'bandpass');
piano_filtered = filter(piano_num, piano_den, piano_mtx);

[violin1_num, violin1_den] = butter(4, [200 10000]./nyquist, 'bandpass');
violin1_filtered = filter(violin1_num, violin1_den, violin1_mtx);

[violin2_num, violin2_den] = butter(4, 6000/nyquist, 'low');
violin2_filtered = filter(violin2_num, violin2_den, violin2_mtx);

[taiko_num, taiko_den] = butter(4, [40 10000]./nyquist, 'bandpass');
taiko_filtered = filter(taiko_num, taiko_den, taiko_mtx);

[choir_num, choir_den] = butter(4, 2000/nyquist, 'low');
choir_filtered = filter(choir_num, choir_den, choir_mtx);

[horn_num, horn_den] = butter(4, 3000/nyquist, 'low');
horn_filtered = filter(horn_num, horn_den, horn_mtx);

% lengths were not the same for some reason
violin1_filtered = violin1_filtered(1:1546006,:);
violin2_filtered = violin2_filtered(1:1546006,:);
taiko_filtered = taiko_filtered(1:1546006,:);
taiko_mtx = taiko_mtx(1:1546006,:);

% I like the original sound for the horn so I didn't use the
% filtered version

my_comp = (1.1)*voice_filtered + (0.9)*piano_filtered + violin1_filtered + ...
    1.05*violin2_filtered + (0.9)*taiko_filtered + choir_filtered + horn_mtx;

%audiowrite('samuel_chen_composition.wav', 0.85*my_comp, 44100);


sound(0.85*my_comp, 44100);

% filtered taiko graph for show
taiko_filtered_ftt = (fft(taiko_filtered))*sample_time;
taiko_filtered_ftt = taiko_filtered_ftt(1:length(taiko_om));

figure(15)
loglog(taiko_om, abs(taiko_filtered_ftt), '.');
grid on;
xlabel('Frequency (Hz)'); ylabel ('Spectrum Magnitude');
title('Filtered Taiko Magnitude Spectra vs Frequency (Hz)');


