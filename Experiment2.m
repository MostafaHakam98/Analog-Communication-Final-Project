%%  MATLAB Analog Project
%   Experiment (2) Frequency Modulation

clc; clear all;

%%  Part (1.1):

%   Use Matlab to read the attached audio file, which has a sampling frequency Fs = 48 KHz.
%   Find the spectrum of this signal (the signal in frequency domain).

fs = 48 * 1000;
ts = 1/fs;

[voice_time, Fs] = audioread('eric.wav');
voice_time = transpose(voice_time);
n = length(voice_time);

t = linspace(0, n*ts, n);
f = linspace(-fs/2, fs/2, n);

figure(1);
subplot(2, 1, 1), plot(t, voice_time, 'r')
title('Audio signal in Time Domain'); xlabel('Time (sec)'); ylabel('Amplitude');

voice_frequency = fft(voice_time, n);
voice_frequency = fftshift(voice_frequency/n);               % Mesh fahmak yad

subplot(2,1,2),plot(f,abs(voice_frequency), 'b')
title('Audio spectrum in Frequency Domain'); xlabel('frequency(hz)'); ylabel('amplitude');

%%  Part (1.2):

%   Using an ideal Filter, remove all frequencies greater than 4 KHz.

figure(2);
low_pass_filter = rectpuls(f, 4000);
subplot(2, 1, 1), plot(f, low_pass_filter, 'g');
title('Ideal Low Pass Filter'); xlabel('Frequency (Hz)'); ylabel('Amplitude');

voice_frequency_filtered = low_pass_filter.*voice_frequency;
subplot(2, 1, 2), plot(f, abs(voice_frequency_filtered), 'b');
title('Audio spectrum in Frequency Domain after LPF'); xlabel('frequency(hz)'); ylabel('amplitude');

%%  Part (1.3):

%   Obtain the filtered signal in time domain, this is a band limited signal of BW=4 KHz.

voice_time_filtered = ifftshift(voice_frequency_filtered);
voice_time_filtered = ifft(voice_time_filtered)*length(voice_frequency_filtered);

figure(3);
subplot(2, 1, 1), plot(f, voice_time_filtered, 'r');
title('Audio signal in Time Domain (After Filter)'); xlabel('Time (sec)'); ylabel('Amplitude');
subplot(2, 1, 2), plot(f, abs(voice_frequency_filtered), 'b');
title('Audio spectrum in Frequency Domain (After Filter)'); xlabel('frequency(hz)'); ylabel('amplitude');

%%  Part (1.4):

%   Sound the filtered audio signal (make sure that there is only a small error in the filtered signal.

sound(real(voice_time_filtered), fs)
pause(8) 

%%  Part (2):

fc = 100 * 1000;
fs = 5 * fc;
ts = 1/fs;

resampled_message = resample(voice_time_filtered, fs, Fs);
cumulative_message = cumsum(resampled_message);

N = length(cumulative_message);
F = linspace(-fs/2, fs/2, N);
t = linspace(0, N*ts, N);

max_sample = max(cumulative_message);

Kf = 73.1;

NBFM_time = cos(2*pi*fc*t + Kf*cumulative_message);
NBFM_frequency = fftshift(fft(NBFM_time, N)/N);

figure(4)

subplot(2, 1, 1), plot(t, NBFM_time, 'r')
title('Audio signal in Time Domain'); xlabel('Time (sec)'); ylabel('Amplitude');

subplot(2, 1, 2), plot(F, NBFM_frequency, 'b');
title('Audio spectrum in Frequency Domain'); xlabel('frequency(hz)'); ylabel('amplitude');

%%  Part (3):

%   What is the condition we needed to achieve NBFM?

%   Beta <<< 1
%   Beta = delta(f)max / fm
%   Beta = Kf.m(t)/(2*pi*fm)
%   Kf.m(t)max/(2*pi*fm) <<< 1
%   Kf <<< 2*pi*fm
%   Kf <<< 8000*pi

%%  Part (4):

%   Demodulate the NBFM signal using a differentiator and an ED. 
%   For the differentiator, you can use the following command: diff. 
%   Assume no noise is introduced.  

NBFM_demodulated = diff(NBFM_time);
envelope_detector = abs(hilbert(NBFM_demodulated));

NBFM_time = 0.1*(envelope_detector - abs(mean (envelope_detector)));
resampled_message = resample(NBFM_time, fs, Fs);

t = linspace(0, length(resampled_message)/Fs, length(resampled_message));


figure(5)
plot(t, resampled_message); 
title('Demodulated Signal in Time Domain'); xlabel('Time (s)'); ylabel('amplitude');


sound(resampled_message, Fs);