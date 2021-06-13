%%  MATLAB Analog Project
%   Experiment (1) Double SideBand Modulation

clc; clear all;

%%  Part (1):

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

%%  Part (2):

%   Using an ideal Filter, remove all frequencies greater than 4 KHz.

figure(2);
low_pass_filter = rectpuls(f, 4000);
subplot(2, 1, 1), plot(f, low_pass_filter, 'g');
title('Ideal Low Pass Filter'); xlabel('Frequency (Hz)'); ylabel('Amplitude');

voice_frequency_filtered = low_pass_filter.*voice_frequency;
subplot(2, 1, 2), plot(f, abs(voice_frequency_filtered), 'b');

%%  Part (3):

%   Obtain the filtered signal in time domain, this is a band limited signal of BW=4 KHz.

voice_time_filtered = ifftshift(voice_frequency_filtered);
voice_time_filtered = ifft(voice_time_filtered)*length(voice_frequency_filtered);

figure(3);
subplot(2, 1, 1), plot(f, voice_time_filtered, 'r');
title('Audio signal in Time Domain (After Filter)'); xlabel('Time (sec)'); ylabel('Amplitude');
subplot(2, 1, 2), plot(f, abs(voice_frequency_filtered), 'b');
title('Audio spectrum in Frequency Domain (After Filter)'); xlabel('frequency(hz)'); ylabel('amplitude');

%%  Part (4):

%   Sound the filtered audio signal (make sure that there is only a small error in the filtered signal.

sound(real(voice_time_filtered), fs)
pause(8) 

%%  Part (5):

%   Modulate the carrier with the filtered signal you obtained
%   You are required to generate both types of modulation (DSB-TC and DSB-SC). 
%   Choose a carrier frequency of 100 KHz. 
%   For the DSB-TC take the DC bias added to message before modulation to be twice the maximum of the message (modulation index =0.5 in this case).
%   You have to sketch the modulated signal of both DSB-TC & DSB-SC in frequency domain.

fc = 100 * 1000;
fs = 5 * fc;

carrier = cos(2*pi*fc*t);

voice_time_SC = voice_time_filtered.*carrier;
voice_frequency_SC = fft(voice_time_SC, n);
voice_frequency_SC = fftshift(voice_frequency_SC/n);

voice_time_TC = (1 + 0.5*voice_time_filtered).*carrier;
voice_frequency_TC = fft(voice_time_TC, n);
voice_frequency_TC = fftshift(voice_frequency_TC/n); 

figure(4)

subplot(4,1,1),plot(t, voice_time_SC, 'r');
title('DSB-SC Modulation in time domain'); xlabel('Time (sec)'); ylabel('Amplitude');

subplot(4,1,2),plot(f, abs(voice_frequency_SC), 'b');
title('DSB-SC Modulation in frequency domain'); xlabel('Frequency (HZ)'); ylabel('Amplitude');

subplot(4,1,3),plot(t, voice_time_TC, 'r');
title('DSB-TC Modulation in time domain'); xlabel('Time (sec)'); ylabel('Amplitude');

subplot(4,1,4),plot(f, abs(voice_frequency_TC), 'b');
title('DSB-TC Modulation in frequency domain'); xlabel('Frequency (HZ)'); ylabel('Amplitude');

%%  Part (6):

%   For both types of modulations (DSB-SC & DSB-TC), use envelop detector to receive the message (assume no noise).

envelope_SC = abs(hilbert(voice_time_SC));
envelope_TC = abs(hilbert(voice_time_TC));

%%  Part (7):

%   Sketch the received signal in time domain, and Play the received signal back.

figure(5)
subplot(2,1,1),plot(t, envelope_SC, 'r');
title('Envelope Detection of DSB-SC'); xlabel('Time (sec)'); ylabel('Amplitude');

subplot(2,1,2),plot(t, envelope_TC, 'r');
title('Envelope Detection of DSB-TC'); xlabel('Time (sec)'); ylabel('Amplitude');

sound(real(envelope_SC), Fs)
pause(8) 
sound(real(envelope_TC), Fs)
pause(8) 

%   Observation: DSB-TC is better used with the envelope detector.

%% Part (8):

%   For DSB-TC only add noise to modulated signal with SNR = 0, 10, and 30 db 
%   Receive them with envelope detector. 
%   Play back the sound file each time after detection and sketch it in time domain.

SNR = [0 10 30];
voice_noise_list = [];

for snr = SNR
    voice_noise = awgn(voice_time_TC, snr, 'measured');
    voice_noise = abs(hilbert(voice_noise));
    voice_noise_list = [voice_noise_list voice_noise];
end

figure(6);

subplot(3,1,1),plot(t, voice_noise_list(1 : n), 'r');
title('DSB-TC Signal at SNR = 0 dB'); xlabel('Time (sec)'); ylabel('Amplitude');

subplot(3,1,2),plot(t, voice_noise_list(1 + n : 2*n), 'r');
title('DSB-TC Signal at SNR = 10 dB'); xlabel('Time (sec)'); ylabel('Amplitude');

subplot(3,1,3),plot(t, voice_noise_list(1 + 2*n : 3*n), 'r');
title('DSB-TC Signal at SNR = 30 dB'); xlabel('Time (sec)'); ylabel('Amplitude');

sound(real(voice_noise_list(1 : n)), Fs)
pause(8)
sound(real(voice_noise_list(1 + n : 2*n)), Fs)
pause(8)
sound(real(voice_noise_list(1 + 2*n: 3*n)), Fs)
pause(8)

%%  Part (9):

%   Use coherent detection to receive the modulated signal with SNR=0, 10, 30 dB 
%   Sound the received signals and plot them in both time and frequency domain.

SNR = [0 10 30];
voice_noise_list = [];

for snr = SNR
    voice_noise = awgn(voice_time_SC, snr, 'measured');
    voice_noise_list = [voice_noise_list voice_noise];
end

figure(7)

voice_noise =  voice_noise_list(1 : n);
voice_noise = voice_noise.* carrier;

voice_noise_frequency = fft(voice_noise, n);
voice_noise_frequency = fftshift(voice_noise_frequency/n);
voice_noise_frequency = low_pass_filter.*voice_noise_frequency;
subplot(2,1,2),plot(f, abs(voice_noise_frequency), 'b');
title('DSB-SC Spectrum at SNR = 0 dB'); xlabel('Frequency (Hz)'); ylabel('Amplitude');

voice_noise = ifftshift(voice_noise_frequency);
voice_noise = ifft(voice_noise)*length(voice_noise);
subplot(2,1,1),plot(t, voice_noise, 'r');
title('DSB-SC Signal at SNR = 0 dB'); xlabel('Time (sec)'); ylabel('Amplitude');

sound(real(voice_noise), Fs)
pause(8)

figure(8)

voice_noise = voice_noise_list(1 + n : 2*n);
voice_noise = voice_noise.* carrier;

voice_noise_frequency = fft(voice_noise, n);
voice_noise_frequency = fftshift(voice_noise_frequency/n);
voice_noise_frequency = low_pass_filter.*voice_noise_frequency;
subplot(2,1,2),plot(f, abs(voice_noise_frequency), 'b');
title('DSB-SC Spectrum at SNR = 10 dB'); xlabel('Frequency (Hz)'); ylabel('Amplitude');

voice_noise = ifftshift(voice_noise_frequency);
voice_noise = ifft(voice_noise)*length(voice_noise);
subplot(2,1,1),plot(t, voice_noise, 'r');
title('DSB-SC Signal at SNR = 10 dB'); xlabel('Time (sec)'); ylabel('Amplitude');

sound(real(voice_noise), Fs)
pause(8)

figure(9)

voice_noise = voice_noise_list(1 + 2*n : 3*n);
voice_noise = voice_noise.* carrier;

voice_noise_frequency = fft(voice_noise, n);
voice_noise_frequency = fftshift(voice_noise_frequency/n);
voice_noise_frequency = low_pass_filter.*voice_noise_frequency;
subplot(2,1,2),plot(f, abs(voice_noise_frequency), 'b');
title('DSB-SC Spectrum at SNR = 30 dB'); xlabel('Frequency (Hz)'); ylabel('Amplitude');

voice_noise = ifftshift(voice_noise_frequency);
voice_noise = ifft(voice_noise)*length(voice_noise);
subplot(2,1,1),plot(t, voice_noise, 'r');
title('DSB-SC Signal at SNR = 30 dB'); xlabel('Time (sec)'); ylabel('Amplitude');

sound(real(voice_noise), Fs)
pause(8)

%%  Part(10)

%   Repeat the coherent detection with frequency error, F=100.1 KHz instead of 100 KHz and Find the error.

SNR = [0 10 30];
voice_noise_list = [];

for snr = SNR
    voice_noise = awgn(voice_time_SC, snr, 'measured');
    voice_noise_list = [voice_noise_list voice_noise];
end

figure(10)

voice_noise =  voice_noise_list(1 : n);
voice_noise = voice_noise.* cos(2*pi*100.1*1000*t);

voice_noise_frequency = fft(voice_noise, n);
voice_noise_frequency = fftshift(voice_noise_frequency/n);
voice_noise_frequency = low_pass_filter.*voice_noise_frequency;
subplot(2,1,2),plot(f, abs(voice_noise_frequency), 'b');
title('DSB-SC Spectrum at SNR = 0 dB with fc = 100.1 KHz'); xlabel('Frequency (Hz)'); ylabel('Amplitude');

voice_noise = ifftshift(voice_noise_frequency);
voice_noise = ifft(voice_noise)*length(voice_noise);
subplot(2,1,1),plot(t, voice_noise, 'r');
title('DSB-SC Signal at SNR = 0 dB with fc = 100.1 KHz'); xlabel('Time (sec)'); ylabel('Amplitude');

sound(real(voice_noise), Fs)
pause(8)

figure(11)

voice_noise = voice_noise_list(1 + n : 2*n);
voice_noise = voice_noise.* cos(2*pi*100.1*1000*t);

voice_noise_frequency = fft(voice_noise, n);
voice_noise_frequency = fftshift(voice_noise_frequency/n);
voice_noise_frequency = low_pass_filter.*voice_noise_frequency;
subplot(2,1,2),plot(f, abs(voice_noise_frequency), 'b');
title('DSB-SC Spectrum at SNR = 10 dB with fc = 100.1 KHz'); xlabel('Frequency (Hz)'); ylabel('Amplitude');

voice_noise = ifftshift(voice_noise_frequency);
voice_noise = ifft(voice_noise)*length(voice_noise);
subplot(2,1,1),plot(t, voice_noise, 'r');
title('DSB-SC Signal at SNR = 10 dB with fc = 100.1 KHz'); xlabel('Time (sec)'); ylabel('Amplitude');

sound(real(voice_noise), Fs)
pause(8)

figure(12)

voice_noise = voice_noise_list(1 + 2*n : 3*n);
voice_noise = voice_noise.* cos(2*pi*100.1*1000*t);

voice_noise_frequency = fft(voice_noise, n);
voice_noise_frequency = fftshift(voice_noise_frequency/n);
voice_noise_frequency = low_pass_filter.*voice_noise_frequency;
subplot(2,1,2),plot(f, abs(voice_noise_frequency), 'b');
title('DSB-SC Spectrum at SNR = 30 dB with fc = 100.1 KHz'); xlabel('Frequency (Hz)'); ylabel('Amplitude');

voice_noise = ifftshift(voice_noise_frequency);
voice_noise = ifft(voice_noise)*length(voice_noise);
subplot(2,1,1),plot(t, voice_noise, 'r');
title('DSB-SC Signal at SNR = 30 dB with fc = 100.1 KHz'); xlabel('Time (sec)'); ylabel('Amplitude');

sound(real(voice_noise), Fs)
pause(8)

%%  Part(10)

%%  Repeat the coherent detection with phase error = 20

SNR = [0 10 30];
voice_noise_list = [];

for snr = SNR
    voice_noise = awgn(voice_time_SC, snr, 'measured');
    voice_noise_list = [voice_noise_list voice_noise];
end

figure(13)

voice_noise =  voice_noise_list(1 : n);
voice_noise = voice_noise.* cos(2*pi*fc*t + (20*pi/180));

voice_noise_frequency = fft(voice_noise, n);
voice_noise_frequency = fftshift(voice_noise_frequency/n);
voice_noise_frequency = low_pass_filter.*voice_noise_frequency;
subplot(2,1,2),plot(f, abs(voice_noise_frequency), 'b');
title('DSB-SC Spectrum at SNR = 0 dB with fc = 100.1 KHz'); xlabel('Frequency (Hz)'); ylabel('Amplitude');

voice_noise = ifftshift(voice_noise_frequency);
voice_noise = ifft(voice_noise)*length(voice_noise);
subplot(2,1,1),plot(t, voice_noise, 'r');
title('DSB-SC Signal at SNR = 0 dB with fc = 100.1 KHz'); xlabel('Time (sec)'); ylabel('Amplitude');

sound(real(voice_noise), Fs)
pause(8)

figure(14)

voice_noise = voice_noise_list(1 + n : 2*n);
voice_noise = voice_noise.* cos(2*pi*fc*t + (20*pi/180));

voice_noise_frequency = fft(voice_noise, n);
voice_noise_frequency = fftshift(voice_noise_frequency/n);
voice_noise_frequency = low_pass_filter.*voice_noise_frequency;
subplot(2,1,2),plot(f, abs(voice_noise_frequency), 'b');
title('DSB-SC Spectrum at SNR = 10 dB with fc = 100.1 KHz'); xlabel('Frequency (Hz)'); ylabel('Amplitude');

voice_noise = ifftshift(voice_noise_frequency);
voice_noise = ifft(voice_noise)*length(voice_noise);
subplot(2,1,1),plot(t, voice_noise, 'r');
title('DSB-SC Signal at SNR = 10 dB with fc = 100.1 KHz'); xlabel('Time (sec)'); ylabel('Amplitude');

sound(real(voice_noise), Fs)
pause(8)

figure(15)

voice_noise = voice_noise_list(1 + 2*n : 3*n);
voice_noise = voice_noise.* cos(2*pi*fc*t + (20*pi/180));

voice_noise_frequency = fft(voice_noise, n);
voice_noise_frequency = fftshift(voice_noise_frequency/n);
voice_noise_frequency = low_pass_filter.*voice_noise_frequency;
subplot(2,1,2),plot(f, abs(voice_noise_frequency), 'b');
title('DSB-SC Spectrum at SNR = 30 dB with fc = 100.1 KHz'); xlabel('Frequency (Hz)'); ylabel('Amplitude');

voice_noise = ifftshift(voice_noise_frequency);
voice_noise = ifft(voice_noise)*length(voice_noise);
subplot(2,1,1),plot(t, voice_noise, 'r');
title('DSB-SC Signal at SNR = 30 dB with fc = 100.1 KHz'); xlabel('Time (sec)'); ylabel('Amplitude');

sound(real(voice_noise), Fs)
pause(8)