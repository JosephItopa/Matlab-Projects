close all;
clc;
clear all;

N  = 1000; % Number of Samples
N1  = 1024;
N2  = 2018;
Fs0 = 150; % Sampling Frequency in Hz % Sampling Time(Ts) = 1/Fs = 1/150 = 0.006 sec
Fs1 = 1000;
Fs2 = 200; 
f1 = 5; % Signal Frequency 1
f2 = 21; % Signal Frequency 2
A1 = 1; % Signal Amplitude 1
A2 = 0.4; % Signal Amplitude 2

% 
t = ((-N/2):(N/2)-1)/Fs0; % Time axis from -2 sec to 2 sec or, t = (-500:499)/250 = (-500:1:499)/250 = -2:0.004:1.996
y = A1.*sin(2*pi*f1*t) + A2.*sin(2*pi*f2*t); % Combination of Sine Wave with Multiple frequencies
% sine wave using sampling frequency(Fs1) = 1000, sampling frequency(Fs2) = 200, number of samples(N1) = 1024, and number of samples(N1) = 2018 
t1_a = ((-N1/2):(N1/2)-1)/Fs1; % number of samples = 1024, Fs1 = 1000
t1_b = ((-N2/2):(N2/2)-1)/Fs1; % number of samples = 2018, Fs1 = 1000
y1_a = A1.*sin(2*pi*f1*t1_a) + A2.*sin(2*pi*f2*t1_a); %sine wave for Fs1 = 1000, and N1 = 1024
y1_b = A1.*sin(2*pi*f1*t1_b) + A2.*sin(2*pi*f2*t1_b); %sine wave for Fs1 = 1000, and N1 = 2018

t2_a = ((-N1/2):(N1/2)-1)/Fs2; % number of samples = 1024, Fs2 = 200
t2_b = ((-N2/2):(N2/2)-1)/Fs2; % number of samples = 2018, Fs2 = 200
y2_a = A1.*sin(2*pi*f1*t2_a) + A2.*sin(2*pi*f2*t2_a); %sine wave for Fs2 = 200, and N1 = 1024
y2_b = A1.*sin(2*pi*f1*t2_b) + A2.*sin(2*pi*f2*t2_b); %sine wave for Fs2 = 200, and N1 = 2018
%% TASK 4.1(a) - FIND THE POWER SPECTRAL DENSITY
% Plot of sine wave

figure('Name','Plot of sine waves');

subplot(3,2,1)
plot(t,y), title('Sine wave'),ylim([-1.5 1.5]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on
% Convolution of hamming window and sine signal 
hammingWindow = hamming(length(t))'; % Hamming Window

sig0 = y.*hammingWindow; % Windowed Signal (Hamming Window)
subplot(3,2,2)
plot(t,sig0),title('Windowed Signal'),ylim([-1.5 1.5]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on

% Autocorrelation function with Hamming window

[r0,lags0] = xcorr(sig0,'biased');
tau0 = lags0/Fs0;

subplot(3,2,3)
plot(tau0,r0),title('ACF with Hamming Window'), xlabel('Time difference \tau (in sec)'), ylabel('Amplitude')
grid on

% Power Spectral Density using Wiener Khintchine Theorem with Hamming window

Rxxdft0 = abs(fftshift(fft(r0)))/N;
freq0 = -Fs0/2:Fs0/length(r0):Fs0/2-(Fs0/length(r0));

subplot(3,2,4)
plot(freq0, Rxxdft0),title({'Power Spectral Density using Wiener Khintchine Theorem ', 'with Hamming window'}),xlim([-50 50]), xlabel('Frequency f (in Hz)'),ylabel('Spectral Power')
grid on


%% FOR SAMPLING FREQUENCY = 1000; NUMBER OF SAMPLES = 1024
% Plot of sine wave

subplot(3,2,1)
plot(t1_a, y1_a),title('Sine wave for sine wave for 1000Hz sampling frequency, and 1024 samples'),ylim([-1.5 1.5]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on
% Convolution of hamming window and sine signal 
hammingWindow0 = hamming(length(t1_a))'; % Hamming Window

sig1 = y1_a.*hammingWindow0; % Windowed Signal (Hamming Window)
subplot(3,2,3)
plot(t1_a, sig1),title('Windowed Signal for 1000Hz sampling frequency, and 1024 samples'),ylim([-1.5 1.5]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on

% Autocorrelation function with Hamming window

[r1,lags1] = xcorr(sig1,'biased');
tau1 = lags1/Fs1;

subplot(3,2,2)
plot(tau1,r1),title('ACF with Hamming Window for 1000Hz sampling frequency and 1024 samples'), xlabel('Time difference \tau (in sec)'), ylabel('Amplitude')
grid on

% Power Spectral Density using Wiener Khintchine Theorem with Hamming window

Rxxdft1 = abs(fftshift(fft(r1)))/N;
freq1 = -Fs1/2:Fs1/length(r1):Fs1/2-(Fs1/length(r1));

subplot(3,2,4)
plot(freq1, Rxxdft1),title({'Power Spectral Density using Wiener Khintchine Theorem for 1000Hz sampling frequency & 1024 samples'}),xlim([-50 50]), xlabel('Frequency f (in Hz)'),ylabel('Spectral Power')
grid on

%% FOR SAMPLING FREQUENCY = 1000; NUMBER OF SAMPLES = 2018
% Plot of sine wave

subplot(3,2,1)
plot(t1_b, y1_b),title('Sine wave for sine wave for 1000Hz sampling frequency, and 2018 samples'),ylim([-1.5 1.5]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on
% Convolution of hamming window and sine signal 
hammingWindow0 = hamming(length(t1_b))'; % Hamming Window

sig1 = y1_b.*hammingWindow0; % Windowed Signal (Hamming Window)
subplot(3,2,3)
plot(t1_b, sig1),title('Windowed Signal for 1000Hz sampling frequency, and 2018 samples'),ylim([-1.5 1.5]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on

% Autocorrelation function with Hamming window

[r1,lags1] = xcorr(sig1,'biased');
tau1 = lags1/Fs1;

subplot(3,2,2)
plot(tau1,r1),title('ACF with Hamming Window for 1000Hz sampling frequency and 2018 samples'), xlabel('Time difference \tau (in sec)'), ylabel('Amplitude')
grid on

% Power Spectral Density using Wiener Khintchine Theorem with Hamming window

Rxxdft1 = abs(fftshift(fft(r1)))/N;
freq1 = -Fs1/2:Fs1/length(r1):Fs1/2-(Fs1/length(r1));

subplot(3,2,4)
plot(freq1, Rxxdft1),title({'Power Spectral Density using Wiener Khintchine Theorem for 1000Hz sampling frequency & 2018 samples'}),xlim([-50 50]), xlabel('Frequency f (in Hz)'),ylabel('Spectral Power')
grid on

%% FOR SAMPLING FREQUENCY = 200; NUMBER OF SAMPLES = 1024

% Plot of sine wave

subplot(3,2,1)
plot(t2_a, y2_a),title('Sine wave for sine wave for 200Hz sampling frequency, and 1024 samples'),ylim([-1.5 1.5]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on
% Convolution of hamming window and sine signal 
hammingWindow0 = hamming(length(t2_a))'; % Hamming Window

sig1 = y2_a.*hammingWindow0; % Windowed Signal (Hamming Window)
subplot(3,2,3)
plot(t2_a, sig1),title('Windowed Signal for 200Hz sampling frequency, and 1024 samples'),ylim([-1.5 1.5]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on

% Autocorrelation function with Hamming window

[r1,lags1] = xcorr(sig1,'biased');
tau1 = lags1/Fs2;

subplot(3,2,2)
plot(tau1,r1),title('ACF with Hamming Window for 200Hz sampling frequency and 1024 samples'), xlabel('Time difference \tau (in sec)'), ylabel('Amplitude')
grid on

% Power Spectral Density using Wiener Khintchine Theorem with Hamming window

Rxxdft1 = abs(fftshift(fft(r1)))/N;
freq1 = -Fs2/2:Fs2/length(r1):Fs2/2-(Fs2/length(r1));

subplot(3,2,4)
plot(freq1, Rxxdft1),title({'Power Spectral Density using Wiener Khintchine Theorem for 200Hz sampling frequency & 1024 samples'}),xlim([-50 50]), xlabel('Frequency f (in Hz)'),ylabel('Spectral Power')
grid on

%%FOR SAMPLING FREQUENCY = 200; NUMBER OF SAMPLES = 2018

% Plot of sine wave

subplot(3,2,1)
plot(t2_b, y2_b),title('Sine wave for sine wave for 200Hz sampling frequency, and 2018 samples'),ylim([-1.5 1.5]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on
% Convolution of hamming window and sine signal 
hammingWindow0 = hamming(length(t2_b))'; % Hamming Window

sig1 = y2_b.*hammingWindow0; % Windowed Signal (Hamming Window)
subplot(3,2,3)
plot(t2_b, sig1),title('Windowed Signal for 200Hz sampling frequency, and 2018 samples'),ylim([-1.5 1.5]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on

% Autocorrelation function with Hamming window

[r1,lags1] = xcorr(sig1,'biased');
tau1 = lags1/Fs2;

subplot(3,2,2)
plot(tau1,r1),title('ACF with Hamming Window for 200Hz sampling frequency and 2018 samples'), xlabel('Time difference \tau (in sec)'), ylabel('Amplitude')
grid on

% Power Spectral Density using Wiener Khintchine Theorem with Hamming window

Rxxdft1 = abs(fftshift(fft(r1)))/N;
freq1 = -Fs2/2:Fs2/length(r1):Fs2/2-(Fs2/length(r1));

subplot(3,2,4)
plot(freq1, Rxxdft1),title({'Power Spectral Density using Wiener Khintchine Theorem for 200Hz sampling frequency & 2018 samples'}),xlim([-50 50]), xlabel('Frequency f (in Hz)'),ylabel('Spectral Power')
grid on

%% TASK 4.1(b) FOR SAMPLING FREQUENCY = 30

Fs = 30;
f1 = 5; % Signal Frequency 1
f2 = 21; % Signal Frequency 2
A1 = 1; % Signal Amplitude 1
A2 = 0.4; % Signal Amplitude 2

t = ((-N/2):(N/2)-1)/Fs; % Time axis 
y = A1.*sin(2*pi*f1*t) + A2.*sin(2*pi*f2*t); % Combination of Sine Wave with Multiple frequencies

% Plot of sine wave

subplot(3,2,1)
plot(t, y),title('Sine wave for sine wave for 30Hz sampling frequency'),ylim([-1.5 1.5]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on
% Convolution of hamming window and sine signal 
hammingWindow0 = hamming(length(t))'; % Hamming Window

sig1 = y.*hammingWindow0; % Windowed Signal (Hamming Window)
subplot(3,2,3)
plot(t, sig1),title('Windowed Signal for 30Hz sampling frequency'),ylim([-1.5 1.5]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on
% Autocorrelation function with Hamming window

[r1,lags1] = xcorr(sig1,'biased');
tau1 = lags1/Fs;

subplot(3,2,2)
plot(tau1,r1),title('ACF with Hamming Window for 30Hz sampling frequency'), xlabel('Time difference \tau (in sec)'), ylabel('Amplitude')
grid on

%Power Spectral Density using Wiener Khintchine Theorem with Hamming window

Rxxdft1 = abs(fftshift(fft(r1)))/N;
freq1 = -Fs/2:Fs/length(r1):Fs/2-(Fs/length(r1));

subplot(3,2,4)
plot(freq1, Rxxdft1),title({'Power Spectral Density using Wiener Khintchine Theorem for 30Hz sampling frequency'}),xlim([-50 50]), xlabel('Frequency f (in Hz)'),ylabel('Spectral Power')
grid on

%% TASK 4.2(a, b)

close all;
clc;
clear all;

alpha = 3.0;
N  = 1024; % Number of Samples
Fs = 100; % Sampling Frequency in Hz % Sampling Time(Ts) = 1/Fs = 1/250 = 0.004 sec
f = 20; % Signal Frequency
t = ((-N/2):(N/2)-1)/Fs; % Time axis 
x = sin(2*pi*f*t); % Sine Wave

noise = alpha.*randn(1,length(t)); % Noise
% plot of original signal
subplot(2,2,1)
plot(t, x),title('Original Signal with no noise'), ylim([-3 3]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on 


% plot of sine signal with noise
y = x + noise; % Noise added to signal
subplot(2,2,2)
plot(t, y),title('Original Signal with noise'), ylim([-3 3]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on 

hammingWindow = hamming(length(t))'; % Hamming Window

figure('Name','Signal evaluation');
subplot(2,2,2)
plot(t, y),title('Signal with Noise'), ylim([-3 3]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on
subplot(2,2,3)
plot(t, hammingWindow), title('Hamming Window'), ylim([-1.5 1.5]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on
% Hamming Window ACF
sig2 = y.*hammingWindow;

[r2, lags2] = xcorr(sig2,'biased');
tau2 = lags2/Fs;

subplot(2,2,4)
plot(tau2, r2),title('ACF using Hamming Window'), xlabel('Time difference \tau (in sec)'), ylabel('Amplitude')
grid on
%% Task 4.2(c, d, e)
%% Hamming window
close all;
clc;
clear all;

alpha = 3.0;
N  = 2048; % Number of Samples
Fs = 100; % Sampling Frequency in Hz % Sampling Time(Ts) = 1/Fs = 1/250 = 0.004 sec
f = 20; % Signal Frequency
t = ((-N/2):(N/2)-1)/Fs; % Time axis 
x = sin(2*pi*f*t); % Sine Wave

noise = alpha.*randn(1,length(t)); % Noise
%%
% plot of original signal
subplot(2,2,1)
plot(t, x),title('Original Signal with no noise'), ylim([-3 3]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on 


% plot of sine signal with noise
y = x + noise; % Noise added to signal
subplot(2,2,2)
plot(t, y),title('Original Signal with noise'), ylim([-3 3]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on 

hammingWindow = hamming(length(t))'; % Hamming Window
rectWindow = rectwin(length(t))';  % Rectangular Window

figure('Name','Signal evaluation');
subplot(2,2,2)
plot(t, y),title('Signal with Noise'), ylim([-3 3]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on
% Hamming Window ACF
sig2 = y.*hammingWindow;

[r2, lags2] = xcorr(sig2,'biased');
tau2 = lags2/Fs;

% Power Spectral Density using Wiener Khintchine Theorem with Hamming window
Rxxdft = abs(fftshift(fft(r2)))/N;
freq = -Fs/2:Fs/length(r2):Fs/2-(Fs/length(r2));

subplot(2,2,3)
plot(freq,Rxxdft),title({'Linear scale:Power Spectral Density using FFT', 'with Hamming window'}),xlim([-50 50]), xlabel('Frequency f (in Hz)'),ylabel('Spectral Power')
grid on
subplot(2,2,4)
plot(freq,10.*log10(Rxxdft)),title({'Logarithm scale: PSD using FFT', 'with Hamming window'}),xlim([-50 50]), xlabel('Frequency f (in Hz)'),ylabel('Spectral Power')
grid on
%% Rectangular window
close all;
clc;
clear all;

alpha = 3.0;
N  = 2048; % Number of Samples
Fs = 100; % Sampling Frequency in Hz % Sampling Time(Ts) = 1/Fs = 1/250 = 0.004 sec
f = 20; % Signal Frequency
t = ((-N/2):(N/2)-1)/Fs; % Time axis 
x = sin(2*pi*f*t); % Sine Wave

noise = alpha.*randn(1,length(t)); % Noise
%%
% plot of original signal
subplot(2,2,1)
plot(t, x),title('Original Signal with no noise'), ylim([-3 3]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on 


% plot of sine signal with noise
y = x + noise; % Noise added to signal
subplot(2,2,2)
plot(t, y),title('Original Signal with noise'), ylim([-3 3]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on 

hammingWindow = hamming(length(t))'; % Hamming Window
rectWindow = rectwin(length(t))';  % Rectangular Window

figure('Name','Signal evaluation');
subplot(2,2,2)
plot(t, y),title('Signal with Noise'), ylim([-3 3]), xlabel('Time (in sec)'), ylabel('Amplitude')
grid on
% Hamming Window ACF
sig2 = y.*rectWindow;

[r2, lags2] = xcorr(sig2,'biased');
tau2 = lags2/Fs;

% Power Spectral Density using Wiener Khintchine Theorem with Hamming window
Rxxdft = abs(fftshift(fft(r2)))/N;
freq = -Fs/2:Fs/length(r2):Fs/2-(Fs/length(r2));

subplot(2,2,3)
plot(freq,Rxxdft),title({'Linear scale:Power Spectral Density using FFT', 'with Rectangular window'}),xlim([-50 50]), xlabel('Frequency f (in Hz)'),ylabel('Spectral Power')
grid on
subplot(2,2,4)
plot(freq,10.*log10(Rxxdft)),title({'Logarithm scale: PSD using FFT', 'with Rectangular window'}),xlim([-50 50]), xlabel('Frequency f (in Hz)'),ylabel('Spectral Power')
grid on