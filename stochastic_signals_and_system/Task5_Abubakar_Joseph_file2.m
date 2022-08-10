clear all;
close all;

N  = 1024; % Number of Samples
Fs = 100; % Sampling Frequency in Hz 
T = 1; % period
w = 0 : pi/T : 6.*pi./T % phase angle
y = (2./((w.^2).*T)).*((1-(2.*cos(w.*T))).^2)

% Power Spectral Density using Wiener Khintchine Theorem with Hamming window
Rxxdft = abs(fftshift(fft(y)));
freq = -Fs/2:Fs/length(y):Fs/2-(Fs/length(y));

subplot(2,1,1);
plot(freq, Rxxdft,'b'),title('Power spectral density'),axis([-1 1 0 2]),xlabel('Phase'),ylabel('PSD');
grid on;