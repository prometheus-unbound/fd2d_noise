
clear all
clc


load sensorData




s1 = s1(1:length(s2));

s1 = repmat(s1,500,1);
s2 = repmat(s2,500,1);

s1_p = [ s1(1:length(s2)); zeros(length(s2)-1,1)];
s2_p = [zeros(length(s2)-1,1); s2];


% Do Fourier Transform
tic
s1_fft = fft(s1_p);
s2_fft = fft(s2_p);

corr_fft = conj(s1_fft) .* s2_fft;
corr = ifft( corr_fft );
toc

tic
[acor,lag] = xcorr(s2,s1);
toc


figure
hold on
plot(corr,'b')
plot(acor,'r')