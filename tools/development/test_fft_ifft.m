
clear all
clc

Fs = 1000;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = 1000;             % Length of signal
t = (0:L-1)*T;        % Time vector


x_time_orig = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);


figure(1)
subplot(2,2,1)
plot(t,x_time_orig,'kx')


% Fourier transform
x_fft = fftshift( fft(x_time_orig,L) );
x_fft2 = fft(x_time_orig);


% compute one sided spectrum
P1 = abs( x_fft / L );
P2 = abs( x_fft2 / L );

% P3 = P2(1:L/2+1);
% P3(2:end-1) = 2*P3(2:end-1);




subplot(2,2,2)
hold on

f = ((-L/2):(L/2)-1) * Fs/L;

plot(f,P1,'b')
plot(f,P2,'r')

title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


% inverse Fourier transform
x_time = ifft( ifftshift( x_fft ) );
x_time_2 = ifft( x_fft2 );


subplot(2,2,3)
hold on
plot(t,x_time,'rx')
plot(t,x_time_2,'bx')




