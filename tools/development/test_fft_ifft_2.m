
clear all
clc


%% user input
dt = 0.09;              % time step
fs = 1/dt;              % sampling frequency

% number of time samples for causal branch, i.e. from 0 to end
nt = 900;

n_noise_sources = 1;
f_peak = [0.5, 1];
bandwidth = [0.2, 0.2];
strength = [1, 0.9];


%% consider time series centered around 0
t = -(nt-1)*dt:dt:(nt-1)*dt;
nt = length(t);         % number of time samples

% set up frequency coordinates, plus some helper variables
f = linspace( -fs/2, fs/2, nt );
nf = length(f);
w = 2*pi*f;
dw = w(2) - w(1);


%% set up two sided power spectral density
psd = zeros( nf, n_noise_sources);
for ns = 1:n_noise_sources
    psd(:,ns) = strength(ns) * exp( -(abs(f)-f_peak(ns)).^2 / bandwidth(ns)^2 );
%     psd(:,ns) = double(f>=0)' .* psd(:,ns);
end


% plot psd
figure(1)
clf
subplot(2,1,1)
hold on
plot(f, psd, 'b')


%% on the fly inverse Fourier transform
ifft_coeff = zeros(nt,nf) + 1i*zeros(nt,nf);
for n = 1:nt
    for k = 1:nf
        ifft_coeff(n,k) = 1/sqrt(2*pi) * exp( 1i*w(k)*t(n) ) * dw;
    end
end

tic
x_time_otf = zeros(nt,1);
for n = 1:nt
    for k = 1:nf
            for ns = 1:n_noise_sources
                x_time_otf(n) = x_time_otf(n) + psd(k,ns) * ifft_coeff(n,k);
            end
    end
end
toc


%% inverse Fourier transform based on ifft
tic
x_time_ifft = zeros( 1, nt );
for ns = 1:n_noise_sources
    x_time_ifft(1,:) = x_time_ifft(1,:) + fftshift( ifft( ifftshift( psd(:,ns) ) * nt ) )';
end
toc


%% plotting
subplot(2,1,2)
hold on

% plot( t, real(x_time_otf), 'r' )
% plot( t, real(x_time_ifft), 'b' )

plot( t, real(x_time_otf)/max(real(x_time_otf)), 'r' )
plot( t, real(x_time_ifft)/max(real(x_time_ifft)), 'b' )

xlim([-10 10])



%% compute amplitude spectrum
x_fft_otf = fftshift( fft( ifftshift( x_time_otf ) ) / nt );
x_fft_ifft = fftshift( fft( ifftshift( x_time_ifft ) ) / nt );

figure(1)
subplot(2,1,1)
hold on

% plot( f_sample, sqrt(x_fft_otf.*conj(x_fft_otf)), 'r' )
% plot( f, sqrt(x_fft_ifft.*conj(x_fft_ifft)), 'r' )

% plot( f, abs(x_fft_otf), 'y' )
% plot( f, abs(x_fft_otf)/max(abs(x_fft_otf)), 'y' )

plot( f, abs(x_fft_ifft), 'g' )
% plot( f, abs(x_fft_ifft)/max(abs(x_fft_ifft)), 'g' )

