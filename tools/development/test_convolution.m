
clear all
clc
close all


[Lx, Lz, nx, nz, dt, nt, order, model_type, source_type, n_basis_fct, fw_nth] = input_parameters();
[f_sample, n_sample, w_sample, dw, freq_samp] = input_interferometry();

t = -(nt-1)*dt:dt:(nt-1)*dt;
nt = length(t);

fft_coeff = zeros(nt,n_sample) + 1i*zeros(nt,n_sample);
ifft_coeff = zeros(nt,n_sample) + 1i*zeros(nt,n_sample);

for n = 1:nt
    
    for k = 1:n_sample
        fft_coeff(n,k) =  1/sqrt(2*pi) * exp(-1i*w_sample(k)*t(n) ) * dt;
        ifft_coeff(n,k) = 1/sqrt(2*pi) * exp( 1i*w_sample(k)*t(n) ) * dw;
    end
    
end


n_noise_sources = 1;
f_peak = [1/11, 1/7];
bandwidth = [0.035, 0.025];
strength = [0.3, 1];


noise_spectrum = zeros(length(f_sample),n_noise_sources);

for ns = 1:n_noise_sources
    noise_spectrum(:,ns) = strength(ns) * exp( -(abs(f_sample)-f_peak(ns)).^2 / bandwidth(ns)^2 );
end

figure(1)
plot( f_sample, noise_spectrum )



source_time_fct = zeros(nt,1);
for n = 1:nt
    
    for ns = 1:n_noise_sources
        for k = 1:n_sample
            source_time_fct(n,1) = source_time_fct(n,1) + noise_spectrum(k,ns) * ifft_coeff(n,k);
        end
    end
    
end

source_time_fct = real(source_time_fct);


Green = zeros(nt,1);
Green(900,1) = 1.0;

test_trace = conv( source_time_fct, Green, 'same' );


figure(2)
hold on

plot( t, source_time_fct, 'r' )
plot( t, test_trace, 'b' )







% % Fourier Transform
% Nsamps = nt;
% Fs = 1/dt;
% 
% f = Fs*(0:Nsamps/2-1)/Nsamps;
% 
% data_fft = abs(fft(source_time_fct));
% data_fft = data_fft( 1:Nsamps/2 );
% 
% figure(3)
% plot(f,data_fft)
% xlim([0 0.2])



