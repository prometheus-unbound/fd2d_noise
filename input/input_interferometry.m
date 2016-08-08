
function [f_sample, n_sample, w_sample, dw, freq_samp] = input_interferometry()


%==========================================================================
% input parameters for Fourier transform
%==========================================================================

f_sample = 0.02:0.002:0.2;


n_sample = length(f_sample);
w_sample = 2 * pi * f_sample;
dw = w_sample(2) - w_sample(1);
freq_samp = 5;


end