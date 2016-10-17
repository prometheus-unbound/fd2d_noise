
%==========================================================================
% input parameters for Fourier transform
%
% [f_sample, n_sample, w_sample, dw, freq_samp] = input_interferometry()
%
% output:
%--------
% f_sample: frequency vector [Hz]
% n_sample: number of frequency samples
% w_sample: frequency vector converted to angular frequency
% dw: step in angular frequency
% freq_samp: calculate Fourier transform every freq_samp time step
%
%==========================================================================


function [f_sample, n_sample, w_sample, dw, freq_samp] = input_interferometry()


    f_sample = 0.04:0.002:0.18;      % in Hz
    freq_samp = 5;


    %- derived quantities -------------------------------------------------
    n_sample = length(f_sample);
    w_sample = 2 * pi * f_sample;
    dw = w_sample(2) - w_sample(1);


end
