
function [f_sample, n_sample, w_sample, dw, freq_samp] = input_interferometry()

%==========================================================================
% input parameters for interferometry
%==========================================================================

%- frequency sampling in Hz -----------------------------------------------
%- The sampling should be evenly spaced for the inverse F transform.
%- The sampling must also be sufficiently dense in order to avoid artefacts
%- on the positive time axis in the time-domain source function. This can 
%- be checked with "/tools/check_correlation_source_function".

%- It is sufficient to consider the positive frequency axis. 

f_sample=0.02:0.002:0.2;

n_sample = length(f_sample);
w_sample = 2*pi*f_sample;
dw = w_sample(2) - w_sample(1);
freq_samp = 1;
