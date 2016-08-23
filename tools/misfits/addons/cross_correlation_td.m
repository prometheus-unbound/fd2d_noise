%==========================================================================
% compute cross correlation function in time domain
%
% [cc,t_cc] = cross_correlation_td(f,g,t)
%
%==========================================================================

function [cc, t_cc] = cross_correlation_td(f, g, t)


    %- initialisations ----------------------------------------------------
    n = length(f);
    dt = t(2) - t(1);
    t_cc = - (n - 1) * dt:dt:(n - 1) * dt;

    %- compute brute-force correlation function ---------------------------
    cc = xcorr(g, f);
   
    
end
