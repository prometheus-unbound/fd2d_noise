%==========================================================================
% compute waveform difference
%
% function [misfit,adstf] = waveform_difference(u,u_0,t)
%
% input:
%--------
% u: synthetic displacement seismogram
% u_0: observed displacement seismogram
% t: time axis
%
% output:
%--------
% misfit
% adstf: adjoint source time function
%
%==========================================================================


function [misfit, adstf] = waveform_difference( u, u_0, win, t )

    dt = abs( t(2) - t(1) );
    
    
    %- compute adjoint misfit and source time function --------------------    
    misfit = 1/2 * sum( win.^2 .* (u - u_0).^2 ) * dt;
    adstf = win.^2 .* (u - u_0) * dt;
    
    
    % time reverse adjoint source time function ---------------------------  
    adstf = fliplr(adstf);
    
        
end