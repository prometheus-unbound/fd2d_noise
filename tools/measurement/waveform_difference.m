%==========================================================================
% waveform difference
%
% function [misfit,adstf] = waveform_difference( u, u_0, du, win, t, deriv_order )
%
% input:
%--------
% u: synthetic displacement seismogram
% u_0: observed displacement seismogram
% du: perturbation
% win: specified window
% t: time axis
% deriv_order: '1st' or '2nd'
%
% output:
%--------
% misfit
% adstf: adjoint source time function
%
%==========================================================================


function [misfit, adstf] = waveform_difference( u, u_0, du, win, t, deriv_order )


    dt = abs( t(2) - t(1) );

    
    if( strcmp( deriv_order, '1st' ) )
    
        %- compute misfit and adjoint source time function ----------------
        misfit = 1/2 * sum( win.^2 .* (u - u_0).^2 ) * dt;
        adstf = win.^2 .* (u - u_0) * dt;
    
    else
        
        misfit = [];
        
        %- compute adjoint source time function ---------------------------
        adstf =  win.^2 .* du * dt;
        
    end
    
    % time reverse adjoint source time function ---------------------------  
    adstf = fliplr(adstf);
    
        
end
