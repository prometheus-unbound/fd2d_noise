%==========================================================================
% waveform difference
%
% function [misfit, adjstf] = waveform_difference( u, u_0, win, t )
%
% input:
%--------
% u: synthetic displacement seismogram
% u_0: observed displacement seismogram
% win: specified window
% t: time axis
%
% output:
%--------
% misfit
% adjstf: adjoint source time function
%
%==========================================================================


function [misfit, adjstf] = waveform_difference( u, u_0, win, t )


    dt = abs( t(2) - t(1) );


    %- compute misfit and adjoint source time function --------------------
    misfit = 1/2 * sum( win.^2 .* (u - u_0).^2 ) * dt;
    adjstf = win.^2 .* (u - u_0) * dt;


    % time reverse adjoint source time function ---------------------------
    adjstf = fliplr(adjstf);


end
