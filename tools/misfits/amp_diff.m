%==========================================================================
% relative amplitude difference
%
% function [misfit, adjstf] = amp_diff( u, u_0, win, t )
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


function [misfit, adjstf] = amp_diff( u, u_0, win, t )


    dt = abs(t(2)-t(1));

    
    %- compute energy in windowed observation -----------------------------
    energy_0 = sum( (win .* u_0).^2 ) * dt;
    
    
    %- compute measurement ------------------------------------------------
    measurement = ( sum( (win .* u).^2 ) * dt - energy_0 ) / ( energy_0 + eps );
    

    %- compute misfit -----------------------------------------------------
    misfit = measurement^2 / 2;
    
    
    %- compute adoint source time function --------------------------------
    adjstf = 2 * measurement * win.^2 .* u * dt / ( energy_0 + eps );
    
    
    % time reverse adjoint source time function ---------------------------
    adjstf = fliplr(adjstf);
    
    
end
