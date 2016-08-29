%==========================================================================
% relative amplitude difference
%
% function [misfit,adsrc] = amp_diff(u,u_0,t)
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


function [misfit, adstf] = amp_diff( u, u_0, du, win, t, deriv_order )


    dt = abs(t(2)-t(1));

    
    %- compute energy in windowed observation -----------------------------
    energy_0 = sum( (win .* u_0).^2 ) * dt;
    
    
    %- compute measurement ------------------------------------------------
    measurement = ( sum( (win .* u).^2 ) * dt - energy_0 ) / ( energy_0 + eps );
    
    
    if( strcmp( deriv_order, '1st' ) )
    
        if( sum(u_0==0) == length(t) )            
            misfit = 1.0;
            adstf = 2 * fliplr(u) / ( sum(u.^2)*dt );            
        else
            
            %- compute misfit ---------------------------------------------
            misfit = measurement^2 / 2;
            
            %- compute adoint source time function ------------------------
            adstf = 2 * measurement * win.^2 .* u / sum( (win .* u_0).^2 );
            
        end
    
    else

        misfit = [];
        
        %- compute adoint source time function ----------------------------
        adstf = 2 * measurement * win.^2 .* du / energy_0 * dt + ...
            4 * win.^2 .* u / energy_0^2 * sum( win.^2 .* u .* du ) * dt^2;
        
    end
    
    
    % time reverse adjoint source time function ---------------------------
    adstf = fliplr(adstf);
    
    
end
