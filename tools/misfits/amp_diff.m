%==========================================================================
% compute relative amplitude difference
%
% function [misfit,adsrc] = amp_diff(u,u_0,t)
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


function [misfit,adstf] = amp_diff(u,u_0,win,t)


    dt = abs(t(2)-t(1));

    
    %- compute adjoint misfit and source time function --------------------
    if sum(u_0==0) == length(t)

        % not sure about that    
        misfit = 1.0;
        adstf = 2 * fliplr(u) / ( sum(u.^2)*dt );

    else

        measurement = ( sum( (win .* u).^2 ) - sum( (win .* u_0).^2 ) ) / sum( (win .* u_0).^2 );
        misfit = measurement^2 / 2;

        adstf = 2 * measurement * win.^2 .* u / sum( (win .* u_0).^2 );

    end
    
    
    % time reverse adjoint source time function ---------------------------
    adstf = fliplr(adstf);
    
    
end