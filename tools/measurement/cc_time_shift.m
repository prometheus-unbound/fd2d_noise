
%==========================================================================
% cross correlation time shift
%
% function [misfit, adjstf] = cc_time_shift( u, u_0, win, t )
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


function [misfit, adjstf] = cc_time_shift(u, u_0, win, t)


    %- compute time shift -------------------------------------------------
    % if sum(u_0==0) == length(t)
    %     
    %     T = 1.0;
    %     
    % else
        
        [cc, t_cc] = cross_correlation_td(win .* u, win .* u_0, t);
        [~, i_max] = max(cc);
        T = t_cc(i_max);
        
        if (abs(T) > 3.5)
            T = 0;
        end
        
    % end


    %- compute misfit -----------------------------------------------------
    misfit = T ^ 2 / 2.0;


    %- compute adjoint source time function -------------------------------
    dt = abs(t(2) - t(1));
    nt = length(t);

    v = zeros(1, nt);
    v(1:nt - 1) = diff(u) / dt;

    adjstf = T * (win .^ 2 .* v) / (sum((win .* v) .^ 2) * dt);


    % time reverse adjoint source time function ---------------------------
    adjstf = fliplr(adjstf);


end


