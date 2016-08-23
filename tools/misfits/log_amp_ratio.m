%==========================================================================
% log amplitude ratio
%
% function [misfit, adjstf] = log_amp_ratio( u, u_0, win_caus, t )
%
% input:
%--------
% u: synthetic displacement seismogram
% u_0: observed displacement seismogram
% win_caus: window around surface wave train on causal branch
% t: time axis
%
% output:
%--------
% misfit
% adjstf: adjoint source time function
%
%==========================================================================


function [misfit, adjstf] = log_amp_ratio(u, u_0, win_caus, t)


    dt = abs(t(2) - t(1));


    %- check if win_caus is on the positive time axis ---------------------
    [~, index] = max(win_caus);
    if (t(index) < 0)
        win_acaus = win_caus;
        win_caus = fliplr(win_acaus);
    else
        win_acaus = fliplr(win_caus);
    end


    %- compute energy for windowed observations and synthetics ------------
    e_caus = trapz((win_caus .* u) .^ 2) * dt;
    e_acaus = trapz((win_acaus .* u) .^ 2) * dt + eps;
    e0_caus = trapz((win_caus .* u_0) .^ 2) * dt;
    e0_acaus = trapz((win_acaus .* u_0) .^ 2) * dt + eps;


    %- compute asymmetry --------------------------------------------------
    A = log(e_caus / e_acaus);
    A0 = log(e0_caus / e0_acaus);


    %- compute misfit -----------------------------------------------------
    misfit = 0.5 * (A - A0) ^ 2;

    
    %- compute adjoint source time function -------------------------------
    de_caus = 2 * win_caus .^ 2 .* u * dt;
    de_acaus = 2 * win_acaus .^ 2 .* u * dt;

    adjstf = (A - A0) * (de_caus / e_caus - de_acaus / e_acaus);


    %- time reverse adjoint source time function --------------------------
    adjstf = fliplr(adjstf);


end

