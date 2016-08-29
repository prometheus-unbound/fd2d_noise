%==========================================================================
% log amplitude ratio
%
% [misfit, adstf] = log_amp_ratio( u, u_0, du, win_caus, t, deriv_order )
%
% input:
%--------
% u: synthetic displacement seismogram
% u_0: observed displacement seismogram
% du: perturbation
% win_caus: window around surface wave train on causal branch
% t: time axis
% deriv_order: '1st' or '2nd'
%
% output:
%--------
% misfit
% adstf: adjoint source time function
%
%==========================================================================


function [misfit, adstf] = log_amp_ratio( u, u_0, du, win_caus, t, deriv_order )


    dt = abs(t(2)-t(1));


    % check if win_caus is on the positive time axis ----------------------
    [~,index] = max(win_caus);
    if( t(index) < 0 )
        win_acaus = win_caus;
        win_caus = fliplr(win_acaus);
    else
        win_acaus = fliplr(win_caus);
    end


    % compute energy for windowed observations and synthetics -------------
    e_caus = trapz( (win_caus .* u).^2 ) * dt;
    e_acaus = trapz( (win_acaus .* u).^2 ) * dt + eps;
    e0_caus = trapz( (win_caus .* u_0).^2 ) * dt;
    e0_acaus = trapz( (win_acaus .* u_0).^2 ) * dt + eps;


    % compute asymmetry ---------------------------------------------------
    A = log( e_caus/e_acaus );
    A0 = log( e0_caus/e0_acaus );

    
    if( strcmp( deriv_order, '1st' ) )
        
        % compute misfit --------------------------------------------------
        misfit = 0.5 * ( A - A0 )^2;
        
        % compute adjoint source time function ----------------------------
        de_caus = 2 * win_caus.^2 .* u * dt;
        de_acaus = 2 * win_acaus.^2 .* u * dt;
        
        if ( sum(u_0==0) == length(t) )
            adstf = de_caus/e_caus - de_acaus/e_acaus;
        else
            adstf = ( A - A0 ) * ( de_caus/e_caus - de_acaus/e_acaus );
        end
        
    else
        
        misfit= [];
        
        % compute adjoint source time function ----------------------------        
        de_caus = trapz( 2 * win_caus.^2 .* u .* du ) * dt;
        de_acaus = trapz( 2 * win_acaus.^2 .* u .* du ) * dt;       
        
        adstf = ( A - A0 ) * ( 2/e_caus * win_caus.^2 - 2/e_acaus * win_acaus.^2 ) .* du * dt ...
            + ( A - A0 ) * ( -2/e_caus^2 * win_caus.^2 .* u .* de_caus + 2/e_acaus^2 * win_acaus.^2 .* u .* de_acaus ) * dt ...
            + ( 2/e_caus * win_caus.^2 .* u - 2/e_acaus * win_acaus.^2 .* u ) .* ( de_caus/e_caus - de_acaus/e_acaus ) * dt;
        
    end

    
    % time reverse adjoint source time function ---------------------------
    adstf = fliplr(adstf);


end
