%==========================================================================
% compute adjoint sources
%
% function [misfit_n,adstf] = make_adjoint_sources(u,u_0,t,veldis,measurement,src,rec)
%
% input:
%-------
% u: synthetic displacement seismograms
% u_0: "observed" displacement seismograms
% t: time axis
% veldis: 'dis' for displacements, 'vel' for velocities
% measurement:  'log_amplitude_ratio'
%               'amplitude_difference'
%               'waveform_difference' 
%               'cc_time_shift'
% src: source position
% rec: receiver positions
%
% output:
%-------
% misfit_n: misfit of each receivers
% adstf: adjoint source time functions for each receiver
%
%
% When u_0, i.e. the observed displacement seismograms, are set to zero, 
% the code performs data-independent measurements. 
% 
%==========================================================================


function [misfit_n, adstf] = make_adjoint_sources( u, u_0, du, t, veldis, measurement, src, rec, deriv_order )

%==========================================================================
%- initialisations --------------------------------------------------------
%==========================================================================

nt = length(t);
dt = abs( t(2) - t(1) );
n_receivers = size( rec, 1 );

%- convert to velocity if wanted ------------------------------------------
if strcmp( veldis, 'vel' )
    
    v = zeros( n_receivers, nt );
    v_0 = zeros( n_receivers, nt );
    
    for k=1:n_receivers
        v(k,1:nt-1) = diff( u(k,:) ) / dt;
        v_0(k,1:nt-1) = diff( u_0(k,:) ) / dt;
    end
   
    u = v;
    u_0 = v_0;
    
end


%==========================================================================
%- march through the various recordings -----------------------------------
%==========================================================================

misfit_n = zeros( n_receivers, 1 );
adstf = zeros( n_receivers, nt );

for n=1:n_receivers
   
    %- select time windows ------------------------------------------------   
    % disp('select left window');
    % [left,~] = ginput(1)
    % disp('select_right_window');
    % [right,~] = ginput(1)

    if strcmp(measurement,'waveform_difference')
        left = t(1);
        right = t(end);
    
    else        
        distance = sqrt( (src(1,1) - rec(n,1)).^2 + (src(1,2) - rec(n,2)).^2 );
        left = distance/4000.0 - 27.0;
        right = distance/4000.0 + 27.0;
        
        if( left < 0 )
            index = find( t==0 );
            left = t(index+1);
        end
        
        if( right > t(end) )
            right = t(end);
        end
        
        % error if time series is not long enough
        if( left > t(end) )
            error('time series too short?')
        end
        
    end


    win = get_window( t, left, right, 'cos_taper' );
    
    
    %- compute misfit and adjoint source time function --------------------    
    if strcmp(measurement,'waveform_difference')
        
        if( strcmp( deriv_order, '1st' ) )
            [misfit_n(n,:), adstf(n,:)] = waveform_difference( u(n,:), u_0(n,:), du(n,:), win, t, deriv_order );
        else
            [~, adstf(n,:)] = waveform_difference( u(n,:), u_0(n,:), du(n,:), win, t, deriv_order );
        end
        
        
    elseif strcmp(measurement,'cc_time_shift')
        
        if( strcmp( deriv_order, '2nd' ) )
            error('\nno Hessian vector products for traveltime measurements\n')
        end
        
        [misfit_n_caus, adstf_caus(1,:)] = cc_time_shift( u(n,:), u_0(n,:), win, t );
        
        [left, right] = swap( -left, -right );
        win = get_window( t, left, right, 'cos_taper' );       
        [misfit_n_acaus, adstf_acaus(1,:)] = cc_time_shift( u(n,:), u_0(n,:), win, t );
        
        if( strcmp( deriv_order, '1st' ) )
            misfit_n(n,:) = misfit_n_caus + misfit_n_acaus;
        end
        
        adstf(n,:) = adstf_caus + adstf_acaus;
        
        
    elseif strcmp(measurement,'amplitude_difference')       
        
        [misfit_n_caus, adstf_caus(1,:)] = amp_diff( u(n,:), u_0(n,:), du(n,:), win, t, deriv_order );
        
        [left, right] = swap( -left, -right );
        win = get_window( t, left, right, 'cos_taper' );       
        [misfit_n_acaus, adstf_acaus(1,:)] = amp_diff( u(n,:), u_0(n,:), du(n,:), win, t, deriv_order );
        
        if( strcmp( deriv_order, '1st' ) )
            misfit_n(n,:) = misfit_n_caus + misfit_n_acaus;
        end
        
        adstf(n,:) = adstf_caus + adstf_acaus;
        
        
    elseif strcmp(measurement,'log_amplitude_ratio')    
        
        win = get_window( t, left, right, 'hann' );
        
        if( strcmp( deriv_order, '1st' ) )
            [misfit_n(n,:), adstf(n,:)] = log_amp_ratio( u(n,:), u_0(n,:), du(n,:), win, t, deriv_order );
        else
            [~, adstf(n,:)] = log_amp_ratio( u(n,:), u_0(n,:), du(n,:), win, t, deriv_order );
        end
        
        
    end
     
    
    %- correct adjoint source time function for velocity measurement ------    
    if strcmp(veldis,'vel')
        adstf( n, 1:nt-1 ) = -diff( adstf(n,:) ) / dt;
    end      

    
end


