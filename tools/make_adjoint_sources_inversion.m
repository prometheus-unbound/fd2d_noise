%==========================================================================
% compute adjoint sources
%
% function [misfit,adstf] = make_adjoint_sources_inversion(u,u_0,t,veldis,measurement,src,rec)
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
% i_ref: number of reference station
% flip_sr: flip source and receiver, important for structure kernel
%
% output:
%-------
% misfit: cumulative misfit for all receivers
% adstf: adjoint source time functions for each receiver
%
%
% When u_0, i.e. the observed displacement seismograms, are set to zero, 
% the code performs data-independent measurements. 
% 
%==========================================================================


function [misfit,adstf] = make_adjoint_sources_inversion(u,u_0,t,veldis,measurement,src,rec)

%==========================================================================
%- initialisations --------------------------------------------------------
%==========================================================================

nt = length(t);
dt = t(2) - t(1);
n_receivers = size(rec,1);

%- convert to velocity if wanted ------------------------------------------
if strcmp(veldis,'vel')
    
    v = zeros(n_receivers,nt);
    v_0 = zeros(n_receivers,nt);
    
    for k=1:n_receivers
        v(k,1:nt-1) = diff(u(k,:))/(t(2)-t(1));
        v_0(k,1:nt-1) = diff(u_0(k,:))/(t(2)-t(1));
    end
   
    u = v;
    u_0 = v_0;
    
end


%==========================================================================
%- march through the various recodings ------------------------------------
%==========================================================================

misfit = 0.0;
adstf = zeros(n_receivers,nt);

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
        
    end


    win = get_window(t,left,right,'cos_taper');
    
    
    %- compute misfit and adjoint source time function --------------------    
    if strcmp(measurement,'waveform_difference')
        
        [misfit_n,adstf(n,:)] = waveform_difference(u(n,:),u_0(n,:),win,t);
        
        
    elseif strcmp(measurement,'cc_time_shift') 
        
        [misfit_n_caus,adstf_caus(1,:)] = cc_time_shift(u(n,:),u_0(n,:),win,t);
        
        [left,right] = swap(-left,-right);
        win = get_window(t,left,right,'cos_taper');       
        [misfit_n_acaus,adstf_acaus(1,:)] = cc_time_shift(u(n,:),u_0(n,:),win,t);
        
        misfit_n = misfit_n_caus + misfit_n_acaus;
        adstf(n,:) = adstf_caus + adstf_acaus;
        
        
    elseif strcmp(measurement,'amplitude_difference')       
        
        [misfit_n_caus,adstf_caus(1,:)] = amp_diff(u(n,:),u_0(n,:),win,t);
        
        [left,right] = swap(-left,-right);
        win = get_window(t,left,right,'cos_taper');       
        [misfit_n_acaus,adstf_acaus(1,:)] = amp_diff(u(n,:),u_0(n,:),win,t);
        
        misfit_n = misfit_n_caus + misfit_n_acaus;
        adstf(n,:) = adstf_caus + adstf_acaus;
        
        
    elseif strcmp(measurement,'log_amplitude_ratio')    
        
        win = get_window(t,left,right,'hann');
        [misfit_n,adstf(n,:)] = log_amp_ratio(u(n,:),u_0(n,:),win,t);
        
    end
    
    
    misfit = misfit + misfit_n;
    
    
    %- correct adjoint source time function for velocity measurement ------    
    if strcmp(veldis,'vel')
        adstf(n,1:nt-1) = -diff(adstf(n,:))/dt;
    end      

    
end


