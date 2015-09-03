function [displacement_seismograms,t] = run_forward(simulation_mode,src,rec,i_ref,flip_sr)

tic

%==========================================================================
% run forward simulation
%
% input:
%--------
% simulation_mode: 'forward', 'forward_green' or 'correlation'
% src: source, i.e. the reference station
% rec: receivers
% i_ref: number of reference station
% flip_sr: flip source and receiver, important for structure kernel
%
% output:
%--------
% u: displacement seismograms
% t: time axis
%
%==========================================================================

cm = cbrewer('div','RdBu',100,'PCHIP');


%==========================================================================
% flip sources and receivers if wanted
%==========================================================================
if( strcmp(flip_sr,'yes') )
    tmp = src;
    src = rec;
    rec = tmp;
end
clear tmp;


%==========================================================================
% initialise simulation
%==========================================================================

%- material and domain ----------------------------------------------------
[Lx,Lz,nx,nz,dt,nt,order,model_type,~,n_basis_fct] = input_parameters();
[X,Z,x,z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);
[mu,rho] = define_material_parameters(nx,nz,model_type); 

output_specs
if (strcmp(make_plots,'yes'))
    plot_model;
end

%- compute indices for receiver locations ---------------------------------
n_receivers = size(rec,1);
rec_id = zeros(n_receivers,2);
for i=1:n_receivers    
    rec_id(i,1) = min( find( min(abs(x-rec(i,1))) == abs(x-rec(i,1)) ) );
    rec_id(i,2) = min( find( min(abs(z-rec(i,2))) == abs(z-rec(i,2)) ) );   
end


%- initialise interferometry ----------------------------------------------
[~, n_sample, w_sample, dw] = input_interferometry();


%- forward simulations ('forward', 'forward_green') -----------------------
if (strcmp(simulation_mode,'forward') || strcmp(simulation_mode,'forward_green'))

    %- time axis
    t = 0:dt:dt*(nt-1);     
    
    
    %- compute indices for source locations
    ns = size(src,1);    
    src_id = zeros(ns,2);
    for i = 1:ns
        src_id(i,1) = min( find( min(abs(x-src(i,1))) == abs(x-src(i,1)) ) );
        src_id(i,2) = min( find( min(abs(z-src(i,2))) == abs(z-src(i,2)) ) );
    end        
    
    
    %- make source time function
    make_source_time_function;
    
    
    %- initialise interferometry ------------------------------------------    
    if strcmp(simulation_mode,'forward_green')    
        
        %- Fourier transform of the forward Greens function
        G_2 = zeros(nx,nz,n_sample) + 1i*zeros(nx,nz,n_sample);
        
        
        % prepare coefficients for Fourier transform
        fft_coeff = zeros(length(t),n_sample) + 1i*zeros(length(t),n_sample);
        for k = 1:n_sample
            fft_coeff(:,k) = exp(-1i*w_sample(k)*t')*dt;
        end
        
    end
    
    
%- forward simulation to compute correlation function ---------------------
elseif strcmp(simulation_mode,'correlation')

    %- time axis
    t = -(nt-1)*dt:dt:(nt-1)*dt; 
    
    
    %- Fourier transform of the correlation velocity field
    C_2 = zeros(nx,nz,n_sample) + 1i*zeros(nx,nz,n_sample);
    
    
    %- Fourier transform of strain field
    C_2_dxv = zeros(nx-1,nz,n_sample) + 1i*zeros(nx-1,nz,n_sample);
    C_2_dzv = zeros(nx,nz-1,n_sample) + 1i*zeros(nx,nz-1,n_sample);
    
    
    %- load frequency-domain Greens function
    if( strcmp(flip_sr,'no') )
        load(['../output/interferometry/G_2_' num2str(i_ref) '.mat']);
    else
        load(['../output/interferometry/G_2_flip_sr_' num2str(i_ref) '.mat']);
    end
    
    
    % prepare coefficients for Fourier transform and its inverse
    fft_coeff = zeros(length(t),n_sample) + 1i*zeros(length(t),n_sample);
    ifft_coeff = zeros(length(t),n_sample) + 1i*zeros(length(t),n_sample);
    for k = 1:n_sample
        G_2(:,:,k) = conj(G_2(:,:,k));
        fft_coeff(:,k) = exp( -1i*w_sample(k)*t' )*dt;
        ifft_coeff(:,k) = dw*exp( 1i*w_sample(k)*t' )/(2*pi);
    end
    
    
    %- initialise noise source locations and spectra
    [noise_source_distribution, noise_spectrum] = make_noise_source(make_plots);
    n_noise_sources = size(noise_spectrum,2);
    
    
    % get integration boundaries for frequency bands -------------------------- 
    if( n_basis_fct ~= 0 )
        [int_limits] = integration_limits(n_sample,n_basis_fct);
    end
    
end
 

%- dynamic fields and absorbing boundary field ----------------------------
v = zeros(nx,nz);
sxy = zeros(nx-1,nz);
szy = zeros(nx,nz-1);


%- initialise seismograms -------------------------------------------------
velocity_seismograms = zeros(n_receivers,nt);


%- initialise absorbing boundary taper a la Cerjan ------------------------
[absbound] = init_absbound();


%==========================================================================
% iterate
%==========================================================================

if (strcmp(make_plots,'yes'))
    figure;
end

t_fft = 0;
t_ifft = 0;


for n = 1:length(t)
    
    %- compute divergence of current stress tensor ------------------------    
    DS = div_s(sxy,szy,dx,dz,nx,nz,order);
    
    
    %- add point sources --------------------------------------------------    
    if (strcmp(simulation_mode,'forward') || strcmp(simulation_mode,'forward_green'))
    
        for i=1:ns
            DS(src_id(i,1),src_id(i,2)) = DS(src_id(i,1),src_id(i,2)) + stf(n);
        end
        
    end
    
    
    %- add sources of the correlation field -------------------------------    
    if( mod(n,5) == 1 && strcmp(simulation_mode,'correlation') && (t(n)<0.0) )
        
        %- transform on the fly to the time domain 
        S = zeros(nx,nz,n_noise_sources) + 1i*zeros(nx,nz,n_noise_sources);

        
%         % noise source for coupled spectra and distributions
%         for ns = 1:n_noise_sources
%             
%             t_ifft_start = tic;
%             
%             for k = 1:n_sample
%                 % S(:,:,ns) = S(:,:,ns) + noise_spectrum(k,ns) * conj(G_2(:,:,k)) * exp( 1i*w_sample(k)*t(n) );
%                 S(:,:,ns) = S(:,:,ns) + noise_spectrum(k,ns) * G_2(:,:,k) * ifft_coeff(n,k);
%             end
%             % S(:,:,ns) = dw*S(:,:,ns)/pi; 
%             
%             t_ifft = t_ifft + toc(t_ifft_start);
%             
%             DS = DS + noise_source_distribution(:,:,ns) .* real(S(:,:,ns));
% 
%         end
        
        
%         % general source version
%         for k = 1:n_sample
%             S(:,:) = S(:,:) + noise_source_distribution(:,:,k) .* G_2(:,:,k) * ifft_coeff(n,k);
%         end
%         DS = DS + real(S);
        
        
        % nicer implementation
        for ns = 1:n_noise_sources
            
            for k = 1:n_sample
                
                if( n_basis_fct ~= 0 )
                    ib = find( k >= int_limits(:,1) & k <= int_limits(:,2) );
                else
                    ib = ns;
                end
                
                S(:,:,ns) = S(:,:,ns) + noise_spectrum(k,ns) * noise_source_distribution(:,:,ib) .* G_2(:,:,k) * ifft_coeff(n,k);
                
            end
            
            DS = DS + real(S(:,:,ns));

        end
        
           
    end
    
    
    %- update velocity field ----------------------------------------------
    v = v + dt*DS./rho;
    
    
    %- apply absorbing boundary taper -------------------------------------    
    v = v.*absbound;
    
    
    %- compute derivatives of current velocity and update stress tensor ---
    strain_dxv = dx_v(v,dx,dz,nx,nz,order);
    strain_dzv = dz_v(v,dx,dz,nx,nz,order);
    
    sxy = sxy + dt*mu(1:nx-1,:) .* strain_dxv;
    szy = szy + dt*mu(:,1:nz-1) .* strain_dzv;
    
    
    %- record velocity seismograms ----------------------------------------    
    for k = 1:n_receivers
        velocity_seismograms(k,n) = v(rec_id(k,1),rec_id(k,2));
    end
    
    
    %- accumulate Fourier transform of the displacement Greens function ---
    if( mod(n,5) == 1 && strcmp(simulation_mode,'forward_green') )
        t_fft_start = tic;           
        
        for k=1:n_sample
            % G_2(:,:,k) = G_2(:,:,k) + v(:,:) * exp(-1i*w_sample(k)*t(n))*dt;
            G_2(:,:,k) = G_2(:,:,k) + v(:,:) * fft_coeff(n,k);
        end
        
        t_fft = t_fft + toc(t_fft_start);
    end
    
    
    %- accumulate Fourier transform of the correlation velocity field -----    
    if( mod(n,5) == 1 && strcmp(simulation_mode,'correlation') )
        t_fft_start = tic; 
        
        for k=1:n_sample
            % C_2(:,:,k) = C_2(:,:,k) + v(:,:) * exp(-1i*w_sample(k)*t(n))*dt;
            C_2(:,:,k) = C_2(:,:,k) + v(:,:) * fft_coeff(n,k);
            
            C_2_dxv(:,:,k) = C_2_dxv(:,:,k) + strain_dxv(:,:) * fft_coeff(n,k);
            C_2_dzv(:,:,k) = C_2_dzv(:,:,k) + strain_dzv(:,:) * fft_coeff(n,k);
        end
        
        t_fft = t_fft + toc(t_fft_start);
    end
    
    
    %- plot velocity field ------------------------------------------------
    if (strcmp(make_plots,'yes'))
        plot_velocity_field;
    end
    
end


%==========================================================================
% output 
%==========================================================================


%- store Fourier transformed velocity Greens function ---------------------
if strcmp(simulation_mode,'forward_green')
    if( strcmp(flip_sr,'no') )
        save(['../output/interferometry/G_2_' num2str(i_ref) '.mat'],'G_2');        
    else
        save(['../output/interferometry/G_2_flip_sr_' num2str(i_ref) '.mat'],'G_2');
    end
end


%- store Fourier transformed correlation velocity field -------------------
if strcmp(simulation_mode,'correlation')
    if( strcmp(flip_sr,'no') )
        save(['../output/interferometry/C_2_' num2str(i_ref) '.mat'],'C_2');
        save(['../output/interferometry/C_2_dxv_' num2str(i_ref) '.mat'],'C_2_dxv');
        save(['../output/interferometry/C_2_dzv_' num2str(i_ref) '.mat'],'C_2_dzv');
    else
        save(['../output/interferometry/C_2_flip_sr_'  num2str(i_ref) '.mat'],'C_2');
        save(['../output/interferometry/C_2_dxv_flip_sr_'  num2str(i_ref) '.mat'],'C_2_dxv');
        save(['../output/interferometry/C_2_dzv_flip_sr_'  num2str(i_ref) '.mat'],'C_2_dzv');
    end
end


%- displacement seismograms -----------------------------------------------
displacement_seismograms = cumsum(velocity_seismograms,2)*dt;

% displacement_seismograms = 0.0*velocity_seismograms;
% index_zero = find( t==0 );
% displacement_seismograms(:,1:index_zero) = fliplr( cumsum(fliplr( velocity_seismograms(:,1:index_zero) ),2) ) * (-dt);
% displacement_seismograms(:,index_zero:end) = cumsum( velocity_seismograms(:,index_zero:end),2) * dt;


%- store the movie if wanted ----------------------------------------------
if strcmp(make_movie,'yes')
    writerObj = VideoWriter(movie_file,'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,M);
    close(writerObj);
end



t_total = toc;

if( strcmp(verbose,'yes') )
    fprintf('\ntotal time: %f\n',t_total)
    fprintf('time fft:   %f\n',t_fft)
    fprintf('percentage: %f\n\n',t_fft/t_total*100)
    fprintf('time ifft:  %f\n',t_ifft)
    fprintf('percentage: %f\n\n',t_ifft/t_total*100)
end

% clean up
close all


end

