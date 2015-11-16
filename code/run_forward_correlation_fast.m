function [displacement_seismograms,t,C_2_dxu,C_2_dzu] = run_forward_correlation_fast(G_2, source_dist, noise_spectrum, mu, rec, mode, df)

%==========================================================================
% run forward correlation
% fast means ready for conversion to mex-files
%
% input:
%--------
% G_2: Green function of reference station
% source_dist: source distribution
% noise_spectrum: spectrum of noise distribution
% mu [N/m^2]
% rec: receivers
% mode: 1 = calculate fourier transform of strain of correlation wavefield
% df: perturbation of rhs, implemented for check of adjoint state
%
% output:
%--------
% correlation recordings
% t: time vector
% C_2_dxu: stress component of correlation wavefield
% C_2_dzu: stress component of correlation wavefield
%
%==========================================================================


%==========================================================================
% initialise simulation
%==========================================================================

%- material and domain ----------------------------------------------------
[Lx,Lz,nx,nz,dt,nt,order,model_type,~,n_basis_fct] = input_parameters();
[~,~,x,z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);
[~,rho] = define_material_parameters(nx,nz,model_type); 
mu = reshape(mu, nx, nz);


%- initialise interferometry ----------------------------------------------
[~,n_sample,w_sample,dw,freq_samp] = input_interferometry();


%- reshape source distribution --------------------------------------------
n_noise_sources = size(noise_spectrum,2);

if( n_basis_fct == 0 )
    noise_source_distribution = reshape(source_dist, nx, nz, n_noise_sources);
else
    noise_source_distribution = reshape(source_dist, nx, nz, n_basis_fct);
end


%- get integration boundaries for frequency bands -------------------------
if( n_basis_fct ~= 0 )
    [int_limits] = integration_limits(n_sample,n_basis_fct);
end


%- compute indices for receiver locations ---------------------------------
n_receivers = size(rec,1);
rec_id = zeros(n_receivers,2);

for i=1:n_receivers    
    rec_id(i,1) = min( find( min(abs(x-rec(i,1))) == abs(x-rec(i,1)) ) );
    rec_id(i,2) = min( find( min(abs(z-rec(i,2))) == abs(z-rec(i,2)) ) );   
end


%==========================================================================
%- forward simulation to compute correlation function ---------------------
%==========================================================================

%- time axis --------------------------------------------------------------
t = -(nt-1)*dt:dt:(nt-1)*dt;
nt = length(t);


%- prepare coefficients for Fourier transform and its inverse -------------
fft_coeff = zeros(length(t),n_sample) + 1i*zeros(length(t),n_sample);
ifft_coeff = zeros(length(t),n_sample) + 1i*zeros(length(t),n_sample);
for k = 1:n_sample
    G_2(:,:,k) = conj(G_2(:,:,k));
    fft_coeff(:,k) = 1/sqrt(2*pi) * exp( -1i*w_sample(k)*t' ) * dt;
    ifft_coeff(:,k) = 1/sqrt(2*pi) * exp( 1i*w_sample(k)*t' ) * dw;
end


%- Fourier transform of the correlation velocity field --------------------
% C_2 = zeros(nx,nz,n_sample) + 1i*zeros(nx,nz,n_sample);

%- Fourier transform of strain field
C_2_dxu = zeros(nx-1,nz,n_sample) + 1i*zeros(nx-1,nz,n_sample);
C_2_dzu = zeros(nx,nz-1,n_sample) + 1i*zeros(nx,nz-1,n_sample);
 

%- dynamic fields and absorbing boundary field ----------------------------
v = zeros(nx,nz);
sxy = zeros(nx-1,nz);
szy = zeros(nx,nz-1);


%- initialise seismograms -------------------------------------------------
displacement_seismograms = zeros(n_receivers,nt);


%- initialise absorbing boundary taper a la Cerjan ------------------------
[absbound] = init_absbound();


%==========================================================================
% iterate
%==========================================================================

%%% TEST TIME DOMAIN VERSION %%%
% G_2 = flip( G_2 , 3 );
% S = zeros(nx,nz,length(t));
% for n = 1:length(t)
%     for ns = 1:n_noise_sources
%         
%         for k = 1:n_sample
%             S(:,:,n) = S(:,:,n) + noise_spectrum(k,ns) * noise_source_distribution(:,:,ns) * ifft_coeff(n,k);
%         end
%         
%     end
% end
% S = real(S);
% 
% for i = 1:size(G_2,1)
%     for j = 1:size(G_2,2)
%         newf(i,j,:) = conv( squeeze(G_2(i,j,:)), squeeze(S(i,j,:)) );
%     end
% end

u = zeros(nx,nz);

for n = 1:length(t)
    
    %- compute divergence of current stress tensor ------------------------    
    DS = div_s(sxy,szy,dx,dz,nx,nz,order);
    
    
    %- add sources of the correlation field -------------------------------    
    if( mod(n,freq_samp)==0 && t(n)<=0.0 )       
        
        %- transform on the fly to the time domain        
        S = zeros(nx,nz,n_noise_sources) + 1i*zeros(nx,nz,n_noise_sources);
        
        
        % % noise source for coupled spectra and distributions
        % for ns = 1:n_noise_sources
        %     
        %     for k=1:n_sample
        %         S(:,:,ns) = S(:,:,ns) + noise_spectrum(k,ns) * G_2(:,:,k) * ifft_coeff(n,k);
        %     end
        %     
        %     DS = DS + noise_source_distribution(:,:,ns) .* real(S(:,:,ns));
        %     
        % end
             
        
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

    
    %%% TEST TIME DOMAIN VERSION %%%
    % DS = DS + newf(:,:,n);
    
    
    %- possible perturbation of rhs to test adjoint state -----------------
    if( size(df,1)==1 && size(df,2)==1 )
        DS = real( DS + repmat(df,nx,nz) );
    else
        DS = real( DS + df(:,:,n) );
    end
    
    
    %- update velocity field ----------------------------------------------
    v = v + dt * DS./rho;
    
    
    %- apply absorbing boundary taper -------------------------------------    
    v = v .* absbound;
    
    
    %- compute derivatives of current velocity and update stress tensor ---
    strain_dxv = dx_v(v,dx,dz,nx,nz,order);
    strain_dzv = dz_v(v,dx,dz,nx,nz,order);
    
    sxy = sxy + dt * mu(1:nx-1,:) .* strain_dxv;
    szy = szy + dt * mu(:,1:nz-1) .* strain_dzv;
    
    
    %- calculate displacement and respective strain -----------------------
    u = u + v * dt;
    strain_dxu = dx_v(u,dx,dz,nx,nz,order);
    strain_dzu = dz_v(u,dx,dz,nx,nz,order);
    
    
    %- record velocity seismograms ----------------------------------------    
    for k = 1:n_receivers
        displacement_seismograms(k,n) = u(rec_id(k,1),rec_id(k,2));
    end
    

    %- accumulate Fourier transform of the correlation displacement field -    
    if( mode==1 && mod(n,freq_samp)==0 )
        
        for k = 1:n_sample            
            C_2_dxu(:,:,k) = C_2_dxu(:,:,k) + strain_dxu * fft_coeff(n,k);
            C_2_dzu(:,:,k) = C_2_dzu(:,:,k) + strain_dzu * fft_coeff(n,k);
        end

    end
    
    
    %%% TEST TIME DOMAIN VERSION %%%
    % C_2_dxu_time(:,:,n) = strain_dxu;
    % C_2_dzu_time(:,:,n) = strain_dzu;
    
end


% time = datetime('now');
% formatOut = 'yyyy-mm-dd-HH-MM-SS';
% save(sprintf('stress_correlation_%s.mat',datestr(time,formatOut,'local')), 'C_2_dxu', 'C_2_dzu')


end

