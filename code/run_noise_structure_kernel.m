function [K_mu, stf_fft] = run_noise_structure_kernel( mu, rho, forward_dxu_time, forward_dzu_time, adstf, adsrc, spectrum, source_dist )

%==========================================================================
% compute sensitivity kernel for mu (and rho)
%
% input:
%--------
% mu [N/m^2]
% rho [kg/m^3]
% C_2_dxu_time & C_2_dzu_time: strain of correlation wavefield
% adstf: adjoint source time functions
% adsrc: adjoint source locations
% spectrum: spectrum of noise distribution
% source_dist: source distribution
%
% output:
%--------
% K_mu: sensitivity kernel for mu
% adjoint_state: necessary for second run
%
%==========================================================================


%==========================================================================
% initialise simulation
%==========================================================================

%- material and domain ----------------------------------------------------
[Lx,Lz,nx,nz,dt,nt,order,~,~,n_basis_fct,fw_nth] = input_parameters();
[~,~,x,z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);
rho = reshape(rho, nx, nz);
mu = reshape(mu, nx, nz);


%- initialise interferometry ----------------------------------------------
[~,n_sample,w_sample,dw,freq_samp] = input_interferometry();


%- reshape source distribution --------------------------------------------
n_noise_sources = size(spectrum,2);

if( n_basis_fct == 0 )
    noise_source_distribution = reshape(source_dist, nx, nz, n_noise_sources);
else
    noise_source_distribution = reshape(source_dist, nx, nz, n_basis_fct);
end


%- get integration boundaries for frequency bands -------------------------
if( n_basis_fct ~= 0 )
    [int_limits] = integration_limits(n_sample,n_basis_fct);
end


%- compute indices for adjoint source locations ---------------------------    
ns_adj = size(adsrc,1);
adsrc_id = zeros(ns_adj,2);

for i=1:ns_adj
    adsrc_id(i,1) = min( find( min(abs(x-adsrc(i,1))) == abs(x-adsrc(i,1))) );
    adsrc_id(i,2) = min( find( min(abs(z-adsrc(i,2))) == abs(z-adsrc(i,2))) );
end


%- time axis --------------------------------------------------------------
t = -(nt-1)*dt:dt:(nt-1)*dt;
n_zero = nt;
nt = length(t);


%- prepare coefficients for Fourier transform and its inverse -------------
n_ftc = floor(nt/freq_samp);
fft_coeff = zeros(n_ftc,n_sample) + 1i*zeros(n_ftc,n_sample);
ifft_coeff = zeros(n_ftc,n_sample) + 1i*zeros(n_ftc,n_sample);

i_ftc = 1;
for n=1:nt
    if( mod(n,freq_samp) == 0 )
        for k = 1:n_sample
            fft_coeff(i_ftc,k) =  1/sqrt(2*pi) * exp(-1i*w_sample(k)*t(n) ) * dt;
            ifft_coeff(i_ftc,k) = 1/sqrt(2*pi) * exp( 1i*w_sample(k)*t(n) ) * dw;
        end
        
        i_ftc = i_ftc + 1;
    end
end


%- allocate dynamic fields ------------------------------------------------
v = zeros(nx,nz);
sxy = zeros(nx-1,nz);
szy = zeros(nx,nz-1);
u = zeros(nx,nz);
stf_fft = zeros(nx,nz,n_sample) + 1i*zeros(nx,nz,n_sample);
K_mu = zeros(nx,nz);


%- initialise absorbing boundary taper a la Cerjan ------------------------
[absbound] = init_absbound();


%==========================================================================
% iterate
%==========================================================================

%- only need first half of second adjoint run - strain of Green function is zero in second half
if( size(adstf,3) ~= 1 )
    nt = n_zero;
end

i_ftc = 1;
i_fw = 1;
for n = 1:nt
    
    
    %- compute divergence of current stress tensor ------------------------    
    DS = div_s(sxy,szy,dx,dz,nx,nz,order);
    
    
    %- add adjoint source time function -----------------------------------
    if( size(adstf,3) == 1 )
        
        for i=1:ns_adj
            DS(adsrc_id(i,1),adsrc_id(i,2)) = DS(adsrc_id(i,1),adsrc_id(i,2)) + real(adstf(i,n));
        end
        
    else
        
        if( mod(n,freq_samp) == 0 && t(n) < 0.0 )
            
            T = zeros(nx,nz) + 1i*zeros(nx,nz);
            for ns = 1:n_noise_sources
                
                if( n_basis_fct ~= 0 )
                    ib = find( k >= int_limits(:,1) & k <= int_limits(:,2) );
                else
                    ib = ns;
                end
                
                for k = 1:n_sample                   
                    T = T + spectrum(k,ns) * noise_source_distribution(:,:,ib) .* conj(adstf(:,:,k)) * ifft_coeff(i_ftc,k);                    
                end
                
            end
            
            i_ftc = i_ftc + 1;
            DS = DS + real(T);
            
        end
        
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
    
       
    %- build up kernel
    if( mod(n, fw_nth) == 0 )
        K_mu(1:nx-1,:) = K_mu(1:nx-1,:) - strain_dxu .* forward_dxu_time(:,:,end-i_fw+1) * fw_nth;
        K_mu(:,1:nz-1) = K_mu(:,1:nz-1) - strain_dzu .* forward_dzu_time(:,:,end-i_fw+1) * fw_nth;
        
        i_fw = i_fw + 1;
    end


    %- save fourier transformed adjoint state for second run --------------
    if( size(adstf,3) == 1 )
        
        if( mod(n,freq_samp) ~= 0 )
            continue
        elseif( mod(n,freq_samp) == 0 && t(n) < 0.0)
            i_ftc = i_ftc + 1;
            continue
        end
        
        for k = 1:n_sample
            stf_fft(:,:,k) = stf_fft(:,:,k) + u * fft_coeff(i_ftc,k);
        end
        i_ftc = i_ftc + 1;
        
    end   
    
    
end

