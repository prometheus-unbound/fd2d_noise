function [K, stf_fft] = run_noise_adjoint( mu, rho, forward_dxu_time, forward_dzu_time, adstf, adsrc, spectrum, source_dist, G_2, mode )

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
% G_2: Fourier transformed displacement Green function of reference station
% mode: integer switch
%       == 0 for source inversion
%       == 1 for structure inversion
%       == 2 for joint inversion
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


%- time axis --------------------------------------------------------------
t = -(nt-1)*dt:dt:(nt-1)*dt;
n_zero = nt;
nt = length(t);


%- initialise interferometry ----------------------------------------------
[~,n_sample,w_sample,dw,freq_samp] = input_interferometry();


%- compute indices for adjoint source locations ---------------------------    
ns_adj = size(adsrc,1);
adsrc_id = zeros(ns_adj,2);

for i=1:ns_adj
    adsrc_id(i,1) = min( find( min(abs(x-adsrc(i,1))) == abs(x-adsrc(i,1))) );
    adsrc_id(i,2) = min( find( min(abs(z-adsrc(i,2))) == abs(z-adsrc(i,2))) );
end


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

%- prepare coefficients for Fourier transform and its inverse -------------
ifft_coeff2 = zeros(nt,n_sample) + 1i*zeros(nt,n_sample);
for k = 1:n_sample
    ifft_coeff2(:,k) = 1/sqrt(2*pi) * exp( 1i*w_sample(k)*t' ) * dw;
end


%- get integration boundaries for frequency bands -------------------------
if( n_basis_fct ~= 0 )
    [int_limits] = integration_limits(n_sample,n_basis_fct);
end


%- initialise absorbing boundary taper a la Cerjan ------------------------
[absbound] = init_absbound();


%- reshape source distribution --------------------------------------------
n_noise_sources = size(spectrum,2);

if( n_basis_fct == 0 )
    noise_source_distribution = reshape(source_dist, nx, nz, n_noise_sources);
else
    noise_source_distribution = reshape(source_dist, nx, nz, n_basis_fct);
end


%- allocate dynamic fields ------------------------------------------------
v = zeros(nx,nz);
sxy = zeros(nx-1,nz);
szy = zeros(nx,nz-1);
u = zeros(nx,nz);
stf_fft = zeros(nx,nz,n_sample) + 1i*zeros(nx,nz,n_sample);


%- allocate source kernel structure ---------------------------------------
K_mu = zeros(nx,nz);

if( n_basis_fct == 0 )
    K_s = zeros(nx,nz,1);
else
    K_s = zeros(nx,nz,n_basis_fct);
end


%==========================================================================
% iterate
%==========================================================================

%- only need first half of second adjoint run
% strain of Green function is zero in second half, i.e. with time reversal
% only first half of second adjoint run is necessary
if( size(adstf,3) ~= 1 )
    nt = n_zero;
end

i_ftc = 1;
i_fw = 1;

for n = 1:nt
    
    
    %- compute divergence of current stress tensor ------------------------    
    DS = div_s(sxy,szy,dx,dz,nx,nz,order);
    
    
    %- add adjoint source time function -----------------------------------
    % size(adstf,3) == 1 means: first run
    if( size(adstf,3) == 1 )
        
        for i=1:ns_adj
            DS(adsrc_id(i,1),adsrc_id(i,2)) = DS(adsrc_id(i,1),adsrc_id(i,2)) + real(adstf(i,n));
        end

    % second run, does not happen with mode == 0
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
    
    
    %- build up source kernel ---------------------------------------------
    % mode == 1 means inversion for structure only, don't need source kernel
    if( mode ~= 1 && size(adstf,3) == 1 && mod(n,freq_samp) == 0 && t(end-n+1) <= 0.0 )
        
        
        if( n_basis_fct == 0 )
            M_tn = zeros(nx,nz,1) + 1i*zeros(nx,nz,1);
        else
            M_tn = zeros(nx,nz,n_basis_fct) + 1i*zeros(nx,nz,n_basis_fct);
        end
        
        
        for k = 1:n_sample
            
            % funny construction with tmp variable needed for mex generation
            if( n_basis_fct == 0 )
                ib = 1;
                tmp = repmat( spectrum(k) * G_2(:,:,k) * ifft_coeff2(end-n+1, k), 1, 1, 1);
            else
                ib = find( k >= int_limits(:,1) & k <= int_limits(:,2) );
                tmp = repmat( spectrum(k) * G_2(:,:,k) * ifft_coeff2(end-n+1, k), 1, 1, n_basis_fct);
            end           
            
            M_tn(:,:,ib) = M_tn(:,:,ib) + tmp(:,:,ib);
            
        end
        
        
        for ib = 1:size(M_tn,3)
            K_s(:,:,ib) = K_s(:,:,ib) + real( M_tn(:,:,ib) .* u );
        end
        
        
    end
    
       
    %- build up structure kernel ------------------------------------------
    % mode == 0 means inversion for source only, don't need structure kernel
    if( mode ~= 0 && mod(n, fw_nth) == 0 )
        K_mu(1:nx-1,:) = K_mu(1:nx-1,:) - strain_dxu .* forward_dxu_time(:,:,end-i_fw+1) * fw_nth;
        K_mu(:,1:nz-1) = K_mu(:,1:nz-1) - strain_dzu .* forward_dzu_time(:,:,end-i_fw+1) * fw_nth;
        
        i_fw = i_fw + 1;
    end


    %- save fourier transformed adjoint state for second run --------------
    % mode == 0 means inversion for source only, don't need second run
    if( mode ~= 0 && size(adstf,3) == 1 )
        
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


%- concatenate kernels ----------------------------------------------------
if( n_basis_fct == 0 )
    K = zeros(nx, nz, 2);
    K(:,:,1) = K_s;
else
    K = zeros(nx, nz, n_basis_fct + 1);
    K(:,:,1:n_basis_fct) = K_s;
end

K(:,:,end) = K_mu;

