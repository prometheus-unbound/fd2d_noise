function [K, stf_fft] = run3_adjoint( structure, noise_source, G_fft, adjsrc, adjstf, wavefield_fwd, mode )

%==========================================================================
% compute sensitivity kernel for mu (and rho)
%
% input:
%--------
% mu [N/m^2]
% rho [kg/m^3]
% u_fwd: forward wavefield
% stf: adjoint source time functions
% src: adjoint source locations
% spectrum: spectrum of noise distribution
% source_dist: source distribution
% G_2: Fourier transformed displacement Green function of reference station
% mode: integer switch
%       == 0 gradient mode: do not Fourier transform adjoint wavefield
%       == 1 gradient mode: Fourier transform adjoint wavefield
%
% output:
%--------
% K_source: sensitivity kernel for the source distribution
% K_mu: sensitivity kernel for mu
% stf_fft: Fourier transformed adjoint_state - necessary for second run
%
%==========================================================================


%- basic configuration ----------------------------------------------------
[Lx, Lz, nx, nz, dt, nt, order, ~, ~, store_fwd_nth] = input_parameters();
[~,~,x,z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);


%- time axis --------------------------------------------------------------
t = -(nt-1)*dt:dt:(nt-1)*dt;
n_zero = nt;
nt = length(t);


%- prepare coefficients for Fourier transform -----------------------------
[~,n_sample,w_sample,dw,freq_samp] = input_interferometry();

fft_coeff = zeros(nt,n_sample) + 1i*zeros(nt,n_sample);
ifft_coeff = zeros(nt,n_sample) + 1i*zeros(nt,n_sample);
for n=1:nt
    for k = 1:n_sample
        fft_coeff(n,k) =  1/sqrt(2*pi) * exp(-1i*w_sample(k)*t(n) ) * dt;
        ifft_coeff(n,k) = 1/sqrt(2*pi) * exp( 1i*w_sample(k)*t(n) ) * dw;
    end
end


%- compute indices for adjoint source locations ---------------------------    
ns_adj = size(adjsrc,1);
src_id = zeros(ns_adj,2);

for i=1:ns_adj
    src_id(i,1) = min( find( min(abs(x-adjsrc(i,1))) == abs(x-adjsrc(i,1))) );
    src_id(i,2) = min( find( min(abs(z-adjsrc(i,2))) == abs(z-adjsrc(i,2))) );
end


%- initialise absorbing boundary taper a la Cerjan ------------------------
[absbound] = init_absbound();


%- allocate dynamic fields ------------------------------------------------
v = zeros(nx,nz);
sxy = zeros(nx-1,nz);
szy = zeros(nx,nz-1);
u = zeros(nx,nz);

if( mode == 1 )
    stf_fft = zeros(nx,nz,n_sample) + 1i*zeros(nx,nz,n_sample);
else
    stf_fft = [];
end


%- allocate source kernel structure ---------------------------------------
K_mu = zeros(nx,nz);
K_source = zeros(nx,nz);


%==========================================================================
% iterate
%==========================================================================

%- only need first half of second adjoint run
if( ( mode == 0 || mode == 1 ) && size(adjstf,3) ~= 1 )
    nt = n_zero;
end

i_fw_in = 1;

for n = 1:nt
    
    
    %- compute divergence of current stress tensor ------------------------    
    DS = div_s(sxy,szy,dx,dz,nx,nz,order); 
    
    
    %- add adjoint source time function -----------------------------------
    if( size(adjstf,3) == 1 && ~isempty(adjstf) )
        
        for i=1:ns_adj
            DS(src_id(i,1),src_id(i,2)) = DS(src_id(i,1),src_id(i,2)) + real(adjstf(i,n));
        end
        
    elseif( size(adjstf,3) ~= 1 && ~isempty(adjstf) )
        
        if( mod(n,freq_samp) == 0 )
            T = zeros(nx,nz) + 1i*zeros(nx,nz);
            
            for k = 1:n_sample
                T = T + noise_source.spectrum(k) * noise_source.distribution .* conj( adjstf(:,:,k) ) * ifft_coeff(n,k);
            end
            
            DS = DS + real(T);
            
        end
        
    end
    
    
    %- update velocity field ----------------------------------------------
    v = v + dt * DS./structure.rho;
    
    
    %- apply absorbing boundary taper -------------------------------------    
    v = v .* absbound;
    
    
    %- compute derivatives of current velocity and update stress tensor ---
    strain_dxv = dx_v(v,dx,dz,nx,nz,order);
    strain_dzv = dz_v(v,dx,dz,nx,nz,order);
    
    sxy = sxy + dt * structure.mu(1:nx-1,:) .* strain_dxv;
    szy = szy + dt * structure.mu(:,1:nz-1) .* strain_dzv;
    
    
    %- calculate displacement and respective strain -----------------------
    u = u + v * dt;
    strain_dxu = dx_v(u,dx,dz,nx,nz,order);
    strain_dzu = dz_v(u,dx,dz,nx,nz,order);  
    
    
    %- build up source kernel ---------------------------------------------
    if( ~isempty(G_fft) && size( G_fft, 3 ) == n_sample && mod(n,freq_samp) == 0 )
        
        M_tn = zeros(nx,nz) + 1i*zeros(nx,nz);
        
        for k = 1:n_sample
            M_tn = M_tn + noise_source.spectrum(k) .* G_fft(:,:,k) * ifft_coeff(n,k);
        end
        
        for ib = 1:size(M_tn,3)
            K_source(:,:,ib) = K_source(:,:,ib) + real( M_tn(:,:,ib) .* u );
        end
        
    end
    
       
    %- build up structure kernel ------------------------------------------
    if( ~isempty( wavefield_fwd ) && size( wavefield_fwd, 3 ) >= i_fw_in && mod( n, store_fwd_nth ) == 0 )
                    
        K_mu(1:nx-1,:) = K_mu(1:nx-1,:) - strain_dxu .* dx_v( wavefield_fwd(:,:,end-i_fw_in+1), dx, dz, nx, nz, order ) * store_fwd_nth;
        K_mu(:,1:nz-1) = K_mu(:,1:nz-1) - strain_dzu .* dz_v( wavefield_fwd(:,:,end-i_fw_in+1), dx, dz, nx, nz, order ) * store_fwd_nth;
        
        i_fw_in = i_fw_in + 1;

    end


    %- save Fourier transformed adjoint state for second run --------------
    if( mode == 1 )
        
        if( mod(n,freq_samp) == 0 )
            for k = 1:n_sample
                stf_fft(:,:,k) = stf_fft(:,:,k) + u * fft_coeff(n,k);
            end
        end
        
    end
    
    
end


%- concatenate kernels ----------------------------------------------------
K = zeros(nx, nz, 2);
K(:,:,1) = K_mu;
K(:,:,2) = K_source;


end
