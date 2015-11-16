function [K_s] = run_noise_source_kernel_fast( G_2, mu, spectrum, stf, adsrc )

%==========================================================================
% run simulation to compute sensitivity kernel for noise power-spectral
% density distribution
% fast means ready for conversion to mex-files
%
% input:
%--------
% G_2: Green function of reference station
% mu [N/m^2]
% stf: adjoint source time function
% adsrc: adjoint source positions
%
% output:
%--------
% K_s: sensitivity kernel
%
%==========================================================================


%==========================================================================
% initialise simulation
%==========================================================================

%- material and domain ----------------------------------------------------
[Lx,Lz,nx,nz,dt,nt,order,model_type,~,n_basis_fct] = input_parameters();
[~,~,x,z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);
[~,rho] = define_material_parameters(nx,nz,model_type); 
mu = reshape(mu,nx,nz);


%- time axis --------------------------------------------------------------    
t = -(nt-1)*dt:dt:(nt-1)*dt;
    

%- compute indices for adjoint source locations ---------------------------
ns = size(adsrc,1);
adsrc_id = zeros(ns,2);

for i=1:ns
    adsrc_id(i,1) = min(find(min(abs(x-adsrc(i,1)))==abs(x-adsrc(i,1))));
    adsrc_id(i,2) = min(find(min(abs(z-adsrc(i,2)))==abs(z-adsrc(i,2))));
end


%- initialise interferometry ----------------------------------------------       
[~,n_sample,w_sample,dw,freq_samp] = input_interferometry();

if( n_basis_fct == 0 )
    n_basis_fct = 1;
end
K_s = zeros(nx,nz,n_basis_fct);


ifft_coeff = zeros(length(t),n_sample) + 1i*zeros(length(t),n_sample);
for k = 1:n_sample
    ifft_coeff(:,k) = 1/sqrt(2*pi) * exp( 1i*w_sample(k)*t' ) * dw;
end


%- get integration boundaries for frequency bands -------------------------
% if( n_basis_fct ~= 0 )
    [int_limits] = integration_limits(n_sample,n_basis_fct);
% end


%- dynamic fields and absorbing boundary field ----------------------------
v = zeros(nx,nz);
sxy = zeros(nx-1,nz);
szy = zeros(nx,nz-1);


%- initialise absorbing boundary taper a la Cerjan ------------------------
[absbound] = init_absbound();


%==========================================================================
% iterate
%==========================================================================

u = zeros(nx,nz);

for n=1:length(t)
    
    %- compute divergence of current stress tensor ------------------------    
    DS = div_s(sxy,szy,dx,dz,nx,nz,order);
    
    
    %- add point sources --------------------------------------------------    
    for i=1:ns
        DS(adsrc_id(i,1),adsrc_id(i,2)) = DS(adsrc_id(i,1),adsrc_id(i,2)) + stf(i,n);
    end
    
    
    %- update velocity field ----------------------------------------------    
    v = v + dt*DS./rho;
    
    
    %- apply absorbing boundary taper -------------------------------------    
    v = v .* absbound;
    
    
    %- compute derivatives of current velocity and update stress tensor ---    
    sxy = sxy + dt * mu(1:nx-1,:) .* dx_v(v,dx,dz,nx,nz,order);
    szy = szy + dt * mu(:,1:nz-1) .* dz_v(v,dx,dz,nx,nz,order);
     
    
    %- accumulate kernel --------------------------------------------------
    u = u + v;
    

    M_tn = zeros(nx,nz,n_basis_fct) + 1i*zeros(nx,nz,n_basis_fct);
    tmp = zeros(nx,nz) + 1i*zeros(nx,nz);
    if( mod(n,freq_samp)==0 && t(end-n+1)<=0.0 )
        
        for k = 1:n_sample
            
            if( n_basis_fct ~= 0 )
                ib = find( k >= int_limits(:,1) & k <= int_limits(:,2) );
            else
                ib = 1;
            end
            
            % funny construction with tmp variable needed for mex generation
            tmp = repmat( spectrum(k) * conj(G_2(:,:,k)) * ifft_coeff(end-n+1, k), 1, 1, n_basis_fct);
            M_tn(:,:,ib) = M_tn(:,:,ib) + tmp(:,:,ib);
            
        end
        
        for ib = 1:n_basis_fct
            K_s(:,:,ib) = K_s(:,:,ib) + real( M_tn(:,:,ib) .* u * dt );
        end
        
    end
    
end

