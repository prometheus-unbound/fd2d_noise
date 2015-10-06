function [G_2] = run_forward_green_fast(mu,src)

%==========================================================================
% run forward simulation
% fast means ready for conversion to mex-files
%
% input:
%--------
% mu [N/m^2]
% src: source, i.e. reference station
%
% output:
%--------
% G_2: displacement Green function for reference station
%
%==========================================================================


%==========================================================================
% initialise simulation
%==========================================================================

%- material and domain ----------------------------------------------------
[Lx,Lz,nx,nz,dt,nt,order,model_type] = input_parameters();
[~,~,x,z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);
[~,rho] = define_material_parameters(nx,nz,model_type);
mu = reshape(mu,nx,nz);


%- initialise interferometry ----------------------------------------------
[~,n_sample,w_sample,~,nt_freq] = input_interferometry();


%- time axis --------------------------------------------------------------
t = 0:dt:dt*(nt-1);


%- compute indices for source locations -----------------------------------
ns = size(src,1);
src_id = zeros(ns,2);
for i = 1:ns
    src_id(i,1) = min( find( min(abs(x-src(i,1))) == abs(x-src(i,1)) ) );
    src_id(i,2) = min( find( min(abs(z-src(i,2))) == abs(z-src(i,2)) ) );
end


%- make source time function ----------------------------------------------
stf = 1.0e9*ones(1,length(t));


%- Fourier transform of the forward Greens function
G_2 = zeros(nx,nz,n_sample) + 1i*zeros(nx,nz,n_sample);


%- prepare coefficients for Fourier transform -----------------------------
fft_coeff = zeros(length(t),n_sample) + 1i*zeros(length(t),n_sample);
for k = 1:n_sample
    fft_coeff(:,k) = 1/sqrt(2*pi) * exp(-1i*w_sample(k)*t') * dt;
end


%- dynamic fields and absorbing boundary field ----------------------------
v = zeros(nx,nz);
sxy = zeros(nx-1,nz);
szy = zeros(nx,nz-1);


%- initialise absorbing boundary taper a la Cerjan ------------------------
[absbound] = init_absbound();


%==========================================================================
% iterate
%==========================================================================

%%% TEST TIME DOMAIN VERSION %%%
% G_2_test = zeros(nx,nz,2*nt-1);

for n = 1:length(t)
    
    %- compute divergence of current stress tensor ------------------------    
    DS = div_s(sxy,szy,dx,dz,nx,nz,order);
    
    
    %- add point sources --------------------------------------------------    
    for i=1:ns
        DS(src_id(i,1),src_id(i,2)) = DS(src_id(i,1),src_id(i,2)) + stf(n);
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
    
    
    %- accumulate Fourier transform of the displacement Greens function ---
    if( mod(n,nt_freq) == 0 )         
        
        for k = 1:n_sample
            G_2(:,:,k) = G_2(:,:,k) + v * fft_coeff(n,k);
        end
        
    end
    
    %%% TEST TIME DOMAIN VERSION %%%
    % G_2_test(:,:,nt-1+n) = v;
    
end


end

