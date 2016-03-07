function [G_2, G_2_dxu_time, G_2_dzu_time] = run_forward1_green( mu, rho, src, mode )

%==========================================================================
% compute Green function for reference station
%
% input:
%--------
% mu [N/m^2]
% rho [kg/m^3]
% src: source position, i.e. the reference station
% mode: integer switch, mode==0 when forward strain is not needed
%
% output:
%--------
% G_2: Fourier transformed displacement Green function of reference station
% G__dxu_time & G__dzu_time: strain of displacemen Green function
%
%==========================================================================


%==========================================================================
% initialise simulation
%==========================================================================

%- material and domain ----------------------------------------------------
[Lx,Lz,nx,nz,dt,nt,order,~,~,~,fw_nth] = input_parameters();
[~,~,x,z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);
rho = reshape(rho,nx,nz);
mu = reshape(mu,nx,nz);


%- initialise interferometry ----------------------------------------------
[~,n_sample,w_sample,~,freq_samp] = input_interferometry();


%- time axis --------------------------------------------------------------
t = 0:dt:(nt-1)*dt;


%- compute indices for source locations -----------------------------------
ns = size(src,1);
src_id = zeros(ns,2);
for i = 1:ns
    src_id(i,1) = min( find( min(abs(x-src(i,1))) == abs(x-src(i,1)) ) );
    src_id(i,2) = min( find( min(abs(z-src(i,2))) == abs(z-src(i,2)) ) );
end


%- make source time function ----------------------------------------------
stf = 1.0e9*ones(1,nt);


%- Fourier transform of the forward Greens function -----------------------
G_2 = zeros(nx,nz,n_sample) + 1i*zeros(nx,nz,n_sample);


%- prepare coefficients for Fourier transform -----------------------------
n_ftc = 0;
for n = nt:(2*nt-1)
    if( mod(n,freq_samp) == 0 )
        n_ftc = n_ftc + 1;
    end
end

fft_coeff = zeros(n_ftc,n_sample) + 1i*zeros(n_ftc,n_sample);
i_ftc = 1;
for n = nt:(2*nt-1)
    if( mod(n,freq_samp) == 0 )
        
        for k = 1:n_sample
            fft_coeff(i_ftc,k) = 1/sqrt(2*pi) * exp(-1i*w_sample(k)*t(n-nt+1)) * dt;
        end
        i_ftc = i_ftc + 1;
        
    end
end


%- allocate dynamic fields ------------------------------------------------
v = zeros(nx,nz);
sxy = zeros(nx-1,nz);
szy = zeros(nx,nz-1);

n_fw = 0;
for n = nt:(2*nt-1)
    if( mod(n,fw_nth) == 0 )
        n_fw = n_fw + 1;
    end
end

if( mode ~= 0 )
    G_2_dxu_time = zeros(nx-1,nz,n_fw,'single');
    G_2_dzu_time = zeros(nx,nz-1,n_fw,'single');
else
    G_2_dxu_time = single(0.0);
    G_2_dzu_time = single(0.0);
end
strain_dxv = zeros(nx-1,nz);
strain_dzv = zeros(nx,nz-1);


%- initialise absorbing boundary taper a la Cerjan ------------------------
[absbound] = init_absbound();


%==========================================================================
% iterate
%==========================================================================

i_ftc = 1;
i_fw = 1;
for n = 1:nt
    
    
    if( mode ~= 0 && mod(n+nt-1, fw_nth) == 0 )
        G_2_dxu_time(:,:,i_fw) = single(strain_dxv);
        G_2_dzu_time(:,:,i_fw) = single(strain_dzv);
        
        i_fw = i_fw + 1;
    end
    
    
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
    if( mod(n+nt-1, freq_samp) == 0 )         
        
        for k = 1:n_sample
            G_2(:,:,k) = G_2(:,:,k) + v * fft_coeff(i_ftc,k);
        end
        i_ftc = i_ftc + 1;
        
    end

    
end


%% return the time reversed Green function
for k = 1:n_sample
    G_2(:,:,k) = conj(G_2(:,:,k));
end


end

