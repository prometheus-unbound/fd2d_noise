function [ G_fft, G_out ] = run_forward1_green( structure, src, mode )

%==========================================================================
% compute Green function for reference station
%
% input:
%--------
% structure: contains mu [N/m^2] and rho [kg/m^3]
% src: source position, i.e. the reference station
% mode: integer switch 
%       == 0 do not save wavefield
%       == 1 save wavefield
%
% output:
%--------
% G_fft: Fourier transformed displacement Green function of reference station
% G_out: Green function wavefield
%
%==========================================================================


%- basic configuration ----------------------------------------------------
[Lx, Lz, nx, nz, dt, nt, order, ~, ~, store_fwd_nth] = input_parameters();
[~,~,x,z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);


%- time axis --------------------------------------------------------------
t = 0:dt:(nt-1)*dt;


%- prepare coefficients for Fourier transform -----------------------------
[~,n_sample,w_sample,~,freq_samp] = input_interferometry();

fft_coeff = zeros(nt,n_sample) + 1i*zeros(nt,n_sample);
for n=1:nt
    for k = 1:n_sample
        fft_coeff(n,k) =  1/sqrt(2*pi) * exp(-1i*w_sample(k)*t(n) ) * dt;
    end
end


%- make source time function ----------------------------------------------
stf = 1.0e9*ones(1,nt);


%- compute indices for source locations -----------------------------------
ns = size(src,1);
src_id = zeros(ns,2);
for i = 1:ns
    src_id(i,1) = min( find( min(abs(x-src(i,1))) == abs(x-src(i,1)) ) );
    src_id(i,2) = min( find( min(abs(z-src(i,2))) == abs(z-src(i,2)) ) );
end


%- initialise absorbing boundary taper a la Cerjan ------------------------
[absbound] = init_absbound();


%- allocate dynamic fields ------------------------------------------------
v = zeros(nx,nz);
sxy = zeros(nx-1,nz);
szy = zeros(nx,nz-1);
G_fft = zeros(nx,nz,n_sample) + 1i*zeros(nx,nz,n_sample);

n_fwd = 0;
for n = nt:(2*nt-1)
    if( mod(n,store_fwd_nth) == 0 )
        n_fwd = n_fwd + 1;
    end
end

if( mode ~= 0 )
    G_out = zeros(nx,nz,n_fwd,'single');
    % G_out = zeros(nx,nz,n_fw);
else
    G_out = [];
end


%==========================================================================
% iterate
%==========================================================================

i_fwd_out = 1;

for n = 1:nt
    
    
    if( mode ~= 0 && mod(n+nt-1, store_fwd_nth) == 0 )
        G_out(:,:,i_fwd_out) = single( v );
        % G_out(:,:,i_fw_out) = v;
        i_fwd_out = i_fwd_out + 1;
    end
        
    
    %- compute divergence of current stress tensor ------------------------    
    DS = div_s(sxy,szy,dx,dz,nx,nz,order); 
    
    
    %- add point sources --------------------------------------------------        
    for i=1:ns
        DS(src_id(i,1),src_id(i,2)) = DS(src_id(i,1),src_id(i,2)) + stf(n);
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
    
    
    %- accumulate Fourier transform of the displacement Greens function ---
    if( mod(n+nt-1, freq_samp) == 0 )         
        
        for k = 1:n_sample
            G_fft(:,:,k) = G_fft(:,:,k) + v * fft_coeff(n,k);
        end
        
    end 
        
    
end


end

