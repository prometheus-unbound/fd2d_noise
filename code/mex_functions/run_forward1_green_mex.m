function [ G_fft, G_out, seismograms ] = run_forward1_green( mu, rho, src, rec, mode, dmu, G_in )

%==========================================================================
% compute Green function for reference station
%
% input:
%--------
% mu [N/m^2]
% rho [kg/m^3]
% src: source position, i.e. the reference station
% rec: receiver
% mode: integer switch 
%       == 0 do not save wavefield
%       == 1 save wavefield
%
% output:
%--------
% G_fft: Fourier transformed displacement Green function of reference station
% G_out: Green function wavefield
% displacement seismograms
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


%- compute indices for receiver locations ---------------------------------
n_receivers = size(rec,1);
rec_id = zeros(n_receivers,2);

for i=1:n_receivers    
    rec_id(i,1) = min( find( min(abs(x-rec(i,1))) == abs(x-rec(i,1)) ) );
    rec_id(i,2) = min( find( min(abs(z-rec(i,2))) == abs(z-rec(i,2)) ) );   
end


%- make source time function ----------------------------------------------
stf = 1.0e9*ones(1,nt);


%- Fourier transform of the forward Greens function -----------------------
G_fft = zeros(nx,nz,n_sample) + 1i*zeros(nx,nz,n_sample);


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
    G_out = zeros(nx,nz,n_fw,'single');
    % G_out = zeros(nx,nz,n_fw);
else
    G_out = single(0.0);
    % G_out = 0.0;
end

% G_out_2 = zeros(nx,nz,2*nt-1,'single');


%- initialise seismograms -------------------------------------------------
seismograms = zeros(n_receivers,nt);


%- initialise absorbing boundary taper a la Cerjan ------------------------
[absbound] = init_absbound();


%==========================================================================
% iterate
%==========================================================================

i_ftc = 1;
i_fw_out = 1;
i_fw_in = 1;

for n = 1:nt
    
    
    if( mode ~= 0 && mod(n+nt-1, fw_nth) == 0 )
        G_out(:,:,i_fw_out) = single( v );
        % G_out(:,:,i_fw_out) = v;
        i_fw_out = i_fw_out + 1;
    end
    
    % G_out_2(:,:,nt-1+n) = single( v );
    
    
    if( ~isempty( dmu ) )
        sxy = sxy + dt * dmu(1:nx-1,:) .* dx_v( G_in(:,:,i_fw_in), dx, dz, nx, nz, order );
        szy = szy + dt * dmu(:,1:nz-1) .* dz_v( G_in(:,:,i_fw_in), dx, dz, nx, nz, order );
        i_fw_in = i_fw_in + 1;
    end
    
    
    %- compute divergence of current stress tensor ------------------------    
    DS = div_s(sxy,szy,dx,dz,nx,nz,order); 
    
    
    %- add point sources --------------------------------------------------    
    if( isempty( dmu ) )
        for i=1:ns
            DS(src_id(i,1),src_id(i,2)) = DS(src_id(i,1),src_id(i,2)) + stf(n);
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
    
    
    %- accumulate Fourier transform of the displacement Greens function ---
    if( mod(n+nt-1, freq_samp) == 0 )         
        
        for k = 1:n_sample
            G_fft(:,:,k) = G_fft(:,:,k) + v * fft_coeff(i_ftc,k);
        end
        i_ftc = i_ftc + 1;
        
    end 
    
    
    %- record seismograms -------------------------------------------------
    for ir = 1:n_receivers
        seismograms(ir,n) = v(rec_id(ir,1), rec_id(ir,2));
    end  
        
    
end


%- return the time reversed Green function --------------------------------
for k = 1:n_sample
    G_fft(:,:,k) = conj(G_fft(:,:,k));
end


end

