function [ K, stf_fft, u_out, seismograms] = run_noise_adjoint( mu, rho, u_fwd, stf, src, rec, spectrum, source_dist, G_fft, mode, dmu, u_in )

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
%       == 10 Hessian mode: do not save nor Fourier transform adjoint wavefield
%       == 11 Hessian mode: only save adjoint wavefield
%       == 12 Hessian mode: only Fourier transform adjoint wavefield
%       == 13 Hessian mode: save adjoint wavefield and Fourier transform adjoint wavefield

%
% output:
%--------
% K_mu: sensitivity kernel for mu
% stf_fft: Fourier transformed adjoint_state - necessary for second run
% u_out: adjoint wavefield
% displacement seismograms
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
ns_adj = size(src,1);
src_id = zeros(ns_adj,2);

for i=1:ns_adj
    src_id(i,1) = min( find( min(abs(x-src(i,1))) == abs(x-src(i,1))) );
    src_id(i,2) = min( find( min(abs(z-src(i,2))) == abs(z-src(i,2))) );
end


%- compute indices for receiver locations ---------------------------------
n_receivers = size(rec,1);
rec_id = zeros(n_receivers,2);

for i=1:n_receivers    
    rec_id(i,1) = min( find( min(abs(x-rec(i,1))) == abs(x-rec(i,1)) ) );
    rec_id(i,2) = min( find( min(abs(z-rec(i,2))) == abs(z-rec(i,2)) ) );   
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
if( ~isempty( source_dist ) )
    n_noise_sources = size(spectrum,2);
    
    if( n_basis_fct == 0 )
        source_distribution = reshape(source_dist, nx, nz, n_noise_sources);
    else
        source_distribution = reshape(source_dist, nx, nz, n_basis_fct);
    end
else
    n_noise_sources = 0;
    source_distribution = 0.0;
end


%- allocate dynamic fields ------------------------------------------------
v = zeros(nx,nz);
sxy = zeros(nx-1,nz);
szy = zeros(nx,nz-1);
u = zeros(nx,nz);
stf_fft = zeros(nx,nz,n_sample) + 1i*zeros(nx,nz,n_sample);

if( mode == 11 || mode == 13 )
    n_fw = floor(nt/fw_nth);
    u_out = zeros(nx,nz,n_fw,'single');
    % u_out = zeros(nx,nz,n_fw);
else
    u_out = single(0.0);
    % u_out = 0.0;
end


%- initialise seismograms -------------------------------------------------
seismograms = zeros(n_receivers,nt);


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
% !!! however, need the full run for hessian vector product
if( ( mode == 0 || mode == 1 ) && size(stf,3) ~= 1 )
    nt = n_zero;
end

i_ftc = 1;
i_fw_in = 1;
i_fw_out = 1;


%%% TEST %%%
if( n_basis_fct == 0 )
    interpol_matrix = ones( n_sample, 1 );
else
    interpol_matrix = zeros( n_sample, n_basis_fct );
end

if( n_basis_fct ~= 0 )
    
    for k = 1:n_sample
        
        ib = find( k >= int_limits(:,1) & k <= int_limits(:,2), 1 );
        indices = int_limits(ib,1) : int_limits(ib,2);
        ni = length(indices);
        
        current = ib;
        
        if( all(ib == 1) )
            below = [];
            above = 2;
        elseif( all(ib == n_basis_fct) )
            below = ib-1;
            above = [];
        else
            below = ib - 1;
            above = ib + 1;
        end
        
        i = find( indices == k, 1 );
        if( all(i < (ni)/2) )
            interpol_matrix( k, current ) = 1 - ((ni)/2-i)/ni;
            if( ~isempty(below) )
                interpol_matrix( k, below ) = ((ni)/2-i)/ni;
            end
        elseif( all(i == (ni)/2) )
            interpol_matrix( k, current ) = 1;
        elseif( all(i > (ni)/2) )
            interpol_matrix( k, current ) = 1 - (i-(ni)/2)/ni;
            if( ~isempty(above) )
                interpol_matrix( k, above ) = (i-(ni)/2)/ni;
            end
        end
        
    end
    
end

distribution = zeros( nx, nz, n_sample );

if( ~isscalar( source_distribution ) )
    for ix = 1:nx
        for iz = 1:nz
            distribution( ix, iz, : ) = interpol_matrix * squeeze( source_distribution( ix, iz, : ) );
        end
    end
end
%%% END TEST %%%


for n = 1:nt
    
    
    %- save adjoint field -------------------------------------------------
    if( ( mode == 11 || mode == 13 ) && mod(n, fw_nth) == 0 )
        u_out(:,:,i_fw_out) = single( u );
        % u_out(:,:,i_fw_out) = u;
        i_fw_out = i_fw_out + 1;
    end

    
    %- compute divergence of current stress tensor ------------------------    
    DS = div_s(sxy,szy,dx,dz,nx,nz,order); 
    
    
    %- add adjoint source time function -----------------------------------
    if( size(stf,3) == 1 && ~isempty(stf) )
        
        for i=1:ns_adj
            DS(src_id(i,1),src_id(i,2)) = DS(src_id(i,1),src_id(i,2)) + real(stf(i,n));
        end
        
    elseif( size(stf,3) ~= 1 && ~isempty(stf) )
        
        if( mod(n,freq_samp) == 0 ) % && t(n) < 0.0 )
            T = zeros(nx,nz) + 1i*zeros(nx,nz);
            for ns = 1:n_noise_sources
                for k = 1:n_sample
                    T = T + spectrum(k,ns) * distribution(:,:,k) .* conj( stf(:,:,k) ) * ifft_coeff(i_ftc,k);
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

    
    if( n == length(t) || isempty(dmu) || size( u_in, 3 ) ~= length(t) )
        sxy = sxy + dt * mu(1:nx-1,:) .* strain_dxv;
        szy = szy + dt * mu(:,1:nz-1) .* strain_dzv;
    else
        sxy = sxy + dt * mu(1:nx-1,:) .* strain_dxv - dmu(1:nx-1,:) .* dx_v( u_in(:,:,n) - u_in(:,:,n+1), dx, dz, nx, nz, order );
        szy = szy + dt * mu(:,1:nz-1) .* strain_dzv - dmu(:,1:nz-1) .* dz_v( u_in(:,:,n) - u_in(:,:,n+1), dx, dz, nx, nz, order );
    end
    
    
    %- calculate displacement and respective strain -----------------------
    u = u + v * dt;
    strain_dxu = dx_v(u,dx,dz,nx,nz,order);
    strain_dzu = dz_v(u,dx,dz,nx,nz,order);  
    
    
    %- build up source kernel ---------------------------------------------
    % if( ~isempty(G_fft) && size( G_fft, 3 ) == n_sample && mod(n,freq_samp) == 0 && t(n) >= 0.0 )
    if( ~isempty(G_fft) && size( G_fft, 3 ) == n_sample && mod(n,freq_samp) == 0 )
        
        if( n_basis_fct == 0 )
            M_tn = zeros(nx,nz,1) + 1i*zeros(nx,nz,1);
        else
            M_tn = zeros(nx,nz,n_basis_fct) + 1i*zeros(nx,nz,n_basis_fct);
        end
            
        for ix = 1:nx
            for iz = 1:nz
                M_tn(ix,iz,:) = interpol_matrix' * ( spectrum .* conj(squeeze(G_fft(ix,iz,:))) .* ifft_coeff2(n, :).' );
            end
        end
        
        for ib = 1:size(M_tn,3)
            K_s(:,:,ib) = K_s(:,:,ib) + real( M_tn(:,:,ib) .* u );
        end
        
    end
    
       
    %- build up structure kernel ------------------------------------------
    if( ~isempty( u_fwd ) && size( u_fwd, 3 ) >= i_fw_in && mod( n, fw_nth ) == 0 )
                    
        K_mu(1:nx-1,:) = K_mu(1:nx-1,:) - strain_dxu .* dx_v( u_fwd(:,:,end-i_fw_in+1), dx, dz, nx, nz, order ) * fw_nth;
        K_mu(:,1:nz-1) = K_mu(:,1:nz-1) - strain_dzu .* dz_v( u_fwd(:,:,end-i_fw_in+1), dx, dz, nx, nz, order ) * fw_nth;
        
        i_fw_in = i_fw_in + 1;

    end


    %- save Fourier transformed adjoint state for second run --------------
    if( ( mode == 1 || mode == 12 || mode == 13 ) )
        
        if( mod(n,freq_samp) == 0 )
            % if( t(n) >= 0.0 )
                for k = 1:n_sample
                    stf_fft(:,:,k) = stf_fft(:,:,k) + u * fft_coeff(i_ftc,k);
                end 
            % end
            i_ftc = i_ftc + 1;
        end
        
    end
    
        
    %- record seismograms -------------------------------------------------
    for ir = 1:n_receivers
        seismograms(ir,n) = u(rec_id(ir,1), rec_id(ir,2));
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

