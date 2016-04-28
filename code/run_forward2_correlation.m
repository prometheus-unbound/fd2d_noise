function [ seismograms, t, C_out ] = run_forward2_correlation( mu, rho, G_fft, spectrum, source_dist, rec, mode, dmu, C_in )

%==========================================================================
% compute correlation wavefield
%
% input:
%--------
% mu [N/m^2]
% rho [kg/m^3]
% G_fft: Fourier transformed Green function of reference station
% spectrum: spectrum of noise distribution
% source_dist: source distribution
% rec: receiverss
% mode: integer switch
%       == 0 do not save wavefield
%       == 1 save wavefield
%
% output:
%--------
% correlation recordings
% t: time vector
% C: correlation wavefield
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
n_noise_sources = size(spectrum, 2);

if( n_basis_fct == 0 )
    source_distribution = reshape(source_dist, nx, nz, n_noise_sources);
else
    source_distribution = reshape(source_dist, nx, nz, n_basis_fct);
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


%- time axis --------------------------------------------------------------
t = -(nt-1)*dt:dt:(nt-1)*dt;
nt = length(t);


%- prepare coefficients for Fourier transform and its inverse -------------
n_ftc = floor(nt/freq_samp);
fft_coeff = zeros(n_ftc,n_sample) + 1i*zeros(n_ftc,n_sample);
ifft_coeff = zeros(n_ftc,n_sample) + 1i*zeros(n_ftc,n_sample);
i_ftc = 1;
for n = 1:nt
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

n_fw = floor(nt/fw_nth);
if( mode ~= 0 )
    C_out = zeros(nx,nz,n_fw,'single');
else
    C_out = single(0.0);
end


%- initialise seismograms -------------------------------------------------
seismograms = zeros(n_receivers,nt);


%- initialise absorbing boundary taper a la Cerjan ------------------------
[absbound] = init_absbound();


%==========================================================================
% iterate
%==========================================================================

i_ftc = 1;
i_fw_out = 1;

for n = 1:nt

    
    %- save correlation wavefield -----------------------------------------
    if( mode ~= 0 && mod(n, fw_nth) == 0 )
        C_out(:,:,i_fw_out) = single( u );
        i_fw_out = i_fw_out + 1;
    end

    
    %- compute divergence of current stress tensor ------------------------    
    DS = div_s(sxy,szy,dx,dz,nx,nz,order);
      
    
    %- add source of the correlation field --------------------------------
    if( mod(n,freq_samp) == 0 && t(n) <= 0.0 )
        
        %- transform on the fly to the time domain        
        S = zeros(nx,nz,n_noise_sources) + 1i*zeros(nx,nz,n_noise_sources);
        
        % calculate source for correlation wavefield
        for ns = 1:n_noise_sources
            
            for k = 1:n_sample
                
                if( n_basis_fct ~= 0 )
                    ib = find( k >= int_limits(:,1) & k <= int_limits(:,2) );
                else
                    ib = ns;
                end
                
                S(:,:,ns) = S(:,:,ns) + spectrum(k,ns) * source_distribution(:,:,ib) .* G_fft(:,:,k) * ifft_coeff(i_ftc,k);
                
            end
            
            DS = DS + real(S(:,:,ns));
            
        end
        
        i_ftc = i_ftc + 1;

    end
       
    
    %- update velocity field ----------------------------------------------
    v = v + dt * DS./rho;
   
    
    %- apply absorbing boundary taper -------------------------------------    
    v = v .* absbound;
    
    
    %- compute derivatives of current velocity and update stress tensor ---
    strain_dxv = dx_v(v,dx,dz,nx,nz,order);
    strain_dzv = dz_v(v,dx,dz,nx,nz,order);
    
    if( n==nt || isempty(dmu) )
        sxy = sxy + dt * mu(1:nx-1,:) .* strain_dxv;
        szy = szy + dt * mu(:,1:nz-1) .* strain_dzv;
    else
        sxy = sxy + dt * mu(1:nx-1,:) .* strain_dxv - dmu(1:nx-1,:) .* dx_v( C_in(:,:,n) - C_in(:,:,n+1), dx, dz, nx, nz, order );
        szy = szy + dt * mu(:,1:nz-1) .* strain_dzv - dmu(:,1:nz-1) .* dz_v( C_in(:,:,n) - C_in(:,:,n+1), dx, dz, nx, nz, order );
    end
    
    
    %- calculate displacement ---------------------------------------------
    u = u + v * dt;
    
    
    %- record seismograms -------------------------------------------------
    for ir = 1:n_receivers
        seismograms(ir,n) = u(rec_id(ir,1), rec_id(ir,2));
    end  
    
    
end


end

