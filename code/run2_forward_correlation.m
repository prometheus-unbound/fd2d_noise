function [ seismograms, C_out ] = run2_forward_correlation( structure, noise_source, G_fft, rec, mode )

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
% rec: receivers
% mode: integer switch
%       == 0 do not save wavefield
%       == 1 save wavefield
%
% output:
%--------
% correlation recordings
% C: correlation wavefield
%
%==========================================================================


%- basic configuration ----------------------------------------------------
[Lx, Lz, nx, nz, dt, nt, order, ~, ~, store_fwd_nth] = input_parameters();
[~,~,x,z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);


%- time axis --------------------------------------------------------------
t = -(nt-1)*dt:dt:(nt-1)*dt;
nt = length(t);


%- prepare coefficients for Fourier transform -----------------------------
[~,n_sample,w_sample,dw,freq_samp] = input_interferometry();

ifft_coeff = zeros(nt,n_sample) + 1i*zeros(nt,n_sample);
for n=1:nt
    for k = 1:n_sample
        ifft_coeff(n,k) = 1/sqrt(2*pi) * exp( 1i*w_sample(k)*t(n) ) * dw;
    end
end


%- compute indices for receiver locations ---------------------------------
n_receivers = size(rec,1);
rec_id = zeros(n_receivers,2);

for i=1:n_receivers    
    rec_id(i,1) = min( find( min(abs(x-rec(i,1))) == abs(x-rec(i,1)) ) );
    rec_id(i,2) = min( find( min(abs(z-rec(i,2))) == abs(z-rec(i,2)) ) );   
end


%- initialise absorbing boundary taper a la Cerjan ------------------------
[absbound] = init_absbound();


%- allocate dynamic fields ------------------------------------------------
v = zeros(nx,nz);
sxy = zeros(nx-1,nz);
szy = zeros(nx,nz-1);
u = zeros(nx,nz);

n_fw = floor(nt/store_fwd_nth);
if( mode ~= 0 )
    C_out = zeros(nx,nz,n_fw,'single');
    % C_out = zeros(nx,nz,n_fw);
else
    C_out = single([]);
end


%- initialise seismograms -------------------------------------------------
seismograms = zeros(n_receivers,nt);


%==========================================================================
% iterate
%==========================================================================

%%% TEST 2 %%%
% for ns = 1:n_noise_sources
%     for k = 1:n_sample
%         distribution(:,:,k) = distribution(:,:,k) + spectrum(k,ns) * distribution(:,:,k);
%     end
% end
% 
% x_time_ifft = fftshift( ifft( ifftshift( distribution ) * nt, [], 3 ) );


% for ix = 1:nx
%     for iz = 1:nz
%         test_G_in(ix, iz, :) = conv( squeeze(test_G_in( ix, iz, : )), squeeze(S( ix, iz, : )), 'same' );
%     end
% end
%%% END TEST 2 %%%


i_fw_out = 1;

for n = 1:nt

    
    %- save correlation wavefield -----------------------------------------
    if( mode ~= 0 && mod(n, store_fwd_nth) == 0 )
        C_out(:,:,i_fw_out) = single( u );
        % C_out(:,:,i_fw_out) = u;
        i_fw_out = i_fw_out + 1;
    end

    
    %- compute divergence of current stress tensor ------------------------    
    DS = div_s(sxy,szy,dx,dz,nx,nz,order);
          
    
    %- add source of the correlation field --------------------------------
    % if( mod(n,freq_samp) == 0 && t(n) <= 0.0 )
    if( mod(n,freq_samp) == 0 )
        
        %- transform on the fly to the time domain
        S = zeros(nx,nz) + 1i*zeros(nx,nz);
        
        % calculate source for correlation wavefield
        for k = 1:n_sample
            S = S + noise_source.spectrum(k) * noise_source.distribution .* conj(G_fft(:,:,k)) * ifft_coeff(n,k);
        end
        
        DS = DS + real(S);
        
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

    
    %- calculate displacement ---------------------------------------------
    u = u + v * dt;
    
    
    %- record seismograms -------------------------------------------------
    for ir = 1:n_receivers
        seismograms(ir,n) = u(rec_id(ir,1), rec_id(ir,2));
    end  
    
    
end


end
