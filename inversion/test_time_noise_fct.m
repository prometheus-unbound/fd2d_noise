

[Lx,Lz,nx,nz,dt,nt,order,model_type,source_type,n_basis_fct] = input_parameters();
t = -(nt-1)*dt:dt:(nt-1)*dt;
[f_sample,n_sample,w_sample,dw,nt_freq] = input_interferometry();


fft_coeff = zeros(length(t),n_sample) + 1i*zeros(length(t),n_sample);
ifft_coeff = zeros(length(t),n_sample) + 1i*zeros(length(t),n_sample);
for k = 1:n_sample
    fft_coeff(:,k) = 1/sqrt(2*pi) * exp( -1i*w_sample(k)*t' ) * dt;
    ifft_coeff(:,k) = 1/sqrt(2*pi) * exp( 1i*w_sample(k)*t' ) * dw;
end


[noise_source_distribution,noise_spectrum] = make_noise_source('no');


S = zeros(nx,nz,length(t)) + 1i * zeros(nx,nz,length(t));
for n = 1:length(t)
    
    n
    
    for k = 1:n_sample
        
        S(:,:,n) = S(:,:,n) + noise_spectrum(k,1) * noise_source_distribution .* conj(G_2(:,:,k)) * ifft_coeff(n,k);
        
    end
    
end
            
S = real(S);


for n = 1:length(t)
   
    mesh(S(:,:,n))
    drawnow
    pause(0.01)
    
end