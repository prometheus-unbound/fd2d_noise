
% load stress field of correlation wavefield
forward_1 = load('stress_correlation_2015-10-05-17-12-42.mat');
forward_2 = load('stress_correlation_2015-10-05-17-13-50.mat');

% load adjoint state
adjoint_1 = load('stress_adjoint_state_2015-10-05-17-13-06.mat');
adjoint_2 = load('stress_adjoint_state_2015-10-05-17-14-17.mat');

% reverse adjoint state
adjoint_1.G_1_dxu_time_flip = flip( adjoint_1.G_1_dxu_time , 3 );
adjoint_1.G_1_dzu_time_flip = flip( adjoint_1.G_1_dzu_time , 3 );

adjoint_2.G_1_dxu_time_flip = flip( adjoint_2.G_1_dxu_time , 3 );
adjoint_2.G_1_dzu_time_flip = flip( adjoint_2.G_1_dzu_time , 3 );



% get frequency configuration
[Lx,Lz,nx,nz,dt,nt,order,model_type,source_type,n_basis_fct] = input_parameters();
t = -(nt-1)*dt:dt:(nt-1)*dt;
[f_sample,n_sample,w_sample,dw,nt_freq] = input_interferometry();


fft_coeff = zeros(length(t),n_sample) + 1i*zeros(length(t),n_sample);
ifft_coeff = zeros(length(t),n_sample) + 1i*zeros(length(t),n_sample);
for k = 1:n_sample
    fft_coeff(:,k) = 1/sqrt(2*pi) * exp( -1i*w_sample(k)*t' ) * dt;
    ifft_coeff(:,k) = 1/sqrt(2*pi) * exp( 1i*w_sample(k)*t' ) * dw;
end

[~,noise_spectrum] = make_noise_source();



% K_mu = zeros(nx,nz);
% for n = 1:length(t)
% 
%     % contribution of run 1
%     K_mu(1:nx-1,:) = K_mu(1:nx-1,:) - forward_1.C_2_dxu_time(:,:,n) .* adjoint_1.G_1_dxu_time_flip(:,:,n);
%     K_mu(:,1:nz-1) = K_mu(:,1:nz-1) - forward_1.C_2_dzu_time(:,:,n) .* adjoint_1.G_1_dzu_time_flip(:,:,n);
%     
%     % contribution of run 2
%     K_mu(1:nx-1,:) = K_mu(1:nx-1,:) - forward_2.C_2_dxu_time(:,:,n) .* adjoint_2.G_1_dxu_time_flip(:,:,n);
%     K_mu(:,1:nz-1) = K_mu(:,1:nz-1) - forward_2.C_2_dzu_time(:,:,n) .* adjoint_2.G_1_dzu_time_flip(:,:,n);
%     
% %     for k = 1:n_sample
% % 
% %         K_mu(1:nx-1,:) = K_mu(1:nx-1,:) + noise_spectrum(k,1) * conj( adjoint_1.G_1_dxv(:,:,k) ) .* forward_1.C_2_dxv(:,:,k);
% %         K_mu(:,1:nz-1) = K_mu(:,1:nz-1) + noise_spectrum(k,1) * conj( adjoint_1.G_1_dzv(:,:,k) ) .* forward_1.C_2_dzv(:,:,k);
% % 
% %         K_mu(1:nx-1,:) = K_mu(1:nx-1,:) + noise_spectrum(k,1) * conj( adjoint_2.G_1_dxv(:,:,k) ) .* forward_2.C_2_dxv(:,:,k);
% %         K_mu(:,1:nz-1) = K_mu(:,1:nz-1) + noise_spectrum(k,1) * conj( adjoint_2.G_1_dzv(:,:,k) ) .* forward_2.C_2_dzv(:,:,k);
% % 
% %     end
% 
% end

K_mu = zeros(nx,nz) + 1i*zeros(nx,nz);
for k=1:n_sample
    
    K_mu(1:nx-1,:) = K_mu(1:nx-1,:) + noise_spectrum(k,1) * conj(adjoint_1.G_1_dxv(:,:,k)) .* forward_1.C_2_dxv(:,:,k);% / w_sample(k)^2;
    K_mu(:,1:nz-1) = K_mu(:,1:nz-1) + noise_spectrum(k,1) * conj(adjoint_1.G_1_dzv(:,:,k)) .* forward_1.C_2_dzv(:,:,k);% / w_sample(k)^2;
    
    K_mu(1:nx-1,:) = K_mu(1:nx-1,:) + noise_spectrum(k,1) * conj(adjoint_2.G_1_dxv(:,:,k)) .* forward_2.C_2_dxv(:,:,k);% / w_sample(k)^2;
    K_mu(:,1:nz-1) = K_mu(:,1:nz-1) + noise_spectrum(k,1) * conj(adjoint_2.G_1_dzv(:,:,k)) .* forward_2.C_2_dzv(:,:,k);% / w_sample(k)^2;
    
end

grad_parameters = real(K_mu);
g = map_gradparameters_to_gradm( m, grad_parameters, usr_par );

