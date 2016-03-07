

clear all
close all
clc


%% INPUT PARAMETERS
[Lx, Lz, nx, nz] = input_parameters();
[X,Z,x,z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);


corr_length = 0.25;
k_max = 2 * pi * dx / (corr_length * Lx);
pert = 1e10 / (corr_length * Lx)^2;


%% TESTING

% step 1: generate a random phase spectrum [-pi, pi]
phase_random = pi * randn( nx, nz );


% step 2: modulate complex phase spectrum with predefined amplitude spectrum to define the Fourier spectrum
amplitude = ones( nx, nz );

k_backup = ones( nx, nz );
for i = 1:nx
    
    for j = 1:nz
        
        kx = 2 * pi * dx * (i-1) / Lx;
        kz = 2 * pi * dz * (j-1) / Lz;
        
        k = sqrt( kx^2 + kz^2 );
        k_backup(i,j) = k;
        
        if( k > k_max )
            amplitude(i,j) = 0.0;
        end
        
    end
    
end

signal_fft = amplitude .* exp( 1i * phase_random );
signal_fft_normal = exp( 1i * phase_random );


% step 3: inverse FFT to obtain space domain representation
signal = real( ifft2( signal_fft ) );
signal_normal = real( ifft2( signal_fft_normal ) );


% step 4: set to zero in absorbing boundary region and smooth out
[absbound] = init_absbound();
signal = signal .* double( absbound == 1 );

[width] = absorb_specs();
absbound = ones(nx,nz);
absbound = absbound .* ( double([X'>2*width]) + cos( (X'-2*width)/width * pi/2 ) .* double([X'<=2*width]) );
absbound = absbound .* ( double([X'<(Lx-2*width)]) + cos( (X'-(Lx-2*width))/width * pi/2 ) .* double([X'>=(Lx-2*width)]) );
absbound = absbound .* ( double([Z'>2*width]) + cos( (Z'-2*width)/width * pi/2 ) .* double([Z'<=2*width]) );
absbound = absbound .* ( double([Z'<(Lz-2*width)]) + cos( (Z'-(Lz-2*width))/width * pi/2 ) .* double([Z'>=(Lz-2*width)]) );
    
signal = absbound .* signal;


% step 5: normalization
rms_signal = rms( rms( signal ) );
max_signal = max( max( abs( signal ) ) );

signal = pert * signal / max_signal;



%% PLOTTING
fig1 = figure(1);
set(fig1,'units','normalized','position',[.1 .6 0.2 0.3])
hold on
mesh( X, Z, signal' )
level = [3 3];
plot3([width,Lx-width],[width,width],level,'k--')
plot3([width,Lx-width],[Lz-width,Lz-width],level,'k--')
plot3([width,width],[width,Lz-width],level,'k--')
plot3([Lx-width,Lx-width],[width,Lz-width],level,'k--')

view([0 90])
axis square
colorbar

% fig2 = figure(2);
% hold on
% set(fig2,'units','normalized','position',[.1 .2 0.2 0.3])
% mesh( X, Z, signal_normal' )
% view([0 90])
% axis image
% colorbar

