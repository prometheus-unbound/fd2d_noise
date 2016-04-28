

clear all
clc


% x_center_ring = 3.0e4;
% z_center_ring = 3.0e4;
% radius = 2.0e4;
% thickness = 5e3;
% taper_strength = 1e7;

% x_center_ring = 1.0e6;
% z_center_ring = 1.0e6;
% radius = 6.9e5;
% thickness = 1.0e5;
% taper_strength = 10e8;

x_center_ring = 1.0e6;
z_center_ring = 1.0e6;
radius = 6.4e5;
thickness = 2.0e5;
taper_strength = 70e8;


[Lx, Lz, nx, nz, ~, ~, ~, model_type] = input_parameters();
[X, Z] = define_computational_domain(Lx,Lz,nx,nz);

R = ( (X-x_center_ring).^2 + (Z-z_center_ring).^2 ).^(1/2);
angle = atan( abs( X-x_center_ring ) ./ abs( Z-z_center_ring ) ) *180/pi;

[k,l] = find(isnan(angle));
angle(k,l) = 0;



% noise_source_distribution = ones(100,100);
noise_source_distribution = exp( -abs( R-radius ).^2 / taper_strength ) .* double(R > (radius-thickness/2) & R < (radius+thickness/2) );
% noise_source_distribution = double(R > (radius-thickness/2) & R < (radius+thickness/2) );

% usr_par.kernel.imfilter = fspecial('gaussian',[20 20], 10);
% noise_source_distribution = imfilter( noise_source_distribution, usr_par.kernel.imfilter, 'circular' );


mesh( X, Z, noise_source_distribution' )
colormap('autumn')
hold on

load ../output/interferometry/array_16_ref_center2.mat

level = 3;
plot3(array(:,1),array(:,2),level+0*array(:,2),'x')

[width] = absorb_specs();
level = [level,level];
plot3([width,Lx-width],[width,width],level,'k--')
plot3([width,Lx-width],[Lz-width,Lz-width],level,'k--')
plot3([width,width],[width,Lz-width],level,'k--')
plot3([Lx-width,Lx-width],[width,Lz-width],level,'k--')

axis square
view([0 90])

