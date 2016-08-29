
function [grad_m] = map_gradparameters_to_gradm( m, grad_parameters, usr_par )
% MAP_GRADPARAMETERS_TO_GRADM function to transform the partial derivatives
% w.r.t. the physical parameters to the partial derivatives
% w.r.t. to the model variable m. This requires to apply the chain rule
% which is implicitly given by the function MAP_M_TO_PARAMETERS.
%
% Input:
% m
% grad_parameters
% usr_par : auxiliary user defined parameters (optional)
%
% Output:
% grad_m
%
% See also MAP_M_TO_PARAMETERS.


temp = 0.0 * grad_parameters;

[absbound] = init_absbound();
[ix,iz] = find( absbound == 1, 1, 'first' );


if( strcmp( usr_par.ring.switch, 'yes' ) )
    
    [Lx, Lz] = input_parameters();
    [X, Z] = define_computational_domain( Lx, Lz, usr_par.config.nx, usr_par.config.nz );
    
    R = ( ( X - usr_par.ring.x_center_ring ).^2 + ( Z - usr_par.ring.z_center_ring ).^2 ).^(1/2);
    
    temp(:,:,1) = grad_parameters(:,:,1) .*  exp( -abs( R - usr_par.ring.radius ).^2 / usr_par.ring.taper_strength ) .* double(R > (usr_par.ring.radius-usr_par.ring.thickness/2) & R < (usr_par.ring.radius+usr_par.ring.thickness/2) );
    temp(ix:end-ix+1, iz:end-iz+1, end ) = imfilter( grad_parameters(ix:end-ix+1, iz:end-iz+1, end), usr_par.kernel.imfilter.structure, 'circular' );
    
else
    
    % temp(ix:end-ix+1, iz:end-iz+1, 1:end-1) = imfilter( grad_parameters(ix:end-ix+1, iz:end-iz+1, 1:end-1), usr_par.kernel.imfilter.source, 'circular' );
    % temp(ix:end-ix+1, iz:end-iz+1, end) = imfilter( grad_parameters(ix:end-ix+1, iz:end-iz+1, end), usr_par.kernel.imfilter.structure, 'circular' );
    
    [Lx, Lz] = input_parameters();
    [~, ~, x, z] = define_computational_domain(Lx, Lz, usr_par.config.nx, usr_par.config.nz);
    
    for i = 1:size(temp,3)-1
        temp(ix:end-ix+1, iz:end-iz+1, i) = gaussblur2d( grad_parameters(ix:end-ix+1, iz:end-iz+1, i), x(ix:end-ix+1), z(iz:end-iz+1), usr_par.kernel.sigma.source );
    end
    temp(ix:end-ix+1, iz:end-iz+1, end) = gaussblur2d( grad_parameters(ix:end-ix+1, iz:end-iz+1, end), x(ix:end-ix+1), z(iz:end-iz+1), usr_par.kernel.sigma.structure );
     
end

temp(:,:,end) = usr_par.initial.mu_0 * temp(:,:,end);

grad_m = reshape(temp, [], 1);


end